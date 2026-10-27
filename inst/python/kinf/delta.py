"""Exact incremental hard NLL for single transpositions (for MCMC / greedy search on large graphs).

A wiring is held as b: for each left inner end in var_a (rows of all variable nodes, concatenated), the right
inner end it is joined to. A move swaps b[i], b[j] for two rows of the same node. Only the walks through those
two rows change, and the set S of tile copies on them is the same before and after the swap. Tile copies on
different walks have no contacts (beyond the constant interchromosomal background), so only copy pairs in
S x S change, hence only expected counts on tile pairs within tiles(S): the update costs O(|S|^2) instead of
O(ntc^2). S is gathered into a fixed-capacity buffer chosen by lax.switch over a ladder of capacities.
State: dict(b, info=(coord, wid, L, iscyc) per end, mu = full symmetric expected map, nll).
"""
import numpy as np
import jax
import jax.numpy as jnp
import model as km
import hic_model as hm

def make_delta(problem, bigdist=1e12, dcut=1e10, theta=2.0, weighted=True, caps=(64, 128, 256, 512, 1024, 2048, 4096)):
    P = problem
    wd = km.make_walkdist_func(P, bigdist)
    walkinfo, var_a = wd.walkinfo, wd.var_a
    K = hm.get_constants()
    uid = P.tc_uid; ntc = len(uid); ntile = P.ntile
    Lu = jnp.asarray(2 * uid); Ru = jnp.asarray(2 * uid + 1); uidj = jnp.asarray(uid)
    oL = jnp.asarray(P.offL); oR = jnp.asarray(P.offR)
    tc_w = jnp.asarray(P.tc_w); tile_of = jnp.asarray(P.tc_tile)
    tw = jnp.asarray(P.tile_w); tcn = jnp.asarray(P.tile_cn)
    hd = 0.5 * P.depth
    area0 = np.median(P.tile_w) ** 2
    wp = jnp.asarray(P.tile_w[:, None] * P.tile_w[None, :] / area0) if weighted else jnp.ones((ntile, ntile))
    x = jnp.asarray(P.x_obs)
    tri = jnp.asarray(np.triu(np.ones((ntile, ntile), bool)))
    caps = [c for c in sorted(caps) if c < ntc] + [ntc]
    caps_arr = jnp.asarray(caps)

    # moves: all transpositions (i, j) of rows within one node block of var_a
    nodes = sorted(P.var_nodes)
    cn = [len(P.var_nodes[n][0]) for n in nodes]
    start = np.concatenate([[0], np.cumsum(cn)])[:-1]
    moves = np.array([(s + a, s + c) for s, k in zip(start, cn) for a in range(k) for c in range(a + 1, k)]).reshape(-1, 2)
    move_node = np.repeat(np.arange(len(nodes)), [k * (k - 1) // 2 for k in cn])
    var_b = {n: np.asarray(km.inner_end(1, P.var_nodes[n][1])) for n in nodes}

    def end_dist(info, e1, e2):
        coord, wid, L, iscyc = info
        d = jnp.abs(coord[e1] - coord[e2])
        d = jnp.where(iscyc[e1], jnp.minimum(d, L[e1] - d), d)
        d = jnp.where(wid[e1] == wid[e2], d, bigdist)
        return jnp.where(e1 == e2, 0.0, d)

    def tc_dist(info, I, J):  # = refine(walkdist(.)) restricted to tile-copy pairs (I, J)
        m = jnp.minimum(jnp.minimum(end_dist(info, Lu[I], Lu[J]) + oL[I] + oL[J], end_dist(info, Ru[I], Ru[J]) + oR[I] + oR[J]),
                        jnp.minimum(end_dist(info, Lu[I], Ru[J]) + oL[I] + oR[J], end_dist(info, Ru[I], Lu[J]) + oR[I] + oL[J]))
        d = jnp.where(uidj[I] == uidj[J], jnp.abs(oL[I] - oL[J]), m)
        return jnp.where(I == J, 0.0, d)

    def full_mu(info):
        a = jnp.arange(ntc)
        return hm.contact_model(tc_dist(info, a[:, None], a[None, :]), tc_w, tile_of, tw, tcn, P.depth, dcut, K)

    def nll_of(mu):
        return -jnp.sum(jnp.where(tri, wp * km.nb_logpmf(x, mu, theta), 0.0))

    @jax.jit
    def init(b):
        info = walkinfo(b); mu = full_mu(info)
        return dict(b=b, info=info, mu=mu, nll=nll_of(mu))

    def branch(c):
        def f(args):
            mask, cnt, info, info2, mu = args
            S = jnp.nonzero(mask, size=c, fill_value=0)[0]
            valid = jnp.arange(c) < cnt
            I, J = S[:, None], S[None, :]
            wS = tc_w[S]
            dM = hm.contact_matrix(tc_dist(info2, I, J), wS, dcut, K) - hm.contact_matrix(tc_dist(info, I, J), wS, dcut, K)
            dM = jnp.where(valid[:, None] & valid[None, :], dM, 0.0)
            tS = jnp.where(valid, tile_of[S], ntile)
            T = jnp.unique(tS, size=c, fill_value=ntile)
            loc = jnp.searchsorted(T, tS)
            dmu = hd * jnp.zeros((c, c)).at[loc[:, None], loc[None, :]].add(dM, mode='drop')
            p, q = T[:, None], T[None, :]
            ok = (p < ntile) & (q < ntile) & (p <= q)
            pc, qc = jnp.minimum(p, ntile - 1), jnp.minimum(q, ntile - 1)
            mo = mu[pc, qc]; xx = x[pc, qc]
            dnll = -jnp.sum(jnp.where(ok, wp[pc, qc] * (km.nb_logpmf(xx, mo + dmu, theta) - km.nb_logpmf(xx, mo, theta)), 0.0))
            return dnll, mu.at[p, q].add(dmu, mode='drop')
        return f
    branches = [branch(c) for c in caps]

    def propose(state, i, j):
        """exact proposal for swapping rows i, j: returns (dnll, new state)"""
        b = state['b']; info = state['info']
        b2 = b.at[i].set(b[j]).at[j].set(b[i])
        info2 = walkinfo(b2)
        wid = info[1]; tcw = wid[Lu]
        mask = (tcw == wid[var_a[i]]) | (tcw == wid[var_a[j]])
        cnt = mask.sum()
        k = jnp.searchsorted(caps_arr, cnt)
        dnll, mu2 = jax.lax.switch(k, branches, (mask, cnt, info, info2, state['mu']))
        return dnll, dict(b=b2, info=info2, mu=mu2, nll=state['nll'] + dnll)

    def affected(state, i, j):  # |S| for diagnostics
        wid = state['info'][1]; tcw = wid[Lu]
        return ((tcw == wid[var_a[i]]) | (tcw == wid[var_a[j]])).sum()

    def b_from_perms(perms):
        return wd.cols({n: jnp.asarray(p) for n, p in perms.items()})

    def perms_from_b(b):
        b = np.asarray(b); out = {}
        for n, s, k in zip(nodes, start, cn):
            pos = {e: c for c, e in enumerate(var_b[n])}
            M = np.zeros((k, k)); M[np.arange(k), [pos[e] for e in b[s:s + k]]] = 1
            out[n] = M
        return out

    return dict(init=init, propose=propose, affected=affected, nll_of=nll_of, full_mu=full_mu,
                moves=moves, move_node=move_node, nodes=nodes, cn=np.array(cn), caps=caps,
                b_from_perms=b_from_perms, perms_from_b=perms_from_b)
