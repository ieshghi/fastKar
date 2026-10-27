"""Hard-wiring forward model + NB loss using the current fastKar Hi-C model (see hic_model.py).

Encoding (from R/prep_data_for_py.R): every copy of every half-node has a unique id (uid, 0-based here).
A left half has its OUTER end at its genomic left end and its INNER end at the right end; a right half the
opposite. Junction copies (external edges) join OUTER ends of two uids; the free variable at a node with CN > 1
is a CN x CN 0/1 permutation joining INNER ends of left copies (rows) to right copies (cols).
End indexing in the distance graph: end 2*uid = genomic left end of the half, 2*uid+1 = genomic right end.
The Sinkhorn-relaxed version of this model is archived in archive/sinkhorn/model.py.
"""
import json
from types import SimpleNamespace
import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
from jax.scipy.special import gammaln
import hic_model as hm

#### problem loading

def load_problem(folder):
    folder = folder.rstrip('/') + '/'
    rd = lambda n: pd.read_csv(folder + n + '.csv', index_col=0)
    nodes = rd('nodes').sort_values('unique.id')
    uid = nodes['unique.id'].values - 1
    assert (uid == np.arange(len(uid))).all()
    nu = len(uid)
    node_of_uid = nodes['node_id'].values.astype(int)
    side_of_uid = (nodes['side'].values == 'right').astype(int)  # 0 left, 1 right
    width_uid = nodes['width'].values.astype(float)

    ie = rd('internal_edges')
    ee = rd('external_edges')
    ext = ee[['n1', 'n2']].values.astype(int) - 1
    inner_fixed = ie[ie['copies'] == 1][['left', 'right']].values.astype(int) - 1  # CN 1: fixed pairing
    var_nodes = {}
    for node, g in ie[ie['copies'] > 1].groupby('node_id'):
        var_nodes[int(node)] = (g['left'].values.astype(int) - 1, g['right'].values.astype(int) - 1)

    tt = rd('tiled_target')
    dt = rd('dedup_tiles')
    ntile = len(tt)
    tc_uid = dt['unique.id'].values.astype(int) - 1
    tc_tile = dt['tile.id'].values.astype(int) - 1
    tc_w = dt['width'].values.astype(float)
    offL = dt['leftdist'].values.astype(float)
    offR = dt['rightdist'].values.astype(float)
    tile_w = tt['width'].values.astype(float)
    tile_cn = np.bincount(tc_tile, minlength=ntile).astype(float)  # graph CN of the tile = number of copies

    gm = rd('gm_dat')
    x_obs = np.zeros((ntile, ntile))
    x_obs[gm['i'].values - 1, gm['j'].values - 1] = np.round(gm['value'].values)  # upper triangle (i <= j)
    assert (gm['i'].values <= gm['j'].values).all()
    depth = float(pd.read_csv(folder + 'depth.csv')['x'].iloc[0])
    purity = float(pd.read_csv(folder + 'purity.csv')['x'].iloc[0])
    assert purity == 1.0, 'purity < 1 not ported'

    return SimpleNamespace(folder=folder, nu=nu, node_of_uid=node_of_uid, side_of_uid=side_of_uid,
                           width_uid=width_uid, ext=ext, inner_fixed=inner_fixed, var_nodes=var_nodes,
                           ntile=ntile, tc_uid=tc_uid, tc_tile=tc_tile, tc_w=tc_w, offL=offL, offR=offR,
                           tile_w=tile_w, tile_cn=tile_cn, x_obs=x_obs, depth=depth, nodes=nodes)

#### graph distances

def outer_end(side, uid):  # genomic end index at which junctions attach
    return 2 * uid + side

def inner_end(side, uid):  # end joined to the other half of the node
    return 2 * uid + (1 - side)

def make_costmat_func(problem, bigdist=1e12):
    """returns cost(perms): (2nu x 2nu) one-step cost matrix over half-node ends.
    internal half edge: width; junction and present permutation links: 0 (R adds nothing); absent: bigdist.
    perms: {node: (CN,CN) 0/1}. Only used by the Floyd-Warshall reference (hard_expected_map_fw)."""
    nu = problem.nu; N = 2 * nu
    side = problem.side_of_uid
    base = np.full((N, N), bigdist)
    np.fill_diagonal(base, 0.0)
    e = np.arange(nu)
    base[2 * e, 2 * e + 1] = problem.width_uid
    base[2 * e + 1, 2 * e] = problem.width_uid
    for n1, n2 in problem.ext:
        a = outer_end(side[n1], n1); b = outer_end(side[n2], n2)
        base[a, b] = base[b, a] = 0.0
    for l, r in problem.inner_fixed:
        a = inner_end(side[l], l); b = inner_end(side[r], r)
        base[a, b] = base[b, a] = 0.0
    nodes = sorted(problem.var_nodes)
    rows, cols = [], []
    for n in nodes:
        l, r = problem.var_nodes[n]
        rr, cc = np.meshgrid(inner_end(0, l), inner_end(1, r), indexing='ij')
        rows.append(rr.ravel()); cols.append(cc.ravel())
    rows = np.concatenate(rows) if rows else np.zeros(0, int)
    cols = np.concatenate(cols) if cols else np.zeros(0, int)
    base = jnp.asarray(base)
    def cost(perms):
        if len(nodes) == 0:
            return base
        p = jnp.concatenate([perms[n].ravel() for n in nodes])
        c = jnp.where(p > 0.5, 0.0, bigdist)
        return base.at[rows, cols].set(c).at[cols, rows].set(c)
    return cost

def floyd_warshall(D):
    """all-pairs shortest path (exact; O(N^3) reference for make_walkdist_func)"""
    N = D.shape[0]
    idx = jnp.arange(N)
    def body(k, D):
        new = D[:, k][:, None] + D[k, :][None, :]
        upd = jnp.minimum(D, new)
        keep = (idx[:, None] == k) | (idx[None, :] == k)
        return jnp.where(keep, D, upd)
    D = jax.lax.fori_loop(0, N, body, D)
    return D.at[jnp.diag_indices(N)].set(0.0)

def make_walkdist_func(problem, bigdist=1e12):
    """exact end-end distances for HARD permutations without Floyd-Warshall: a hard wiring gives every end at
    most one link partner, so the genome is a set of disjoint walks (paths and cycles). Distance = length
    along the shared walk (shorter way round on a cycle), bigdist between walks. Walks are found by pointer
    jumping (log2 N rounds of gathers), so this is jittable and vmappable.
    State a = 'enter half-copy a//2 at end a, leave at end a^1'; succ(a) = link partner of a^1 (or sink N)."""
    nu = problem.nu; N = 2 * nu
    lp0 = np.full(N, -1)
    def link(a, b):
        assert a != b and lp0[a] < 0 and lp0[b] < 0, 'end with more than one link'
        lp0[a] = b; lp0[b] = a
    side = problem.side_of_uid
    for n1, n2 in problem.ext:
        link(outer_end(side[n1], n1), outer_end(side[n2], n2))
    for l, r in problem.inner_fixed:
        link(inner_end(side[l], l), inner_end(side[r], r))
    nodes = sorted(problem.var_nodes)
    # rows: inner ends of left copies, cols: inner ends of right copies, concatenated over nodes (one scatter)
    var_a = jnp.asarray(np.concatenate([inner_end(0, problem.var_nodes[n][0]) for n in nodes] or [np.zeros(0, int)]))
    var_b = {n: jnp.asarray(inner_end(1, problem.var_nodes[n][1])) for n in nodes}
    lp0 = jnp.asarray(lp0)
    w = jnp.asarray(np.append(np.repeat(problem.width_uid, 2), 0.0))  # width of the half crossed in state a; sink 0
    K = int(np.ceil(np.log2(N + 1))) + 1
    idx = jnp.arange(N + 1)
    def jump(s, val, comb):  # K rounds of pointer doubling: val(a) <- comb(val(a), val(s(a))), s <- s(s)
        def body(_, c):
            s, v = c
            return s[s], comb(v, jax.tree.map(lambda x: x[s], v))
        return jax.lax.fori_loop(0, K, body, (s, val))
    def cols(perms):  # right inner end joined to each row of var_a
        return jnp.concatenate([var_b[n][jnp.argmax(perms[n], axis=1)] for n in nodes]) if nodes else jnp.zeros(0, int)
    def walkinfo(b):
        """b: right inner end joined to each left inner end in var_a. Returns per end: coord (bp along its walk),
        wid (walk id), L (cycle length), iscyc."""
        lp = lp0.at[var_a].set(b).at[b].set(var_a) if nodes else lp0
        nxt = lp[jnp.arange(N) ^ 1]
        succ = jnp.append(jnp.where(nxt >= 0, nxt, N), N)
        # pass 1: cycle detection and the minimum state of each cycle
        sK, mn = jump(succ, jnp.where(idx < N, idx, N + 1), jnp.minimum)
        cyc = sK != N
        # break each cycle in front of its minimum state, then every state lies on a path
        succ2 = jnp.where(cyc & (succ == mn), N, succ)
        # pass 2: R(a) = width from state a to the tail inclusive; T(a) = tail state of a's path (T(sink) = N)
        _, (R, T) = jump(succ2, (w, idx), lambda v, vs: (v[0] + vs[0], jnp.where(vs[1] == N, v[1], vs[1])))
        R, T, cyc = R[:N], T[:N], cyc[:N]
        circ = jax.ops.segment_max(R, T, num_segments=N)[T]  # cycle length (= R at the head)
        # per half-copy pick the direction with the smaller tail (consistent along each walk)
        s0 = 2 * jnp.arange(nu)
        c = jnp.where(T[s0] <= T[s0 + 1], s0, s0 + 1)  # canonical state per uid
        ce = c ^ 1
        coord = jnp.zeros(N).at[c].set(R[c]).at[ce].set(R[c] - w[c])  # entry end at R, exit end at R - width
        wid = jnp.zeros(N, int).at[c].set(T[c]).at[ce].set(T[c])
        L = jnp.zeros(N).at[c].set(circ[c]).at[ce].set(circ[c])
        iscyc = jnp.zeros(N, bool).at[c].set(cyc[c]).at[ce].set(cyc[c])
        return coord, wid, L, iscyc
    def walkdist(perms):
        coord, wid, L, iscyc = walkinfo(cols(perms))
        d = jnp.abs(coord[:, None] - coord[None, :])
        d = jnp.where(iscyc[:, None], jnp.minimum(d, L[:, None] - d), d)
        D = jnp.where(wid[:, None] == wid[None, :], d, bigdist)
        return D.at[jnp.diag_indices(N)].set(0.0)
    walkdist.walkinfo, walkdist.cols, walkdist.var_a = walkinfo, cols, var_a
    return walkdist

def make_refine_func(problem):
    """end-end distances -> tile-copy x tile-copy distances.
    For copies in different half-copies: min over the 4 end pairings of offset + end-end distance.
    Within the same half copy: |difference of offsets|. (Circular wrap inside one half is never shorter:
    a circle contains whole nodes so L >= 2 * half width.)"""
    uid = problem.tc_uid
    Lu = jnp.asarray(2 * uid); Ru = jnp.asarray(2 * uid + 1)
    oL = jnp.asarray(problem.offL); oR = jnp.asarray(problem.offR)
    same = jnp.asarray(uid[:, None] == uid[None, :])
    delta = jnp.abs(oL[:, None] - oL[None, :])
    n = len(uid)
    eye = jnp.eye(n, dtype=bool)
    def refine(Dend):
        LL = Dend[Lu][:, Lu] + oL[:, None] + oL[None, :]
        RR = Dend[Ru][:, Ru] + oR[:, None] + oR[None, :]
        LR = Dend[Lu][:, Ru] + oL[:, None] + oR[None, :]
        RL = Dend[Ru][:, Lu] + oR[:, None] + oL[None, :]
        st = jnp.stack([LL, RR, LR, RL])
        d = jnp.where(same, delta, jnp.min(st, axis=0))
        return jnp.where(eye, 0.0, d)
    return refine

#### likelihood

def nb_logpmf(x, mu, theta):
    mu_s = jnp.where(mu > 0, mu, 1.0)
    p = theta / (theta + mu_s)
    lp = gammaln(x + theta) - gammaln(theta) - gammaln(x + 1) + theta * jnp.log(p) + x * jnp.log1p(-p)
    return jnp.where(mu > 0, lp, 0.0)  # R compdats: logprob[y == 0] = 0

def make_nll_func(problem, theta=2.0, area0=None, weighted=True):
    """compmaps conventions: upper triangle (i <= j, incl. diagonal), weight widthprod/area0,
    area0 = median(tile width)^2, observed counts rounded. weighted=False gives every pixel weight 1, which is what
    R's compmaps() does on gM objects (gM drops the widthprod column; see scripts/ref_nll.R)."""
    w = problem.tile_w
    if area0 is None:
        area0 = np.median(w) ** 2
    tri = jnp.asarray(np.triu(np.ones((problem.ntile,) * 2, bool)))
    wp = jnp.asarray(w[:, None] * w[None, :]) / area0 if weighted else jnp.ones((problem.ntile,) * 2)
    x = jnp.asarray(problem.x_obs)
    def nll(Y):
        return -jnp.sum(jnp.where(tri, wp * nb_logpmf(x, Y, theta), 0.0))
    return nll

#### put it together

def make_model(problem, bigdist=1e12, dcut=1e10, theta=2.0, area0=None, weighted=True):
    """returns SimpleNamespace with
      hard_expected_map(perms)     full symmetric tile x tile expected counts (R map = upper triangle) for hard 0/1
                                   permutations {node: (CN,CN)}, by walk traversal
      hard_expected_map_fw(perms)  same by Floyd-Warshall (O(N^3) reference)
      walkdist(perms)              end x end distances by walk traversal
      nll(Y)                       weighted NB NLL (compmaps convention)
      loss_perms(perms)            nll(hard_expected_map(perms))"""
    costfn = make_costmat_func(problem, bigdist)
    refine = make_refine_func(problem)
    nll = make_nll_func(problem, theta, area0, weighted)
    w = jnp.asarray(problem.tc_w); tile_of = jnp.asarray(problem.tc_tile)
    tw = jnp.asarray(problem.tile_w); tcn = jnp.asarray(problem.tile_cn)
    K = hm.get_constants()

    def from_enddist(Dend):
        return hm.contact_model(refine(Dend), w, tile_of, tw, tcn, problem.depth, dcut, K)

    hard_expected_map_fw = lambda perms: from_enddist(floyd_warshall(costfn({n: jnp.asarray(p, dtype=jnp.float64) for n, p in perms.items()})))
    walkdist = make_walkdist_func(problem, bigdist)
    hard_expected_map = lambda perms: from_enddist(walkdist(perms))  # equals hard_expected_map_fw
    loss_perms = lambda perms: nll(hard_expected_map(perms))
    return SimpleNamespace(hard_expected_map=hard_expected_map, hard_expected_map_fw=hard_expected_map_fw,
                           walkdist=walkdist, loss_perms=loss_perms, nll=nll, problem=problem)

def upper_triangle(Y):
    i, j = np.triu_indices(Y.shape[0])
    return i, j, np.asarray(Y)[i, j]
