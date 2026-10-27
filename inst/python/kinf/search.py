"""MCMC search on the exact hard NLL with delta updates (delta.py): greedy (T = 0) or parallel tempering (PT).

Proposal: a node uniformly at random, then a uniform transposition of its rows (symmetric). Each chain runs
blocks of `block` steps on the device (lax.scan); PT swaps temperatures between adjacent rungs after every block
(even/odd rungs alternating). greedy = the same chains all at T = 0 (accept only dNLL < 0).
Chains start from independent uniformly random wirings.
"""
import time
import numpy as np
import jax
import jax.numpy as jnp
import model
import delta
import optimize as opt

def run(P, method='pt', chains=16, steps=20000, block=500, tmin=1.0, tmax=1000.0, seed=0, D=None,
        nll_ref=0.0, trace=None, log=print, log_every=10, topk=1, history=None, history_every=0):
    """P: model.load_problem(...). Returns dict(b = best wiring (delta encoding), nll, perms (hard permutations
    {node: (CN,CN)}), summary, top). nll_ref is subtracted from every reported NLL (e.g. the truth NLL -> gaps).
    trace: open text file for the per-block curve (steps, best, current NLL per rung, all minus nll_ref).
    history, history_every > 0: open text file for the NLL of every chain every history_every proposals (and at
    step 0), with the temperature the chain had then: columns steps, chain, T, nll (minus nll_ref).
    topk > 1: each chain also keeps the 2*topk lowest-NLL states it visited (any rung), counting states whose NLLs
    agree to 1e-9 (relative) once (copy relabelings of one karyotype have the same NLL); top = the topk best of
    these over all chains as a list of dict(nll, perms), distinct as karyotypes (dedup in topk_karyotypes())."""
    if method not in ('pt', 'greedy'):
        raise ValueError(f'method must be pt or greedy, not {method}')
    C = chains
    temps = np.zeros(C) if method == 'greedy' else np.geomspace(tmin, tmax, C)
    nb = max(steps // block, 1)
    if len(P.var_nodes) == 0:  # nothing to search: every node has CN 1, the wiring is determined
        nll = float(model.make_model(P).loss_perms({}))
        summ = dict(method=method, seed=seed, chains=C, steps=0, block=block, temps=list(temps), best_nll=nll,
                    best_rel=nll - nll_ref, n_var_nodes=0, n_moves=0, wall_s=0.0)
        return dict(b=np.zeros(0, int), nll=nll, perms={}, summary=summ, top=[dict(nll=nll, perms={})])
    D = D or delta.make_delta(P)

    moves = jnp.asarray(D['moves']); mnode = D['move_node']
    n_nodes = len(D['nodes'])
    off = jnp.asarray(np.searchsorted(mnode, np.arange(n_nodes))); nmv = jnp.asarray(np.bincount(mnode, minlength=n_nodes))

    K = 2 * topk if topk > 1 else 0  # per-chain buffer of distinct low-NLL states

    def buf_insert(bn, bbuf, nll, b):
        dup = jnp.any(jnp.abs(bn - nll) <= 1e-9 * jnp.maximum(1.0, jnp.abs(nll)))
        w = jnp.argmax(bn)
        ins = (nll < bn[w]) & ~dup
        return bn.at[w].set(jnp.where(ins, nll, bn[w])), bbuf.at[w].set(jnp.where(ins, b, bbuf[w]))

    @jax.jit
    def run_block(state, best, key, T, buf):
        def step(c, k):
            st, bnll, bb, acc, bf = c
            k1, k2, k3 = jax.random.split(k, 3)
            nd = jax.random.randint(k1, (), 0, n_nodes)
            i, j = moves[off[nd] + jax.random.randint(k2, (), 0, nmv[nd])]
            dn, st2 = D['propose'](st, i, j)
            ok = (dn < 0) | (jnp.log(jax.random.uniform(k3)) < -dn / jnp.maximum(T, 1e-300))
            st = jax.tree.map(lambda u, v: jnp.where(ok, u, v), st2, st)
            better = st['nll'] < bnll
            if K:
                bf = buf_insert(bf[0], bf[1], st['nll'], st['b'])
            return (st, jnp.where(better, st['nll'], bnll), jnp.where(better, st['b'], bb), acc + ok, bf), st['nll']
        (state, bnll, bb, acc, buf), ys = jax.lax.scan(step, (state, best[0], best[1], 0, buf),
                                                        jax.random.split(key, block))
        return state, (bnll, bb), acc, buf, ys

    key = jax.random.PRNGKey(seed)
    r0 = opt.random_perms(P, jax.random.fold_in(key, 999), C)
    states = [D['init'](D['b_from_perms']({n: v[c] for n, v in r0.items()})) for c in range(C)]
    bests = [(s['nll'], s['b']) for s in states]
    bufs = [buf_insert(jnp.full(K, jnp.inf), jnp.zeros((K,) + s['b'].shape, s['b'].dtype), s['nll'], s['b'])
            if K else () for s in states]
    rung = np.arange(C)  # rung[c] = temperature index of chain c
    hist = history is not None and history_every > 0
    if hist:
        history.write('steps\tchain\tT\tnll\n')
        for c in range(C):
            history.write(f"0\t{c}\t{temps[rung[c]]:g}\t{float(states[c]['nll']) - nll_ref:.10g}\n")
    if trace is not None:
        trace.write('steps\tbest\t' + '\t'.join(f'T{t:g}' for t in temps) + '\n')
    acc_tot = np.zeros(C); sw_try = np.zeros(max(C - 1, 1)); sw_acc = np.zeros(max(C - 1, 1))
    t0 = time.time()
    for blk in range(nb):
        for c in range(C):
            key, k = jax.random.split(key)
            states[c], bests[c], acc, bufs[c], ys = run_block(states[c], bests[c], k, temps[rung[c]], bufs[c])
            acc_tot[rung[c]] += float(acc)
            if hist:
                ix = np.arange(history_every - 1, block, history_every)
                yv = np.asarray(ys)[ix] - nll_ref
                history.write(''.join(f"{blk * block + i + 1}\t{c}\t{temps[rung[c]]:g}\t{v:.10g}\n"
                                      for i, v in zip(ix, yv)))
        E = np.array([float(s['nll']) for s in states])
        if method == 'pt':
            at = np.argsort(rung)  # chain at each rung
            for r in range(blk % 2, C - 1, 2):
                ca, cb = at[r], at[r + 1]
                d = (1 / temps[r] - 1 / temps[r + 1]) * (E[ca] - E[cb])
                sw_try[r] += 1
                if np.log(np.random.default_rng([seed, blk, r]).random()) < d:
                    rung[ca], rung[cb] = r + 1, r; sw_acc[r] += 1
        bg = min(float(b[0]) for b in bests) - nll_ref
        at = np.argsort(rung)
        if trace is not None:
            trace.write(f'{(blk + 1) * block}\t{bg:.6g}\t' + '\t'.join(f'{E[at[r]] - nll_ref:.6g}' for r in range(C)) + '\n')
            trace.flush()
        if log and (blk % log_every == 0 or blk == nb - 1):
            log(f'block {blk + 1}/{nb} best {bg:.1f} rungs (min/mid/max T) {np.round(E[at] - nll_ref)[[0, C // 2, C - 1]]} '
                f'{time.time() - t0:.0f}s')
    wall = time.time() - t0
    ib = int(np.argmin([float(b[0]) for b in bests]))
    b = np.asarray(bests[ib][1]); nll = float(bests[ib][0])
    perms = {n: np.asarray(M) for n, M in D['perms_from_b'](b).items()}
    top = [dict(nll=nll, perms=perms)]
    if K:
        cand = sorted(((float(v), np.asarray(bb)) for bf in bufs for v, bb in zip(bf[0], bf[1]) if np.isfinite(v)),
                      key=lambda t: t[0])
        top = [dict(nll=v, perms={n: np.asarray(M) for n, M in D['perms_from_b'](bb).items()}) for v, bb in cand]
        top = topk_karyotypes(P, [dict(nll=nll, perms=perms)] + top, topk)
    summ = dict(method=method, seed=seed, chains=C, steps=nb * block, block=block, temps=list(temps), best_nll=nll,
                best_rel=nll - nll_ref, n_var_nodes=n_nodes, n_moves=len(D['moves']),
                acc_rate_per_rung=list(acc_tot / (nb * block)), swap_rate=list(sw_acc / np.maximum(sw_try, 1)),
                wall_s=wall, ms_per_proposal=wall / (nb * block * C) * 1e3)
    return dict(b=b, nll=nll, perms=perms, summary=summ, top=top)

def topk_karyotypes(P, cands, k):
    """the k lowest-NLL entries of cands (dicts with nll, perms) that are distinct karyotypes; adds walks/circular."""
    import walks
    out, seen = [], set()
    for c in sorted(cands, key=lambda c: c['nll']):
        wk, circ = walks.perms_to_walks(P, c['perms'])
        key = tuple(sorted(walks.canonical(w, x) for w, x in zip(wk, circ)))
        if key in seen:
            continue
        seen.add(key)
        out.append(dict(c, walks=wk, circular=circ))
        if len(out) == k:
            break
    return out
