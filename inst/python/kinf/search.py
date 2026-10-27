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
        nll_ref=0.0, trace=None, log=print, log_every=10):
    """P: model.load_problem(...). Returns dict(b = best wiring (delta encoding), nll, perms (hard permutations
    {node: (CN,CN)}), summary). nll_ref is subtracted from every reported NLL (e.g. the truth NLL -> gaps).
    trace: open text file for the per-block curve (steps, best, current NLL per rung, all minus nll_ref)."""
    if method not in ('pt', 'greedy'):
        raise ValueError(f'method must be pt or greedy, not {method}')
    C = chains
    temps = np.zeros(C) if method == 'greedy' else np.geomspace(tmin, tmax, C)
    nb = max(steps // block, 1)
    if len(P.var_nodes) == 0:  # nothing to search: every node has CN 1, the wiring is determined
        nll = float(model.make_model(P).loss_perms({}))
        summ = dict(method=method, seed=seed, chains=C, steps=0, block=block, temps=list(temps), best_nll=nll,
                    best_rel=nll - nll_ref, n_var_nodes=0, n_moves=0, wall_s=0.0)
        return dict(b=np.zeros(0, int), nll=nll, perms={}, summary=summ)
    D = D or delta.make_delta(P)

    moves = jnp.asarray(D['moves']); mnode = D['move_node']
    n_nodes = len(D['nodes'])
    off = jnp.asarray(np.searchsorted(mnode, np.arange(n_nodes))); nmv = jnp.asarray(np.bincount(mnode, minlength=n_nodes))

    @jax.jit
    def run_block(state, best, key, T):
        def step(c, k):
            st, bnll, bb, acc = c
            k1, k2, k3 = jax.random.split(k, 3)
            nd = jax.random.randint(k1, (), 0, n_nodes)
            i, j = moves[off[nd] + jax.random.randint(k2, (), 0, nmv[nd])]
            dn, st2 = D['propose'](st, i, j)
            ok = (dn < 0) | (jnp.log(jax.random.uniform(k3)) < -dn / jnp.maximum(T, 1e-300))
            st = jax.tree.map(lambda u, v: jnp.where(ok, u, v), st2, st)
            better = st['nll'] < bnll
            return (st, jnp.where(better, st['nll'], bnll), jnp.where(better, st['b'], bb), acc + ok), None
        (state, bnll, bb, acc), _ = jax.lax.scan(step, (state, best[0], best[1], 0), jax.random.split(key, block))
        return state, (bnll, bb), acc

    key = jax.random.PRNGKey(seed)
    r0 = opt.random_perms(P, jax.random.fold_in(key, 999), C)
    states = [D['init'](D['b_from_perms']({n: v[c] for n, v in r0.items()})) for c in range(C)]
    bests = [(s['nll'], s['b']) for s in states]
    rung = np.arange(C)  # rung[c] = temperature index of chain c
    if trace is not None:
        trace.write('steps\tbest\t' + '\t'.join(f'T{t:g}' for t in temps) + '\n')
    acc_tot = np.zeros(C); sw_try = np.zeros(max(C - 1, 1)); sw_acc = np.zeros(max(C - 1, 1))
    t0 = time.time()
    for blk in range(nb):
        for c in range(C):
            key, k = jax.random.split(key)
            states[c], bests[c], acc = run_block(states[c], bests[c], k, temps[rung[c]])
            acc_tot[rung[c]] += float(acc)
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
    summ = dict(method=method, seed=seed, chains=C, steps=nb * block, block=block, temps=list(temps), best_nll=nll,
                best_rel=nll - nll_ref, n_var_nodes=n_nodes, n_moves=len(D['moves']),
                acc_rate_per_rung=list(acc_tot / (nb * block)), swap_rate=list(sw_acc / np.maximum(sw_try, 1)),
                wall_s=wall, ms_per_proposal=wall / (nb * block * C) * 1e3)
    return dict(b=b, nll=nll, perms=perms, summary=summ)
