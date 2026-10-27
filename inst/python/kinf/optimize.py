"""Hard-wiring search utilities on the exact NLL (model.py): random wirings, greedy transposition search ("polish"),
random baselines. The Sinkhorn optimiser that used to live here is archived in archive/sinkhorn/optimize.py.
"""
import time
import numpy as np
import jax
import jax.numpy as jnp

def make_hard_nll(m):
    return jax.jit(lambda perms: m.nll(m.hard_expected_map(perms)))

def polish(perms, hard_nll, max_sweeps=50):
    """greedy hard local search: at each node try every transposition of the pairing, keep improvements,
    sweep until no single transposition improves the NLL."""
    perms = {n: np.array(M) for n, M in perms.items()}
    best = float(hard_nll(perms)); n_acc = 0
    for sweep in range(max_sweeps):
        improved = False
        for n in sorted(perms):
            k = perms[n].shape[0]
            for a in range(k):
                for b in range(a + 1, k):
                    trial = dict(perms)
                    M = perms[n].copy(); M[[a, b]] = M[[b, a]]
                    trial[n] = M
                    v = float(hard_nll(trial))
                    if v < best - 1e-9:
                        perms, best, improved = trial, v, True
                        n_acc += 1
        if not improved:
            break
    return perms, best, dict(sweeps=sweep + 1, accepted=n_acc)

def random_perms(problem, key, k):
    """k random hard wirings, stacked: {node: (k, CN, CN)}"""
    out = {}
    keys = jax.random.split(key, len(problem.var_nodes))
    for kk, n in zip(keys, sorted(problem.var_nodes)):
        cn = len(problem.var_nodes[n][0])
        idx = jax.vmap(lambda q: jax.random.permutation(q, cn))(jax.random.split(kk, k))
        out[n] = jax.nn.one_hot(idx, cn, dtype=jnp.float64)
    return out

def random_baseline(problem, m, key, k=1000, batch=250):
    """best hard NLL among k uniformly random wirings (same search budget idea as fastKar's sampling)"""
    f = jax.jit(jax.vmap(lambda p: m.nll(m.hard_expected_map(p))))
    best, best_perms = np.inf, None
    for s in range(0, k, batch):
        P = random_perms(problem, jax.random.fold_in(key, s), min(batch, k - s))
        v = np.asarray(f(P))
        i = int(np.argmin(v))
        if v[i] < best:
            best, best_perms = float(v[i]), {n: np.asarray(M[i]) for n, M in P.items()}
    return best, best_perms

def random_climb(problem, m, key, n_starts=8, hard_nll=None):
    """baseline without gradients: greedy transposition search (polish) from n_starts random wirings"""
    hard_nll = hard_nll or make_hard_nll(m)
    P = random_perms(problem, key, n_starts)
    t0 = time.time()
    res = [polish({n: np.asarray(M[s]) for n, M in P.items()}, hard_nll) for s in range(n_starts)]
    v = np.array([r[1] for r in res])
    b = int(np.argmin(v))
    return dict(perms=res[b][0], nll=v[b], nlls=v, accepted=[r[2]['accepted'] for r in res], t=time.time() - t0)
