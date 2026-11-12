"""optimize.py on the selftest problem (agreement with the old engine was verified before it was removed, commit 3798ff7); incremental NLL; scores; plain-score runs."""
import os, sys
import numpy as np
import jax, jax.numpy as jnp
jax.config.update('jax_enable_x64', True)
here = os.path.dirname(os.path.abspath(__file__)); kin = os.path.dirname(here)
sys.path.insert(0, kin)
import problem as pr, optimize as op, scores as sc

folder = os.path.join(kin, 'selftest')
g = pr.Graph.load(folder); h = pr.HiC.load(folder, g)
def check(name, ok): print(('PASS ' if ok else 'FAIL ') + name, flush=True)
quiet = lambda s: None

nl = pr.nll_term(g, h)

# 2. incremental vs full NLL after 300 random moves
b = pr.random_b(g, jax.random.PRNGKey(5), 1)[0]
aux = nl.init(op.kar(g, b)); ko = op.kar(g, b); rng = np.random.default_rng(0)
upd = jax.jit(lambda aux, ko, kn, i, j: nl.update(aux, ko, kn, i, j))
for _ in range(300):
    i, j = g.moves[rng.integers(len(g.moves))]
    b2 = b.at[i].set(b[j]).at[j].set(b[i]); kn = op.kar(g, b2)
    aux = upd(aux, ko, kn, i, j); b, ko = b2, kn
full = float(pr.nll(g, h, b)); inc = float(nl.value(aux))
check(f'incremental NLL after 300 moves (rel {abs(full-inc)/abs(full):.1e})', abs(full - inc) <= 1e-8 * abs(full))

# 3. scores
bs = pr.random_b(g, jax.random.PRNGKey(7), 10)
amp = [int(n) for n in g.vnodes[:3]]
fw, fc, fe = jax.jit(sc.n_walks(g)), jax.jit(sc.n_circular(g)), jax.jit(sc.amp_entropy(g, amp))
ok1 = ok2 = ok3 = True
for b in bs:
    wk, circ = pr.perms_to_walks(g, pr.perms_from_b(g, b)); k = op.kar(g, b)
    ok1 &= int(fw(k)) == len(wk); ok2 &= int(fc(k)) == int(sum(circ))
    c = np.array([sum(abs(int(x)) in amp for x in w) for w in wk], float)
    p = c / c.sum(); p = p[p > 0]
    ok3 &= abs(float(fe(k)) + np.sum(p * np.log(p))) < 1e-10
check('n_walks', ok1); check('n_circular', ok2); check('amp_entropy', ok3)

# 4. plain scores
for vm in (False, True):
    r = op.optimize(g, [sc.n_circular(g)], T=0.5, chains=4, steps=100, block=50, seed=2, vmap=vm, log=quiet)
    check(f'n_circular run vmap={vm} (best {r["best"]["E"]})', np.isfinite(r['E']).all() and np.isfinite(r['best']['E']))
r = op.optimize(g, [sc.n_circular(g)], energy=sc.target(0, sc.schedule([0, 200], [0, 3])), T=0.3, chains=4, steps=200, block=50,
                seed=3, history_every=25, log=quiet)
check('target/schedule run', np.isfinite(r['E']).all() and r['history']['E'].shape[0] == 4 * (1 + 8) and r['t_end'] == 200)
r2 = op.optimize(g, [nl], T=lambda t: jnp.geomspace(1, 100, 4) * (1 + 0 * t), chains=4, steps=100, block=50, seed=4, log=quiet)
r3 = op.optimize(g, [nl], T=np.geomspace(1, 100, 4), chains=4, steps=100, block=50, seed=4, init=r2, log=quiet)
check('callable T and continuation', np.isfinite(r2['E']).all() and r3['t_end'] == 200)
