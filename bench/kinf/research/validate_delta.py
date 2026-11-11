"""Delta NLL (delta.py) vs full recompute. Usage: python bench/kinf/research/validate_delta.py FOLDER [n_moves=200]
From a random wiring, proposes n random transpositions, accepts half of them, and checks every dnll against
the full hard NLL difference, then the carried map/NLL against a fresh evaluation. Prints timing and |S|."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../inst/python/kinf'))
import jax; jax.config.update('jax_enable_x64', True)
import jax.numpy as jnp
import numpy as np
import model, delta, optimize as opt

folder = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
P = model.load_problem(folder); m = model.make_model(P); hn = opt.make_hard_nll(m)
D = delta.make_delta(P)
r0 = opt.random_perms(P, jax.random.PRNGKey(0), 1)
st = D['init'](D['b_from_perms']({k: v[0] for k, v in r0.items()}))
print('init nll delta-model', float(st['nll']), 'full', float(hn({k: v[0] for k, v in r0.items()})), flush=True)
prop = jax.jit(D['propose']); aff = jax.jit(D['affected'])
rng = np.random.default_rng(1); mv = D['moves']
errs, sizes, ts = [], [], []
for t in range(n):
    i, j = mv[rng.integers(len(mv))]
    sizes.append(int(aff(st, i, j)))
    t0 = time.time(); dn, st2 = prop(st, i, j); dn = float(dn); ts.append(time.time() - t0)
    full_new = float(hn(D['perms_from_b'](st2['b'])))
    errs.append(abs((float(st['nll']) + dn) - full_new))
    if rng.random() < 0.5:
        st = st2
fresh = D['init'](st['b'])
print(f'moves {n}: max |nll+dnll - full| = {max(errs):.3e}  (nll scale {float(st["nll"]):.4g})')
print(f'carried vs fresh: nll diff {abs(float(st["nll"]) - float(fresh["nll"])):.3e}, map max rel '
      f'{float(jnp.max(jnp.abs(st["mu"] - fresh["mu"]) / jnp.maximum(jnp.abs(fresh["mu"]), 1e-300))):.3e}')
ts = np.array(ts[20:]) * 1e3  # skip compile
print(f'|S| median {np.median(sizes):.0f}, 90% {np.percentile(sizes, 90):.0f}, max {max(sizes)}; ms/proposal median {np.median(ts):.2f} mean {ts.mean():.2f}')
