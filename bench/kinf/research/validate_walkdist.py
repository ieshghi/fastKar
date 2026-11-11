"""Walk-traversal hard map (model.hard_expected_map) vs Floyd-Warshall hard map (hard_expected_map_fw).
Usage: python bench/kinf/research/validate_walkdist.py ROOT [folders...]   (default: every problem folder under ROOT)
Checks the truth and 20 random wirings per problem; prints max relative map difference, NLL difference, timing."""
import sys, os, glob, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../inst/python/kinf'))
import jax; jax.config.update('jax_enable_x64', True)
import jax.numpy as jnp
import numpy as np
import model, walks, optimize as opt

root = sys.argv[1].rstrip('/') + '/'
folders = sys.argv[2:] or sorted(os.path.relpath(os.path.dirname(f), root) for f in glob.glob(root + '**/truth_gwalks.json', recursive=True))
worst = 0.0
for f in folders:
    P = model.load_problem(root + f)
    m = model.make_model(P)
    walk = jax.jit(jax.vmap(m.hard_expected_map)); fw = jax.jit(jax.vmap(m.hard_expected_map_fw))
    sn, circ = walks.read_gwalks(root + f + '/truth_gwalks.json')
    T = walks.walks_to_perms(P, sn, circ)
    R = opt.random_perms(P, jax.random.PRNGKey(0), 20)
    B = {n: jnp.concatenate([jnp.asarray(T[n])[None], R[n]]) for n in R}
    a = np.asarray(walk(B)); b = np.asarray(fw(B))
    rel = float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))
    dn = float(np.max(np.abs(np.asarray(jax.vmap(m.nll)(a)) - np.asarray(jax.vmap(m.nll)(b)))))
    t = {}
    for name, g in (('walk', walk), ('fw', fw)):
        g(B).block_until_ready(); t0 = time.time(); g(B).block_until_ready(); t[name] = (time.time() - t0) / 21 * 1e3
    worst = max(worst, rel)
    print(f'{f}: N_ends={2 * P.nu} max_rel={rel:.2e} max_dNLL={dn:.2e} ms/eval walk={t["walk"]:.2f} fw={t["fw"]:.2f}', flush=True)
print('WORST max_rel', worst)
