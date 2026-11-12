"""problem.py reproduces regression.json (frozen from the old modules); walkdist == Floyd-Warshall; Graph.from_tables == Graph.load."""
import os, sys, json, shutil, tempfile
import numpy as np, pandas as pd
import jax, jax.numpy as jnp
jax.config.update('jax_enable_x64', True)
here = os.path.dirname(os.path.abspath(__file__)); kin = os.path.dirname(here)
sys.path.insert(0, kin)
import problem as pr

folder = os.path.join(kin, 'selftest')
R = json.load(open(os.path.join(here, 'regression.json')))
res = []
def check(name, ok):
    res.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name)
close = lambda a, b, tol=1e-9: np.allclose(a, b, rtol=tol, atol=0)

g = pr.Graph.load(folder); h = pr.HiC.load(folder, g)
sn, circ = pr.read_gwalks(folder + '/truth_gwalks.json')
b_of = lambda perms: pr.b_from_perms(g, perms)
check('truth NLL', close(pr.nll(g, h, b_of(pr.walks_to_perms(g, sn, circ))), R['truth_nll']))
bs = [jnp.asarray(b) for b in R['b']]
check('random wirings valid (check_b)', all(len(pr.check_b(g, b)) == len(g.var_a) for b in bs))
rb = pr.random_b(g, jax.random.PRNGKey(0), 10)
check('random_b == frozen wirings', (np.asarray(rb) == np.array(R['b'])).all())
check('10 full NLLs', close([float(pr.nll(g, h, b)) for b in bs], R['nll']))
check('perms_from_b round trip', all((np.asarray(pr.b_from_perms(g, pr.perms_from_b(g, b))) == np.asarray(b)).all() for b in bs))
T = pr.make_nll_term(g, h)
st = T['init'](bs[0]); check('init NLL', close(float(st['nll']), R['init_nll']))
seq = []
for k in R['moves']:
    i, j = g.moves[k]; _, st = T['propose'](st, int(i), int(j)); seq.append(float(st['nll']))
check('300 incremental NLLs', close(seq, R['move_nll']))
check('incremental == full at the end', close(float(st['nll']), float(pr.nll(g, h, st['b']))))
kars = []
for b in bs:
    w, c = pr.perms_to_walks(g, pr.perms_from_b(g, b))
    kars.append(sorted([list(map(int, k[1])), bool(k[0])] for k in (pr.canonical(a, cc) for a, cc in zip(w, c))))
check('canonical karyotypes', kars == R['canonical'])
w0, c0 = pr.perms_to_walks(g, pr.walks_to_perms(g, sn, circ))
check('truth walks round trip same_karyotype', pr.same_karyotype(w0, c0, sn, circ))
h2 = pr.HiC.load(folder, g)
check('rescale_depth ratio', close(pr.rescale_depth(g, h2, 64, 0), R['rescale_ratio']))

# (b) walk traversal vs Floyd-Warshall
ok = True
for b in bs:
    Dw = pr.walkdist(g, b); Df = pr.floyd_warshall(pr.fw_costmat(g, pr.perms_from_b(g, b)))
    ok &= bool(jnp.all(jnp.abs(Dw - Df) <= 1e-9 * jnp.maximum(1.0, jnp.abs(Df))))
check('walkdist == Floyd-Warshall (10 wirings)', ok)
check('expected_map == expected_map_fw', close(pr.expected_map(g, h, bs[0]), pr.expected_map_fw(g, h, pr.perms_from_b(g, bs[0]))))

# (c) from_tables vs load
nd = pd.read_csv(folder + '/nodes.csv', index_col=0).sort_values('unique.id')
nodes = pd.DataFrame([dict(node_id=n, cn=int((q['side'] == 'left').sum()), width=q[q['side'] == 'left']['width'].iloc[0] + q[q['side'] == 'right']['width'].iloc[0]) for n, q in nd.groupby('node_id', sort=True)])
J = [dict(n1=int(g.node_of_uid[a]), side1='right' if g.side_of_uid[a] else 'left',
          n2=int(g.node_of_uid[c]), side2='right' if g.side_of_uid[c] else 'left', cn=1) for a, c in g.ext]
g2 = pr.Graph.from_tables(nodes, pd.DataFrame(J))
same = all(np.array_equal(getattr(g, k), getattr(g2, k)) for k in
           ['nu', 'node_of_uid', 'side_of_uid', 'width_uid', 'ext', 'inner_fixed', 'cn', 'start', 'var_a', 'lp0', 'moves', 'move_node'])
same &= g.var_nodes.keys() == g2.var_nodes.keys() and all(np.array_equal(g.var_nodes[n][i], g2.var_nodes[n][i]) for n in g.var_nodes for i in (0, 1))
check('from_tables == load (Graph arrays)', same)

# (d) graph-only folder
with tempfile.TemporaryDirectory() as d:
    for f in ['nodes', 'internal_edges', 'external_edges']:
        shutil.copy(f'{folder}/{f}.csv', d)
    g3 = pr.Graph.load(d)
    check('Graph.load on 3-CSV folder', np.array_equal(g3.ext, g.ext) and np.array_equal(g3.var_a, g.var_a))
print(f'{sum(res)} PASS, {len(res) - sum(res)} FAIL')
