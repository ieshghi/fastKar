"""Freeze the numbers of the OLD modules (model, delta, walks, optimize, infer) on selftest/ -> regression.json."""
import os, sys, json, copy
import numpy as np
import jax
jax.config.update('jax_enable_x64', True)
here = os.path.dirname(os.path.abspath(__file__)); kin = os.path.dirname(here)
sys.path.insert(0, kin)
import model, delta, walks, optimize as opt, infer

folder = os.path.join(kin, 'selftest')
P = model.load_problem(folder); m = model.make_model(P); D = delta.make_delta(P)
sn, circ = walks.read_gwalks(folder + '/truth_gwalks.json')
out = dict(truth_nll=float(m.loss_perms(walks.walks_to_perms(P, sn, circ))))
rp = opt.random_perms(P, jax.random.PRNGKey(0), 10)
bs, nlls, kars = [], [], []
for c in range(10):
    perms = {n: v[c] for n, v in rp.items()}
    bs.append([int(x) for x in D['b_from_perms'](perms)])
    nlls.append(float(m.loss_perms(perms)))
    w, cc = walks.perms_to_walks(P, perms)
    kars.append(sorted([list(map(int, k[1])), bool(k[0])] for k in (walks.canonical(a, b) for a, b in zip(w, cc))))
out.update(b=bs, nll=nlls, canonical=kars)
rng = np.random.RandomState(1)
mv = [int(x) for x in rng.randint(0, len(D['moves']), 300)]
st = D['init'](jax.numpy.asarray(bs[0])); seq = []
for k in mv:
    i, j = D['moves'][k]
    _, st = D['propose'](st, int(i), int(j)); seq.append(float(st['nll']))
out.update(moves=mv, move_nll=seq, init_nll=float(D['init'](jax.numpy.asarray(bs[0]))['nll']))
P2 = model.load_problem(folder)
out['rescale_ratio'] = infer.rescale_depth(P2, 64, 0)
json.dump(out, open(os.path.join(here, 'regression.json'), 'w'))
print('truth', out['truth_nll'], 'ratio', out['rescale_ratio'], 'nmoves', len(mv), 'nll0', nlls[0])
