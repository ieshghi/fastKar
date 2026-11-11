"""Benchmark driver for search.py on simulated problems (needs FOLDER/truth_gwalks.json): greedy (T = 0) vs PT.
For general use (no truth, gwalks.json output) see infer.py.
Usage: python bench/kinf/research/wg_search.py FOLDER TAG --method pt|greedy [--chains 16 --steps 20000 --block 500
       --tmin 1 --tmax 1000 --seed 0]
Writes <folder>/../results/<TAG>/<problem>_<method>_s<seed>.{tsv,json}: curve (proposals per chain, best gap,
per-rung current gap) and final summary (best gap to the truth NLL, acceptance, timing), plus _bestb.npy.
"""
import sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../inst/python/kinf'))
import jax; jax.config.update('jax_enable_x64', True)
import numpy as np
import model, delta, walks, search

ap = argparse.ArgumentParser()
ap.add_argument('folder'); ap.add_argument('tag')
ap.add_argument('--method', default='pt'); ap.add_argument('--chains', type=int, default=16)
ap.add_argument('--steps', type=int, default=20000); ap.add_argument('--block', type=int, default=500)
ap.add_argument('--tmin', type=float, default=1.0); ap.add_argument('--tmax', type=float, default=1000.0)
ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args()
folder = a.folder.rstrip('/')
P = model.load_problem(folder); D = delta.make_delta(P)
sn, circ = walks.read_gwalks(folder + '/truth_gwalks.json')
nll_truth = float(D['init'](D['b_from_perms'](walks.walks_to_perms(P, sn, circ)))['nll'])

out = os.path.join(os.path.dirname(folder), 'results', a.tag); os.makedirs(out, exist_ok=True)
base = os.path.join(out, f'{os.path.basename(folder)}_{a.method}_s{a.seed}')

with open(base + '.tsv', 'w') as f:  # columns steps, best, T<t> (gaps to the truth NLL)
    res = search.run(P, method=a.method, chains=a.chains, steps=a.steps, block=a.block, tmin=a.tmin, tmax=a.tmax,
                     seed=a.seed, D=D, nll_ref=nll_truth, trace=f, log=lambda s: print(s, flush=True))
s = res['summary']
summ = dict(folder=folder, method=a.method, seed=a.seed, chains=a.chains, steps=a.steps, temps=s['temps'],
            nll_truth=nll_truth, best_gap=s['best_rel'], acc_rate_per_rung=s['acc_rate_per_rung'],
            swap_rate=s['swap_rate'], wall_s=s['wall_s'], ms_per_proposal=s['ms_per_proposal'])
json.dump(summ, open(base + '.json', 'w'), indent=1)
np.save(base + '_bestb.npy', res['b'])
print(json.dumps({k: v for k, v in summ.items() if k not in ('temps',)}, default=float))
