"""Summarise wg_search.py runs: re-score each best wiring with the full model (checks delta drift), compare it to
the truth (exact karyotype; recall of the multiset of walk adjacencies, i.e. consecutive signed-node pairs), plot
best gap vs proposals. Usage: python bench/kinf/research/wg_analyze.py RESULTS_DIR PLOT_PDF"""
import sys, os, glob, json
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../inst/python/kinf'))
import jax; jax.config.update('jax_enable_x64', True)
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import model, delta, walks, optimize as opt

def adjacencies(ws, circ):
    c = Counter()
    for w, cy in zip(ws, circ):
        w = list(w); pairs = list(zip(w, w[1:])) + ([(w[-1], w[0])] if cy and len(w) > 1 else [])
        c.update(min((x, y), (-y, -x)) for x, y in pairs)
    return c

res, cache = [], {}
for j in sorted(glob.glob(sys.argv[1].rstrip('/') + '/*.json')):
    s = json.load(open(j)); f = s['folder']
    if f not in cache:
        P = model.load_problem(f); m = model.make_model(P)
        sn, circ = walks.read_gwalks(f + '/truth_gwalks.json')
        cache[f] = (P, m, opt.make_hard_nll(m), delta.make_delta(P), sn, circ)
    P, m, hn, D, sn, circ = cache[f]
    perms = D['perms_from_b'](np.load(j.replace('.json', '_bestb.npy')))
    gap_full = float(hn(perms)) - float(hn(walks.walks_to_perms(P, sn, circ)))
    w2, c2 = walks.perms_to_walks(P, perms)
    A, B = adjacencies(w2, c2), adjacencies(sn, circ)
    res.append(dict(truth=os.path.basename(f), method=s['method'], seed=s['seed'], gap=s['best_gap'], gap_full=gap_full,
                    exact=walks.same_karyotype(w2, c2, sn, circ), adj_recall=sum((A & B).values()) / sum(B.values()),
                    n_walks=len(w2), n_walks_truth=len(sn), n_cyc=int(sum(c2)), n_cyc_truth=int(sum(circ)),
                    swap_rate_min=min(s['swap_rate']) if s['method'] == 'pt' else np.nan))
d = pd.DataFrame(res)
pd.set_option('display.width', 200)
print(d.round(3).to_string(index=False))
d.to_csv(sys.argv[1].rstrip('/') + '/summary.tsv', sep='\t', index=False)

truths = sorted(d.truth.unique())
fig, ax = plt.subplots(1, len(truths), figsize=(4.2 * len(truths), 3.6), sharey=True)
col = {'greedy': '#c0504d', 'pt': '#2f6db5'}
for a_, t in zip(np.atleast_1d(ax), truths):
    for tsv in sorted(glob.glob(sys.argv[1].rstrip('/') + f'/{t}_*.tsv')):
        meth = 'pt' if '_pt_' in tsv else 'greedy'; seed = tsv.split('_s')[-1].split('.')[0]
        c = pd.read_csv(tsv, sep='\t')
        a_.plot(c.steps, c['best_gap' if 'best_gap' in c else 'best'].clip(lower=1), color=col[meth], ls='-' if seed == '0' else '--', label=f'{meth} s{seed}')
    a_.set_yscale('log'); a_.set_xscale('log'); a_.set_title(f'SK-OV-3 whole genome, {t}'); a_.set_xlabel('steps per chain (16 chains)')
np.atleast_1d(ax)[0].set_ylabel('best NLL - truth NLL'); np.atleast_1d(ax)[0].legend(frameon=False, fontsize=8)
plt.tight_layout(); plt.savefig(sys.argv[2]); print('plot', sys.argv[2])
