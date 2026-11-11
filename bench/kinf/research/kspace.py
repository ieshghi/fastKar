"""Estimate the number of distinct karyotypes K(G) of exported problems from uniformly random wirings.
Draws n random wirings (one uniform permutation per variable node), canonicalises the karyotype (walk
rotation/reversal, copy relabeling), and reports distinct count, singletons and Chao1.
Usage: python bench/kinf/research/kspace.py [-n 20000] folder [folder ...]   (prints one tab-separated line per folder)
"""
import sys, os, math, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../inst/python/kinf'))
from collections import Counter
import numpy as np
import model, walks

def kspace(folder, n=20000, seed=0):
    P = model.load_problem(folder)
    rng = np.random.default_rng(seed)
    cnt = Counter()
    for _ in range(n):
        perms = {k: np.eye(len(l))[rng.permutation(len(l))] for k, (l, r) in P.var_nodes.items()}
        w, c = walks.perms_to_walks(P, perms)
        cnt[tuple(sorted(walks.canonical(a, b) for a, b in zip(w, c)))] += 1
    f1 = sum(v == 1 for v in cnt.values()); f2 = sum(v == 2 for v in cnt.values())
    cns = [len(l) for l, r in P.var_nodes.values()]
    return dict(folder=os.path.basename(folder.rstrip('/')), n_var_nodes=len(cns), max_cn=max(cns, default=1),
                n_halfcopies=P.nu, n_tilecopies=len(P.tc_uid), n_tiles=P.ntile,
                log10_wirings=sum(math.log10(math.factorial(k)) for k in cns), draws=n, distinct=len(cnt),
                singletons=f1, chao1=len(cnt) + f1 * (f1 - 1) / (2 * (f2 + 1)))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-n', type=int, default=20000)
    ap.add_argument('folders', nargs='+')
    a = ap.parse_args()
    for i, f in enumerate(a.folders):
        r = kspace(f, a.n)
        if i == 0:
            print('\t'.join(r))
        print('\t'.join(str(v) for v in r.values()), flush=True)
