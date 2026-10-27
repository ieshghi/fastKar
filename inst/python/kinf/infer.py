"""Command-line entry point: infer a karyotype (wiring of a genome graph) from Hi-C by search on the exact NLL.

Input: a problem folder written by R/prep_data_for_py.R (or from R, use infer_karyotype() in R/infer_karyotype.R,
which writes the folder, calls this script and reads the result back as a gWalk).
Output (OUT): gwalks.json (best karyotype: signed snode ids on the gGraph passed to prep_data_for_py, circular
flags; read in R with read_kinf_gwalks()), summary.json, trace.tsv (best and per-rung NLL per block), best_b.npy.

Usage: python infer.py PROBLEM OUT [--method pt|greedy] [--chains 16] [--steps 20000] [--block 500]
       [--tmin 1] [--tmax 1000] [--seed 0] [--device auto|gpu|cpu] [--truth]
Run with PYTHONPATH unset (env -u PYTHONPATH); entries from PYTHONPATH are dropped from sys.path here as well.
"""
import os, sys
_pp = [os.path.abspath(p) for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p]
sys.path[:] = [p for p in sys.path if os.path.abspath(p or '.') not in _pp]  # foreign PYTHONPATH breaks numpy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import argparse, json, time

def parse():
    ap = argparse.ArgumentParser(description='Karyotype inference from Hi-C (exact-NLL search).')
    ap.add_argument('problem', help='problem folder from prep_data_for_py()')
    ap.add_argument('out', help='output folder (created)')
    ap.add_argument('--method', default='pt', choices=['pt', 'greedy'],
                    help='pt: parallel tempering (default); greedy: all chains at T = 0')
    ap.add_argument('--chains', type=int, default=16, help='number of chains (= PT rungs) [16]')
    ap.add_argument('--steps', type=int, default=20000, help='proposals per chain, rounded down to whole blocks [20000]')
    ap.add_argument('--block', type=int, default=500, help='proposals per chain between PT swaps [500]')
    ap.add_argument('--tmin', type=float, default=1.0, help='lowest PT temperature [1]')
    ap.add_argument('--tmax', type=float, default=1000.0, help='highest PT temperature [1000]')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--device', default='auto', choices=['auto', 'gpu', 'cpu'],
                    help='auto: GPU if available, else CPU with a warning; gpu: fail without a GPU [auto]')
    ap.add_argument('--truth', action='store_true',
                    help='also score PROBLEM/truth_gwalks.json (simulations): NLL gap and exact recovery')
    return ap.parse_args()

def main():
    a = parse()
    if a.device == 'cpu':
        os.environ['JAX_PLATFORMS'] = 'cpu'
    os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')  # don't grab 75% of a shared GPU
    import jax; jax.config.update('jax_enable_x64', True)
    import numpy as np
    import model, walks, search
    log = lambda s: print(s, flush=True)

    devs = jax.devices()
    log(f'devices: {devs}')
    if a.device != 'cpu' and devs[0].platform != 'gpu':
        if a.device == 'gpu':
            sys.exit('error: --device gpu but JAX found no GPU (CUDA init failed or CUDA_VISIBLE_DEVICES empty)')
        log('WARNING: no GPU found, running on CPU (~20x slower)')

    folder = a.problem.rstrip('/')
    P = model.load_problem(folder)
    log(f'problem {folder}: {P.ntile} tiles, {len(P.tc_uid)} tile copies, {len(P.var_nodes)} variable nodes')
    nll_truth = None
    if a.truth:
        sn, circ = walks.read_gwalks(folder + '/truth_gwalks.json')
        nll_truth = float(model.make_model(P).loss_perms(walks.walks_to_perms(P, sn, circ)))
        log(f'truth NLL {nll_truth:.6g} (reported values below are gaps to it)')

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, 'trace.tsv'), 'w') as tr:
        res = search.run(P, method=a.method, chains=a.chains, steps=a.steps, block=a.block, tmin=a.tmin,
                         tmax=a.tmax, seed=a.seed, nll_ref=nll_truth or 0.0, trace=tr, log=log)

    wk, circular = walks.perms_to_walks(P, res['perms'])
    with open(os.path.join(a.out, 'gwalks.json'), 'w') as f:
        json.dump({'snode.ids': [[int(x) for x in w] for w in wk], 'circular': [bool(c) for c in circular]}, f)
    np.save(os.path.join(a.out, 'best_b.npy'), res['b'])
    summ = dict(problem=os.path.abspath(folder), device=str(devs[0].platform), n_walks=len(wk),
                n_circular=int(sum(circular)), **res['summary'])
    if nll_truth is not None:
        summ.update(nll_truth=nll_truth, gap_to_truth=res['nll'] - nll_truth,
                    exact_recovery=bool(walks.same_karyotype(wk, circular, sn, circ)))
    json.dump(summ, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    log('done: ' + json.dumps({k: v for k, v in summ.items() if k in
        ('best_nll', 'gap_to_truth', 'exact_recovery', 'n_walks', 'n_circular', 'wall_s', 'ms_per_proposal')}))

if __name__ == '__main__':
    main()
