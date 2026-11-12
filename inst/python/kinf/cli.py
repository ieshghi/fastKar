"""Command-line entry point: infer a karyotype (wiring of a genome graph) from Hi-C by search on the exact NLL.

Input: a problem folder written by R/prep_data_for_py.R (or from R, use infer_karyotype() in R/infer_karyotype.R,
which writes the folder, calls this script and reads the result back as a gWalk).
Output (OUT): gwalks.json (best karyotype: signed snode ids on the gGraph passed to prep_data_for_py, circular
flags; read in R with read_kinf_gwalks()), summary.json, trace.tsv (best and per-rung NLL per block), best_b.npy;
with --topk K > 1 also topk.json (the K best distinct karyotypes visited: list of {nll, snode.ids, circular});
with --history N > 0 also history.tsv (NLL of every chain every N proposals: steps, chain, T, nll; minus the truth
NLL with --truth), for convergence plots.

Usage: python cli.py PROBLEM OUT [--method pt|greedy] [--chains 16] [--steps 20000] [--block 500]
       [--tmin 1] [--tmax 1000] [--seed 0] [--device auto|gpu|cpu] [--truth] [--topk 1] [--rescale 0] [--history 0]
       [--init FILE] [--init-random 0] [--spec FILE]
--spec FILE: Python file defining make(g, hic) -> dict(terms=..., energy=... optional) (see optimize.optimize); replaces the
Hi-C NLL, the problem folder then needs only nodes/internal_edges/external_edges (hic = None); incompatible with --rescale/--truth.
CPU use is capped to FASTKAR_KINF_CORES CPUs (default 8; 0 = no cap) before JAX starts.
Run with PYTHONPATH unset (env -u PYTHONPATH); entries from PYTHONPATH are dropped from sys.path here as well.
"""
import os, sys
def _cap_cpus(n=None):
    """narrow CPU affinity to a block of FASTKAR_KINF_CORES CPUs (chosen by pid) so XLA's thread pools stay under ulimit -u;
    an affinity already at or below the cap (taskset, infer_karyotype(cores =)) is left alone. Before the first JAX computation."""
    n = int(os.environ.get('FASTKAR_KINF_CORES', 8)) if n is None else int(n)
    if n <= 0 or not hasattr(os, 'sched_getaffinity'):
        return
    cur = sorted(os.sched_getaffinity(0))
    if len(cur) > n:
        k = os.getpid() % (len(cur) // n)
        os.sched_setaffinity(0, cur[k * n:(k + 1) * n])
_cap_cpus()
_pp = [os.path.abspath(p) for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p]
sys.path[:] = [p for p in sys.path if os.path.abspath(p or '.') not in _pp]  # foreign PYTHONPATH breaks numpy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import argparse, json, importlib.util

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
    ap.add_argument('--topk', type=int, default=1,
                    help='also write the K best distinct karyotypes visited by any chain to topk.json [1: off]')
    ap.add_argument('--init', metavar='FILE',
                    help='warm start: best_b.npy (or a stack of b rows, .npy) or a gwalks.json / topk.json from a '
                         'previous run on the same graph and target region (any pix.size); chains take the given '
                         'wirings in order, round robin')
    ap.add_argument('--init-random', type=int, default=0, metavar='N',
                    help='with --init: leave the last N chains on random starts (control) [0]')
    ap.add_argument('--history', type=int, default=0, metavar='N',
                    help='write the NLL of every chain every N proposals to history.tsv [0: off]')
    ap.add_argument('--rescale', type=int, default=0, metavar='N',
                    help='rescale depth by observed total / median expected total over N random wirings, as '
                         'fastKar bestfit_search(rescale = TRUE) does over its sampled candidates [0: off]')
    ap.add_argument('--spec', metavar='FILE',
                    help='Python file defining make(g, hic) -> dict(terms=..., energy=... optional) instead of the Hi-C NLL; '
                         'the problem folder then needs no Hi-C files; not with --rescale/--truth')
    ap.add_argument('--truth', action='store_true',
                    help='also score PROBLEM/truth_gwalks.json (simulations): NLL gap and exact recovery')
    return ap.parse_args()

def main():
    a = parse()
    if a.spec and (a.rescale > 0 or a.truth):
        sys.exit('error: --spec cannot be combined with --rescale or --truth')
    if a.device == 'cpu':
        os.environ['JAX_PLATFORMS'] = 'cpu'
    os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')  # don't grab 75% of a shared GPU
    import jax; jax.config.update('jax_enable_x64', True)
    import numpy as np
    import problem as pr, optimize as opt
    log = lambda s: print(s, flush=True)

    devs = jax.devices()
    log(f'devices: {devs}')
    if a.device != 'cpu' and devs[0].platform != 'gpu':
        if a.device == 'gpu':
            sys.exit('error: --device gpu but JAX found no GPU (CUDA init failed or CUDA_VISIBLE_DEVICES empty)')
        log('WARNING: no GPU found, running on CPU (~20x slower)')

    folder = a.problem.rstrip('/')
    g = pr.Graph.load(folder)
    hic = None if a.spec else pr.HiC.load(folder, g)
    log(f'problem {folder}: {len(g.nodes)} node ends, {len(g.vnodes)} variable nodes' + ('' if hic is None else f', {hic.nbin} data bins (scored), {hic.ntile} tiles, {len(hic.tc_uid)} tile copies'))
    depth_ratio = 1.0
    if a.rescale > 0:
        depth_ratio = pr.rescale_depth(g, hic, a.rescale, a.seed)
        log(f'depth rescaled by {depth_ratio:.4g} (median over {a.rescale} random wirings): {hic.depth:.4g}')
    nll_truth = None
    if a.truth:
        sn, circ = pr.read_gwalks(folder + '/truth_gwalks.json')
        nll_truth = float(pr.nll(g, hic, pr.b_from_perms(g, pr.walks_to_perms(g, sn, circ))))
        log(f'truth NLL {nll_truth:.6g} (reported values below are gaps to it)')
    ref = nll_truth or 0.0
    if a.spec:
        s = importlib.util.spec_from_file_location('kinf_spec', a.spec); m = importlib.util.module_from_spec(s)
        s.loader.exec_module(m)
        spec = m.make(g, hic)
    else:
        spec = dict(terms=[pr.nll_term(g, hic)])

    init = load_init(g, a.init) if a.init else None
    if init:
        try:
            init = [pr.check_b(g, b) for b in init]
        except ValueError as e:
            sys.exit(f'error: --init {a.init}: {e}')
    os.makedirs(a.out, exist_ok=True)
    T = 0.0 if a.method == 'greedy' else np.geomspace(a.tmin, a.tmax, a.chains)
    res = opt.optimize(g, spec['terms'], spec.get('energy'), T=T, chains=a.chains, steps=a.steps, block=a.block,
                       seed=a.seed, init=init, init_random=a.init_random, history_every=a.history, topk=a.topk, log=log)
    best, tr, H = res['best'], res['trace'], res['history']
    temps = np.zeros(a.chains) if a.method == 'greedy' else T
    with open(os.path.join(a.out, 'trace.tsv'), 'w') as f:
        if tr:
            f.write('steps\tbest\t' + '\t'.join(f'T{t:g}' for t in temps) + '\n')
            for s, b, cur in zip(tr['steps'], tr['best'], tr['cur']):
                f.write(f'{s}\t{b - ref:.6g}\t' + '\t'.join(f'{x - ref:.6g}' for x in cur) + '\n')
    if a.history > 0:
        with open(os.path.join(a.out, 'history.tsv'), 'w') as f:
            if H:
                f.write('steps\tchain\tT\tnll\n')
                f.write(''.join(f'{s}\t{c}\t{t:g}\t{e - ref:.10g}\n' for s, c, t, e in zip(H['step'], H['chain'], H['T'], H['E'])))

    wk, circular = pr.perms_to_walks(g, pr.perms_from_b(g, best['b']))
    pr.write_gwalks(os.path.join(a.out, 'gwalks.json'), wk, circular)
    np.save(os.path.join(a.out, 'best_b.npy'), best['b'])
    if a.topk > 1:
        top = pr.distinct_karyotypes(g, res['top'], a.topk)
        with open(os.path.join(a.out, 'topk.json'), 'w') as f:
            json.dump([{'nll': t['E'], 'snode.ids': [[int(x) for x in w] for w in t['walks']],
                        'circular': [bool(c) for c in t['circular']]} for t in top], f)
        log(f'top {len(top)} karyotypes, NLL - best: ' + ' '.join(f"{t['E'] - best['E']:.3g}" for t in top))
    r = res['summary']
    summ = dict(problem=os.path.abspath(folder), device=str(devs[0].platform), depth=None if hic is None else hic.depth,
                depth_ratio=depth_ratio, n_walks=len(wk), n_circular=int(sum(circular)),
                method=a.method, seed=a.seed, chains=r['chains'], steps=r['steps'], block=r['block'],
                temps=[float(t) for t in temps], best_nll=best['E'], best_rel=best['E'] - ref,
                n_var_nodes=r['n_var_nodes'], n_moves=r['n_moves'])
    summ.update({k: r[k] for k in ('acc_rate_per_rung', 'swap_rate', 'wall_s', 'ms_per_proposal') if k in r})
    if nll_truth is not None:
        summ.update(nll_truth=nll_truth, gap_to_truth=best['E'] - nll_truth,
                    exact_recovery=bool(pr.same_karyotype(wk, circular, sn, circ)))
    json.dump(summ, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    log('done: ' + json.dumps({k: v for k, v in summ.items() if k in
        ('best_nll', 'gap_to_truth', 'exact_recovery', 'n_walks', 'n_circular', 'wall_s', 'ms_per_proposal')}))

def load_init(g, path):
    """warm-start wirings from a file: .npy = one b (1-D) or a stack (2-D); .json = gwalks.json ({snode.ids,
    circular}) or topk.json (a list of them), converted with problem.walks_to_perms."""
    import numpy as np
    import problem as pr
    if path.endswith('.npy'):
        b = np.load(path)
        return [b] if b.ndim == 1 else list(b)
    x = json.load(open(path))
    out = []
    for k in (x if isinstance(x, list) else [x]):
        try:
            out.append(pr.b_from_perms(g, pr.walks_to_perms(g, [[int(v) for v in w] for w in k['snode.ids']],
                                                            [bool(c) for c in k['circular']])))
        except (ValueError, AssertionError, KeyError) as e:
            sys.exit(f'error: --init {path}: walks do not fit this problem ({e}); same graph and target region?')
    return out

if __name__ == '__main__':
    main()
