"""kinf problem: the graph (wiring space), wirings and walks, the walk structure every score sees (Kar), the Hi-C data and
forward model, and the exact NLL (full, and incremental as a Term for optimize.py).

Encoding (R/backwards.R prep_data_for_py): every copy of every half-node has a unique id (uid, 0-based). A left half has its
OUTER end at its genomic left end and its INNER end at the right end; a right half the opposite. Junction copies (ext)
join OUTER ends of two uids; at a node with CN > 1 the free variable is a CN x CN 0/1 permutation joining INNER ends of left
copies (rows) to right copies (cols). End index: 2*uid = genomic left end of the half, 2*uid+1 = right end.
A wiring is b: for each row of var_a (left inner end of every variable node, concatenated), the right inner end it is joined to.
"""
import json
import os

def cap_cpus(n=None):
    """narrow CPU affinity to a block of FASTKAR_KINF_CORES CPUs (chosen by pid) so XLA's thread pools stay under ulimit -u;
    an affinity already at or below the cap (taskset, infer_karyotype(cores =)) is left alone. Must run before the first JAX
    computation (the backend starts lazily)."""
    n = int(os.environ.get('FASTKAR_KINF_CORES', 8)) if n is None else int(n)
    if n <= 0 or not hasattr(os, 'sched_getaffinity'):
        return
    cur = sorted(os.sched_getaffinity(0))
    if len(cur) > n:
        k = os.getpid() % (len(cur) // n)
        os.sched_setaffinity(0, cur[k * n:(k + 1) * n])
cap_cpus()  # on import: every user of the package (cli.py, scripts, tests) is capped

from collections import defaultdict
from typing import Any, Callable, NamedTuple
import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp
from jax.scipy.special import gammaln

#### graph

def outer_end(side, uid):  # genomic end index at which junctions attach
    return 2 * uid + side

def inner_end(side, uid):  # end joined to the other half of the node
    return 2 * uid + (1 - side)

class Graph:
    """nu, node_of_uid, side_of_uid (0 left, 1 right), width_uid, ext (junction uid pairs), inner_fixed (CN 1 links),
    var_nodes {node: (left uids, right uids)}, nodes (annotation table), and derived: vnodes (sorted variable nodes),
    cn, start (block offsets in b), var_a, var_b {node: right inner ends}, lp0 (fixed link partner per end, -1 none),
    moves (transpositions of rows within a node block of b), move_node."""
    def __init__(self, nodes, ext, inner_fixed, var_nodes):
        self.nodes = nodes
        self.nu = len(nodes)
        self.node_of_uid = nodes['node_id'].values.astype(int)
        self.side_of_uid = (nodes['side'].values == 'right').astype(int)
        self.width_uid = nodes['width'].values.astype(float)
        self.ext, self.inner_fixed, self.var_nodes = ext, inner_fixed, var_nodes
        self.vnodes = sorted(var_nodes)
        self.cn = np.array([len(var_nodes[n][0]) for n in self.vnodes], dtype=int)
        self.start = np.concatenate([[0], np.cumsum(self.cn)])[:-1].astype(int)
        self.var_a = np.concatenate([inner_end(0, var_nodes[n][0]) for n in self.vnodes] or [np.zeros(0, int)])
        self.var_b = {n: inner_end(1, var_nodes[n][1]) for n in self.vnodes}
        self.moves = np.array([(s + a, s + c) for s, k in zip(self.start, self.cn) for a in range(k)
                               for c in range(a + 1, k)], dtype=int).reshape(-1, 2)
        self.move_node = np.repeat(np.arange(len(self.vnodes)), self.cn * (self.cn - 1) // 2)
        N = 2 * self.nu; side = self.side_of_uid
        lp0 = np.full(N, -1)
        def link(a, b):
            assert a != b and lp0[a] < 0 and lp0[b] < 0, 'end with more than one link'
            lp0[a] = b; lp0[b] = a
        for n1, n2 in ext:
            link(outer_end(side[n1], n1), outer_end(side[n2], n2))
        for l, r in inner_fixed:
            link(inner_end(side[l], l), inner_end(side[r], r))
        self.lp0 = lp0

    @classmethod
    def _from_parts(cls, nodes, ext, ie):
        ext = np.asarray(ext).reshape(-1, 2).astype(int)
        inner_fixed = ie[ie['copies'] == 1][['left', 'right']].values.astype(int)
        var_nodes = {int(n): (g['left'].values.astype(int), g['right'].values.astype(int))
                     for n, g in ie[ie['copies'] > 1].groupby('node_id')}
        return cls(nodes, ext, inner_fixed, var_nodes)

    @classmethod
    def load(cls, folder):
        """needs only nodes.csv, internal_edges.csv, external_edges.csv of the problem folder"""
        folder = folder.rstrip('/') + '/'
        rd = lambda n: pd.read_csv(folder + n + '.csv', index_col=0)
        nodes = rd('nodes').sort_values('unique.id')
        assert (nodes['unique.id'].values - 1 == np.arange(len(nodes))).all()
        ie = rd('internal_edges'); ee = rd('external_edges')
        ie = ie.assign(left=ie['left'] - 1, right=ie['right'] - 1)
        return cls._from_parts(nodes, ee[['n1', 'n2']].values.astype(int) - 1, ie)

    @classmethod
    def from_tables(cls, nodes, junctions):
        """nodes: DataFrame node_id, cn, width (+ optional seqnames, start, end); junctions: DataFrame n1, side1, n2, side2, cn
        (side = 'left'/'right' of the node the junction attaches to). Same uid / junction-copy assignment as prep_data_for_py:
        node by node, left copies then right copies (width: left floor(w/2), right the rest); junction rows repeated cn times,
        the k-th use of (node, side) takes its k-th copy, all n1 ends of all rows before all n2 ends."""
        rows = []
        for _, r in nodes.iterrows():
            wl = int(r['width']) // 2
            for side, w in (('left', wl), ('right', int(r['width']) - wl)):
                for _ in range(int(r['cn'])):
                    d = dict(node_id=int(r['node_id']), side=side, width=float(w))
                    if 'seqnames' in r and 'start' in r and 'end' in r:
                        d.update(seqnames=r['seqnames'], start=r['start'] if side == 'left' else r['end'] - w + 1,
                                 end=r['start'] + w - 1 if side == 'left' else r['end'])
                    rows.append(d)
        nd = pd.DataFrame(rows); nd['unique.id'] = np.arange(1, len(nd) + 1)
        uids = defaultdict(list)
        for u, (n, s) in enumerate(zip(nd['node_id'], nd['side'])):
            uids[(n, s)].append(u)
        jr = junctions.loc[junctions.index.repeat(junctions['cn'].astype(int))]
        occ = defaultdict(int); ends = {}
        for role in ('1', '2'):
            for e, (n, s) in enumerate(zip(jr['n' + role], jr['side' + role])):
                k = (int(n), s)
                assert occ[k] < len(uids[k]), f'more junction copies than half-node copies at {k}'
                ends[(role, e)] = uids[k][occ[k]]; occ[k] += 1
        ext = np.array([(ends[('1', e)], ends[('2', e)]) for e in range(len(jr))], dtype=int)
        left = nd[nd['side'] == 'left']['node_id']; l_u = np.nonzero(nd['side'].values == 'left')[0]
        r_u = np.nonzero(nd['side'].values == 'right')[0]
        ie = pd.DataFrame(dict(node_id=left.values, left=l_u, right=r_u))
        ie['copies'] = ie.groupby('node_id')['node_id'].transform('size')
        return cls._from_parts(nd, ext, ie)

#### wirings

def b_from_perms(g, perms):
    return jnp.concatenate([jnp.asarray(g.var_b[n])[jnp.argmax(jnp.asarray(perms[n]), axis=1)] for n in g.vnodes]) \
        if g.vnodes else jnp.zeros(0, int)

def perms_from_b(g, b):
    b = np.asarray(b); out = {}
    for n, s, k in zip(g.vnodes, g.start, g.cn):
        pos = {e: c for c, e in enumerate(g.var_b[n])}
        M = np.zeros((k, k)); M[np.arange(k), [pos[e] for e in b[s:s + k]]] = 1
        out[n] = M
    return out

def check_b(g, b):
    """b as a numpy int array, or ValueError if it is not a wiring of g"""
    b = np.asarray(b)
    if b.ndim != 1 or len(b) != g.cn.sum():
        raise ValueError(f'init wiring has shape {b.shape}; this problem needs {g.cn.sum()} entries '
                         f'({len(g.vnodes)} variable nodes): was it made for the same graph and target region?')
    for n, s, k in zip(g.vnodes, g.start, g.cn):
        if not np.array_equal(np.sort(b[s:s + k]), np.sort(g.var_b[n])):
            raise ValueError(f'init wiring is not a permutation at node {n}: same graph and target region?')
    return b

def random_b(g, key, k):
    """k uniformly random wirings, (k, len(var_a))"""
    keys = jax.random.split(key, len(g.vnodes))
    blocks = [jnp.asarray(g.var_b[n])[jax.vmap(lambda q: jax.random.permutation(q, c))(jax.random.split(kk, k))]
              for kk, n, c in zip(keys, g.vnodes, g.cn)]
    return jnp.concatenate(blocks, axis=1) if blocks else jnp.zeros((k, 0), int)

#### walks

def read_gwalks(path):
    d = json.load(open(path))
    sn = [[w] if isinstance(w, int) else list(w) for w in d['snode.ids']]  # jsonlite writes length-one vectors as numbers
    circ = d['circular']
    circ = [circ] if isinstance(circ, bool) else list(circ)
    return sn, circ

def write_gwalks(path, walks, circular):
    with open(path, 'w') as f:
        json.dump({'snode.ids': [[int(x) for x in w] for w in walks], 'circular': [bool(c) for c in circular]}, f)

def walks_to_perms(g, snode_ids, circular):
    """hard permutation matrices {node: (CN,CN) 0/1, rows = left copies, cols = right copies} (only nodes with CN>1).
    Copies of a junction (and of a node) are exchangeable; unused copies are assigned greedily in walk order."""
    side = g.side_of_uid; node = g.node_of_uid
    pool = defaultdict(list)  # junction type -> list of (n1, n2) uid pairs still unused
    for n1, n2 in g.ext:
        key = tuple(sorted([(int(node[n1]), int(side[n1])), (int(node[n2]), int(side[n2]))]))
        pool[key].append((int(n1), int(n2)))
    used = set()
    visits = [[[abs(int(s)), 1 if s > 0 else -1, None, None] for s in w] for w in snode_ids]  # [node, sign, entry_uid, exit_uid]
    def junction(va, vb):
        a_side = 1 if va[1] > 0 else 0       # exit through the right (+) / left (-) outer end
        b_side = 0 if vb[1] > 0 else 1       # enter through the left (+) / right (-) outer end
        la, lb = (va[0], a_side), (vb[0], b_side)
        key = tuple(sorted([la, lb]))
        if not pool[key]:
            raise ValueError(f'no unused junction copy of type {key}')
        n1, n2 = pool[key].pop()
        if (int(node[n1]), int(side[n1])) != la:
            n1, n2 = n2, n1
        assert (int(node[n1]), int(side[n1])) == la and (int(node[n2]), int(side[n2])) == lb
        va[3] = n1; vb[2] = n2
        used.add(n1); used.add(n2)
    for vs, c in zip(visits, circular):
        for k in range(len(vs) - 1):
            junction(vs[k], vs[k + 1])
        if c:
            junction(vs[-1], vs[0])
    left_over = sum(len(v) for v in pool.values())
    if left_over:
        raise ValueError(f'{left_over} junction copies unused by the karyotype')
    free = defaultdict(list)  # (node, side) -> uids with no junction
    for u in range(g.nu):
        if u not in used:
            free[(int(node[u]), int(side[u]))].append(u)
    ext_uids = set(int(x) for x in g.ext.ravel())
    free = {k: [u for u in v if u not in ext_uids] for k, v in free.items()}
    link = defaultdict(list)  # node -> list of (left_uid, right_uid)
    for vs in visits:
        for v in vs:
            n, s = v[0], v[1]
            for idx, sd in ((2, 0 if s > 0 else 1), (3, 1 if s > 0 else 0)):
                if v[idx] is None:
                    v[idx] = free[(n, sd)].pop()
                    used.add(v[idx])
            lu, ru = (v[2], v[3]) if s > 0 else (v[3], v[2])
            assert side[lu] == 0 and side[ru] == 1 and node[lu] == n and node[ru] == n
            link[n].append((lu, ru))
    assert len(used) == g.nu, 'not every half-node copy is used by the karyotype'
    perms = {}
    for n, (lefts, rights) in g.var_nodes.items():
        M = np.zeros((len(lefts), len(lefts)))
        li = {int(u): i for i, u in enumerate(lefts)}; ri = {int(u): i for i, u in enumerate(rights)}
        assert len(link[n]) == len(lefts)
        for lu, ru in link[n]:
            M[li[lu], ri[ru]] += 1
        assert (M.sum(0) == 1).all() and (M.sum(1) == 1).all()
        perms[n] = M
    for l, r in g.inner_fixed:  # CN 1 nodes: only one possible link
        assert (int(l), int(r)) in set(link[int(node[l])])
    return perms

def perms_to_walks(g, perms):
    """walks (lists of signed node ids) and circular flags from hard permutations; linear walks start at an open end"""
    side = g.side_of_uid; node = g.node_of_uid
    outer = {}
    for n1, n2 in g.ext:
        outer[int(n1)] = int(n2); outer[int(n2)] = int(n1)
    inner = {}
    for l, r in g.inner_fixed:
        inner[int(l)] = int(r); inner[int(r)] = int(l)
    for n, (lefts, rights) in g.var_nodes.items():
        M = np.asarray(perms[n])
        for i, j in zip(*np.nonzero(M > 0.5)):
            inner[int(lefts[i])] = int(rights[j]); inner[int(rights[j])] = int(lefts[i])
    visited = set()
    walks, circ = [], []
    def trace(u0):
        out = []; u = u0
        while True:
            out.append(int(node[u]) * (1 if side[u] == 0 else -1))
            visited.add(u); x = inner[u]; visited.add(x)
            if x not in outer: return out, False
            u = outer[x]
            if u == u0: return out, True
    for u in range(g.nu):
        if u not in outer and u not in visited:
            w, c = trace(u); walks.append(w); circ.append(c)
    for u in range(g.nu):
        if u not in visited:
            w, c = trace(u); walks.append(w); circ.append(c)
    return walks, circ

def canonical(walk, circular):
    rev = [-s for s in reversed(walk)]
    if not circular:
        return (False, min(tuple(walk), tuple(rev)))
    cands = []
    for w in (walk, rev):
        for k in range(len(w)):
            cands.append(tuple(w[k:] + w[:k]))
    return (True, min(cands))

def same_karyotype(w1, c1, w2, c2):
    return sorted(canonical(w, c) for w, c in zip(w1, c1)) == sorted(canonical(w, c) for w, c in zip(w2, c2))

#### walk structure: what every score sees

def walkinfo(g, b):
    """exact walk structure of a hard wiring b, without Floyd-Warshall: every end has at most one link partner, so the genome
    is a set of disjoint walks (paths and cycles), found by pointer jumping (log2 N gather rounds; jittable, vmappable).
    State a = 'enter half-copy a//2 at end a, leave at end a^1'; succ(a) = link partner of a^1 (or sink N).
    Returns per end: coord (bp along its walk), wid (walk id), L (cycle length), iscyc."""
    nu = g.nu; N = 2 * nu
    w = jnp.asarray(np.append(np.repeat(g.width_uid, 2), 0.0))  # width of the half crossed in state a; sink 0
    K = int(np.ceil(np.log2(N + 1))) + 1
    idx = jnp.arange(N + 1)
    def jump(s, val, comb):  # K rounds of pointer doubling: val(a) <- comb(val(a), val(s(a))), s <- s(s)
        def body(_, c):
            s, v = c
            return s[s], comb(v, jax.tree.map(lambda x: x[s], v))
        return jax.lax.fori_loop(0, K, body, (s, val))
    lp0 = jnp.asarray(g.lp0); var_a = jnp.asarray(g.var_a)
    lp = lp0.at[var_a].set(b).at[b].set(var_a) if g.vnodes else lp0
    nxt = lp[jnp.arange(N) ^ 1]
    succ = jnp.append(jnp.where(nxt >= 0, nxt, N), N)
    # pass 1: cycle detection and the minimum state of each cycle
    sK, mn = jump(succ, jnp.where(idx < N, idx, N + 1), jnp.minimum)
    cyc = sK != N
    # break each cycle in front of its minimum state, then every state lies on a path
    succ2 = jnp.where(cyc & (succ == mn), N, succ)
    # pass 2: R(a) = width from state a to the tail inclusive; T(a) = tail state of a's path (T(sink) = N)
    _, (R, T) = jump(succ2, (w, idx), lambda v, vs: (v[0] + vs[0], jnp.where(vs[1] == N, v[1], vs[1])))
    R, T, cyc = R[:N], T[:N], cyc[:N]
    circ = jax.ops.segment_max(R, T, num_segments=N)[T]  # cycle length (= R at the head)
    # per half-copy pick the direction with the smaller tail (consistent along each walk)
    s0 = 2 * jnp.arange(nu)
    c = jnp.where(T[s0] <= T[s0 + 1], s0, s0 + 1)  # canonical state per uid
    ce = c ^ 1
    coord = jnp.zeros(N).at[c].set(R[c]).at[ce].set(R[c] - w[c])  # entry end at R, exit end at R - width
    wid = jnp.zeros(N, int).at[c].set(T[c]).at[ce].set(T[c])
    L = jnp.zeros(N).at[c].set(circ[c]).at[ce].set(circ[c])
    iscyc = jnp.zeros(N, bool).at[c].set(cyc[c]).at[ce].set(cyc[c])
    return coord, wid, L, iscyc

class Kar(NamedTuple):
    """a wiring b and its walk structure (per end: coord, wid, L, iscyc; see walkinfo)"""
    b: Any; coord: Any; wid: Any; L: Any; iscyc: Any

def kar(g, b):
    return Kar(b, *walkinfo(g, b))

class Term(NamedTuple):
    """a score with an incremental form (optimize.py): init(kar) -> aux; value(aux) -> scalar;
    update(aux, kar_old, kar_new, i, j) -> aux after swapping rows i, j of b"""
    init: Callable; value: Callable; update: Callable

#### Hi-C forward model (port of fastKar R/forward.R + src/forward_kernel.cpp; constants from data/fastkar_model.json)

_MODEL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'fastkar_model.json')
_CACHE = {}

def get_constants(path=_MODEL_PATH):
    if path not in _CACHE:
        m = json.load(open(path))
        s = m['distance_spline']
        assert s['method'] == 3, 'only R fmm splines are supported'
        _CACHE[path] = dict(
            x=np.asarray(s['x'], dtype=np.float64), y=np.asarray(s['y'], dtype=np.float64),
            b=np.asarray(s['b'], dtype=np.float64), c=np.asarray(s['c'], dtype=np.float64),
            d=np.asarray(s['d'], dtype=np.float64),
            lmin=float(m['lmin']), slope=float(m['slope']), maxval=float(m['maxval']),
            diag_a=float(m['diag_spline']['a']), diag_b=float(m['diag_spline']['b']),
            rho=float(m['interchrom_density']), raw=m)
        # knots at log(step * k) for integer k (the fastKar spline is fitted on a bp grid): a table from
        # floor(d / step) to the knot index replaces the binary search
        x = _CACHE[path]['x']; step = np.exp(x[0]); k = np.exp(x) / step
        if np.allclose(k, np.round(k), rtol=1e-9, atol=0):
            with np.errstate(divide='ignore'):
                grid = np.log(step * np.arange(int(np.round(k[-1])) + 1))
            _CACHE[path]['knot_step'] = float(step)
            _CACHE[path]['knot_table'] = np.searchsorted(x, grid, side='right') - 1
    return _CACHE[path]

def _spline(u, K):
    """R splinefun (fmm) evaluation at u; outside the knots the first/last cubic is extended (as R does)."""
    x = jnp.asarray(K['x'])
    n = x.shape[0]
    if 'knot_step' in K:  # = searchsorted(x, u, 'right') - 1, without the binary search
        tab = jnp.asarray(K['knot_table'])
        q = jnp.clip(jnp.floor(jnp.exp(jnp.minimum(u, 100.0)) / K['knot_step']), 0, tab.shape[0] - 1).astype(jnp.int64)
        i = jnp.clip(tab[q], 0, n - 1)
        for _ in range(2):  # fix rounding at knots (floor and table can each be one off)
            i = i - (x[i] > u) + ((i + 1 < n) & (x[jnp.minimum(i + 1, n - 1)] <= u))
            i = jnp.clip(i, 0, n - 1)
    else:
        i = jnp.clip(jnp.searchsorted(x, u, side='right') - 1, 0, n - 1)
    dx = u - x[i]
    return jnp.asarray(K['y'])[i] + dx * (jnp.asarray(K['b'])[i] + dx * (jnp.asarray(K['c'])[i] + dx * jnp.asarray(K['d'])[i]))

def distance_density(d, K=None):
    """contact density per bp^2 at distance d > 0 (bp): spline above exp(lmin), power law below, 1/d tail above maxval"""
    K = K or get_constants()
    d = jnp.asarray(d)
    ld = jnp.log(d)
    lmin = K['lmin']
    out = _spline(jnp.maximum(ld, lmin), K) + K['slope'] * jnp.minimum(ld - lmin, 0.0)
    dens = jnp.exp(out)
    dmax = K['maxval']
    dens_max = jnp.exp(_spline(jnp.log(dmax), K))
    return jnp.where(d > dmax, dens_max * dmax / d, dens)

def diag_contacts(widthprod, K=None):
    """contacts in the diagonal pixel of a tile: exp(diag_spline(log(widthprod)))"""
    K = K or get_constants()
    return jnp.exp(K['diag_a'] + K['diag_b'] * jnp.log(widthprod))

def pair_contacts(d, widthprod, dcut=1e10, dfloor=1e-3, K=None):
    """widthprod * distance_density(d), times the smooth gate 1/(1+(d/dcut)^4) that removes pairs on different walks"""
    d = jnp.maximum(d, dfloor)
    gate = 1.0 / (1.0 + (d / dcut) ** 4)
    return widthprod * distance_density(d, K) * gate

def contact_matrix(D, w, dcut=1e10, K=None):
    """per tile-copy pair contacts M[u,v] (symmetric, no mult/depth). D: (ntc,ntc) distances in bp, w: (ntc,) widths."""
    wp = w[:, None] * w[None, :]
    eye = jnp.eye(w.shape[0], dtype=bool)
    Ds = jnp.where(eye, 1.0, D)  # keep the unused branch finite
    off = pair_contacts(Ds, wp, dcut, K=K)
    return jnp.where(eye, diag_contacts(wp, K), off)

def contact_model(D, w, tile_of_copy, w_tile, cn_tile, depth, dcut=1e10, K=None):
    """expected counts per reference tile pair (full symmetric ntile x ntile; the R map is its upper triangle):
    value(p,q) = depth/2 * [ sum_{u in p, v in q} M[u,v] + w_p w_q rho cn_p cn_q ]."""
    K = K or get_constants()
    M = contact_matrix(D, w, dcut, K)
    ntile = w_tile.shape[0]
    seg = lambda X: jax.ops.segment_sum(X, tile_of_copy, num_segments=ntile)
    R = seg(seg(M).T).T
    bg = w_tile[:, None] * w_tile[None, :] * K['rho'] * (cn_tile[:, None] * cn_tile[None, :])
    return 0.5 * depth * (R + bg)

#### Hi-C data + likelihood

class HiC:
    """tiles (ntile, tc_uid/tc_tile/tc_w/offL/offR per tile copy, tile_w, tile_cn), data bins (nbin, bin_w, bin_of_tile),
    observed counts x_obs (upper triangle), depth, purity"""
    @classmethod
    def load(cls, folder, graph):
        h = cls(); folder = folder.rstrip('/') + '/'
        rd = lambda n: pd.read_csv(folder + n + '.csv', index_col=0)
        tt = rd('tiled_target'); dt = rd('dedup_tiles')
        h.ntile = len(tt)
        h.tc_uid = dt['unique.id'].values.astype(int) - 1
        h.tc_tile = dt['tile.id'].values.astype(int) - 1
        h.tc_w = dt['width'].values.astype(float)
        h.offL = dt['leftdist'].values.astype(float)
        h.offR = dt['rightdist'].values.astype(float)
        h.tile_w = tt['width'].values.astype(float)
        h.tile_cn = np.bincount(h.tc_tile, minlength=h.ntile).astype(float)  # graph CN of the tile = number of copies
        # scoring is on the data bins: every tile lies in one bin, the tile map is summed into bins (to_bins)
        if not os.path.exists(folder + 'bins.csv') or 'bin.id' not in tt.columns:
            raise ValueError(f'{folder}: no bins.csv / bin.id (exported before 2026-11-11, when kinf scored on tiles): '
                             're-export with fastKar::prep_data_for_py()')
        bins = rd('bins')
        h.nbin = len(bins)
        h.bin_w = bins['width'].values.astype(float)
        h.bin_of_tile = tt['bin.id'].values.astype(int) - 1
        gm = rd('gm_dat')
        assert (gm['i'].values <= gm['j'].values).all()
        h.x_obs = np.zeros((h.nbin, h.nbin))
        h.x_obs[gm['i'].values - 1, gm['j'].values - 1] = np.round(gm['value'].values)  # upper triangle; compdats rounds too
        h.depth = float(pd.read_csv(folder + 'depth.csv')['x'].iloc[0])
        h.purity = float(pd.read_csv(folder + 'purity.csv')['x'].iloc[0])
        assert h.purity == 1.0, 'purity < 1 not ported'
        return h

def nb_logpmf(x, mu, theta):
    mu_s = jnp.where(mu > 0, mu, 1.0)
    p = theta / (theta + mu_s)
    lp = gammaln(x + theta) - gammaln(theta) - gammaln(x + 1) + theta * jnp.log(p) + x * jnp.log1p(-p)
    return jnp.where(mu > 0, lp, 0.0)  # R compdats: logprob[y == 0] = 0

def to_bins(hic, Y):
    """tile map Y (full symmetric ntile x ntile, R map = upper triangle) -> full symmetric nbin x nbin map whose upper triangle
    is rebin_matrix() of Y's: off-diagonal bin pairs sum all tile pairs (p in A, q in B); a diagonal bin sums its tile pairs
    p <= q once (= (sum_{p,q in A} Y + trace_A) / 2)."""
    bot = jnp.asarray(hic.bin_of_tile); nb = hic.nbin
    seg = lambda X: jax.ops.segment_sum(X, bot, num_segments=nb)
    B = seg(seg(Y).T).T
    d = (jnp.diag(B) + seg(jnp.diag(Y))) / 2
    return B.at[jnp.arange(nb), jnp.arange(nb)].set(d)

def tc_dist(hic, kar, I, J, bigdist=1e12):
    """distances between tile copies I and J (index arrays, broadcast) on the wiring kar: min over the 4 end pairings of offset +
    end-end distance along the shared walk (shorter way round a cycle; bigdist between walks); |difference of offsets| within
    one half copy; 0 for I == J"""
    uid = jnp.asarray(hic.tc_uid); oL = jnp.asarray(hic.offL); oR = jnp.asarray(hic.offR)
    def ed(e1, e2):
        d = jnp.abs(kar.coord[e1] - kar.coord[e2])
        d = jnp.where(kar.iscyc[e1], jnp.minimum(d, kar.L[e1] - d), d)
        d = jnp.where(kar.wid[e1] == kar.wid[e2], d, bigdist)
        return jnp.where(e1 == e2, 0.0, d)
    Li, Ri, Lj, Rj = 2 * uid[I], 2 * uid[I] + 1, 2 * uid[J], 2 * uid[J] + 1
    m = jnp.minimum(jnp.minimum(ed(Li, Lj) + oL[I] + oL[J], ed(Ri, Rj) + oR[I] + oR[J]),
                    jnp.minimum(ed(Li, Rj) + oL[I] + oR[J], ed(Ri, Lj) + oR[I] + oL[J]))
    d = jnp.where(uid[I] == uid[J], jnp.abs(oL[I] - oL[J]), m)
    return jnp.where(I == J, 0.0, d)

def expected_map(hic, kar, dcut=1e10):
    """expected counts, full symmetric tile x tile map (R map = upper triangle), for the wiring kar"""
    a = jnp.arange(len(hic.tc_uid))
    return contact_model(tc_dist(hic, kar, a[:, None], a[None, :]), jnp.asarray(hic.tc_w), jnp.asarray(hic.tc_tile),
                         jnp.asarray(hic.tile_w), jnp.asarray(hic.tile_cn), hic.depth, dcut)

def nll_of_bins(hic, mu, theta=2.0, weighted=True):
    """compmaps conventions: NB NLL of the bin map mu over the upper triangle (incl. diagonal), pixel weight widthprod / area0
    (area0 = median(bin width)^2; weighted=False: 1), observed counts rounded (HiC.load)"""
    w = hic.bin_w; n = hic.nbin
    wp = jnp.asarray(w[:, None] * w[None, :]) / np.median(w) ** 2 if weighted else jnp.ones((n, n))
    tri = jnp.asarray(np.triu(np.ones((n, n), bool)))
    return -jnp.sum(jnp.where(tri, wp * nb_logpmf(jnp.asarray(hic.x_obs), mu, theta), 0.0))

def nll(g, hic, b, theta=2.0, weighted=True):
    """full NLL of the wiring b on the data bins"""
    return nll_of_bins(hic, to_bins(hic, expected_map(hic, kar(g, b))), theta, weighted)

def rescale_depth(g, hic, n, seed, batch=64):
    """multiply hic.depth (in place) by sum(observed) / median(sum(expected)) over n uniformly random wirings, both over the
    upper triangle of the data-bin map; returns the ratio. The expected map is linear in depth."""
    if len(g.vnodes) == 0:
        n = 1
    tot = jax.jit(jax.vmap(lambda b: jnp.sum(jnp.triu(to_bins(hic, expected_map(hic, kar(g, b)))))))
    r = random_b(g, jax.random.PRNGKey(seed + 7919), n)
    totals = np.concatenate([np.asarray(tot(r[i:i + batch])) for i in range(0, n, batch)])
    ratio = float(np.sum(hic.x_obs) / np.median(totals))
    hic.depth *= ratio
    return ratio

def nll_term(g, hic, dcut=1e10, theta=2.0, weighted=True, caps=(64, 128, 256, 512, 1024, 2048, 4096)):
    """the exact NLL as an incremental Term; aux = (mu = full symmetric expected BIN map (diagonal = the bin's upper-triangle
    sum), nll). A move swaps b[i], b[j] for two rows of one node. Only the walks through those rows change, and the set S of
    tile copies on them is the same before and after. Copies on different walks have no contacts (beyond the constant
    interchromosomal background), so only copy pairs in S x S change, hence only tile pairs within tiles(S): the update costs
    O(|S|^2) instead of O(ntc^2). S is gathered into a fixed-capacity buffer chosen by lax.switch over a ladder of capacities
    (so the update cannot be vmapped efficiently); the change of the tile map is summed into the bins of tiles(S) (to_bins)."""
    K = get_constants()
    var_a = jnp.asarray(g.var_a)
    ntc = len(hic.tc_uid); ntile = hic.ntile; nbin = hic.nbin; hd = 0.5 * hic.depth
    Lu = jnp.asarray(2 * hic.tc_uid); tc_w = jnp.asarray(hic.tc_w); tile_of = jnp.asarray(hic.tc_tile)
    w = hic.bin_w
    wp = jnp.asarray(w[:, None] * w[None, :] / np.median(w) ** 2) if weighted else jnp.ones((nbin, nbin))
    x = jnp.asarray(hic.x_obs)
    bin_ext = jnp.asarray(np.append(hic.bin_of_tile, nbin))  # padded tile index ntile -> padded bin nbin
    caps = [c for c in sorted(caps) if c < ntc] + [ntc]
    caps_arr = jnp.asarray(caps)

    def branch(c):
        def f(args):
            mask, cnt, ko, kn, mu = args
            S = jnp.nonzero(mask, size=c, fill_value=0)[0]
            valid = jnp.arange(c) < cnt
            I, J = S[:, None], S[None, :]
            wS = tc_w[S]
            dM = contact_matrix(tc_dist(hic, kn, I, J), wS, dcut, K) - contact_matrix(tc_dist(hic, ko, I, J), wS, dcut, K)
            dM = jnp.where(valid[:, None] & valid[None, :], dM, 0.0)
            tS = jnp.where(valid, tile_of[S], ntile)
            T = jnp.unique(tS, size=c, fill_value=ntile)
            loc = jnp.searchsorted(T, tS)
            dmu = hd * jnp.zeros((c, c)).at[loc[:, None], loc[None, :]].add(dM, mode='drop')  # tile map change on T x T
            # sum into bins: off-diagonal bin pairs take every tile pair, a diagonal bin its pairs p <= q once
            tb = bin_ext[T]
            U = jnp.unique(tb, size=c, fill_value=nbin)
            lb = jnp.searchsorted(U, tb)
            dB = jnp.zeros((c, c)).at[lb[:, None], lb[None, :]].add(dmu, mode='drop')
            dd = (jnp.diag(dB) + jnp.zeros(c).at[lb].add(jnp.diag(dmu), mode='drop')) / 2
            dB = dB.at[jnp.arange(c), jnp.arange(c)].set(dd)
            p, q = U[:, None], U[None, :]
            ok = (p < nbin) & (q < nbin) & (p <= q)
            pc, qc = jnp.minimum(p, nbin - 1), jnp.minimum(q, nbin - 1)
            mo = mu[pc, qc]; xx = x[pc, qc]
            dnll = -jnp.sum(jnp.where(ok, wp[pc, qc] * (nb_logpmf(xx, mo + dB, theta) - nb_logpmf(xx, mo, theta)), 0.0))
            return dnll, mu.at[p, q].add(dB, mode='drop')
        return f
    branches = [branch(c) for c in caps]

    def init(kar):
        mu = to_bins(hic, expected_map(hic, kar, dcut))
        return mu, nll_of_bins(hic, mu, theta, weighted)
    def update(aux, ko, kn, i, j):
        tcw = ko.wid[Lu]
        mask = (tcw == ko.wid[var_a[i]]) | (tcw == ko.wid[var_a[j]])
        cnt = mask.sum()
        dnll, mu2 = jax.lax.switch(jnp.searchsorted(caps_arr, cnt), branches, (mask, cnt, ko, kn, aux[0]))
        return mu2, aux[1] + dnll
    return Term(init, lambda aux: aux[1], update)


def distinct_karyotypes(g, top, k):
    """the k lowest-E entries of top (pairs (E, b)) that are distinct karyotypes: list of dict(E, b, walks, circular)"""
    out, seen = [], set()
    for E, b in sorted(top, key=lambda t: t[0]):
        wk, circ = perms_to_walks(g, perms_from_b(g, b))
        key = tuple(sorted(canonical(w, x) for w, x in zip(wk, circ)))
        if key in seen:
            continue
        seen.add(key)
        out.append(dict(E=float(E), b=np.asarray(b), walks=wk, circular=circ))
        if len(out) == k:
            break
    return out
