"""Karyotype <-> hard permutation conversions in the encoding of model.py / R/prep_data_for_py.R."""
import json
from collections import defaultdict
import numpy as np

SIDE = {0: 'left', 1: 'right'}

def read_gwalks(path):
    d = json.load(open(path))
    sn = d['snode.ids']
    # R's jsonlite writes length-one vectors as bare numbers
    sn = [[w] if isinstance(w, int) else list(w) for w in sn]
    circ = d['circular']
    circ = [circ] if isinstance(circ, bool) else list(circ)
    return sn, circ

def walks_to_perms(problem, snode_ids, circular):
    """hard permutation matrices {node: (CN,CN) 0/1, rows = left copies, cols = right copies} (only nodes with CN>1).
    Copies of a junction (and of a node) are exchangeable; unused copies are assigned greedily in walk order."""
    P = problem
    side = P.side_of_uid; node = P.node_of_uid
    pool = defaultdict(list)  # junction type -> list of (n1, n2) uid pairs still unused
    for n1, n2 in P.ext:
        key = tuple(sorted([(int(node[n1]), int(side[n1])), (int(node[n2]), int(side[n2]))]))
        pool[key].append((int(n1), int(n2)))
    used = set()
    visits = []  # per walk: list of [node, sign, entry_uid, exit_uid]
    for w in snode_ids:
        visits.append([[abs(int(s)), 1 if s > 0 else -1, None, None] for s in w])
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
    for u in range(P.nu):
        if u not in used:
            free[(int(node[u]), int(side[u]))].append(u)
    # remaining open ends (walk ends): need a uid without a junction
    ext_uids = set(int(x) for x in P.ext.ravel())
    free = {k: [u for u in v if u not in ext_uids] for k, v in free.items()}
    link = defaultdict(list)  # node -> list of (left_uid, right_uid)
    for vs in visits:
        for v in vs:
            n, s = v[0], v[1]
            for idx, (sd) in ((2, 0 if s > 0 else 1), (3, 1 if s > 0 else 0)):
                if v[idx] is None:
                    v[idx] = free[(n, sd)].pop()
                    used.add(v[idx])
            lu, ru = (v[2], v[3]) if s > 0 else (v[3], v[2])
            assert side[lu] == 0 and side[ru] == 1 and node[lu] == n and node[ru] == n
            link[n].append((lu, ru))
    assert len(used) == P.nu, 'not every half-node copy is used by the karyotype'
    perms = {}
    for n, (lefts, rights) in P.var_nodes.items():
        M = np.zeros((len(lefts), len(lefts)))
        li = {int(u): i for i, u in enumerate(lefts)}; ri = {int(u): i for i, u in enumerate(rights)}
        assert len(link[n]) == len(lefts)
        for lu, ru in link[n]:
            M[li[lu], ri[ru]] += 1
        assert (M.sum(0) == 1).all() and (M.sum(1) == 1).all()
        perms[n] = M
    for l, r in P.inner_fixed:  # CN 1 nodes: only one possible link
        assert (int(l), int(r)) in set(link[int(node[l])])
    return perms

def perms_to_walks(problem, perms):
    """walks (list of signed node id lists) and circular flags from hard permutations. Each walk is
    reported once, starting at an open end for linear walks."""
    P = problem
    side = P.side_of_uid; node = P.node_of_uid
    outer = {}
    for n1, n2 in P.ext:
        outer[int(n1)] = int(n2); outer[int(n2)] = int(n1)
    inner = {}
    for l, r in P.inner_fixed:
        inner[int(l)] = int(r); inner[int(r)] = int(l)
    for n, (lefts, rights) in P.var_nodes.items():
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
    for u in range(P.nu):
        if u not in outer and u not in visited:
            w, c = trace(u); walks.append(w); circ.append(c)
    for u in range(P.nu):
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
    a = sorted(canonical(w, c) for w, c in zip(w1, c1))
    b = sorted(canonical(w, c) for w, c in zip(w2, c2))
    return a == b
