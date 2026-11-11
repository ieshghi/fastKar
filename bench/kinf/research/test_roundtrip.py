import sys, os, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../inst/python/kinf'))
import model, walks
V2 = '/gpfs/commons/groups/imielinski_lab/projects/karyotyping/karyotype-inference/temp_data/v2/'
prob = model.load_problem(V2 + 'sim_t1_d4')
nok = 0
for t in range(1, 34):
    sn, circ = walks.read_gwalks(V2 + f'sim_t{t}_d4/truth_gwalks.json')
    perms = walks.walks_to_perms(prob, sn, circ)
    w2, c2 = walks.perms_to_walks(prob, perms)
    assert walks.same_karyotype(sn, circ, w2, c2), (t, sn, w2)
    nok += 1
print('round trip OK for', nok, 'truths')

# random hard permutations -> walks -> perms -> walks must give the same karyotype
import numpy as np
rng = np.random.default_rng(0)
for k in range(200):
    perms = {n: np.eye(len(l))[rng.permutation(len(l))] for n, (l, r) in prob.var_nodes.items()}
    w1, c1 = walks.perms_to_walks(prob, perms)
    p2 = walks.walks_to_perms(prob, w1, c1)
    w2, c2 = walks.perms_to_walks(prob, p2)
    assert walks.same_karyotype(w1, c1, w2, c2)
    # same karyotype => same hard map
print('random perm round trips OK (200), e.g. n walks', len(w1), 'circular', sum(c1))
