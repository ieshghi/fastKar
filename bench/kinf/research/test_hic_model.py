import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../inst/python/kinf'))
import jax; jax.config.update('jax_enable_x64', True)
import numpy as np, jax.numpy as jnp, hic_model as hm
K = hm.get_constants(); chk = K['raw']['check']
d = np.array(chk['dist']); got = np.array(hm.distance_density(d))
print('density rel err', np.max(np.abs(got/np.array(chk['dens'])-1)))
assert np.max(np.abs(got/np.array(chk['dens'])-1)) < 1e-10
wp = np.array(chk['wp']); g = np.array(hm.diag_contacts(wp)); print('diag rel err', np.max(np.abs(g/np.array(chk['diag'])-1)))
assert np.max(np.abs(g/np.array(chk['diag'])-1)) < 1e-10
# tail and small-d branch vs R (checked via R later in validate); gradients finite
f = lambda x: hm.pair_contacts(x, 1e9)
for x in [0.5, 1., 5e3, 1e4, 1e6, 1e8, 3e8, 1e12]:
    gr = jax.grad(f)(jnp.float64(x)); assert np.isfinite(gr), x
print('grads finite; OK')
