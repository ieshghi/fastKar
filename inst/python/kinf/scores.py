"""Plain scores f(kar) -> scalar (jittable; static arrays precomputed outside JAX) and energy helpers for optimize.py."""
import numpy as np
import jax
import jax.numpy as jnp


def walk_sum(kar, per_end_values):
    """sum of per-end values over the ends of each walk: array of length 2*nu indexed by walk id (tail state; unused ids are 0)"""
    v = jnp.asarray(per_end_values)
    return jax.ops.segment_sum(v, kar.wid, num_segments=len(v))

def n_walks(g):
    N = 2 * g.nu
    return lambda kar: jnp.zeros(N).at[kar.wid].set(1.0).sum()

def n_circular(g):
    N = 2 * g.nu
    return lambda kar: jnp.zeros(N).at[kar.wid].max(kar.iscyc.astype(float)).sum()

def amp_entropy(g, amp_nodes):
    """entropy of the distribution of amplicon node visits over walks (a visit = left + right half copy, 1/2 each)"""
    N = 2 * g.nu
    w = np.zeros(N); w[0::2] = np.where(np.isin(g.node_of_uid, np.asarray(amp_nodes)), 0.5, 0.0)  # one end per half copy
    w = jnp.asarray(w)
    def f(kar):
        c = walk_sum(kar, w)
        p = c / c.sum()
        return -jnp.sum(jnp.where(p > 0, p * jnp.log(jnp.where(p > 0, p, 1.0)), 0.0))
    return f

def schedule(ts, vs):
    """piecewise-linear t -> value, constant outside [ts[0], ts[-1]]"""
    ts, vs = jnp.asarray(ts, float), jnp.asarray(vs, float)
    return lambda t: jnp.interp(t, ts, vs)

def target(k, sched):
    """energy(values, t) = |values[k] - sched(t)|"""
    return lambda values, t: jnp.abs(values[k] - sched(t))
