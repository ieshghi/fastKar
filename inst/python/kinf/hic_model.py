"""Port of the fastKar (R) Hi-C forward model: contact densities and per-tile-copy contact assembly.

Everything here mirrors ~/git/fastKar/R/forward.R (distance_density, simulate_map) and
src/forward_kernel.cpp (forward_walks_cpp, forward_assemble_cpp). Constants come from
data/fastkar_model.json (written by scripts/export_model.R).
"""
import json
import os
import numpy as np
import jax
import jax.numpy as jnp

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
    """contact density per bp^2 at distance d > 0 (bp). Exact port: spline above exp(lmin), power law below,
    1/d tail above 1e8 (the tail is applied here so the function matches what simulate_map() computes)."""
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
    """contacts in the diagonal pixel of a tile: exp(diag_spline(log(widthprod))) (linear in log-log)"""
    K = K or get_constants()
    return jnp.exp(K['diag_a'] + K['diag_b'] * jnp.log(widthprod))

def pair_contacts(d, widthprod, dcut=1e10, dfloor=1e-3, K=None):
    """R's f(a,b) for two distinct rows at distance d: widthprod * distance_density(d).
    Pairs on different components of the graph sit at ~bigdist; the smooth gate 1/(1+(d/dcut)^4) removes them
    (R gives them no contact). The gate is >1-1e-12 for d <= 1e7."""
    d = jnp.maximum(d, dfloor)
    gate = 1.0 / (1.0 + (d / dcut) ** 4)
    return widthprod * distance_density(d, K) * gate

def contact_matrix(D, w, dcut=1e10, K=None):
    """per tile-copy pair contacts M[u,v] (symmetric, no mult/depth). D: (ntc,ntc) distances in bp, w: (ntc,) widths.
    Off-diagonal = widthprod*density(D); diagonal (u == v) = diag_contacts(w^2)."""
    wp = w[:, None] * w[None, :]
    eye = jnp.eye(w.shape[0], dtype=bool)
    Ds = jnp.where(eye, 1.0, D)  # keep the unused branch finite
    off = pair_contacts(Ds, wp, dcut, K=K)
    return jnp.where(eye, diag_contacts(wp, K), off)

def contact_model(D, w, tile_of_copy, w_tile, cn_tile, depth, dcut=1e10, K=None):
    """expected counts per reference tile pair (full symmetric ntile x ntile; the R map is its upper triangle).
    value(p,q) = depth/2 * [ sum_{u in p, v in q} M[u,v] + w_p w_q rho cn_p cn_q ].
    Summing M over ordered copy pairs reproduces R's rules: a==b once, a!=b in the same tile p==q twice
    (mult 2), p!=q once."""
    K = K or get_constants()
    M = contact_matrix(D, w, dcut, K)
    ntile = w_tile.shape[0]
    # sum copy pairs into tile pairs (= T.T @ M @ T with T one-hot, without the O(ntile ncopy^2) float64 matmuls)
    seg = lambda X: jax.ops.segment_sum(X, tile_of_copy, num_segments=ntile)
    R = seg(seg(M).T).T
    bg = w_tile[:, None] * w_tile[None, :] * K['rho'] * (cn_tile[:, None] * cn_tile[None, :])
    return 0.5 * depth * (R + bg)
