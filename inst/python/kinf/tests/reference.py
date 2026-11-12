"""Slow reference implementations for the tests (not used by the engine): the expected map from all-pairs end distances, by
walk traversal (walkdist) and by Floyd-Warshall on the one-step cost matrix (fw_costmat), then refined to tile copies."""
import numpy as np
import jax
import jax.numpy as jnp
import problem as pr

def walkdist(g, b, bigdist=1e12):
    """end x end distances: along the shared walk (shorter way round on a cycle), bigdist between walks"""
    coord, wid, L, iscyc = pr.walkinfo(g, b)
    d = jnp.abs(coord[:, None] - coord[None, :])
    d = jnp.where(iscyc[:, None], jnp.minimum(d, L[:, None] - d), d)
    D = jnp.where(wid[:, None] == wid[None, :], d, bigdist)
    return D.at[jnp.diag_indices(2 * g.nu)].set(0.0)

def fw_costmat(g, perms, bigdist=1e12):
    """(2nu x 2nu) one-step cost over ends: half width inside a half; junction and present permutation links 0; absent bigdist"""
    N = 2 * g.nu; side = g.side_of_uid
    base = np.full((N, N), bigdist)
    np.fill_diagonal(base, 0.0)
    e = np.arange(g.nu)
    base[2 * e, 2 * e + 1] = g.width_uid
    base[2 * e + 1, 2 * e] = g.width_uid
    for n1, n2 in g.ext:
        a = pr.outer_end(side[n1], n1); c = pr.outer_end(side[n2], n2)
        base[a, c] = base[c, a] = 0.0
    for l, r in g.inner_fixed:
        a = pr.inner_end(side[l], l); c = pr.inner_end(side[r], r)
        base[a, c] = base[c, a] = 0.0
    for n in g.vnodes:
        l, r = g.var_nodes[n]
        rr, cc = np.meshgrid(pr.inner_end(0, l), pr.inner_end(1, r), indexing='ij')
        v = np.where(np.asarray(perms[n]).ravel() > 0.5, 0.0, bigdist)
        base[rr.ravel(), cc.ravel()] = v; base[cc.ravel(), rr.ravel()] = v
    return jnp.asarray(base)

def floyd_warshall(D):
    """all-pairs shortest path, O(N^3)"""
    N = D.shape[0]
    idx = jnp.arange(N)
    def body(k, D):
        upd = jnp.minimum(D, D[:, k][:, None] + D[k, :][None, :])
        return jnp.where((idx[:, None] == k) | (idx[None, :] == k), D, upd)
    return jax.lax.fori_loop(0, N, body, D).at[jnp.diag_indices(N)].set(0.0)

def map_from_enddist(hic, Dend, dcut=1e10):
    """end-end distances -> tile-copy distances (min over the 4 end pairings of offset + end distance; |offset difference|
    within one half copy) -> expected tile map"""
    uid = hic.tc_uid
    Lu = jnp.asarray(2 * uid); Ru = jnp.asarray(2 * uid + 1)
    oL = jnp.asarray(hic.offL); oR = jnp.asarray(hic.offR)
    LL = Dend[Lu][:, Lu] + oL[:, None] + oL[None, :]
    RR = Dend[Ru][:, Ru] + oR[:, None] + oR[None, :]
    LR = Dend[Lu][:, Ru] + oL[:, None] + oR[None, :]
    RL = Dend[Ru][:, Lu] + oR[:, None] + oL[None, :]
    d = jnp.where(jnp.asarray(uid[:, None] == uid[None, :]), jnp.abs(oL[:, None] - oL[None, :]),
                  jnp.min(jnp.stack([LL, RR, LR, RL]), axis=0))
    d = jnp.where(jnp.eye(len(uid), dtype=bool), 0.0, d)
    return pr.contact_model(d, jnp.asarray(hic.tc_w), jnp.asarray(hic.tc_tile), jnp.asarray(hic.tile_w),
                            jnp.asarray(hic.tile_cn), hic.depth, dcut)
