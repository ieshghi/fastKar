"""Generic optimizer of kinf: MCMC (Metropolis, optional parallel tempering / annealing) over wirings b of a genome graph.

A score is a Term(init, value, update) or a plain jittable f(kar) -> scalar. Every score sees a Kar (b + walk structure,
built once per proposal). Time dependence lives only in energy(values, t). See optimize() below.
"""
import time
from typing import NamedTuple, Callable, Any
import numpy as np
import jax
import jax.numpy as jnp
import problem as pr


class Kar(NamedTuple):
    b: Any; coord: Any; wid: Any; L: Any; iscyc: Any


class Term(NamedTuple):
    init: Callable    # init(kar) -> aux
    value: Callable   # value(aux) -> scalar
    update: Callable  # update(aux, kar_old, kar_new, i, j) -> aux (i, j = swapped rows of b)


def kar(g, b):
    return Kar(b, *pr.walkinfo(g, b))


def as_term(f):
    if isinstance(f, Term):
        return f
    return Term(init=f, value=lambda aux: aux, update=lambda aux, ko, kn, i, j: f(kn))


def sum_energy(values, t):
    return sum(values)


def optimize(g, terms, energy=None, T=1.0, chains=16, steps=20000, block=500, seed=0, init=None, init_random=0,
             t0=None, history_every=0, topk=1, vmap=False, log=print):
    """Metropolis search over wirings of g minimising energy(values, t); values = tuple of the terms' values.
    terms: Term or plain f(kar) -> scalar (jittable), a list of them. energy: (values, t) -> scalar, jittable, t = global
    step (traced); default sum(values). Each step the energy of the current and the proposed state are both computed at
    the current t. Accept if dE < 0 or log u < -dE / T. At every block boundary every term is re-initialised from scratch.
    T: scalar (all chains); 1-D array of length chains (ladder = parallel tempering: adjacent-rung swaps after every block,
    even/odd alternating); or a jittable callable t -> (chains,) array (annealing; swaps use T at the block end). T = 0 is greedy.
    init: None (random), list/array of b (round robin; last init_random chains random) or a previous Result (continue; use a
    new seed; t0 defaults to its t_end). Per-chain loop over one jitted block function; vmap=True vmaps it over chains.
    Returns a dict (Result): b (chains, n), rung, E, values (nterms, chains), t_end, best (E, b, values), best_E/best_b/
    best_vals per chain, top (list of (E, b) sorted: best + the 2*topk lowest-E states per chain, see
    problem.distinct_karyotypes), history (dict of arrays: step, chain, T, E, v0..), trace (dict: steps, best, cur (nblk, chains
    by rung)), summary."""
    C = chains
    terms = [as_term(f) for f in (terms if isinstance(terms, (list, tuple)) and not isinstance(terms, Term) else [terms])]
    nt = len(terms)
    energy = energy or sum_energy
    callable_T = callable(T)
    Tarr = None if callable_T else np.broadcast_to(np.asarray(T, float), (C,)).copy()
    ladder = callable_T or np.ndim(T) == 1
    prev = init if isinstance(init, dict) else None
    if t0 is None:
        t0 = int(prev['t_end']) if prev else 0
    Tat = lambda t: np.asarray(T(t), float) if callable_T else Tarr  # host: T per rung at step t
    Tf = (lambda t, r: jnp.asarray(T(t))[r]) if callable_T else (lambda t, r: jnp.asarray(Tarr)[r])

    start = jax.jit(lambda b: tuple(tm.init(kar(g, b)) for tm in terms))
    vals_of = lambda aux: tuple(tm.value(a) for tm, a in zip(terms, aux))
    if len(g.vnodes) == 0:
        b = jnp.zeros(0, int); v = vals_of(start(b)); E = float(energy(v, t0))
        summ = dict(chains=C, steps=0, block=block, best_E=E, n_var_nodes=0, n_moves=0, wall_s=0.0)
        return dict(b=np.zeros((1, 0), int), rung=np.zeros(1, int), E=np.array([E]), values=np.array([[float(x)] for x in v]),
                    t_end=t0, best=dict(E=E, b=np.zeros(0, int), values=[float(x) for x in v]), top=[(E, np.zeros(0, int))],
                    history={}, trace={}, summary=summ)

    nb = max(steps // block, 1)
    moves = jnp.asarray(g.moves); mnode = g.move_node; n_nodes = len(g.vnodes)
    off = jnp.asarray(np.searchsorted(mnode, np.arange(n_nodes))); nmv = jnp.asarray(np.bincount(mnode, minlength=n_nodes))
    K = 2 * topk if topk > 1 else 0

    def buf_insert(bn, bb, E, b):
        dup = jnp.any(jnp.abs(bn - E) <= 1e-9 * jnp.maximum(1.0, jnp.abs(E)))
        w = jnp.argmax(bn)
        ins = (E < bn[w]) & ~dup
        return bn.at[w].set(jnp.where(ins, E, bn[w])), bb.at[w].set(jnp.where(ins, b, bb[w]))

    def block_fn(S, key, rung, tstart):
        ko = kar(g, S['b'])
        def step(c, x):
            ko, aux, bE, bb, bv, acc, bufE, bufb = c
            k, tt = x
            k1, k2, k3 = jax.random.split(k, 3)
            nd = jax.random.randint(k1, (), 0, n_nodes)
            i, j = moves[off[nd] + jax.random.randint(k2, (), 0, nmv[nd])]
            b = ko.b
            kn = kar(g, b.at[i].set(b[j]).at[j].set(b[i]))
            aux2 = tuple(tm.update(a, ko, kn, i, j) for tm, a in zip(terms, aux))
            vo, vn = vals_of(aux), vals_of(aux2)
            Eo, En = energy(vo, tt), energy(vn, tt)
            dE = En - Eo
            ok = (dE < 0) | (jnp.log(jax.random.uniform(k3)) < -dE / jnp.maximum(Tf(tt, rung), 1e-300))
            sel = lambda u, v: jax.tree.map(lambda p, q: jnp.where(ok, p, q), u, v)
            ko, aux = sel(kn, ko), sel(aux2, aux)
            E, v = jnp.where(ok, En, Eo), sel(vn, vo)
            better = E < bE
            bE, bb = jnp.where(better, E, bE), jnp.where(better, ko.b, bb)
            bv = tuple(jnp.where(better, x, y) for x, y in zip(v, bv))
            if K:
                bufE, bufb = buf_insert(bufE, bufb, E, ko.b)
            return (ko, aux, bE, bb, bv, acc + ok, bufE, bufb), (E, jnp.stack(v))
        c0 = (ko, S["aux"], S["bE"], S["bb"], S["bv"], jnp.int32(0), S['bufE'], S['bufb'])
        ks = jax.random.split(key, block)
        (ko, aux, bE, bb, bv, acc, bufE, bufb), ys = jax.lax.scan(step, c0, (ks, tstart + jnp.arange(block)))
        S2 = dict(b=ko.b, aux=aux, bE=bE, bb=bb, bv=bv, bufE=bufE, bufb=bufb)
        return S2, acc, ys

    bj = jax.vmap(block_fn, in_axes=(0, 0, 0, None)) if vmap else None
    bj = jax.jit(bj if vmap else block_fn)

    @jax.jit
    def reinit(S, t):  # drift control: aux from scratch; returns state, energy, values
        aux = tuple(tm.init(kar(g, S['b'])) for tm in terms)
        v = vals_of(aux)
        return dict(S, aux=aux), energy(v, t), jnp.stack(v)

    # initial wirings
    key = jax.random.PRNGKey(seed)
    b0 = pr.random_b(g, jax.random.fold_in(key, 999), C)
    rung = np.arange(C)
    if prev is not None:
        b0 = jnp.asarray(prev['b'], b0.dtype); rung = np.asarray(prev['rung']).copy()
        assert len(b0) == C, 'previous Result has a different number of chains'
        log(f'continuing {C} chains from t = {t0}')
    elif init is not None:
        bi = [pr.check_b(g, b) for b in init]
        n_init = max(C - init_random, 0)
        b0 = jnp.stack([jnp.asarray(bi[c % len(bi)], b0.dtype) if c < n_init else b0[c] for c in range(C)])
        log(f'warm start: {min(n_init, C)} chains from {len(bi)} given wirings, {C - min(n_init, C)} random')
    nK = max(K, 1)
    def mk(b):
        aux = start(b); v = vals_of(aux); E = energy(v, t0)
        bufE, bufb = jnp.full(nK, jnp.inf), jnp.zeros((nK,) + b.shape, b.dtype)
        if K:
            bufE, bufb = buf_insert(bufE, bufb, E, b)
        return dict(b=b, aux=aux, bE=E, bb=b, bv=v, bufE=bufE, bufb=bufb)
    states = [mk(b0[c]) for c in range(C)]
    stack = lambda L: jax.tree.map(lambda *x: jnp.stack(x), *L)
    unstack = lambda S, n: [jax.tree.map(lambda x: x[c], S) for c in range(n)]
    E = np.array([float(energy(vals_of(s['aux']), t0)) for s in states])
    swaps = ladder and C > 1 and np.any(Tat(t0) > 0)

    hist = history_every > 0
    H = []
    def hrows(step, ch, Tc, Ec, vc):
        H.append(np.concatenate([[step, ch, Tc, Ec], vc]))
    if hist:
        for c in range(C):
            hrows(t0, c, Tat(t0)[rung[c]], E[c], [float(x) for x in vals_of(states[c]['aux'])])
    tr_steps, tr_best, tr_cur = [], [], []
    acc_tot = np.zeros(C); sw_try = np.zeros(max(C - 1, 1)); sw_acc = np.zeros(max(C - 1, 1))
    wall0 = time.time()
    for blk in range(nb):
        tb = t0 + blk * block
        ks = []
        for c in range(C):
            key, k = jax.random.split(key); ks.append(k)
        if vmap:
            S2, acc, ys = bj(stack(states), jnp.stack(ks), jnp.asarray(rung), tb)
            states = unstack(S2, C); acc = np.asarray(acc); ysE = np.asarray(ys[0]); ysV = np.asarray(ys[1])
            outs = [(acc[c], ysE[c], ysV[c]) for c in range(C)]
        else:
            outs = []
            for c in range(C):
                states[c], a, ys = bj(states[c], ks[c], int(rung[c]), tb)
                outs.append((float(a), np.asarray(ys[0]), np.asarray(ys[1])))
        for c in range(C):
            acc_tot[rung[c]] += float(outs[c][0])
            if hist:
                ix = np.arange(history_every - 1, block, history_every)
                for i in ix:
                    t = tb + i + 1
                    hrows(t, c, Tat(t)[rung[c]], outs[c][1][i], outs[c][2][i])
        te = t0 + (blk + 1) * block
        rs = [reinit(states[c], te) for c in range(C)]
        states = [r[0] for r in rs]; E = np.array([float(r[1]) for r in rs])
        if swaps:
            temps = Tat(te); at = np.argsort(rung)
            for r in range(blk % 2, C - 1, 2):
                ca, cb = at[r], at[r + 1]
                with np.errstate(all='ignore'):
                    d = (1 / temps[r] - 1 / temps[r + 1]) * (E[ca] - E[cb])
                sw_try[r] += 1
                if np.log(np.random.default_rng([seed, blk, r]).random()) < d:
                    rung[ca], rung[cb] = r + 1, r; sw_acc[r] += 1
        bg = min(float(s['bE']) for s in states); at = np.argsort(rung)
        tr_steps.append(te); tr_best.append(bg); tr_cur.append(E[at])
        if log and (blk % 10 == 0 or blk == nb - 1):
            log(f'block {blk + 1}/{nb} best {bg:.1f} rungs (min/mid/max T) {np.round(E[at])[[0, C // 2, C - 1]]} '
                f'{time.time() - wall0:.0f}s')
    wall = time.time() - wall0
    ibE = np.array([float(s['bE']) for s in states]); ib = int(np.argmin(ibE))
    best = dict(E=float(ibE[ib]), b=np.asarray(states[ib]['bb']), values=[float(x) for x in states[ib]['bv']])
    top = [(best['E'], best['b'])]
    if K:
        cand = sorted(((float(v), np.asarray(bb)) for s in states for v, bb in zip(s['bufE'], s['bufb']) if np.isfinite(v)),
                      key=lambda t: t[0])
        top = top + cand
    history = {}
    if hist:
        Hh = np.array(H); history = dict(step=Hh[:, 0].astype(int), chain=Hh[:, 1].astype(int), T=Hh[:, 2], E=Hh[:, 3],
                                         values=Hh[:, 4:])
    summ = dict(chains=C, steps=nb * block, block=block, seed=seed, best_E=best['E'], n_var_nodes=n_nodes, n_moves=len(g.moves),
                acc_rate_per_rung=list(acc_tot / (nb * block)), swap_rate=list(sw_acc / np.maximum(sw_try, 1)),
                wall_s=wall, ms_per_proposal=wall / (nb * block * C) * 1e3)
    return dict(b=np.stack([np.asarray(s['b']) for s in states]), rung=rung, E=E,
                values=np.array([[float(x) for x in vals_of(s['aux'])] for s in states]).T, t_end=t0 + nb * block,
                best=best, best_E=ibE, best_b=np.stack([np.asarray(s['bb']) for s in states]),
                best_vals=np.array([[float(x) for x in s['bv']] for s in states]).T, top=top, history=history,
                trace=dict(steps=np.array(tr_steps), best=np.array(tr_best), cur=np.array(tr_cur)), summary=summ)
