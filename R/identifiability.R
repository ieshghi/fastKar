# Identifiability analysis utilities
#
# Local identifiability (manuscript Fig. 4): for one event, sample a truth
# karyotype and nearby candidates (local swaps at footprint nodes), simulate
# noisy assay observations of the truth, and record P(MLE = truth) per assay
# (run_graph_local_identifiability()). identifiability_context() sets up the
# buffered footprint, Hi-C resolution and edit-distance settings it uses.
#
# The earlier whole-graph "edit-distance ball" analysis (run_identifiability(),
# score_identifiability() and helpers) was removed on 2026-10-23 and archived in
# ~/projects/karyotyping/archive/fastKar/.

identifiability_context <- function(graph,
                                    footprint = NULL,
                                    graph_mode = c('footprint', 'full'),
                                    wholegenome = NULL,
                                    hic_binsperft = 100,
                                    footprint_buffer_fraction = 1/5,
                                    footprint_buffer_min = 1e4,
                                    edit_mode = c('fragment', 'footprint', 'node')) {
    graph_mode = match.arg(graph_mode)
    edit_mode = match.arg(edit_mode)
    if (is.null(footprint)) {
        footprint = graph$footprint
    }
    if (!is.null(wholegenome)) {
        footprint = footprint %&% wholegenome
    }
    ft0 = streduce(footprint)
    if (!is.numeric(footprint_buffer_fraction) || length(footprint_buffer_fraction) != 1 ||
        !is.finite(footprint_buffer_fraction) || footprint_buffer_fraction < 0) {
        stop('footprint_buffer_fraction must be a non-negative finite scalar')
    }
    if (!is.numeric(footprint_buffer_min) || length(footprint_buffer_min) != 1 ||
        !is.finite(footprint_buffer_min) || footprint_buffer_min < 0) {
        stop('footprint_buffer_min must be a non-negative finite scalar')
    }
    buffer = max(sum(width(ft0)) * footprint_buffer_fraction, footprint_buffer_min)
    ft = streduce(ft0 + buffer)
    target_width = sum(width(ft))
    hic_resolution = 1e3 * ceiling(target_width / (1e3 * hic_binsperft))
    graph_cut = loosefix(graph$copy$disjoin(ft) %&% ft)
    if (graph_mode == 'footprint') {
        analysis_graph = graph_cut
    } else {
        analysis_graph = graph
    }
    use.nodes = unique((analysis_graph$gr %&% ft)$node.id)
    ## How candidates are compared to the truth (edit_dist_cpp on use.nodes):
    ## 'node': whole walks, each node weighted by its full width (original).
    ## 'footprint': whole walks, each node weighted by its bp inside the
    ##   buffered footprint, the region the assays observe.
    ## 'fragment': as 'footprint', but walks are cut wherever they leave the
    ##   buffered footprint and compared as in-region fragments, so karyotypes
    ##   that differ only in how in-region pieces are linked through sequence
    ##   outside it (e.g. which copies share a chromosome) count as identical.
    edit_widths = if (edit_mode != 'node') .identifiability_edit_widths(analysis_graph, ft) else NULL
    split_ends = if (edit_mode == 'fragment') .identifiability_split_ends(analysis_graph, ft) else NULL
    list(
        graph = analysis_graph,
        target_region = ft,
        pix.size = hic_resolution,
        use.nodes = use.nodes,
        edit_mode = edit_mode,
        edit_widths = edit_widths,
        split_ends = split_ends,
        event_size = .identifiability_event_size(analysis_graph, use.nodes, edit_widths),
        footprint0 = ft0,
        target_width = target_width,
        footprint_buffer = buffer,
        footprint_buffer_fraction = footprint_buffer_fraction,
        footprint_buffer_min = footprint_buffer_min
    )
}

## bp of each node inside region `ft`, named by node.id (nodes outside omitted).
.identifiability_edit_widths <- function(graph, ft) {
    ngr = graph$nodes$gr
    ov = gr.findoverlaps(ngr, ft)
    w = data.table(node.id = ngr$node.id[ov$query.id], w = as.numeric(width(ov)))[
      , .(w = sum(w)), by = node.id]
    setNames(w$w, w$node.id)
}

## Whether each node's left and right end lies in region `ft`, as logicals
## named by node.id (nodes not overlapping ft omitted).
.identifiability_split_ends <- function(graph, ft) {
    ngr = graph$nodes$gr
    ngr = ngr[ngr %^% ft]
    left = gr.start(ngr, ignore.strand = TRUE) %^% ft
    right = gr.end(ngr, ignore.strand = TRUE) %^% ft
    list(left = setNames(left, ngr$node.id), right = setNames(right, ngr$node.id))
}

## Copy-number-weighted size (bp) of the event nodes, sum(width * cn), using the
## edit-distance widths when given. This is the edit distance from a karyotype
## to one with the event deleted, so it puts edit distances on a per-event scale.
.identifiability_event_size <- function(graph, use.nodes, edit_widths = NULL) {
    ## Row filter computed outside `[`: an `%in%` inside i can be optimised
    ## through a stale data.table index on graphs read back from .rds and
    ## silently return no rows.
    nd = graph$nodes$dt
    keep = nd$node.id %in% use.nodes
    nd = nd[keep, .(node.id, width = as.numeric(width), cn)]
    if (!is.null(edit_widths)) {
        nd[, width := edit_widths[as.character(node.id)]]
        nd[is.na(width), width := 0]
    }
    sum(nd$width * nd$cn, na.rm = TRUE)
}

## mclapply() turns a worker that dies (e.g. OOM-killed) into a NULL element
## and a worker error into a try-error, with only a warning; unlist() then
## silently drops or stringifies them. Fail loudly instead.
.mclapply_strict <- function(X, FUN, mc.cores = 1) {
    out = mclapply(X, FUN, mc.cores = mc.cores)
    bad = vapply(out, function(x) is.null(x) || inherits(x, 'try-error'), logical(1))
    if (length(out) != length(X) || any(bad)) {
        errs = out[vapply(out, inherits, logical(1), 'try-error')]
        stop(sum(bad), ' of ', length(X), ' parallel tasks failed',
             if (length(errs)) paste0(': ', conditionMessage(attr(errs[[1]], 'condition')))
             else ' (worker killed, likely out of memory)')
    }
    out
}

## Progress line for batch logs: time, step and the parent's resident memory.
.identifiability_log <- function(label, ...) {
    rss = tryCatch(as.numeric(gsub('[^0-9]', '', grep('^VmRSS', readLines('/proc/self/status'), value = TRUE))) / 2^20,
                   error = function(e) NA_real_)
    message(format(Sys.time(), '%H:%M:%S'), ' [', label, '] ', paste0(...), sprintf(' (RSS %.1f GB)', rss))
}

## Long-read word table restricted to reads touching the event (non-empty
## words). Walks that contain no event node only produce the empty word, so they
## are skipped before reads_fromwalk(); in full-graph mode that is most of the
## genome.
.identifiability_event_reads <- function(gw, readL, use.nodes = NULL) {
    walks = list(graph = gw$graph, snode.id = gw$snode.id, circular = gw$circular)
    if (!is.null(use.nodes)) {
        hit = vapply(walks$snode.id, function(w) any(abs(w) %in% use.nodes), logical(1))
        if (!any(hit)) {
            return(data.table(words = character(), nums = numeric()))
        }
        walks$snode.id = walks$snode.id[hit]
        walks$circular = walks$circular[hit]
    }
    reads_fromwalk(walks, readL, use.nodes = use.nodes)[words != '']
}

## -----------------------------------------------------------------------------
## Local identifiability: P(MLE = truth) for one truth karyotype among up to
## n_neighbors other sampled karyotypes, estimated from n_obs noisy
## observations of the truth per assay. Candidate expectations (Hi-C maps,
## long-read word tables) are computed once and reused for every observation.
## -----------------------------------------------------------------------------

## Footprint signature of a walk set (anything with $snode.id and $circular):
## the in-footprint pieces edit_dist_cpp() compares (fragments for
## edit_mode = 'fragment', else walks stripped to footprint nodes), each in a
## canonical orientation, sorted and hashed. Equal signatures imply edit
## distance 0, so candidates sharing one would tie on every assay.
.local_footprint_signature <- function(walk, ctx) {
    ids = ctx$graph$nodes$dt$node.id
    max_id = max(abs(ids))
    penalty = numeric(max_id)
    w = if (!is.null(ctx$edit_widths)) ctx$edit_widths[as.character(ids)] else ctx$graph$nodes$dt$width
    w[is.na(w) | !(ids %in% ctx$use.nodes)] = 0
    penalty[abs(ids)] = -w
    s = sort_snodes(walk$snode.id, walk$circular)
    if (!is.null(ctx$split_ends)) {
        left.in = right.in = logical(max_id)
        left.in[as.integer(names(ctx$split_ends$left))] = ctx$split_ends$left
        right.in[as.integer(names(ctx$split_ends$right))] = ctx$split_ends$right
        pieces = split_walk_fragments(s$nodelist, s$circ, penalty, left.in, right.in)
    } else {
        pieces = lapply(s$nodelist, function(x) x[penalty[abs(x)] != 0])
    }
    pieces = vapply(pieces[lengths(pieces) > 0], function(v) {
        min(paste(v, collapse = ','), paste(-rev(v), collapse = ','))
    }, character(1))
    digest::digest(sort(pieces), algo = 'xxhash64')
}

## A random truth plus up to n_neighbors karyotypes near it. The truth is a
## uniformly random wiring (as in sample.gwalks()); n_walkers independent
## local.sampling() chains each make n_swaps swaps from it, only at nodes inside
## the footprint (use.nodes). Every karyotype visited is a candidate; those with
## edit distance 0 to the truth are dropped, as are all but one of each set
## identical inside the footprint (same .local_footprint_signature()), and a
## random n_neighbors of the rest are kept. Events with few swappable nodes give
## fewer candidates: then every distinct local karyotype is compared.
## kars[[1]] is the truth; swaps is the fewest swaps at which each was reached.
sample_local_candidates <- function(ctx, n_neighbors = 99, n_walkers = 200, n_swaps = 3,
                                    mc.cores = 1) {
    wiring = gg.to.wiring(ctx$graph)
    w = wiring$internal.edges
    w[, right := if (.N > 1) sample(right) else right, by = n]
    truth_raw = traverse_graph_cpp(w, wiring$loose.ends)
    as_gw = function(x) gW(graph = ctx$graph, snode.id = x$snode.id, circular = x$circular)
    truth = as_gw(truth_raw)
    frozen = setdiff(unique(w$n), ctx$use.nodes)
    if (!nrow(w[cn > 1 & !(n %in% frozen)])) {
        return(list(kars = list(truth), edit = 0, swaps = 0L))
    }
    walks = .mclapply_strict(seq_len(n_walkers), function(i) {
        local.sampling(ctx$graph, nsteps = n_swaps + 1, nwalk = 1, frozen.nodes = frozen,
                       starter_edges = data.table::copy(w), return.gw = FALSE)[[1]][-1]
    }, mc.cores = mc.cores)
    walks = unlist(walks, recursive = FALSE)
    swaps = rep(seq_len(n_swaps), n_walkers)
    hash = vapply(walks, function(x) x$hash, character(1))
    ## fewest swaps first, so dedup keeps the closest visit of each karyotype
    ord = order(swaps)
    ord = ord[!duplicated(hash[ord]) & hash[ord] != truth$hash]
    walks = walks[ord]; swaps = swaps[ord]
    stats = .mclapply_strict(walks, function(x) {
        list(edit = edit_dist_cpp(truth_raw, x, graph = ctx$graph, use.nodes = ctx$use.nodes,
                                  node.widths = ctx$edit_widths, split.ends = ctx$split_ends),
             sig = .local_footprint_signature(x, ctx))
    }, mc.cores = mc.cores)
    edit = vapply(stats, `[[`, numeric(1), 'edit')
    sig = vapply(stats, `[[`, character(1), 'sig')
    keep = which(edit > 0 & !duplicated(sig))
    keep = keep[sample.int(length(keep), min(length(keep), n_neighbors))]
    list(kars = c(list(truth), lapply(walks[keep], as_gw)),
         edit = c(0, edit[keep]), swaps = c(0L, swaps[keep]))
}

## Hi-C NLL of n_obs noisy maps of kars[[1]] (the truth) under every candidate:
## an n_obs x length(kars) matrix. Candidates are simulated once, scaled to
## each map's depth by the first candidate's total, and scored by the negative
## binomial likelihood summed over pixels, vectorised over maps. Each pixel is weighted by its area
## width_i * width_j / pix.size^2 (as compdats() does): the tiles are cut at node boundaries, and a sliver
## tile should not count as much as a full pixel.
local_hic_nll_matrix <- function(kars, ctx, n_obs = 100, depth = 10, theta = 2) {
    true_expected = forward_simulate(kars[[1]], target_region = ctx$target_region,
                                     pix.size = ctx$pix.size, depth = depth)
    obs = make_noisydat(true_expected, nsamp = n_obs, theta = theta)
    gr = true_expected$gr
    template = make_template_dat(gr)[, .(i, j, id, widthprod)]
    setkey(template, id)
    weight = template$widthprod / ctx$pix.size^2
    X = round(vapply(obs, function(d) d[template, on = 'id']$value, numeric(nrow(template))))
    X = matrix(X, nrow = nrow(template))

    prepped = prep_for_sim(kars[[1]], ctx$target_region, ctx$pix.size)
    Y = vapply(kars, function(gw) {
        s = simulate_walks(gw, prepped$tiled.target, prepped$widthdt,
                           if.sum = TRUE, depth = depth / 2, gm.out = FALSE)$dat
        s = s[, .(i, j, value)][template, on = c('i', 'j')]
        s[is.na(value), value := 0]
        s[order(id)]$value
    }, numeric(nrow(template)))
    Y = matrix(Y, nrow = nrow(template))
    ratio = colSums(X) / sum(Y[, 1])

    out = matrix(NA_real_, n_obs, length(kars))
    for (k in seq_along(kars)) {
        nz = Y[, k] > 0
        mu = outer(Y[nz, k], ratio)
        lp = matrix(dnbinom(X[nz, , drop = FALSE], mu = mu, size = theta, log = TRUE), nrow = sum(nz))
        out[, k] = -colSums(lp * weight[nz], na.rm = TRUE)
    }
    out
}

## Long-read NLL of n_obs read sets drawn from kars[[1]] (the truth) under
## every candidate: reads restricted to those touching the event, words
## hashed, background smoothing over the shared vocabulary. n_obs x length(kars) matrix.
local_longread_nll_matrix <- function(kars, ctx, readL, n_obs = 100, depth = 10,
                                      background = 1e-10) {
    hash_words = digest::getVDigest(algo = 'xxhash64')
    tabs = lapply(kars, function(gw) {
        rt = .identifiability_event_reads(gw, readL, ctx$use.nodes)
        if (!nrow(rt)) return(numeric())
        rt = data.table(h = hash_words(rt$words, serialize = FALSE), n = rt$nums)[, .(n = sum(n)), by = h]
        setNames(rt$n, rt$h)
    })
    out = matrix(Inf, n_obs, length(kars))
    truth = tabs[[1]]
    if (!length(truth)) return(out)
    V = length(unique(unlist(lapply(tabs, names), use.names = FALSE)))
    context = ctx$target_width + 2 * readL
    n_reads = max(1L, as.integer(round(depth / (readL / context))))
    for (m in seq_len(n_obs)) {
        oc = table(sample(names(truth), size = n_reads, prob = truth, replace = TRUE))
        for (k in seq_along(tabs)) {
            tk = tabs[[k]]
            if (!length(tk)) next
            cnt = tk[names(oc)]
            cnt[is.na(cnt)] = 0
            out[m, k] = -sum(as.numeric(oc) * log((cnt + background) / (sum(tk) + background * V)))
        }
    }
    attr(out, 'n_reads') = n_reads
    out
}

## Per observation: the truth's (column 1) NLL rank, the number of candidates
## tied for the best NLL, and p_pick = P(MLE = truth) with ties broken
## uniformly (1/k if the truth ties k-way for best, else 0).
.local_pick_summary <- function(nll, tol = 1e-8) {
    rbindlist(lapply(seq_len(nrow(nll)), function(m) {
        x = nll[m, ]
        best = x <= min(x) + tol * max(abs(min(x)), 1)
        data.table(obs_id = m,
                   truth_rank = 1L + sum(x[-1] < x[1] - tol * max(abs(x[1]), 1)),
                   n_best = sum(best),
                   p_pick = best[1] / sum(best))
    }))
}

run_graph_local_identifiability <- function(graph,
                                            footprint = NULL,
                                            event = NULL,
                                            seed = NULL,
                                            n_neighbors = 99,
                                            n_walkers = 200,
                                            n_swaps = 3,
                                            n_obs = 100,
                                            graph_mode = 'full',
                                            modalities = c('hic', 'longread'),
                                            longread_lengths = c(1e5, 1e7),
                                            depth = 10,
                                            theta = 2,
                                            hic_binsperft = 300,
                                            footprint_buffer_fraction = 1/5,
                                            footprint_buffer_min = 1e4,
                                            edit_mode = 'fragment',
                                            longread_background = 1e-10,
                                            wholegenome = NULL,
                                            mc.cores = 1) {
    modalities = match.arg(modalities, choices = c('hic', 'longread'), several.ok = TRUE)
    ctx = identifiability_context(graph = graph, footprint = footprint, graph_mode = graph_mode,
                                  wholegenome = wholegenome, hic_binsperft = hic_binsperft,
                                  footprint_buffer_fraction = footprint_buffer_fraction,
                                  footprint_buffer_min = footprint_buffer_min, edit_mode = edit_mode)
    label = if (!is.null(event$event_id)) as.character(event$event_id) else 'event'
    .identifiability_log(label, 'start: ', length(ctx$graph$nodes), ' nodes, target ',
                         ctx$target_width, ' bp, pix ', ctx$pix.size)
    if (!is.null(seed)) set.seed(seed)
    cand = sample_local_candidates(ctx, n_neighbors = n_neighbors, n_walkers = n_walkers,
                                   n_swaps = n_swaps, mc.cores = mc.cores)
    .identifiability_log(label, length(cand$kars) - 1L, ' candidates besides the truth')

    nll = list()
    if ('hic' %in% modalities && length(cand$kars) > 1) {
        nll[['hic']] = local_hic_nll_matrix(cand$kars, ctx, n_obs = n_obs, depth = depth, theta = theta)
        .identifiability_log(label, 'hic scored')
    }
    if ('longread' %in% modalities && length(cand$kars) > 1) {
        for (readL in longread_lengths) {
            nll[[paste0('longread_', readL)]] = local_longread_nll_matrix(
                cand$kars, ctx, readL, n_obs = n_obs, depth = depth, background = longread_background)
            .identifiability_log(label, 'long reads ', readL, ' scored')
        }
    }

    picks = rbindlist(lapply(names(nll), function(a) {
        dt = .local_pick_summary(nll[[a]])
        dt[, `:=`(modality = sub('_.*', '', a),
                  readL = if (a == 'hic') NA_real_ else as.numeric(sub('.*_', '', a)),
                  n_reads = if (a == 'hic') NA_integer_ else attr(nll[[a]], 'n_reads'))]
        dt
    }), fill = TRUE)
    ## Only one karyotype (inside the footprint) is possible, so every assay
    ## finds the truth: one row per assay with p_pick = 1.
    unique_walk = length(cand$kars) == 1
    if (unique_walk) {
        picks = data.table(obs_id = NA_integer_, truth_rank = 1L, n_best = 1L, p_pick = 1,
                           modality = c(if ('hic' %in% modalities) 'hic',
                                        if ('longread' %in% modalities) rep('longread', length(longread_lengths))),
                           readL = c(if ('hic' %in% modalities) NA_real_,
                                     if ('longread' %in% modalities) longread_lengths),
                           n_reads = NA_integer_)
    }
    meta = data.table(n_cand = length(cand$kars), unique_walk = unique_walk,
                      event_size = ctx$event_size, target_width = ctx$target_width)
    if (!is.null(event)) meta = cbind(as.data.table(event), meta)
    if (nrow(picks)) picks = cbind(meta[rep(1L, nrow(picks))], picks)
    .identifiability_log(label, 'done')
    list(picks = picks, meta = meta, edit = cand$edit, swaps = cand$swaps, nll = nll)
}
