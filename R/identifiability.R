# Identifiability analysis utilities
#
# These functions implement the karyotype identifiability analysis: sample a
# truth karyotype, simulate noisy assay data, naively sample one candidate set
# (scoring Hi-C as in bestfit_search() when requested), rank that candidate set by each
# assay's likelihood, and record each candidate's edit distance to the truth.
# The candidate edit distances double as the truth-centred null distribution;
# scoring is done downstream from `ranked`.

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

## Event sizes for run_identifiability() results saved before `event_size` was
## added to `ranked`: one row per event, computed from the stored contexts.
identifiability_event_sizes <- function(result) {
    unique(rbindlist(lapply(result$raw_output, function(ro) {
        ctx = ro$context
        size = if (!is.null(ctx$event_size)) ctx$event_size
               else .identifiability_event_size(ctx$graph, ctx$use.nodes)
        data.table(event_id = ro$ranked$event_id[1], event_size = size)
    })), by = 'event_id')
}

.sample_identifiability_gwalks <- function(graph,
                                         n,
                                         mc.cores = 1,
                                         min_n = 1L,
                                         max_tries = 5L) {
    out = list()
    seen = character()
    for (try_i in seq_len(max_tries)) {
        need = max(n - length(out), min_n)
        gw = sample.gwalks(graph, need, verbose = FALSE, mc.cores = mc.cores)
        if (length(gw)) {
            h = vapply(gw, function(x) x$hash, character(1))
            keep = !(h %in% seen)
            out = c(out, gw[keep])
            seen = c(seen, h[keep])
        }
        if (length(out) >= min_n && length(out) >= n) break
    }
    if (length(out) > n) out = out[seq_len(n)]
    if (length(out) < min_n) {
        stop('Could not sample any valid karyotypes from graph after ', max_tries, ' tries')
    }
    out
}

.identifiability_add_metadata <- function(dt,
                                          event,
                                          ctx,
                                          truth,
                                          rep_id,
                                          nsample,
                                          graph_mode,
                                          modality,
                                          readL = NA_real_) {
    if (!is.null(event)) {
        for (nm in names(event)) {
            val = event[[nm]]
            if (length(val) == 1 && !(nm %in% names(dt))) {
                dt[, (nm) := val]
            }
        }
    }
    dt[, rep_id := rep_id]
    dt[, nsample := nsample]
    dt[, graph_mode := graph_mode]
    dt[, modality := modality]
    dt[, readL := readL]
    dt[, pix.size := ctx$pix.size]
    dt[, target_width := ctx$target_width]
    dt[, event_size := ctx$event_size]
    dt[, edit_mode := ctx$edit_mode]
    dt[, footprint_buffer := ctx$footprint_buffer]
    dt[, footprint_buffer_fraction := ctx$footprint_buffer_fraction]
    dt[, footprint_buffer_min := ctx$footprint_buffer_min]
    dt[, truth_hash := truth$true_gw$hash]
    dt
}

make_identifiability_truth <- function(ctx,
                                       depth = 10,
                                       theta = 2,
                                       simulate_hic = TRUE,
                                       mc.cores = 1) {
    true_gw = .sample_identifiability_gwalks(ctx$graph, 1, mc.cores = mc.cores, min_n = 1L)[[1]]
    if (!simulate_hic) {
        return(list(true_gw = true_gw, true_expected = NULL, true_hic = NULL))
    }
    true_expected = forward_simulate(
        true_gw,
        target_region = ctx$target_region,
        pix.size = ctx$pix.size,
        depth = depth
    )
    true_hic = make_noisymap(true_expected, nsamp = 1, theta = theta, mc.cores = 1)[[1]]
    list(true_gw = true_gw, true_expected = true_expected, true_hic = true_hic)
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

## Hi-C NLL of the truth's noisy map under each candidate. Equivalent to
## bestfit_search(topk = -1)$nlls, but the tiling (shared by all candidates on
## the same graph) is prepared once and each worker returns only a scalar, so
## no simulated maps are kept.
identifiability_hic_nlls <- function(kars, truth, ctx, depth = 10, theta = 2, mc.cores = 1) {
    prepped = prep_for_sim(kars[[1]], ctx$target_region, ctx$pix.size)
    obs = truth$true_hic
    obs_dat = rebin_matrix(obs, obs$gr)$dat
    template = make_template_dat(obs$gr)[, .(i, j, id)]
    area0 = median(width(obs$gr))^2
    sim_dat = function(gw) {
        s = simulate_walks(gw, prepped$tiled.target, prepped$widthdt,
                           if.sum = TRUE, depth = depth / 2, gm.out = FALSE)$dat
        s = merge.data.table(template, s[, .(i, j, value)], by = c('i', 'j'), all.x = TRUE)
        s[is.na(value), value := 0]
        s
    }
    ## As in bestfit_search(), scale all candidates by the first one's depth ratio.
    ratio = sum(obs$value) / sum(sim_dat(kars[[1]])$value)
    unlist(.mclapply_strict(kars, function(gw) {
        s = sim_dat(gw)[, value := value * ratio]
        compdats(copy(obs_dat), s, theta = theta, ifsum = TRUE, area0 = area0)
    }, mc.cores = mc.cores))
}

## Progress line for batch logs: time, step and the parent's resident memory.
.identifiability_log <- function(label, ...) {
    rss = tryCatch(as.numeric(gsub('[^0-9]', '', grep('^VmRSS', readLines('/proc/self/status'), value = TRUE))) / 2^20,
                   error = function(e) NA_real_)
    message(format(Sys.time(), '%H:%M:%S'), ' [', label, '] ', paste0(...), sprintf(' (RSS %.1f GB)', rss))
}

## Sample one candidate set for a truth and compute each candidate's edit
## distance to the truth (and Hi-C NLL when requested). Other modalities re-rank
## this same candidate set, so edit distances are computed once per truth and
## assays are compared on identical candidates.
identifiability_candidates <- function(ctx,
                                       truth,
                                       nsample,
                                       depth = 10,
                                       theta = 2,
                                       run_hic = TRUE,
                                       mc.cores = 1,
                                       label = '') {
    kars = sample.gwalks(ctx$graph, nsample, verbose = FALSE, mc.cores = mc.cores)
    .identifiability_log(label, 'sampled ', length(kars), ' candidates')
    hic_nll = NULL
    if (run_hic) {
        hic_nll = identifiability_hic_nlls(kars, truth, ctx, depth = depth,
                                           theta = theta, mc.cores = mc.cores)
        .identifiability_log(label, 'hic scored')
    }
    edit = unlist(.mclapply_strict(kars, function(gw) {
        edit_dist_cpp(truth$true_gw, gw, graph = ctx$graph, use.nodes = ctx$use.nodes,
                      node.widths = ctx$edit_widths, split.ends = ctx$split_ends)
    }, mc.cores = mc.cores))
    list(
        kars = kars,
        hic_nll = hic_nll,
        candidates = data.table(
            cand_id = seq_along(kars),
            hash = vapply(kars, function(gw) gw$hash, character(1)),
            edit = edit
        )
    )
}

## Rank the candidate set by one modality's NLL. Returns one row per candidate;
## scoring is left to downstream analysis.
.identifiability_rank <- function(candidates, nll, n_reads = NA_integer_) {
    ranked = copy(candidates)
    ranked[, nll := nll]
    ranked[, n_reads := n_reads]
    setorder(ranked, nll, cand_id)
    ranked[, rank := seq_len(.N)]
    setcolorder(ranked, c('rank', 'cand_id', 'nll', 'edit', 'hash'))
    ranked
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

simulate_identifiability_longread_words <- function(gw,
                                                    readL,
                                                    depth = 10,
                                                    use.nodes = NULL,
                                                    context_width = NULL) {
    ## Reads that touch no event node collapse to the empty word. Drop them so
    ## the observed reads, whose count is scaled to the event context below,
    ## are drawn from reads overlapping the event; in full-graph mode they are
    ## otherwise ~99% of the read mass and the observation carries no signal.
    reads = .identifiability_event_reads(gw, readL, use.nodes)
    if (!nrow(reads)) {
        return(character())
    }

    if (is.null(context_width)) {
        context_width = sum(width(gw$graph$footprint))
    }
    context = context_width + 2 * readL
    n_reads = max(1L, as.integer(round(depth / (readL / context))))

    sample(reads$words, size = n_reads, prob = reads$nums, replace = TRUE)
}

identifiability_longread_candidate_nlls <- function(candidates,
                                                    obs_words,
                                                    readL,
                                                    use.nodes = NULL,
                                                    background = 1e-10,
                                                    mc.cores = 1) {
    if (!length(obs_words)) {
        return(rep(Inf, length(candidates)))
    }

    obs_counts = table(obs_words)
    obs_vocab = names(obs_counts)

    ## Condition on reads overlapping the event, matching the observed reads.
    ## Workers return only what the NLL needs (the candidate's counts of the
    ## observed words, its total, and its vocabulary for the shared V), not the
    ## read tables. The vocabulary is returned as 64-bit word hashes: for
    ## high-copy amplicons, long reads span hundreds of nodes and the words are
    ## kb-long strings that would otherwise accumulate in (and never be released
    ## from) the parent, which every later fork then copies.
    hash_words = digest::getVDigest(algo = 'xxhash64')
    summ = .mclapply_strict(candidates, function(gw) {
        rt = .identifiability_event_reads(gw, readL, use.nodes)
        counts = setNames(rt$nums, rt$words)[obs_vocab]
        counts[is.na(counts)] = 0
        list(obs = unname(counts), total = sum(rt$nums),
             words = if (nrow(rt)) hash_words(rt$words, serialize = FALSE) else character())
    }, mc.cores = mc.cores)

    V = length(unique(c(hash_words(obs_vocab, serialize = FALSE),
                        unlist(lapply(summ, `[[`, 'words'), use.names = FALSE))))
    for (i in seq_along(summ)) summ[[i]]$words = NULL
    invisible(gc())

    vapply(summ, function(x) {
        if (x$total == 0) {
            return(Inf)
        }
        p_obs = (x$obs + background) / (x$total + background * V)
        -sum(as.numeric(obs_counts) * log(p_obs))
    }, numeric(1))
}

run_graph_identifiability <- function(graph,
                                      footprint = NULL,
                                      event = NULL,
                                      nsamples = c(25, 100, 500),
                                      truth_reps = 1,
                                      rep_ids = NULL,
                                      seeds = NULL,
                                      graph_mode = 'footprint',
                                      modalities = c('hic', 'longread'),
                                      longread_lengths = c(1e5, 1e7),
                                      depth = 10,
                                      theta = 2,
                                      hic_binsperft = 100,
                                      footprint_buffer_fraction = 1/5,
                                      footprint_buffer_min = 1e4,
                                      edit_mode = 'fragment',
                                      longread_background = 1e-10,
                                      wholegenome = NULL,
                                      mc.cores = 1) {
    modalities = match.arg(modalities, choices = c('hic', 'longread'), several.ok = TRUE)
    run_hic = 'hic' %in% modalities
    run_longread = 'longread' %in% modalities && length(longread_lengths) > 0
    ## rep_ids names the truths to run (default 1..truth_reps), so one truth of
    ## an event can run on its own as a shard. seeds, if given, has one seed per
    ## rep_id and is set before that truth is sampled, so a truth does not depend
    ## on which other truths ran in the same job.
    if (is.null(rep_ids)) rep_ids = seq_len(truth_reps)
    rep_ids = as.integer(rep_ids)
    if (!is.null(seeds) && length(seeds) != length(rep_ids)) {
        stop('seeds must have one entry per rep_id')
    }
    ctx = identifiability_context(
        graph = graph,
        footprint = footprint,
        graph_mode = graph_mode,
        wholegenome = wholegenome,
        hic_binsperft = hic_binsperft,
        footprint_buffer_fraction = footprint_buffer_fraction,
        footprint_buffer_min = footprint_buffer_min,
        edit_mode = edit_mode
    )
    ranked = list()
    label = if (!is.null(event$event_id)) as.character(event$event_id) else 'event'
    .identifiability_log(label, 'start: ', length(ctx$graph$nodes), ' nodes, target ',
                         ctx$target_width, ' bp, pix ', ctx$pix.size)
    for (k in seq_along(rep_ids)) {
        rep_id = rep_ids[k]
        invisible(gc())
        if (!is.null(seeds)) set.seed(seeds[k])
        truth = make_identifiability_truth(ctx, depth = depth, theta = theta,
                                           simulate_hic = run_hic)
        for (nsample in nsamples) {
            cand = identifiability_candidates(ctx, truth, nsample, depth = depth,
                                              theta = theta, run_hic = run_hic,
                                              mc.cores = mc.cores,
                                              label = paste0(label, ' rep ', rep_id))
            .identifiability_log(paste0(label, ' rep ', rep_id), 'edit distances done')
            if (run_hic) {
                r = .identifiability_rank(cand$candidates, cand$hic_nll)
                ranked[[length(ranked) + 1L]] = .identifiability_add_metadata(
                    r, event, ctx, truth, rep_id, nsample, graph_mode, 'hic', NA_real_)
            }
            if (run_longread) {
                for (readL in longread_lengths) {
                    obs_words = simulate_identifiability_longread_words(
                        truth$true_gw,
                        readL = readL,
                        depth = depth,
                        use.nodes = ctx$use.nodes,
                        context_width = ctx$target_width
                    )
                    nll = identifiability_longread_candidate_nlls(
                        cand$kars,
                        obs_words = obs_words,
                        readL = readL,
                        use.nodes = ctx$use.nodes,
                        background = longread_background,
                        mc.cores = mc.cores
                    )
                    r = .identifiability_rank(cand$candidates, nll, length(obs_words))
                    ranked[[length(ranked) + 1L]] = .identifiability_add_metadata(
                        r, event, ctx, truth, rep_id, nsample, graph_mode, 'longread', readL)
                    .identifiability_log(paste0(label, ' rep ', rep_id), 'long reads ', readL, ' scored')
                }
            }
        }
    }
    list(
        ranked = rbindlist(ranked, fill = TRUE),
        context = ctx
    )
}

run_identifiability <- function(graphs,
                                footprints = NULL,
                                events = NULL,
                                nsamples = c(25, 100, 500),
                                truth_reps = 1,
                                rep_ids = NULL,
                                seeds = NULL,
                                graph_mode = 'footprint',
                                modalities = c('hic', 'longread'),
                                longread_lengths = c(1e5, 1e7),
                                depth = 10,
                                theta = 2,
                                hic_binsperft = 100,
                                footprint_buffer_fraction = 1/5,
                                footprint_buffer_min = 1e4,
                                edit_mode = 'fragment',
                                longread_background = 1e-10,
                                wholegenome = NULL,
                                core_budget = 40L,
                                outer_cores = NULL) {
    if (!is.list(graphs)) {
        graphs = list(graphs)
    }
    n = length(graphs)
    if (is.null(footprints)) {
        footprints = vector('list', n)
    }
    if (!is.list(footprints)) {
        footprints = rep(list(footprints), n)
    }
    if (is.null(events)) {
        events = vector('list', n)
    } else if (is.data.frame(events) || is.data.table(events)) {
        events = split(as.data.table(events), seq_len(nrow(events)))
    } else if (!is.list(events)) {
        events = rep(list(events), n)
    }
    ## rep_ids/seeds: NULL, or a list with one integer vector per graph (see
    ## run_graph_identifiability).
    if (!is.null(rep_ids) && !is.list(rep_ids)) rep_ids = rep(list(rep_ids), n)
    if (!is.null(seeds) && !is.list(seeds)) seeds = rep(list(seeds), n)
    stopifnot(length(footprints) == n, length(events) == n,
              is.null(rep_ids) || length(rep_ids) == n,
              is.null(seeds) || length(seeds) == n)
    if (is.null(outer_cores)) {
        ## Events can vary substantially in runtime, so default to processing
        ## events serially and spending the full core budget inside each event.
        ## Users can still request event-level parallelism explicitly via
        ## outer_cores, with mc.preschedule=FALSE below for load balancing.
        outer_cores = 1L
    }
    outer_cores = min(n, max(1L, as.integer(outer_cores)))
    inner_cores = max(1L, floor(as.integer(core_budget) / outer_cores))
    if (outer_cores * inner_cores > core_budget) {
        stop('Nested core configuration exceeds core_budget')
    }
    message('Running ', n, ' graphs with outer_cores=', outer_cores,
            ', inner_cores=', inner_cores,
            ', max concurrent cores=', outer_cores * inner_cores)
    t0 = Sys.time()
    run_one = function(i) {
        run_graph_identifiability(
            graph = graphs[[i]],
            footprint = footprints[[i]],
            event = events[[i]],
            nsamples = nsamples,
            truth_reps = truth_reps,
            rep_ids = rep_ids[[i]],
            seeds = seeds[[i]],
            graph_mode = graph_mode,
            modalities = modalities,
            longread_lengths = longread_lengths,
            depth = depth,
            theta = theta,
            hic_binsperft = hic_binsperft,
            footprint_buffer_fraction = footprint_buffer_fraction,
            footprint_buffer_min = footprint_buffer_min,
            edit_mode = edit_mode,
            longread_background = longread_background,
            wholegenome = wholegenome,
            mc.cores = inner_cores
        )
    }
    if (outer_cores > 1L) {
        out = pbmclapply(seq_len(n), run_one, mc.cores = outer_cores, mc.preschedule = FALSE)
    } else {
        out = lapply(seq_len(n), run_one)
    }
    runtime = difftime(Sys.time(), t0, units = 'mins')
    ranked = rbindlist(lapply(out, `[[`, 'ranked'), fill = TRUE)
    ranked[, runtime_min := as.numeric(runtime)]
    list(
        ranked = ranked,
        raw_output = out,
        config = list(
            n_graphs = n,
            nsamples = nsamples,
            truth_reps = truth_reps,
            graph_mode = graph_mode,
            modalities = modalities,
            longread_lengths = longread_lengths,
            footprint_buffer_fraction = footprint_buffer_fraction,
            footprint_buffer_min = footprint_buffer_min,
            edit_mode = edit_mode,
            core_budget = core_budget,
            outer_cores = outer_cores,
            inner_cores = inner_cores,
            runtime_min = as.numeric(runtime)
        )
    )
}

## Score one assay's ranking of one candidate set. The plausible set T is all
## candidates whose NLL is within `delta` nats of the minimum (plus a relative
## `tol` for floating-point ties); the assay is modelled as picking uniformly
## from T. For each ball radius, p_ball_<cutoff> is the probability that the
## pick lies within that edit distance of the truth, and prior_ball_<cutoff>
## the probability that a random candidate does.
.identifiability_score_one <- function(nll, edit, radius, delta = 1, tol = 1e-8) {
    m = min(nll)
    if (is.finite(m)) {
        plausible = !is.na(nll) & nll - m <= delta + tol * max(abs(m), 1)
    } else {
        plausible = rep(TRUE, length(nll))
    }
    c(list(n_candidates = length(edit),
           n_plausible = sum(plausible),
           frac_plausible = mean(plausible),
           min_edit = min(edit)),
      setNames(lapply(radius, function(r) mean(edit[plausible] <= r)), paste0('p_ball_', names(radius))),
      setNames(lapply(radius, function(r) mean(edit <= r)), paste0('prior_ball_', names(radius))))
}

## Ball radii for one candidate set. ball = 'size': cutoff x event size, the
## copy-number-weighted bp of the event nodes (ranked needs an `event_size`
## column). ball = 'quantile': the cutoff-quantile of the candidates' edit
## distances (0 = the closest candidates); with ties at the radius the ball
## holds more than that fraction of candidates. Names are 100 x cutoff.
.identifiability_ball_radius <- function(edit, event_size, ball, cutoffs) {
    r = if (ball == 'size') cutoffs * event_size[1]
        else quantile(edit, cutoffs, type = 1, names = FALSE)
    setNames(r, sub('\\.', '_', format(100 * cutoffs, trim = TRUE, drop0trailing = TRUE)))
}

## Score the `ranked` output of run_identifiability(). One row per truth
## (event x rep_id x nsample) x assay on the full candidate set, with the
## probability that the assay's pick (and a random candidate) lands within
## each edit-distance ball around the truth. `ball` picks how radii are set
## (see .identifiability_ball_radius) and `cutoffs` the radii. `delta` is the
## NLL tolerance (nats) defining the plausible set; 0 = exact ties. With
## `subsample_n`, also rescores random candidate subsets of each size
## (`subsample_reps` draws, averaged), using the same candidate subsets for
## every assay, to check convergence in the number of sampled candidates.
score_identifiability <- function(ranked,
                                  ball = c('size', 'quantile'),
                                  cutoffs = c(0.05, 0.1, 0.25),
                                  delta = 1,
                                  subsample_n = NULL,
                                  subsample_reps = 50,
                                  tol = 1e-8,
                                  seed = 1,
                                  mc.cores = 1) {
    ball = match.arg(ball)
    dt = as.data.table(ranked)
    if (ball == 'size' && !('event_size' %in% names(dt))) {
        stop("ball = 'size' needs an event_size column in ranked (see identifiability_event_sizes())")
    }
    if (!('event_size' %in% names(dt))) dt[, event_size := NA_real_]
    cand_cols = c('rank', 'cand_id', 'nll', 'edit', 'hash', 'n_reads', 'runtime_min')
    unit_cols = setdiff(names(dt), cand_cols)
    assay_cols = c('modality', 'readL')
    set_cols = setdiff(unit_cols, assay_cols)
    score = function(nll, edit, event_size) {
        .identifiability_score_one(nll, edit,
                                   .identifiability_ball_radius(edit, event_size, ball, cutoffs),
                                   delta, tol)
    }

    full = dt[, score(nll, edit, event_size), by = unit_cols]
    full[, `:=`(subsampled = FALSE, subsample_n = n_candidates, subsample_reps = 1L)]
    if (is.null(subsample_n) || !length(subsample_n)) {
        return(full)
    }
    metric_cols = setdiff(names(full), c(unit_cols, 'subsampled', 'subsample_n', 'subsample_reps'))

    dt[, set_grp := .GRP, by = set_cols]
    groups = split(dt, by = 'set_grp', keep.by = FALSE)
    set.seed(seed)
    seeds = sample.int(.Machine$integer.max, length(groups))
    sub = rbindlist(pbmclapply(seq_along(groups), function(g) {
        set.seed(seeds[g])
        gdt = groups[[g]]
        ids = unique(gdt$cand_id)
        rbindlist(lapply(subsample_n[subsample_n < length(ids)], function(s) {
            draws = rbindlist(lapply(seq_len(subsample_reps), function(r) {
                data.table(sub_rep = r, cand_id = sample(ids, s))
            }))
            sc = merge(gdt, draws, by = 'cand_id', allow.cartesian = TRUE)[
              , score(nll, edit, event_size), by = c(unit_cols, 'sub_rep')]
            sc[, c(lapply(.SD, mean, na.rm = TRUE), list(subsample_reps = .N)),
               by = unit_cols, .SDcols = metric_cols][
              , `:=`(subsampled = TRUE, subsample_n = s)]
        }), fill = TRUE)
    }, mc.cores = mc.cores), fill = TRUE)

    out = rbindlist(list(full, sub), fill = TRUE, use.names = TRUE)
    setorderv(out, c(set_cols, assay_cols, 'subsample_n'))
    out[]
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
## an n_obs x length(kars) matrix. Per map this is identifiability_hic_nlls()
## (candidates simulated once, scaled to the map's depth by the first
## candidate's total, negative binomial likelihood summed over pixels),
## vectorised over maps. Pixels are unweighted, as there: the observed map
## reaches compdats() without a widthprod column, so every pixel gets area0.
local_hic_nll_matrix <- function(kars, ctx, n_obs = 100, depth = 10, theta = 2) {
    true_expected = forward_simulate(kars[[1]], target_region = ctx$target_region,
                                     pix.size = ctx$pix.size, depth = depth)
    obs = make_noisydat(true_expected, nsamp = n_obs, theta = theta)
    gr = true_expected$gr
    template = make_template_dat(gr)[, .(i, j, id)]
    setkey(template, id)
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
        out[, k] = -colSums(lp, na.rm = TRUE)
    }
    out
}

## Long-read NLL of n_obs read sets drawn from kars[[1]] (the truth) under
## every candidate, as identifiability_longread_candidate_nlls(): reads
## restricted to those touching the event, words hashed, background smoothing
## over the shared vocabulary. n_obs x length(kars) matrix.
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
