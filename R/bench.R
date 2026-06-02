# fastKar — benchmark helpers
#
# Tracks per-phase runtime for traverse and hash so we can compare baseline
# (traverse_graph_cpp + hash_snodelist) against new implementations as we
# land them. Verifies correctness via partition-equivalence: two configs are
# considered equivalent iff they group samples into the same equivalence
# classes (regardless of how they spell their hash strings).

#' Standard synthetic benchmark battery.
#'
#' Returns a named list of gGraph objects. `gimme_bfb` is stochastic (random
#' BFB breakpoints), so we seed before generating each graph to make sizes
#' reproducible across runs — important for tracking speedups over time.
#' @export
bench_synthetic_graphs <- function(seed = 42L) {
  set.seed(seed); g_bfb_small  <- gimme_bfb(1e6, 8,  cycles = 5)$graph
  set.seed(seed); g_bfb_medium <- gimme_bfb(1e6, 16, cycles = 8)$graph
  list(
    bfb_small    = g_bfb_small,
    bfb_medium   = g_bfb_medium,
    repdup_cis1  = makerepdup(1e6, geometry = "cis1")$graph,
    repdup_trans = makerepdup(1e6, geometry = "trans")$graph
  )
}

#' Load real-data graphs from bench/graphs/*.rds (if present).
#'
#' Each .rds may deserialize to:
#'   * a gGraph directly
#'   * an object with a `$graph` slot that is a gGraph
#'   * a list whose elements are gGraphs (each becomes its own bench entry,
#'     named `<basename>__<elemname-or-index>`)
#' Empty list if the directory is missing or empty.
#' @export
bench_real_graphs <- function(dir = "bench/graphs") {
  if (!dir.exists(dir)) return(list())
  paths <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(paths)) return(list())

  out <- list()
  for (p in paths) {
    base <- tools::file_path_sans_ext(basename(p))
    obj  <- readRDS(p)

    if (inherits(obj, "gGraph")) {
      out[[base]] <- obj
    } else if (is.list(obj) && !is.null(obj$graph) && inherits(obj$graph, "gGraph")) {
      out[[base]] <- obj$graph
    } else if (is.list(obj)) {
      gg_idx <- which(vapply(obj, inherits, logical(1), what = "gGraph"))
      if (!length(gg_idx)) {
        stop("File ", p, " contained no gGraph objects")
      }
      labels <- if (!is.null(names(obj))) names(obj)[gg_idx] else as.character(gg_idx)
      labels[is.na(labels) | labels == ""] <- as.character(gg_idx[is.na(labels) | labels == ""])
      for (k in seq_along(gg_idx)) {
        out[[paste0(base, "__", labels[k])]] <- obj[[gg_idx[k]]]
      }
    } else {
      stop("File ", p, " did not contain a gGraph, an object with $graph, or a list of gGraphs")
    }
  }
  out
}

#' Default benchmark configs: baseline, traverse-only, traverse+hash, batched.
#'
#' A config is a list with:
#'   * traverse     : per-sample function `(internal.edges, loose.ends) -> walks`
#'                    (when batch = FALSE, the default)
#'                    OR
#'                    batch function `(internal.edges, perms_list, loose.ends)
#'                    -> list of walks` (when batch = TRUE)
#'   * hash         : per-walk function `(snode.id, circular) -> chr`
#'   * batch        : optional logical (default FALSE)
#' @export
bench_default_configs <- function() {
  list(
    v1           = list(traverse = traverse_graph_cpp,
                        hash     = hash_snodelist),                   # baseline
    v2t          = list(traverse = traverse_graph_v2_cpp,
                        hash     = hash_snodelist),                   # v2 traversal, R hash
    v2t_h2       = list(traverse = traverse_graph_v2_cpp,
                        hash     = hash_karyotype_cpp),               # v2 traversal + v2 hash
    v2t_batch_h2 = list(traverse = traverse_graph_v2_batch_cpp,
                        hash     = hash_karyotype_cpp,
                        batch    = TRUE)                              # batched + v2 hash
  )
}

# Internal: partition fingerprint — for each element, the position of its
# first occurrence in the vector. Two hash vectors that induce the same
# equivalence relation will produce identical match(h, h) regardless of
# what their hash strings spell.
.bench_partition <- function(h) match(h, h)

#' Compare configs (each a (traverse, hash) pair) on a single graph.
#'
#' Verifies that all configs produce the same partition of the N samples
#' into equivalence classes; aborts on mismatch.
#'
#' @param gg          a gGraph
#' @param N           number of permutation samples to time
#' @param configs     named list of list(traverse, hash); see bench_default_configs()
#' @param graph_name  string label for result rows
#' @param seed        random seed for permutation generation
#' @export
bench_one_graph <- function(gg,
                            N          = 1000,
                            configs    = bench_default_configs(),
                            graph_name = NA_character_,
                            seed       = 1) {

  wiring         <- gg.to.wiring(gg)
  internal.edges <- wiring$internal.edges
  loose.ends     <- wiring$loose.ends

  set.seed(seed)
  perms <- lapply(seq_len(N), function(i) {
    internal.edges[, if (.N > 1) sample(right, .N) else right, by = n]$V1
  })

  out          <- list()
  ref_part     <- NULL
  ref_name     <- names(configs)[1]

  for (cfg_name in names(configs)) {
    cfg <- configs[[cfg_name]]
    is_batch <- isTRUE(cfg$batch)

    gc(verbose = FALSE)
    t_traverse <- system.time({
      walks <- if (is_batch) {
        cfg$traverse(internal.edges, perms, loose.ends)
      } else {
        lapply(perms,
               function(p) cfg$traverse(internal.edges[, right := p], loose.ends))
      }
    })

    gc(verbose = FALSE)
    t_hash <- system.time({
      hashes <- vapply(walks,
                       function(w) cfg$hash(w$snode.id, w$circular),
                       character(1))
    })

    part <- .bench_partition(hashes)
    if (is.null(ref_part)) {
      ref_part <- part
    } else if (!identical(part, ref_part)) {
      stop(sprintf(
        "Sample partitions differ between '%s' and '%s' on '%s' (N=%d): %d vs %d unique classes.",
        ref_name, cfg_name, graph_name, N,
        length(unique(ref_part)), length(unique(part))))
    }

    elapsed_t <- unname(t_traverse["elapsed"])
    elapsed_h <- unname(t_hash["elapsed"])
    out[[length(out) + 1L]] <- data.table::data.table(
      graph   = graph_name,
      N       = N,
      version = cfg_name,
      phase   = c("traverse", "hash", "total"),
      seconds = c(elapsed_t, elapsed_h, elapsed_t + elapsed_h)
    )
  }

  data.table::rbindlist(out)
}

#' Run the full benchmark battery and append results to history.
#'
#' @param Ns           vector of sample sizes per graph
#' @param configs      named list of (traverse, hash) configs
#' @param history_path rds file collecting results across runs
#' @param graphs       named list of graphs; if NULL, uses synthetic+real
#' @export
bench_all <- function(Ns           = c(1000, 10000),
                      configs      = bench_default_configs(),
                      history_path = "bench/history.rds",
                      graphs       = NULL) {

  if (is.null(graphs)) {
    graphs <- c(bench_synthetic_graphs(), bench_real_graphs())
  }
  if (!length(graphs)) stop("No graphs to benchmark.")

  results <- list()
  for (gname in names(graphs)) {
    for (N in Ns) {
      message(sprintf("[bench] %s N=%d", gname, N))
      r <- bench_one_graph(graphs[[gname]], N = N, configs = configs,
                           graph_name = gname)
      results[[length(results) + 1L]] <- r
    }
  }
  results <- data.table::rbindlist(results)

  results[, date       := Sys.time()]
  results[, git_commit := tryCatch(system("git rev-parse HEAD", intern = TRUE),
                                   error = function(e) NA_character_)]
  results[, n_configs  := length(configs)]

  if (!dir.exists(dirname(history_path))) {
    dir.create(dirname(history_path), recursive = TRUE)
  }
  prior    <- if (file.exists(history_path)) readRDS(history_path) else NULL
  combined <- data.table::rbindlist(list(prior, results), fill = TRUE)
  saveRDS(combined, history_path)

  results
}

#' Wide summary: rows = (graph, N), columns = version_phase.
#' @export
bench_summary <- function(results) {
  data.table::dcast(results, graph + N ~ version + phase,
                    value.var = "seconds")
}
