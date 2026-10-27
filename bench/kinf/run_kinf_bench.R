# run_kinf_bench.R: convergence benchmark of the kinf search (infer_karyotype) on small CCLE events with
# simulated truths. A machine check, not a hardness study: 5 cut events (bench/kinf/graphs, from
# make_kinf_bench_graphs.R) x 3 random truths (sample.gwalks, hashed seeds) x depths 1 and 10, noisy maps
# make_noisymap(forward_simulate(truth), 2); PT and greedy, 16 chains x 5000 steps each.
# Success = best NLL <= truth NLL + 1e-6 (a negative gap: another karyotype fits the data better, i.e. the truth
# is not identifiable at that depth).
# Usage (from the repo root; needs install_kinf()):
#   Rscript bench/kinf/run_kinf_bench.R pilot [gpu]      one problem (first event, truth 1, depth 10, PT)
#   Rscript bench/kinf/run_kinf_bench.R run [gpus]       all problems; gpus e.g. 0,1,2,3 (one run per GPU at a time)
#   Rscript bench/kinf/run_kinf_bench.R plot [outdir]    figure from results/ (default ~/public_html/karyotype_inference)
# Writes bench/kinf/results/{summary,traces,curves}.tsv (pilot: pilot_*.tsv); curves = best NLL - truth NLL so far
# over all chains every 10 proposals from the random starts (from infer_karyotype(history_every = 10)), kept at
# its change points (+ last step). Raw run folders go to tempdir() and are
# deleted.
suppressMessages(devtools::load_all('.', quiet = TRUE))
data.table::setDTthreads(1)
bdir = 'bench/kinf'
args = commandArgs(trailingOnly = TRUE)
stage = if (length(args)) args[1] else 'pilot'

seed_of = function(...) strtoi(substr(digest::digest(paste(...), algo = 'xxhash32', serialize = FALSE), 1, 7), 16L)
load_event = function(f){
  x = readRDS(f)
  x$graph = gG(nodes = x$nodes, edges = x$edges)
  x
}
depths = c(1, 10); ntruth = 3; chains = 16; steps = 5000  # pilot: truth NLL reached by step 1000 (A2058_928)

design = function(){
  fs = sort(list.files(file.path(bdir, 'graphs'), pattern = '\\.rds$', full.names = TRUE))
  d = CJ(file = fs, truth = seq_len(ntruth), depth = depths, method = c('pt', 'greedy'), sorted = FALSE)
  d[, event_id := sub('\\.rds$', '', basename(file))]
  d[, run_id := paste(event_id, paste0('t', truth), paste0('d', depth), method, sep = '_')]
  d[]
}

run_one = function(r, gpu, cores){
  x = load_event(r$file)
  set.seed(seed_of('truth', r$event_id, r$truth))
  truth = sample.gwalks(x$graph, 1, return.gw = TRUE, verbose = FALSE)[[1]]
  set.seed(seed_of('noise', r$event_id, r$truth, r$depth))
  gm = make_noisymap(forward_simulate(truth, target_region = x$ft, depth = r$depth, pix.size = x$pix.size), 2)[[1]]
  wd = tempfile(paste0('kinfbench_', r$run_id, '_'))
  t0 = Sys.time()
  res = infer_karyotype(x$graph, gm, depth = r$depth, pix.size = x$pix.size, target_region = x$ft,
                        method = r$method, chains = chains, steps = steps, seed = r$truth, truth = truth,
                        gpu = gpu, cores = cores, workdir = wd, history_every = 10, verbose = FALSE)
  tot = as.numeric(difftime(Sys.time(), t0, units = 'secs'))
  s = res$summary
  tr = res$trace[, .(steps, best)]
  cv = res$history[, .(nll = min(nll)), keyby = steps][, .(steps, best = cummin(nll))]
  cv = cv[c(TRUE, diff(best) != 0) | steps == max(steps)]  # step function: change points + last step
  hit = tr[best <= 1e-6]$steps
  unlink(wd, recursive = TRUE)
  list(summary = data.table(run_id = r$run_id, event_id = r$event_id, type = x$type, truth = r$truth,
                            depth = r$depth, method = r$method, n_var_nodes = s$n_var_nodes, chains = s$chains,
                            steps = s$steps, nll_truth = s$nll_truth, best_nll = s$best_nll,
                            gap = s$gap_to_truth, success = s$gap_to_truth <= 1e-6,
                            exact_recovery = s$exact_recovery, steps_to_truth = if (length(hit)) min(hit) else NA_real_,
                            search_s = s$wall_s, total_s = tot, ms_per_proposal = s$ms_per_proposal, device = s$device),
       trace = cbind(data.table(run_id = r$run_id), tr), curve = cbind(data.table(run_id = r$run_id), cv))
}

run_all = function(d, gpus, prefix){
  slot = function(k) list(gpu = gpus[(k - 1) %% length(gpus) + 1], cores = 8 * ((k - 1) %% length(gpus)) + 0:7 + 64)
  out = parallel::mclapply(seq_len(nrow(d)), function(k){
    s = slot(k)
    tryCatch(run_one(d[k], s$gpu, s$cores), error = function(e) list(error = paste(d$run_id[k], conditionMessage(e))))
  }, mc.cores = length(gpus), mc.preschedule = TRUE)
  err = unlist(lapply(out, `[[`, 'error'))
  if (length(err)) message('errors:\n', paste(err, collapse = '\n'))
  ok = out[vapply(out, function(o) is.null(o$error), logical(1))]
  summ = rbindlist(lapply(ok, `[[`, 'summary')); trc = rbindlist(lapply(ok, `[[`, 'trace'))
  fwrite(rbindlist(lapply(ok, `[[`, 'curve')), file.path(bdir, 'results', paste0(prefix, 'curves.tsv')), sep = '\t')
  fwrite(summ, file.path(bdir, 'results', paste0(prefix, 'summary.tsv')), sep = '\t')
  fwrite(trc, file.path(bdir, 'results', paste0(prefix, 'traces.tsv')), sep = '\t')
  print(summ[, .(run_id, n_var_nodes, gap = signif(gap, 3), success, exact_recovery, steps_to_truth,
                 search_s = round(search_s, 1), device)], nrows = 200)
  invisible(summ)
}

if (stage == 'pilot'){
  gpu = if (length(args) > 1) as.integer(args[2]) else 0L
  d = design()[truth == 1 & depth == 10 & method == 'pt'][1]
  run_all(d, gpu, 'pilot_')
}

if (stage == 'run'){
  gpus = if (length(args) > 1) as.integer(strsplit(args[2], ',')[[1]]) else 0L
  d = design()
  message(nrow(d), ' runs on GPUs ', paste(gpus, collapse = ','))
  summ = run_all(d, gpus, '')
  print(summ[, .(runs = .N, success = sum(success), exact = sum(exact_recovery), median_s = median(search_s)),
             by = .(method, depth)])
}

if (stage == 'plot'){
  out = if (length(args) > 1) args[2] else '~/public_html/karyotype_inference'
  summ = fread(file.path(bdir, 'results', 'summary.tsv')); trc = fread(file.path(bdir, 'results', 'curves.tsv'))
  trc = merge(trc, summ[, .(run_id, event_id, type, truth, depth, method)], by = 'run_id')
  trc[, panel := paste0(sub('^G32831\\.', '', event_id), ' (', type, ')')]
  trc[, depth_lab := factor(paste0(depth, 'x'), levels = paste0(sort(unique(depth)), 'x'))]
  symlog = scales::trans_new('symlog', function(x) sign(x) * log10(1 + abs(x) / 1e-3),
                             function(y) sign(y) * 1e-3 * (10^abs(y) - 1))
  brk = c(-10, -0.1, 0, 0.1, 10, 1e3, 1e5)
  p = ggplot(trc, aes(steps, best, group = run_id, colour = depth_lab, linetype = method)) +
    geom_hline(yintercept = 0, colour = 'grey50', linewidth = 0.3) +
    ggplot2::geom_step(linewidth = 0.6, alpha = 0.85) +
    facet_wrap(~panel, ncol = 3) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 10, base = 10),
                                breaks = c(0, 10, 100, 1000, 5000)) +
    ggplot2::scale_y_continuous(trans = symlog, breaks = brk, labels = format(brk, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    ggplot2::scale_colour_manual(values = c('1x' = '#2a7ab0', '10x' = '#c0562f')) +
    ggplot2::scale_linetype_manual(values = c(pt = 'solid', greedy = 'dashed')) +
    labs(x = 'proposals per chain', y = 'best NLL - truth NLL (symlog)', colour = 'depth', linetype = 'search',
         title = 'kinf search convergence on CCLE cut events (3 simulated truths each)') +
    theme_classic(base_size = 10) + theme(strip.background = element_blank())
  dir.create(path.expand(out), recursive = TRUE, showWarnings = FALSE)
  for (ext in c('pdf', 'png')) ggplot2::ggsave(file.path(path.expand(out), paste0('kinf_bench.', ext)), p,
                                               width = 10, height = 6, dpi = 150)
  tab = summ[, .(runs = .N, success = sum(success), exact = sum(exact_recovery),
                 median_gap = signif(median(gap), 3), median_search_s = round(median(search_s), 1)),
             by = .(event_id, depth, method)]
  fwrite(tab, file.path(bdir, 'results', 'summary_by_event.tsv'), sep = '\t')
  print(tab)
}
