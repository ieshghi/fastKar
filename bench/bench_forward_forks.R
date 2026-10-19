## F5: does forking over walks still pay once the walk kernel is in C++?
## Times the kernel on all walks in one process vs split over 4 forked workers (summed), and
## forward_simulate with mc.cores 1 vs 4 for the oracle and fastKarOptim.
##   Rscript bench/bench_forward_forks.R [optim_path]
args = commandArgs(TRUE)
optim_path = normalizePath(if (length(args)) args[1] else '~/git/fastKar_optim')
suppressMessages({library(data.table); setDTthreads(1); library(parallel)})
suppressMessages(pkgload::load_all(Sys.getenv('FASTKAR_ORACLE', '~/git/fastKar_oracle'), export_all=FALSE, attach=FALSE, quiet=TRUE))
suppressMessages(pkgload::load_all(optim_path, export_all=FALSE, attach=FALSE, quiet=TRUE))
ns = asNamespace('fastKarOptim'); ref = function(n) get(n, envir=asNamespace('fastKar'))
fd = file.path(optim_path, 'tests', 'testthat', 'fixtures')
lf = function(name){fx = readRDS(file.path(fd,paste0(name,'.rds'))); fx$gw = lapply(fx$kars, function(k){gw=gGnome::gW(graph=fx$gg, snode.id=k$snode.id, circular=k$circular); if(!is.null(k$cn)) gw$set(cn=k$cn); gw}); fx}
for (nm in c('ccle_dm','ccle_chromothripsis')) { fx = lf(nm); gw = fx$gw[[1]]
  p = ns$prep_for_sim(gw, fx$ft, 1e4); tt = p$tiled.target
  sn = gw$snode.id; hit = ns$walks_hit_target(sn, tt$node.id); sn = sn[hit]; circ = gw$circular[hit]; cn = rep(1, length(sn))
  rows = ns$walk_rows(sn, tt, p$widthdt); ka = ns$forward_kernel_args(0)
  run = function(w) { r = seq.int(rows$walk_ptr[w[1]]+1, length.out=rows$walk_ptr[max(w)+1]-rows$walk_ptr[w[1]])
    ns$forward_walks_cpp(as.integer(rows$walk_ptr[c(w, max(w)+1)] - rows$walk_ptr[w[1]]), rows$orig[r], rows$width[r], rows$int_width, as.logical(circ[w]), cn[w], nrow(tt), ka$z, ka$lmin, ka$slope, ka$maxval, ka$maxval_density, 0, 1, ka$diagf) }
  t1 = system.time(a1 <- run(seq_along(sn)))[3]
  t4 = system.time({ ch = split(seq_along(sn), cut(seq_along(sn), 4, labels=FALSE)); parts = mclapply(ch, run, mc.cores=4); a4 = Reduce(`+`, parts) })[3]
  cat(sprintf('%s 1e4: %d walks, %d tiles; kernel 1 core %.2fs, 4 forks %.2fs, maxrel %.1e\n', nm, length(sn), nrow(tt), t1, t4, max(abs(a1-a4)/pmax(a1,1e-300))))
  tw = system.time(ns$simulate_walks(gw, tt, p$widthdt, gm.out=FALSE))[3]
  cat(sprintf('  simulate_walks total %.2fs\n', tw))
}
fx = lf('ccle_dm_local'); 
for (mc in c(1,4)) cat('oracle dm_local 1e5 mc.cores', mc, system.time(ref('forward_simulate')(fx$gw[[1]], fx$ft, 1e5, mc.cores=mc))[3], '\n')
for (mc in c(1,4)) cat('optim  dm_local 1e5 mc.cores', mc, system.time(ns$forward_simulate(fx$gw[[1]], fx$ft, 1e5, mc.cores=mc))[3], '\n')
