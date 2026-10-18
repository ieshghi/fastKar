## Benchmark the fastKar workhorse functions: reference fastKar (oracle, pinned 0a5497e)
## vs fastKarOptim (this worktree), on the fixtures in tests/testthat/fixtures.
##
##   Rscript bench/bench_workhorse.R [--funcs f1,f2] [--cases c1,c2] [--pkgs oracle,optim]
##                                   [--reps 3] [--tag label] [--optim PATH] [--out FILE]
##
## Per (pkg, case, bench) it records the median elapsed time over reps, the peak RSS increase
## (VmHWM reset through /proc/self/clear_refs before each run) and R's max heap use (gc reset).
## Rows are appended to bench/history_workhorse.tsv with the optim commit and the tag.
suppressMessages({ library(data.table); setDTthreads(1) })

args = commandArgs(trailingOnly = TRUE)
opt_arg = function(flag, default) {
    k = match(flag, args)
    if (is.na(k)) default else args[k + 1]
}
optim_path = normalizePath(opt_arg('--optim', '~/git/fastKar_optim'))
oracle_path = normalizePath(Sys.getenv('FASTKAR_ORACLE', '~/git/fastKar_oracle'))
reps = as.integer(opt_arg('--reps', '3'))
tag = opt_arg('--tag', '')
pkgs = strsplit(opt_arg('--pkgs', 'oracle,optim'), ',')[[1]]
out_file = opt_arg('--out', file.path(optim_path, 'bench', 'history_workhorse.tsv'))
fixture_dir = file.path(optim_path, 'tests', 'testthat', 'fixtures')

suppressMessages({
    if ('oracle' %in% pkgs) pkgload::load_all(oracle_path, export_all = FALSE, attach = FALSE, quiet = TRUE)
    if ('optim' %in% pkgs) pkgload::load_all(optim_path, export_all = FALSE, attach = FALSE, quiet = TRUE)
})
ns = c(oracle = 'fastKar', optim = 'fastKarOptim')
commit = tryCatch(system2('git', c('-C', optim_path, 'rev-parse', '--short', 'HEAD'), stdout = TRUE), error = function(e) NA)
dirty = length(tryCatch(system2('git', c('-C', optim_path, 'status', '--porcelain', '--untracked-files=no'), stdout = TRUE), error = function(e) character())) > 0

load_fixture = function(name) {
    fx = readRDS(file.path(fixture_dir, paste0(name, '.rds')))
    fx$gw = lapply(fx$kars, function(k) {
        gw = gGnome::gW(graph = fx$gg, snode.id = k$snode.id, circular = k$circular)
        if (!is.null(k$cn)) gw$set(cn = k$cn)
        gw
    })
    fx
}

rss_mb = function(field) {
    l = grep(paste0('^', field), readLines('/proc/self/status'), value = TRUE)
    as.numeric(gsub('[^0-9]', '', l)) / 1024
}
measure = function(expr_fn) {
    invisible(gc(reset = TRUE))
    try(writeLines('5', '/proc/self/clear_refs'), silent = TRUE)
    base = rss_mb('VmRSS')
    t = system.time(res <- expr_fn())[['elapsed']]
    g = gc()
    list(time = t, rss = rss_mb('VmHWM') - base, heap = sum(g[, ncol(g)]), res = res)
}

## Benchmarks: each takes (f = function getter for one package, fx = fixture, pix) and returns a
## closure that runs the timed work. Setup (noisy observations etc.) is done outside the closure
## with the same seed for both packages.
benches = list(
    forward_1e5 = list(pix = 1e5, run = function(f, fx, pix) function()
        f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10)),
    forward_1e4 = list(pix = 1e4, run = function(f, fx, pix) function()
        f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10)),
    candidates_nll_1e4 = list(pix = 1e4, run = function(f, fx, pix) {
        set.seed(1)
        obs = f('make_noisymap')(f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10), 1)[[1]]
        function() vapply(fx$gw, function(gw)
            f('compmaps')(obs, f('forward_simulate')(gw, target_region = fx$ft, pix.size = pix, depth = 10), ifsum = TRUE), 1)
    }),
    compmaps_1e4 = list(pix = 1e4, run = function(f, fx, pix) {
        set.seed(1)
        sim = f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10)
        obs = f('make_noisymap')(sim, 1)[[1]]
        function() f('compmaps')(obs, sim, ifsum = TRUE)
    }),
    noisydat_1e4 = list(pix = 1e4, run = function(f, fx, pix) {
        sim = f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10)
        function() { set.seed(1); f('make_noisydat')(sim, 10) }
    }),
    liktest_1e5 = list(pix = 1e5, run = function(f, fx, pix) {
        if (length(fx$gw) < 2) return(NULL)
        a = f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10)
        b = f('forward_simulate')(fx$gw[[2]], target_region = fx$ft, pix.size = pix, depth = 10)
        function() { set.seed(1); f('liktest_separable')(a, b, nsamp = 20) }
    }),
    get_dists_1e5 = list(pix = 1e5, run = function(f, fx, pix) {
        if (length(fx$gw) < 3) return(NULL)
        function() f('get_dists')(fx$gw[1:min(4, length(fx$gw))], target_region = fx$ft, pix.size = pix, depth = 10)
    }),
    wiring = list(pix = NA, run = function(f, fx, pix) function() f('gg.to.wiring')(fx$gg)),
    local_sampling = list(pix = NA, run = function(f, fx, pix) function() {
        set.seed(1); f('local.sampling')(fx$gg, nsteps = 4, nwalk = 20, onlyhash = TRUE)
    }),
    reads_1e5 = list(pix = 1e5, run = function(f, fx, pix) function() f('reads_fromwalk')(fx$gw[[1]], pix)),
    reads_1e7 = list(pix = 1e7, run = function(f, fx, pix) function() f('reads_fromwalk')(fx$gw[[1]], pix)),
    longread_kl_1e5 = list(pix = 1e5, run = function(f, fx, pix) {
        if (length(fx$gw) < 2) return(NULL)
        function() f('longread_kl')(fx$gw[[1]], fx$gw[[2]], readL = pix, depth = 10, use.nodes = NULL)
    }),
    compmaps_kl_1e4 = list(pix = 1e4, run = function(f, fx, pix) {
        a = f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10)
        b = if (length(fx$gw) > 1) f('forward_simulate')(fx$gw[[2]], target_region = fx$ft, pix.size = pix, depth = 10) else a
        function() f('compmaps')(a, b, return_kl = TRUE, area0 = pix^2, ifsum = TRUE)
    }),
    hic_cov_1e4 = list(pix = 1e4, run = function(f, fx, pix) {
        sim = f('forward_simulate')(fx$gw[[1]], target_region = fx$ft, pix.size = pix, depth = 10)
        function() f('hic_cov')(sim)
    }),
    liktest_lr_1e5 = list(pix = 1e5, run = function(f, fx, pix) {
        if (length(fx$gw) < 2) return(NULL)
        function() { set.seed(1); f('liktest_separable_lr')(fx$gw[[1]], fx$gw[[2]], readL = pix, depth = 10, nsamp = 20) }
    }),
    group_klcomp_1e5 = list(pix = 1e5, run = function(f, fx, pix) {
        if (length(fx$gw) < 3) return(NULL)
        ## gWalk$hash looks hash_karyotype_cpp up from the global env (packages are not attached here)
        assign('hash_karyotype_cpp', f('hash_karyotype_cpp'), envir = globalenv())
        function() f('group_klcomp')(fx$gw[1:min(4, length(fx$gw))], target_region = fx$ft, hic_resolution = pix, depth = 10)
    }),
    traverse = list(pix = NA, run = function(f, fx, pix) {
        w = f('gg.to.wiring')(fx$gg)
        function() f('traverse_graph_cpp')(w$internal.edges, w$loose.ends)
    }),
    sample_neighborhood = list(pix = NA, run = function(f, fx, pix) function() {
        set.seed(1); f('sample.neighborhood')(fx$gg, 20, return.gw = FALSE)
    }),
    edit_dist = list(pix = NA, run = function(f, fx, pix) {
        if (length(fx$gw) < 2) return(NULL)
        function() vapply(2:length(fx$gw), function(k) f('edit_dist_cpp')(fx$gw[[1]], fx$gw[[k]]), 1)
    })
)

default_cases = c('toy_bfb', 'ccle_pyrgo', 'ccle_dm', 'ccle_dm_local', 'ccle_chromothripsis', 'ccle_tyfonas', 'ccle_tyfonas_local')
cases = strsplit(opt_arg('--cases', paste(default_cases, collapse = ',')), ',')[[1]]
funcs = strsplit(opt_arg('--funcs', paste(names(benches), collapse = ',')), ',')[[1]]
stopifnot(all(funcs %in% names(benches)))

rows = list()
for (cs in cases) {
    fx = load_fixture(cs)
    for (bn in funcs) {
        b = benches[[bn]]
        for (pk in pkgs) {
            f = function(name) get(name, envir = asNamespace(ns[[pk]]))
            run = tryCatch(b$run(f, fx, b$pix), error = function(e) { message(cs, ' ', bn, ' ', pk, ' setup failed: ', conditionMessage(e)); NULL })
            if (is.null(run)) next
            ms = lapply(seq_len(reps), function(r) tryCatch(measure(run), error = function(e) { message(cs, ' ', bn, ' ', pk, ' failed: ', conditionMessage(e)); NULL }))
            ms = Filter(Negate(is.null), ms)
            if (!length(ms)) next
            row = data.table(date = format(Sys.time(), '%Y-%m-%d %H:%M'), commit = commit, dirty = dirty, tag = tag,
                             pkg = pk, case = cs, bench = bn, reps = length(ms),
                             time_s = median(vapply(ms, `[[`, 1, 'time')),
                             peak_rss_mb = max(vapply(ms, `[[`, 1, 'rss')),
                             r_heap_mb = max(vapply(ms, `[[`, 1, 'heap')))
            message(sprintf('%-26s %-20s %-7s %8.3f s  rss +%7.0f MB  heap %7.0f MB', cs, bn, pk, row$time_s, row$peak_rss_mb, row$r_heap_mb))
            rows[[length(rows) + 1]] = row
        }
    }
}
res = rbindlist(rows)
fwrite(res, out_file, sep = '\t', append = file.exists(out_file))
if (all(c('oracle', 'optim') %in% res$pkg)) {
    w = dcast(res, case + bench ~ pkg, value.var = c('time_s', 'peak_rss_mb'))
    w[, speedup := time_s_oracle / time_s_optim]
    w[, mem_ratio := (peak_rss_mb_oracle + 1) / (peak_rss_mb_optim + 1)]
    print(w[, .(case, bench, time_s_oracle, time_s_optim, speedup = round(speedup, 2),
                rss_oracle = round(peak_rss_mb_oracle), rss_optim = round(peak_rss_mb_optim), mem_ratio = round(mem_ratio, 2))])
}
