## The Python search end to end (needs a kinf environment: install_kinf()). Small problems, a few seconds each.

test_that('toy problems with a known truth are recovered exactly (incl. a circular walk)', {
    skip_if_no_loosefix()
    skip_if_not(kinf_available(), 'kinf Python not configured (install_kinf())')
    data.table::setDTthreads(1)
    for (nm in c('toy_repdup', 'ccle_pyrgo_local')) {
        x = if (nm == 'toy_repdup') kinf_toy(nm, pix = 5e4, depth = 30) else kinf_toy(nm, k = 'circ', pix = NULL, depth = 30)
        res = infer_karyotype(x$gg, x$gm, depth = x$depth, pix.size = x$pix, target_region = x$ft, chains = 4,
                              steps = 2000, block = 250, truth = x$truth, topk = 3, verbose = FALSE)
        expect_s3_class(res$gw, 'gWalk')
        expect_true(isTRUE(res$summary$exact_recovery), info = nm)
        expect_lte(res$summary$gap_to_truth, 1e-6)
        expect_equal(sum(res$gw$circular), sum(x$truth$circular))
        ## top-k: best first, distinct karyotypes, the first is the returned best
        expect_true(length(res$kars) >= 1 && length(res$kars) <= 3)
        expect_false(is.unsorted(res$nlls))
        expect_equal(res$nlls[1], res$summary$best_nll)
        expect_equal(lapply(res$kars[[1]]$snode.id, as.integer), lapply(res$gw$snode.id, as.integer))
    }
})

test_that('the Python truth NLL equals the R weighted NB NLL of the exported maps', {
    skip_if_no_loosefix()
    skip_if_not(kinf_available(), 'kinf Python not configured (install_kinf())')
    data.table::setDTthreads(1)
    for (nm in c('toy_repdup', 'ccle_pyrgo_local')) {
        x = if (nm == 'toy_repdup') kinf_toy(nm, pix = 5e4) else kinf_toy(nm, k = 'circ', pix = NULL)
        ex = kinf_export(x)
        wd = tempfile('kinf_nll_'); dir.create(wd)
        opt('.kinf_launch')(problem = ex$dir, workdir = wd, python = kinf_python(),
                            args = list(chains = 2, steps = 100, block = 50), truth = TRUE, verbose = FALSE)
        s = jsonlite::read_json(file.path(wd, 'out', 'summary.json'), simplifyVector = TRUE)
        expect_equal(s$nll_truth, kinf_ref_nll(ex$dir), tolerance = 1e-9, info = nm)
    }
})

test_that('rescale leaves the depth unchanged up to the wiring spread and reports the ratio', {
    skip_if_no_loosefix()
    skip_if_not(kinf_available(), 'kinf Python not configured (install_kinf())')
    data.table::setDTthreads(1)
    x = kinf_toy(pix = 5e4)
    res = infer_karyotype(x$gg, x$gm, depth = x$depth * 2, pix.size = x$pix, target_region = x$ft, chains = 2,
                          steps = 500, block = 250, rescale = 64, verbose = FALSE)
    expect_equal(res$summary$depth, res$summary$depth_ratio * x$depth * 2)
    expect_equal(res$summary$depth_ratio, 0.5, tolerance = 0.1)
})

test_that('the slurm backend writes a job script, submits it and collects the result (fake sbatch)', {
    skip_if_no_loosefix()
    skip_if_not(kinf_available(), 'kinf Python not configured (install_kinf())')
    data.table::setDTthreads(1)
    ## fake sbatch/sacct/squeue first on PATH: sbatch runs the job script here and records its arguments;
    ## sacct/squeue report nothing (so kinf_collect() decides from the result files)
    bin = tempfile('fakeslurm_'); dir.create(bin)
    writeLines(c('#!/bin/bash', 'echo "$@" > "$(dirname "$0")/sbatch_args"', 'out=""; prev=""',
                 'for a in "$@"; do [ "$prev" = "--output" ] && out="$a"; prev="$a"; done',
                 'script="${@: -1}"; bash "$script" > "${out//%j/4242}" 2>&1; echo 4242'), file.path(bin, 'sbatch'))
    writeLines(c('#!/bin/bash', 'exit 0'), file.path(bin, 'sacct'))
    writeLines(c('#!/bin/bash', 'exit 0'), file.path(bin, 'squeue'))
    Sys.chmod(file.path(bin, c('sbatch', 'sacct', 'squeue')), '755')
    old = Sys.getenv('PATH'); on.exit(Sys.setenv(PATH = old), add = TRUE)
    Sys.setenv(PATH = paste(bin, old, sep = ':'))
    x = kinf_toy(pix = 5e4, depth = 30)
    wd = tempfile('kinf_slurm_')
    gpu = identical(kinf_config()$build, 'cuda12')  # device 'auto' becomes --device gpu in a job
    ## the workdir under tempdir() also exercises the node-local warning
    expect_warning(job <- infer_karyotype(x$gg, x$gm, depth = x$depth, pix.size = x$pix, target_region = x$ft, chains = 4,
                          steps = 1000, block = 250, truth = x$truth, backend = 'slurm', wait = FALSE,
                          device = if (gpu) 'auto' else 'cpu', workdir = wd, slurm = list(partition = 'gpu', account = 'acct1',
                                                                       extra = '--qos=test'), verbose = FALSE), 'node-local')
    expect_s3_class(job, 'kinf_job')
    expect_equal(job$job_id, '4242')
    args = readLines(file.path(bin, 'sbatch_args'))
    expect_match(args, '--partition=gpu', fixed = TRUE)
    expect_match(args, '--gres=gpu:1', fixed = TRUE)
    expect_match(args, '--account=acct1', fixed = TRUE)
    expect_match(args, '--qos=test', fixed = TRUE)
    expect_false(grepl('--wait', args, fixed = TRUE))
    sh = readLines(file.path(wd, 'kinf_job.sh'))
    expect_true(any(grepl('^unset PYTHONPATH', sh)))
    expect_true(any(grepl(paste0("'--device' '", if (gpu) 'gpu' else 'cpu', "'"), sh, fixed = TRUE)))
    expect_true(file.exists(file.path(wd, 'slurm-4242.out')))
    res = kinf_collect(job)
    expect_s3_class(res$gw, 'gWalk')
    expect_true(isTRUE(res$summary$exact_recovery))
})

test_that('history_every records every chain and agrees with the per-block trace', {
    skip_if_no_loosefix()
    skip_if_not(kinf_available(), 'kinf Python not configured (install_kinf())')
    data.table::setDTthreads(1)
    x = kinf_toy(pix = 5e4)
    res = infer_karyotype(x$gg, x$gm, depth = x$depth, pix.size = x$pix, target_region = x$ft, chains = 3,
                          steps = 500, block = 250, truth = x$truth, history_every = 25, verbose = FALSE)
    h = res$history
    expect_equal(names(h), c('steps', 'chain', 'T', 'nll'))
    expect_equal(sort(unique(h$steps)), seq(0, 500, by = 25))
    expect_true(all(h[, .N, by = steps]$N == 3))
    ## at block ends the chains' NLLs are the trace's per-rung NLLs
    for (b in res$trace$steps) {
        tv = sort(unlist(res$trace[steps == b, -(1:2)]))
        expect_equal(sort(h[steps == b]$nll), unname(tv), tolerance = 1e-5)
    }
})

test_that('init warm-starts chains from a previous (coarser) run or a gWalk, and rejects foreign wirings', {
    skip_if_no_loosefix()
    skip_if_not(kinf_available(), 'kinf Python not configured (install_kinf())')
    data.table::setDTthreads(1)
    x = kinf_toy(pix = 1e5, depth = 30)
    coarse = infer_karyotype(x$gg, x$gm, depth = x$depth, pix.size = 1e5, target_region = x$ft, chains = 4,
                             steps = 2000, block = 250, truth = x$truth, verbose = FALSE)
    expect_true(isTRUE(coarse$summary$exact_recovery))
    ## fine pass on the same graph + target: chains 1-3 start from the coarse best, chain 4 random
    xf = kinf_toy(pix = 5e4, depth = 30)
    fine = infer_karyotype(xf$gg, xf$gm, depth = xf$depth, pix.size = 5e4, target_region = xf$ft, chains = 4,
                           steps = 250, block = 250, method = 'greedy', truth = xf$truth, init = coarse,
                           init_random = 1, history_every = 250, verbose = FALSE)
    h0 = fine$history[steps == 0][order(chain)]
    expect_equal(h0$nll[1:3], rep(0, 3), tolerance = 1e-6)
    expect_gt(h0$nll[4], 1)
    expect_true(isTRUE(fine$summary$exact_recovery))
    ## a gWalk on gg as init
    g = infer_karyotype(xf$gg, xf$gm, depth = xf$depth, pix.size = 5e4, target_region = xf$ft, chains = 2,
                        steps = 250, block = 250, method = 'greedy', truth = xf$truth, init = xf$truth,
                        history_every = 250, verbose = FALSE)
    expect_equal(g$history[steps == 0]$nll, c(0, 0), tolerance = 1e-6)
    ## a wiring of another graph fails with a clear message
    other = kinf_toy('ccle_pyrgo_local', k = 'circ', pix = NULL)
    expect_error(infer_karyotype(xf$gg, xf$gm, depth = xf$depth, pix.size = 5e4, target_region = xf$ft,
                                 chains = 2, steps = 250, block = 250, init = other$truth, verbose = FALSE),
                 'do not fit this problem')
})

test_that('Python and R agree when the footprint covers only part of the graph (outside nodes wired, not tiled)', {
    skip_if_no_loosefix()
    skip_if_not(kinf_available(), 'kinf Python not configured (install_kinf())')
    data.table::setDTthreads(1)
    ## whole-graph fixtures whose footprint cuts through nodes (ccle_dm: chr2 node 16.2-53.1 Mb covered at both ends)
    for (nm in c('ccle_pyrgo', 'ccle_dm')) {
        x = kinf_toy(nm, pix = NULL, depth = 30)
        ex = kinf_export(x)
        expect_gt(length(x$gg$nodes), length(x$gg$nodes %&% x$ft))
        wd = tempfile('kinf_sub_'); dir.create(wd)
        opt('.kinf_launch')(problem = ex$dir, workdir = wd, python = kinf_python(),
                            args = list(chains = 2, steps = 100, block = 50), truth = TRUE, verbose = FALSE)
        s = jsonlite::read_json(file.path(wd, 'out', 'summary.json'), simplifyVector = TRUE)
        expect_equal(s$nll_truth, kinf_ref_nll(ex$dir), tolerance = 1e-9, info = nm)
    }
})
