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
