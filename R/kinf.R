# kinf.R: R front end to the Python karyotype search shipped in inst/python/kinf (JAX; exact-NLL parallel
# tempering over the wirings of a genome graph). Setup once with install_kinf(); then
#   res = infer_karyotype(gg, hic = 'sample.mcool', pix.size = 1e6)
# runs prep_data_for_py() -> infer.py -> read_kinf_gwalks(), locally or as a Slurm job.

.kinf_uv_version = '0.12.24'
.kinf_cache = new.env(parent = emptyenv())

#' Folder of the shipped Python code (infer.py and modules); works installed and under devtools::load_all().
kinf_dir = function(){
  d = system.file('python', 'kinf', package = 'fastKar')
  if (!nzchar(d)) stop('inst/python/kinf not found in the fastKar installation')
  d
}

.kinf_config_path = function() file.path(tools::R_user_dir('fastKar', 'config'), 'kinf.json')

#' The kinf configuration written by install_kinf() (python, envdir, build, standalone, slurm defaults), or an
#' empty list.
kinf_config = function(){
  f = .kinf_config_path()
  if (!file.exists(f)) return(list())
  jsonlite::read_json(f, simplifyVector = TRUE)
}

.kinf_write_config = function(cfg){
  f = .kinf_config_path()
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, null = 'null'), f)
  invisible(f)
}

#' Python interpreter for kinf. Lookup order: argument > options(fastKar.kinf.python) > env var
#' FASTKAR_KINF_PYTHON > config file (install_kinf()).
kinf_python = function(python = NULL, mustWork = TRUE){
  cand = c(python, getOption('fastKar.kinf.python'), Sys.getenv('FASTKAR_KINF_PYTHON'), kinf_config()$python)
  cand = cand[!is.na(cand) & nzchar(cand)]
  if (!length(cand)){
    if (mustWork) stop('no Python configured for kinf: run install_kinf() (or set options(fastKar.kinf.python = ...))')
    return(NA_character_)
  }
  if (mustWork && !file.exists(cand[1])) stop('kinf python not found: ', cand[1], ' (rerun install_kinf()?)')
  cand[1]
}

# Environment settings for every kinf Python subprocess: drop PYTHONPATH/PYTHONHOME (the easybuild module
# PYTHONPATH breaks numpy), and LD_LIBRARY_PATH too when the interpreter is the standalone uv Python (it needs
# none; stray module libraries can shadow the CUDA wheels).
.kinf_envs = function(python, extra = NULL){
  cfg = kinf_config()
  standalone = isTRUE(cfg$standalone) && identical(normalizePath(python, mustWork = FALSE),
                                                   normalizePath(cfg$python %||% '', mustWork = FALSE))
  c('-u', 'PYTHONPATH', '-u', 'PYTHONHOME', if (standalone) c('-u', 'LD_LIBRARY_PATH'),
    'XLA_PYTHON_CLIENT_PREALLOCATE=false', extra)
}

`%||%` = function(a, b) if (is.null(a)) b else a

#' TRUE when a kinf Python is configured and imports jax, numpy and pandas (cached per interpreter).
kinf_available = function(python = NULL){
  py = kinf_python(python, mustWork = FALSE)
  if (is.na(py) || !file.exists(py)) return(FALSE)
  if (is.null(.kinf_cache[[py]])){
    out = suppressWarnings(system2('env', c(.kinf_envs(py), shQuote(py), '-c', shQuote('import jax, numpy, pandas')),
                                   stdout = FALSE, stderr = FALSE))
    .kinf_cache[[py]] = identical(as.integer(out), 0L)
  }
  .kinf_cache[[py]]
}

# Default Slurm settings for kinf jobs: the built-in defaults, overridden by the config file, overridden by slurm.
.kinf_slurm = function(slurm = list()){
  def = list(partition = 'gpu', gres = 'gpu:1', time = '2:00:00', mem = '32G', cpus = 8, account = NULL,
             extra = NULL)
  cfg = kinf_config()$slurm
  if (length(cfg)) def = utils::modifyList(def, as.list(cfg))
  utils::modifyList(def, as.list(slurm))
}

.kinf_fstype = function(path){
  d = path
  while (!dir.exists(d) && d != dirname(d)) d = dirname(d)
  out = tryCatch(suppressWarnings(system2('stat', c('-f', '-c', '%T', shQuote(d)), stdout = TRUE, stderr = FALSE)),
                 error = function(e) NA_character_)
  if (length(out)) out[1] else NA_character_
}

#' Set up the Python environment for karyotype inference (one time).
#'
#' Builds a dedicated virtual environment with the pinned JAX/numpy/pandas versions of
#' inst/python/kinf/requirements-<build>.txt and records it in the kinf config file
#' (tools::R_user_dir('fastKar', 'config')/kinf.json), then runs kinf_check(). By default the interpreter is a
#' standalone Python fetched by uv (relocatable; runs without environment modules, e.g. in Slurm jobs); uv itself
#' is downloaded into envdir/bin when it isn't on PATH. Nothing is installed into the system Python or R libraries.
#'
#' @param envdir folder for the environment (~1 GB cpu, ~3-4 GB cuda12). It must be visible from the compute
#'   nodes that will run the search (a shared filesystem, not /tmp). Default tools::R_user_dir('fastKar',
#'   'data')/kinf.
#' @param build 'cuda12' (NVIDIA GPUs; needs only the driver, not a system CUDA), 'cpu', or 'auto' (cuda12 on
#'   Linux x86_64 when nvidia-smi works here, else cpu; pick cuda12 explicitly when installing from a login node
#'   for GPU compute nodes).
#' @param python with uv = FALSE: an existing Python >= 3.11 to build the venv from (default: python3 on PATH).
#' @param uv use a uv-managed standalone Python (recommended).
#' @param test_gpu_job also run the self test as a short Slurm GPU job (uses your allocation).
#' @param slurm Slurm defaults to store in the config (partition, gres, time, mem, cpus, account, extra).
#' @param force rebuild the environment from scratch.
#' @return the config list, invisibly.
install_kinf = function(envdir = NULL, build = c('auto', 'cuda12', 'cpu'), python = NULL, uv = TRUE,
                        test_gpu_job = FALSE, slurm = list(), force = FALSE){
  build = match.arg(build)
  if (is.null(envdir)) envdir = file.path(tools::R_user_dir('fastKar', 'data'), 'kinf')
  dir.create(envdir, recursive = TRUE, showWarnings = FALSE)
  envdir = normalizePath(envdir)
  fs = .kinf_fstype(envdir)
  if (grepl('^(/tmp|/var/tmp|/dev/shm)', envdir) || fs %in% c('tmpfs', 'ramfs'))
    warning('envdir ', envdir, ' looks node-local (', fs, '); Slurm jobs on other nodes will not see it')
  sys = Sys.info()
  if (build == 'auto'){
    gpu_here = nzchar(Sys.which('nvidia-smi')) &&
      identical(suppressWarnings(system2('nvidia-smi', '-L', stdout = FALSE, stderr = FALSE)), 0L)
    build = if (sys[['sysname']] == 'Linux' && sys[['machine']] == 'x86_64' && gpu_here) 'cuda12' else 'cpu'
    message('build = ', build, if (build == 'cpu') ' (no GPU seen here; use build = "cuda12" for GPU compute nodes)')
  }
  if (build == 'cuda12' && !(sys[['sysname']] == 'Linux' && sys[['machine']] == 'x86_64'))
    stop('build = "cuda12" needs Linux x86_64')
  req = file.path(kinf_dir(), paste0('requirements-', build, '.txt'))
  venv = file.path(envdir, 'venv')
  vpy = file.path(venv, 'bin', 'python')
  if (force && dir.exists(venv)) unlink(venv, recursive = TRUE)
  run = function(cmd, args, envs = NULL){
    message('$ ', cmd, ' ', paste(args, collapse = ' '))
    st = system2('env', c('-u', 'PYTHONPATH', '-u', 'PYTHONHOME', envs, shQuote(cmd), shQuote(args)))
    if (!identical(as.integer(st), 0L)) stop('command failed (exit ', st, '): ', cmd, ' ', paste(args, collapse = ' '))
  }

  if (uv){
    uvbin = .kinf_get_uv(envdir)
    uvenv = c('-u', 'VIRTUAL_ENV', '-u', 'CONDA_PREFIX',  # env: unsets before assignments
              paste0('UV_PYTHON_INSTALL_DIR=', file.path(envdir, 'python')),
              paste0('UV_CACHE_DIR=', file.path(envdir, 'uv-cache')), 'UV_PYTHON_PREFERENCE=only-managed')
    run(uvbin, c('python', 'install', '3.11'), uvenv)
    if (!file.exists(vpy)) run(uvbin, c('venv', '--python', '3.11', venv), uvenv)
    run(uvbin, c('pip', 'install', '--python', vpy, '-r', req), uvenv)
    run(uvbin, c('cache', 'clean'), uvenv)
  } else {
    if (is.null(python)) python = Sys.which('python3')
    if (!nzchar(python)) stop('no python3 found; give python = or use uv = TRUE')
    ver = system2('env', c('-u', 'PYTHONPATH', shQuote(python), '-c',
                           shQuote('import sys; print(sys.version_info >= (3, 11))')), stdout = TRUE)
    if (!identical(ver, 'True')) stop(python, ' is older than Python 3.11')
    if (!file.exists(vpy)) run(python, c('-m', 'venv', venv))
    run(vpy, c('-m', 'pip', 'install', '--upgrade', 'pip'))
    run(vpy, c('-m', 'pip', 'install', '-r', req))
    # a venv on a module-system Python needs its libpython from LD_LIBRARY_PATH; jobs without the module fail
    bare = suppressWarnings(system2('env', c('-i', shQuote(vpy), '-c', shQuote('import sys')), stdout = FALSE,
                                    stderr = FALSE))
    if (!identical(as.integer(bare), 0L))
      warning(vpy, ' does not start in an empty environment (it needs a module or LD_LIBRARY_PATH): Slurm jobs ',
              'must load the same module. uv = TRUE avoids this.')
  }

  cfg = kinf_config()
  cfg$python = vpy; cfg$envdir = envdir; cfg$build = build; cfg$standalone = uv
  cfg$slurm = utils::modifyList(as.list(cfg$slurm), as.list(slurm))
  if (!length(cfg$slurm)) cfg$slurm = NULL
  f = .kinf_write_config(cfg)
  message('kinf config written to ', f)
  if (exists(vpy, envir = .kinf_cache)) rm(list = vpy, envir = .kinf_cache)
  kinf_check(python = vpy)
  if (test_gpu_job) kinf_check(python = vpy, backend = 'slurm')
  invisible(cfg)
}

# uv binary: from PATH, or the pinned release downloaded into envdir/bin (sha256 checked).
.kinf_get_uv = function(envdir){
  local = file.path(envdir, 'bin', 'uv')
  if (file.exists(local)) return(local)
  onpath = Sys.which('uv')
  if (nzchar(onpath)) return(unname(onpath))
  sys = Sys.info()
  arch = c(x86_64 = 'x86_64', aarch64 = 'aarch64', arm64 = 'aarch64')[sys[['machine']]]
  os = c(Linux = 'unknown-linux-gnu', Darwin = 'apple-darwin')[sys[['sysname']]]
  if (is.na(arch) || is.na(os)) stop('no uv build for ', sys[['sysname']], '/', sys[['machine']], '; use uv = FALSE')
  name = paste0('uv-', arch, '-', os)
  url = sprintf('https://github.com/astral-sh/uv/releases/download/%s/%s.tar.gz', .kinf_uv_version, name)
  tmp = tempfile('uv_'); dir.create(tmp)
  tgz = file.path(tmp, 'uv.tar.gz'); shaf = file.path(tmp, 'uv.sha256')
  message('downloading uv ', .kinf_uv_version, ' from ', url)
  utils::download.file(url, tgz, mode = 'wb', quiet = TRUE)
  utils::download.file(paste0(url, '.sha256'), shaf, quiet = TRUE)
  want = strsplit(readLines(shaf, warn = FALSE)[1], '[[:space:]]+')[[1]][1]
  got = strsplit(system2(if (nzchar(Sys.which('sha256sum'))) 'sha256sum' else 'shasum',
                         c(if (!nzchar(Sys.which('sha256sum'))) c('-a', '256'), shQuote(tgz)), stdout = TRUE),
                 '[[:space:]]+')[[1]][1]
  if (!identical(tolower(want), tolower(got))) stop('uv download checksum mismatch')
  utils::untar(tgz, exdir = tmp)
  dir.create(dirname(local), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(tmp, name, 'uv'), local, overwrite = TRUE)
  Sys.chmod(local, '755')
  unlink(tmp, recursive = TRUE)
  local
}

#' Check the kinf setup: JAX devices, then a short search on the shipped self-test problem (a small simulated
#' graph with a known truth) that must reach the truth NLL.
#'
#' @param python interpreter (default: kinf_python()).
#' @param backend 'local' or 'slurm' (submits the test as a short GPU job with the Slurm defaults).
#' @param device passed to infer.py ('auto', 'gpu', 'cpu').
#' @param slurm Slurm overrides for backend = 'slurm'.
#' @return list(ok, device, summary), invisibly.
kinf_check = function(python = NULL, backend = c('local', 'slurm'), device = c('auto', 'gpu', 'cpu'),
                      slurm = list(), gpu = NULL, cores = NULL){
  backend = match.arg(backend); device = match.arg(device)
  py = kinf_python(python)
  if (backend == 'local'){
    # minimal environment (as a batch job might see it): the interpreter must start and find jax on its own
    probe = "import jax; jax.config.update('jax_enable_x64', True); import numpy, pandas; print(jax.__version__, jax.devices())"
    out = suppressWarnings(system2('env', c('-i', paste0('HOME=', Sys.getenv('HOME')), 'PATH=/usr/bin:/bin',
                                            'XLA_PYTHON_CLIENT_PREALLOCATE=false',
                                            if (!is.null(gpu)) paste0('CUDA_VISIBLE_DEVICES=', gpu),
                                            shQuote(py), '-c', shQuote(probe)), stdout = TRUE, stderr = TRUE))
    st = attr(out, 'status')
    if (!is.null(st) && st != 0) stop('kinf python failed to import jax in a clean environment:\n',
                                      paste(utils::tail(out, 15), collapse = '\n'))
    message('jax ', utils::tail(out, 1))
  }
  workdir = if (backend == 'slurm') .kinf_runs_dir() else tempfile('kinf_check_')
  if (backend == 'slurm') workdir = tempfile('kinf_check_', tmpdir = workdir)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  job = .kinf_launch(problem = file.path(kinf_dir(), 'selftest'), workdir = workdir, python = py,
                     args = list(method = 'pt', chains = 4, steps = 2000, block = 250, tmin = 1, tmax = 1000,
                                 seed = 0, device = if (backend == 'slurm') 'gpu' else device),
                     truth = TRUE, backend = backend, slurm = slurm, wait = TRUE, gpu = gpu, cores = cores,
                     verbose = FALSE)
  summ = jsonlite::read_json(file.path(job$workdir, 'out', 'summary.json'), simplifyVector = TRUE)
  ok = isTRUE(summ$exact_recovery) || abs(summ$gap_to_truth) < 1e-6
  message(sprintf('kinf self test (%s, device %s): %s; gap to truth NLL %.3g, exact recovery %s, %.1f s',
                  backend, summ$device, if (ok) 'OK' else 'FAILED', summ$gap_to_truth, summ$exact_recovery,
                  summ$wall_s))
  if (!ok) warning('kinf self test did not reach the truth; see ', job$workdir)
  invisible(list(ok = ok, device = summ$device, summary = summ, workdir = job$workdir))
}

.kinf_runs_dir = function(){
  cfg = kinf_config()
  d = cfg$workroot %||% (if (!is.null(cfg$envdir)) file.path(cfg$envdir, 'runs'))
  if (is.null(d)) stop('no shared run folder: give workdir = (on a filesystem the compute nodes see)')
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Run (local) or submit (slurm) infer.py on an exported problem folder; output goes to workdir/out.
.kinf_launch = function(problem, workdir, python, args, truth = FALSE, backend = 'local', slurm = list(),
                        wait = TRUE, gpu = NULL, cores = NULL, verbose = TRUE){
  odir = file.path(workdir, 'out')
  cmd = c(python, file.path(kinf_dir(), 'infer.py'), problem, odir)
  for (nm in names(args)) cmd = c(cmd, paste0('--', nm), format(args[[nm]], scientific = FALSE))
  if (truth) cmd = c(cmd, '--truth')
  if (backend == 'local'){
    if (!is.null(cores)) cmd = c('taskset', '-c', paste(cores, collapse = ','), cmd)
    envs = .kinf_envs(python, if (!is.null(gpu)) paste0('CUDA_VISIBLE_DEVICES=', paste(gpu, collapse = ',')))
    if (verbose) message('running: ', paste(cmd, collapse = ' '))
    log = file.path(workdir, 'infer.log')
    status = system2('env', c(envs, shQuote(cmd)), stdout = if (verbose) '' else log,
                     stderr = if (verbose) '' else log)
    if (!identical(as.integer(status), 0L))
      stop('infer.py failed (exit ', status, '); problem folder kept in ', problem,
           if (!verbose) paste0('\n', paste(utils::tail(readLines(log), 20), collapse = '\n')))
    return(structure(list(job_id = NA, workdir = workdir, backend = 'local'), class = 'kinf_job'))
  }
  if (!nzchar(Sys.which('sbatch'))) stop('sbatch not found: backend = "slurm" needs a Slurm submit host')
  s = .kinf_slurm(slurm)
  envs = .kinf_envs(python)
  unsets = envs[which(envs == '-u') + 1]
  script = file.path(workdir, 'kinf_job.sh')
  writeLines(c('#!/bin/bash', 'set -eo pipefail', paste('unset', paste(unsets, collapse = ' ')),
               'export XLA_PYTHON_CLIENT_PREALLOCATE=false',
               'echo "host $(hostname) CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"',
               paste(shQuote(cmd), collapse = ' ')), script)
  Sys.chmod(script, '755')
  sargs = c('--parsable', '--job-name', 'kinf', '--output', file.path(workdir, 'slurm-%j.out'),
            if (!is.null(s$partition)) paste0('--partition=', s$partition),
            if (!is.null(s$gres)) paste0('--gres=', s$gres),
            if (!is.null(s$time)) paste0('--time=', s$time),
            if (!is.null(s$mem)) paste0('--mem=', s$mem),
            if (!is.null(s$cpus)) paste0('--cpus-per-task=', s$cpus),
            if (!is.null(s$account)) paste0('--account=', s$account),
            s$extra, if (wait) '--wait', script)
  if (verbose) message('sbatch ', paste(sargs, collapse = ' '))
  out = suppressWarnings(system2('sbatch', shQuote(sargs), stdout = TRUE, stderr = TRUE))
  id = utils::tail(grep('^[0-9]+', out, value = TRUE), 1)
  id = sub(';.*', '', id)
  job = structure(list(job_id = if (length(id)) id else NA, workdir = workdir, backend = 'slurm'),
                  class = 'kinf_job')
  st = attr(out, 'status')
  if (!length(id)) stop('sbatch failed:\n', paste(out, collapse = '\n'))
  if (verbose) message('submitted Slurm job ', id)
  if (wait && !is.null(st) && st != 0)
    stop('Slurm job ', id, ' failed (sbatch exit ', st, '):\n', .kinf_logtail(job))
  if (wait && !file.exists(file.path(odir, 'summary.json')))
    stop('Slurm job ', id, ' finished without results:\n', .kinf_logtail(job))
  job
}

.kinf_logtail = function(job, n = 25){
  f = file.path(job$workdir, paste0('slurm-', job$job_id, '.out'))
  if (!file.exists(f)) return('(no slurm log)')
  paste(utils::tail(readLines(f, warn = FALSE), n), collapse = '\n')
}

# Stop unless gg can be exported: copy numbers present and integer, and on every node side the junction copies
# are at most the node copies (loose ends are fine: they become walk ends).
.kinf_check_graph = function(gg){
  nd = gg$nodes$dt; ed = gg$edges$dt
  if (is.null(nd$cn) || anyNA(nd$cn) || (nrow(ed) && (is.null(ed$cn) || anyNA(ed$cn))))
    stop('gg needs node and junction copy numbers (cn) without NAs; run loosefix(gg) / balance() first')
  if (any(nd$cn != round(nd$cn)) || (nrow(ed) && any(ed$cn != round(ed$cn)))) stop('gg copy numbers must be integers')
  if (nrow(ed)){
    inc = rbind(ed[, .(node.id = n1, side = n1.side, cn)], ed[, .(node.id = n2, side = n2.side, cn)])
    inc = inc[, .(ecn = sum(cn)), by = .(node.id, side)]
    inc[, ncn := nd$cn[match(node.id, nd$node.id)]]
    bad = inc[ecn > ncn]
    if (nrow(bad)) stop('gg is not balanced: on ', nrow(bad), ' node sides (nodes ',
                        paste(utils::head(unique(bad$node.id), 10), collapse = ','),
                        ') the junction copies exceed the node copies; fix the graph (e.g. loosefix(), balance())')
  }
  invisible(TRUE)
}

#' Infer a karyotype (wiring of gg) from Hi-C by exact-NLL search (parallel tempering by default).
#'
#' Exports gg and the map with prep_data_for_py(), runs inst/python/kinf/infer.py (locally or as a Slurm job)
#' and reads the best wiring back as a gWalk on gg. Set up the Python side once with install_kinf().
#'
#' @param gg gGraph with node copy numbers (gg$nodes$dt$cn) and junction copy numbers (gg$edges$dt$cn).
#' @param hic gMatrix of observed counts, or the path to an .mcool (read at the finest resolution <= pix.size).
#' @param depth Hi-C depth (fastKar convention). Required for a gMatrix; for an .mcool, estimated with
#'   estimate.depthratio(hic, ploidy = ploidy) when NULL.
#' @param pix.size tile size (bp) of the fit.
#' @param target_region GRanges to fit; default: the whole graph (streduce(gg$nodes$gr)).
#' @param ploidy used for the depth estimate from an .mcool.
#' @param purity only 1 is supported.
#' @param method 'pt' (parallel tempering) or 'greedy' (all chains at T = 0).
#' @param chains,steps,block,tmin,tmax,seed search settings: chains (= PT rungs), proposals per chain, proposals
#'   between temperature swaps, geometric ladder tmin..tmax, random seed.
#' @param device 'auto' (GPU if available, else CPU with a warning), 'gpu' (fail without one), 'cpu'.
#' @param gpu local backend: GPU index (sets CUDA_VISIBLE_DEVICES); NULL leaves it unset.
#' @param cores local backend: CPU cores to pin the Python process to, e.g. 0:7 (each JAX process otherwise
#'   starts ~1000 threads, and the per-user thread limit can make parallel runs fail). NULL: no pinning.
#' @param backend 'local' (run here) or 'slurm' (submit a GPU job with sbatch).
#' @param slurm Slurm settings overriding the config defaults: partition, gres, time, mem, cpus, account, extra
#'   (character vector of further sbatch arguments).
#' @param wait slurm backend: wait for the job (sbatch --wait) and return the result; FALSE returns a kinf_job
#'   handle for kinf_collect().
#' @param workdir folder for the exported problem (workdir/problem) and the search output (workdir/out). Default
#'   a temp folder (local) or a new folder under the config's shared run folder (slurm).
#' @param truth optional gWalk on gg (simulations): also scores it (gap to its NLL, exact recovery).
#' @param topk also return the topk best distinct karyotypes visited by any chain (as bestfit_search(topk =)):
#'   each chain keeps its 2*topk lowest-NLL states (NLLs equal to 1e-9 relative count once), merged over chains
#'   and deduplicated as karyotypes.
#' @param rescale as bestfit_search(rescale = TRUE): multiply depth by observed total / median expected total
#'   over random wirings (the analogue of its sampled candidates) before the search. TRUE = 256 random wirings,
#'   or give the number. The ratio and the depth used are in summary (depth_ratio, depth).
#' @param history_every record the NLL of every chain every history_every proposals (and at step 0) for
#'   convergence plots; returned as history (data.table steps, chain, T = the chain's temperature then, nll; minus
#'   the truth NLL when truth is given). 0: off. The best-so-far curve is cummin over chains.
#' @param python kinf interpreter (default kinf_python()).
#' @return list(gw = gWalk on gg, summary = list from summary.json, trace = data.table of the search curve
#'   (best and per-rung NLL per block), workdir); with topk > 1 also kars (list of gWalks, best first) and nlls;
#'   with history_every > 0 also history;
#'   with backend = 'slurm', wait = FALSE a kinf_job handle.
infer_karyotype = function(gg, hic, depth = NULL, pix.size = 1e6, target_region = NULL, ploidy = 2, purity = 1,
                           method = c('pt', 'greedy'), chains = 16, steps = 20000, block = 500, tmin = 1,
                           tmax = 1000, seed = 0, device = c('auto', 'gpu', 'cpu'), gpu = NULL, cores = NULL,
                           backend = c('local', 'slurm'), slurm = list(), wait = TRUE, workdir = NULL,
                           truth = NULL, topk = 1, rescale = FALSE, history_every = 0, python = NULL,
                           verbose = TRUE){
  method = match.arg(method); device = match.arg(device); backend = match.arg(backend)
  if (purity != 1) stop('purity != 1 is not supported by the Python model yet')
  py = kinf_python(python)
  .kinf_check_graph(gg)
  if (!is.null(truth) && !inherits(truth, 'gWalk')) stop('truth must be a gWalk on gg')
  if (is.null(target_region)) target_region = streduce(gg$nodes$gr)
  if (is.null(workdir))
    workdir = if (backend == 'slurm') tempfile('kinf_', tmpdir = .kinf_runs_dir()) else tempfile('kinf_')
  if (backend == 'slurm' && grepl('^(/tmp|/var/tmp|/dev/shm)', workdir))
    warning('workdir ', workdir, ' is probably node-local; the Slurm job may not see it')

  if (is.character(hic)){
    if (tools::file_ext(hic) != 'mcool') stop('hic: give a gMatrix or an .mcool path')
    avail = as.integer(names(rhdf5::h5read(hic, '/resolutions')))
    if (!any(avail <= pix.size)) stop('no mcool resolution <= pix.size (available: ', paste(sort(avail), collapse = ', '), ')')
    res = max(avail[avail <= pix.size])
    if (verbose) message('reading ', hic, ' at ', res, ' bp')
    if (is.null(depth)) depth = estimate.depthratio(hic, ploidy = ploidy)
    hic = cooler(hic, res = res, gr = target_region)
  } else if (is.null(depth)) stop('depth is required when hic is a gMatrix')
  if (verbose) message('depth ', signif(depth, 4))

  cpu_likely = device == 'cpu' || (backend == 'local' && device == 'auto' && identical(kinf_config()$build, 'cpu'))
  if (cpu_likely && chains * steps >= 2e5)
    message('note: ', chains, ' x ', steps, ' proposals on CPU; expect roughly 20x the GPU time')

  workdir = normalizePath(workdir, mustWork = FALSE)
  pdir = file.path(workdir, 'problem')
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  prep_data_for_py(gg = gg, gm = hic, depth = depth, purity = purity, pix.size = pix.size,
                   target_region = target_region, write.to.folder = paste0(pdir, '/'), truth = truth)
  job = .kinf_launch(problem = pdir, workdir = workdir, python = py,
                     args = list(method = method, chains = chains, steps = steps, block = block, tmin = tmin,
                                 tmax = tmax, seed = seed, topk = topk, history = as.integer(history_every),
                                 rescale = if (isTRUE(rescale)) 256L else as.integer(rescale),
                                 device = if (backend == 'slurm' && device == 'auto') 'gpu' else device),
                     truth = !is.null(truth), backend = backend, slurm = slurm, wait = wait, gpu = gpu,
                     cores = cores, verbose = verbose)
  job$graph = gg
  if (backend == 'slurm' && !wait) return(job)
  .kinf_read_results(job)
}

.kinf_read_results = function(job){
  odir = file.path(job$workdir, 'out')
  out = list(gw = read_kinf_gwalks(file.path(odir, 'gwalks.json'), job$graph),
             summary = jsonlite::read_json(file.path(odir, 'summary.json'), simplifyVector = TRUE),
             trace = data.table::fread(file.path(odir, 'trace.tsv')),
             workdir = job$workdir)
  if (file.exists(file.path(odir, 'history.tsv'))) out$history = data.table::fread(file.path(odir, 'history.tsv'))
  if (file.exists(file.path(odir, 'topk.json'))) out = c(out, read_kinf_topk(file.path(odir, 'topk.json'), job$graph))
  out
}

#' Read a topk.json written by infer.py --topk: list(kars = gWalks on graph (or plain walk lists when graph is
#' NULL), best first, nlls).
read_kinf_topk = function(path, graph = NULL){
  x = jsonlite::read_json(path, simplifyVector = FALSE)
  kars = lapply(x, function(k){
    snode.id = lapply(k$snode.ids, function(w) as.integer(unlist(w)))
    circular = as.logical(unlist(k$circular))
    if (is.null(graph)) list(snode.id = snode.id, circular = circular)
    else gW(graph = graph, snode.id = snode.id, circular = circular)
  })
  list(kars = kars, nlls = vapply(x, function(k) as.numeric(k$nll), numeric(1)))
}

#' Collect a Slurm kinf job submitted with infer_karyotype(..., backend = 'slurm', wait = FALSE).
#'
#' @param job the kinf_job handle.
#' @return the infer_karyotype() result when the job completed; otherwise reports the state and returns the
#'   handle invisibly (pending/running) or stops with the tail of the Slurm log (failed).
kinf_collect = function(job){
  if (!inherits(job, 'kinf_job')) stop('job must be a kinf_job handle from infer_karyotype(wait = FALSE)')
  done = file.exists(file.path(job$workdir, 'out', 'summary.json'))
  state = NA_character_
  if (!is.na(job$job_id) && nzchar(Sys.which('sacct'))){
    st = suppressWarnings(system2('sacct', c('-j', job$job_id, '-X', '-n', '-P', '-o', 'State'), stdout = TRUE,
                                  stderr = FALSE))
    if (length(st) && nzchar(st[1])) state = trimws(sub(' .*', '', st[1]))
  }
  if (is.na(state) && !is.na(job$job_id) && nzchar(Sys.which('squeue'))){
    st = suppressWarnings(system2('squeue', c('-j', job$job_id, '-h', '-o', '%T'), stdout = TRUE, stderr = FALSE))
    if (length(st) && nzchar(st[1])) state = trimws(st[1])
  }
  if (state %in% c('PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'REQUEUED', 'SUSPENDED')){
    message('job ', job$job_id, ' is ', state)
    return(invisible(job))
  }
  if (done && (is.na(state) || state == 'COMPLETED')) return(.kinf_read_results(job))
  stop('job ', job$job_id, ' ended as ', state, ' without results:\n', .kinf_logtail(job))
}

#' @export
print.kinf_job = function(x, ...){
  cat('kinf job', x$job_id, '(', x$backend, ') in', x$workdir, '\n')
  invisible(x)
}

#' Read a gwalks.json written by infer.py (signed snode ids + circular flags) as a gWalk on graph
#' (the gGraph that was exported), or as a plain list when graph is NULL.
read_kinf_gwalks = function(path, graph = NULL){
  x = jsonlite::read_json(path, simplifyVector = FALSE)  # keep length-1 walks as lists
  snode.id = lapply(x$snode.ids, function(w) as.integer(unlist(w)))
  circular = as.logical(unlist(x$circular))
  if (is.null(graph)) return(list(snode.id = snode.id, circular = circular))
  gW(graph = graph, snode.id = snode.id, circular = circular)
}

# The model check values stored in inst/python/kinf/data/fastkar_model.json, recomputed from the current
# package (used by scripts/export_kinf_model.R and the sync test).
.kinf_model_check = function(){
  fg = splineobj$diag_spline
  dist = c(1,10,100,1e3,5e3,1e4,3e4,1e5,1e6,1e7,1e8)
  wp = c(1e2,1e4,1e6,1e8,2.5e9,1e10,1e11)
  list(dist = dist, dens = distance_density(dist), wp = wp, diag = exp(unname(fg(log(wp)))))
}
