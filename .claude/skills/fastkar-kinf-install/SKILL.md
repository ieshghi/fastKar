---
name: fastkar-kinf-install
description: Set up fastKar's karyotype inference (kinf; Python/JAX search behind infer_karyotype()) on a new machine or HPC cluster. Use when the user wants to install or configure kinf, install_kinf(), GPU/CUDA or Slurm settings for infer_karyotype(), or fastKar's R dependencies.
---

# Install fastKar karyotype inference (kinf)

The installer logic lives in R: `install_kinf()` builds a dedicated Python environment (standalone uv Python +
pinned JAX/numpy/pandas from `inst/python/kinf/requirements-*.txt`) and writes the config
`tools::R_user_dir('fastKar', 'config')/kinf.json`. Your job is to survey the machine, choose the arguments,
run it, and report. Keep it thin; don't reimplement the installer.

## 1. Survey (read-only; report what you find)
- OS / arch: `uname -sm` (cuda12 build: Linux x86_64 only).
- GPU: `nvidia-smi` (driver version; the CUDA 12 wheels need a recent NVIDIA driver, no system CUDA).
  Login nodes often have no GPU: check the compute nodes via the scheduler instead.
- Slurm: `sinfo -o '%P %G %l %a'` (GPU partitions, gres names such as `gpu:l40s:1`, time limits);
  `sacctmgr show assoc user=$USER format=account,partition` (is `--account` needed?).
- Filesystems: the env and run folders must be visible from compute nodes. `stat -f -c %T <dir>` (nfs/gpfs/lustre =
  shared; tmpfs/xfs under /tmp = node-local); check quota (`quota -s`, `mmlsquota`, `lfs quota`). cuda12 needs
  ~4-6 GB, cpu ~1 GB. Prefer a project or scratch space over a small home quota.
- R: `R --version` (fastKar needs R >= 4.3), `Rscript -e '.libPaths()'`.
- Environment pitfalls: a set `PYTHONPATH` or module-system Python libraries in `LD_LIBRARY_PATH` (kinf unsets
  them for its own processes; mention them anyway), conda/venv activation, `~/.Rprofile` side effects.
- Internet: `curl -sI https://pypi.org/simple/jax/` and `https://github.com` (uv and wheels are downloaded).
  Without internet: `uv = FALSE` with an existing Python >= 3.11 and a local wheel mirror.

## 2. Decide (ask the user when unclear)
- `build`: `'cuda12'` if the search will run on NVIDIA GPU nodes (even if installing from a GPU-less login
  node), else `'cpu'` (~20x slower; fine for small graphs).
- `envdir`: a shared directory with enough quota, e.g. `/scratch/<lab>/envs/fastkar_kinf`.
- Slurm defaults: `slurm = list(partition =, gres =, time =, mem = '32G', cpus = 8, account =)`.

## 3. fastKar's R dependencies (if fastKar isn't installed yet)
Remotes in DESCRIPTION: `mskilab-org/gGnome`, `gUtils`, `GxG`, `gTrack`, `skitools`. Ask which library to install
into (e.g. a project library via `.libPaths()`), then
`remotes::install_github(c('mskilab-org/gUtils', 'mskilab-org/gGnome', 'mskilab-org/GxG', 'mskilab-org/gTrack',
'mskilab-org/skitools'), lib = LIB)` and `remotes::install_github('ieshghi/fastKar', lib = LIB)` (or
`R CMD INSTALL -l LIB .` from a clone). If the user's startup files load development packages, run installs with
`R_PROFILE_USER=/dev/null`.

## 4. Install and check
```r
library(fastKar)
install_kinf(envdir = '<shared dir>', build = 'cuda12',
             slurm = list(partition = 'gpu', gres = 'gpu:1', account = NULL))
kinf_check()                      # jax devices + self test: must print "OK"
kinf_check(backend = 'slurm')     # optional: same test as a short GPU job (ask first: uses their allocation)
```
`install_kinf(test_gpu_job = TRUE)` runs the Slurm check right after installing. Rerunning `install_kinf()` is
safe; `force = TRUE` rebuilds the env.

## 5. Report
Tell the user: build, envdir (size), the Python path and config file written, the device seen by
`kinf_check()` (local and, if run, on the GPU node), and anything unresolved (quota, account, no GPU visible).

## Guardrails
- Use `install_kinf()`. Never pip-install into a system/base/conda-base Python or into the user's R libraries.
- Ask before: submitting cluster jobs, creating multi-GB environments, editing `~/.bashrc`, `~/.Rprofile` or
  other startup files, or choosing where to install R packages.
- Keep the pinned versions in `requirements-*.txt`; don't upgrade JAX ad hoc.
- If `kinf_check()` reports `device cpu` on a GPU machine: check `nvidia-smi`, `CUDA_VISIBLE_DEVICES`, and that
  the env was built with `build = 'cuda12'`.

This skill is picked up automatically when Claude Code runs inside a clone of fastKar. If fastKar was installed
with `install_github`, point your agent to
https://github.com/ieshghi/fastKar/blob/master/.claude/skills/fastkar-kinf-install/SKILL.md.
