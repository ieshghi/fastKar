# make_kinf_selftest.R: write the kinf self-test problem inst/python/kinf/selftest/ (used by kinf_check()):
# the toy_repdup test fixture (6 nodes, CN 2-3), its first karyotype as the truth, a noisy simulated map at
# depth 10 on 100 kb tiles, exported with prep_data_for_py(..., truth =).
# Usage (from the repo root): Rscript scripts/make_kinf_selftest.R
suppressMessages(devtools::load_all('.', quiet = TRUE))
data.table::setDTthreads(1)
fx = readRDS('tests/testthat/fixtures/toy_repdup.rds')
truth = gW(graph = fx$gg, snode.id = fx$kars[[1]]$snode.id, circular = fx$kars[[1]]$circular)
depth = 10; pix = 1e5
set.seed(20261027)
gm = make_noisymap(forward_simulate(truth, target_region = fx$ft, depth = depth, pix.size = pix), 2)[[1]]
out = 'inst/python/kinf/selftest/'
unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
prep_data_for_py(gg = fx$gg, gm = gm, depth = depth, pix.size = pix, target_region = fx$ft, write.to.folder = out,
                 truth = truth)
cat(list.files(out), '\n')
