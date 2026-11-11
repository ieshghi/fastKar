# export_wg.R
# Whole-genome SK-OV-3 problems: CCLE complex.rds graph (SK-OV-3 is not in nyu_pairs), whole genome as target,
# tiles of pix.size. Truths sampled from the whole graph with sample.gwalks, simulated at depth 1 and 10 as in v3.
# Outputs under <project>/temp_data/wg_skov3_<pix>/: sim_t<i>_d<j>/, realdata/, manifest.csv
# Usage: Rscript bench/kinf/research/export_wg.R [pix.size=1e6] [ntruth=3]   (run from the fastKar repo root)
suppressMessages(devtools::load_all('~/git/fastKar'))
args = commandArgs(trailingOnly=TRUE)
pix = if (length(args) >= 1) as.numeric(args[1]) else 1e6
ntruth = if (length(args) >= 2) as.integer(args[2]) else 3L
root = '/gpfs/commons/groups/imielinski_lab/projects/karyotyping/karyotype-inference/'
outdir = paste0(root, 'temp_data/wg_skov3_', format(pix, scientific=FALSE), '/')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
mk = function(f){dir.create(f, recursive=TRUE, showWarnings=FALSE); paste0(f, '/')}
seed_of = function(...) strtoi(substr(digest::digest(paste(...), algo='xxhash32', serialize=FALSE), 1, 7), 16L)

pr = as.data.table(readRDS('~/projects/karyotyping/db/ccle.pairs.rds')); setkey(pr, pair)
p = pr['G16474.SK-OV-3']
gg = readRDS(p$complex)
if (!is.null(gg$copy)) gg = gg$copy
gg = loosefix(gg)
ft = streduce(gg$nodes$gr)
cat('nodes', length(gg$nodes), 'sum cn', sum(gg$nodes$dt$cn), 'target', sum(width(ft))/1e6, 'Mb\n')

depths = c(1, 10)
man = list()
for (i in 1:ntruth){
  tseed = seed_of('wgtruth', i); set.seed(tseed)
  t0 = Sys.time()
  truth = sample.gwalks(gg, 1, return.gw=TRUE, verbose=FALSE)[[1]]
  cat('truth', i, 'walks', length(truth), 'sampled in', format(Sys.time() - t0), '\n')
  for (j in seq_along(depths)){
    f = paste0('sim_t', i, '_d', j); sd = seed_of('wgnoise', i, j); set.seed(sd)
    t0 = Sys.time()
    gm = make_noisymap(forward_simulate(truth, target_region=ft, depth=depths[j], pix.size=pix), 2)[[1]]
    prep_data_for_py(gg=gg, gm=gm, depth=depths[j], pix.size=pix, target_region=ft, write.to.folder=mk(paste0(outdir, f)), truth=truth)
    cat(f, 'done in', format(Sys.time() - t0), '\n')
    man[[length(man)+1]] = data.table(folder=f, truth=i, truth_seed=tseed, depth=depths[j], seed=sd, pix_size=pix)
  }
}
# observed Hi-C at the export resolution (finest available <= pix), rebinned by prep_data_for_py
mc = p$mcool; ploidy = p$ploidy; if (is.null(ploidy) || is.na(ploidy)) ploidy = 2
avail = as.integer(names(rhdf5::h5read(mc, '/resolutions'))); res_use = max(avail[avail <= pix])
gm = cooler(mc, res=res_use, gr=ft)
dep = estimate.depthratio(mc, ploidy=ploidy)
prep_data_for_py(gg=gg, gm=gm, depth=dep, pix.size=pix, target_region=ft, write.to.folder=mk(paste0(outdir, 'realdata')))
man[[length(man)+1]] = data.table(folder='realdata', depth=dep, pix_size=pix, mcool_res=res_use, ploidy=ploidy)
fwrite(rbindlist(man, fill=TRUE), paste0(outdir, 'manifest.csv'))
cat('done\n')
