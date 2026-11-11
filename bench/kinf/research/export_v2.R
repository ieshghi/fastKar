# export_v2.R
# Re-exports karyotype-inference benchmark problems with the fixed (count-conserving) rebin of the
# observed map (fastKar prep_data_for_py, uses fastKar::rebin_matrix instead of gm$disjoin()$agg()).
# Graph = wsamp[[1]]$graph from sampled_walks_chr3.rds (all 100 truths share it; checked).
# Outputs under <project>/temp_data/v2/:
#   realdata/            SK-OV-3 observed Hi-C (res 1e4 cooler), region, pix.size 5e4
#   sim_t<i>_d<j>/       truth i=1..33 x depth fraction j=1..4 of c(0.03,0.1,0.3,1)*SK-OV-3 depth;
#                        data = make_noisymap(forward_simulate(truth,...),2)[[1]]; plus truth_expected.csv
#                        (on the export tiling, i,j = rows of tiled_target.csv), truth_gwalks.json
#   manifest.csv         folder, truth index, depth fraction, depth, seed
# Usage: Rscript bench/kinf/research/export_v2.R   (run from the fastKar repo root)
suppressMessages(devtools::load_all('~/git/fastKar'))
root = '/gpfs/commons/groups/imielinski_lab/projects/karyotyping/karyotype-inference/'
outdir = paste0(root,'temp_data/v2/')
pairs = readRDS('/gpfs/commons/home/sclarke/Projects/jabba_hic/db/pairs.rds')
region = parse.gr("1:34603609-37500000,2:121488406-123511096,5:138827707-140552898,20:41565779-43466165")
wsamp = readRDS(paste0(root,'sampled_walks_chr3.rds'))
skov.bench = wsamp[[1]]$graph
skov.hic = cooler(pairs[pair=='SK-OV-3']$mcool, res=as.integer(1e4), gr=region)
skov.depth = estimate.depthratio(pairs[pair=='SK-OV-3']$mcool)
cat('depth',skov.depth,'\n')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
mk = function(f){dir.create(f,showWarnings=FALSE);paste0(f,'/')}
prep_data_for_py(gg=skov.bench,gm=skov.hic,depth=skov.depth,pix.size=5e4,target_region=region,write.to.folder=mk(paste0(outdir,'realdata')))
fracs = c(0.03,0.1,0.3,1)
man = CJ(truth=1:33,dfrac_idx=1:4)
man[,`:=`(folder=paste0('sim_t',truth,'_d',dfrac_idx),depth_fraction=fracs[dfrac_idx],depth=skov.depth*fracs[dfrac_idx],seed=1000L*truth+dfrac_idx)]
res = mclapply(1:nrow(man),function(r){
  m = man[r]
  set.seed(m$seed)
  truth = wsamp[[m$truth]]
  gm = make_noisymap(forward_simulate(truth,target_region=region,depth=m$depth,pix.size=5e4),2)[[1]]
  prep_data_for_py(gg=skov.bench,gm=gm,depth=m$depth,pix.size=5e4,target_region=region,write.to.folder=mk(paste0(outdir,m$folder)),truth=truth)
  sum(gm$dat$value[gm$dat$i<=gm$dat$j],na.rm=TRUE)
},mc.cores=12)
print(table(sapply(res,class)))
man[,obs_total:=unlist(res)]
fwrite(man[,.(folder,truth,depth_fraction,depth,seed)],paste0(outdir,'manifest.csv'))
