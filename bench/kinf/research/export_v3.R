# export_v3.R
# Harder benchmark (v3): CCLE complex-rearrangement events with a LARGE karyotype space K(G).
# Events/footprints/graphs follow the fastKar Fig. 4 local identifiability run
# (~/projects/karyotyping/data/figure4_local_identifiability_flow/ccle_v3_local_swap_event_assets.rds);
# per-event graph prep = identifiability_context(): buffered footprint (0.2, min 1e4) -> target_region,
# graph_cut = loosefix(graph$copy$disjoin(ft) %&% ft) -> the graph (all truths are sample.gwalks(graph_cut)).
# pix.size: buffered footprint tiled into ~150 tiles (multiple of 1 kb, >= 5 kb).
# Stages (Rscript bench/kinf/research/export_v3.R <stage>, run from the fastKar repo root; outputs under <root>/temp_data/v3/):
#   cand     candidate table (sizes of cut graphs; mcool lines only; types chromothripsis/tyfonas/chromoplexy/bfb/
#            pyrgo/rigma/tic/cpxdm/dm)  -> candidates.csv
#   probe    probe folder per kept candidate (one truth, depth 10) -> probe/<event_id>/ ; then run
#            python bench/kinf/research/kspace.py -n 20000 on them -> probe_kspace.tsv
#   export   (after choosing selected_events.txt) for each selected event: sim_t<i>_d<j>/ (5 truths x depth
#            c(1,10); noisy = make_noisymap(forward_simulate(...),2)[[1]]), realdata/, manifest.csv, events.csv
suppressMessages(devtools::load_all('~/git/fastKar'))
root = '/gpfs/commons/groups/imielinski_lab/projects/karyotyping/karyotype-inference/'
outdir = paste0(root,'temp_data/v3/')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
fdir = '~/projects/karyotyping/data/figure4_local_identifiability_flow/'
mk = function(f){dir.create(f,recursive=TRUE,showWarnings=FALSE);paste0(f,'/')}
stage = commandArgs(trailingOnly=TRUE)[1]
types = c('chromothripsis','tyfonas','chromoplexy','bfb','pyrgo','rigma','tic','cpxdm','dm')

load_event_graph = function(gg_path, footprint_str){
  gg = readRDS(normalizePath(path.expand(gg_path),mustWork=TRUE))
  if (!is.null(gg$copy)) gg = gg$copy
  fp = streduce(parse.gr(strsplit(footprint_str,';',fixed=TRUE)[[1]]))
  ft0 = streduce(fp)
  buffer = max(sum(width(ft0))*0.2,1e4)
  ft = streduce(ft0+buffer)
  gcut = loosefix(gg$copy$disjoin(ft) %&% ft)
  list(graph=gcut,ft=ft,ft0=ft0,target_width=sum(width(ft)))
}
pick_pix = function(w){ max(5e3, round(w/150/1e3)*1e3) }

if (stage=='cand'){
  ev = as.data.table(readRDS(paste0(fdir,'ccle_v3_local_swap_event_assets.rds')))
  pr = as.data.table(readRDS('~/projects/karyotyping/db/ccle.pairs.rds'))
  okp = pr[!is.na(mcool) & mcool!='' & file.exists(mcool)]$pair
  ev = ev[type %in% types & pair %in% okp]
  cat(nrow(ev),'events on',length(unique(ev$pair)),'mcool lines\n')
  res = mclapply(split(ev,by='pair'),function(d){
    rbindlist(lapply(1:nrow(d),function(r){
      e = d[r]
      tryCatch({
        x = load_event_graph(e$gg,e$footprint)
        n = x$graph$nodes$dt
        data.table(event_id=e$event_id,pair=e$pair,type=e$type,n_nodes=nrow(n),n_cn_gt1=sum(n$cn>1),max_cn=max(n$cn),sum_cn=sum(n$cn),
                   target_width=x$target_width,pix_size=pick_pix(x$target_width),err=NA_character_)
      },error=function(err) data.table(event_id=e$event_id,pair=e$pair,type=e$type,err=conditionMessage(err)),fill=TRUE)
    }),fill=TRUE)
  },mc.cores=12)
  cand = rbindlist(res,fill=TRUE)
  fwrite(cand,paste0(outdir,'candidates.csv'))
  print(table(cand$type,!is.na(cand$err)))
  print(cand[is.na(err)][order(type,sum_cn)][,.(event_id,type,n_nodes,n_cn_gt1,max_cn,sum_cn,target_width,pix_size)],nrows=500)
}

ev_all = function(){
  ev = as.data.table(readRDS(paste0(fdir,'ccle_v3_local_swap_event_assets.rds')))
  setkey(ev,event_id); ev
}
seed_of = function(...) strtoi(substr(digest::digest(paste(...),algo='xxhash32',serialize=FALSE),1,7),16L)

# one event: write a folder for a given truth seed / depth (noisy map, truth=)
export_one = function(x, truth, depth, seed, pix, folder){
  set.seed(seed)
  gm = make_noisymap(forward_simulate(truth,target_region=x$ft,depth=depth,pix.size=pix),2)[[1]]
  prep_data_for_py(gg=x$graph,gm=gm,depth=depth,pix.size=pix,target_region=x$ft,write.to.folder=mk(folder),truth=truth)
  invisible(NULL)
}

if (stage=='probe'){
  cand = fread(paste0(outdir,'candidates.csv'))[sum_cn<=150 & max_cn<=8]
  # spread over types: all chromothripsis/dm, then the largest sum_cn per type
  quota = c(chromothripsis=3,chromoplexy=8,pyrgo=9,rigma=8,tic=12,bfb=2,dm=2,tyfonas=2,cpxdm=2)
  setorder(cand,type,-sum_cn)
  if (!identical(commandArgs(trailingOnly=TRUE)[2],'all')) cand = cand[,head(.SD,quota[type[1]]),by=type]  # 'all': probe every kept candidate
  cand = cand[!file.exists(paste0(outdir,'probe/',event_id,'/truth_expected.csv'))]
  cat(nrow(cand),'probes\n'); print(table(cand$type))
  ev = ev_all()
  res = mclapply(cand$event_id,function(id){
    tryCatch({
      e = ev[id]; x = load_event_graph(e$gg,e$footprint)
      pix = pick_pix(x$target_width)
      sd = seed_of('probe',id); set.seed(sd)
      truth = sample.gwalks(x$graph,1,return.gw=TRUE,verbose=FALSE)[[1]]
      export_one(x,truth,10,sd,pix,paste0(outdir,'probe/',id))
      'ok'
    },error=function(err) paste('ERR',conditionMessage(err)))
  },mc.cores=12)
  print(table(unlist(res))); print(unlist(res)[unlist(res)!='ok'])
}

if (stage=='export'){
  sel = readLines(paste0(outdir,'selected_events.txt'))
  ev = ev_all()
  pr = as.data.table(readRDS('~/projects/karyotyping/db/ccle.pairs.rds')); setkey(pr,pair)
  depths = c(1,10); ntruth = 5L
  res = mclapply(sel,function(id){
    e = ev[id]; x = load_event_graph(e$gg,e$footprint); pix = pick_pix(x$target_width)
    man = list()
    for (i in 1:ntruth){
      tseed = seed_of('truth',id,i); set.seed(tseed)
      truth = sample.gwalks(x$graph,1,return.gw=TRUE,verbose=FALSE)[[1]]
      for (j in seq_along(depths)){
        f = paste0('sim_t',i,'_d',j); sd = seed_of('noise',id,i,j)
        export_one(x,truth,depths[j],sd,pix,paste0(outdir,id,'/',f))
        man[[length(man)+1]] = data.table(event_id=id,pair=e$pair,type=e$type,folder=paste0(id,'/',f),truth=i,truth_seed=tseed,depth=depths[j],seed=sd,pix_size=pix)
      }
    }
    # observed Hi-C: finest available resolution (1 kb) <= pix.size, rebinned to the export tiles by prep_data_for_py
    mc = pr[e$pair]$mcool; ploidy = pr[e$pair]$ploidy; if (is.na(ploidy)) ploidy = 2
    avail = as.integer(names(rhdf5::h5read(mc,'/resolutions'))); res_use = min(avail[avail<=pix])
    gm = cooler(mc,res=res_use,gr=x$ft)
    dep = estimate.depthratio(mc,ploidy=ploidy)
    prep_data_for_py(gg=x$graph,gm=gm,depth=dep,pix.size=pix,target_region=x$ft,write.to.folder=mk(paste0(outdir,id,'/realdata')))
    man[[length(man)+1]] = data.table(event_id=id,pair=e$pair,type=e$type,folder=paste0(id,'/realdata'),truth=NA_integer_,truth_seed=NA_integer_,depth=dep,seed=NA_integer_,pix_size=pix,mcool_res=res_use,ploidy=ploidy)
    rbindlist(man,fill=TRUE)
  },mc.cores=5)
  bad = sapply(res,function(r) inherits(r,'try-error') || !is.data.table(r)); print(bad); if (any(bad)) print(res[bad])
  man = rbindlist(res[!bad],fill=TRUE)
  fwrite(man,paste0(outdir,'manifest.csv'))
  # events.csv: sizes + kspace probe stats + Fig. 4 P(MLE = truth) per assay
  ps = fread(paste0(outdir,'probe_summary.csv'))[event_id %in% sel]
  s = as.data.table(readRDS(paste0(fdir,'ccle_v3_local_swap_event_summary.rds')))[event_id %in% sel]
  pw = dcast(s,event_id~assay,value.var='p_pick'); setnames(pw,names(pw)[-1],paste0('fig4_p_',gsub('[^A-Za-z0-9]+','_',names(pw)[-1])))
  nc = unique(s[,.(event_id,fig4_n_cand=n_cand,fig4_chance=chance)])
  fwrite(merge(merge(ps,pw,by='event_id'),nc,by='event_id'),paste0(outdir,'events.csv'))
  print(man)
}
