#' Export a gGraph + Hi-C map as a problem folder for the Python karyotype search (kinf)
#'
#' Splits every node of gg into a left and a right half, gives each copy of each half-node its
#' own id, fixes the junction copies to half-node copies, tiles the data bins at pix.size and at node boundaries
#' (the tiles forward_simulate() uses; each tile copy sits on the half-node copy containing its centre), and
#' rebins the observed map onto the data bins with rebin_matrix() (count-conserving). Nodes outside
#' the target are wired like the others but get no tiles: they enter the fit only through the walk distances
#' between target tiles. The only free
#' choice left is, at each node with CN > 1, which left-half copy continues into which right-half copy.
#' Used by infer_karyotype(); the folder is read by inst/python/kinf/problem.py (Graph.load, HiC.load).
#'
#' @param gg gGraph with node and junction copy numbers (cn).
#' @param gm gMatrix of observed counts.
#' @param depth Hi-C depth (fastKar convention).
#' @param purity only 1 is supported by the Python model.
#' @param pix.size simulation tile size (bp); tiles are cut inside the data bins, which are what is scored.
#' @param target_region GRanges of the Hi-C fit region; the data bins of gm overlapping it are scored (bins.csv). Default
#'   gg$footprint. Nodes outside it are still wired.
#' @param write.to.folder folder path ending in '/'; NULL: return the tables only.
#' @param truth optional gWalk on gg (simulations): also writes truth_gwalks.json and truth_expected.csv (the
#'   noiseless forward_simulate() map on the export's tiles, summed into the data bins; i, j index rows of bins.csv).
#'   gm_dat.csv holds the observed counts on the data bins (i, j index bins.csv rows); tiled_target.csv has each
#'   tile's bin.id.
#' @return named list of the exported tables.
prep_data_for_py <- function(gg,gm,depth,purity=1,pix.size=1e5,target_region=NULL,write.to.folder=NULL,truth=NULL){ 
  if (is.null(target_region)){
    target_region = gg$footprint
  }
  #pull data from graph
  nodesdt = gg$nodes$dt[,.(start,end,seqnames,snode.id,cn)]
  nodesgr = gg$nodes$gr
  edgesdt = gg$edges$dt[,.(cn,n1,n2,n1.side,n2.side,type)]

  #every node of the graph is wired (walks through nodes outside the target set the distances between target tiles), but only the target is tiled and scored
  targetnodes = gr2dt(nodesgr)

  #split nodes in half: our model of walks on the genome graph requires that we optimize over all possible wirings of left-junctions to right-junctions at a given node. To this end we consider the left- and right- halves of a node as distinct.
  nodesdt.split = targetnodes[rep(1:.N,each=2)][,side:=ifelse(mod(.I,2)==1,'left','right')][,.(start,end,seqnames,node.id,cn,side)][,width:=ifelse(side=='left',floor((end-start+1)/2),ceil((end-start+1)/2))]
  nodesdt.split[side=='left',end:=start+width-1]
  nodesdt.split[side=='right',start:=end-width+1]

  #now we consider each copy of each half-node as unique, and give it a unique ID. This is so we can keep track of junction mappings
  dedup.nodes = nodesdt.split[rep(1:.N,cn)][,.(start,end,seqnames,node_id = node.id,side)][,unique.id:=.I]

  #Now we have to tile the target region which we are going to simulate and fit to data, making sure our tiles match up with edges of nodes

  # the fit is scored on the data's own bins (2026-11-11, user: no assumption about how counts are spread inside a bin;
  # reload the data at a finer resolution if needed): bins = the observed bins overlapping the target. The model is
  # simulated on tiles = bins cut at pix.size and at node boundaries (exactly the tiles forward_simulate() uses; not at
  # half-node boundaries, 2026-11-12), so every tile lies inside one bin, and the expected map is summed into bins before
  # scoring (as bestfit_search() simulates on the observed bins and compmaps() rebins onto them).
  bins = gm$gr[gm$gr %^% target_region]
  bins = GRanges(as.character(seqnames(bins)),ranges(bins))   # no seqinfo: a map's seqlengths can differ from the graph's
  tiled.target = gr.merge(gr.tile(bins,pix.size),dt2gr(targetnodes[cn>0,.(seqnames,start,end)]))[,c()]
  tb = findOverlaps(GRanges(as.character(seqnames(tiled.target)),ranges(tiled.target)),bins,type='within')
  if (length(tb) != length(tiled.target) || any(duplicated(queryHits(tb)))){stop('prep_data_for_py: a tile is not inside exactly one data bin')}
  tiled.target$bin.id = subjectHits(tb)[order(queryHits(tb))]

  # every copy of every tile: the forward model treats a tile as a point at its centre (forward_kernel.cpp: mid), so a
  # tile copy belongs to the copy of the half-node that contains its centre (one copy per node copy), even when the tile
  # spans the node's midpoint; leftdist / rightdist = distance from the centre to that half's left / right end (+ 1/2),
  # through which kinf measures distances along the walk.
  tc = gr2dt(tiled.target)[,.(seqnames,start,end,width)][,tile.id:=.I]
  ctr = tc[,.(seqnames,start=floor((start+end)/2),end=floor((start+end)/2),tile.id)]
  ov = findOverlaps(GRanges(as.character(ctr$seqnames),IRanges(ctr$start,ctr$end)),
                    GRanges(as.character(dedup.nodes$seqnames),IRanges(dedup.nodes$start,dedup.nodes$end)))
  dedup.tiles = data.table(tile.id=ctr$tile.id[queryHits(ov)],unique.id=dedup.nodes$unique.id[subjectHits(ov)])
  dedup.tiles[,`:=`(width=tc$width[tile.id],start=tc$start[tile.id],end=tc$end[tile.id])]
  dedup.tiles[,leftdist:=start-dedup.nodes[unique.id]$start + width/2]
  dedup.tiles[,rightdist:=dedup.nodes[unique.id]$end - end + width/2]
  dedup.tiles = dedup.tiles[order(tile.id,unique.id),.(width,unique.id,tile.id,leftdist,rightdist)]

  # now we actually set up the machinery for permutations at a given node. keep track of every unique junction
  dedup.edges = edgesdt[rep(1:.N,cn)][,cn:=NULL][,edge_id:=.I]
  edges_long = data.table::melt(dedup.edges,
                     id.vars = c("edge_id", "type"),
                     measure.vars = list(c("n1", "n2"), c("n1.side", "n2.side")),
                     variable.name = "role", 
                     value.name = c("node_id", "side"))
  edges_long[, occ := seq_len(.N), by = .(node_id, side)]
  dedup.nodes[, occ := seq_len(.N), by = .(node_id, side)]
  edges_long = merge.data.table(edges_long,
                      dedup.nodes[, .(node_id, side, occ, unique.id)],
                      by = c("node_id", "side", "occ"),
                      all.x = TRUE)
  edges_unique = dcast(edges_long, edge_id + type ~ role, value.var = "unique.id")
  setnames(edges_unique, c("1", "2"), c("unique.id.n1", "unique.id.n2"))

  # two data tables. External edges are fixed junctions connecting two nodes, while internal edges are the edges connecting the two halves of a given node, which we can optimize over by permuting them appropriately
  external_edges = merge.data.table(dedup.edges, edges_unique, by = c("edge_id", "type"))[,.(edge_id,n1 = unique.id.n1, n2 = unique.id.n2)]
  internal_edges = cbind(dedup.nodes[side=='left'][,.(node_id,left=unique.id)],dedup.nodes[side=='right'][,.(right=unique.id)])[,copies:=.N,by=node_id]

  # We do a final merge of the nodes
  dedup.nodes = merge.data.table(dedup.nodes,nodesdt.split[,.(node_id=node.id,side,width)],by=c('node_id','side'),all.x=TRUE)[,.(node_id,side,unique.id,width)]

  # observed counts on the data bins, as they are (no split, no rounding); i, j index bins.csv rows
  gm_rebin = rebin_matrix(gm,bins)
  
  # put in a named list and write to files for Python to read
  data.out = list(tiled_target = gr2dt(tiled.target),dedup_tiles = dedup.tiles, nodes = dedup.nodes, internal_edges = internal_edges, external_edges = external_edges,bins = gr2dt(bins),gm_dat = gm_rebin$dat,depth=depth,purity=purity)
  if (!is.null(write.to.folder)){
    lapply(1:length(data.out),function(i){
      thisname = names(data.out)[i]
      write.csv(data.out[[i]],paste0(write.to.folder,thisname,'.csv'))
    })
  }
  # optional ground truth (a gWalk): noiseless R map at the same depth, simulated on the EXPORT's own tiling
  # (rows of tiled_target.csv, tile.id = row index), then summed into the data bins: truth_expected.csv i,j index bins.csv rows
  if (!is.null(truth) && !is.null(write.to.folder)){
    nodesgr.t = truth$graph$gr[,c('node.id','cn')]
    nodesgr.t = nodesgr.t[strand(nodesgr.t)=='+']
    tt.gr = tiled.target[,c()]
    ov = gr.findoverlaps(tt.gr,nodesgr.t)
    ov = ov[!duplicated(ov$query.id)]
    tt.dt = gr2dt(tt.gr)[,.(start,end,seqnames,width)][,tile.id:=.I]
    tt.dt[ov$query.id,`:=`(node.id=nodesgr.t$node.id[ov$subject.id],cn=nodesgr.t$cn[ov$subject.id])]
    tt.dt[is.na(node.id),cn:=0]
    tt.dt[,compartment:='A']
    prepped.t = list(tiled.target=tt.dt[,.(start,end,seqnames,tile.id,width,node.id,cn,compartment)],widthdt=gr2dt(nodesgr.t)[,.(nodewidth=width,node.id,start,end)])
    truth_sim = forward_simulate(truth,pix.size=pix.size,depth=depth,purity=purity,prepped=prepped.t)
    # summed into the data bins (tiles nest in bins, so the rebin is an exact sum): i, j index bins.csv rows
    truth_sim = rebin_matrix(gM(gr=tiled.target[,c()],dat=truth_sim$dat),bins)
    write.csv(truth_sim$dat[,.(i,j,value)],paste0(write.to.folder,'truth_expected.csv'),row.names=FALSE)
    writeLines(jsonlite::toJSON(list(snode.ids = lapply(truth$snode.id,function(x){as.integer(x)}), circular = as.logical(truth$circular)),auto_unbox=FALSE),paste0(write.to.folder,'truth_gwalks.json'))
  }
  return(data.out)
}

bestfit_search = function(graph,hic,nsample,depth,pix.size=NULL,target_region=NULL,topk=1,mc.cores=1,keep.circular=T,purity=1,verbose=F,rescale=F){
	if (is.null(target_region)){
		target_region = graph$footprint
	}
	if (is.null(pix.size)){
		pix.size = median(width(hic$gr))
	}
	gw = sample.gwalks(graph,nsample,verbose=F,mc.cores=mc.cores,keep.circular=keep.circular)
	# simulate on the observed bins that overlap target_region: every simulated tile then lies inside one
	# observed bin, so compmaps() sums tiles exactly instead of splitting them across bins
	sim_region = target_region
	if (!is.null(target_region)){
		obs_bins = hic$gr[,c()]
		obs_bins = obs_bins[obs_bins %^% target_region]
		if (length(obs_bins)){sim_region = obs_bins}
	}
	# the tiling depends only on the graph, target and resolution: prepare it once for all candidates
	prepped = NULL
	if (length(gw) && !is.null(sim_region)){
		g1 = gw[[1]]$graph
		if (all(vapply(gw, function(w){identical(w$graph, g1)}, logical(1)))){
			prepped = prep_for_sim(gw[[1]],sim_region,pix.size)
		}
	}
	if (verbose==F){
		sims = mclapply(gw,function(w){forward_simulate(w,target_region=sim_region,pix.size=pix.size,depth=depth,purity=purity,prepped=prepped)},mc.cores=mc.cores)
	}else{
		sims = pbmclapply(gw,function(w){forward_simulate(w,target_region=sim_region,pix.size=pix.size,depth=depth,purity=purity,prepped=prepped)},mc.cores=mc.cores)
	}
	# The genome-wide depth estimate can be off locally, so rescale all candidates by one shared factor
	# to the observed total. Use the median candidate total so that no single sample sets it.
	#ratio = sum(hic$value)/median(unlist(lapply(sims,function(s){sum(s$value)})))
	#neglogliks = unlist(mclapply(sims,function(s){compmaps(hic,s*ratio,theta=2,ifsum=T)},mc.cores=mc.cores))
	# rescale = TRUE: one shared ratio = observed total / median candidate total, both on the observed bins
	ratio = 1
	if (rescale){
		ratio = sim_depth_ratio(hic,sims,mc.cores=mc.cores)
		sims = lapply(sims,function(s){s*ratio})
	}
	neglogliks = unlist(mclapply(sims,function(s){compmaps(hic,s,theta=2,ifsum=T)},mc.cores=mc.cores))
	if (topk==-1){
		best = order(neglogliks)
	}else{
		best = order(neglogliks)[1:topk]
	}
	return(list(true_hic = hic,kars=gw[best],nlls=neglogliks[best],sims = sims[best],purity=purity,ratio=ratio))
}

#' Shared depth correction for a set of simulated maps: observed total / median simulated total,
#' with every map rebinned to the observed bins first (as compmaps() does).
sim_depth_ratio = function(hic,sims,mc.cores=1){
	obs_total = sum(rebin_matrix(hic,hic$gr)$dat$value,na.rm=T)
	sim_totals = unlist(mclapply(sims,function(s){sum(rebin_matrix(s,hic$gr)$dat$value,na.rm=T)},mc.cores=mc.cores))
	obs_total/median(sim_totals)
}
