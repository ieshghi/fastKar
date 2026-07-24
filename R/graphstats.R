#' @useDynLib fastKar, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#goes from a ggraph to a "wiring", which gives all the internal edges of the graph (going from left side of a node to right side)
#along with the loose node ids and a reference data table with the new node ids (all copies of each node are de-duplicated)
gg.to.wiring = function(gg){ 
  if(!(('loose.cn.left' %in% colnames(gg$nodes$dt)) & ('loose.cn.right' %in% colnames(gg$nodes$dt)))){
  	gg = loosefix(gg)
  }
  nodesdt = gg$nodes$dt[,.(start,end,seqnames,snode.id,cn,loose.cn.left,loose.cn.right)]
  nodesgr = gg$nodes$gr[,c('snode.id','cn')]
  edgesdt = gg$edges$dt[,.(cn,n1,n2,n1.side,n2.side,type)]
  #make loose edges from nodesdt, then combine them will all other edges, then split all edges as left- and right-edges coming off of nodes
  nodesdt[is.na(loose.cn.right),loose.cn.right:=cn]
  nodesdt[is.na(loose.cn.left),loose.cn.left:=cn] #adjusting CN at chromosome ends
  # separate left and right loose ends and label them appropriately as edges going to node "0", before adding them to the total set of external edges
  left.looseedges = data.table::copy(nodesdt)[loose.cn.left>0][,.(n2=snode.id,cn = loose.cn.left,n2.side='left',n1=0,n1.side='right')]
  right.looseedges = data.table::copy(nodesdt)[loose.cn.right>0][,.(n1=snode.id,cn = loose.cn.right,n2.side='left',n2=0,n1.side='right')]
  external.edges = rbind(edgesdt,rbind(left.looseedges,right.looseedges)[,type:='LOO'][,.(cn,n1,n2,n1.side,n2.side,type)])
  dedup.edges = external.edges[cn>0][rep(1:.N,cn)][,.(n1,n2,n1.side,n2.side,type)] #separate all copies of all edges
  split.edgetable = melt.data.table(dedup.edges[,.(n1=paste0(n1,substr(n1.side,1,1)),n2=paste0(n2,substr(n2.side,1,1)),subid=.I,type)],id=c('subid','type'))[,.(subid,n=value,type)] #separate all edges and label them with a unique ID "subid"
  left.split.edgetable = split.edgetable[grepl('l',n)][,n:=as.integer(substr(n,1,nchar(n)-1))][n>0] %>% setkeyv('n') #setkey important so we have all left and right edges ordered appropriately.
  right.split.edgetable = split.edgetable[grepl('r',n)][,n:=as.integer(substr(n,1,nchar(n)-1))][n>0] %>% setkeyv('n')
  # each copy of each node has two edges coming out of it, one on the left one on the right. Each of those edges has a unique edge id which we call "subid"
  internal.edges = cbind(left.split.edgetable[,.(n,left=subid)],right.split.edgetable$subid)[,.(n,left,right=V2)][,cn:=.N,by=n] #each unique mapping of right-edged to left-edges (reorderings of the third column) corresponds to a unique walk decomposition of the graph
  # where are the loose edges?
  loose.subids = split.edgetable[type=='LOO']$subid %>% unique
  #split the nodes in the graph to half-nodes. This needs to be done only one time. Then, aggregate Hi-C data to those half-nodes.
  internal.edges[,id:=.I]
  return(list(internal.edges = internal.edges,loose.ends = loose.subids,gg=gg))
}

#sample gWalks from graph gg, take N samples and return all unique permutations among them
#' @import pbmcapply
#' @import digest
sample.gwalks = function(gg,N=1,mc.cores=1,chunksize = 1e3,return.gw=T,remove.dups=T,verbose=T,onlyhash=F,keep.circular=T,frozen.nodes = NULL){
  wiring = gg.to.wiring(gg)
  internal.edges = wiring$internal.edges
  loose.ends = wiring$loose.ends

  hash_fn  =  hash_karyotype_cpp

  if(verbose){
  message('Sampling permutations')
  }
  n_groups   <- split(seq_len(nrow(internal.edges)), internal.edges$n)
  n_vals = as.integer(names(n_groups))
  right_vec0 <- internal.edges$right
  shuffle_edges <- function(edges) {
    out <- right_vec0
    for (i in 1:length(n_groups)){
            g = n_groups[[i]]
            if (length(g) > 1 & !(n_vals[i] %in% frozen.nodes)) out[g] <- sample(out[g])
    }
    out
  }
  shuffle_chunk <- function(K, edges) {
    replicate(K, shuffle_edges(edges), simplify = FALSE)
  }
  chunks <- N %/% chunksize
  if (chunks > 1){
	  if(verbose){
  	perms <- do.call('c',pbmclapply(seq_len(chunks), function(i) shuffle_chunk(chunksize, internal.edges),mc.cores = mc.cores))
  }else{
  	perms <- do.call('c',mclapply(seq_len(chunks), function(i) shuffle_chunk(chunksize, internal.edges),mc.cores = mc.cores))
	  }
  }else{
	if(verbose){
  		perms <- pbmclapply(seq_len(N), function(i) shuffle_edges(internal.edges),mc.cores = mc.cores)
	}else{
  		perms <- mclapply(seq_len(N), function(i) shuffle_edges(internal.edges),mc.cores = mc.cores)
	}
  }

  if(verbose & remove.dups){
  	message('Only keeping unique permutations')
  }
  hashes = unlist(mclapply(perms,digest,mc.cores=mc.cores))
  dt = data.table(hash=hashes,idx=1:N)[,id:=1:.N,by=hash]
  if(remove.dups){uniqueperms = perms[dt[id==1]$idx]}else{uniqueperms = perms}
  permchunks = split(uniqueperms, ceiling(seq_along(uniqueperms)/chunksize))

  if(verbose){
  	message('Generating walks')
  }
  makewalk_chunk <- function(permchunk) {
    ws = traverse_graph_v2_batch_cpp(internal.edges, permchunk, loose.ends)
    if (return.gw){
      ws = lapply(ws,function(w){gW(graph=gg,snode.id=w$snode.id,circular=w$circular)})
    }
    if (onlyhash){
	return(lapply(ws,function(w){hash_fn(w$snode.id,w$circular)}))
    } else{
    	return(ws)
    }
  }
  if(verbose){
  	walks.out <- do.call('c',pbmclapply(permchunks, makewalk_chunk,mc.cores = mc.cores))
  }else{
  	walks.out <- do.call('c',mclapply(permchunks, makewalk_chunk,mc.cores = mc.cores))
  }
  if (onlyhash){
  		hash.dt = data.table(hash=unlist(walks.out))[,idx:=.I][,id:=1:.N,by=hash]
  		if (remove.dups){return(unique(hash.dt$hash))}else{return(hash.dt)}
  }else{if (remove.dups){
  	if (return.gw){
  	        hashes = do.call('c',mclapply(walks.out,function(w){w$hash},mc.cores=mc.cores))}
  	else{
  	        hashes = do.call('c',mclapply(walks.out,function(w){hash_fn(w$snode.id,w$circular)},mc.cores=mc.cores))}
  	hash.dt = data.table(hash=hashes)[,idx:=.I][,id:=1:.N,by=hash]
  	walks.out = walks.out[hash.dt[id==1]$idx]
  }
  if (!keep.circular){
	  keep = unlist(lapply(walks.out,function(w){sum(w$circular)==0}))
	  walks.out = walks.out[keep]
  }
  return(walks.out)
  }
}

#samples karyotype space using markov chains always starting at the same walk, initialized randomly. 
local.sampling = function(gg,nsteps,nwalk,onlyhash=F,starter_edges=NULL,return.edges=F,return.gw=F){
  wiring = gg.to.wiring(gg)
  shuffle_edges = function(edges) {
        new_right = edges[, if (.N > 1) sample(right, .N) else right, by = n]$V1
  	return(edges[,right:=new_right])
  }
  if (is.null(starter_edges)){
  	internal.edges = shuffle_edges(wiring$internal.edges)
  } else {
	internal.edges = starter_edges
  }
  loose.ends = wiring$loose.ends
  hashhist = c()
  permute.node = function(edges) {
  	if (nrow(edges[cn>1])==0){return(edges)}
	pivot.node = sample(edges[cn>1]$n,1)
	new_edges = data.table::copy(edges)
	edges.to.permute = edges[n==pivot.node]$right
	inds = sample(seq_along(edges.to.permute),2)
	edges.to.permute[c(inds[1],inds[2])] <- edges.to.permute[c(inds[2],inds[1])]
  	new_edges[n==pivot.node,right:=edges.to.permute]
	return(new_edges)
  }
  gw0 = traverse_graph_cpp(internal.edges,loose.ends)
  gw0$hash = hash_snodelist(gw0$snode.id,gw0$circular)
  walkhist = lapply(1:nwalk,function(i){list(gw0)}) #initialize nwalk walkers at the same point
  hashhist = lapply(1:nwalk,function(i){gw0$hash}) #initialize hashes
  edges = lapply(1:nwalk,function(i){internal.edges}) #initialize edge table
  for (i in seq_len(nwalk)){
  for (j in seq_len(nsteps-1)){
	newedges = permute.node(edges[[i]])	
	newwalk = traverse_graph_cpp(newedges,loose.ends)
	newhash = hash_snodelist(newwalk$snode.id,newwalk$circular)
	edges[[i]] = newedges
	newwalk$hash = newhash
	walkhist[[i]][[j+1]] = newwalk 
	hashhist[[i]] = c(hashhist[[i]],newhash)	
  }}
  if (return.edges){
	hashhist = list(hashes=hashhist,starter_edges = internal.edges)
	walkhist = list(walks=walkhist,starter_edges = internal.edges)
  }
  if (onlyhash){
	  return(hashhist)
  } else{
	  return(walkhist)
  }
}

#starts a markov chain at a random location and samples for a given length
markov.gwalk = function(gg,len,self.avoid = F,attempts = 10,return.gw=F,seed=NULL){
  if (is.null(seed)){set.seed(sample(1:1e5,1))
  } else{set.seed(seed)}
  wiring = gg.to.wiring(gg)
  internal.edges = wiring$internal.edges
  loose.ends = wiring$loose.ends
  hashhist = c()
  walkhist = list()
  permute.node = function(edges) {
	pivot.node = sample(edges[cn>1]$n,1)
	new_edges = data.table::copy(edges)
	edges.to.permute = edges[n==pivot.node]$right
	inds = sample(seq_along(edges.to.permute),2)
	edges.to.permute[c(inds[1],inds[2])] <- edges.to.permute[c(inds[2],inds[1])]
  	new_edges[n==pivot.node,right:=edges.to.permute]
	return(new_edges)
  }
  shuffle_edges = function(edges) {
        new_right = edges[, if (.N > 1) sample(right, .N) else right, by = n]$V1
  	return(edges[,right:=new_right])
  }
  gw0 = traverse_graph_cpp(shuffle_edges(internal.edges),loose.ends)
  gw0$hash = hash_snodelist(gw0$snode.id,gw0$circular)
  walkhist[[1]] = gw0
  hashhist = c(hashhist,gw0$hash)
  for (i in seq_len(len-1)){
	valid = F
  	killwalk = T #gets switched to F if we find the next step
  	for (j in seq_len(attempts)){
		if (!valid){
			edge_try = permute.node(internal.edges)	
			walk_try = traverse_graph_cpp(edge_try,loose.ends)
			hash_try = hash_snodelist(walk_try$snode.id,walk_try$circular)
			if (self.avoid){
				valid = !(hash_try %in% hashhist)
			} else{valid=T}
		} else{
			internal.edges = edge_try
			walk_try$hash = hash_try
			walkhist[[i+1]] = walk_try
			hashhist = c(hashhist,hash_try)	
			killwalk = F
			break
		}
	} #give a search ending condition, in case we get stuck return a shorter walk
	if (killwalk){
		break
	}
  }
  if (return.gw){
	  gwhist = lapply(walkhist,function(w){gW(graph=gg,snode.id=w$snode.id,circular=w$circular)})
	  return(gwhist)
  }else{return(walkhist)}
}

#runs markov.gwalk in parallel starting at random locations. different from local.sampling because that function starts all markov chains at the same location!
multi.markov = function(gg,N,len,self.avoid=F,attempts=10,mc.cores=1,return.gw=F,seed=NULL){
	if (is.null(seed)){set.seed(sample(1:1e5,1))
	}else{set.seed(seed)}
	seeds = sample(1:1e5,N)
	mclapply(1:N,function(i){
		markov.gwalk(gg=gg,len=len,self.avoid=self.avoid,attempts=attempts,return.gw=return.gw,seed=seeds[i])  
  },mc.cores=mc.cores)
}

booth_rotate = function(x) { #an implementation of Booth's algorithm to disambiguate circular walk hashes. See https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
	booth = function(s) {
		n = length(s)
		if (n == 0) return(1L)
		s2 = c(s, s)
		i = 1L
		j = 2L
		k = 0L
		while (i <= n && j <= n && k < n) {
			a = s2[i + k]
			b = s2[j + k]
			if (a == b) {
				k = k + 1L
			} else if (a > b) {
				# rotation at i is worse than rotation at j -> skip i's prefix
				i = i + k + 1L
				if (i <= j) i = j + 1L
				k = 0L
			} else {
				# rotation at j is worse -> skip j's prefix
				j = j + k + 1L
				if (j <= i) j = i + 1L
				k = 0L
			}
		}
		pos = min(i, j)
		# ensure returned index is in 1..n (not > n)
		if (pos > n) pos = pos - n
		return(pos)
	}
	start = booth(x)
	n = length(x)
	if (start==1){
		return(x)
	}else{
		return(x[c(start:n, 1:(start-1))])
	}
}

hash_snodelist = function(snode.id,circular){
	sorted = sort_snodes(snode.id,circular)
	snode.id = sorted$nodelist
	circular = sorted$circ
	circ = ifelse(rep(circular,2),'C','L')
	nodepcomp = c(snode.id,lapply(snode.id,function(s){-rev(s)}))
	nodestring = lapply(1:length(nodepcomp),function(i){
		return(paste0(toString(nodepcomp[[i]]),circ[i]))
	})
	return(toString(sort(do.call('c',nodestring))))
}

sort_snodes = function(nodelist,circ=NULL) {
	# Canonicalize circular walks to their Booth lex-min rotation in place;
	# we still need the raw walk later to compute Booth(rc(W)) properly.
	if (sum(circ)>0){
		circ_walks = nodelist[circ]
		circ_walks_rot = lapply(circ_walks,function(w){
			return(booth_rotate(w))
	})
		nodelist[circ] = circ_walks_rot
	}
	# For each walk pick the lex-min over its rotation+RC orbit.
	# Linear:   compare W and rc(W).
	# Circular: compare Booth(W) and Booth(rc(W)). The earlier code compared
	#   Booth(W) to rc(Booth(W)), which is wrong because rc(Booth(W)) is not
	#   generally the Booth rotation of rc(W) -- they differ by a cyclic
	#   shift. That left RC-equivalent circular karyotypes presented in
	#   different rotations with different canonical forms.
	choose_compl = lapply(seq_along(nodelist), function(i) {
		x = nodelist[[i]]
		rc = -rev(x)
		if (!is.null(circ) && length(circ) >= i && isTRUE(circ[i])) {
			rc = booth_rotate(rc)
		}
		if (paste(x, collapse = ",") <= paste(rc, collapse = ",")) {
			x
		} else {
			rc
		}
	})
	ord <- order(sapply(choose_compl, paste, collapse = ","))
	sorted_nodes = choose_compl[ord]
	if (!is.null(circ)){
		sorted_circ = circ[ord]
		return(list(nodelist=sorted_nodes,circ=sorted_circ))
	}else{
		return(sorted_nodes)
	}
}

#should introduce a threshold width for removing del/dups
smoothdeldups = function(ggraph,res=NULL){
	#dup_junctions = ggraph$edgesdt[!is.na(ggraph$edgesdt$dup)][n1==n2]
	dup_junctions = ggraph$edgesdt[n1==n2 & n1.side!=n2.side]
	del_junctions = ggraph$edgesdt[!is.na(ggraph$edgesdt$del)][abs(n1-n2)==2]
	del_junctions[,this.n:=round((n1+n2)/2)]
	nodesgr = ggraph$nodes$gr
	cnvec = nodesgr$cn
	cnvec[dup_junctions$n1] = cnvec[dup_junctions$n1]-dup_junctions$cn
	cnvec[del_junctions$this.n] = cnvec[del_junctions$this.n]+del_junctions$cn
	nodesgr$cn = cnvec
	nodesgr$fix = T
	nagraph = gG(nodes = ggraph$nodes$gr,edges = ggraph$edgesdt)
	nagraph$nodes$mark(cn=NA)
	nagraph$edges$mark(cn=NA)
	if (nrow(nagraph$edges$dt)==0){
		return(NULL)
	}else{
		nodeldups = loosefix(balance(nagraph,marginal=nodesgr,verbose=F)[,cn>0][cn>0])
		nodeldups$simplify()
		return(loosefix(nodeldups))
	}
}

sample_and_collapse = function(ggraph,samplesize,min.wid=3e7,bands=NULL,mc.cores=1){
	message('Sampling walks')
	walkset = sample.gwalks(ggraph,samplesize,return.gw=F,mc.cores=mc.cores)
  	gr = ggraph$nodes$gr[,c('node.id')]
	message('Collapsing to bands')
	if (is.null(bands)){
		merged.dt = gr2dt(gr)[,.(node.id,band.id=node.id,width,seqnames,start,end)]
	} else if (inherits(bands,'GRanges')){
 		merged.dt = gr2dt(gr.merge(gr,bands)[,c('query.id','subject.id')])[,.(node.id=query.id,band.id=subject.id,width=width,seqnames,start,end)]
	} else if (inherits(bands,'character')){
		bands.td = gTrack::karyogram(file = bands)
		bands = bands.td@data
		bands.gr = gr.nochr(grl.unlist(do.call(`GRangesList`, bands)))
		bands.gr = bands.gr %Q% (seqnames %in% c(1:22,'X','Y'))
		bands.gr = dt2gr(gr2dt(bands.gr)[,start:=start+1],seqlengths=seqlengths(bands.gr))
 		merged.dt = gr2dt(gr.merge(gr,bands.gr)[,c('query.id','subject.id')])[,.(node.id=query.id,band.id=subject.id,width=width,seqnames,start,end)]
	} else {
		error('Bands provided must either be GRanges or a path to cytobands file')
	}
	collapsed.walks = collapse_gwalklist(walkset,merged.dt,min.wid,mc.cores)

	message('Hashing and counting samples')
	coll_hashes = sapply(collapsed.walks, `[[`, "hash")
	coll_hash.dt = data.table(collapsed.hash=unlist(coll_hashes))[,walkset.id:=.I][,collapsed.id:=as.integer(factor(coll_hashes))]
	coll_hash.dt[,count:=.N,by=collapsed.id]
	setkeyv(coll_hash.dt,'collapsed.id')
	coll_hash.dt[,instance:=1:.N,by=collapsed.id]
	coll_idx = coll_hash.dt[instance==1]$walkset.id
	unique.collapsed = coll_hash.dt[instance==1,.(walkset.id,collapsed.id)]
	collapsed.walks = collapsed.walks[coll_idx]

	return(list(graph=ggraph,walkset=walkset,collapsed.walks=collapsed.walks,bands=dt2gr(merged.dt),hash.dt=coll_hash.dt))
}

#' @import pbapply
collapse_gwalklist <- function(gwlist,merged.dt,min.wid,mc.cores=1){
  merged.dt[,row_id:=.I]
  collapsed.list = pbmclapply(gwlist, function(gw) {
      walknodes = gw$snode.id
      circular = gw$circular
      collapsed.walk = lapply(walknodes,function(x){
        walk = merge.data.table(data.table(step=1:length(x),snode.id=x,node.id = abs(x),strand=sign(x)),merged.dt[,.(node.id,row_id,width,seqnames)],by='node.id',allow.cartesian=T)[,.(step,strand,snode.id,width,row_id,seqnames)]
	walk[strand==-1,row_id:=-rev(row_id),by=step]
	setkeyv(walk,'step')
	walk[,breakpt:=((row_id-shift(row_id)!=1)|(seqnames!=shift(seqnames)))][,breakpt:=ifelse(is.na(breakpt),T,breakpt)]
	walk[,run_id:=cumsum(breakpt)]
	walk[,small:=sum(width)<min.wid,by=run_id]
	walk[,runofruns:=cumsum(abs(rev(circdiff(rev(small),-1))))] #lol
	walk[,keep:=(!small | sum(width)>min.wid),by=runofruns]
	walk_kept = walk[keep==T][,.(row_id,seqnames,width,ufo=small)]
	walk_kept[,breakpt:=((row_id-shift(row_id)!=1)|(seqnames!=shift(seqnames)))][,breakpt:=ifelse(is.na(breakpt),T,breakpt)]
	walk_kept[,run_id:=cumsum(breakpt)]
	walk_kept[,row_out:=ifelse(ufo,NA,row_id)]
	outdat = walk_kept$row_out
	outdat = outdat[c(TRUE, !(is.na(outdat[-1]) & is.na(outdat[-length(outdat)])))]
	return(outdat)
      })
      return(list(snode.id = collapsed.walk,circular = circular,hash=hash_snodelist(collapsed.walk,circular)))
      }, mc.cores=mc.cores)
  return(collapsed.list)
}

to_gwalk = function(walklist,gr,mc.cores=1){
	grl = mclapply(walklist$snode.id,function(nl){
		nl = nl[!is.na(nl)]
		this.gr = gr[abs(nl)]
		strand(this.gr) = ifelse(nl>0,'+','-')
		if(length(this.gr)>0){return(this.gr)}else{return(NULL)}
        },mc.cores=mc.cores)
	keep = do.call('c',lapply(grl,function(x){!is.null(x)}))
	grl = do.call('GRangesList',grl[keep])
	return(gW(grl=grl,circular=walklist$circular[keep])$disjoin())
}

reads_fromwalk = function(walk,readL){
	gr = walk$graph$gr[,c('node.id')]
	reads = do.call('rbind',lapply(1:length(walk$snode.id),function(i){
		snodes = walk$snode.id[[i]]
		widths = width(gr[abs(snodes)])
		grdt = data.table(snode.id=snodes)
		grdt[,end:=cumsum(widths)]
		grdt[,start:=end-widths+1]
		total_length = max(grdt$end)
		if (walk$circular[i]) {
		    effective_readL = min(readL, total_length)
		    domain_lo = 1
		    domain_hi = total_length
		    grdt[, wrap := 0]
		    add.gr = copy(grdt)[,`:=`(start = start + total_length,end   = end   + total_length,wrap  = 1)]
		    max_read_end = total_length + effective_readL - 1
		    add.gr = add.gr[start <= max_read_end]
		    add.gr[end > max_read_end, end := max_read_end]
		    grdt <- rbindlist(list(grdt, add.gr),use.names = TRUE)
		} else {
		    effective_readL = readL
		    domain_lo = min(grdt$start) - effective_readL + 1
		    domain_hi = max(grdt$end)
		    grdt[, wrap := 0]
		}
    		setorder(grdt, start, end)
    		grdt[, node_order := .I]
    		cuts = sort(unique(c(domain_lo,domain_hi + 1,grdt$start -effective_readL + 1, grdt$end + 1)))
    		cuts = cuts[cuts >= domain_lo & cuts <= domain_hi + 1]
    		# Each row represents an interval of starts with an identical overlap word.
    		bins = data.table(read_start = head(cuts, -1),next_start = tail(cuts, -1))
    		bins = bins[read_start < next_start]
    		bins[, `:=`(bin_id = .I,start_max = next_start - 1,multiplicity = next_start - read_start,read_end = read_start + effective_readL - 1)]
    		# Find all nodes overlapped by the representative read from each bin.
    		nodes = grdt[, .(start,end,snode.id,node_order)]
    		setkey(nodes, start, end)
    		hits = foverlaps(x = bins,y = nodes, by.x = c("read_start", "read_end"),by.y = c("start", "end"),type = "any",nomatch = 0L)
    		setorder(hits, bin_id, node_order)
    		reads_by_bin = hits[,.(read = list(snode.id),multiplicity = first(multiplicity)),by = bin_id]
		rev_reads_by_bin = reads_by_bin[,.(bin_id,read=lapply(reads_by_bin$read,function(x){-rev(x)}),multiplicity)]
		all_reads = rbind(reads_by_bin,rev_reads_by_bin)[,row:=.I]
		readnames = unlist(lapply(all_reads$read,function(x){paste0(as.character(x),collapse='|')}))
		all_reads$words = readnames
    		out = all_reads[,.(nums = sum(multiplicity)),by = words]
    		setorder(out, -nums, words)
		out
      }))
      reads[,nums:=sum(nums),by=words]
      return(unique(reads))
}

longread_kl = function(walk_x, walk_y, graph=NULL,readL=1e4, depth = 1, background = 1e-5,mc.cores=1) {
	if (is.null(walk_x$graph)){
		if(is.null(graph)){
			error('Must provide either a gWalk object or a graph as input to function')
		}
		walk_x$graph = graph
		walk_y$graph = graph
	}else{
		graph = walk_x$graph
	}
	liktest_separable_lr(walk_x,walk_y,readL=readL,depth=depth,background=background,mc.cores=mc.cores,return.kl=T)
}

hic_kl = function(walk_x, walk_y, target_region=NULL, graph=NULL,pix.size=1e6, depth = 1,theta=2,mask=NULL) {
	if (is.null(walk_x$graph)){
		if(is.null(graph)){
			error('Must provide either a gWalk object or a graph as input to function')
		}
		walk_x$graph = graph
		walk_y$graph = graph
	}else{
		graph = walk_x$graph
	}
	if (is.null(target_region)){target_region=graph$footprint}
	hic_x = forward_simulate(walk_x,target_region = target_region,pix.size=pix.size,depth=depth)
	hic_y = forward_simulate(walk_y,target_region = target_region,pix.size=pix.size,depth=depth)
	kl = compmaps(hic_x,hic_y,theta=theta,return_kl=T,area0=pix.size^2,ifsum=T,mask=mask)
	#kl = kl_nb(hic_x$value,hic_y$value,r=theta)
	return(kl)
}

alignscore = function(x,y,gap_penalty=-1){
	if (is.null(x)){
		return(sum(gap_penalty[as.character(abs(y))]))
	} else if (is.null(y)){
		return(sum(gap_penalty[as.character(abs(x))]))
	}
	xc = as.character(abs(x))
	yc = as.character(abs(y))
	if (length(gap_penalty)==1){
		nodes = unique(c(xc,yc))
		gap_penalty = rep(gap_penalty,length(nodes))
		names(gap_penalty) = nodes
	}
	gx = gap_penalty[xc]
	gy = gap_penalty[yc]
	n = length(x)
	m = length(y)
	S = matrix(0, n + 1, m + 1)
	S[,1] = cumsum(c(0,gx))
	S[1,] = cumsum(c(0,gy))
	for (i in 2:(n+1)){
		for (j in 2:(m+1)){
			match = S[i-1,j-1] + ifelse(x[i-1]==y[j-1],0,gx[i-1]+gy[j-1])
			del = S[i-1,j] + gx[i-1]
			ins = S[i,j-1] + gy[j-1]
			S[i,j] = max(match,del,ins)
		}
	}

	return(S[n + 1, m + 1])
}

alignscore_compl = function(x,y,gp){
	s1 = alignscore(x,y,gp)
	if (!is.null(y)){
		s2 = alignscore(x,-rev(y),gp)
		return(max(s1,s2))}
	else{return(s1)}
}

edit_dist_cpp = function(gwx,gwy,graph=NULL,thresh=0,return_all = F,constpenalty=F){
	if (is.null(graph)){
		graph = gwx$graph
		if (is.null(graph)){
			error('Must provide a graph object or have a graph as an element of gwx')
		}
	}

	nodedt = graph$nodes$dt[,.(width,node.id)]
	widthvec = setNames(graph$nodes$dt$width,graph$nodes$dt$node.id)
	ids = as.integer(names(widthvec))
	max_id = max(abs(ids))
	penalty = numeric(max_id)
	if (constpenalty){
		penalty[abs(ids)] = -1
	}else{
		penalty[abs(ids)] = -widthvec
	}

	sn_x = sort_snodes(gwx$snode.id,gwx$circular)$nodelist
	sn_y = sort_snodes(gwy$snode.id,gwy$circular)$nodelist
	x_totlen = vapply(sn_x,function(ids) sum(widthvec[as.character(abs(ids))]),numeric(1))
	y_totlen = vapply(sn_y,function(ids) sum(widthvec[as.character(abs(ids))]),numeric(1))

	sn_x = sn_x[x_totlen > thresh]
	sn_y = sn_y[y_totlen > thresh]

	n = length(sn_x)
	m = length(sn_y)
	if (n<m){
		sn_x = c(sn_x,rep(list(integer(0)),m-n))
	}else if (m<n){
		sn_y = c(sn_y,rep(list(integer(0)),n-m))
	}
	n = max(n,m)

	costmat = compute_cost_matrix_cpp(sn_x, sn_y, penalty)
	assignment = clue::solve_LSAP(pmax(costmat,0),maximum=F)
	optscores = costmat[cbind(seq_along(assignment),assignment)]
	if (return_all){return(list(assignment,optscores))}else{return(sum(optscores[optscores > thresh]))}
}

edit_dist = function(gwx,gwy,graph=NULL,thresh=0,return_all = F){
	if (is.null(graph)){
		graph = gwx$graph
		if (is.null(graph)){
			error('Must provide a graph object or have a graph as an element of gwx')
		}
	}

	nodedt = graph$nodes$dt[,.(width,node.id)]
	widthvec = setNames(graph$nodes$dt$width,graph$nodes$dt$node.id)
	gap_penalties = -widthvec

	sn_x = sort_snodes(gwx$snode.id,gwx$circular)$nodelist
	sn_y = sort_snodes(gwy$snode.id,gwy$circular)$nodelist
	x_totlen = vapply(sn_x,function(ids) sum(widthvec[as.character(abs(ids))]),numeric(1))
	y_totlen = vapply(sn_y,function(ids) sum(widthvec[as.character(abs(ids))]),numeric(1))

	sn_x = sn_x[x_totlen > thresh]
	sn_y = sn_y[y_totlen > thresh]

	n = length(sn_x)
	m = length(sn_y)
	if (n<m){
		sn_x = c(sn_x,rep(list(NULL),m-n))
	}else if (m<n){
		sn_y = c(sn_y,rep(list(NULL),n-m))
	}
	n = max(n,m)

	comppairs = CJ(i=1:n,j=1:n)
	costs = -unlist(lapply(1:nrow(comppairs),function(x){alignscore_compl(sn_x[[comppairs[x]$i]],sn_y[[comppairs[x]$j]],gap_penalties)}))
	comppairs[,cost:=costs]
	costmat = data.matrix(tidyr::pivot_wider(comppairs,names_from=j,values_from=cost)[,-1])
	assignment = clue::solve_LSAP(costmat,maximum=F)

	optscores = costmat[cbind(seq_along(assignment),assignment)]
	if (return_all){return(list(assignment,optscores))}else{return(sum(optscores[optscores > thresh]))}
}

#just a wrapper function to calculate all distance pairs
#if you want to skip calculating Hi-C (or long-read) distances, set pix.size=0 (or readL=0). To avoid edit distances, set edit_thresh=NA
get_dists = function(gw,graph=NULL,target_region=NULL,pix.size=0,readL=0,edit_thresh=NA,depth=1,mc.cores=1){
	if (is.null(target_region)){target_region = gw[[1]]$footprint}
	comppairs = CJ(i=seq_along(gw),j=seq_along(gw))[i>j]
	if (pix.size>0){
	message('Calculating Hi-C distances')
	hic_dists = unlist(pbmclapply(1:nrow(comppairs),function(x){
		hic_kl(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,pix.size=pix.size,depth=depth,theta=2,target_region=target_region)},mc.cores=mc.cores))
	comppairs[,hic:=hic_dists]
	}
	if (!is.na(edit_thresh)){
	message('Calculating edit distances')
	edit_dists = unlist(pbmclapply(1:nrow(comppairs),function(x){
		edit_dist_cpp(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,thresh=edit_thresh)
				   },mc.cores=mc.cores))
	comppairs[,edit:=edit_dists]
	}
	if (readL>0){
	message('Calculating long-read distances')
	lr_dists = unlist(pbmclapply(1:nrow(comppairs),function(x){
		longread_kl(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,readL=readL,depth=depth)
			   },mc.cores=mc.cores))
	comppairs[,longread:=lr_dists]
	}
	return(comppairs)
}
	
# does something similar to walks() but only returns linear walks on a given graph by sampling and keeping only unique walks. Doesn't return linear walks that exist on no possible decomposition, unlike walks()
get_unique_walks = function(gg,N=1,mode = c('all','lin','circular'),frozen.nodes = NULL,mc.cores=1){
        walksets = sample.gwalks(gg,N,mc.cores=mc.cores,frozen.nodes = frozen.nodes,verbose=F,return.gw=F)
	snodeslist = unlist(mclapply(walksets,function(gw){gw$snode.id},mc.cores=mc.cores),recursive=F)
	circlist = unlist(mclapply(walksets,function(gw){gw$circular},mc.cores=mc.cores))
	if (mode=='lin'){
		keep = circlist==F
	}else if(mode=='circular'){
		keep = circlist==T
	}else{
		keep = rep(T,length(circlist))
	}
	snodeslist = snodeslist[keep]
	circlist = circlist[keep]
	all_walks = gW(graph=gg,snode.id=snodeslist,circular=circlist)
	if (mode=='lin'){
		all_walks = all_walks[all_walks$dt$circular==F]
	}else if(mode=='circular'){
		all_walks = all_walks[all_walks$dt$circular==T]
	}
	hashvec = unlist(mclapply(1:length(circlist),function(x){hash_karyotype_cpp(snodeslist[x],circlist[x])},mc.cores=mc.cores))
	#hashvec = unlist(lapply(1:length(all_walks),function(i){all_walks[i]$hash}))
	all_hashes = data.table(hash=hashvec)[,idx:=.I][,instance:=1:.N,by=hash]
	return(all_walks[all_hashes[instance==1]$idx])
}

get_all_pair_distances = function(gw,depth,pix.size,mc.cores=1){
	comppairs = CJ(i=seq_along(gw),j=seq_along(gw))[i>j]
	hic_dists = unlist(pbmclapply(1:nrow(comppairs),function(x){
			hic_kl(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],pix.size=pix.size,depth=depth,theta=2)
				   },mc.cores=mc.cores))
	return(hic_dists)}


library(data.table)

paste_loose_ends_timed = function(gg,seed=NULL,maxtime=NULL){
	if (is.null(maxtime)){return(paste_loose_ends(gg,seed))}
	result = tryCatch(
	  R.utils::withTimeout({paste_loose_ends(gg,seed)},timeout = maxtime,onTimeout = "error"),
	  TimeoutException = function(e) {return(gg)}
	)
}
paste_loose_ends <- function(gg,seed=NULL){
    work = gr2dt(gg$loose)[cn>0 & terminal==F][,row:=.I][,.(strand,node.cn,cn,index,node.id,orientation,row)]

    leftover = work[0,]
    leftover[, unpaired.cn := integer()]

    if (sum(work$cn) %% 2 == 1) {
        if (!is.null(seed))
            set.seed(seed)
        ## Weighted by CN = uniform sampling over loose-end copies.
        k = sample.int(nrow(work),1,prob = work$cn)
        leftover = work[k]
        leftover[,unpaired.cn:=1]
        work[k, cn := cn - 1]
        work = work[cn > 0]
    }
    #helper functions
    canonical = function(r) {sort(as.integer(r[r > 0]))}
    state_key = function(r) {paste0("s:", paste(r, collapse = ","))}
    lower_bound = function(r) {(length(r) + 1L) %/% 2L}## One junction type can remove at most two active ends.
    greedy_upper_bound = function(r) {
        r = canonical(r)
        n.used = 0
        if (sum(r) %% 2 != 0)stop("Internal error: residual CN sum is odd.")
        while (length(r) >= 2) {
            ## Best easy case: completely paste two equal CNs.
            duplicated.cn = unique(r[duplicated(r)])
            if (length(duplicated.cn)) {
                v = max(duplicated.cn)
                ii = which(r == v)[1:2]
                r = r[-ii]
            } else {
                ## Otherwise completely consume the smallest end
                ## against the largest.
                r[length(r)] = r[length(r)] - r[1]
                r = canonical(r[-1])
            }
            n.used = n.used + 1L
        }
        ## The final residual, if present, must be even and can
        ## be absorbed by one fold-back.
        if (length(r) == 1L) {
            if (r[1L] %% 2L != 0L)
                stop("Internal error: odd residual CN.")
            n.used = n.used + 1L
        }
        n.used
    }
    make_moves = function(r) {
	    r = canonical(r)
	    if (!length(r))
	        return(list())
	    n = length(r)
	    ## Always choose the largest residual-CN end.
	    i = n
	    ri = r[i]
	    moves = list()
	    seen = new.env(hash = TRUE, parent = emptyenv())
	    ## Add a move unless another move has already produced
	    ## the same canonical child state.
	    push_move = function(j, w, foldback, child) {
	        child = canonical(child)
	        key = state_key(child)
	        if (exists(key, envir = seen, inherits = FALSE)){return(NULL)}
	        assign(key, TRUE, envir = seen)
	        list(i = i,j = j,w = as.integer(w),foldback = foldback,child = child)
	    }
	    if (n >= 2L) {
	        ## Only one representative partner is needed for each
	        ## distinct residual CN.
	        partner.cn = unique(r[seq_len(n - 1L)])
	        for (v in partner.cn) {
	            j = which(r == v)[1L]
	            ## Since i is the maximum residual CN:
	            ## min(ri, v) == v.
	            for (w in seq_len(v)) {
	                child = r
	                child[i] = child[i] - w
	                child[j] = child[j] - w
	                move = push_move(j = j,w = w,foldback = FALSE,child = child)
	                if (!is.null(move)){moves[[length(moves) + 1L]] = move}
	            }
	        }
	    }
	    if (ri >= 2L) {
	        for (w in seq_len(ri %/% 2L)) {
	            child = r
	            child[i] = child[i] - 2L * w
	            move = push_move(j = i,w = w,foldback = TRUE,child = child)
	            if (!is.null(move)){moves[[length(moves) + 1L]] = move}
	        }
	    }
	    if (!length(moves)){stop("Internal error: no moves from state ",paste(r, collapse = ",")," (sum = ", sum(r), ")")}
	    score = unlist(lapply(moves,function(a){1L + greedy_upper_bound(a$child)}))
	    active = unlist(lapply(moves,function(a){length(a$child)}))
	    moves[order(score, active)]
	}
        memo <- new.env(hash = TRUE, parent = emptyenv())
	cost.memo = new.env(hash = TRUE,parent = emptyenv())
	move.memo <- new.env(hash = TRUE,parent = emptyenv())
	solve_cost <- function(r) {
	    r <- canonical(r)
	    if (!length(r)){return(0)}
	    if (sum(r) %% 2 != 0)stop("Internal error: odd residual state ",paste(r, collapse = ","))
	    key = state_key(r)
	    if (exists(key,envir = cost.memo,inherits = FALSE)) {
	        return(get(key,envir = cost.memo,inherits = FALSE))
	    }
	    ## A single active end is necessarily handled by
	    ## one fold-back.
	    if (length(r) == 1L) {
	        stopifnot(r[1L] %% 2L == 0L)
	        move <- list(i = 1L,j = 1L,w = r[1L] %/% 2L,foldback = TRUE,child = integer())
	        assign(key,1L,envir = cost.memo)
	        assign(key,move,envir = move.memo)
	        return(1L)
	    }
	    lb = lower_bound(r)
	    best = Inf
	    best.move = NULL
	    for (move in make_moves(r)) {
	        ## Once we have a feasible incumbent, prune states
	        ## that cannot possibly improve it.
	        if (is.finite(best) && 1 + lower_bound(move$child) >= best) {next}
	        cost = 1 + solve_cost(move$child)
	        if (cost < best) {
	            best = cost
	            best.move = move
	            ## The theoretical lower bound has been reached.
	            if (best == lb){break}
	        }
	    }
	    if (is.null(best.move))stop("Internal error: failed to solve state ",paste(r, collapse = ","))
	    best = as.integer(best)
	    assign(key,best,envir = cost.memo)
	    assign(key,best.move,envir = move.memo)
	    best
	}
    initial.state = canonical(work$cn)
    optimal.cost = solve_cost(initial.state)
    residual = work$cn
    moves = vector("list", optimal.cost)
    for (step in seq_len(optimal.cost)){
	active = which(residual>0)
	active = active[order(residual[active],active)]
	r = residual[active]
	chosen = get(state_key(r),envir=move.memo)
	ii = active[chosen$i]
	jj = active[chosen$j]
	if (ii==jj){
		residual[ii] = residual[ii]-2*chosen$w
	}else{
		residual[ii] = residual[ii]-chosen$w
		residual[jj] = residual[jj]-chosen$w
	}
	aa = min(ii,jj)
	bb = max(ii,jj)
	moves[[step]] = data.table(i=work$row[aa],j=work$row[bb],
		loose1=work$index[aa],loose2=work$index[bb],
		cn=chosen$w,foldback=aa==bb)
    }
    if (any(residual != 0L)){stop("Internal error: reconstruction left residual CN.")}
    junctions = if (length(moves)) {rbindlist(moves)
    } else {
        data.table(i = integer(),j = integer(),loose1 = work$index[0],loose2 = work$index[0],cn = integer(),foldback = logical())
    }
    repeated = junctions[,.N,by = .(i, j)][N > 1L]
    if (nrow(repeated)){stop("Internal error: reconstruction repeated a junction pair.")}
    #list(junctions = junctions,leftover = leftover,n.unique.junctions = nrow(junctions),n.states = length(ls(memo, all.names = TRUE)))
    newgg = gg$copy
    loosedt = gr2dt(gg$loose)[cn>0 & terminal==F][,row:=.I][,.(row,node.id,orientation,cn)]
    newedges = junctions[,.(cn,
	n1=loosedt$node.id[i],
	n2=loosedt$node.id[j],
	n1.side=loosedt$orientation[i],
	n2.side=loosedt$orientation[j],
	type='ALT')]
    newgg = gG(nodes=gg$nodes$gr,edges=rbind(gg$edges$dt,newedges,fill=T))
    newgg = loosefix(newgg)
    newgg$set(y.field = "cn")
    return(newgg)
}

boil = function(gg,ft,N,k_return = 1,verbose=F,mc.cores=mc.cores){
	if (is(ft,'character')){ft = streduce(parse.gr(ft))}
	amplicon_nodes = (gg$nodes$gr %&% ft)$node.id %>% unique 
	amplicon_context_nodes = (gg$nodes$gr %&% streduce(ft + 1e3))$node.id %>% unique 
	freeze.nodes = gg$nodes$dt[!(node.id %in% amplicon_context_nodes)]$node.id #we don't want to freeze nodes directly adjacent to the amplicon, hence ft + 1e3
	if(verbose){message('Sampling walks from graph')}
	walks = get_unique_walks(gg,N,mode='circular',frozen.nodes = freeze.nodes,mc.cores=mc.cores) %&% ft
	#get max CN of walks in the graph
	if(verbose){message('Scoring walks')}
	walknodes = walks$nodesdt[,.(walk.id,node.id=abs(snode.id),walk.iid)]
	walknodes[,cn:=.N,by=c('node.id','walk.id')]
	walknodes = merge.data.table(walknodes,gg$nodes$dt[,.(node.id,graph.cn=cn,width)],by=c('node.id'),all.x=T,cartesian=T)
	walknodes[is.na(graph.cn),graph.cn:=0]
	walknodes[,maxN:=graph.cn%/%cn]
	walknodes[,max.walk.cn.nodes:=min(maxN),by=walk.id]
	walknodes[,walk.cn:=max.walk.cn.nodes*cn]
	walkedges = walks$edgesdt[,.(walk.id,edge.id=abs(sedge.id),walk.iid)]
	walkedges[,cn:=.N,by=c('edge.id','walk.id')]
	walkedges = merge.data.table(walkedges,gg$edges$dt[,.(edge.id,graph.cn=cn)],by=c('edge.id'),all.x=T,allow.cartesian=T)
	walkedges[is.na(graph.cn),graph.cn:=0]
	walkedges[,maxN:=graph.cn%/%cn]
	walkedges[,max.walk.cn.edges:=min(maxN),by=walk.id]
	walknodes = unique(merge.data.table(walknodes,walkedges[,.(walk.id,max.walk.cn.edges)],by='walk.id',all.x=T,allow.cartesian=T))
	walknodes[,max.walk.cn:=min(max.walk.cn.nodes,max.walk.cn.edges),by='walk.id']
	#what fraction of the amplicon does the walk account for
	amplicon_weight = sum(gg$nodes$dt[amplicon_nodes]$cn * gg$nodes$dt[amplicon_nodes]$width)
	walknodes[,ampfrac := max.walk.cn*sum(width*(node.id %in% amplicon_nodes)) / amplicon_weight,by=walk.id]
	walkdt = unique(walknodes[,.(walk.id,ampfrac,score=max.walk.cn*ampfrac,max.walk.cn)])
	#what is the entropy of the walk
	#walknodes[,node.entropy := -max.walk.cn*cn/walk.cn*log(cn/walk.cn)]
	#walknodes[,walk.entropy := sum(node.entropy),by=walk.id]
	#score walks combining both these things.. entropy * amplicon fraction perhaps?
	#NOTE: entropy calculation seems bad so far, so we use walk.cn * ampfrac as a score instead
	ord = order(walkdt$score,decreasing=T)
	sorted_walks = walks[walkdt[ord]$walk.id]
	#make an ecDNA solution from the kth best walk
	if(verbose){message(paste0('Genrating top-',k_return,' solutions'))}
	k_solns = mclapply(1:k_return,function(k){
		walk_to_peel = sorted_walks[k]
		peel_cn=walkdt[ord[k]]$max.walk.cn
		#calculate what the remaining node and edge CN in the graph will be
		edges_sub = merge.data.table(unique(walk_to_peel$edgesdt[,.(edge.id=abs(sedge.id))][,.(edge.id,sub_cn = .N*peel_cn),by=edge.id])[,.(edge.id,sub_cn)],gg$edges$dt[,.(edge.id,graph.cn=cn)],by='edge.id')[,.(edge.id,remaining_cn=graph.cn-sub_cn)]
		nodes_sub = merge.data.table(unique(walk_to_peel$nodesdt[,.(node.id=abs(snode.id))][,.(node.id,sub_cn = .N*peel_cn),by=node.id])[,.(node.id,sub_cn)],gg$nodes$dt[,.(node.id,graph.cn=cn)],by='node.id')[,.(node.id,remaining_cn=graph.cn-sub_cn)]
		#instantiate new graph with updated CNs, then sample from it to get the final walks to concatenate to the walks
		gg_out = gg$copy
		gg_out$nodes[nodes_sub$node.id]$mark(cn = nodes_sub$remaining_cn)
		gg_out$edges[edges_sub$edge.id]$mark(cn = edges_sub$remaining_cn)
		gg_out = loosefix(gg_out[,cn>0])
		remaining.walks = sample.gwalks(gg_out,1,verbose=F)[[1]]
		nr = length(remaining.walks)
		combined.gw = gW(graph=gg,snode.id=c(remaining.walks$snode.id,rep(walk_to_peel$snode.id,peel_cn)),circular=c(remaining.walks$circular,rep(T,peel_cn)))
		return(combined.gw)
	},mc.cores=mc.cores)
	return(k_solns)
}

#calculate entropy of the given nodes as distributed in the given gwalk object
amp.entropy = function(gw,amp.nodes){
	amp_nodesdt = gw$nodesdt[,.(walk.id,node.id=abs(snode.id))][node.id %in% amp.nodes]
	amp_nodesdt[,walk.amp := .N,by=walk.id]
	amp_nodesdt[,ampfrac:=walk.amp/nrow(.SD)]
	walk.dist = unique(amp_nodesdt[,.(walk.id,ampfrac)])
	-sum(walk.dist$ampfrac*log(walk.dist$ampfrac))
}

squeeze = function(gg,ft,N,k_return=1,verbose=F,mc.cores=1){
	amplicon_nodes = (gg$nodes$gr %&% ft)$node.id %>% unique 
	amplicon_context_nodes = (gg$nodes$gr %&% streduce(ft + 1e3))$node.id %>% unique 
	freeze.nodes = gg$nodes$dt[!(node.id %in% amplicon_context_nodes)]$node.id #we don't want to freeze nodes directly adjacent to the amplicon, hence ft + 1e3
	if (verbose){message('Sampling walks from graph')}
	gwl = sample.gwalks(gg,N,frozen.nodes = freeze.nodes,verbose=F,mc.cores=mc.cores)
	if (verbose){message('Scoring walks')}
	entropies = unlist(lapply(gwl,function(gw){amp.entropy(gw,amplicon_nodes)}))
	top_walks = gwl[order(entropies)[1:k_return]]
	if(verbose){message(paste0('Genrating top-',k_return,' solutions'))}
	mclapply(top_walks,function(gw){
			 if (sum(gw$circular)>0){embedloops(gw)
			 }else{gw}
			 },mc.cores=mc.cores)
}
