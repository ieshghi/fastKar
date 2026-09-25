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
local.sampling = function(gg,nsteps,nwalk,frozen.nodes=NULL,onlyhash=F,starter_edges=NULL,return.edges=F,return.gw=F){
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
	pivot.node = sample(setdiff(edges[cn>1]$n,frozen.nodes),1)
	new_edges = data.table::copy(edges)
	edges.to.permute = edges[n==pivot.node]$right
	inds = sample(seq_along(edges.to.permute),2)
	edges.to.permute[c(inds[1],inds[2])] <- edges.to.permute[c(inds[2],inds[1])]
  	new_edges[n==pivot.node,right:=edges.to.permute]
	return(new_edges)
  }
  gw0 = traverse_graph_cpp(internal.edges,loose.ends)
  if (return.gw){
	gw0 = gW(graph=gg,snode.id=gw0$snode.id,circular=gw0$circular)
  }else{
  	gw0$hash = hash_karyotype_cpp(gw0$snode.id,gw0$circular)
  }
  walkhist = lapply(1:nwalk,function(i){list(gw0)}) #initialize nwalk walkers at the same point
  hashhist = lapply(1:nwalk,function(i){gw0$hash}) #initialize hashes
  edges = lapply(1:nwalk,function(i){internal.edges}) #initialize edge table
  for (i in seq_len(nwalk)){
  for (j in seq_len(nsteps-1)){
	newedges = permute.node(edges[[i]])	
	newwalk = traverse_graph_cpp(newedges,loose.ends)
	if (return.gw){
		newwalk = gW(graph=gg,snode.id=newwalk$snode.id,circular=newwalk$circular)
		newhash = newwalk$hash
	}
	else{
		newhash = hash_karyotype_cpp(newwalk$snode.id,newwalk$circular)
		newwalk$hash = newhash
	}
	edges[[i]] = newedges
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

#faster version of local.sampling for sampling the neighborhood of one walk. It only samples walks that are exactly one step away from the starting point
sample.neighborhood = function(gg,n_neighbor,starter_edges=NULL,return.edges=F,return.gw=T,mc.cores=1){
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
  if (return.gw){
	gw0 = gW(graph=gg,snode.id=gw0$snode.id,circular=gw0$circular)
  }else{
  	gw0$hash = hash_karyotype_cpp(gw0$snode.id,gw0$circular)
  }
  new_perms = mclapply(1:n_neighbor,function(x){permute.node(copy(internal.edges))},mc.cores=mc.cores)
  neighbor_walks = mclapply(new_perms,function(p){
			newwalk = traverse_graph_cpp(p,loose.ends)
			if (return.gw){
				newwalk = gW(graph=gg,snode.id=newwalk$snode.id,circular=newwalk$circular)
				newhash = newwalk$hash
			}
			else{
				newhash = hash_karyotype_cpp(newwalk$snode.id,newwalk$circular)
				newwalk$hash = newhash
			}
			return(newwalk)
		},mc.cores=mc.cores)
  if (return.edges){
	return(list(orig_walk = gw0,neighbor_walks=neighbor_walks,starter_edges = internal.edges))
  } else{
	return(list(orig_walk = gw0,neighbor_walks=neighbor_walks))
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
  gw0$hash = hash_karyotype_cpp(gw0$snode.id,gw0$circular)
  walkhist[[1]] = gw0
  hashhist = c(hashhist,gw0$hash)
  for (i in seq_len(len-1)){
	valid = F
  	killwalk = T #gets switched to F if we find the next step
  	for (j in seq_len(attempts)){
		if (!valid){
			edge_try = permute.node(internal.edges)	
			walk_try = traverse_graph_cpp(edge_try,loose.ends)
			hash_try = hash_karyotype_cpp(walk_try$snode.id,walk_try$circular)
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

sort_snodes = function(nodelist,circ=NULL) {
	# Fast C++ canonicalization using the same Booth-rotation implementation
	# and lexicographic ordering as hash_karyotype_cpp().
	if (is.null(circ)) {
		return(sort_snodes_cpp(nodelist, rep(FALSE, length(nodelist))))
	}
	return(sort_snodes_cpp(nodelist, circ))
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

reads_fromwalk = function(walk,readL,minsize=0,use.nodes=NULL){
	gr = walk$graph$gr[,c('node.id')]
	tinynodes = gr[width(gr) < minsize]$node.id
	if (!is.null(use.nodes)){
		tinynodes = c(tinynodes,gr$node.id[!(gr$node.id %in% use.nodes)])
	}
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
		readnames = unlist(lapply(all_reads$read,function(x){
						  x = x[!(abs(x) %in% tinynodes)]
						  paste0(as.character(x),collapse='|')
						}))
		all_reads$words = readnames
    		out = all_reads[,.(nums = sum(multiplicity)),by = words]
    		setorder(out, -nums, words)
		out
      }))
      reads[,nums:=sum(nums),by=words]
      return(unique(reads))
}

longread_kl = function(walk_x, walk_y, graph=NULL,readL=1e4, depth = 1, background = 1e-5,use.nodes=use.nodes,mc.cores=1) {
	if (is.null(walk_x$graph)){
		if(is.null(graph)){
			error('Must provide either a gWalk object or a graph as input to function')
		}
		walk_x$graph = graph
		walk_y$graph = graph
	}else{
		graph = walk_x$graph
	}
	liktest_separable_lr(walk_x,walk_y,readL=readL,depth=depth,background=background,mc.cores=mc.cores,use.nodes=use.nodes,return.kl=T)
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

edit_dist_cpp = function(gwx,gwy,graph=NULL,return_all = F,constpenalty=F,use.nodes=NULL){
	if (is.null(graph)){
		graph = gwx$graph
		if (is.null(graph)){
			error('Must provide a graph object or have a graph as an element of gwx')
		}
	}
	nodedt = graph$nodes$dt[,.(width,node.id)]
	if (!is.null(use.nodes)){
		nodedt[!(node.id %in% use.nodes)]$width=0
	}
	widthvec = setNames(nodedt$width,nodedt$node.id)
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
	if (return_all){return(list(assignment,optscores))}else{return(sum(optscores))}
}

#just a wrapper function to calculate all distance pairs
#if you want to skip calculating Hi-C (or long-read) distances, set pix.size=0 (or readL=0). To avoid edit distances, set edit_thresh=NA
get_dists = function(gw,graph=NULL,target_region=NULL,pix.size=0,readL=0,edit_thresh=NA,use.nodes=NULL,depth=1,mc.cores=1){
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
		edit_dist_cpp(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,use.nodes=use.nodes)
				   },mc.cores=mc.cores))
	comppairs[,edit:=edit_dists]
	}
	if (readL>0){
	message('Calculating long-read distances')
	lr_dists = unlist(pbmclapply(1:nrow(comppairs),function(x){
		longread_kl(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,readL=readL,depth=depth,use.nodes=use.nodes)
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
	  TimeoutException = function(e) {stop(e)}
	)
}
paste_loose_ends <- function(gg, seed = NULL) {
    work <- gr2dt(gg$loose)[cn > 0 & terminal == FALSE][, row := .I][, .(strand, node.cn, cn, index, node.id, orientation, row)]
    if (!nrow(work))
        return(gg$copy)
    work[, cn := as.integer(cn)]
    ## If the total CN is odd, one copy necessarily remains unpaired.
    leftover <- work[0]
    leftover[, unpaired.cn := integer()]
    if (sum(work$cn) %% 2L == 1L) {
        if (!is.null(seed))
            set.seed(seed)
        ## Sampling proportional to CN is equivalent to selecting
        ## a uniformly random loose-end copy.
        k <- sample.int(nrow(work), 1L, prob = work$cn)
        leftover <- copy(work[k])
        leftover[, unpaired.cn := 1L]
        work[k, cn := cn - 1L]
    }
    residual <- work$cn
    moves <- list()
    add_move <- function(i, j, cn) {
        moves[[length(moves) + 1L]] <<- data.table(i = i,j = j,loose1 = work$index[i],loose2 = work$index[j],cn = as.integer(cn),foldback = i == j)
    }
    while (any(residual > 0L)) {
        active <- which(residual > 0L)
        ## Only one end remains: close it with a foldback.
        if (length(active) == 1L) {
            i <- active[1L]
            if (residual[i] %% 2L != 0L) {
                stop(
                    "Cannot consume the final loose end: residual CN is odd."
                )
            }
            add_move(i, i, residual[i] %/% 2L)
            residual[i] <- 0L
            next
        }
        active.cn <- residual[active]
        ## Prefer completely matching two equal-CN ends.
        duplicated.cn <- unique(active.cn[duplicated(active.cn)])
        if (length(duplicated.cn)) {
            ## Prefer the largest exact match.
            target.cn <- max(duplicated.cn)
            pair <- active[active.cn == target.cn][1:2]
            i <- pair[1L]
            j <- pair[2L]
            w <- target.cn
        } else {
            ## Otherwise pair the two largest residual ends.
            pair <- active[order(residual[active], decreasing = TRUE)][1:2]
            i <- pair[1L]
            j <- pair[2L]
            w <- min(residual[i], residual[j])
        }
        add_move(i, j, w)
        residual[i] <- residual[i] - w
        residual[j] <- residual[j] - w
    }
    junctions <- if (length(moves)) {rbindlist(moves)
    } else {
        data.table(i = integer(),j = integer(),loose1 = work$index[0],loose2 = work$index[0],cn = integer(),foldback = logical())
    }
    ## Every non-foldback move exhausts at least one endpoint,
    ## so the same pair should never be generated twice.
    repeated <- junctions[, .N, by = .(pmin(i, j), pmax(i, j))][N > 1L]
    if (nrow(repeated))
        stop("Internal error: a junction pair was generated more than once.")
    newedges <- junctions[, .(cn,n1 = work$node.id[i],n2= work$node.id[j],n1.side = work$orientation[i],n2.side = work$orientation[j],type    = "ALT")]
    newgg <- gG(nodes = gg$nodes$gr,edges = rbind(gg$edges$dt, newedges, fill = TRUE))
    newgg <- loosefix(newgg)
    newgg$set(y.field = "cn")
    ## Attach these if you want to inspect what the heuristic did.
    newgg
}


