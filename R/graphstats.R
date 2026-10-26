#' @useDynLib fastKar, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#goes from a ggraph to a "wiring", which gives all the internal edges of the graph (going from left side of a node to right side)
#along with the loose node ids and a reference data table with the new node ids (all copies of each node are de-duplicated)
## Output is identical to the reference (fastKar 0a5497e) implementation, which encoded each edge
## end as the string paste0(node, substr(side, 1, 1)), melted, and recovered left/right ends with
## grepl/substr. Here the ends are kept as integer node ids plus a side code:
##  * edge copies: external edges (graph edges, then left loose ends, then right loose ends, as in
##    the reference rbind) with cn > 0, each repeated cn times; subid = copy index;
##  * the melted long table lists every copy's n1 end, then every copy's n2 end;
##  * left (right) table = long rows whose side starts with 'l' ('r') and whose node is > 0,
##    stably sorted by node (as setkeyv did);
##  * loose.ends = subids of copies of type 'LOO', ascending (= unique() over the melted order).
gg.to.wiring = function(gg){ 
  if(!(('loose.cn.left' %in% colnames(gg$nodes$dt)) & ('loose.cn.right' %in% colnames(gg$nodes$dt)))){
  	gg = loosefix(gg)
  }
  nodesdt = gg$nodes$dt[,.(start,end,seqnames,snode.id,cn,loose.cn.left,loose.cn.right)]
  edgesdt = gg$edges$dt[,.(cn,n1,n2,n1.side,n2.side,type)]
  ## same NA fill as the reference (data.table coercion rules included)
  nodesdt[is.na(loose.cn.right),loose.cn.right:=cn]
  nodesdt[is.na(loose.cn.left),loose.cn.left:=cn] #adjusting CN at chromosome ends
  lsel = which(nodesdt$loose.cn.left > 0)
  rsel = which(nodesdt$loose.cn.right > 0)
  nL = length(lsel); nR = length(rsel)
  ## left loose ends: n1 = 0 (right side), n2 = node (left side)
  ## right loose ends: n1 = node (right side), n2 = 0 (left side)
  ext.cn = c(edgesdt$cn, nodesdt$loose.cn.left[lsel], nodesdt$loose.cn.right[rsel])
  ext.n1 = c(edgesdt$n1, rep(0, nL), nodesdt$snode.id[rsel])
  ext.n2 = c(edgesdt$n2, nodesdt$snode.id[lsel], rep(0, nR))
  ext.s1 = c(substr(edgesdt$n1.side, 1, 1), rep('r', nL + nR))
  ext.s2 = c(substr(edgesdt$n2.side, 1, 1), rep('l', nL + nR))
  ext.loo = c(as.character(edgesdt$type), rep("LOO", nL + nR)) == "LOO"
  keep = which(ext.cn > 0)
  copies = keep[rep(seq_along(keep), ext.cn[keep])] #separate all copies of all edges
  ncp = length(copies)
  node.long = suppressWarnings(c(as.integer(ext.n1[copies]), as.integer(ext.n2[copies])))
  side.long = c(ext.s1[copies], ext.s2[copies])
  subid.long = c(seq_len(ncp), seq_len(ncp))
  lrows = which(side.long == 'l' & node.long > 0)
  rrows = which(side.long == 'r' & node.long > 0)
  lrows = lrows[order(node.long[lrows], method = 'radix')]
  rrows = rrows[order(node.long[rrows], method = 'radix')]
  # each copy of each node has two edges coming out of it, one on the left one on the right. Each of those edges has a unique edge id which we call "subid"
  internal.edges = data.table(n = node.long[lrows], left = subid.long[lrows], right = subid.long[rrows])
  setattr(internal.edges, 'sorted', 'n')
  internal.edges[,cn:=.N,by=n] #each unique mapping of right-edged to left-edges (reorderings of the third column) corresponds to a unique walk decomposition of the graph
  # where are the loose edges?
  loose.subids = which(ext.loo[copies])
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
	  # full chunks plus a last partial one, so exactly N permutations are drawn (multiples of chunksize are unchanged)
	  sizes = c(rep(chunksize, chunks), if (N %% chunksize) N %% chunksize)
	  if(verbose){
  	perms <- do.call('c',pbmclapply(sizes, function(K) shuffle_chunk(K, internal.edges),mc.cores = mc.cores))
  }else{
  	perms <- do.call('c',mclapply(sizes, function(K) shuffle_chunk(K, internal.edges),mc.cores = mc.cores))
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

## Integer view of a wiring table (internal.edges) for the samplers below. Only `right` changes
## between karyotypes, so the samplers carry a `right` vector instead of copying the table, and
## look up the rows of a node through `rows` (rows of each distinct n, in row order: exactly what
## the reference selected with edges[n == x]).
.wiring_state = function(edges, loose.ends) {
  n = edges$n
  un = unique(n)
  list(n = n, cn = edges$cn, right = edges$right, un = un,
       rows = unname(split(seq_along(n), factor(match(n, un), levels = seq_along(un)))),
       n_int = as.integer(n), left_int = as.integer(edges$left), id_int = as.integer(edges$id),
       loose = as.numeric(loose.ends))
}
.wiring_rows = function(st, node) st$rows[[match(node, st$un)]]
## traverse_graph_cpp(edges with this right column, loose.ends)
.wiring_traverse = function(st, right) traverse_graph_fast_cpp(st$n_int, st$left_int, as.integer(right), st$id_int, st$loose)
## the reference permute.node swap: same sample() call, same assignment
.wiring_swap = function(right, rows) {
	edges.to.permute = right[rows]
	inds = sample(seq_along(edges.to.permute),2)
	edges.to.permute[c(inds[1],inds[2])] <- edges.to.permute[c(inds[2],inds[1])]
	right[rows] = edges.to.permute
	right
}
## the reference initial shuffle (in place on `edges`), kept verbatim for its RNG stream
.wiring_shuffle = function(edges) {
        new_right = edges[, if (.N > 1) sample(right, .N) else right, by = n]$V1
  	return(edges[,right:=new_right])
}
## internal.edges of a wiring: the wiring's own table when it was computed here (the reference
## shuffles that in place), a copy when the caller passed one in, so the caller's object is kept.
.wiring_edges = function(wiring, supplied) {
  if (supplied) data.table::copy(wiring$internal.edges) else wiring$internal.edges
}

#samples karyotype space using markov chains always starting at the same walk, initialized randomly. 
## wiring: optional precomputed gg.to.wiring(gg) (not modified). Same RNG draws and outputs as the
## reference: the by-group shuffle of `right` (unless starter_edges), then per step
## sample.int(#pivots, 1) and sample(seq_along(rows), 2).
local.sampling = function(gg,nsteps,nwalk,frozen.nodes=NULL,onlyhash=F,starter_edges=NULL,return.edges=F,return.gw=F,wiring=NULL){
  supplied = !is.null(wiring)
  if (!supplied) wiring = gg.to.wiring(gg)
  if (is.null(starter_edges)){
  	internal.edges = .wiring_shuffle(.wiring_edges(wiring, supplied))
  } else {
	internal.edges = starter_edges
  }
  loose.ends = wiring$loose.ends
  st = .wiring_state(internal.edges, loose.ends)
  ## n and cn never change between steps, so the pivots (and their rows) are fixed
  pivots = setdiff(st$n[which(st$cn > 1)], frozen.nodes)
  pivot.rows = lapply(pivots, function(p) .wiring_rows(st, p))
  permute.node = function(right) {
  	if (length(pivots)==0){return(right)}
	k = sample.int(length(pivots),1) #sample(x,1) on a single integer x would draw from 1:x
	.wiring_swap(right, pivot.rows[[k]])
  }
  ## the reference's first step evaluates edges[n == pivot] on internal.edges itself, which makes
  ## data.table add an auto-index when the table is not keyed by n; repeat it so starter_edges
  ## comes back with the same attributes
  if (nwalk >= 1 && nsteps >= 2 && length(pivots) > 0) invisible(internal.edges[n == pivots[1L]])
  gw0 = .wiring_traverse(st, st$right)
  if (return.gw){
	gw0 = gW(graph=gg,snode.id=gw0$snode.id,circular=gw0$circular)
  }else{
  	gw0$hash = hash_karyotype_cpp(gw0$snode.id,gw0$circular)
  }
  walkhist = lapply(1:nwalk,function(i){list(gw0)}) #initialize nwalk walkers at the same point
  hashhist = lapply(1:nwalk,function(i){gw0$hash}) #initialize hashes
  for (i in seq_len(nwalk)){
  	steps = seq_len(nsteps-1)
	if (!length(steps)) next
	right = st$right
	wl = c(walkhist[[i]], vector('list', length(steps)))
	hl = c(hashhist[[i]], character(length(steps)))
  	for (j in steps){
		right = permute.node(right)
		newwalk = .wiring_traverse(st, right)
		if (return.gw){
			newwalk = gW(graph=gg,snode.id=newwalk$snode.id,circular=newwalk$circular)
			newhash = newwalk$hash
		}
		else{
			newhash = hash_karyotype_cpp(newwalk$snode.id,newwalk$circular)
			newwalk$hash = newhash
		}
		wl[[j+1]] = newwalk 
		hl[j+1] = newhash
  	}
	walkhist[[i]] = wl
	hashhist[[i]] = hl
  }
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
## wiring: optional precomputed gg.to.wiring(gg) (not modified). Draws as in the reference:
## sample(edges[cn > 1]$n, 1) (pivot weighted by CN), then sample(seq_along(rows), 2).
sample.neighborhood = function(gg,n_neighbor,starter_edges=NULL,return.edges=F,return.gw=T,mc.cores=1,wiring=NULL){
  supplied = !is.null(wiring)
  if (!supplied) wiring = gg.to.wiring(gg)
  if (is.null(starter_edges)){
  	internal.edges = .wiring_shuffle(.wiring_edges(wiring, supplied))
  } else {
	internal.edges = starter_edges
  }
  loose.ends = wiring$loose.ends
  st = .wiring_state(internal.edges, loose.ends)
  cand = st$n[which(st$cn > 1)]
  permute.node = function(right) {
  	if (length(cand)==0){return(right)}
	pivot.node = sample(cand,1)
	.wiring_swap(right, .wiring_rows(st, pivot.node))
  }
  gw0 = .wiring_traverse(st, st$right)
  if (return.gw){
	gw0 = gW(graph=gg,snode.id=gw0$snode.id,circular=gw0$circular)
  }else{
  	gw0$hash = hash_karyotype_cpp(gw0$snode.id,gw0$circular)
  }
  new_perms = mclapply(1:n_neighbor,function(x){permute.node(st$right)},mc.cores=mc.cores)
  neighbor_walks = mclapply(new_perms,function(p){
			newwalk = .wiring_traverse(st, p)
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
## wiring: optional precomputed gg.to.wiring(gg) (not modified).
markov.gwalk = function(gg,len,self.avoid = F,attempts = 10,return.gw=F,seed=NULL,wiring=NULL){
  if (is.null(seed)){set.seed(sample(1:1e5,1))
  } else{set.seed(seed)}
  supplied = !is.null(wiring)
  if (!supplied) wiring = gg.to.wiring(gg)
  internal.edges = .wiring_edges(wiring, supplied)
  loose.ends = wiring$loose.ends
  hashhist = c()
  walkhist = list()
  internal.edges = .wiring_shuffle(internal.edges)
  st = .wiring_state(internal.edges, loose.ends)
  cand = st$n[which(st$cn > 1)]
  permute.node = function(right) {
	pivot.node = sample(cand,1)
	.wiring_swap(right, .wiring_rows(st, pivot.node))
  }
  right = st$right
  gw0 = .wiring_traverse(st, right)
  gw0$hash = hash_karyotype_cpp(gw0$snode.id,gw0$circular)
  walkhist[[1]] = gw0
  hashhist = c(hashhist,gw0$hash)
  for (i in seq_len(len-1)){
	valid = F
  	killwalk = T #gets switched to F if we find the next step
  	for (j in seq_len(attempts)){
		if (!valid){
			edge_try = permute.node(right)	
			walk_try = .wiring_traverse(st, edge_try)
			hash_try = hash_karyotype_cpp(walk_try$snode.id,walk_try$circular)
			if (self.avoid){
				valid = !(hash_try %in% hashhist)
			} else{valid=T}
		} else{
			right = edge_try
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
## The wiring is computed once and shared by all chains (gg.to.wiring draws no random numbers).
multi.markov = function(gg,N,len,self.avoid=F,attempts=10,mc.cores=1,return.gw=F,seed=NULL,wiring=NULL){
	if (is.null(seed)){set.seed(sample(1:1e5,1))
	}else{set.seed(seed)}
	seeds = sample(1:1e5,N)
	if (is.null(wiring)) wiring = gg.to.wiring(gg)
	mclapply(1:N,function(i){
		markov.gwalk(gg=gg,len=len,self.avoid=self.avoid,attempts=attempts,return.gw=return.gw,seed=seeds[i],wiring=wiring)  
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

## Long-read words of a walk set: one row per distinct word (node ids of the reads' overlapped
## nodes joined by '|', both orientations, nodes narrower than minsize or outside use.nodes
## dropped) with nums = number of read start positions giving it, summed over walks.
## Identical (words, nums and row order) to the reference foverlaps implementation; the sweep is
## in reads_fromwalk_cpp (src/graphstats_reads.cpp).
reads_fromwalk = function(walk,readL,minsize=0,use.nodes=NULL){
	if (length(.reads_cache$entries)){
		hit = .reads_cache_get(walk,readL,minsize,use.nodes)
		if (!is.null(hit)) return(data.table::copy(hit))
	}
	gr = walk$graph$gr[,c('node.id')]
	tinynodes = gr[width(gr) < minsize]$node.id
	if (!is.null(use.nodes)){
		tinynodes = c(tinynodes,gr$node.id[!(gr$node.id %in% use.nodes)])
	}
	sn = walk$snode.id
	lens = lengths(sn)
	circ = as.logical(walk$circular)
	if (!(is.numeric(readL) && length(readL) == 1 && is.finite(readL) && readL >= 1 && readL == round(readL)) ||
	    !length(sn) || any(lens == 0) || length(circ) != length(sn) || anyNA(circ)){
		return(.reads_fromwalk_ref(walk,readL,minsize=minsize,use.nodes=use.nodes))
	}
	x = unlist(sn, use.names = FALSE)
	ax = abs(x)
	## tokens per walk, so each walk's ids print as as.character() prints that walk's vector
	fwd = unlist(lapply(sn, function(v) as.character(v)), use.names = FALSE)
	rtok = unlist(lapply(sn, function(v) as.character(-v)), use.names = FALSE)
	res = reads_fromwalk_cpp(width(gr)[ax], as.integer(lens), circ, as.numeric(readL), fwd, rtok, !(ax %in% tinynodes))
	data.table(words = res$words, nums = res$nums)
}

## Reference (fastKar 0a5497e) implementation, kept verbatim as the fallback of reads_fromwalk()
## for read lengths that are not integer-valued (where foverlaps' end adjustment on doubles could
## matter) and for degenerate walks.
.reads_fromwalk_ref = function(walk,readL,minsize=0,use.nodes=NULL){
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
			stop('Must provide either a gWalk object or a graph as input to function')
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
			stop('Must provide either a gWalk object or a graph as input to function')
		}
		walk_x$graph = graph
		walk_y$graph = graph
	}else{
		graph = walk_x$graph
	}
	if (is.null(target_region)){target_region=graph$footprint}
	hic_x = forward_simulate(walk_x,target_region = target_region,pix.size=pix.size,depth=depth)
	hic_y = forward_simulate(walk_y,target_region = target_region,pix.size=pix.size,depth=depth)
	kl = .hic_kl_maps(hic_x,hic_y,theta=theta,pix.size=pix.size,mask=mask)
	#kl = kl_nb(hic_x$value,hic_y$value,r=theta)
	return(kl)
}
## the comparison step of hic_kl, on already simulated maps
.hic_kl_maps = function(hic_x,hic_y,theta=2,pix.size,mask=NULL){
	compmaps(hic_x,hic_y,theta=theta,return_kl=T,area0=pix.size^2,ifsum=T,mask=mask)
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

## Cut sorted walks into fragments of nonzero-penalty nodes that stay inside a
## region. left.in/right.in (indexed by node id) say whether each node's left
## and right end lies in the region. Consecutive kept nodes share a fragment
## only if no dropped node lies between them and the walk leaves the first and
## enters the second through in-region ends. For a circular walk, the last
## fragment continues into the first when the wrap-around step qualifies.
split_walk_fragments = function(nodelist,circ,penalty,left.in,right.in){
	exit.in = function(s) ifelse(s > 0, right.in[abs(s)], left.in[abs(s)])
	entry.in = function(s) ifelse(s > 0, left.in[abs(s)], right.in[abs(s)])
	## pieces of each walk are collected per walk and concatenated once at the end (the reference
	## grew `out = c(out, pieces)` inside the loop, quadratic in the number of walks)
	res = vector('list', length(nodelist))
	for (k in seq_along(nodelist)){
		v = nodelist[[k]]
		idx = which(penalty[abs(v)] != 0)
		if (!length(idx)) next
		kv = v[idx]
		nk = length(kv)
		joined = if (nk > 1) diff(idx) == 1 & exit.in(kv[-nk]) & entry.in(kv[-1]) else logical(0)
		pieces = unname(split(kv, cumsum(c(TRUE, !joined))))
		np = length(pieces)
		if (isTRUE(circ[k]) && np > 1 && idx[1] == 1 && idx[nk] == length(v) &&
		    exit.in(kv[nk]) && entry.in(kv[1])){
			pieces[[1]] = c(pieces[[np]], pieces[[1]])
			pieces = pieces[-np]
		}
		res[[k]] = pieces
	}
	out = unlist(res, recursive = FALSE, use.names = FALSE)
	if (is.null(out)) list() else out
}

## node.widths: optional per-node weights (named by node.id) replacing node widths;
## nodes not named get weight 0.
## split.ends: optional list(left = , right = ) of logicals named by node.id,
## whether each node's left/right end lies in the comparison region. When given,
## walks are compared as fragments inside that region (see split_walk_fragments)
## rather than as whole walks. Both NULL reproduces the original behaviour.
edit_dist_cpp = function(gwx,gwy,graph=NULL,return_all = F,constpenalty=F,use.nodes=NULL,node.widths=NULL,split.ends=NULL){
	if (is.null(graph)){
		graph = gwx$graph
		if (is.null(graph)){
			stop('Must provide a graph object or have a graph as an element of gwx')
		}
	}
	penalty = .edit_penalty(graph,constpenalty=constpenalty,use.nodes=use.nodes,node.widths=node.widths)
	sx = sort_snodes(gwx$snode.id,gwx$circular)
	sy = sort_snodes(gwy$snode.id,gwy$circular)
	.edit_pair(.edit_pieces(sx,penalty,split.ends),.edit_pieces(sy,penalty,split.ends),penalty,return_all)
}

## edit_dist_cpp in three steps, so callers comparing many walks (get_dists) can sort each walk
## once. Same code as the reference body, split at the walk-independent / per-walk boundaries.
.edit_penalty = function(graph,constpenalty=F,use.nodes=NULL,node.widths=NULL){
	nodedt = graph$nodes$dt[,.(width,node.id)]
	if (!is.null(node.widths)){
		w = node.widths[as.character(nodedt$node.id)]
		w[is.na(w)] = 0
		nodedt$width = as.numeric(w)
	}
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
	penalty
}
## sorted walks (sort_snodes output) -> the pieces that are aligned
.edit_pieces = function(sx,penalty,split.ends=NULL){
	if (!is.null(split.ends)){
		max_id = length(penalty)
		left.in = right.in = logical(max_id)
		left.in[as.integer(names(split.ends$left))] = split.ends$left
		right.in[as.integer(names(split.ends$right))] = split.ends$right
		split_walk_fragments(sx$nodelist,sx$circ,penalty,left.in,right.in)
	}else{
		## Nodes with zero edit penalty cannot affect alignment scores.  In
		## footprint-local comparisons on full graphs, most nodes are often outside
		## use.nodes and have zero penalty, so strip them before the pairwise DP cost
		## matrix calculation.
		lapply(sx$nodelist,function(x){x[penalty[abs(x)]!=0]})
	}
}
.edit_pair = function(sn_x,sn_y,penalty,return_all=F){
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

## Scoped cache of reads_fromwalk() results, keyed by the walk object (identical(): for gWalks,
## the same R6 object) and the arguments. Used by get_dists() so longread_kl() on every pair
## reuses each walk's reads instead of recomputing them; outside .with_reads_cache() it is empty.
.reads_cache = new.env(parent = emptyenv())
.reads_cache_get = function(walk, readL, minsize, use.nodes){
	for (e in .reads_cache$entries){
		if (identical(e$walk, walk) && identical(e$readL, readL) && identical(e$minsize, minsize) &&
		    identical(e$use.nodes, use.nodes)) return(e$reads)
	}
	NULL
}
.with_reads_cache = function(walks, readL, minsize, use.nodes, expr){
	old = .reads_cache$entries
	on.exit(.reads_cache$entries <- old, add = TRUE)
	.reads_cache$entries = lapply(walks, function(w){
		list(walk = w, readL = readL, minsize = minsize, use.nodes = use.nodes,
		     reads = reads_fromwalk(w, readL, minsize = minsize, use.nodes = use.nodes))
	})
	expr
}
.all_have_graph = function(gw) all(vapply(gw, function(w) !is.null(w$graph), logical(1)))

#just a wrapper function to calculate all distance pairs
#if you want to skip calculating Hi-C (or long-read) distances, set pix.size=0 (or readL=0).
## Same columns as the reference, which re-simulated both karyotypes of every pair (n(n-1) calls
## to forward_simulate) and re-sorted / re-read both walks of every pair. Here, when every walk
## carries its graph (so hic_kl/longread_kl would not substitute `graph`), each karyotype is
## simulated once, sorted once and read once, and the pairwise comparisons reuse those.
get_dists = function(gw,graph=NULL,target_region=NULL,pix.size=0,readL=0,use.nodes=NULL,depth=1,mc.cores=1){
	if (is.null(target_region)){target_region = gw[[1]]$footprint}
	comppairs = CJ(i=seq_along(gw),j=seq_along(gw))[i>j]
	fast = nrow(comppairs) > 0 && .all_have_graph(gw)
	if (pix.size>0){
	if (fast){
		sims = mclapply(seq_along(gw),function(k){
			forward_simulate(gw[[k]],target_region = target_region,pix.size=pix.size,depth=depth)},mc.cores=mc.cores)
		hic_dists = unlist(mclapply(1:nrow(comppairs),function(x){
			.hic_kl_maps(sims[[comppairs$i[x]]],sims[[comppairs$j[x]]],theta=2,pix.size=pix.size)},mc.cores=mc.cores))
		rm(sims)
	}else{
	hic_dists = unlist(mclapply(1:nrow(comppairs),function(x){
		hic_kl(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,pix.size=pix.size,depth=depth,theta=2,target_region=target_region)},mc.cores=mc.cores))
	}
	comppairs[,hic:=hic_dists]
	}
	if (fast){
		## edit_dist_cpp(gw_i, gw_j, graph): penalty from `graph`, or from gw_i's graph
		pens = if (!is.null(graph)) rep(list(.edit_penalty(graph,use.nodes=use.nodes)),length(gw)) else
			lapply(gw,function(w) .edit_penalty(w$graph,use.nodes=use.nodes))
		sorted = lapply(gw,function(w) sort_snodes(w$snode.id,w$circular))
		edit_dists = unlist(mclapply(1:nrow(comppairs),function(x){
			i = comppairs$i[x]; j = comppairs$j[x]; pen = pens[[i]]
			.edit_pair(.edit_pieces(sorted[[i]],pen),.edit_pieces(sorted[[j]],pen),pen)
					   },mc.cores=mc.cores))
	}else{
	edit_dists = unlist(mclapply(1:nrow(comppairs),function(x){
		edit_dist_cpp(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,use.nodes=use.nodes)
				   },mc.cores=mc.cores))
	}
	comppairs[,edit:=edit_dists]
	if (readL>0){
	lr_pairs = function(){
		unlist(mclapply(1:nrow(comppairs),function(x){
		longread_kl(gw[[comppairs[x]$i]],gw[[comppairs[x]$j]],graph=graph,readL=readL,depth=depth,use.nodes=use.nodes)
			   },mc.cores=mc.cores))
	}
	## longread_kl -> liktest_separable_lr -> longread_probdist -> reads_fromwalk(w, readL, 0, use.nodes)
	lr_dists = if (fast) .with_reads_cache(gw,readL,0,use.nodes,lr_pairs()) else lr_pairs()
	comppairs[,longread:=lr_dists]
	}
	return(comppairs)
}
	
# does something similar to walks() but only returns linear walks on a given graph by sampling and keeping only unique walks. Doesn't return linear walks that exist on no possible decomposition, unlike walks()
get_unique_walks = function(gg,N=1,mode = c('all','lin','circular'),frozen.nodes = NULL,mc.cores=1){
        walksets = sample.gwalks(gg,N,mc.cores=mc.cores,frozen.nodes = frozen.nodes,verbose=F,return.gw=F)
	snodeslist = unlist(lapply(walksets,'[[','snode.id'),recursive=F)
	circlist = unname(unlist(lapply(walksets,'[[','circular')))
	if (mode=='lin'){
		keep = circlist==F
	}else if(mode=='circular'){
		keep = circlist==T
	}else{
		keep = rep(T,length(circlist))
	}
	snodeslist = snodeslist[keep]
	circlist = circlist[keep]
	if (length(snodeslist)){
		hashvec = unlist(mclapply(1:length(circlist),function(x){hash_karyotype_cpp(snodeslist[x],circlist[x])},mc.cores=mc.cores))
		unq = which(!duplicated(hashvec))
		all_walks = gW(graph=gg,snode.id=snodeslist[unq],circular=circlist[unq])
		return(all_walks)
	}else{
		return(gW())
	}
}

## Same values as the reference, which called hic_kl on every pair (target = gw_i's graph
## footprint). When every walk has a graph and all footprints are identical, each karyotype is
## simulated once on that footprint and the pairs reuse the maps.
get_all_pair_distances = function(gw,depth,pix.size,mc.cores=1){
	comppairs = CJ(i=seq_along(gw),j=seq_along(gw))[i>j]
	fast = nrow(comppairs) > 0 && .all_have_graph(gw)
	if (fast){
		fps = lapply(gw,function(w) w$graph$footprint)
		fast = all(vapply(fps,identical,logical(1),fps[[1]]))
	}
	if (fast){
		sims = mclapply(seq_along(gw),function(k){
			forward_simulate(gw[[k]],target_region = fps[[1]],pix.size=pix.size,depth=depth)},mc.cores=mc.cores)
		hic_dists = unlist(pbmclapply(1:nrow(comppairs),function(x){
				.hic_kl_maps(sims[[comppairs$i[x]]],sims[[comppairs$j[x]]],theta=2,pix.size=pix.size)
					   },mc.cores=mc.cores))
		return(hic_dists)
	}
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
## foldback = FALSE: when only one loose end has residual CN left, leave it loose
## (walks end there) instead of closing it with a foldback onto itself.
## region: if given (GRanges), only loose ends on nodes overlapping it are pasted, to each other;
## all other loose ends stay loose. NULL pastes genome-wide.
paste_loose_ends <- function(gg, seed = NULL, foldback = TRUE, region = NULL) {
    if (is.null(region)) {
        lo <- gg$loose
    } else {
        ids <- unique((gg$nodes$gr %&% region)$node.id)
        if (!length(ids))
            return(gg$copy)
        lo <- gg$nodes[ids]$loose
    }
    if (!length(lo))
        return(gg$copy)
    work <- gr2dt(lo)[cn > 0 & terminal == FALSE][, row := .I][, .(strand, node.cn, cn, index, node.id, orientation, row)]
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
            if (!foldback) {
                residual[i] <- 0L
                next
            }
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

## EagleC SV calls (the *_combined.txt table: chr1 pos1 chr2 pos2 then the probability of each orientation
## ++ +- -+ --), one row per call and orientation with probability >= min_p, highest first.
## EagleC's orientation is the opposite of gGnome's strand (as in sclarke's read_eagle, srctools R/EagleC.R):
## EagleC '+' = gGnome '-' = right side of a node (at its end), EagleC '-' = gGnome '+' = left side (at its start).
## gGnome junctions and loose ends use the same strand <-> side convention.
read_eaglec = function(file, min_p = 0.1){
	d = fread(file)
	setnames(d, c('chr1','pos1','chr2','pos2','++','+-','-+','--'))
	d[, call_id := .I]
	d[, (5:8) := lapply(.SD, as.numeric), .SDcols = 5:8]
	d = data.table::melt(d, id.vars = c('call_id','chr1','pos1','chr2','pos2'), variable.name = 'orient', value.name = 'p')[p >= min_p]
	flip = c('+' = '-', '-' = '+')
	d[, orient := as.character(orient)]
	d[, `:=`(chr1 = gsub('chr','',chr1), chr2 = gsub('chr','',chr2),
		 str1 = unname(flip[substr(orient,1,1)]), str2 = unname(flip[substr(orient,2,2)]))]
	d[order(-p)]
}

## index of the nearest end (ends: seqnames, pos, strand) within tol of each breakpoint, same strand; NA if none
.eaglec_near = function(chr, bp, str, ends, tol){
	vapply(seq_along(chr), function(i){
		d = abs(ends$pos - bp[i])
		ok = which(ends$seqnames == chr[i] & ends$strand == str[i] & d <= tol)
		if (length(ok)) ok[which.min(d[ok])] else NA_integer_}, 0L)
}

## Join loose ends with Hi-C junction calls (read_eaglec()) instead of the heuristic in paste_loose_ends().
## A call is used if both its breakpoints lie within tol of a loose end with the matching strand and at least
## one of the two loose ends is on a node overlapping region (NULL: anywhere). Calls are used in order of
## probability, each with CN = min of the two ends' residual loose CN (both ends on one loose end: a foldback,
## CN = residual %/% 2). Every other loose end stays loose (a walk endpoint). Calls whose partner is not a
## loose end are not used (that would need a new breakpoint and CN at the partner, i.e. a graph refit).
## Returns the new graph (always freshly built, so it never shares an edge table with gg), or with
## return_all = TRUE list(gg, calls = calls touching region loose ends with the partner type, pasted, loose).
paste_loose_ends_eaglec = function(gg, calls, region = NULL, tol = 1e4, return_all = FALSE){
	lo = gr2dt(gg$loose)
	lo = if (nrow(lo)) lo[cn > 0 & terminal == FALSE] else lo
	lo = lo[, .(seqnames = as.character(seqnames), pos = start, strand = as.character(strand), node.id, orientation, cn)]
	lo[, inreg := if (is.null(region)) TRUE else node.id %in% (gg$nodes$gr %&% region)$node.id]
	calls = copy(calls)[, `:=`(lo1 = .eaglec_near(chr1, pos1, str1, lo, tol), lo2 = .eaglec_near(chr2, pos2, str2, lo, tol))]
	calls[, `:=`(inreg1 = lo$inreg[lo1] %in% TRUE, inreg2 = lo$inreg[lo2] %in% TRUE)]
	hit = calls[inreg1 | inreg2]
	hit[, partner := ifelse(is.na(lo1) | is.na(lo2), 'no loose end',
		ifelse(lo1 == lo2, 'same loose end (foldback)', ifelse(inreg1 & inreg2, 'loose end in region', 'loose end outside')))]
	resid = as.integer(lo$cn)
	moves = list()
	for (r in which(!is.na(hit$lo1) & !is.na(hit$lo2))){
		i = hit$lo1[r]; j = hit$lo2[r]
		w = if (i == j) resid[i] %/% 2L else min(resid[i], resid[j])
		if (w < 1) next
		resid[i] = resid[i] - w; resid[j] = resid[j] - w
		moves[[length(moves) + 1]] = data.table(call_id = hit$call_id[r], cn = w, n1 = lo$node.id[i], n2 = lo$node.id[j],
							n1.side = lo$orientation[i], n2.side = lo$orientation[j])
	}
	pasted = rbindlist(moves)
	edges = copy(gg$edges$dt)
	if (nrow(pasted)){
		newedges = pasted[, .(cn = sum(cn), type = 'ALT'), by = .(n1, n2, n1.side, n2.side)] #several calls can join the same two ends
		edges = rbind(edges, newedges, fill = TRUE)
	}
	newgg = loosefix(gG(nodes = gg$nodes$gr, edges = edges))
	newgg$set(y.field = 'cn')
	if (!return_all) return(newgg)
	lo[, cn_left := resid]
	list(gg = newgg, calls = hit, pasted = pasted, loose = lo)
}


