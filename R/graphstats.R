#TODO for all inference methods: make some hashing method so that user can look at all the sampled walk sets and their respective NLLs

#goes from a ggraph to a "wiring", which gives all the internal edges of the graph (going from left side of a node to right side)
#along with the loose node ids and a reference data table with the new node ids (all copies of each node are de-duplicated)
gg.to.wiring <- function(gg){ 
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
  left.looseedges = copy(nodesdt)[loose.cn.left>0][,.(n2=snode.id,cn = loose.cn.left,n2.side='left',n1=0,n1.side='right')]
  right.looseedges = copy(nodesdt)[loose.cn.right>0][,.(n1=snode.id,cn = loose.cn.right,n2.side='left',n2=0,n1.side='right')]
  external.edges = rbind(edgesdt,rbind(left.looseedges,right.looseedges)[,type:='LOO'][,.(cn,n1,n2,n1.side,n2.side,type)])
  #browser()
  dedup.edges = external.edges[rep(1:.N,cn)][,.(n1,n2,n1.side,n2.side,type)] #separate all copies of all edges
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

# given a wiring, with internal edges and loose ends, traverse graph and return paths and cycles
cppFunction('List traverse_graph_cpp(DataFrame A, NumericVector loose_ends) {
  // Extract columns from the data frame
  IntegerVector n = A["n"];
  IntegerVector left = A["left"];
  IntegerVector right = A["right"];
  IntegerVector id = A["id"];
  
  int nrow = n.size();
  
  // Data structures for tracking
  std::set<int> visited_edges;
  std::set<int> visited_rows;
  List traversed_paths;
  
  // Make a copy of loose_ends to modify
  std::set<int> remaining_loose_ends;
  for (int i = 0; i < loose_ends.size(); i++) {
    remaining_loose_ends.insert(loose_ends[i]);
  }
  
  // Remove loose ends that have already been visited
  for (auto it = remaining_loose_ends.begin(); it != remaining_loose_ends.end();) {
    if (visited_edges.find(*it) != visited_edges.end()) {
      it = remaining_loose_ends.erase(it);
    } else {
      ++it;
    }
  }
  
  // First part: traverse paths starting from loose ends
  while (!remaining_loose_ends.empty()) {
    int start_edge = *remaining_loose_ends.begin();  // Pick an unvisited loose end
    std::vector<int> path;
    std::vector<int> nodepath;
    int current_edge = start_edge;

    
    while (true) {
      path.push_back(current_edge);
      visited_edges.insert(current_edge);
      
      // Find the row where current_edge appears in "left" or "right"
      int row_idx = -1;
      for (int i = 0; i < nrow; i++) {
        if ((left[i] == current_edge || right[i] == current_edge) && 
            visited_rows.find(id[i]) == visited_rows.end()) {
          row_idx = i;
          break;
        }
      }
      
      if (row_idx == -1) break;  // No valid row found
      
      visited_rows.insert(id[row_idx]);
      
      // Determine the next edge
      int next_edge;
      if (left[row_idx] == current_edge) {
        nodepath.push_back(n[row_idx]);
        next_edge = right[row_idx];
      } else {
        nodepath.push_back(-n[row_idx]);
        next_edge = left[row_idx];
      }
      
      // Check if next_edge is a loose end and weve already traversed something
      if (!path.empty() && 
          remaining_loose_ends.find(next_edge) != remaining_loose_ends.end()) {
        break;
      }
      
      // Find the next row containing next_edge
      int next_row_idx = -1;
      for (int i = 0; i < nrow; i++) {
        if ((left[i] == next_edge || right[i] == next_edge) && 
            visited_rows.find(id[i]) == visited_rows.end()) {
          next_row_idx = i;
          break;
        }
      }
      
      if (next_row_idx == -1) {
        // No more paths, add the final edge
        visited_edges.insert(next_edge);
        path.push_back(next_edge);
        break;
      }
      
      current_edge = next_edge;  // Move to the next edge
    }
    
    // Add the completed path
    if(nodepath.size()>0){
        traversed_paths.push_back(nodepath);
    }
    
    // Update unvisited loose ends
    remaining_loose_ends.clear();
    for (int i = 0; i < loose_ends.size(); i++) {
      if (visited_edges.find(loose_ends[i]) == visited_edges.end()) {
        remaining_loose_ends.insert(loose_ends[i]);
      }
    }
  }
  
  // Second part: detect cycles
  List cycles;
  
  // Find remaining unvisited rows
  std::vector<int> remaining_row_indices;
  for (int i = 0; i < nrow; i++) {
    if (visited_rows.find(id[i]) == visited_rows.end()) {
      remaining_row_indices.push_back(i);
    }
  }
  
  while (!remaining_row_indices.empty()) {
    int start_idx = remaining_row_indices[0];
    int start_edge = right[start_idx];  // Start from right edge
    
    visited_rows.insert(id[start_idx]);
    
    std::vector<int> path;
    std::vector<int> nodepath;
    nodepath.push_back(n[start_idx]);  // Start with the node value
    
    int current_edge = start_edge;
    
    while (true) {
      path.push_back(current_edge);
      
      // Find the row containing this edge
      int row_idx = -1;
      for (int i : remaining_row_indices) {
        if ((left[i] == current_edge || right[i] == current_edge) && 
            visited_rows.find(id[i]) == visited_rows.end()) {
          row_idx = i;
          break;
        }
      }
      
      if (row_idx == -1) {
        // No more valid rows, cycle is complete
        cycles.push_back(nodepath);
        break;
      }
      
      visited_rows.insert(id[row_idx]);
      
      // Determine the next edge
      int next_edge;
      if (left[row_idx] == current_edge) {
        nodepath.push_back(n[row_idx]);
        next_edge = right[row_idx];
      } else {
        nodepath.push_back(-n[row_idx]);
        next_edge = left[row_idx];
      }
      
      current_edge = next_edge;  // Move to next edge
    }
    
    // Update remaining row indices
    remaining_row_indices.clear();
    for (int i = 0; i < nrow; i++) {
      if (visited_rows.find(id[i]) == visited_rows.end()) {
        remaining_row_indices.push_back(i);
      }
    }
  }

  int n_paths = traversed_paths.size();
  int n_cycles = cycles.size();
  int n_total = n_paths + n_cycles;
  List combined(n_total);
  LogicalVector circular(n_total);

  for (int i = 0; i < n_paths; i++) {
      combined[i] = traversed_paths[i];
      circular[i] = false;
  }

  for (int i = 0; i < n_cycles; i++) {
      combined[n_paths + i] = cycles[i];
      circular[n_paths + i] = true;
  }
  
  return List::create(
      _["snode.id"] = combined,
      _["circular"] = circular
  );
}
')

#sample gWalks from graph gg, take N samples and return all unique permutations among them

#' @import pbmcapply
#' @import digest
sample.gwalks = function(gg,N=1,mc.cores=1,chunksize = 1e3,return.gw=T,remove.dups=T){ 
  wiring = gg.to.wiring(gg)
  internal.edges = wiring$internal.edges
  loose.ends = wiring$loose.ends
  message('Sampling permutations')

  shuffle_edges <- function(edges) {
    new_right <- edges[, if (.N > 1) sample(right, .N) else right, by = n]$V1
  }
  shuffle_chunk <- function(K, edges) {
    replicate(K, shuffle_edges(edges), simplify = FALSE)
  }
  chunks <- N %/% chunksize
  if (chunks > 1){
  	perms <- do.call('c',pbmclapply(seq_len(chunks), function(i) shuffle_chunk(chunksize, internal.edges),mc.cores = mc.cores))
  }else{
  	perms <- pbmclapply(seq_len(N), function(i) shuffle_edges(internal.edges),mc.cores = mc.cores)
  }
  message('Hashing and keeping unique permutations')
  hashes = unlist(lapply(perms,digest))

  message('Generating walks')
  dt = data.table(hash=hashes,idx=1:N)[,id:=1:.N,by=hash]
  uniqueperms = perms[dt[id==1]$idx]
  permchunks = split(uniqueperms, ceiling(seq_along(uniqueperms)/chunksize))
  makewalk_chunk <- function(permchunk) {
    ws = lapply(permchunk,function(p){traverse_graph_cpp(internal.edges[,right:=p],loose.ends)})
    if (return.gw){
      ws = lapply(ws,function(w){gW(graph=gg,snode.id=w$snode.id,circular=w$circular)})
    }
    return(ws)
  }
  walks.out <- do.call('c',pbmclapply(permchunks, makewalk_chunk,mc.cores = mc.cores))
  if (remove.dups){
  	if (return.gw){
  	        hashes = do.call('c',lapply(walks.out,function(w){w$hash}))}
  	else{
  	        hashes = do.call('c',lapply(walks.out,function(w){hash_snodelist(w$snode.id,w$circular)}))}
  	hash.dt = data.table(hash=hashes)[,idx:=.I][,id:=1:.N,by=hash]
  	walks.out = walks.out[hash.dt[id==1]$idx]
  } 
  return(walks.out)
}


hash_snodelist = function(snode.id,circular){
  circ = ifelse(rep(circular,2),'C','L')
  nodepcomp = c(snode.id,lapply(snode.id,function(s){-rev(s)}))
  nodestring = lapply(1:length(nodepcomp),function(i){
    	return(paste0(toString(nodepcomp[[i]]),circ[i]))
  })
  return(toString(sort(do.call('c',nodestring))))
}

sort_snodes = function(nodelist,arr=NULL) {
  choose_compl = lapply(nodelist, function(x) {
    rc = -rev(x)
    if (paste(x, collapse = ",") <= paste(rc, collapse = ",")) {
      x
    } else {
      rc
    }
  })
  ord <- order(sapply(choose_compl, paste, collapse = ","))
  sorted_nodes = choose_compl[ord]
  if (!is.null(arr)){
    sorted_arr = arr[ord]
    return(list(nodelist=sorted_nodes,arr=sorted_arr))
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

sample_and_collapse = function(ggraph,samplesize,min.wid=3e7,ifsmooth=F,bands=NULL,region=NULL,mc.cores=1){
	if (is.null(region)){
		region = si2gr(hg_seqlengths(chr=FALSE)) %Q% (seqnames %in% c(1:22,'X','Y'))
	}
	if (ifsmooth){
		message('Smoothing dels and dups')
		finalgraph = smoothdeldups(ggraph %&% region)
	}else{ finalgraph = ggraph$copy
	}
	message('Sampling walks')
	if (is.null(finalgraph)){
		return(list(card.estimate=1,coverage=1,collapsed.walks=NULL,bands=bands))
	}else if (nrow(finalgraph$edges$dt) == 0 ){
		return(list(card.estimate=1,coverage=1,collapsed.walks=NULL,bands=bands))
	}else{
		walkset = sample.gwalks(finalgraph,samplesize,return.gw=F,mc.cores=mc.cores)
		message('Collapsing to bands')
		if (is.null(bands)){
			#bands = gr.tile(finalgraph$footprint,min.wid)
			datafolder = '~/projects/karyotyping/data'
			bands = readRDS(paste0(datafolder,'/coarsebands.rds'))
		} else if (inherits(bands,'GRanges')){
		} else {
			error('Bands provided must either be GRanges or a scale to tile the genome')
		}
		collapsed.walks = collapse_gwalklist(walkset,bands.gr = bands,graph = finalgraph,min.wid=min.wid,mc.cores=mc.cores)
		message('Hashing and counting samples')
		coll_hashes = sapply(collapsed.walks, `[[`, "hash")
	}
	coll_hash.dt = data.table(collapsed.hash=unlist(coll_hashes))[,walkset.id:=.I][,collapsed.id:=as.integer(factor(coll_hashes))]
	coll_hash.dt[,count:=.N,by=collapsed.id]
	setkeyv(coll_hash.dt,'collapsed.id')
	coll_hash.dt[,instance:=1:.N,by=collapsed.id]
	coll_idx = coll_hash.dt[instance==1]$walkset.id
	unique.collapsed = coll_hash.dt[instance==1,.(walkset.id,collapsed.id)]
	collapsed.walks = collapsed.walks[coll_idx]
	return(list(graph=finalgraph,walkset=walkset,bands=NULL,collapsed.walks=collapsed.walks,unique.collapsed=unique.collapsed,bands=bands,hash.dt=coll_hash.dt))
}


#' @import pbapply
collapse_gwalklist <- function(gwlist,bands.gr,graph,min.wid,mc.cores=1,chunksize=1e3){
  gr = graph$nodes$gr[,c('snode.id')]
  grdt = gr2dt(gr[width(gr)>=min.wid])
  keepnodes = grdt$snode.id
#  bands.gr$band.width = width(bands.gr)
#  merged.dt = gr2dt(gr.merge(gr,bands.gr)[,c('query.id','subject.id','band.width')])[,.(node.id=query.id,band.id=subject.id,width=width,seqnames,start,end,band.width)]
  gw.chunks <- split(gwlist, ceiling(1:length(gwlist) / chunksize))
  collapsed.list = pblapply(gw.chunks, function(chunk) {
    collapsed.chunk = mclapply(chunk, function(gw){
      walknodes = gw$snode.id
      circular = gw$circular
#      sorted = sort_snodes(walknodes,arr=circular)
#      walknodes_sorted = sorted$nodelist
#      circular = sorted$arr
#      collapsed.walk = lapply(walknodes_sorted,function(x){
      collapsed.walk = lapply(walknodes,function(x){
        walk = merge.data.table(data.table(step=1:length(x),snode.id=x,node.id = abs(x)),grdt[,.(node.id=snode.id,width)],by='node.id',allow.cartesian=T)[,.(step,snode.id,width)]
	setkeyv(walk,'step')
	if (nrow(walk)){
		return(walk$snode.id)
	}else{
		return(NULL)
	}
#	walk = x[abs(x) %in% keepnodes]
#	if (total_walkwidth > min.wid){ #only keep walks longer than min.wid
#        	walk[strand<0,step.id:=rev(step.id),by=order] #need to add concensus strand functionality
#		setkeyv(walk,c('order','step.id'))
#		wdt = walk[,.(band.id,width,band.width,strand)][,cluster:=rleid(band.id)]
#		wdt[,bandfrac:=sum(width)/band.width,by=cluster]
#		wdt[,strandmean:=sum(strand*width)/band.width,by=cluster]
#        	if (nrow(walk)){
#		  bands.dt = unique(wdt[,.(band.id,bandfrac,strandmean)])
#        	  bandsvec = bands.dt$band.id
#		  bandscopy = round(bands.dt$bandfrac)
#		  strandvec = ifelse(bands.dt$strandmean>=0,1,-1)
#        	  bandsvec = rep(bandsvec,bandscopy)
#        	  strandvec = rep(strandvec,bandscopy)
#		  if (length(strandvec*bandsvec)>0){
#        	  return(strandvec*bandsvec)}
#		  else{return(NULL)}
#        	}else{
#        	  return(NULL)
#        	}
#	} else{return(NULL)}
      })
      keep = unlist(lapply(collapsed.walk,function(w){!is.null(w)}))
      collapsed.walk = collapsed.walk[keep]
      circular = circular[keep]
#      sorted = sort_snodes(collapsed.walk,arr=circular)
      return(list(snode.id = collapsed.walk,circular = circular,hash=hash_snodelist(collapsed.walk,circular)))
      }, mc.cores=mc.cores)
      return(collapsed.chunk)
      })
      return(do.call('c',collapsed.list))
}

to_gwalk = function(walklist,gr,mc.cores=1){
	grl = mclapply(walklist$snode.id,function(nl){
		this.gr = gr[abs(nl)]
		strand(this.gr) = ifelse(nl>0,'+','-')
		return(this.gr)
      },mc.cores=mc.cores)
	grl = do.call('GRangesList',grl)
	return(gW(grl=grl,circular=walklist$circular)$disjoin())
}

reads_fromwalk = function(walk,readL,min.res = NULL){
	gr = walk$graph$gr[,c('node.id')]
	reads = do.call('rbind',lapply(1:length(walk$snode.id),function(i){
		snodes = walk$snode.id[[i]]
		widths = width(gr[abs(snodes)])
		if (!is.null(min.res)){
			keep = widths >= min.res
			snodes = snodes[keep]
			widths = widths[keep]
			if (sum(keep)==0){
				return(NULL)
			}
		}
		grdt = data.table(snode.id=snodes)
		grdt[,end:=cumsum(widths)]
		grdt[,start:=end-widths+1]
		grdt = grdt[,.(start,end,snode.id)]
		maxlen = max(grdt$end)
		lin_maxlen = copy(maxlen)
		if (walk$circular[i]){ #check if this is a circular walk
			add.gr = data.table::copy(grdt)
			add.gr = add.gr[start < readL][,start:=start+maxlen][,end:=min(end+maxlen+1,maxlen+readL)] #add at the end however much is equal to the hanging piece
			grdt = rbind(grdt,add.gr) #add at the end however much is equal to the hanging piece
			maxlen = max(grdt$end)
		}
		breakpts = c(grdt$start)
		reads = rbind(data.table(start = breakpts)[,end:=start+readL],data.table(end=breakpts)[,start:=end-readL])[,.(start,end)]
		reads[end>maxlen,end:=maxlen]
		reads[start<=0,start:=1]
		reads = unique(reads[end>start & start <= lin_maxlen])
		setkeyv(reads,c('start','end'))
		reads[,num:=diff(c(start,lin_maxlen+1))][,id:=.I]
		setkeyv(reads,c('start','end'))
		setkeyv(grdt,c('start','end'))
		ovdt = foverlaps(x=reads[,.(start,end,num,id)],y=grdt[,.(start,end,snode.id)],by.x=c('start','end'),by.y=c('start','end'),type='any',nomatch=0L)[,copy:=1:.N,by=id]
		nodeid_list = ovdt[, .(snode.id = list(snode.id)), by = id][order(id), snode.id]
		nums = ovdt[copy==1]$num
		nodeid_list = c(nodeid_list,lapply(nodeid_list,function(n){-rev(n)}))
		nodeid_list = sapply(nodeid_list,paste,collapse='')
		return(data.table(words=nodeid_list,nums=rep(nums,2)))
	}))
	reads[,nums:=sum(nums),by=words]
	return(unique(reads))
}

longread_kl <- function(walk_x, walk_y, graph=NULL,readL=1e4, depth = NULL, min.res=NULL, background = 1e-5,mc.cores=1) {
	if (is.null(walk_x$graph)){
		if(is.null(graph)){
			error('Must provide either a gWalk object or a graph as input to function')
		}
		walk_x$graph = graph
		walk_y$graph = graph
	}else{
		graph = walk_x$graph
	}
	x_list = reads_fromwalk(walk_x,readL,min.res=min.res)
	y_list = reads_fromwalk(walk_y,readL,min.res=min.res)
	all_words = union(x_list$words,y_list$words)
	c1 = setNames(x_list$nums, x_list$words)[all_words]
	c2 = setNames(y_list$nums, y_list$words)[all_words]
	c1[is.na(c1)] <- 0
	c2[is.na(c2)] <- 0
	p = (c1 + background) / (sum(c1) + background * length(all_words))
	q = (c2 + background) / (sum(c2) + background * length(all_words))
	#m = 0.5 * (p + q)
	klval = sum(p*log(p/q))
  	#KL = function(x, y) sum(x * log(x / y))
  	#JS = 0.5 * KL(p, m) + 0.5 * KL(q, m)
	if (is.null(depth)){
		return(klval)
	} else{
		domain = walk_x$graph$footprint
		covperread = min(readL/sum(width(domain)),1)
		N = depth/covperread
		return(N*klval)
	}
}

hic_kl <- function(walk_x, walk_y, graph=NULL,pix.size=1e6, depth = 1,theta=2) {
	if (is.null(walk_x$graph)){
		if(is.null(graph)){
			error('Must provide either a gWalk object or a graph as input to function')
		}
		walk_x$graph = graph
		walk_y$graph = graph
	}else{
		graph = walk_x$graph
	}
	hic_x = forward_simulate(walk_x,target_region = graph$footprint,pix.size=pix.size,depth=depth)
	hic_y = forward_simulate(walk_y,target_region = graph$footprint,pix.size=pix.size,depth=depth)
	kl = kl_nb(hic_x$value,hic_y$value,r=theta)
	return(kl)
}
