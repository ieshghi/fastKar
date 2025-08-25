#TODO for all inference methods: make some hashing method so that user can look at all the sampled walk sets and their respective NLLs

#goes from a ggraph to a "wiring", which gives all the internal edges of the graph (going from left side of a node to right side)
#along with the loose node ids and a reference data table with the new node ids (all copies of each node are de-duplicated)
gg.to.wiring <- function(gg){ 
  if (!(('loose.cn.left' %in% colnames(gg$dt)) & ('loose.cn.right' %in% colnames(gg$dt)))){
    gg = balance(gg$copy)
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

sample.gwalks = function(gg,N=1,mc.cores=1,chunksize = 1e3,return.gw=T){ 
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
  return(walks.out)
}


hash_snodelist = function(snode.id,circular){
  circ = ifelse(rep(circular,2),'C','L')
  nodepcomp = c(snode.id,lapply(snode.id,function(s){-rev(s)}))
  nodestring = lapply(1:length(nodepcomp),function(i){
    paste0(toString(nodepcomp[[i]]),circ[i])
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

collapse_gwalklist <- function(gwlist,bands.gr,graph,min.wid=1e3,mc.cores=1,chunksize=1e3){
  gr = graph$nodes$gr[,c('snode.id')]
  merged.dt = gr2dt(gr.merge(gr,bands.gr)[,c('query.id','subject.id')])[,.(node.id=query.id,band.id=subject.id,width=width,seqnames,start,end)]

  gw.chunks <- split(gwlist, ceiling(1:length(gwlist) / chunksize))
  collapsed.list = pblapply(gw.chunks, function(chunk) {
    collapsed.chunk = mclapply(chunk, function(gw){
      walknodes = gw$snode.id
      circular = gw$circular
      sorted = sort_snodes(walknodes,arr=circular)
      walknodes_sorted = sorted$nodelist
      circular = sorted$arr
      collapsed.walk = lapply(walknodes_sorted,function(x){
        walk = merge.data.table(data.table(order = 1:length(x),node.id = abs(x), strand = sign(x)),merged.dt,by='node.id')[,step.id:=.I]
        walk = walk[width>min.wid]
        # we should do something more sophisticated to simulate blurriness!!!
        if (nrow(walk)){
          walk[strand<0,step.id:=rev(step.id),by=order]
          setkeyv(walk,c('order','step.id'))
          bandsvec = walk$band.id
          nodup_inds = c(TRUE, diff(bandsvec) != 0)
          bandsvec_nodup = bandsvec[nodup_inds]
          strandvec_nodup = walk$strand[nodup_inds]
          return(strandvec_nodup*bandsvec_nodup)
        }else{
          return(NULL)
        }
      })
      keep = do.call('c',lapply(collapsed.walk,function(x){!is.null(x)}))
      sorted = sort_snodes(collapsed.walk[keep],arr=circular[keep])
      return(list(band.id = sorted$nodelist,circular = sorted$arr,hash=hash_snodelist(sorted$nodelist,sorted$arr)))
      }, mc.cores=mc.cores)
      return(collapsed.chunk)
      })
      return(do.call('c',collapsed.list))
}

#should introduce a threshold width for removing del/dups
smoothdeldups = function(ggraph){
	dup_junctions = ggraph$edgesdt[!is.na(ggraph$edgesdt$dup)][n1==n2]
	del_junctions = ggraph$edgesdt[!is.na(ggraph$edgesdt$del)][abs(n1-n2)==2]
	del_junctions[,this.n:=round((n1+n2)/2)]
	nodesgr = ggraph$nodes$gr
	cnvec = nodesgr$cn
	cnvec[dup_junctions$n1] = cnvec[dup_junctions$n1]-dup_junctions$cn
	cnvec[del_junctions$this.n] = cnvec[del_junctions$this.n]+del_junctions$cn
	nodesgr$cn = cnvec
	nagraph = ggraph$copy
	nagraph$nodes$mark(cn=NA)
	nagraph$edges$mark(cn=NA)
	nodeldups = balance(nagraph,marginal=nodesgr,verbose=0)[,cn>0]
	return(loosefix(nodeldups)$simplify())
}

cardinality_estimate = function(ggraph,samplesize,min.wid=5e6,region=NULL,datafolder=NULL,bandsfile = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg38.cytoband.txt",mc.cores=1){
	if (is.null(region)){
		region = si2gr(hg_seqlengths(chr=FALSE)) %Q% (seqnames %in% c(1:22,'X','Y'))
	}
	finalgraph = smoothdeldups(ggraph %&% region)

	walkset = sample.gwalks(finalgraph,samplesize,return.gw=F,mc.cores=mc.cores)
	hashes = pbmclapply(walkset,function(gw){hash_snodelist(gw$snode.id,gw$circular)},mc.cores=mc.cores)
	hash.dt = data.table(hashes=unlist(hashes))
	hash.dt[,count:=.N,by=hashes]
	estimate_card = chao1(unique(hash.dt)$count)
	
	#do these walks collapse
	idx = hash.dt[,.I[1],by=hashes]$V1
	bands.td = gTrack::karyogram(file = bandsfile)
	bands = bands.td@data
	bands = gr.nochr(grl.unlist(do.call(`GRangesList`, bands)))
	bands = bands %Q% (seqnames %in% c(1:22,'X','Y'))
	bands = dt2gr(gr2dt(bands)[,start:=start+1],seqlengths=seqlengths(bands))
	
	collapsed.walks = collapse_gwalklist(walkset[idx],bands,finalgraph,min.wid=min.wid,mc.cores=mc.cores,chunksize=1e3)
	
	coll_hashes = sapply(collapsed.walks, `[[`, "hash")
	coll_hash.dt = data.table(hashes=unlist(coll_hashes))
	coll_idx = coll_hash.dt[,.I[1],by=hashes]$V1
	coll_hash.dt[,count:=.N,by=hashes]
	coll_hash.dt[,countcount:=.N,by=count]
	coll_hash.dt = unique(coll_hash.dt)
	counttable = unique(coll_hash.dt[,.(count,countcount)])
	coverage = 1 - (counttable[count==1]$countcount)/(counttable$countcount %>% sum)
	setkeyv(coll_hash.dt,'count')
	estimate_coll_card = chao1(coll_hash.dt$count)

	return(list(card.estimate = estimate_coll_card,coverage=coverage,chao_reliability = reliability,collapsed.walks=collapsed.walks,bands=bands))
}
