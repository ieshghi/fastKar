# this file contains the tools necessary for performing inference on graphs

#goes from a ggraph to a "wiring", which gives all the internal edges of the graph (going from left side of a node to right side)
#along with the loose node ids and a reference data table with the new node ids (all copies of each node are de-duplicated)

gg.to.wiring <- function(gg){ 
  nodesdt = gg$nodes$dt[,.(start,end,seqnames,snode.id,cn,loose.cn.left,loose.cn.right)]
  nodesgr = gg$nodes$gr[,c('snode.id','cn')]
  edgesdt = gg$edges$dt[,.(cn,n1,n2,n1.side,n2.side,type)]
  #make loose edges from nodesdt, then combine them will all other edges, then split all edges as left- and right-edges coming off of nodes
  nodesdt[is.na(loose.cn.right),loose.cn.right:=cn]
  nodesdt[is.na(loose.cn.left),loose.cn.left:=cn] #adjusting CN at chromosome ends
  # separate left and right loose ends and label them appropriately as edges going to node "0", before adding them to the total set of external edges
  left.looseedges = nodesdt[loose.cn.left>0][,.(n2=snode.id,cn = loose.cn.left,n2.side='left',n1=0,n1.side='right')]
  right.looseedges = nodesdt[loose.cn.right>0][,.(n1=snode.id,cn = loose.cn.right,n2.side='left',n2=0,n1.side='right')]
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
traverse_graph <- function(A, loose.ends) {
  visited_edges <- c()
  visited_rows <- c()
  traversed_paths <- list()
  # Only process loose ends that haven't been traversed
  remaining_loose_ends <- setdiff(loose.ends, visited_edges)
  while (length(remaining_loose_ends) > 0) {
    start_edge <- remaining_loose_ends[1]  # Pick an unvisited loose end
    path <- c()
    nodepath <- c()
    current_edge <- start_edge
    while (!is.na(current_edge) && !(current_edge %in% remaining_loose_ends && length(path) > 0)) {
      path <- c(path, current_edge)
      visited_edges <- c(visited_edges, current_edge)
      
      # Find the row where current_edge appears in 'left' or 'right'
      row <- A[(left == current_edge | right == current_edge) & !(id %in% visited_rows)]
      visited_rows <- c(visited_rows,row$id)
      # Determine the next edge
      if (row$left == current_edge) {
        nodepath <- c(nodepath,row$n)
        next_edge <- row$right
      } else {
        nodepath <- c(nodepath,-row$n)
        next_edge <- row$left
      }
      # Find the next occurrence of next_edge (excluding current row)
      next_row <- A[(left == next_edge | right == next_edge) & !(id %in% visited_rows)]
      if (nrow(next_row) == 0){ 
        visited_edges = c(visited_edges,next_edge) 
        path = c(path,next_edge) 
        break  # No more paths
      } 
      current_edge <- next_edge  # Move to the next edge
    }
    traversed_paths = append(traversed_paths,list(nodepath))
    remaining_loose_ends <- setdiff(loose.ends, visited_edges)  # Update unvisited loose ends
  }
  #now we detect cycles
  remaining_edges <- A[!(id %in% visited_rows)]
  cycles <- list()
  while (nrow(remaining_edges) > 0) {
    start_edge <- remaining_edges$right[1]  # Pick an unvisited edge
    visited_rows = c(visited_rows,remaining_edges$id[1])
    path <- c()
    nodepath = c(remaining_edges$n[1])
    current_edge <- start_edge
    while (!is.na(current_edge)) {
      path <- c(path, current_edge)
      # Find the row containing this edge
      row <- remaining_edges[(left == current_edge | right == current_edge) & !(id %in% visited_rows)]
      if (nrow(row) == 0) {
        cycles <- append(cycles, list(nodepath))
        break
      }
      visited_rows <- c(visited_rows,row$id)
      # Determine the next edge
      if (row$left == current_edge) {
        nodepath <- c(nodepath,row$n)
        next_edge <- row$right
      } else {
        nodepath <- c(nodepath,-row$n)
        next_edge <- row$left
      }
      # If we revisit an edge, a cycle is found
      current_edge <- next_edge  # Move to next edge
    }
    # Remove visited edges from remaining_edges
    remaining_edges <- remaining_edges[!(id %in% visited_rows)]
  }
  return(list(paths=traversed_paths,cycles=cycles))
} # ---> to be made into a Rcpp function for faster evaluation

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
  
  // Construct and return the result list
  List result;
  result["paths"] = traversed_paths;
  result["cycles"] = cycles;
  return result;
}
')

walks.from.edges <- function(wiring,shuffle=0,return.gw = FALSE,ifcpp=TRUE){ #wrapper function
  internal.edges = wiring$internal.edges
  loose.ends = wiring$loose.ends
  gg = wiring$gg
  edges = copy(internal.edges) 
  if (shuffle==1){
    edges[,right:=ifelse(cn>1,sample(right,unique(cn)),right),by=n] #shuffle rewiring
  } else if (shuffle==2){
    nswap = sample(edges[cn>1]$n,1)
    edges[n==nswap,sright:=sample(right,unique(cn))]
    edges[,right:=ifelse(is.na(sright),right,sright)][,sright:=NULL]
  } else if (shuffle < 0){
    nswap = -shuffle
    edges[n==nswap,sright:=sample(right,unique(cn))]
    edges[,right:=ifelse(is.na(sright),right,sright)][,sright:=NULL]
  }
  if (ifcpp){
    walks.out = traverse_graph_cpp(internal.edges,loose.ends)
  } else{
    walks.out = traverse_graph(internal.edges,loose.ends)
  }
  if (return.gw){
    return(c(gW(graph=gg,snode.id=walks.out$paths),gW(graph=gg,snode.id=walks.out$cycles,circular=TRUE)))
  }else{
    circular = c(rep(F,length(walks.out$paths)),rep(T,length(walks.out$cycles)))
    snode.id = c(walks.out$paths,walks.out$cycles)
    return(list(graph=gg,snode.id=snode.id,circular=circular))
  }
  return(walks.out)
}