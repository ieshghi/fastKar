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
  nodesdt.split = nodesdt[cn>0][rep(1:.N,each=2)][,lr:=ifelse(mod(.I,2)==1,'l','r')][,.(start,end,seqnames,snode.id,cn,lr)][,width:=ifelse(lr=='l',floor((end-start+1)/2),ceil((end-start+1)/2))]
  nodesdt.split[lr=='l',end:=start+width-1]
  nodesdt.split[lr=='r',start:=end-width+1]
  ref.splitnodes = nodesdt.split[,.(start,end,nodeid=paste0(abs(snode.id),lr),width,refid=.I)]
  internal.edges[,id:=.I]
  return(list(internal.edges = internal.edges,loose.ends = loose.subids,gr.splitnodes=ref.splitnodes,gg=gg))
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
