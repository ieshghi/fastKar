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
  loose.subids = split.edgetable[type=='LOO']$subid
  #split the nodes in the graph to half-nodes. This needs to be done only one time. Then, aggregate Hi-C data to those half-nodes.
  nodesdt.split = nodesdt[cn>0][rep(1:.N,each=2)][,lr:=ifelse(mod(.I,2)==1,'l','r')][,.(start,end,seqnames,snode.id,cn,lr)][,width:=ifelse(lr=='l',floor((end-start+1)/2),ceil((end-start+1)/2))]
  nodesdt.split[lr=='l',end:=start+width-1]
  nodesdt.split[lr=='r',start:=end-width+1]
  ref.splitnodes = nodesdt.split[,.(nodeid=paste0(abs(snode.id),lr),width,refid=.I)]
  return(list(internal.edges,loose.subids,ref.splitnodes,gg))
}

# given a wiring (a mapping of left-edges to right-edges), there should be a unique set of walks which traverse that graph.
paths_from_edges_R <- function(iedge, jedge, nodes, num_edges, looseends) {
  # Number of edges and rows
  n <- num_edges
  nrows <- length(nodes)
  # Create adjacency list representation
  adj_list <- vector("list", n)
  # Fill adjacency list from edge list
  for (k in seq_len(nrows)) {
    u <- iedge[k]  # 1-based indexing
    v <- jedge[k]
    adj_list[[u]] <- append(adj_list[[u]], list(c(v, k)))
    adj_list[[v]] <- append(adj_list[[v]], list(c(u, -k)))  # Sign to track strand
  }
  # Result lists to store paths and cycles
  paths <- list()
  cycles <- list()
  visited_edge <- rep(FALSE, n)  # Track visited edges
  # Function to traverse paths
  traverse_path <- function(start) {
    path <- c()
    current <- start
    visited_edge[current] <- TRUE
    prev <- -1
    while (TRUE) {
      found_next <- FALSE
      for (neighbor in adj_list[[current]]) {
        next_edge <- neighbor[1]
        row_id <- neighbor[2]

        if (!visited_edge[next_edge]) {
          visited_edge[next_edge] <- TRUE
          path <- c(path, if (row_id > 0) nodes[row_id] else -nodes[-row_id])
          prev <- current
          current <- next_edge
          found_next <- TRUE
          break
        }
      }
      if (!found_next) break  # End of path
    }
    return(if (length(path) > 0) path else NULL)
  }
  # Traverse paths from looseends
  for (start in looseends) {
    path <- traverse_path(start)
    if (!is.null(path)) {
      paths <- append(paths, list(path))
    }
  }
  # Detect cycles for remaining unvisited edges
  for (k in seq_len(nrows)) {
    if (!visited_edge[iedge[k]]) {
      cycle <- c()
      current <- iedge[k]
      visited_edge[current] <- TRUE
      prev <- -1

      while (TRUE) {
        found_next <- FALSE

        for (neighbor in adj_list[[current]]) {
          next_edge <- neighbor[1]
          row_id <- neighbor[2]

          if (!visited_edge[next_edge]) {
            visited_edge[next_edge] <- TRUE
            cycle <- c(cycle, if (row_id > 0) nodes[row_id] else -nodes[-row_id])
            prev <- current
            current <- next_edge
            found_next <- TRUE
            break
          }
        }
        if (!found_next || (length(cycle) > 1 && current == iedge[k])) {
          break  # Cycle closed
        }
      }
      if (length(cycle) > 0) {
        cycles <- append(cycles, list(cycle))
      }
    }
  }
  return(list(paths = paths, cycles = cycles))
}
