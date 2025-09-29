prep_data_for_py <- function(gg,gm,depth,pix.size=1e5,target_region=NULL,write.to.folder=NULL){ 
  if (is.null(target_region)){
    target_region = gg$footprint
  }
  #pull data from graph
  nodesdt = gg$nodes$dt[,.(start,end,seqnames,snode.id,cn)]
  nodesgr = gg$nodes$gr
  edgesdt = gg$edges$dt[,.(cn,n1,n2,n1.side,n2.side,type)]

  #make sure edges of nodes and target region match
  targetnodes = gr2dt(gr.merge(target_region,nodesgr))

  #split nodes in half: our model of walks on the genome graph requires that we optimize over all possible wirings of left-junctions to right-junctions at a given node. To this end we consider the left- and right- halves of a node as distinct.
  nodesdt.split = targetnodes[rep(1:.N,each=2)][,side:=ifelse(mod(.I,2)==1,'left','right')][,.(start,end,seqnames,node.id,cn,side)][,width:=ifelse(side=='left',floor((end-start+1)/2),ceil((end-start+1)/2))]
  nodesdt.split[side=='left',end:=start+width-1]
  nodesdt.split[side=='right',start:=end-width+1]

  #now we consider each copy of each half-node as unique, and give it a unique ID. This is so we can keep track of junction mappings
  dedup.nodes = nodesdt.split[rep(1:.N,cn)][,.(start,end,seqnames,node_id = node.id,side)][,unique.id:=.I]

  #Now we have to tile the target region which we are going to simulate and fit to data, making sure our tiles match up with edges of nodes
  tiled.target = gr.merge(gr.tile(target_region,pix.size),dt2gr(nodesdt.split))[,c()]

  #we want to consider every "copy" of every tile uniquely
  dedup.tiles = gr2dt(gr.merge(tiled.target,dt2gr(dedup.nodes)))[,tile.id:=query.id]

  #When we simulate later, we will be calculating distances between node-halves because that is the finest degree of information given by the genome graph, but then we will extrapolate from this the distances between pixels. So here we calculate the offset between a given pixel and the left/right edges of the half-node in which it is.
  dedup.tiles[,leftdist:=start-dedup.nodes[unique.id]$start + width/2]
  dedup.tiles[,rightdist:=dedup.nodes[unique.id]$end - end + width/2]
  dedup.tiles = unique(dedup.tiles[,.(width,unique.id,tile.id,leftdist,rightdist)])

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

  # rebin Hi-C data to match simulation tiles
  gm_rebin = gm$disjoin(tiled.target)$agg(tiled.target)
  
  # put in a named list and write to files for Python to read
  data.out = list(tiled_target = gr2dt(tiled.target),dedup_tiles = dedup.tiles, nodes = dedup.nodes, internal_edges = internal_edges, external_edges = external_edges,gm_dat = gm_rebin$dat,depth=depth)
  if (!is.null(write.to.folder)){
    lapply(1:length(data.out),function(i){
      thisname = names(data.out)[i]
      write.csv(data.out[[i]],paste0(write.to.folder,thisname,'.csv'))
    })
  }
  return(data.out)
}
