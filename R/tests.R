# this file contains some basic tests for the fastKar package

test_fastKar_forward <- function(){
    tiles = gr.tile(parse.gr('1:1:1e6'),3e5)
    nodeseq = c(1,2,-4,-3)
    walk = makewalk(tiles,nodeseq)
    somewalks$dt[,cn:=1]
    gmats = run_analysis(somewalks,pix.size=1e4,mc.cores=3)
    return(gmats)
}

# simple reciprocal dup
makerepdup <- function(regsize,geometry){
  target_chroms = c(1,2)
  locus1 = GRanges(target_chroms[1],ranges=IRanges(1,regsize*3))
  locus2 = GRanges(target_chroms[2],ranges=IRanges(1,regsize*3))
  gr = c(gr.tile(locus1,regsize),gr.tile(locus2,regsize))[,c()]
  gr$label = c('A','B','C','D','E','F')
  gr$cn=1
  names(gr) = gr$label
  if (geometry=='cis1'){
    gr1 = gr[c('A','B','C')]
    gr2 = gr[c('D','E','F')]
    gr3 = gr[c('A','B','E','B','C')]
    walk = gW(grl=GRangesList(gr1,gr2,gr2,gr3))
    walk$graph$edges$mark(cn=1)
    walk$graph$nodes$mark(cn=1)
    walk$disjoin()
  }else if(geometry=='cis2'){
    gr1 = gr[c('A','B','C')]
    gr2 = gr[c('D','E','F')]
    gr3 = gr[c('D','E','B','E','F')]
    walk = gW(grl=GRangesList(gr1,gr1,gr2,gr3))
    walk$graph$edges$mark(cn=1)
    walk$graph$nodes$mark(cn=1)
    walk$disjoin()
  }else if(geometry=='trans'){
    gr1 = gr[c('A','B','C')]
    gr2 = gr[c('D','E','F')]
    gr3 = gr[c('D','E','B','C')]
    gr4 = gr[c('A','B','E','F')]
    walk = gW(grl=GRangesList(gr1,gr2,gr3,gr4))
    walk$graph$edges$mark(cn=1)
    walk$graph$nodes$mark(cn=1)
    walk$disjoin()
  }
  return(walk)
}

makewalk <- function(tiles,nodewalk,background = FALSE){
  segments = gG(nodes=tiles)$gr[,c()]
  segments$node.id = 1:length(segments)
  if (is.list(nodewalk)){
    walk.gr = lapply(nodewalk,function(w){segments[w]})
  }else{
    walk.gr = segments[nodewalk]
  }
  if (background){
    norm.gr = segments %Q% (strand=='+')
    walk = gW(grl = GRangesList(walk.gr,norm.gr))
    walk$set(cn=c(1,as.numeric(background)))
  }else{
    walk = gW(grl = GRangesList(walk.gr))
    walk$set(cn=1)
  }
  walk$nodes$mark(cn=1)
  walk$edges$mark(cn=1)
  walk$disjoin()
  return(walk)
}