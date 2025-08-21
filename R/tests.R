gimme_bfb <- function(w,n,cycles=5){
  tiles = gr.tile(parse.gr('1'),w)[1:n]
  return(makebfb(tiles,cycles=cycles,background=F)[1])
}

makebfb  <- function(segments,cycles=5,background=TRUE){
  segs = gG(nodes=unique(gr.stripstrand(segments)))$gr[,c()]
  segs$node.id = 1:(length(segs))
  binids = (segs %Q% (strand=='+'))$node.id
  #assume these are bins at the end of a chomosome
  for(this.cycle in 1:cycles){
    nbin = length(binids)
    binlabels = 1:nbin
    breakpt = sample(binlabels[2:(length(binlabels)-1)],1)
    if (breakpt <= nbin/2){
      keepbins = binlabels[binlabels>=breakpt]
      right_seg = binids[keepbins]
      left_seg = -rev(binids[keepbins])
    }else{
      keepbins = binlabels[binlabels<=breakpt]
      left_seg = binids[keepbins]
      right_seg = -rev(binids[keepbins])
    }
    binids = c(left_seg,right_seg)
  }
  walk.gr = segs[abs(binids)]
  strand(walk.gr) = ifelse(sign(binids)>0,'+','-')
  norm.gr = segs %Q% (strand=='+')
  if(background){
    nodewalk = gW(grl=GRangesList(walk.gr,norm.gr))
  }else{
    nodewalk = gW(grl=GRangesList(walk.gr))
  }
  nodewalk$set(cn=1)
  nodewalk$nodes$mark(cn=1)
  nodewalk$edges$mark(cn=1)
  nodewalk$disjoin()
  return(nodewalk)
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

makewalk <- function(tiles,nodewalk,circ=FALSE,background = FALSE,stranded=TRUE){
  segments = gG(nodes=tiles)$gr[,c()]
  segments$node.id = 1:length(segments)
  if (is.list(nodewalk)){
    if (length(circ)==1){
      circ = rep(circ,length(nodewalk))
    } else if(length(circ)!=length(nodewalk)){
      stop('Declare circular array of length 1 or equal to number of walks')
    }
    walk.gr = lapply(nodewalk,function(w){
      s = rep('+',length(w))
      wa = abs(w)
      s[w<0]='-'
      out = segments[abs(w)]
      strand(out) = s
      return(out)
      })
  }else{
      s = rep('+',length(nodewalk))
      wa = abs(nodewalk)
      s[nodewalk<0]='-'
      out = segments[abs(nodewalk)]
      strand(out) = s
      walk.gr = out
  }
  if (background){
    circ = c(circ,rep(FALSE,background))
    norm.gr = segments %Q% (strand=='+')
    walk = gW(grl = GRangesList(walk.gr,norm.gr),circular=circ)
    walk$set(cn=c(1,as.numeric(background)))
  }else{
    walk = gW(grl = GRangesList(walk.gr),circular=circ)
    walk$set(cn=1)
  }
  walk$nodes$mark(cn=1)
  walk$edges$mark(cn=1)
  walk$disjoin()
  return(walk)
}
