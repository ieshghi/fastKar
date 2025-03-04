# this file contains some basic tests for the fastKar package

test_fastKar_forward <- function(){
    somewalks = makerepdup(1e5,'cis1')
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
