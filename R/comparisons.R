test.walks.with.hic <- function(walkset,hic.data,resolution=1e5,mc.cores=1,return='scores'){
    if(!is.list(walkset)){
        stop('Give me multiple walks with the same footprint in a list to compare!')
    }
    firstwalk = walkset[[1]]
    target=firstwalk$footprint
    depth.est = estimate.depthratio(hic.data,resolution)
    predictions = mclapply(walkset,function(w){
        return(run_analysis(w,target_region=target,if.comps=F,pix.size=resolution,mc.cores=1,verbose=F,if.sum=T,depth=depth.est,model=0))
    },mc.cores=mc.cores)
    rebin.data = (hic.data$disjoin(predictions[[1]]$gr))$agg(predictions[[1]]$gr) #make sure data is aggregated on the same GRanges as the predictions

    if (return=='scores'){
        scores = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=TRUE,theta=3)},mc.cores=mc.cores)
        return(scores)
    } else if (return=='scoremaps'){
        scoremaps = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=FALSE,theta=3)},mc.cores=mc.cores)
        return(scoremaps)
    } else if (return=='predictions'){
        return(predictions)
    }
}

compmaps <- function(map_test,map_true,theta=0,ifsum=FALSE,ifscale = FALSE,if.diag=TRUE){ 
    #returns the negative log likelihood of test map given true map and negative-binomial noise
    template = make_template_dat(map_true$gr)
    dat_test = merge.data.table(template[,.(i,j,id)],map_test$dat[,.(i,j,value)],all.x=TRUE,by=c('i','j'))
    dat_true = merge.data.table(template[,.(i,j,id)],map_true$dat[,.(i,j,value)],all.x=TRUE,by=c('i','j'))
    dat_test[is.na(value),value:=0]
    dat_true[is.na(value),value:=0]
    comp_dat = compdats(dat_test,dat_true,theta,ifscale,ifsum=FALSE,if.diag=if.diag)
    if (ifsum){
        return(sum(comp_dat$value))
    }else{
        return(gM(gr=map_true$gr,dat=comp_dat))
    }
}

ensemble_compmaps <- function(map_test,map_true,n=100,theta=0,ifscale=FALSE,if.diag=TRUE,ifsum=FALSE,mc.cores=1){
    gr_test = map_test$gr
    dat_true = map_true$dat
    test.ensemble = make_noisydat(map_test,num.copies=n,theta=theta)
    comps = mclapply(test.ensemble,function(m){compdats(m,dat_true,theta,ifscale=ifscale,ifsum=ifsum,if.diag=if.diag)},mc.cores=mc.cores)
    if(ifsum){
        return(mean(unlist(comps)))
    }else{
        return(gM(gr=gr_test,dat=sum_matrices(comps))/n)
    }
}

estimate.depthratio <- function(hic.data,res,ifplot=FALSE){ #put here hic data in a non-rearranged place. Estimates the depth ratio between this data and the reference data
  lookup.data = fastKar::small_lookup
  dat = hic.data$dat[i!=j]
  dat[,d:=res*(j-i)]
  dat[,source:='skov3']
  dat[,tot_density:=mean(value/res^2,na.rm=TRUE),by=d]
  dat = unique(dat[d>=res & d<=1e7,.(d,tot_density,source)])
  est = dat[d==1e5]$tot_density/mean(lookup.data[d==1e5]$tot_density)
  return(est)
}