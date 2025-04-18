f1score_comparemaps <- function(null_map,hyp_map,theta=3,significance = 0.05,return_samples=FALSE){
  num_samples = 10/significance
  samples = lapply(list(null_map,hyp_map),function(m){make_noisydat(m,num_samples,theta)}) #samples to classify}
  likelihooddiff = lapply(samples,function(s){
                   unlist(lapply(s,function(ss){
                    return(compdats(ss,null_map$dat,theta)-compdats(ss,hyp_map$dat,theta))
                    }))
                   })
  likcutoff = unname(quantile(likelihooddiff[[1]],1-significance)) #bootstrap likelihood ratio cutoff by sampling from the null
  if (return_samples){
    return(list(likelihooddiff[[1]],likelihooddiff[[2]],likcutoff))
  }else{
  reject_null = likelihooddiff[[2]] > likcutoff

  tps = sum(reject_null==TRUE) #since hyp_samples are non-null, rejecting the null in those samples is a true positive
  fns = sum(reject_null==FALSE) #for the same reason as above, keeping the null in those samples is a false negative
  fps = significance*num_samples #by design, false positive rate is the significance level we choose times the number of samples
  f1 = 2*tps/(2*tps + fps + fns)
  return(f1)
  }
}

test.walks.with.hic <- function(walkset,hic.data,resolution=1e5,mc.cores=1,target_region=NULL,return='scores'){
    if(!is.list(walkset)){
        stop('Give me multiple walks with the same footprint in a list to compare!')
    }
    if(is.null(target_region)) {
        firstwalk = walkset[[1]]
        target_region=firstwalk$footprint
    }
    depth.est = estimate.depthratio(hic.data,resolution)
    predictions = mclapply(walkset,function(w){
        return(forward_simulate(w,target_region=target_region,if.comps=F,pix.size=resolution,mc.cores=1,verbose=F,if.sum=T,depth=depth.est,model=0))
    },mc.cores=mc.cores)
    rebin.data = (hic.data$disjoin(predictions[[1]]$gr))$agg(predictions[[1]]$gr) #make sure data is aggregated on the same GRanges as the predictions

    if (return=='scores'){
        scores = mclapply(predictions,function(pred){compmaps(rebin.data,pred,ifsum=TRUE,theta=3)},mc.cores=mc.cores)
        return(scores)
    } else if (return=='scoremaps'){
        scoremaps = mclapply(predictions,function(pred){compmaps(rebin.data,pred,ifsum=FALSE,theta=3)},mc.cores=mc.cores)
        return(scoremaps)
    } else if (return=='predictions'){
        return(predictions)
    } else if (return == 'sp') {
        scores = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=TRUE,theta=3)},mc.cores=mc.cores)
        return(list(scores = scores, predictions = predictions))
    } else if (return == 'all') {
        scores = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=TRUE,theta=3)},mc.cores=mc.cores)
        scoremaps = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=FALSE,theta=3)},mc.cores=mc.cores)
        return(list(scores = scores, predictions = predictions, scoremaps = scoremaps, hic.data = rebin.data))
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

estimate.depthratio <- function(hic.data,res=1e5){ #put here hic data in a non-rearranged place. Estimates the depth ratio between this data and the reference data
  lookup.data = fastKar::small_lookup
  dat = hic.data$dat[i!=j]
  dat[,d:=res*(j-i)]
  dat[,source:='skov3']
  dat[,tot_density:=mean(value/res^2,na.rm=TRUE),by=d]
  dat = unique(dat[d>=res & d<=1e7,.(d,tot_density,source)])
  est = dat[d==1e5]$tot_density/mean(lookup.data[d==1e5]$tot_density) #gets ratio between counts in reference data (1500x) and provided data
  return(1500*est)
}

compdats <- function(dat_test,dat_true,theta=0,ifscale=FALSE,ifsum=TRUE,checkinds=TRUE,if.diag=TRUE){
    if(checkinds){
        setkeyv(dat_test,c('i','j'))
        setkeyv(dat_true,c('i','j'))
        dat_test[,id:=.I]
        dat_true[,id:=.I]
    }
    if (ifscale){
        dat_test_scale = (dat_test$value)*mean(dat_true$value,na.rm=TRUE)/mean(dat_test$value,na.rm=TRUE) #re-scale using means
        dat_test$value = dat_test_scale
    }
    combdat = merge.data.table(dat_test[,.(i,j,id,value)],dat_true[,.(id,value)],by='id')
    if(if.diag==FALSE){
        combdat = combdat[i!=j]
    }
    if (theta>0){
        combdat[,logprob:=dnbinom(round(value.x),mu = value.y,size = theta,log=TRUE)]
    }else{
        combdat[,logprob:=dpois(round(value.x),lambda = value.y,log=TRUE)]
    }
    combdat[value.y==0,logprob:=0]
    if (ifsum==TRUE){
        return(-sum(combdat$logprob,na.rm=TRUE))
    }else{
        return(combdat[,.(value=-logprob,i,j,id)])
    }
}

make_noisymap <- function(map_in,theta=0,num.copies=1,mc.cores=20){
    dats = make_noisydat(map_in,num.copies,theta)
    out.maps = mclapply(dats,function(i){gM(gr=map_in$gr,dat=i)},mc.cores=mc.cores)
    return(out.maps)
}

make_noisydat <- function(map_in,num.copies=1,theta=0,backlambda = 0){ #samples from negative binomial distribution, unless theta=0 then uses poisson
    mapdat = map_in$dat
    meanvals = mapdat$value
    if (theta==0){
        newval = rpois(length(meanvals)*num.copies,rep(meanvals,num.copies))
    }else {
        newval = rnbinom(length(meanvals)*num.copies,size=theta,mu=rep(meanvals,num.copies))
    }
    if (backlambda > 0){
        newval = newval + rpois(length(newval),backlambda)
    }
    multimap = do.call(rbind,rep(list(mapdat),num.copies))
    multimap$map.ids = lapply(1:num.copies,function(i){rep(i,nrow(mapdat))}) %>% unlist
    multimap$value = newval
    out.dats = split(multimap,by='map.ids')
    #out.gms = lapply(out.dats,function(i){gM(gr=map_in$gr,dat=i)})
    return(out.dats)
}