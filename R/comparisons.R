f1score_comparemaps <- function(null_map,hyp_map,theta=2,significance = 0.05,return_samples=FALSE,mask=NULL){
  if (!is.null(mask)){
      gr = null_map$gr %&% mask
      bad.inds = gr$tile.id
      nulldat = null_map$dat
      hypdat = hyp_map$dat
      nulldat[i %in% bad.inds,value:=0]
      nulldat[j %in% bad.inds,value:=0]
      hypdat[i %in% bad.inds,value:=0]
      hypdat[j %in% bad.inds,value:=0]
      null_map = gM(gr = null_map$gr,dat=nulldat)
      hyp_map = gM(gr = hyp_map$gr,dat=dat)
  }
  num_samples = 4/significance
  samples = lapply(list(null_map,hyp_map),function(m){make_noisydat(m,num_samples,theta)}) #samples to classify}
  likelihooddiff = lapply(samples,function(s){
                   unlist(lapply(s,function(ss){
                    return(compdats(ss,null_map$dat,theta)-compdats(ss,hyp_map$dat,theta))
                    }))
                   })
  likcutoff = unname(quantile(likelihooddiff[[1]],1-significance)) #bootstrap likelihood ratio cutoff by sampling from the null
  reject_null = likelihooddiff[[2]] > likcutoff
  tps = sum(reject_null==TRUE) #since hyp_samples are non-null, rejecting the null in those samples is a true positive
  fns = sum(reject_null==FALSE) #for the same reason as above, keeping the null in those samples is a false negative
  fps = significance*num_samples #by design, false positive rate is the significance level we choose times the number of samples
  f1 = 2*tps/(2*tps + fps + fns)
  if (return_samples){
    return(list(likelihooddiff[[1]],likelihooddiff[[2]],likcutoff))
  } else{
  return(f1)
  }
}

kl_nb <- function(mu1, mu2, r=2,tinyval=1e-10,ifsum=T) { #get KL divergence between two sets of negative binomials
  mu1 <- as.numeric(pmax(mu1,tinyval))
  mu2 <- as.numeric(pmax(mu2,tinyval))
  stopifnot(length(mu1) == length(mu2))
  kl <- mu1*log(mu1/mu2) + (r+mu1)*log((r+mu2)/(r+mu1))
  kl[kl<0] = 0
  if(ifsum){
  	return(sum(kl))
  }else{
	return(kl)
  }
}

klsymm_nb <- function(mu1, mu2, r=2,tinyval=1e-10) { #get KL divergence between two sets of negative binomials
  mu1 <- as.numeric(pmax(mu1,tinyval))
  mu2 <- as.numeric(pmax(mu2,tinyval))
  stopifnot(length(mu1) == length(mu2))
  
  kl <- (mu1-mu2)*log((mu1/mu2)*((r+mu2)/(r+mu1)))
  kl[kl<0] = 0
  return(sum(kl))
}

test.walks.with.hic <- function(walkset,hic.data,resolution=1e5,mc.cores=1,depth.est=1,target_region=NULL,if.diag=TRUE,return='scores',mask=NULL){
    if(!is.list(walkset)){
        stop('Give me multiple walks with the same footprint in a list to compare!')
    }
    if(is.null(target_region)) {
        firstwalk = walkset[[1]]
        target_region=firstwalk$footprint
    }
    predictions = mclapply(walkset,function(w){
        return(forward_simulate(w,target_region=target_region,if.comps=F,pix.size=resolution,mc.cores=1,if.sum=T,depth=depth.est,model=0))
    },mc.cores=mc.cores)
    rebin.data = (hic.data$disjoin(predictions[[1]]$gr))$agg(predictions[[1]]$gr) #make sure data is aggregated on the same GRanges as the predictions

    if (return=='scores'){
        scores = mclapply(predictions,function(pred){compmaps(rebin.data,pred,ifsum=TRUE,theta=2,if.diag=if.diag,mask=mask)},mc.cores=mc.cores)
        return(scores)
    } else if (return=='scoremaps'){
        scoremaps = mclapply(predictions,function(pred){compmaps(rebin.data,pred,ifsum=FALSE,theta=2,if.diag=if.diag,mask=mask)},mc.cores=mc.cores)
        return(scoremaps)
    } else if (return=='predictions'){
        return(predictions)
    } else if (return == 'sp') {
        scores = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=TRUE,theta=2,if.diag=if.diag,mask=mask)},mc.cores=mc.cores)
        return(list(scores = scores, predictions = predictions, hic.data = rebin.data))
    } else if (return == 'all') {
        scores = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=TRUE,theta=2,if.diag=if.diag,mask=mask)},mc.cores=mc.cores)
        scoremaps = mclapply(predictions,function(pred){compmaps(pred,rebin.data,ifsum=FALSE,theta=2,if.diag=if.diag,mask=mask)},mc.cores=mc.cores)
        return(list(scores = scores, predictions = predictions, scoremaps = scoremaps, hic.data = rebin.data))
    }
}

compmaps <- function(map_test,map_true,theta=0,ifsum=FALSE,ifscale = FALSE,if.diag=TRUE,mask=NULL,return_kl=F,area0=NULL){ 
    if (!is.null(mask)){
        gr = map_true$gr %&% mask
        bad.inds = gr$tile.id
        trudat = map_true$dat
        trudat[i %in% bad.inds,value:=0]
        trudat[j %in% bad.inds,value:=0]
        map_true = gM(gr=map_true$gr,dat=trudat)
    }
   if (is.null(area0)){
	gr = map_true$gr
        medianwid = median(width(gr))
	area0 = medianwid^2
   }
    template = make_template_dat(map_true$gr)
    #returns the negative log likelihood of test map given true map and negative-binomial noise
    dat_test = merge.data.table(template[,.(i,j,widthprod,id)],map_test$dat[,.(i,j,value)],all.x=TRUE,by=c('i','j'))
    dat_true = merge.data.table(template[,.(i,j,widthprod,id)],map_true$dat[,.(i,j,value)],all.x=TRUE,by=c('i','j'))
    dat_test[is.na(value),value:=0]
    dat_true[is.na(value),value:=0]
    comp_dat = compdats(dat_test,dat_true,theta,ifscale,ifsum=FALSE,if.diag=if.diag,return_kl=return_kl,area0=area0)
    if (ifsum){
        return(sum(comp_dat$value))
    }else{
        return(gM(gr=map_true$gr,dat=comp_dat))
    }
}

estimate.depthratio <- function(filepath,mode='hic',res=1e6,ploidy=2,if.chr=FALSE){ #put here hic data for a whole genome
    wholegenome = gr.tile(si2gr(hg_seqlengths()[1:24]),res)
    if (if.chr){
        wholegenome = gr.chr(wholegenome)
    }
    if (mode=='hic'){
        depthest.hic = straw(filepath,res=as.integer(res),gr=wholegenome)
    }else if (mode=='mcool'){
        depthest.hic = cooler(filepath,res=as.integer(res),gr=wholegenome)
    }else if (mode=='rds'){
        depthest.hic = readRDS(filepath)
    }
    return(depthest.hic$value%>%sum * 300 / 3e9 * 2/ploidy)
}

compdats <- function(dat_test,dat_true,theta=0,ifscale=FALSE,ifsum=TRUE,checkinds=TRUE,if.diag=TRUE,return_kl=F,area0=1e8){
    if (is.null(dat_test$widthprod)){
        dat_test$widthprod= area0 #if the data.table doesn't have an area column, assume all pixels are the same size
    }
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
    combdat = merge.data.table(dat_test[,.(i,j,widthprod,id,value)],dat_true[,.(id,value)],by='id')
    if(if.diag==FALSE){
        combdat = combdat[i!=j]
    }
    if (theta>0){
	if (return_kl){
		combdat[,kldiv:=kl_nb(value.x,value.y,r=theta,ifsum=F)]
	}else{
        	combdat[,logprob:=dnbinom(round(value.x),mu = value.y,size = theta,log=TRUE)]
	}
    }else{
        combdat[,logprob:=dpois(round(value.x),lambda = value.y,log=TRUE)]
    }
    combdat[value.y==0,logprob:=0]
    if (ifsum==TRUE){
        return(-sum((combdat$widthprod/area0)*combdat$logprob,na.rm=TRUE))
    }else{
	if (return_kl){
        	return(combdat[,.(value=widthprod/area0*kldiv,i,j,id)])
	}else{
        	return(combdat[,.(value=-logprob*widthprod/area0,i,j,id)])
	}
    }
}

make_noisymap <- function(map_in,theta=0,num.copies=1,mc.cores=20){
    dats = make_noisydat(map_in,num.copies,theta)
    out.maps = mclapply(dats,function(i){gM(gr=map_in$gr,dat=i)},mc.cores=mc.cores)
    return(out.maps)
}

make_noisydat <- function(map_in,num.copies=1,theta=0,backlambda = 0){ #samples from negative binomial distribution, unless theta=0 then uses poisson
    mapdat = map_in$dat
    template = make_template_dat(map_in$gr)
    mapdat = merge.data.table(template[,.(i,j,widthprod,id)],mapdat[,.(i,j,value)],all.x=TRUE,by=c('i','j'))
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

hic_cov = function(gm){ #get coverage from Hi-C map, at the granges in that map
	gr = gm$gr
	gr$tile.id=1:length(gr)
	temp = make_template_dat(gr)
	grdt = gr2dt(gr)
	dat = merge.data.table(gm$dat,temp[,.(i,j)],by=c('i','j'),all.y=T)
	dat[is.na(value),value:=0]
	out = rbind(
	  dat[, .(i, value)],
	  dat[i != j, .(i = j, value)]
	)[, .(sum.value = sum(value)), by = i]
	grdt = merge.data.table(grdt,out[,.(tile.id=i,cov=sum.value)],by='tile.id')
	return(dt2gr(grdt,seqlengths = seqlengths(gm$gr)))
}

callloops = function(hic,gw,depth,fdr_cut = 0.01,resolution = NULL,sim.dat = NULL,reg = NULL,qqplot_path = NULL,mode='peaks'){ #call Hi-C loops 
	if (is.null(reg)){
		reg = streduce(hic$gr)
	}else{
		hic.gr2 = hic$gr %&% reg
		hic = hic$disjoin(hic.gr2)$agg(hic.gr2)
	}
	if (is.null(resolution)){
		resolution = width(hic$gr[1])
	}
	if (is.null(sim.dat)){
		sim.dat = forward_simulate(gw,target_region=reg,pix.size=resolution,depth=depth)
	}
	sim.dat = sim.dat$disjoin(hic$gr)$agg(hic$gr)
	compdat = merge.data.table(sim.dat$dat[,.(i,j,expected=value)],hic$dat[,.(i,j,measured=value)],by=c('i','j'),all.x=T)
	compdat[is.na(measured),measured:=0]
	if(mode=='peaks'){
		compdat[,nbinom_p:=pnbinom(measured-1,size=2,mu=expected,lower.tail=F)]
	}else if(mode=='all'){
		compdat[,nbinom_p:=dnbinom(measured,size=2,mu=expected)]
	}
	compdat[,fdr:=p.adjust(nbinom_p,'BH')]
	compdat[,signif:=fdr<fdr_cut]
	compdat_out = compdat[signif==T,.(i,j,measured,expected,logoe=log(measured/expected),neglogp=-log(nbinom_p))]
	compdat_out[measured==0 & expected==0,logoe:=0]
	outlier.gm = list(gr=sim.dat$gr,dat=compdat_out[,.(i,j,logoe,neglogp)])
	simulated = data.table(value=rnbinom(nrow(compdat),mu=compdat$expected,size=2),expected=compdat$expected)
	simulated[,nbinom_p:=pnbinom(value-1,size=2,mu=expected,lower.tail=F)]
	if (!is.null(qqplot_path)){
		ppdf(qq_pval(compdat$nbinom_p,exp=simulated$nbinom_p,max.x=10,max.y=10),qqplot_path,width=5,height=4) 
	}
	return(outlier.gm)
}

