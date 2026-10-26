## ---------------------------------------------------------------------------------------------
## Internal helpers for the fast paths below. Every fast path reproduces the reference
## (merge/setkey based) implementation exactly; the reference code is kept as the fallback for
## inputs the fast path does not cover (unaligned or duplicated pixels, unusual column types).
## ---------------------------------------------------------------------------------------------

## The pixel table of a gMatrix without the defensive copy that gMatrix$dat makes on every
## access (8.5 M rows at 4131 tiles). Read-only: callers must never modify the result.
.gm_rawdat = function(gm){
	if (inherits(gm, 'R6')){
		p = tryCatch(gm$.__enclos_env__$private, error = function(e) NULL)
		if (!is.null(p) && is.data.table(p$pdat)) return(p$pdat)
	}
	gm$dat
}

## tile widths as make_template_dat() sees them, or NULL if it would not accept the bins
.tile_widths = function(gr){
	if (inherits(gr, 'GRanges')) return(as.numeric(width(gr)))
	if (inherits(gr, 'data.table') && !is.null(gr[['width']]) && is.numeric(gr[['width']])) return(as.numeric(gr[['width']]))
	NULL
}

## TRUE iff d's i, j, id columns are identical() to make_template_dat()'s for n tiles
.is_template = function(d, n){
	is.list(d) && !is.null(d[['i']]) && !is.null(d[['j']]) && !is.null(d[['id']]) &&
		cmp_template_match(d[['i']], d[['j']], d[['id']], as.integer(n))
}

## TRUE iff compdats(a, b) would merge a and b row by row after setkey(., 'id'): identical
## i/j/id columns with strictly increasing (so already id-sorted, unique, non-NA) ids
.aligned_pair = function(a, b){
	cols = c('i', 'j', 'id', 'value')
	all(cols %in% names(a)) && all(cols %in% names(b)) &&
		identical(a[['i']], b[['i']]) && identical(a[['j']], b[['j']]) && identical(a[['id']], b[['id']]) &&
		isFALSE(is.unsorted(a[['id']], strictly = TRUE))
}

## Per-pixel compdats() values on aligned vectors, exactly as the reference computes them on the
## merged table: x = observed (value.x), y = expected (value.y), wp = widthprod (vector or scalar).
## Returns the value vector over the rows kept by if.diag (keep = NULL means all rows).
.compdats_values = function(x, y, wp, theta, return_kl, area0){
	if (theta>0){
		if (return_kl){
			kldiv = kl_nb(x, y, r=theta, ifsum=F)
		}else{
			logprob = dnbinom(round(x), mu = y, size = theta, log=TRUE)
		}
	}else{
		logprob = dpois(round(x), lambda = y, log=TRUE)
	}
	if (return_kl){
		return(wp/area0*kldiv)
	}
	logprob[which(y==0)] = 0
	-logprob*wp/area0
}

## Fast compdats() core for an aligned pair (see .aligned_pair). Returns the per-pixel values
## and the row subset (keep = NULL: all rows) that compdats(ifsum = FALSE) would return.
.compdats_aligned = function(dat_test, dat_true, theta, if.diag, return_kl, area0){
	x = dat_test[['value']]
	y = dat_true[['value']]
	wp = if (is.null(dat_test$widthprod)) area0 else dat_test[['widthprod']]
	keep = NULL
	if (if.diag==FALSE){
		keep = which(dat_test[['i']] != dat_test[['j']])
		x = x[keep]; y = y[keep]
		if (length(wp) > 1) wp = wp[keep]
	}
	list(value = .compdats_values(x, y, wp, theta, return_kl, area0), keep = keep)
}

## can the fast core stand in for the reference merge path?
.compdats_fast_ok = function(theta, return_kl, area0){
	is.numeric(area0) && length(area0) == 1 && is.numeric(theta) && length(theta) == 1 && !is.na(theta) &&
		!(return_kl && !(theta > 0)) ## reference errors (no kldiv column) for KL with theta <= 0
}

## the ifsum = FALSE table of compdats(): value, i, j, id over the kept rows, keyed by id
.compdats_table = function(dat_test, res){
	sel = function(v) if (is.null(res$keep)) v else v[res$keep]
	out = data.table(value = res$value, i = sel(dat_test[['i']]), j = sel(dat_test[['j']]), id = sel(dat_test[['id']]))
	setattr(out, 'sorted', 'id')
	out
}

## Template layout and mean values that make_noisydat() draws from, without the merge.
## NULL (and no RNG use) when the inputs need the reference path.
.noisy_prep = function(map_in, nsamp){
	if (!is.numeric(nsamp) || length(nsamp) != 1 || is.na(nsamp) || nsamp < 1 || nsamp != round(nsamp)) return(NULL)
	w = .tile_widths(map_in$gr)
	if (is.null(w) || !length(w)) return(NULL)
	n = length(w)
	d = .gm_rawdat(map_in)
	if (!is.list(d) || !all(c('i', 'j', 'value') %in% names(d))) return(NULL)
	v = d[['value']]
	if (!(is.double(v) || is.integer(v) || is.logical(v))) return(NULL)
	tc = cmp_template_cols(w, TRUE)
	P = length(tc$i)
	if (length(v) == P && identical(d[['i']], tc$i) && identical(d[['j']], tc$j)){
		meanvals = as.double(v)
	}else{
		if (!is.integer(d[['i']]) || !is.integer(d[['j']])) return(NULL)
		if (anyDuplicated(as.data.table(list(i = d[['i']], j = d[['j']])))) return(NULL) ## merge would add rows
		idx = cmp_pix_index(d[['i']], d[['j']], n)
		ok = which(!is.na(idx))
		meanvals = rep(NA_real_, P)
		meanvals[idx[ok]] = as.double(v[ok])
	}
	meanvals[is.na(meanvals)] = 0
	list(i = tc$i, j = tc$j, widthprod = tc$widthprod, id = seq_len(P) + 0L, meanvals = meanvals, P = P) # id materialized (not ALTREP)
}

## the single vectorized draw of make_noisydat() (same RNG stream)
.noisy_draw = function(prep, nsamp, theta){
	meanvals = prep$meanvals
	if (theta==0){
		rpois(length(meanvals)*nsamp,rep(meanvals,nsamp))
	}else {
		rnbinom(length(meanvals)*nsamp,size=theta,mu=rep(meanvals,nsamp))
	}
}

## dense template-ordered values of a pixel table (missing and NA -> 0), or NULL if the table
## has non-integer indices, duplicated pixels or non-double values
.dense_values = function(d, n){
	if (!is.list(d) || !all(c('i', 'j', 'value') %in% names(d)) || !is.double(d[['value']])) return(NULL)
	P = n * (n + 1) / 2
	if (.is_template(d, n)){
		v = d[['value']]
	}else{
		if (!is.integer(d[['i']]) || !is.integer(d[['j']])) return(NULL)
		if (anyDuplicated(as.data.table(list(i = d[['i']], j = d[['j']])))) return(NULL)
		idx = cmp_pix_index(d[['i']], d[['j']], as.integer(n))
		ok = which(!is.na(idx))
		v = numeric(P)
		v[idx[ok]] = d[['value']][ok]
	}
	v[is.na(v)] = 0
	v
}

## gM on the template of target_gr with d's values filled in (0 elsewhere): what
## rebin_matrix's template join + gM() returns. NULL if d needs the reference join.
.fill_template_gm = function(target_gr, d){
	if (!is.integer(d[['i']]) || !is.integer(d[['j']])) return(NULL)
	v = d[['value']]
	if (!(is.double(v) || is.integer(v) || is.logical(v))) return(NULL)
	n = length(target_gr)
	tc = cmp_template_cols(numeric(n), FALSE)
	idx = cmp_pix_index(d[['i']], d[['j']], as.integer(n))
	ok = which(!is.na(idx))
	val = numeric(length(tc$i))
	val[idx[ok]] = as.double(v[ok])
	gM(gr=target_gr,dat=data.table(i = tc$i, j = tc$j, value = val))
}

## ---------------------------------------------------------------------------------------------

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
  ## fast path: draw the samples as make_noisydat() would, and score them against both maps'
  ## values directly (compdats(ss, map$dat, theta) on aligned vectors), without building tables
  maps = list(null_map,hyp_map)
  preps = lapply(maps,function(m){.noisy_prep(m,num_samples)})
  raws = lapply(maps,.gm_rawdat)
  fast = !any(vapply(preps,is.null,logical(1))) && .compdats_fast_ok(theta,FALSE,1e8) &&
	  preps[[1]]$P == preps[[2]]$P && all(vapply(raws,function(d){.is_template(d,length(.tile_widths(maps[[1]]$gr)))},logical(1)))
  if (fast){
    draws = lapply(preps,function(p){.noisy_draw(p,num_samples,theta)}) #samples to classify
    P = preps[[1]]$P
    likelihooddiff = lapply(seq_along(draws),function(s){
	    wp = preps[[s]]$widthprod
	    unlist(setNames(lapply(seq_len(num_samples),function(k){
		    x = draws[[s]][((k-1)*P+1):(k*P)]
		    return(sum(.compdats_values(x,raws[[1]][['value']],wp,theta,FALSE,1e8),na.rm=TRUE)-
			   sum(.compdats_values(x,raws[[2]][['value']],wp,theta,FALSE,1e8),na.rm=TRUE))
	    }),as.character(seq_len(num_samples))))
    })
  }else{
  samples = lapply(list(null_map,hyp_map),function(m){make_noisydat(m,num_samples,theta)}) #samples to classify}
  likelihooddiff = lapply(samples,function(s){
                   unlist(lapply(s,function(ss){
                    return(compdats(ss,null_map$dat,theta)-compdats(ss,hyp_map$dat,theta))
                    }))
                   })
  }
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

call_gwalk_longreads = function(longreads, gwa, gwb, readL, depth){
	longreads_aswords = unlist(lapply(longreads,function(r){paste0(r,collapse='|')}))
	pq = longread_probdist(list(gwa,gwb),readL)
	p_a = pq[[1]]
	p_b = pq[[2]]

	loglikratio = sum(log(p_a[longreads_aswords]) - log(p_b[longreads_aswords]))
	if (loglikratio > 0){
		called_gw = gwa
	}else if (loglikratio < 0){
		called_gw = gwb
	}else{
		called_gw = sample(list(gwa,gwb))
	}
	simdata = liktest_separable_lr(gwa,gwb,readL=readL,depth=depth)
	return(list(called_gw = called_gw, loglikratio = loglikratio, sim_error_rate = simdata$error_rate, gwa_llr = simdata$gwa_ratios, gwb_llr = simdata$gwb_ratios))
}

longread_probdist = function(gwlist,readL,minsize=0,background=1e-10,use.nodes=NULL,mc.cores=1){
	if (is(gwlist,'gWalk')){gwlist = list(gwlist)}
	reads_list = mclapply(gwlist,function(x){reads_fromwalk(x,readL,minsize=minsize,use.nodes=use.nodes)},mc.cores=mc.cores)
	.probdist_from_reads(reads_list,background=background,mc.cores=mc.cores)
}

## longread_probdist() given the reads_fromwalk() tables: smoothed word frequencies over the
## shared vocabulary. Words are looked up with match() rather than by name; "" and NA words
## never match a name, so they get count 0 exactly as the named lookup gives them.
.probdist_from_reads = function(reads_list,background=1e-10,mc.cores=1){
	wordlist = lapply(reads_list,function(x){x$words})
	all_words = unique(unname(unlist(wordlist)))#do.call('union',lapply(reads_list,function(x){x$words}))
	nomatch = is.na(all_words) | all_words == ''
	mclapply(reads_list,function(x){
			m = match(all_words,x$words)
			m[nomatch] = NA_integer_
			cx = unname(x$nums)[m]
			names(cx) = all_words
			cx[is.na(cx)] = 0
			return((cx + background)/(sum(cx) + background * length(all_words)))
		   },mc.cores=mc.cores)
}

liktest_separable_lr = function(gwa,gwb,readL,depth=1,nsamp = 20,mc.cores=1,background=1e-10,use.nodes=NULL,return.kl = F,return.all=F){
	context = sum(width(gwa$graph$footprint)) + 2*readL
	covperread = readL/context
	N = depth/covperread

	pq = longread_probdist(list(gwa,gwb),readL,background=background,use.nodes=use.nodes)
	.liktest_lr_core(pq,N,nsamp,return.kl,return.all)
}

## the sampling part of liktest_separable_lr() given longread_probdist() output
.liktest_lr_core = function(pq,N,nsamp,return.kl=F,return.all=F){
	p = pq[[1]]
	q = pq[[2]]

	if(return.kl){return(N*sum(p*log(p/q)))}

	nm = names(p)
	if (length(p) >= 1 && !is.null(nm) && identical(nm, names(q))){
		## sample(names(p), N, prob = p, replace = T) is names(p)[sample.int(length(p), N, T, p)], and
		## p[names] / q[names] with unique names is positional, except that "" and NA never match
		lr = log(p) - log(q)
		lr[is.na(nm) | nm == ''] = NA_real_
		np = length(p)
		p_samples = lapply(1:nsamp,function(x){sample.int(np,size = N,replace=T,prob=p)})
		q_samples = lapply(1:nsamp,function(x){sample.int(np,size = N,replace=T,prob=q)})
		loglikratio = function(sample){sum(lr[sample])}
	}else{
		p_samples = lapply(1:nsamp,function(x){sample(names(p),size = N,prob=p,replace=T)})
		q_samples = lapply(1:nsamp,function(x){sample(names(q),size = N,prob=q,replace=T)})
		loglikratio = function(sample){sum(log(p[sample]) - log(q[sample]))}
	}

	p_ratios = unlist(lapply(p_samples,loglikratio))
	q_ratios = unlist(lapply(q_samples,loglikratio))
	error_rate = (sum(p_ratios <= 0) + sum(q_ratios >= 0)) / (2*nsamp)
	if (return.all){return(list(error_rate = error_rate, gwa_ratios = p_ratios, gwb_ratios = q_ratios))}
	return(error_rate)
}

liktest_separable = function(map_a,map_b,theta=2,nsamp = 20,mc.cores=1,error.thresh = NULL){
	## fast path: the samples are only scored against map_a and map_b (compdats with its defaults,
	## theta = 2 and area0 = 1e8), so draw them as make_noisydat() does and score the value vectors
	prep_a = .noisy_prep(map_a,nsamp)
	prep_b = .noisy_prep(map_b,nsamp)
	raw_a = .gm_rawdat(map_a)
	raw_b = .gm_rawdat(map_b)
	if (!is.null(prep_a) && !is.null(prep_b) && prep_a$P == prep_b$P &&
	    .is_template(raw_a,length(.tile_widths(map_a$gr))) && .is_template(raw_b,length(.tile_widths(map_a$gr)))){
		draws_a = .noisy_draw(prep_a,nsamp,theta)
		draws_b = .noisy_draw(prep_b,nsamp,theta)
		P = prep_a$P
		ya = raw_a[['value']]
		yb = raw_b[['value']]
		nll = function(x,y,wp){sum(.compdats_values(x,y,wp,2,FALSE,1e8),na.rm=TRUE)}
		ratios = function(draws,wp){
			unlist(mclapply(seq_len(nsamp),function(k){
				x = draws[((k-1)*P+1):(k*P)]
				-nll(x,ya,wp)+nll(x,yb,wp)
			},mc.cores=mc.cores))
		}
		a_ratios = ratios(draws_a,prep_a$widthprod)
		b_ratios = ratios(draws_b,prep_b$widthprod)
	}else{
	a_samples = make_noisydat(map_a,nsamp=nsamp,theta=theta)
	b_samples = make_noisydat(map_b,nsamp=nsamp,theta=theta)
	#c,maskzeros=Tompmaps is the negative log likelihood
	loglikratio = function(sample){-compdats(sample,map_a$dat)+compdats(sample,map_b$dat)}

	a_ratios = unlist(mclapply(a_samples,loglikratio,mc.cores=mc.cores))
	b_ratios = unlist(mclapply(b_samples,loglikratio,mc.cores=mc.cores))
	}
	error_rate = (sum(a_ratios <= 0) + sum(b_ratios >= 0)) / (2*nsamp)
	return(error_rate)
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

compmaps = function(map_exp,map_theory,theta=2,ifsum=F,if.diag=T,ifscale=F,mask=NULL,return_kl=F,area0=NULL,maskzeros=F){
   if (is.null(area0)){
	gr = map_theory$gr
        medianwid = median(width(gr))
	area0 = medianwid^2
   }
   map_theory = rebin_matrix(map_theory,map_exp$gr)
   map_exp = rebin_matrix(map_exp,map_exp$gr)
   ## Fast path: both maps on the same template (the usual case after rebin_matrix). This is
   ## compdats(map_exp$dat, map_theory$dat, ...) on the aligned value vectors, without the two
   ## table copies, the setkey, the all.equal check and the merge.
   dat_exp = .gm_rawdat(map_exp)
   dat_theory = .gm_rawdat(map_theory)
   if (.compdats_fast_ok(theta,return_kl,area0) && .aligned_pair(dat_exp,dat_theory)){
       res = .compdats_aligned(dat_exp,dat_theory,theta,if.diag,return_kl,area0)
       value = res$value
       if (!is.null(mask)){
           bad.inds = (map_exp$gr %&% mask)$tile.id
           ii = dat_exp[['i']]; jj = dat_exp[['j']]
           if (!is.null(res$keep)){ii = ii[res$keep]; jj = jj[res$keep]}
           value[(ii %in% bad.inds) & (jj %in% bad.inds)] = 0
           res$value = value
       }
       if (ifsum){
           return(sum(value))
       }else{
           return(gM(gr=map_exp$gr,dat=.compdats_table(dat_exp,res)))
       }
   }
   comp_dat = compdats(map_exp$dat,map_theory$dat,theta,ifsum=FALSE,if.diag=if.diag,return_kl=return_kl,area0=area0)
   if (!is.null(mask)){
       bad.inds = (map_exp$gr %&% mask)$tile.id
       comp_dat[(i %in% bad.inds) & (j %in% bad.inds),value:=0]
   }
   if (ifsum){
       return(sum(comp_dat$value))
   }else{
       return(gM(gr=map_exp$gr,dat=comp_dat))
   }
}

estimate.depthratio <- function(filepath,res=1e6,ploidy=2,if.chr=FALSE,purity=1){ #put here hic data for a whole genome
    if (!is.numeric(purity) || length(purity) != 1 || is.na(purity) || purity < 0 || purity > 1){
        stop('purity must be a single numeric value between 0 and 1')
    }
    if (!is.numeric(ploidy) || length(ploidy) != 1 || is.na(ploidy) || ploidy <= 0){
        stop('ploidy must be a single positive numeric value')
    }
    mode = tools::file_ext(filepath)
    if (mode=='mcool'){ #chromosome lengths from the file itself, so any build works (e.g. hg38 Pore-C; hg19 files give the hg19 lengths)
        chroms = rhdf5::h5read(filepath, paste0('/resolutions/', as.integer(res), '/chroms'))
        main = grepl('^(chr)?([0-9]+|X|Y)$', chroms$name)
        wholegenome = gr.tile(si2gr(setNames(as.numeric(chroms$length[main]), sub('^chr', '', chroms$name[main]))),res)
    }else{
        wholegenome = gr.tile(si2gr(hg_seqlengths()[1:24]),res)
    }
    if (if.chr){
        wholegenome = gr.chr(wholegenome)
    }
    if (mode=='hic'){
        depthest.hic = straw(filepath,res=as.integer(res),gr=wholegenome)
    }else if (mode=='mcool'){
        depthest.hic = cooler(filepath,res=as.integer(res),gr=wholegenome)
    }else if (mode=='rds'){
        depthest.hic = readRDS(filepath)
    }else{
        stop('Unrecognized Hi-C file mode: ',mode)
    }
    sample.ploidy = purity*ploidy + (1-purity)*2
    return(depthest.hic$value%>%sum * 300 / 3e9 * 2/sample.ploidy)
}

#returns negative log likelihood of test data given true data
compdats = function(dat_test,dat_true,theta=2,ifsum=TRUE,if.diag=TRUE,return_kl=F,area0=1e8){
    setkey(dat_test,'id')
    setkey(dat_true,'id')
    ## fast path: same template on both sides (identical i/j/id after the sort), so the merge
    ## below is a row-by-row pairing and the all.equal() check passes
    if (.compdats_fast_ok(theta,return_kl,area0) && .aligned_pair(dat_test,dat_true)){
        res = .compdats_aligned(dat_test,dat_true,theta,if.diag,return_kl,area0)
        if (ifsum==TRUE){
            return(sum(res$value,na.rm=TRUE))
        }else{
            return(.compdats_table(dat_test,res))
        }
    }
    if (is.null(dat_test$widthprod)){
        dat_test$widthprod= area0 #if the data.table doesn't have an area column, assume all pixels are the same size
    }
    if (all.equal(dat_test[,.(i,j,id)],dat_true[,.(i,j,id)])!=T){
	    stop('Equalize gMatrices before comparing data')}
    combdat = merge.data.table(dat_test[,.(i,j,id,widthprod,value)],dat_true[,.(id,value)],by=c('id'))
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
    if (return_kl){
	    combdat[,value:=widthprod/area0*kldiv]
    }else{
	    combdat[,value:=-logprob*widthprod/area0]
    }
    if (ifsum==TRUE){
        return(sum(combdat$value,na.rm=TRUE))
    }else{
	return(combdat[,.(value,i,j,id)])
    }
}

make_noisymap <- function(map_in,nsamp=1,theta=2,mc.cores=1){
    dats = make_noisydat(map_in,nsamp,theta)
    out.maps = mclapply(dats,function(i){gM(gr=map_in$gr,dat=i)},mc.cores=mc.cores)
    return(out.maps)
}

make_noisydat = function(map_in,nsamp=1,theta=2){ #samples from negative binomial distribution, unless theta=0 then uses poisson
    prep = .noisy_prep(map_in,nsamp)
    if (is.null(prep)){return(.make_noisydat_ref(map_in,nsamp,theta))}
    ## same single vectorized draw as the reference, then one table per sample sliced out of it
    ## (instead of rbind-ing nsamp copies of the template and split()-ing them again)
    newval = .noisy_draw(prep,nsamp,theta)
    P = prep$P
    out.dats = lapply(seq_len(nsamp),function(k){
        fresh = function(v){if (k==1) v else copy(v)} #each table owns its columns (one copy each, no data.table() copies)
        setDT(list(i=fresh(prep$i),j=fresh(prep$j),widthprod=fresh(prep$widthprod),id=fresh(prep$id),
                   value=newval[((k-1)*P+1):(k*P)],map.ids=rep(k,P)))
    })
    names(out.dats) = as.character(seq_len(nsamp))
    return(out.dats)
}

## reference make_noisydat(), for inputs the fast path does not take
.make_noisydat_ref = function(map_in,nsamp=1,theta=2){
    mapdat = map_in$dat
    template = make_template_dat(map_in$gr)
    mapdat = merge.data.table(template[,.(i,j,widthprod,id)],mapdat[,.(i,j,value)],all.x=TRUE,by=c('i','j'))
    meanvals = mapdat$value
    meanvals[is.na(meanvals)] = 0
    if (theta==0){
        newval = rpois(length(meanvals)*nsamp,rep(meanvals,nsamp))
    }else {
        newval = rnbinom(length(meanvals)*nsamp,size=theta,mu=rep(meanvals,nsamp))
    }
    multimap = do.call(rbind,rep(list(mapdat),nsamp))
    multimap$map.ids = lapply(1:nsamp,function(i){rep(i,nrow(mapdat))}) %>% unlist
    multimap$value = newval
    out.dats = split(multimap,by='map.ids')
    #out.gms = lapply(out.dats,function(i){gM(gr=map_in$gr,dat=i)})
    return(out.dats)
}

hic_cov = function(gm){ #get coverage from Hi-C map, at the granges in that map
	gr = gm$gr
	gr$tile.id=1:length(gr)
	grdt = gr2dt(gr)
	## fast path: per-tile sums over the dense template-ordered values, no template merge
	v = .dense_values(.gm_rawdat(gm),length(gr))
	if (!is.null(v)){
		out = data.table(i = seq_len(length(gr)), sum.value = cmp_tile_cov(v,length(gr)))
	}else{
	temp = make_template_dat(gr)
	dat = merge.data.table(gm$dat,temp[,.(i,j)],by=c('i','j'),all.y=T)
	dat[is.na(value),value:=0]
	out = rbind(
	  dat[, .(i, value)],
	  dat[i != j, .(i = j, value)]
	)[, .(sum.value = sum(value)), by = i]
	}
	grdt = merge.data.table(grdt,out[,.(tile.id=i,cov=sum.value)],by='tile.id')
	return(dt2gr(grdt,seqlengths = seqlengths(gm$gr)))
}

callloops = function(hic,gw,depth,fdr_cut = 0.01,resolution = NULL,sim.dat = NULL,reg = NULL,qqplot_path = NULL,mode='peaks'){ #call Hi-C loops 
	if (is.null(reg)){
		reg = streduce(hic$gr)
	}else{
		hic.gr2 = hic$gr %&% reg
		hic = hic$disjoin(hic.gr2)$agg(hic.gr2) #FIX THIS: should rebin sims, not true data
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

group_klcomp = function(gw1, gw2=NULL,target_region = NULL, hic_resolution = NULL,depth=10,mc.cores=1){
	if (is.null(target_region)){ft = gw1[[1]]$footprint
	}else{ft = target_region}
	symm=F
	if (is.null(gw2)){
		gw2=gw1
		symm=T
	} #compare among selves
	if (is.null(hic_resolution)){
		w = sum(width(ft))
		hic_resolution = 1e3*round(w/1e5) #about 100 bins across the footprint,rounded off to the nearest kb
	}
	comppairs = CJ(i=seq_along(gw1),j=seq_along(gw2))
	if (symm){comppairs = comppairs[j>i]}
	## Each pair's KL is hic_kl(gw1[[i]], gw2[[j]], target_region, pix.size = hic_resolution,
	## depth, theta = 2): both walks simulated on target_region (or, if NULL, on the footprint of
	## gw1[[i]]'s graph), then compmaps(KL). Simulate every (walk, target) once and reuse the maps.
	same_hash = vapply(seq_len(nrow(comppairs)),function(x){gw1[[comppairs$i[x]]]$hash==gw2[[comppairs$j[x]]]$hash},logical(1))
	todo = which(!same_hash)
	has_graph = all(vapply(c(gw1[unique(comppairs$i[todo])],gw2[unique(comppairs$j[todo])]),function(g){!is.null(g$graph)},logical(1)))
	if (!length(todo) || !has_graph){
		kldiv = mclapply(1:nrow(comppairs),function(x){
			i = comppairs[x]$i
			j = comppairs[x]$j
			if (gw1[[i]]$hash==gw2[[j]]$hash){return(0)}
			hic_kl(gw1[[i]],gw2[[j]],target_region = target_region,pix.size=hic_resolution,depth=depth,theta=2)
		},mc.cores=mc.cores)
		comppairs$kldiv = unlist(kldiv)
		return(comppairs)
	}
	## target of each pair, as an index into a list of distinct targets
	if (is.null(target_region)){
		targets = list()
		tix = integer(nrow(comppairs))
		for (x in todo){
			fp = gw1[[comppairs$i[x]]]$graph$footprint
			k = Position(function(t){identical(t,fp)},targets)
			if (is.na(k)){targets[[length(targets)+1]] = fp; k = length(targets)}
			tix[x] = k
		}
	}else{
		targets = list(target_region)
		tix = rep(1L,nrow(comppairs))
	}
	## the simulations needed: (set, walk index, target); in the symmetric case both sides are gw1
	side2 = if (symm) 'a' else 'b'
	jobs = unique(rbind(data.table(set='a',w=comppairs$i[todo],t=tix[todo]),
			    data.table(set=side2,w=comppairs$j[todo],t=tix[todo])))
	jobs[,key:=paste(set,w,t)]
	sims = mclapply(seq_len(nrow(jobs)),function(k){
		gw = if (jobs$set[k]=='a') gw1[[jobs$w[k]]] else gw2[[jobs$w[k]]]
		forward_simulate(gw,target_region = targets[[jobs$t[k]]],pix.size=hic_resolution,depth=depth)
	},mc.cores=mc.cores)
	names(sims) = jobs$key
	kldiv = mclapply(1:nrow(comppairs),function(x){
		if (same_hash[x]){return(0)}
		hic_x = sims[[paste('a',comppairs$i[x],tix[x])]]
		hic_y = sims[[paste(side2,comppairs$j[x],tix[x])]]
		compmaps(hic_x,hic_y,theta=2,return_kl=T,area0=hic_resolution^2,ifsum=T,mask=NULL)
	},mc.cores=mc.cores)
	comppairs$kldiv = unlist(kldiv)
	return(comppairs)
}

#what is the minimal readlength required to tell the difference between two karyotypes?
minimal_readlength = function(gw1,gw2){
	if (gw1==gw2){return(NA)}
	count_kmers = function(walks, k) {
	  as.data.table(count_kmers_cpp(snode_ids = walks$snode.id,circular = walks$circular,k = as.integer(k)))}
	maxl = max(c(gw1$dt$length,gw2$dt$length))
	for (k in 3:maxl){
		k1 = count_kmers(gw1,k) %>% setkey('kmer')
		k2 = count_kmers(gw2,k) %>% setkey('kmer')
		compk = merge.data.table(k1[,.(count,kmer)],k2[,.(count,kmer)],by='kmer',all=T)
		compk[is.na(count.x),count.x:=0]
		compk[is.na(count.y),count.y:=0]
		diffk = compk[count.x!=count.y]
		if (nrow(diffk)){break}
	}
	keepcols = c(paste0('node',c(2:(k-1))))
	midnodes = unique(rbind(k1[diffk$kmer,nomatch=0],k2[diffk$kmer,nomatch=0]))[,..keepcols]
	width_lookup = setNames(gw1$graph$nodes$dt$width, gw1$graph$nodes$dt$node.id)
	midnodes[, width_sum := rowSums(matrix(width_lookup[as.character(abs(as.matrix(.SD)))],nrow = .N),na.rm = TRUE), .SDcols = keepcols]
	return(min(midnodes$width_sum))
}

.same_gr_bins = function(gr1, gr2){
	length(gr1)==length(gr2) &&
		identical(as.character(seqnames(gr1)), as.character(seqnames(gr2))) &&
		identical(start(gr1), start(gr2)) &&
		identical(end(gr1), end(gr2))
}

rebin_matrix = function(gm,target_gr,writezeros=T){
	## Fast path: if the gMatrix is already defined on the requested bins, avoid
	## the relatively expensive disjoin()/agg() calls.  Still fill missing pixels
	## with explicit zeroes when needed so compdats() sees identical templates.
	## The pixel table is read in place (gm$dat copies it on every access) and checked
	## against the template layout without building make_template_dat().
	if (.same_gr_bins(gm$gr,target_gr)){
		d = .gm_rawdat(gm)
		if (.is_template(d,length(target_gr))){
			return(gm)
		}
		if (all(c('i','j','value') %in% names(d)) && !anyDuplicated(as.data.table(list(i=d[['i']],j=d[['j']])))){
			out = .fill_template_gm(target_gr,d)
			if (!is.null(out)){return(out)}
			outdt = make_template_dat(target_gr)
			setkey(outdt,i,j)
			outdt[gm$dat, on = .(i, j), value := i.value]
			return(gM(gr=target_gr,dat=outdt))
		}
	}
	## Different bins: split counts by overlap (not gm$disjoin()$agg(), which copies each pixel's
	## full value to every piece it is cut into and then sums the pieces, inflating the map)
	d2 = .rebin_counts(gm,target_gr)
	out = .fill_template_gm(target_gr,d2)
	if (!is.null(out)){return(out)}
	outdt = make_template_dat(target_gr)
	setkey(outdt,i,j)
	outdt[d2, on = .(i, j), value := i.value]
	return(gM(gr=target_gr,dat=outdt))
}

## Count-conserving rebin of gm's upper-triangle pixels onto target_gr. Each source bin is cut by
## the target bins, with f = overlap width / source bin width per piece, and a source pixel's value
## goes to the target pixels of its pieces in proportion to f_i * f_j (counts assumed uniform within
## a pixel). Ordered piece pairs are summed, so a split diagonal pixel gives 2 f_a f_b to its (a,b)
## target pixel. Source bins nested in target bins are not cut: coarsening onto them is an exact sum.
## Parts of the source outside target_gr are dropped. Returns data.table(i, j, value) with i <= j.
.rebin_counts = function(gm,target_gr){
	src = gm$gr
	d = .gm_rawdat(gm)
	ov = gr2dt(gr.findoverlaps(src,target_gr))
	if (!nrow(ov) || !nrow(d)){
		return(data.table(i=integer(),j=integer(),value=numeric()))
	}
	ov = ov[,.(s=query.id,t=subject.id,f=width/as.numeric(width(src))[query.id])]
	pix = data.table(i=as.integer(d[['i']]),j=as.integer(d[['j']]),value=as.double(d[['value']]))
	pix = merge(pix,ov[,.(i=s,ti=t,fi=f)],by='i',allow.cartesian=TRUE)
	pix = merge(pix,ov[,.(j=s,tj=t,fj=f)],by='j',allow.cartesian=TRUE)
	pix[,.(value=sum(value*fi*fj,na.rm=TRUE)),by=.(i=pmin(ti,tj),j=pmax(ti,tj))]
}
