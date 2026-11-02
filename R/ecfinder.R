fastKar_ecfinder = function(gg,ft,hic.res,true_hic_path,ec_fraction_threshold = 0.6,n_sample=100,n_ec = 10,n_circle_samples=1000,n_hsr_samples=5000,figure_path = NULL,mc.cores=1,return_all=F,purity=1,rescale=F,score_region=NULL,seed=NULL,depth=NULL,sep_D_min=0.9,res_ladder=c(5e3,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6),keep.acen=TRUE,centromeres=NULL,hypotheses=c('independent','paired'),n_boil=3,call_q=0.05){
	library(ggforce)
	# true_hic_path = NULL: no Hi-C. Only asks whether the graph allows ecDNA and how separated the simulated ecDNA and HSR
	# score distributions are (sep_D, sep_margin), at each Hi-C depth in 'depth' (required; X coverage as from
	# estimate.depthratio). Bins tile the window at the res_ladder resolution closest to hic.res; no rescale, no observed
	# score, no plots; call 'None'. Hypotheses are simulated once at depth[1] and scaled (simulated maps are linear in depth)
	no_hic = is.null(true_hic_path)
	if (no_hic && is.null(depth)){stop('fastKar_ecfinder: depth is required when true_hic_path is NULL')}
	# seed: reproducible run (L'Ecuyer-CMRG streams, so forked mclapply workers are seeded too; same mc.cores needed)
	if (!is.null(seed)){
		old_kind = RNGkind("L'Ecuyer-CMRG")
		set.seed(seed)
		on.exit(RNGkind(old_kind[1]),add=TRUE)
	}
	# return_all = TRUE: also return the graph, hypotheses and observed/simulated maps (element 'hyp'), e.g. for figures
	# sep_D_min: a Hi-C call (ecDNA or CI) needs the simulated ecDNA and HSR score distributions to be separated, KS p < 0.01
	# and KS D > sep_D_min (the p-value alone passes heavy overlap at 100 samples per side, e.g. D = 0.5); otherwise 'None'.
	# 0.9 since 2026-11-02 (user; was >= 0.95): COLO320-HSR with the contig fold-back, D 0.906, observed at the HSR mode
	# hypotheses: 'independent' (default): ecDNA hypotheses boiled from the whole graph (boil_greedy), HSR hypotheses from
	# squeeze(), generated separately. 'paired': each ecDNA hypothesis is one of the squeeze() HSR hypotheses with its
	# amplicon-touching walks boiled into circles (boil_paired(), n_boil boils each), identical outside the amplicon; noisy
	# samples are then scored without their whole group. Kept as an option to compare on more lines (user, 2026-11-02).
	# call_q: separated scores call ecDNA if the observed score is above the 1 - call_q quantile of the HSR samples, CI if
	# below the call_q quantile of the ecDNA samples, else None (was: above the ecDNA minimum / below the HSR maximum)
	# keep.acen / centromeres: passed to squeeze(). keep.acen = TRUE (default; user, 2026-11-02): acentric HSR hypotheses
	# allowed (with the paired ecDNA hypotheses, the centric filter flipped COLO320-HSR to ecDNA in 2 of 3 seeds; centric
	# filtering is kept for ranking HSR hypotheses later). keep.acen = FALSE: HSR hypotheses must carry the amplicon on a
	# centric walk, and an event with no such karyotype is called 'None' ('No centric HSR solutions')
	# chromosome lengths from the graph, not DEFAULT_GENOME / DEFAULT_BSGENOME (see graph_seqlengths())
	sl = GenomeInfoDb::seqlengths(graph_seqlengths(gg$nodes$gr, gg))
	sl = sl[grepl('^(chr)?([0-9]+|X|Y)$', names(sl)) & !is.na(sl)]
	wholegenome = GRanges(names(sl), IRanges(1, sl))
	if (is(ft,'character')){ft = parse.gr(ft)}
	ft = streduce(graph_seqlengths(ft, gg))
	ft_context = streduce(ft + sum(width(ft))/10) #add 10% to the footprint for context
	# score_region: Hi-C window to score (and plot) instead of the footprint + 10%, e.g. a few Mb around the amplified
	# gene of a large multi-chromosome event, so bins are finer; hypotheses are still sampled on the whole footprint
	if (!is.null(score_region)){
		if (is(score_region,'character')){score_region = parse.gr(score_region)}
		ft_context = streduce(score_region)
	}
	event_nodes = (gg$nodes$gr %&% ft)$node.id
	message('Sampling and generating solutions')
	hypotheses = match.arg(hypotheses)
	amp_nodes = unique(event_nodes)
	# diagnostics (never filters); centromeres only needed for the centric fractions, toy graphs may not have a table
	cen_nodes = tryCatch({
		cen = if (is.null(centromeres)) fastKar_centromeres(gg) else centromeres
		unique((gg$nodes$gr %&% cen)$node.id)
	},error=function(e){NULL})
	diagnostics = function(hsr_solns,ecdna_solns,ec_group,pair_outside_match){
		hsr_cov = vapply(hsr_solns,amp_linear_coverage,numeric(1),ft=ft)
		list(n_hsr=length(hsr_solns),n_groups=if (is.null(ec_group)) NA_integer_ else length(unique(ec_group)),
		     pair_match_frac=if (length(pair_outside_match)) mean(pair_outside_match) else NA_real_,
		     hsr_top_cov=if (length(hsr_cov)) max(hsr_cov) else NA_real_,
		     hsr_centric_frac=if (is.null(cen_nodes) || !length(hsr_solns)) NA_real_ else mean(vapply(hsr_solns,amp_centric_frac,numeric(1),ft=ft,cen_nodes=cen_nodes)),
		     ec_centric_frac=if (is.null(cen_nodes) || !length(ecdna_solns)) NA_real_ else mean(vapply(ecdna_solns,amp_centric_frac,numeric(1),ft=ft,cen_nodes=cen_nodes)),
		     hsr_telo_frac=if (length(hsr_solns)) mean(vapply(hsr_solns,amp_telo_frac,numeric(1),ft=ft)) else NA_real_)
	}
	ec_group = NULL; pair_outside_match = logical(0); pair_outside_frac_bp = numeric(0)
	if (hypotheses == 'independent'){
		# ecDNA and HSR hypotheses generated separately: ecDNA by boiling the whole graph (boil_greedy), HSR by squeeze()
		ecdna_solns = boil_greedy(gg=gg,ft=ft,N=n_circle_samples,n_solns_try=n_ec,mc.cores=mc.cores) #karyotypes sampled to find circles
		ecdna_circ = unname(unlist(lapply(ecdna_solns,function(x){amp_circular_fraction(x,ft)})))
		ecdna_all = ecdna_solns
		ecdna_solns = ecdna_solns[ecdna_circ > ec_fraction_threshold]
		if (!length(ecdna_solns)){
			out = c(list(call='CI', reason='No good ecDNA solutions', max_ecdna_circ=max(c(0,ecdna_circ),na.rm=TRUE)),diagnostics(list(),list(),NULL,logical(0)))
			if (return_all){out$hyp = list(gg=gg,ft=ft,ecdna_candidates=ecdna_all,ecdna_circ=ecdna_circ)}
			return(out)
		}
		hsr_solns = squeeze(gg=gg,ft=ft,N=n_hsr_samples,k_return=n_ec,verbose = F,mc.cores=mc.cores,keep.acen=keep.acen,centromeres=centromeres) #karyotypes sampled to pick the HSR solutions from
		if (!length(hsr_solns)){
			out = list(call='None', reason='No centric HSR solutions', max_ecdna_circ=max(c(0,ecdna_circ),na.rm=TRUE))
			if (return_all){out$hyp = list(gg=gg,ft=ft,ecdna_candidates=ecdna_all,ecdna_circ=ecdna_circ)}
			return(out)
		}
	} else{
		# paired: each ecDNA hypothesis is an HSR hypothesis with its amplicon-touching walks boiled into circles, everything
		# outside the amplicon unchanged (docs/ecdna_paired_hypotheses_plan.md in the context repo)
		hsr_solns = squeeze(gg=gg,ft=ft,N=n_hsr_samples,k_return=n_ec,verbose = F,mc.cores=mc.cores,keep.acen=keep.acen,centromeres=centromeres)
		if (!length(hsr_solns)){
			out = list(call='None', reason='No centric HSR solutions', max_ecdna_circ=NA_real_)
			if (return_all){out$hyp = list(gg=gg,ft=ft)}
			return(out)
		}
		circles = find_circles(gg,ft,N=n_circle_samples,mc.cores=mc.cores) #once, on the full graph; score_circles drops those that don't fit a hypothesis' CN
		boils = mclapply(seq_along(hsr_solns),function(k){boil_paired(hsr_solns[[k]],circles,ft,n_boil=n_boil)},mc.cores=mc.cores)
		ecdna_all = unlist(boils,recursive=FALSE)
		if (is.null(ecdna_all)){ecdna_all = list()}
		ec_group_all = rep(seq_along(boils),base::lengths(boils)) #index of the HSR hypothesis each ecDNA candidate was boiled from
		ecdna_circ = unname(unlist(lapply(ecdna_all,function(x){amp_circular_fraction(x,ft)})))
		if (is.null(ecdna_circ)){ecdna_circ = numeric(0)}
		keep = ecdna_circ > ec_fraction_threshold
		ecdna_solns = ecdna_all[keep]
		ec_group = ec_group_all[keep]
		# pairing check: ecDNA hypothesis vs its HSR partner, outside the amplicon
		hsr_sig = lapply(hsr_solns,outside_signature,amp_nodes=amp_nodes)
		pair_mm = lapply(seq_along(ecdna_solns),function(j){outside_mismatch(outside_signature(ecdna_solns[[j]],amp_nodes),hsr_sig[[ec_group[j]]])})
		pair_outside_match = vapply(pair_mm,function(m){m$match},logical(1))
		pair_outside_frac_bp = vapply(pair_mm,function(m){m$frac_bp},numeric(1))
		if (!length(ecdna_solns)){
			out = c(list(call='CI', reason='No good ecDNA solutions', max_ecdna_circ=max(c(0,ecdna_circ),na.rm=TRUE)),
				diagnostics(hsr_solns,list(),ec_group,logical(0)))
			if (return_all){out$hyp = list(gg=gg,ft=ft,ecdna_candidates=ecdna_all,ecdna_circ=ecdna_circ,hsr_solns=hsr_solns,ec_group=ec_group_all,
						      hsr_cov=vapply(hsr_solns,amp_linear_coverage,numeric(1),ft=ft))}
			return(out)
		}
	}
	diag = diagnostics(hsr_solns,ecdna_solns,ec_group,pair_outside_match)
	hsr_cov = vapply(hsr_solns,amp_linear_coverage,numeric(1),ft=ft)
	n_ec = length(ecdna_solns)
	n_hsr = length(hsr_solns)
	leave_out_mode = pick_leave_out_mode(ec_group,n_ec,n_hsr)
	random_walks = sample.gwalks(gg,n_sample,mc.cores=mc.cores,verbose = F)
	random_circ = unname(unlist(mclapply(
		random_walks,function(x){amp_circular_fraction(x,ft)},mc.cores=mc.cores)))
	# estimate depth of hi-c data and pick resolution to simulate at
	ploidy = sum(gg$gr$cn*width(gg$gr))/sum(width(gg$gr))
	if (hic.res<0){ #if hic.res < 0 this is asking for a number of bins to cover the footprint
		hic.res = sum(width(ft_context))/abs(hic.res)
	}
	if (no_hic){
		depths = depth
		depth = depths[1]
		hic.res = res_ladder[which.min(abs(log10(res_ladder) - log10(hic.res)))]
		hic.gr = (gr.tile(wholegenome,hic.res) %&% ft_context)[,c()] #genome-aligned bins, as a cooler query returns
	} else{
	is.chr = hic_is_chr(true_hic_path) #e.g. the NYU cell line mcools name chromosomes 'chr1', ...; graphs don't
	depth = estimate.depthratio(true_hic_path,ploidy=ploidy,purity=purity,if.chr=is.chr)
	hictype = tools::file_ext(true_hic_path)
	# now find available resolution that is closest to requested resolution
	avail_res = hic_res(true_hic_path)
	dists = abs(log10(avail_res) - log10(hic.res))
	hic.res = avail_res[dists==min(dists)]
	# load hic data
	if (hictype=='mcool'){
		if (is.chr){
			true_hic = cooler(true_hic_path,gr=gr.chr(ft_context),res=hic.res)
			true_hic = gM(gr=gr.nochr(true_hic$gr),dat=true_hic$dat)
		} else{true_hic = cooler(true_hic_path,gr=ft_context,res=hic.res)}
	} else{true_hic = straw(true_hic_path,gr=ft_context,res=hic.res)}
	hic.gr = true_hic$gr
	true_hic = rebin_matrix(true_hic,hic.gr) #make sure the matrix matches future matrices (basically just adds in 0s)
	}
	message('Simulating Hi-C for hypotheses and random samples')
	# simulate on the data bins themselves, so every simulated tile lies inside one bin and the rebin is an exact sum.
	# all solutions are walks on gg, so the bins are tiled once and reused (walks on another graph are tiled anew)
	prepped = prep_for_sim(list(graph=gg),hic.gr,hic.res)
	node_cn = gg$nodes$dt$cn
	simfun=function(x){
		p = if (identical(x$graph$nodes$dt$cn,node_cn)) prepped else NULL
		rebin_matrix(forward_simulate(x,target_region = hic.gr,
					      pix.size=hic.res,depth=depth,purity=purity,prepped=p),hic.gr)
	}
	hsr_sims = mclapply(hsr_solns,simfun,mc.cores=mc.cores)
	ecdna_sims = mclapply(ecdna_solns,simfun,mc.cores=mc.cores)
	random_sims = mclapply(random_walks,simfun,mc.cores=mc.cores)
	# rescale = TRUE: the genome-wide depth estimate can be off locally, so scale every sim by one shared
	# ratio = observed total / median total over all hypothesis and random sims (as in bestfit_search)
	ratio = 1
	if (rescale && !no_hic){
		ratio = sim_depth_ratio(true_hic,c(hsr_sims,ecdna_sims,random_sims),mc.cores=mc.cores)
		hsr_sims = lapply(hsr_sims,function(s){s*ratio})
		ecdna_sims = lapply(ecdna_sims,function(s){s*ratio})
		random_sims = lapply(random_sims,function(s){s*ratio})
	}
	# noise and null scores for one set of sims: noisy samples of every solution scored against all solutions
	score_null = function(hsr_sims,ecdna_sims,random_sims){
		#simulate noise, 1 sample per random walk solution and n_sample for each ecDNA/HSR solution
		n_ec_sample = n_sample %/% n_ec
		n_hsr_sample = n_sample %/% n_hsr
		message('Simulating noise')
		hsr_noisy_sims = unlist(mclapply(hsr_sims,function(x){make_noisydat(x,n_hsr_sample)},
						 mc.cores=mc.cores),recursive=F)
		ec_noisy_sims = unlist(mclapply(ecdna_sims,function(x){make_noisydat(x,n_ec_sample)},
						mc.cores=mc.cores),recursive=F)
		random_noisy_sims = mclapply(random_sims,function(x){make_noisydat(x)[[1]]},mc.cores=mc.cores)
		#
		message('Scoring karyotypes')
		area0 = median(width(ecdna_sims[[1]]$gr))^2
		# pixel areas from the bins: noisy samples carry widthprod, the observed gMatrix table does not
		bin_widths = width(hic.gr)
		# ecDNA score is mean NLL over HSR solutions - mean NLL over ecDNA solutions (mean, not best-fit min: as accurate or
		# better in the 2026-10-19 K-sweep on simulated truths, with fewer CI truths called ecDNA at large n_ec)
		# pixel tables pulled out once: gm$dat copies the table on every access
		hsr_dats = lapply(hsr_sims,function(sim){sim$dat})
		ecdna_dats = lapply(ecdna_sims,function(sim){sim$dat})
		nll = function(dats,data,leave_out=NULL){
			nlls = unlist(lapply(dats,function(d){compdats(data,d,area0=area0,widths=bin_widths)}))
			if (length(leave_out)){nlls = nlls[-leave_out]} #guard: x[-integer(0)] would drop everything
			mean(nlls)}
		combined_llr = function(x,leave_hsr=NULL,leave_ec=NULL){nll(hsr_dats,x,leave_hsr)-nll(ecdna_dats,x,leave_ec)}
		# noisy samples are scored without the hypotheses that generated them, since the real data's karyotype is not
		# among them either. 'group' (ecDNA hypotheses from > 1 HSR hypothesis): without HSR k and every ecDNA hypothesis
		# boiled from it; 'solution': without the one solution (each side has > 1; switching on at > 3 left the null
		# thresholds too far out for 2-3 solutions, so those calls mostly went to None in the K-sweep); 'none' otherwise
		ec_src = rep(seq_len(n_ec),each=n_ec_sample)
		hsr_src = rep(seq_len(n_hsr),each=n_hsr_sample)
		# score all samples and true Hi-C
		random_scores = mclapply(random_noisy_sims,combined_llr,mc.cores=mc.cores) %>% unlist
		ec_scores = mclapply(seq_along(ec_noisy_sims),function(k){
			lo = leave_out_ix('ec',ec_src[k],ec_group,leave_out_mode)
				combined_llr(ec_noisy_sims[[k]],leave_hsr=lo$hsr,leave_ec=lo$ec)},mc.cores=mc.cores) %>% unlist
		hsr_scores = mclapply(seq_along(hsr_noisy_sims),function(k){
			lo = leave_out_ix('hsr',hsr_src[k],ec_group,leave_out_mode)
				combined_llr(hsr_noisy_sims[[k]],leave_hsr=lo$hsr,leave_ec=lo$ec)},mc.cores=mc.cores) %>% unlist
		#
		dt = rbind(data.table(score=random_scores,type='Random',circle = random_circ),
			   data.table(score=ec_scores,type='ecDNA'),
			   data.table(score=hsr_scores,type='HSR'),
			   fill=T)
		list(dt=dt,combined_llr=combined_llr,hsr_dats=hsr_dats,ecdna_dats=ecdna_dats,area0=area0,bin_widths=bin_widths)
	}
	if (no_hic){ #one set of null scores per depth; sims scaled from depth[1]
		dt = rbindlist(lapply(depths,function(d){
			f = d/depth
			sc = score_null(lapply(hsr_sims,function(s){s*f}),lapply(ecdna_sims,function(s){s*f}),lapply(random_sims,function(s){s*f}))
			cbind(depth=d,sc$dt)
		}))
		sep_dt = dt[,c(list(sep_pval=suppressWarnings(ks.test(score[type=='HSR'],score[type=='ecDNA'])$p.value)),
			       sep_stats(score[type=='ecDNA'],score[type=='HSR'])),by=depth]
		diag$n_hsr = NULL #reported once, below
		out = c(list(call='None',reason='No Hi-C',hypotheses=hypotheses,n_ec=n_ec,n_hsr=n_hsr,leave_out=leave_out_mode != 'none',leave_out_mode=leave_out_mode,
			   max_ecdna_circ=max(ecdna_circ,na.rm=TRUE),hic.res=hic.res),diag,list(score_dt=dt,sep_dt=sep_dt))
		if (return_all){
			out$hyp = list(gg=gg,ft=ft,ft_context=ft_context,hic.res=hic.res,depth=depth,ecdna_solns=ecdna_solns,hsr_solns=hsr_solns,
				       ecdna_sims=ecdna_sims,hsr_sims=hsr_sims,
				       ec_group=ec_group,pair_outside_match=pair_outside_match,pair_outside_frac_bp=pair_outside_frac_bp,hsr_cov=hsr_cov)
		}
		return(out)
	}
	sc = score_null(hsr_sims,ecdna_sims,random_sims)
	dt = sc$dt; combined_llr = sc$combined_llr; hsr_dats = sc$hsr_dats; ecdna_dats = sc$ecdna_dats; area0 = sc$area0; bin_widths = sc$bin_widths
	true_dat = true_hic$dat
	true_llr = combined_llr(true_dat)
	ecdna_nll = unlist(lapply(ecdna_dats,function(d){compdats(true_dat,d,area0=area0,widths=bin_widths)}))
	best_ecdna_ix = which.min(ecdna_nll)
	hsr_nll = unlist(lapply(hsr_dats,function(d){compdats(true_dat,d,area0=area0,widths=bin_widths)}))
	best_hsr_ix = which.min(hsr_nll)
	likrat_gm = compmaps(true_hic,hsr_sims[[best_hsr_ix]]) - compmaps(true_hic,ecdna_sims[[best_ecdna_ix]])
	likrat_max = max(abs(likrat_gm$value))
	sum_ratio = sum(true_hic$value)/sum(hsr_sims[[best_hsr_ix]]$value)
	cmax = 1/5*max(hsr_sims[[best_hsr_ix]]$value) #pick some shared maximum color for Hi-C
	# make plots. First
	ppdf(plot(ggplot(dt,aes(x=type,y=score))+
		  geom_sina()+
		  geom_hline(yintercept=true_llr,linetype='dashed',color='red')
	  ),width=5,height=4,paste0(figure_path,'/true_hic_topology_call'))
	#
	ppdf(plot(c(gg$gtrack(y.field='cn'),
		    true_hic$gtrack(name='True Hi-C',cmap.max=cmax),
		    hsr_sims[[best_hsr_ix]]$gtrack(name='HSR Hi-C',cmap.max=cmax),
		    ecdna_sims[[best_ecdna_ix]]$gtrack(name='ecDNA Hi-C',cmap.max=cmax),
		    likrat_gm$gtrack(name='Loglik-diff EC - HSR',
				     colormap=c('blue','white','red'),
				     clim=c(-likrat_max,likrat_max))
		    ),ft_context),
	     width=4,height=20,paste0(figure_path,'/true_hic_vs_training'))
	#
	ppdf(plot(c(hsr_solns[[1]]$gtrack(name='HSR'),
		    ecdna_solns[[1]]$gtrack(name='ecDNA'),
		    random_walks[[1]]$gtrack(name='Random'),
		    gg$gtrack(name='gGraph',y.field='cn')),
		  ft_context),
	     width=5,height=20,paste0(figure_path,'/training_examples'))
	#
	hsr_max = max(dt[type=='HSR']$score)
	ec_min = min(dt[type=='ecDNA']$score)
	# call thresholds: quantiles of the simulated scores, not their extremes (user, 2026-11-02). With min/max one outlying
	# ecDNA sample below the observed score flipped COLO320-HSR to ecDNA in several seeds although the observed score sat
	# inside or left of the HSR distribution
	hsr_hi = unname(quantile(dt[type=='HSR']$score,1-call_q))
	ec_lo = unname(quantile(dt[type=='ecDNA']$score,call_q)) #as in ecfinder_call()
	ppdf(plot(ggplot(dt[type=='Random'],aes(x=score,y=circle))+geom_point(size=1) + 
	  geom_vline(xintercept=hsr_hi,color='red',linetype='dashed') + 
	  geom_vline(xintercept=ec_lo,color='green',linetype='dashed') + 
	  labs(x='ecDNA - HSR score',y='amplicon circ. fraction')),
	     width=5,height=4,paste0(figure_path,'/circular_dependence'))
	#
	dt_r = dt[type=='Random']
	circ_cor = cor(dt_r$score,dt_r$circle)
	ec_hsr_sep = ks.test(dt[type=='HSR']$score,dt[type=='ecDNA']$score)$p.value
	ss = sep_stats(dt[type=='ecDNA']$score,dt[type=='HSR']$score)
	separated = ec_hsr_sep < 1e-2 & ss$sep_D > sep_D_min
	call = ecfinder_call(true_llr,dt[type=='ecDNA']$score,dt[type=='HSR']$score,separated,call_q)
	diag$n_hsr = NULL
	out = c(list(call=call,sep_pval = ec_hsr_sep,sep_D=ss$sep_D,sep_D_min=sep_D_min,sep_margin=ss$sep_margin,hypotheses=hypotheses,call_q=call_q,ec_lo=ec_lo,hsr_hi=hsr_hi,max_ecdna_circ=max(ecdna_circ,na.rm=TRUE),circ_cor = circ_cor,ec_score_true = true_llr,score_dt=dt,sum_ratio=sum_ratio,depth_ratio=ratio,n_ec=n_ec,n_hsr=n_hsr,leave_out=leave_out_mode != 'none',leave_out_mode=leave_out_mode),diag)
	if (return_all){
		out$hyp = list(gg=gg,ft=ft,ft_context=ft_context,hic.res=hic.res,true_hic=true_hic,
			       ecdna_solns=ecdna_solns,hsr_solns=hsr_solns,ecdna_sims=ecdna_sims,hsr_sims=hsr_sims,
			       ecdna_nll=ecdna_nll,hsr_nll=hsr_nll,best_ecdna_ix=best_ecdna_ix,best_hsr_ix=best_hsr_ix,likrat_gm=likrat_gm,
			       ec_group=ec_group,pair_outside_match=pair_outside_match,pair_outside_frac_bp=pair_outside_frac_bp,hsr_cov=hsr_cov)
	}
	return(out)
}


# Hi-C call from the observed score and the simulated ecDNA / HSR scores (separated: the KS gate passed): 'ecDNA' above
# the 1 - call_q quantile of the HSR scores, 'CI' below the call_q quantile of the ecDNA scores, else 'None'
ecfinder_call = function(obs,ec,hsr,separated,call_q=0.05){
	if (!separated){return('None')}
	if (obs > quantile(hsr,1-call_q)){'ecDNA'} else if (obs < quantile(ec,call_q)){'CI'} else{'None'}
}

# separation of the simulated ecDNA and HSR score distributions. sep_D: Kolmogorov-Smirnov D, the largest gap between
# the two ECDFs (1 = no overlap; the best single score threshold misclassifies (1-D)/2 of samples at equal priors).
# sep_margin: (ecDNA 5th percentile - HSR 95th percentile) / mean IQR, > 0 when the bulks are apart (D saturates at 1)
sep_stats = function(ec,hsr){
	list(sep_D=unname(suppressWarnings(ks.test(ec,hsr)$statistic)),
	     sep_margin=unname((quantile(ec,0.05) - quantile(hsr,0.95))/((IQR(ec) + IQR(hsr))/2)))
}

find_circles = function(gg,ft,N,context.padding = 1e3,verbose=F,mc.cores=1){
	if (is(ft,'character')){ft = streduce(graph_seqlengths(parse.gr(ft), gg))}
	amplicon_context_nodes = (gg$nodes$gr %&% streduce(ft + context.padding))$node.id %>% unique
	freeze.nodes = gg$nodes$dt[!(node.id %in% amplicon_context_nodes)]$node.id #we don't want to freeze nodes directly adjacent to the amplicon
	if(verbose){message('Sampling circular walks from graph')}
	circles = get_unique_walks(gg,N,mode='circular',frozen.nodes = freeze.nodes,mc.cores=mc.cores)
	if (length(circles)==0){return(circles)}
	circles.ft = tryCatch({
		circles %&% ft
	}, error=function(e){
		ft.nodes = (gg$nodes$gr %&% ft)$node.id %>% unique
		keep = unique(circles$nodesdt[abs(snode.id) %in% ft.nodes,walk.id])
		circles[keep]
	})
	return(circles.ft)
}

score_circles = function(circles,gg,ft){
	if (is(ft,'character')){ft = streduce(parse.gr(ft))}
	if (length(circles)==0){
		return(data.table(walk.id=integer(),ampfrac=numeric(),score=numeric(),max.walk.cn=integer()))
	}
	amplicon_nodes = (gg$nodes$gr %&% ft)$node.id %>% unique
	amplicon_weight = sum(gg$nodes$dt[amplicon_nodes]$cn * gg$nodes$dt[amplicon_nodes]$width)
	if (!is.finite(amplicon_weight) || amplicon_weight <= 0){
		stop('Amplicon has zero graph copy-weight; cannot score circles')
	}
	node.cn = gg$nodes$dt$cn
	names(node.cn) = gg$nodes$dt$node.id
	node.width = gg$nodes$dt$width
	names(node.width) = gg$nodes$dt$node.id
	edge.cn = gg$edges$dt$cn
	names(edge.cn) = gg$edges$dt$edge.id
	walknodes = circles$nodesdt[,.(walk.id,node.id=abs(snode.id),walk.iid)]
	walknodes[,cn:=.N,by=c('node.id','walk.id')]
	walknodes[,graph.cn:=node.cn[as.character(node.id)]]
	walknodes[,width:=node.width[as.character(node.id)]]
	walknodes[is.na(graph.cn),graph.cn:=0]
	walknodes[,maxN:=graph.cn%/%cn]
	walknodes[,max.walk.cn.nodes:=min(maxN),by=walk.id]
	walknodes[,walk.cn:=max.walk.cn.nodes*cn]
	walkedges = circles$edgesdt[,.(walk.id,edge.id=abs(sedge.id),walk.iid)]
	walkedges[,cn:=.N,by=c('edge.id','walk.id')]
	walkedges[,graph.cn:=edge.cn[as.character(edge.id)]]
	walkedges[is.na(graph.cn),graph.cn:=0]
	walkedges[,maxN:=graph.cn%/%cn]
	walkedge_scores = walkedges[,.(max.walk.cn.edges=min(maxN)),by=walk.id]
	walknodes = merge.data.table(walknodes,walkedge_scores,by='walk.id',all.x=T)
	walknodes[,max.walk.cn:=min(max.walk.cn.nodes,max.walk.cn.edges),by='walk.id']
	walknodes[,ampfrac := max.walk.cn*sum(width*(node.id %in% amplicon_nodes)) / amplicon_weight,by=walk.id]
	walkdt = unique(walknodes[,.(walk.id,ampfrac,score=max.walk.cn*ampfrac,max.walk.cn)])
	walkdt[,prob:=score/sum(score)]
	return(walkdt[order(score,decreasing=T)])
}

# greedy circle peel: repeatedly subtract a circle (the best-scoring one if deterministic, else sampled by score) at its
# max copy number from a working copy of gg. Returns the peeled circles (snode.id lists, one per copy) and what is left
# of the graph (loosefix-ed, so it can be decomposed with sample.gwalks)
greedypeel = function(gg,circles,ft,deterministic=TRUE){
	peeled = list()
	mycirc = circles$copy
	# not gg$copy: for graphs read from rds the copy shares gg's edge table, so subtract_walk's
	# edge $mark would strip copy number from gg itself (and from everything sampled from it later)
	gg.cur = gG(nodes = gg$nodes$gr, edges = copy(gg$edges$dt))
	scores = score_circles(mycirc,gg.cur,ft)[max.walk.cn>0 & score>0]
	while (nrow(scores)){
		if (deterministic | (nrow(scores)==1)){
			best = scores[1]
		}else{
			best = scores[sample(1:nrow(scores),1,prob=scores$prob)]
		}
		walk = mycirc[best$walk.id]
		mycirc = mycirc[!(walk.id == best$walk.id)]
		peeled = c(peeled,rep(walk$snode.id,best$max.walk.cn))
		gg.cur = subtract_walk(gg.cur,walk,best$max.walk.cn)
		scores = score_circles(mycirc,gg.cur,ft)[max.walk.cn>0 & score>0]
	}
	list(circles=peeled,rest=loosefix(gg.cur))
}

# sub-graph of gw's graph (same nodes and edges, so ids match and circles, walks on that graph, can be scored and
# subtracted on it directly) with node and edge CN = traversal counts of the walks of gw that touch the amplicon (0
# elsewhere). Edges are copied (copy(gg$edges$dt)), so no $mark reaches the input graph's tables.
touch_subgraph = function(gw,amp){
	gg = gw$graph
	touch = vapply(gw$snode.id,function(s){any(abs(s) %in% amp)},logical(1))
	sub = gG(nodes = gg$nodes$gr, edges = copy(gg$edges$dt))
	sub$nodes$mark(cn = tabulate(abs(unlist(gw$snode.id[touch])),length(gg$nodes)))
	edge.cn = table(abs(unlist(gw$sedge.id[touch])))[as.character(sub$edges$dt$edge.id)]
	sub$edges$mark(cn = ifelse(is.na(edge.cn),0L,as.integer(edge.cn)))
	sub
}

# ecDNA hypotheses paired with the HSR hypothesis gw: all walks touching the amplicon are boiled together (greedypeel on
# touch_subgraph(), first deterministic, then sampled; n_boil times), the rest is decomposed again with every node
# outside the amplicon + 1 kb frozen (as squeeze() does), and the other walks of gw are kept as is. Walks are on
# gw$graph, so amp_circular_fraction() gets the parent graph's amplicon CN. Hypotheses are deduplicated by hash, keeping
# the first. Empty list if there is no circle or no boil peels one.
boil_paired = function(gw,circles,ft,n_boil=3){
	gg = gw$graph
	amp = unique((gg$nodes$gr %&% ft)$node.id)
	touch = vapply(gw$snode.id,function(s){any(abs(s) %in% amp)},logical(1))
	if (!any(touch) || !length(circles)){return(list())}
	sub = touch_subgraph(gw,amp)
	frozen = sub$nodes$dt[!(node.id %in% (gg$nodes$gr %&% streduce(ft + 1e3))$node.id)]$node.id
	oth_snode = gw$snode.id[!touch]
	oth_circ = gw$circular[!touch]
	boils = lapply(seq_len(n_boil),function(i){
		p = greedypeel(sub,circles,ft,deterministic = i==1)
		if (!length(p$circles)){return(NULL)}
		rest = sample.gwalks(p$rest,1,frozen.nodes=frozen,verbose=F)[[1]]
		gW(graph=gg,snode.id=c(oth_snode,rest$snode.id,p$circles),
		   circular=c(oth_circ,rest$circular,rep(TRUE,length(p$circles))))
	})
	boils = Filter(Negate(is.null),boils)
	boils[!duplicated(unlist(lapply(boils,'[[','hash')))] #keeps the deterministic boil
}

# maximal runs of one walk's signed node ids outside the amplicon, as list(key, bp). A run's key carries its boundary
# types ('E' walk end, 'A' amplicon) and is orientation-normalised: the run and its reverse complement (reversed, ids
# negated, boundaries swapped) give the smaller of the two strings. ends: boundary type at the walk's two ends.
walk_runs = function(s,in_amp,node_width,ends='E'){
	n = length(s)
	if (!n){return(list(key=character(0),bp=numeric(0)))}
	grp = cumsum(c(TRUE,in_amp[-1] != in_amp[-n]))
	runs = split(seq_len(n),grp)
	runs = runs[vapply(runs,function(ix){!in_amp[ix[1]]},logical(1))]
	key = vapply(runs,function(ix){
		l = if (ix[1]==1) ends else 'A'
		r = if (ix[length(ix)]==n) ends else 'A'
		k1 = paste0(l,'|',paste(s[ix],collapse=','),'|',r)
		k2 = paste0(r,'|',paste(rev(-s[ix]),collapse=','),'|',l)
		sort(c(k1,k2),method='radix')[1]
	},character(1))
	list(key=unname(key),bp=unname(vapply(runs,function(ix){sum(node_width[abs(s[ix])])},numeric(1))))
}

# the wiring of a hypothesis outside the amplicon: the multiset of maximal runs of non-amplicon snode.ids in its walks
# (between amplicon visits or walk ends), as a sorted character vector of orientation-normalised keys ('E|1,2,3|A'; a
# circular walk without amplicon node is one canonical 'C|...' key over rotations and orientations), with attribute 'bp'
# giving the bp of each run (sum of node widths on gw$graph). Two hypotheses that differ only inside the amplicon have
# equal signatures.
outside_signature = function(gw,amp_nodes){
	node_width = as.numeric(width(gw$graph$nodes$gr))
	keys = character(0); bps = numeric(0)
	for (w in seq_along(gw$snode.id)){
		s = gw$snode.id[[w]]
		in_amp = abs(s) %in% amp_nodes
		if (gw$circular[w]){
			if (any(in_amp)){
				i = which(in_amp)[1]
				s = c(s[i:length(s)],s[seq_len(i-1)]) #rotate to start at an amplicon node, so no run wraps around
				r = walk_runs(s,abs(s) %in% amp_nodes,node_width,ends='A')
			} else{
				n = length(s)
				rots = function(v){vapply(seq_len(n),function(i){paste(v[c(i:n,seq_len(i-1))],collapse=',')},character(1))}
				cand = c(rots(s),rots(rev(-s)))
				r = list(key=paste0('C|',sort(cand,method='radix')[1]),bp=sum(node_width[abs(s)]))
			}
		} else{r = walk_runs(s,in_amp,node_width,ends='E')}
		keys = c(keys,r$key); bps = c(bps,r$bp)
	}
	o = order(keys,method='radix')
	structure(keys[o],bp=bps[o])
}

# compare two outside_signature()s: match = identical multisets; frac_bp = bp in runs present in different numbers,
# sum(|count_a - count_b| * bp) / (total bp_a + total bp_b), 0 if both are empty
outside_mismatch = function(sig_a,sig_b){
	bp = c(attr(sig_a,'bp'),attr(sig_b,'bp'))
	names(bp) = c(sig_a,sig_b)
	bp = bp[!duplicated(names(bp))]
	ca = table(factor(as.character(sig_a),levels=names(bp)))
	cb = table(factor(as.character(sig_b),levels=names(bp)))
	tot = sum(attr(sig_a,'bp')) + sum(attr(sig_b,'bp'))
	list(match = identical(sort(as.character(sig_a),method='radix'),sort(as.character(sig_b),method='radix')),
	     frac_bp = if (tot > 0) sum(abs(as.vector(ca) - as.vector(cb)) * bp) / tot else 0)
}

# copy-weighted amplicon fraction (amp_walk_fraction) carried by the linear walks that visit a centromere node
amp_centric_frac = function(gw,ft,cen_nodes){
	has_cen = vapply(gw$snode.id,function(s){any(abs(s) %in% cen_nodes)},logical(1))
	sum(amp_walk_fraction(gw,ft)[!gw$circular & has_cen])
}

# copy-weighted amplicon fraction carried by the linear walks whose two ends are both canonical-chromosome telomeres
# (not loose ends). The first end is the left side of node s1 if s1 > 0, else the right side of |s1|; the last end is
# the right side of the last node if > 0, else the left side of its |id|. A left side is a telomere if the node starts at
# the lowest start of its chromosome among the graph's nodes (normally 1), a right side if it ends at the highest end.
amp_telo_frac = function(gw,ft){
	nd = gw$graph$nodes$gr
	chr = as.character(seqnames(nd))
	canon = grepl('^(chr)?([0-9]+|X|Y)$',chr)
	left_telo = canon & start(nd) == ave(start(nd),chr,FUN=min)
	right_telo = canon & end(nd) == ave(end(nd),chr,FUN=max)
	telo = function(id,left){if (left) left_telo[abs(id)] else right_telo[abs(id)]}
	ok = vapply(seq_along(gw$snode.id),function(w){
		s = gw$snode.id[[w]]
		if (gw$circular[w] || !length(s)){return(FALSE)}
		a = s[1]; z = s[length(s)]
		telo(a,a > 0) && telo(z,z < 0)
	},logical(1))
	sum(amp_walk_fraction(gw,ft)[ok])
}

# leave-out mode of the scoring: 'group' if the ecDNA hypotheses come from > 1 HSR hypothesis (a noisy sample is scored
# without its whole group), 'solution' if each side has > 1 hypothesis, else 'none'
pick_leave_out_mode = function(ec_group,n_ec,n_hsr){
	if (length(unique(ec_group)) > 1){'group'} else if (n_ec > 1 & n_hsr > 1){'solution'} else{'none'}
}

# hypotheses to leave out when scoring a noisy sample generated by hypothesis k of src_type ('hsr' or 'ec'): list(hsr, ec)
# of integer indices (integer(0) = leave nothing out). ec_group: HSR index of each ecDNA hypothesis.
leave_out_ix = function(src_type,k,ec_group,mode){
	none = list(hsr=integer(0),ec=integer(0))
	if (mode == 'none'){return(none)}
	if (mode == 'solution'){
		if (src_type == 'hsr'){none$hsr = k} else{none$ec = k}
		return(none)
	}
	g = if (src_type == 'hsr') k else ec_group[k]
	list(hsr=g,ec=which(ec_group == g))
}

subtract_walk = function(gg,walk,peel_cn=1){
	if (length(walk)!=1){stop('subtract_walk expects a single-walk gWalk')}
	if (length(peel_cn)!=1 || is.na(peel_cn) || peel_cn < 0){stop('peel_cn must be a non-negative scalar')}
	if (peel_cn==0){return(gg)}
	edge.cn = gg$edges$dt$cn
	names(edge.cn) = gg$edges$dt$edge.id
	node.cn = gg$nodes$dt$cn
	names(node.cn) = gg$nodes$dt$node.id
	edges_sub = walk$edgesdt[,.(edge.id=abs(sedge.id))][,.(sub_cn=.N*peel_cn),by=edge.id]
	edges_sub[,graph.cn:=edge.cn[as.character(edge.id)]]
	edges_sub[,remaining_cn:=graph.cn-sub_cn]
	nodes_sub = walk$nodesdt[,.(node.id=abs(snode.id))][,.(sub_cn=.N*peel_cn),by=node.id]
	nodes_sub[,graph.cn:=node.cn[as.character(node.id)]]
	nodes_sub[,remaining_cn:=graph.cn-sub_cn]
	if (any(edges_sub$remaining_cn < 0) || any(nodes_sub$remaining_cn < 0)){
		stop('subtract_walk would produce negative copy number')
	}
	gg$nodes[nodes_sub$node.id]$mark(cn = nodes_sub$remaining_cn)
	gg$edges[edges_sub$edge.id]$mark(cn = edges_sub$remaining_cn)
	#gg = loosefix(gg)
	return(gg)
}

# ecDNA karyotypes boiled from the whole graph (greedypeel, then the rest decomposed); fastKar_ecfinder boils each HSR
# hypothesis instead (boil_paired)
boil_greedy = function(gg,ft,N,n_solns_try=1,mc.cores=1){
	circles = find_circles(gg,ft,N=N,verbose=FALSE,mc.cores=mc.cores)
	if(!length(circles)){return(sample.gwalks(gg,n_solns_try,verbose=F))}
	boil = function(deterministic){
		p = greedypeel(gg,circles,ft,deterministic)
		rest.gw = sample.gwalks(p$rest,1,verbose=F)[[1]]
		gW(graph=gg,snode.id=c(rest.gw$snode.id,p$circles),circular=c(rest.gw$circular,rep(T,length(p$circles))))
	}
	output = boil(deterministic=T)
	if (n_solns_try>1){
		output_rand = mclapply(1:n_solns_try,function(x){boil(deterministic=F)},mc.cores=mc.cores)
		output = c(list(output),output_rand)
		return(output[!duplicated(unlist(lapply(output,'[[','hash')))])
	}else{
		return(output)
	}
}

#calculate entropy of the given nodes as distributed in the given gwalk object
amp.entropy = function(gw,amp.nodes){
	amp_nodesdt = gw$nodesdt[,.(walk.id,node.id=abs(snode.id))][node.id %in% amp.nodes]
	amp_nodesdt[,walk.amp := .N,by=walk.id]
	amp_nodesdt[,ampfrac:=walk.amp/nrow(.SD)]
	walk.dist = unique(amp_nodesdt[,.(walk.id,ampfrac)])
	-sum(walk.dist$ampfrac*log(walk.dist$ampfrac))
}

#amp.entropy() of a raw karyotype as sample.gwalks(return.gw=F) gives it (list(snode.id, circular)), without building a gWalk
amp.entropy.raw = function(kar,amp.nodes){
	nodes = abs(unlist(kar$snode.id))
	walk = rep(seq_along(kar$snode.id),lengths(kar$snode.id))
	walk.amp = tabulate(walk[nodes %in% amp.nodes],length(kar$snode.id))
	ampfrac = walk.amp[walk.amp>0]/sum(walk.amp)
	-sum(ampfrac*log(ampfrac))
}

# HSR hypotheses: the k_return sampled karyotypes with the lowest amplicon entropy (amplicon on few walks), non-amplicon
# nodes frozen to the default wiring; circular walks touching the amplicon are embedded into linear ones (embedloops).
# keep.acen = TRUE (default since 2026-11-02, user): all karyotypes eligible. keep.acen = FALSE: only karyotypes whose amplicon-heaviest walk (copy-weighted amplicon bp)
# carries a centromere are eligible, since an acentric integrated array would be lost at division (on COLO320-HSR the
# best-fit unrestricted HSR put 68% of the amplicon on an acentric palindrome). centromeres: GRanges, default
# fastKar_centromeres(gg). No eligible karyotype: warning and an empty list.
squeeze = function(gg,ft,N,k_return=1,verbose=F,mc.cores=1,keep.acen=TRUE,centromeres=NULL){
	amplicon_nodes = (gg$nodes$gr %&% ft)$node.id %>% unique 
	amplicon_context_nodes = (gg$nodes$gr %&% streduce(ft + 1e3))$node.id %>% unique 
	freeze.nodes = gg$nodes$dt[!(node.id %in% amplicon_context_nodes)]$node.id #we don't want to freeze nodes directly adjacent to the amplicon, hence ft + 1e3
	if (verbose){message('Sampling walks from graph')}
	# raw karyotypes: gWalks for all N cost ~10x the memory and ~6x the time (NCIH-2170, N = 1000), and only k_return are kept
	kars = sample.gwalks(gg,N,frozen.nodes = freeze.nodes,verbose=F,mc.cores=mc.cores,return.gw=F)
	if (!keep.acen){
		if (is.null(centromeres)){centromeres = fastKar_centromeres(gg)}
		cen_nodes = unique((gg$nodes$gr %&% centromeres)$node.id)
		node_width = as.numeric(width(gg$nodes$gr))
		centric = vapply(kars,function(k){amp_walk_centric(k,amplicon_nodes,cen_nodes,node_width)},logical(1))
		if (verbose){message(sum(centric),' of ',length(kars),' karyotypes have a centric amplicon walk')}
		if (!any(centric)){
			warning('squeeze: none of ',length(kars),' sampled karyotypes has its amplicon on a centric walk')
			return(list())
		}
		kars = kars[centric]
	}
	if (verbose){message('Scoring walks')}
	entropies = unlist(lapply(kars,function(k){amp.entropy.raw(k,amplicon_nodes)}))
	top_walks = lapply(kars[head(order(entropies),k_return)],function(k){gW(graph=gg,snode.id=k$snode.id,circular=k$circular)})
	rm(kars)
	if(verbose){message(paste0('Genrating top-',k_return,' solutions'))}
	mclapply(top_walks,function(gw){
			 if (sum((gw %&% ft)$circular)>0){
				 tryCatch({embedloops(gw)
				 }, error = function(msg){
				 return(gw)})
			 }else{gw}
			 },mc.cores=mc.cores)
}

# TRUE if the walk of a raw karyotype (list(snode.id, circular), as sample.gwalks(return.gw = FALSE) gives it) that carries
# the most amplicon bp (copy-weighted: a node traversed twice counts twice) visits a centromere node. Ties: the first walk.
amp_walk_centric = function(kar,amp_nodes,cen_nodes,node_width){
	nodes = abs(unlist(kar$snode.id))
	walk = rep(seq_along(kar$snode.id),base::lengths(kar$snode.id))
	in_amp = nodes %in% amp_nodes
	if (!any(in_amp)){return(FALSE)}
	amp_bp = vapply(split(node_width[nodes[in_amp]],factor(walk[in_amp],levels=seq_along(kar$snode.id))),sum,numeric(1))
	any(nodes[walk == which.max(amp_bp)] %in% cen_nodes)
}

# Centromeres (UCSC cytoband acen bands, inst/extdata/centromeres_<build>.tsv) as GRanges named like gg's chromosomes
# ('chr8' or '8'). build: 'hg19' or 'hg38'; NULL infers it from gg's chr1 length (graph_seqlengths(); JaBbA rounds
# lengths to kb, so +- 1e4), and errors if chr1 matches neither.
fastKar_centromeres = function(gg,build=NULL){
	sl = GenomeInfoDb::seqlengths(graph_seqlengths(gg$nodes$gr,gg))
	chr1 = sl[names(sl) %in% c('1','chr1')]
	if (is.null(build)){
		ref = c(hg19 = 249250621, hg38 = 248956422)
		hit = names(ref)[length(chr1) == 1 & !is.na(chr1) & abs(ref - chr1[1]) <= 1e4]
		if (length(hit) != 1){stop('fastKar_centromeres: cannot tell hg19 from hg38 by chr1 length (', paste(chr1, collapse = ','), '); pass build or centromeres')}
		build = hit
	}
	f = system.file('extdata',paste0('centromeres_',build,'.tsv'),package='fastKar')
	if (!nzchar(f)){stop('fastKar_centromeres: no centromere table for build ', build)}
	d = data.table::fread(f,skip='seqnames')
	if (any(grepl('^chr',names(sl)))){d[, seqnames := paste0('chr',seqnames)]}
	GRanges(d$seqnames,IRanges(d$start,d$end))
}

amp_circular_fraction = function(gw,ft){
	gw = gw %&% ft
	if (sum(gw$circular)==0){return(0)}
	gg = gw$graph
	event_nodes = unique((gg$gr %&% ft)$node.id)
	amplicon_size = sum(width(gg$nodes$gr[event_nodes])*gg$nodes$dt[event_nodes]$cn)
	circnodes = abs(unlist(gw$snode.id[which(gw$circular==T)]))
	if (!is.null(circnodes)){
		circnodes = circnodes[circnodes %in% event_nodes]
	}
	if (length(circnodes)){return(sum(width(gg$nodes$gr[circnodes]))/amplicon_size)
	}else{return(0)}
}

# copy-weighted fraction of the amplicon carried by each walk of gw: sum(width x traversals of amplicon nodes in the
# walk) / sum(width x graph CN of amplicon nodes), one value per walk
amp_walk_fraction = function(gw,ft){
	gg = gw$graph
	event_nodes = unique((gg$gr %&% ft)$node.id)
	amplicon_size = sum(width(gg$nodes$gr[event_nodes])*gg$nodes$dt[event_nodes]$cn)
	nodes = abs(unlist(gw$snode.id))
	walk = rep(seq_along(gw$snode.id),lengths(gw$snode.id))
	keep = nodes %in% event_nodes
	w = tapply(width(gg$nodes$gr[nodes[keep]]),factor(walk[keep],levels=seq_along(gw)),sum)
	as.vector(ifelse(is.na(w),0,w))/amplicon_size
}

# HSR coverage: amp_walk_fraction() of the best linear walk (0 if there is none)
amp_linear_coverage = function(gw,ft){
	if (all(gw$circular)){return(0)}
	max(amp_walk_fraction(gw,ft)[!gw$circular])
}

hic_res = function(path) {
	path = normalizePath(path)
	pathtype = tools::file_ext(path)
	if (pathtype=='hic'){
		reses = strawr::readHicBpResolutions(path) %>% sort() %>% signif(., digits = 5)
	}else if (pathtype=='mcool'){
		cooler_details.dt = rhdf5::h5ls(path) %>% as.data.table()
    		suppressWarnings(cooler_details.dt[,name := as.integer(name)])
    		reses= cooler_details.dt[!is.na(name),]$name %>% sort()
	}
	return(reses)
}

# TRUE if an mcool names its chromosomes 'chr1', ... (.hic files are assumed not to)
hic_is_chr = function(path) {
	if (tools::file_ext(path) != 'mcool') return(FALSE)
	chroms = rhdf5::h5read(path, paste0('/resolutions/', hic_res(path)[1], '/chroms/name'))
	any(grepl('^chr', chroms))
}

## Hi-C scoring window around a focus locus (e.g. a FISH probe gene) inside a large amplicon, for fastKar_ecfinder(score_region=).
## core: the contiguous run of nodes at CN > amp_mult x ploidy that overlaps focus (adjacent amplified nodes merged), +- pad;
## plus +- partner_pad around the far breakpoints of ALT junctions with an end in the core and CN > amp_mult x ploidy.
## Only canonical chromosomes (1-22, X, Y, with or without 'chr') are kept: unplaced/alt contigs (e.g. 16_KI270728v1_random,
## a high-CN junction partner in COLO320) are not scored.
## Returns the reduced window (GRanges), or NULL if focus is not on an amplified node.
focus_score_window = function(gg, focus, amp_mult = 2, pad = 3e5, partner_pad = 1e6){
	if (is(focus,'character')){focus = parse.gr(focus)}
	nd = gg$nodes$dt
	ploidy = nd[!is.na(cn), sum(cn*as.numeric(width))/sum(as.numeric(width))]
	amp = gg$nodes$gr[which(nd$cn > amp_mult*ploidy)]
	if (!length(amp)) return(NULL)
	segs = reduce(amp[, c()]) #adjacent amplified nodes become one segment
	core = segs[segs %^% focus]
	if (!length(core)) return(NULL)
	core.nodes = (gg$nodes$gr %&% core)$node.id
	alt = gg$edges[type == 'ALT' & cn > amp_mult*ploidy]
	alt = alt[alt$dt$n1 %in% core.nodes | alt$dt$n2 %in% core.nodes]
	win = core + pad
	if (length(alt)){
		bps = unlist(alt$junctions$grl)
		bps = bps[!(bps %^% (core + pad))]
		if (length(bps)){win = c(win[, c()], bps[, c()] + partner_pad)}
	}
	win = win[grepl('^(chr)?([0-9]+|X|Y)$', as.character(seqnames(win)))]
	streduce(win[, c()])
}
