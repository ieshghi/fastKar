fastKar_ecfinder = function(gg,ft,hic.res,true_hic_path,ec_fraction_threshold = 0.6,n_sample=100,n_ec = 10,figure_path = NULL,mc.cores=1,return_all=F,purity=1){
	library(ggforce)
	wholegenome = si2gr(hg_seqlengths(chr=FALSE)) %Q% (seqnames %in% c(1:22,'X','Y'))
	if (is(ft,'character')){ft = streduce(parse.gr(ft))}
	ft_context = streduce(ft + sum(width(ft))/10) #add 10% to the footprint for context
	event_nodes = (gg$nodes$gr %&% ft)$node.id
	message('Sampling and generating solutions')
	# generate HSR and ecDNA solutions
	ecdna_solns = boil_greedy(gg=gg,ft=ft,N=100,n_solns_try=n_ec,mc.cores=mc.cores)
	ecdna_circ = unname(unlist(lapply(
		ecdna_solns,function(x){amp_circular_fraction(x,ft)})))
	ecdna_solns = ecdna_solns[ecdna_circ > ec_fraction_threshold]
	if (!length(ecdna_solns)){return(list(call='CI', reason='No good ecDNA solutions'))}
	hsr_solns = squeeze(gg=gg,ft=ft,N=1000,k_return=1,verbose = F,mc.cores=mc.cores)
	n_ec = length(ecdna_solns)
	n_hsr = length(hsr_solns)
	random_walks = sample.gwalks(gg,n_sample,mc.cores=mc.cores,verbose = F)
	random_circ = unname(unlist(mclapply(
		random_walks,function(x){amp_circular_fraction(x,ft)},mc.cores=mc.cores)))
	# estimate depth of hi-c data and pick resolution to simulate at
	ploidy = sum(gg$gr$cn*width(gg$gr))/sum(width(gg$gr))
	depth = estimate.depthratio(true_hic_path,ploidy=ploidy,purity=purity)
	hictype = tools::file_ext(true_hic_path)
	if (hic.res<0){ #if hic.res < 0 this is asking for a number of bins to cover the footprint
		hic.res = sum(width(ft_context))/abs(hic.res)
	}
	# now find available resolution that is closest to requested resolution
	avail_res = hic_res(true_hic_path)
	dists = abs(log10(avail_res) - log10(hic.res))
	hic.res = avail_res[dists==min(dists)]
	# load hic data
	if (hictype=='mcool'){true_hic = cooler(true_hic_path,gr=ft_context,res=hic.res)
	} else{true_hic = straw(true_hic_path,gr=ft_context,res=hic.res)}
	hic.gr = true_hic$gr
	true_hic = rebin_matrix(true_hic,hic.gr) #make sure the matrix matches future matrices (basically just adds in 0s)
	message('Simulating Hi-C for hypotheses and random samples')
	# simulate Hi-C for each of the generated solutions, rebin to the data GRanges
	simfun=function(x){
		rebin_matrix(forward_simulate(x,target_region = ft_context,
					      pix.size=hic.res,depth=depth,purity=purity),hic.gr)
	}
	hsr_sims = mclapply(hsr_solns,simfun,mc.cores=mc.cores)
	ecdna_sims = mclapply(ecdna_solns,simfun,mc.cores=mc.cores)
	random_sims = mclapply(random_walks,simfun,mc.cores=mc.cores)
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
	sum_ratio_ec = mean(unlist(lapply(ec_noisy_sims,function(x){sum(x$value)})))/sum(true_hic$value)
	sum_ratio_hsr = mean(unlist(lapply(hsr_noisy_sims,function(x){sum(x$value)})))/sum(true_hic$value)
	scaling_ratio = mean(c(sum_ratio_ec,sum_ratio_hsr))
	true_hic_scaled = true_hic*scaling_ratio
	#
	message('Scoring karyotypes')
	area0 = median(width(ecdna_sims[[1]]$gr))^2
	# ecDNA score is mean likelihood of ecDNA - likelihood of HSR
	nll = function(sims,data){
		min(unlist(lapply(sims,function(sim){compdats(data,sim$dat,area0=area0)})))}
	combined_llr = function(x){nll(hsr_sims,x)-nll(ecdna_sims,x)} 
	# score all samples and true Hi-C
	random_scores = mclapply(random_noisy_sims,combined_llr,mc.cores=mc.cores) %>% unlist
	ec_scores = mclapply(ec_noisy_sims,combined_llr,mc.cores=mc.cores) %>% unlist
	hsr_scores = mclapply(hsr_noisy_sims,combined_llr,mc.cores=mc.cores) %>% unlist
	true_llr = combined_llr(true_hic_scaled$dat)
	#
	dt = rbind(data.table(score=random_scores,type='Random',circle = random_circ),
		   data.table(score=ec_scores,type='ecDNA'),
		   data.table(score=hsr_scores,type='HSR'),
		   fill=T)
	ecdna_nll = unlist(lapply(ecdna_sims,function(sim){compdats(true_hic_scaled$dat,sim$dat,area0=area0)}))
	best_ecdna_ix = which.min(ecdna_nll)
	hsr_nll = unlist(lapply(hsr_sims,function(sim){compdats(true_hic_scaled$dat,sim$dat,area0=area0)}))
	best_hsr_ix = which.min(hsr_nll)
	likrat_gm = compmaps(true_hic,hsr_sims[[best_hsr_ix]]) - compmaps(true_hic,ecdna_sims[[best_ecdna_ix]])
	likrat_max = max(abs(likrat_gm$value))
	# make plots. First
	ppdf(plot(ggplot(dt,aes(x=type,y=score))+
		  geom_sina()+
		  geom_hline(yintercept=true_llr,linetype='dashed',color='red')
	  ),width=5,height=4,paste0(figure_path,'/true_hic_topology_call'))
	#
	ppdf(plot(c(gg$gtrack(y.field='cn'),
		    true_hic_scaled$gtrack(name='True Hi-C'),
		    hsr_sims[[best_hsr_ix]]$gtrack(name='HSR Hi-C'),
		    ecdna_sims[[best_ecdna_ix]]$gtrack(name='ecDNA Hi-C'),
		    likrat_gm$gtrack(name='Loglik-diff EC - HSR',
				     colormap=c('blue','white','red'),
				     clim=c(-likrat_max,likrat_max))
		    ),ft_context),
	     width=5,height=20,paste0(figure_path,'/true_hic_vs_training'))
	#
	ppdf(plot(c(hsr_solns[[1]]$gtrack(name='HSR'),
		    ecdna_solns[[1]]$gtrack(name='ecDNA'),
		    random_walks[[1]]$gtrack(name='Random'),
		    gg$gtrack(name='gGraph',y.field='cn')),
		  ft_context),
	     width=7,height=25,paste0(figure_path,'/training_examples'))
	#
	hsr_max = max(dt[type=='HSR']$score)
	ec_min = min(dt[type=='ecDNA']$score)
	ppdf(plot(ggplot(dt[type=='Random'],aes(x=score,y=circle))+geom_point(size=1) + 
	  geom_vline(xintercept=hsr_max,color='red',linetype='dashed') + 
	  geom_vline(xintercept=ec_min,color='green',linetype='dashed') + 
	  labs(x='ecDNA - HSR score',y='amplicon circ. fraction')),
	     width=5,height=4,paste0(figure_path,'/circular_dependence'))
	#
	dt_r = dt[type=='Random']
	circ_cor = cor(dt_r$score,dt_r$circle)
	ec_hsr_sep = ks.test(dt[type=='HSR']$score,dt[type=='ecDNA']$score)$p.value
	if (ec_hsr_sep < 1e-2 & true_llr > ec_min){
		return(list(call='ecDNA',sep_pval = ec_hsr_sep,circ_cor = circ_cor,ec_score_true = true_llr,score_dt=dt,scaling_ratio=scaling_ratio))
	}else if (ec_hsr_sep < 1e-2 & true_llr < hsr_max){
		return(list(call='CI',sep_pval = ec_hsr_sep,circ_cor = circ_cor,ec_score_true = true_llr,score_dt=dt,scaling_ratio=scaling_ratio))
	}else{
		return(list(call='None',sep_pval = ec_hsr_sep,circ_cor = circ_cor,ec_score_true = true_llr,score_dt=dt,scaling_ratio=scaling_ratio))
	}
}


find_circles = function(gg,ft,N,context.padding = 1e3,verbose=F,mc.cores=1){
	if (is(ft,'character')){ft = streduce(parse.gr(ft))}
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

boil_greedy = function(gg,ft,N,n_solns_try=1,mc.cores=1){
	peeledcircles = list()
	circles = find_circles(gg,ft,N=N,verbose=FALSE,mc.cores=mc.cores)
	if(!length(circles)){return(sample.gwalks(gg,n_solns_try,verbose=F))}
	greedypeel = function(gg,circles,deterministic){	
		mycirc = circles$copy
		gg.cur = gg$copy
		scores = score_circles(mycirc,gg.cur,ft)[max.walk.cn>0 & score>0]
		while (nrow(scores)){
			if (deterministic | (nrow(scores)==1)){
				best = scores[1]
			}else{
				best = scores[sample(1:nrow(scores),1,prob=scores$prob)]
			}
			walk = mycirc[best$walk.id]
			mycirc = mycirc[!(walk.id == best$walk.id)]
			peeledcircles = c(peeledcircles,rep(walk$snode.id,best$max.walk.cn))
			gg.cur = subtract_walk(gg.cur,walk,best$max.walk.cn)
			scores = score_circles(mycirc,gg.cur,ft)[max.walk.cn>0 & score>0]
		}
		peeled_gw = gW(graph=gg,snode.id = peeledcircles,circular = rep(T,length(peeledcircles)))
		rest.gw = sample.gwalks(loosefix(gg.cur),1,verbose=F)[[1]]
		return(c(gW(graph=gg,snode.id=rest.gw$snode.id,circular=rest.gw$circular),peeled_gw))
	}
	output = mclapply(1:n_solns_try,function(x){greedypeel(gg,circles,deterministic=n_solns_try==1)},mc.cores=mc.cores)
	return(output[!duplicated(unlist(lapply(output,'[[','hash')))])
}

#calculate entropy of the given nodes as distributed in the given gwalk object
amp.entropy = function(gw,amp.nodes){
	amp_nodesdt = gw$nodesdt[,.(walk.id,node.id=abs(snode.id))][node.id %in% amp.nodes]
	amp_nodesdt[,walk.amp := .N,by=walk.id]
	amp_nodesdt[,ampfrac:=walk.amp/nrow(.SD)]
	walk.dist = unique(amp_nodesdt[,.(walk.id,ampfrac)])
	-sum(walk.dist$ampfrac*log(walk.dist$ampfrac))
}

squeeze = function(gg,ft,N,k_return=1,verbose=F,mc.cores=1){
	amplicon_nodes = (gg$nodes$gr %&% ft)$node.id %>% unique 
	amplicon_context_nodes = (gg$nodes$gr %&% streduce(ft + 1e3))$node.id %>% unique 
	freeze.nodes = gg$nodes$dt[!(node.id %in% amplicon_context_nodes)]$node.id #we don't want to freeze nodes directly adjacent to the amplicon, hence ft + 1e3
	if (verbose){message('Sampling walks from graph')}
	gwl = sample.gwalks(gg,N,frozen.nodes = freeze.nodes,verbose=F,mc.cores=mc.cores)
	if (verbose){message('Scoring walks')}
	entropies = unlist(lapply(gwl,function(gw){amp.entropy(gw,amplicon_nodes)}))
	top_walks = gwl[order(entropies)[1:k_return]]
	if(verbose){message(paste0('Genrating top-',k_return,' solutions'))}
	mclapply(top_walks,function(gw){
			 if (sum((gw %&% ft)$circular)>0){
				 tryCatch({embedloops(gw)
				 }, error = function(msg){
				 return(gw)})
			 }else{gw}
			 },mc.cores=mc.cores)
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
