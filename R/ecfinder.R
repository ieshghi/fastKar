fastKar_ecfinder = function(gg,ft,n_sample=100,n_ec = 10,sim.resolution = 1e4,true_hic_path = NULL,hic.res=NULL,figure_path = NULL,mc.cores=1,return_all=F){
	library(ggforce)
	checked_mclapply = function(X, FUN, ..., mc.cores = 1, step = 'mclapply'){
		out = parallel::mclapply(X, FUN, ..., mc.cores = mc.cores)
		bad = which(vapply(out, inherits, logical(1), what = 'try-error'))
		if (length(bad)){
			msg = paste(vapply(out[bad], function(x) paste(as.character(x), collapse = '\n'), character(1)), collapse = '\n---\n')
			stop(sprintf('%s failed in %d/%d forked workers. First failing indices: %s\n%s', step, length(bad), length(out), paste(head(bad, 10), collapse = ','), msg), call. = FALSE)
		}
		return(out)
	}
	wholegenome = si2gr(hg_seqlengths(chr=FALSE)) %Q% (seqnames %in% c(1:22,'X','Y'))
	if (length(grep('chr',as.character(seqnames(ft))))){wholegenome = gr.chr(wholegenome)}
	gg = loosefix(gg %&% wholegenome)
	ft = ft %&% wholegenome
	#message('Pasting graph loose ends')
	#gg = paste_loose_ends_timed(gg$copy,maxtime=120)
	if (is(ft,'character')){ft = streduce(parse.gr(ft))}
	ft_context = streduce(ft + sum(width(ft))/10) #add 10% to the footprint for context
	event_nodes = (gg$nodes$gr %&% ft)$node.id
	message('Sampling and generating solutions')
	hsr_solns = squeeze(gg=gg,ft=ft,N=1000,k_return=n_ec,verbose = F,mc.cores=mc.cores)
	ecdna_solns = boil(gg=gg,ft=ft,N=100,k_return=n_ec,verbose = F,mc.cores=mc.cores)
	random_walks = sample.gwalks(gg,n_sample,mc.cores=mc.cores,verbose = F)
	#
	resolution = sim.resolution
	if (!is.null(true_hic_path)){
		hictype = tools::file_ext(true_hic_path)
		depth = estimate.depthratio(true_hic_path,mode=hictype)
	}else{
		depth = 10
	}
	#
	message('Simulating Hi-C for hypotheses and random samples')
	simfun = function(x){forward_simulate(x,target_region = ft_context,pix.size=resolution,depth=depth)}
	hsr_sims = checked_mclapply(hsr_solns,simfun,mc.cores=mc.cores,step='forward_simulate HSR')
	ecdna_sims = checked_mclapply(ecdna_solns,simfun,mc.cores=mc.cores,step='forward_simulate ecDNA')
	random_sims = checked_mclapply(random_walks,simfun,mc.cores=mc.cores,step='forward_simulate random')
	if(!is.null(true_hic_path)){
		if (is.null(hic.res)){
			hic.res=resolution
		}
		if (hictype=='mcool'){
			true_hic = cooler(true_hic_path,gr=ft_context,res=hic.res)
		} else{
			true_hic = straw(true_hic_path,gr=ft_context,res=hic.res)
		}
		hic.gr = true_hic$gr
		message('Re-binning to data coordinates')
		hsr_sims = mclapply(hsr_sims,function(s){rebin_matrix(s,hic.gr)},mc.cores=mc.cores)
		ecdna_sims = mclapply(ecdna_sims,function(s){rebin_matrix(s,hic.gr)},mc.cores=mc.cores)
		random_sims = mclapply(random_sims,function(s){rebin_matrix(s,hic.gr)},mc.cores=mc.cores)
	}
	if (length(ecdna_sims) < n_ec | length(hsr_sims) < n_ec){
		n_ec = min(length(ecdna_sims),length(hsr_sims))
	}
	ecdna_train = ecdna_sims[1:n_ec]
	ecdna_test = ecdna_train
	hsr_train = hsr_sims[1:n_ec]
	hsr_test = hsr_train
	message('Simulating noise')
	n_hsr_sample = n_sample %/% length(hsr_test)
	n_ec_sample = n_sample %/% length(ecdna_test)
	samplefun = function(x,n){make_noisydat(x,nsamp=n,theta=2)}
	hsr_testing_sims = unlist(checked_mclapply(hsr_test,function(x){samplefun(x,n_hsr_sample)},mc.cores=mc.cores,step='make_noisydat HSR'),recursive=F)
	ec_testing_sims = unlist(checked_mclapply(ecdna_test,function(x){samplefun(x,n_ec_sample)},mc.cores=mc.cores,step='make_noisydat EC'),recursive=F)
	random_sims_noise = checked_mclapply(random_sims,function(x){make_noisydat(x)[[1]]},mc.cores=mc.cores,step='make_noisydat random')
	#
	message('Calculating circular fraction')
	amplicon_size = sum(width(gg$nodes$gr[event_nodes])*gg$nodes$dt[event_nodes]$cn)
	circfunc = function(gw){
		circnodes = unlist(gw$snode.id[which(gw$circular==T)])
		if (!is.null(circnodes)){
			circnodes = circnodes[circnodes %in% event_nodes]
			circnodes = abs(circnodes)
		}
		if (length(circnodes)){circ_size = sum(width(gg$nodes$gr[circnodes]))
		}else{circ_size = 0}
		return(circ_size / amplicon_size)}
	random_circ = unname(unlist(checked_mclapply(random_walks,circfunc,mc.cores=mc.cores,step='circular fraction random')))
	#
	message('Scoring karyotypes')
	area0 = median(width(ecdna_sims[[1]]$gr))^2
	nll = function(sims,data){sum(unlist(lapply(sims,function(sim){compdats(data,sim$dat,area0=area0)})))}
	combined_llr = function(x){nll(hsr_train,x)-nll(ecdna_train,x)} #likelihood of ecDNA - likelihood of HSR
	#
	random_scores = checked_mclapply(random_sims_noise,combined_llr,mc.cores=mc.cores,step='score random') %>% unlist
	ec_scores = checked_mclapply(ec_testing_sims,combined_llr,mc.cores=mc.cores,step='score ecDNA') %>% unlist
	hsr_scores = checked_mclapply(hsr_testing_sims,combined_llr,mc.cores=mc.cores,step='score HSR') %>% unlist
	#
	dt = rbind(data.table(score=random_scores,type='Random',circle = random_circ),
		   data.table(score=ec_scores,type='ecDNA'),
		   data.table(score=hsr_scores,type='HSR'),fill=T)
	#
	ppdf(plot(c(hsr_solns[[1]]$gtrack(name='HSR'),ecdna_solns[[1]]$gtrack(name='ecDNA'),random_walks[[1]]$gtrack(name='Random'),gg$gtrack(name='gGraph',y.field='cn')),ft_context),width=7,height=25,paste0(figure_path,'/training_examples'))
	ppdf(plot(ggplot(dt,aes(x=type,y=score))+geom_sina()),width=5,height=4,paste0(figure_path,'/testing_sinaplots'))
	hsr_max = max(dt[type=='HSR']$score)
	ec_min = min(dt[type=='ecDNA']$score)
	ppdf(plot(ggplot(dt[type=='Random'],aes(x=score,y=circle))+geom_point(size=1) + 
	  geom_vline(xintercept=hsr_max,color='red',linetype='dashed') + 
	  geom_vline(xintercept=ec_min,color='green',linetype='dashed') + 
	  labs(x='ecDNA - HSR score',y='amplicon circ. fraction')),width=5,height=4,paste0(figure_path,'/circular_dependence'))
	#
	if(!is.null(true_hic_path)){
		llr_gm = function(sims,data_gm){sum(unlist(lapply(sims,function(sim){compmaps(data_gm,sim,ifsum=T)})))}
		combined_llr_gm = function(x){llr_gm(hsr_train,x)-llr_gm(ecdna_train,x)}
		true_llr = combined_llr_gm(true_hic)
		ppdf(plot(ggplot(dt,aes(x=type,y=score))+geom_sina() + geom_hline(yintercept=true_llr,linetype='dashed',color='red')),width=5,height=4,paste0(figure_path,'/true_hic_topology_call'))
		likrat_gm = compmaps(true_hic,hsr_train[[1]]) - compmaps(true_hic,ecdna_train[[1]])
		likrat_max = max(abs(likrat_gm$value))
		ppdf(plot(c(true_hic$gtrack(name='True Hi-C'),hsr_sims[[1]]$gtrack(name='HSR Hi-C'),ecdna_sims[[1]]$gtrack(name='ecDNA Hi-C'),likrat_gm$gtrack(name='Loglik-diff EC - HSR',colormap=c('blue','white','red'),clim=c(-likrat_max,likrat_max))),ft_context),width=5,height=20,paste0(figure_path,'/true_hic_vs_training'))
	}else{true_llr = NA}
	dt_r = dt[type=='Random']
	if (is.numeric(dt_r$score) & is.numeric(dt_r$circle)){
		circ_cor = cor(dt_r$score,dt_r$circle)
		ec_hsr_sep = ks.test(dt[type=='HSR']$score,dt[type=='ecDNA']$score)$p.value
	}else{
		circ_cor = NA
		ec_hsr_sep = NA
	}
	return(list(sep_pval = ec_hsr_sep,circ_cor = circ_cor,ec_score_true = true_llr, ec_min=ec_min,hsr_max=hsr_max))
}

boil = function(gg,ft,N,k_return = 1,verbose=F,mc.cores=mc.cores){
	if (is(ft,'character')){ft = streduce(parse.gr(ft))}
	amplicon_nodes = (gg$nodes$gr %&% ft)$node.id %>% unique 
	amplicon_context_nodes = (gg$nodes$gr %&% streduce(ft + 1e3))$node.id %>% unique 
	freeze.nodes = gg$nodes$dt[!(node.id %in% amplicon_context_nodes)]$node.id #we don't want to freeze nodes directly adjacent to the amplicon, hence ft + 1e3
	if(verbose){message('Sampling walks from graph')}
	walks = get_unique_walks(gg,N,mode='circular',frozen.nodes = freeze.nodes,mc.cores=mc.cores) %&% ft
	#get max CN of walks in the graph
	if(verbose){message('Scoring walks')}
	walknodes = walks$nodesdt[,.(walk.id,node.id=abs(snode.id),walk.iid)]
	walknodes[,cn:=.N,by=c('node.id','walk.id')]
	walknodes = merge.data.table(walknodes,gg$nodes$dt[,.(node.id,graph.cn=cn,width)],by=c('node.id'),all.x=T,cartesian=T)
	walknodes[is.na(graph.cn),graph.cn:=0]
	walknodes[,maxN:=graph.cn%/%cn]
	walknodes[,max.walk.cn.nodes:=min(maxN),by=walk.id]
	walknodes[,walk.cn:=max.walk.cn.nodes*cn]
	walkedges = walks$edgesdt[,.(walk.id,edge.id=abs(sedge.id),walk.iid)]
	walkedges[,cn:=.N,by=c('edge.id','walk.id')]
	walkedges = merge.data.table(walkedges,gg$edges$dt[,.(edge.id,graph.cn=cn)],by=c('edge.id'),all.x=T,allow.cartesian=T)
	walkedges[is.na(graph.cn),graph.cn:=0]
	walkedges[,maxN:=graph.cn%/%cn]
	walkedges[,max.walk.cn.edges:=min(maxN),by=walk.id]
	walknodes = unique(merge.data.table(walknodes,walkedges[,.(walk.id,max.walk.cn.edges)],by='walk.id',all.x=T,allow.cartesian=T))
	walknodes[,max.walk.cn:=min(max.walk.cn.nodes,max.walk.cn.edges),by='walk.id']
	#what fraction of the amplicon does the walk account for
	amplicon_weight = sum(gg$nodes$dt[amplicon_nodes]$cn * gg$nodes$dt[amplicon_nodes]$width)
	walknodes[,ampfrac := max.walk.cn*sum(width*(node.id %in% amplicon_nodes)) / amplicon_weight,by=walk.id]
	walkdt = unique(walknodes[,.(walk.id,ampfrac,score=max.walk.cn*ampfrac,max.walk.cn)])
	#what is the entropy of the walk
	#walknodes[,node.entropy := -max.walk.cn*cn/walk.cn*log(cn/walk.cn)]
	#walknodes[,walk.entropy := sum(node.entropy),by=walk.id]
	#score walks combining both these things.. entropy * amplicon fraction perhaps?
	#NOTE: entropy calculation seems bad so far, so we use walk.cn * ampfrac as a score instead
	ord = order(walkdt$score,decreasing=T)
	sorted_walks = walks[walkdt[ord]$walk.id]
	#make an ecDNA solution from the kth best walk
	if(verbose){message(paste0('Genrating top-',k_return,' solutions'))}
	if (length(sorted_walks) < k_return){k_return = length(sorted_walks)}
	k_solns = mclapply(1:k_return,function(k){
		walk_to_peel = sorted_walks[k]
		peel_cn=walkdt[ord[k]]$max.walk.cn
		#calculate what the remaining node and edge CN in the graph will be
		edges_sub = merge.data.table(unique(walk_to_peel$edgesdt[,.(edge.id=abs(sedge.id))][,.(edge.id,sub_cn = .N*peel_cn),by=edge.id])[,.(edge.id,sub_cn)],gg$edges$dt[,.(edge.id,graph.cn=cn)],by='edge.id')[,.(edge.id,remaining_cn=graph.cn-sub_cn)]
		nodes_sub = merge.data.table(unique(walk_to_peel$nodesdt[,.(node.id=abs(snode.id))][,.(node.id,sub_cn = .N*peel_cn),by=node.id])[,.(node.id,sub_cn)],gg$nodes$dt[,.(node.id,graph.cn=cn)],by='node.id')[,.(node.id,remaining_cn=graph.cn-sub_cn)]
		#instantiate new graph with updated CNs, then sample from it to get the final walks to concatenate to the walks
		gg_out = gg$copy
		gg_out$nodes[nodes_sub$node.id]$mark(cn = nodes_sub$remaining_cn)
		gg_out$edges[edges_sub$edge.id]$mark(cn = edges_sub$remaining_cn)
		gg_out = loosefix(gg_out[,cn>0])
		remaining.walks = sample.gwalks(gg_out,1,verbose=F)[[1]]
		nr = length(remaining.walks)
		combined.gw = gW(graph=gg,snode.id=c(remaining.walks$snode.id,rep(walk_to_peel$snode.id,peel_cn)),circular=c(remaining.walks$circular,rep(T,peel_cn)))
		return(combined.gw)
	},mc.cores=mc.cores)
	return(k_solns)
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
