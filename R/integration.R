fastKar_ecdna = function(gg,mcool,amp,mask=NULL,amp_cn=8,min_amp_bp=5e5,pad=1e6,chrs=c(1:22,'X'),res=1e6,alpha=0.05){
	# ecDNA test (user, 2026-11-10): does the amplicon contact the whole genome more, per copy, than ordinary sequence does, on most chromosomes? One scale (res tiles from the mcool),
	# no fit, no pseudocount, no trimming. amp: the amplicon footprint (GRanges, e.g. the JaBbA event footprint).
	# 1. tiles: the mcool's res bins on chrs, dropping tiles with >= res / 2 masked. amplicon tiles: >= min_amp_bp inside CN >= amp_cn nodes of amp (if none, the one tile with the most such bp).
	#    genome tiles: every node overlapping the tile has 1 <= CN < amp_cn (CN.bp > 0), and the tile is > pad from any CN >= amp_cn node
	# 2. expected contacts between two tiles on different chromosomes: CN.bp x CN.bp (CN.bp = sum over the tile's graph nodes of CN x overlap bp)
	# 3. per genome tile g (trans pairs only): r_g = log2[(OA_g / EA_g) / (OG_g / EG_g)], OA / EA = reads / expected with all amplicon tiles, OG / EG = same with all other genome tiles.
	#    OA = 0 gives r = -Inf (below 0); tiles with EA = 0 (all amplicon tiles on their chromosome) or OG = 0 are dropped
	# 4. per chromosome: median r_g; call: one-sided sign test across chromosomes, n_chr_pos of n_chr medians > 0 ~ Binomial(n_chr, 1/2); ecDNA if p < alpha
	# Returns list(call, p, n_chr_pos, n_chr, med (median r_g), score (descriptive: log2[(sum OA / sum EA) / (sum OG / sum EG)]), n_amp_tiles, tiles, chr)
	nodes = gg$nodes$gr; nd = gg$nodes$dt; nodes = GRanges(gsub('^chr','',as.character(seqnames(nodes))),ranges(nodes)); keep = as.character(seqnames(nodes)) %in% as.character(chrs)
	nodes = nodes[keep]; cn = nd$cn[keep]; hi = reduce(nodes[cn>=amp_cn])
	amp = reduce(GRanges(gsub('^chr','',as.character(seqnames(amp))),ranges(amp)))
	G = paste0('resolutions/',format(res,scientific=FALSE)); on.exit(rhdf5::h5closeAll())
	ch = gsub('^chr','',as.character(rhdf5::h5read(mcool,paste0(G,'/bins/chrom')))); st = as.vector(rhdf5::h5read(mcool,paste0(G,'/bins/start'))); en = as.vector(rhdf5::h5read(mcool,paste0(G,'/bins/end')))
	tl = GRanges(ch,IRanges(st+1,en)); n = length(tl)
	ovbp = function(R){a = numeric(n); if (!length(R)) return(a); o = findOverlaps(tl,R); if (length(o)){s = tapply(as.numeric(width(pintersect(tl[queryHits(o)],R[subjectHits(o)]))),queryHits(o),sum); a[as.integer(names(s))] = s}; a}
	mbp = if (is.null(mask)) numeric(n) else ovbp(reduce(GRanges(gsub('^chr','',as.character(seqnames(mask))),ranges(mask))))
	o = findOverlaps(tl,nodes); w = as.numeric(width(pintersect(tl[queryHits(o)],nodes[subjectHits(o)])))
	cnbp = numeric(n); s = tapply(cn[subjectHits(o)]*w,queryHits(o),sum); cnbp[as.integer(names(s))] = s
	mxcn = numeric(n); s = tapply(cn[subjectHits(o)],queryHits(o),max); mxcn[as.integer(names(s))] = s
	ok = ch %in% as.character(chrs) & mbp<res/2; abp = ovbp(GenomicRanges::intersect(hi,amp))*ok
	isA = abp>=min_amp_bp; if (!any(isA)){if (max(abp)==0) stop('no CN >= amp_cn sequence inside amp'); isA[which.max(abp)] = TRUE}
	isG = ok & !isA & mxcn<amp_cn & cnbp>0 & !overlapsAny(tl,hi+pad)
	px = data.table(i=as.vector(rhdf5::h5read(mcool,paste0(G,'/pixels/bin1_id')))+1L,j=as.vector(rhdf5::h5read(mcool,paste0(G,'/pixels/bin2_id')))+1L,
		v=as.numeric(rhdf5::h5read(mcool,paste0(G,'/pixels/count'))))[ch[i]!=ch[j]]
	OA = OG = numeric(n)
	x = px[isA[i] & isG[j],.(v=sum(v)),by=j]; OA[x$j] = OA[x$j]+x$v; x = px[isG[i] & isA[j],.(v=sum(v)),by=i]; OA[x$i] = OA[x$i]+x$v
	x = px[isG[i] & isG[j]]; x = rbind(x[,.(k=i,v)],x[,.(k=j,v)])[,.(v=sum(v)),by=k]; OG[x$k] = OG[x$k]+x$v
	chi = match(ch,unique(ch)); TA = tapply(cnbp*isA,chi,sum); TG = tapply(cnbp*isG,chi,sum)
	tiles = data.table(seqnames=ch,start=st,end=en,cnbp,OA,OG,EA=cnbp*(sum(TA)-TA[chi]),EG=cnbp*(sum(TG)-TG[chi]))[isG][EA>0 & OG>0][,r:=log2((OA/EA)/(OG/EG))][]
	cm = tiles[,.(med=median(r),n_tiles=.N),by=seqnames]; npos = sum(cm$med>0); p = binom.test(npos,nrow(cm),0.5,alternative='greater')$p.value
	list(call=if (p<alpha) 'ecDNA' else 'not ecDNA',p=p,n_chr_pos=npos,n_chr=nrow(cm),med=median(tiles$r),score=log2((sum(tiles$OA)/sum(tiles$EA))/(sum(tiles$OG)/sum(tiles$EG))),
		n_amp_tiles=sum(isA),tiles=tiles,chr=cm)
}

fastKar_integration_tiles = function(gg,mcool,target,mask=NULL,amp_cn=8,flank=5e4,pad=1e6,chrs=c(1:22,'X'),res=1e6,alpha=0.05,nb=1,depth=NULL){
	# Integration-site caller in one function (user, 2026-11-10), one scale (the mcool's res tiles). One NB fit (as the previous caller).
	# target: the requested amplicon gene (GRanges, e.g. MYCN); mask: GRanges or path to an rds of one
	# 1. tiles: the mcool's res bins on chrs, dropping tiles with >= res / 2 masked. amplicon tile (user: the target only, as the previous caller): the one tile containing the target midpoint;
	#    the highest CN among nodes within flank of the midpoint must be >= amp_cn, otherwise the amplicon is not called (error). Other CN >= amp_cn sequence is neither amplicon nor genome.
	#    genome tiles: every node overlapping the tile has 1 <= CN < amp_cn (CN.bp > 0), and the tile is > pad from any CN >= amp_cn node. CN.bp = sum over the tile's nodes of CN x overlap bp
	# 2. per genome tile g: obs = reads with the amplicon tile (trans and cis); q_g = (OG / EG) / median over tiles, OG = reads with genome tiles on other chromosomes, EG = CN.bp(g) x their summed CN.bp
	#    (the tile's own sensitivity, from its genome contacts)
	# 3. background: mu_g = q_g x (r x E_g + cis_g), E_g = CN.bp(g) x CN.bp(amplicon tile) x interchrom_density; cis_g (g on the target's chromosome) = (depth / 2) x ploidy x
	#    bp(a) x bp(g) x distance_density(|mid(a) - mid(g)|) (fastKar prediction for ploidy unrearranged copies, as the previous caller; depth = estimate.depthratio unless given,
	#    ploidy = bp-weighted graph CN). NB(mu_g, theta), (r, theta) by plain MLE on all genome tiles
	# 4. sites: tiles within nb tiles of a breakpoint (CN change between adjacent nodes, ALT junction end, loose end with CN > 0) with p = P(X >= obs) < alpha / number of such tiles (Bonferroni).
	#    call: 'Integrated' if any site, otherwise 'no site'. Adjacent sites merged into loci
	# Returns list(call, sites, loci, r, theta, cut, depth, ploidy, amp_cn, amp_tile, tiles)
	if (is.character(mask)) mask = readRDS(mask)
	nochr = function(g) GRanges(gsub('^chr','',as.character(seqnames(g))),ranges(g))
	nd = as.data.table(gg$nodes$dt); nodes = nochr(gg$nodes$gr); keep = as.character(seqnames(nodes)) %in% as.character(chrs)
	cn = nd$cn[keep]; nodes = nodes[keep]; hi = reduce(nodes[cn>=amp_cn]); ploidy = sum(cn*as.numeric(width(nodes)))/sum(as.numeric(width(nodes)))
	G = paste0('resolutions/',format(res,scientific=FALSE)); on.exit(rhdf5::h5closeAll())
	ch = gsub('^chr','',as.character(rhdf5::h5read(mcool,paste0(G,'/bins/chrom')))); st = as.vector(rhdf5::h5read(mcool,paste0(G,'/bins/start'))); en = as.vector(rhdf5::h5read(mcool,paste0(G,'/bins/end')))
	tl = GRanges(ch,IRanges(st+1,en)); n = length(tl); bp = as.numeric(en-st); mid = (st+en)/2
	ovbp = function(R){a = numeric(n); if (!length(R)) return(a); o = findOverlaps(tl,R); if (length(o)){s = tapply(as.numeric(width(pintersect(tl[queryHits(o)],R[subjectHits(o)]))),queryHits(o),sum); a[as.integer(names(s))] = s}; a}
	mbp = if (is.null(mask)) numeric(n) else ovbp(reduce(nochr(mask)))
	o = findOverlaps(tl,nodes); w = as.numeric(width(pintersect(tl[queryHits(o)],nodes[subjectHits(o)])))
	cnbp = numeric(n); s = tapply(cn[subjectHits(o)]*w,queryHits(o),sum); cnbp[as.integer(names(s))] = s
	mxcn = numeric(n); s = tapply(cn[subjectHits(o)],queryHits(o),max); mxcn[as.integer(names(s))] = s
	ok = ch %in% as.character(chrs) & mbp<res/2
	target = nochr(target); m = round((min(start(target))+max(end(target)))/2); tg = GRanges(as.character(seqnames(target))[1],IRanges(m-flank,m+flank-1))
	acn = max(c(-Inf,cn[overlapsAny(nodes,tg)])); if (acn<amp_cn) stop(sprintf('fastKar_integration_tiles: target CN %s < amp_cn = %s, not called',acn,amp_cn))
	isA = overlapsAny(tl,GRanges(as.character(seqnames(tg)),IRanges(m,m))) & ok; if (!any(isA)) stop('fastKar_integration_tiles: the target tile is masked or not on chrs')
	isG = ok & !isA & mxcn<amp_cn & cnbp>0 & !overlapsAny(tl,hi+pad)
	px = data.table(i=as.vector(rhdf5::h5read(mcool,paste0(G,'/pixels/bin1_id')))+1L,j=as.vector(rhdf5::h5read(mcool,paste0(G,'/pixels/bin2_id')))+1L,
		v=as.numeric(rhdf5::h5read(mcool,paste0(G,'/pixels/count'))))[(isA[i] & isG[j]) | (isG[i] & isA[j]) | (isG[i] & isG[j] & ch[i]!=ch[j])]
	OA = OG = numeric(n)
	x = rbind(px[isA[i] & isG[j],.(g=j,v)],px[isG[i] & isA[j],.(g=i,v)])[,.(v=sum(v)),by=g]; OA[x$g] = x$v
	x = px[isG[i] & isG[j]]; x = rbind(x[,.(k=i,v)],x[,.(k=j,v)])[,.(v=sum(v)),by=k]; OG[x$k] = x$v
	chi = match(ch,unique(ch)); TG = tapply(cnbp*isG,chi,sum)
	t = data.table(id=seq_len(n),seqnames=ch,start=st,end=en,cnbp,obs=OA,OG,EG=cnbp*(sum(TG)-TG[chi]))[isG][OG>0]
	acnbp = sum(cnbp[isA]) # amplicon tile CN.bp (outside the data.table: inside j, cnbp is the genome-tile column)
	t[,q:=(OG/EG)/median(OG/EG)][,E:=cnbp*acnbp*fastKar::interchrom_density][,cis:=0]
	ia = which(isA); hc = intersect(unique(t$seqnames),unique(ch[ia]))
	if (length(hc)){
		if (is.null(depth)) depth = estimate.depthratio(mcool,ploidy=ploidy,if.chr=hic_is_chr(mcool))
		for (cc in hc){a = ia[ch[ia]==cc]; k = which(t$seqnames==cc)
			t[k,cis:=(depth/2)*ploidy*bp[t$id[k]]*sapply(t$id[k],function(g) sum(bp[a]*distance_density(abs(mid[a]-mid[g]))))]}
	}
	nll = function(par) -sum(dnbinom(t$obs,size=exp(par[2]),mu=t$q*(exp(par[1])*t$E+t$cis),log=TRUE))
	f = optim(c(log(sum(t$obs)/sum(t$q*t$E)),log(2)),nll); rr = exp(f$par[1]); th = exp(f$par[2])
	t[,mu:=q*(rr*E+cis)][,p:=pnbinom(obs-1,size=th,mu=mu,lower.tail=FALSE)]
	nd2 = data.table(seqnames=as.character(seqnames(nochr(gg$nodes$gr))),start=start(gg$nodes$gr),end=end(gg$nodes$gr),cn=nd$cn,
		loose.left=if ('loose.left' %in% names(nd)) nd$loose.left else FALSE,loose.right=if ('loose.right' %in% names(nd)) nd$loose.right else FALSE,
		loose.cn.left=if ('loose.cn.left' %in% names(nd)) nd$loose.cn.left else NA_real_,loose.cn.right=if ('loose.cn.right' %in% names(nd)) nd$loose.cn.right else NA_real_)[order(seqnames,start)]
	b = rbind(nd2[,.(pos=end[-.N],chg=diff(cn)!=0),by=seqnames][chg==TRUE,.(seqnames,pos)],
		nd2[loose.left %in% TRUE & !is.na(loose.cn.left) & loose.cn.left>0,.(seqnames,pos=start)],
		nd2[loose.right %in% TRUE & !is.na(loose.cn.right) & loose.cn.right>0,.(seqnames,pos=end)])
	if (isTRUE(nrow(gg$edges$dt)>0) && any(gg$edges$dt$type=='ALT')){jb = unlist(gg$junctions[type=='ALT']$grl); b = rbind(b,data.table(seqnames=gsub('^chr','',as.character(seqnames(jb))),pos=start(jb)))}
	t[,tested:=if (nrow(b)) overlapsAny(GRanges(seqnames,IRanges(start+1,end)),GRanges(b$seqnames,IRanges(b$pos,width=1))+nb*res) else FALSE]
	cut = alpha/max(1,sum(t$tested)); t[,site:=tested & p<cut]
	sites = t[site==TRUE][order(p),.(seqnames,start,end,obs,mu,cis,p)]
	loci = if (nrow(sites)){sg = GRanges(sites$seqnames,IRanges(sites$start+1,sites$end)); g = reduce(sg); o = findOverlaps(sg,g)
		data.table(locus=subjectHits(o),obs=sites$obs,mu=sites$mu,p=sites$p)[,.(n_tiles=.N,obs=sum(obs),mu=sum(mu),min_p=min(p)),by=locus][,`:=`(seqnames=as.character(seqnames(g))[locus],start=start(g)[locus]-1,end=end(g)[locus])][order(min_p),.(seqnames,start,end,n_tiles,obs,mu,min_p)]} else
		data.table(seqnames=character(),start=numeric(),end=numeric(),n_tiles=integer(),obs=numeric(),mu=numeric(),min_p=numeric())
	list(call=if (nrow(sites)) 'Integrated' else 'no site',sites=sites,loci=loci,r=rr,theta=th,cut=cut,depth=if (is.null(depth)) NA else depth,ploidy=ploidy,amp_cn=acn,amp_tile=tl[isA],tiles=t[,!'id'])
}

# mcool helpers, moved from R/ecfinder.R (archived 2026-11-10); used by fastKar_integration_tiles() for estimate.depthratio(if.chr)
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

hic_is_chr = function(path) {
	if (tools::file_ext(path) != 'mcool') return(FALSE)
	chroms = rhdf5::h5read(path, paste0('/resolutions/', hic_res(path)[1], '/chroms/name'))
	any(grepl('^chr', chroms))
}

## ---- Karyotype hypotheses (restored 2026-11-11 from the archived R/ecfinder.R, ~/projects/karyotyping/archive/fastKar/ecfinder_removed_20261110.R) ----
## boil_greedy(): ecDNA karyotypes (circles peeled greedily from the graph, rest decomposed); squeeze(): HSR karyotypes (lowest amplicon entropy).
## Helpers: find_circles, score_circles, greedypeel, subtract_walk, amp.entropy, amp.entropy.raw, amp_walk_centric, fastKar_centromeres.

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
	if (is(ft,'character')){ft = streduce(parse.gr(ft))}
	amp_nodes = unique((gg.cur$nodes$gr %&% ft)$node.id)
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
		# amplicon fully peeled (100% circular; easy to reach when ft is a single node): done. score_circles would stop on it
		if (sum(gg.cur$nodes$dt[amp_nodes]$cn) <= 0){break}
		scores = score_circles(mycirc,gg.cur,ft)[max.walk.cn>0 & score>0]
	}
	list(circles=peeled,rest=loosefix(gg.cur))
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
			 # circular walk on the amplicon (node test, not gw %&% ft: gUtils::grl.in's aggregate(formula = ) fails on R >= 4.4)
			 if (any(gw$circular & vapply(gw$snode.id,function(x) any(abs(x) %in% amplicon_nodes),logical(1)))){
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
