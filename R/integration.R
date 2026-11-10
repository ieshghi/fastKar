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
