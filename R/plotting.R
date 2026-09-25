
#make a gwalk in cancer coordinates, colored by chromosome
cancercoord_gw = function(gw){
	keepseqs = unique(as.character(seqnames(gw$footprint)))
	sl = seqlengths(gw)[keepseqs]
	seq_colors = colorspace::qualitative_hcl(length(sl))
	seq_colors = setNames(seq_colors,names(sl))
	fake_seqlengths = c(10000000000)
	names(fake_seqlengths)='1'
	cancercoord = function(gr){
		dt = gr2dt(gr[,c('node.id')])
		dt[,seqlabel:=seqnames]
		dt[,fracpos:=((start+end)/2)/sl[seqnames]]
		dt[,col:=colorspace::lighten(seq_colors[seqnames],amount=fracpos)]
		dt[,seqnames:=1]
		dt[,strand:='+']
		dt[,end:=cumsum(width)]
		dt[,start:=end-width+1]
		dt[,oldnode:=node.id]
		dt$node.id=NULL
		return(dt2gr(dt,seqlengths=fake_seqlengths))
	}
	cancer_list = lapply(gw$grl,function(gr){cancercoord(gr)})
	return(gW(grl=GRangesList(cancer_list)))
}
#takes a list of matrices with shared row / column names and plots them together 
plot_matrices<- function(mat_list,names = NULL,shared_scale = TRUE,palette = 'magma',axis_text_size = 7,title_size = 10) {
    stopifnot(is.list(mat_list))
    stopifnot(all(sapply(mat_list, is.matrix)))
    N <- length(mat_list)
    if (is.null(names)) {names <- paste0("Hypothesis ", seq_len(N))}
    dims <- lapply(mat_list, dim)
    if (!all(sapply(dims, identical, dims[[1]]))) {stop("All matrices must have the same dimensions.")}
    dt <- rbindlist(lapply(seq_along(mat_list), function(i) {
        m <- mat_list[[i]]
        rn <- rownames(m)
        cn <- colnames(m)
        if (is.null(rn)) rn <- paste0("E", seq_len(nrow(m)))
        if (is.null(cn)) cn <- paste0("P", seq_len(ncol(m)))
        out <- as.data.table(as.table(m))
        setnames(out, c("enhancer", "promoter", "score"))
        out[, enhancer := factor(enhancer, levels = rev(rn))]
        out[, promoter := factor(promoter, levels = cn)]
        out[, rearrangement := factor(names[i], levels = names)]
        out
    }))
    p <- ggplot(dt, aes(x = promoter, y = enhancer, fill = score)) +
        geom_tile() +
        facet_wrap(~ rearrangement, nrow = 1) +
        coord_equal() +
        theme_classic(base_size = 10) +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = title_size),
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = axis_text_size
            ),
            axis.text.y = element_text(size = axis_text_size)
        ) +
        labs(
            x = "Promoter",
            y = "Enhancer",
            fill = "E-P score"
        )
    if (shared_scale) {
        p <- p +
            scale_fill_viridis(
		option=palette,
                limits = range(dt$score, na.rm = TRUE)
            )
    } else {
        p <- p +
            scale_fill_viridis(option=palette)
    }
    return(p)
}

#make a nice gtrack object for a list of gwalks
nice_gwgt = function(gw,colorfield='node.id',ywid=0.35){
dodo.call('c',lapply(gw,function(gwi){gwi$gtrack(gr.colorfield=colorfield,lwd=1,height=5,ywid=ywid,labels.suppress=T,border='black')}))
}
