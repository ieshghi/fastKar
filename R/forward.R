# In this file I define some functions used in the prediction of 3D structure at rearrangements.

#' Simulate the Hi-C data of a set of gWalks in reference coordinates using the fastKar framework
#' 
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param target_region a gRange corresponding to the region of interest for the simulation. If left blank, assume the footprint of the walks.
#' @param pix.size resolution of the simulation tiles, default is 100kb. If set to 0, each tile corresponds to half a node in the gGraph.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @param mc.cores accepted for compatibility (defaults to 1). The walk kernel runs single-threaded in C++, where forking over walks costs more than it saves (bench/bench_forward_forks.R); results do not depend on it. Parallelize over candidates instead.
#' @param if.sum boolean deciding whether to sum the simulated Hi-C signal from all walks or whether to return a list of Hi-C maps (defaults to TRUE)
#' @param depth depth of the simulated Hi-C data, in X coverage (defaults to 1)
#' @param purity tumor purity/fraction of the sample explained by the supplied walks (defaults to 1). If less than 1, the remaining fraction is simulated as a diploid reference genome over the same context and added to the tumor simulation.
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @param if.interchr boolean, sets whether to simulate interchromosomal contacts (default is TRUE). If if.sum = FALSE, then this is set to FALSE. 
#' @param gm.out boolean, chooses whether to return a gMatrix or a data.table (defaults to TRUE, returning a gMatrtix)
#' @param prepped optional result of prep_for_sim() for this graph, target_region and pix.size, reused instead of tiling again (target_region, pix.size and if.comps are then not used for tiling)
#' @return either a gMatrix (gm.out = if.sum = T), a list of gMatrices (gm.out = T, if.sum = F), a data.table (gm.out = F, if.sum = T), or a list of data.tables (gm.out = if.sum = F)
forward_simulate <- function(walks,target_region = NULL,pix.size=1e5,if.comps=FALSE,mc.cores=1,if.sum=TRUE,depth=1,purity=1,model=0,if.interchr=T,gm.out=T,prepped=NULL){
    if (!is.numeric(purity) || length(purity) != 1 || is.na(purity) || purity < 0 || purity > 1){
        stop('purity must be a single numeric value between 0 and 1')
    }
    if (purity < 1 && !if.sum){
        stop('purity < 1 is currently supported only when if.sum = TRUE')
    }
    haploid.depth = depth/2
    if (is.null(prepped)){
        if (is.null(target_region)){target_region = walks$footprint}
        prepped.data = prep_for_sim(walks,target_region,pix.size,if.comps)
    } else{
        prepped.data = prepped
    }
    if (purity == 1){
        return(simulate_walks(walks,prepped.data$tiled.target,prepped.data$widthdt,if.comps,mc.cores,if.sum,haploid.depth,model,if.interchr,gm.out))
    }
    tumor.sim = simulate_walks(walks,prepped.data$tiled.target,prepped.data$widthdt,if.comps,mc.cores,TRUE,haploid.depth*purity,model,if.interchr,F)
    normal.sim = simulate_diploid_map(prepped.data$tiled.target,if.comps=if.comps,depth=haploid.depth*(1-purity),model=model,if.interchr=if.interchr,gm.out=F)
    if (identical(normal.sim$dat$id, tumor.sim$dat$id) && !is.unsorted(tumor.sim$dat$id, strictly = TRUE)){
        # both maps on the same template, in pixel id order: add the value vectors
        combined.dat = normal.sim$dat[,.(i,j,id)]
        combined.dat[,value:=normal.sim$dat$value + tumor.sim$dat$value]
        setattr(combined.dat,'sorted','id')
    } else{
        combined.dat = merge.data.table(normal.sim$dat[,.(i,j,id,value)],
                                        tumor.sim$dat[,.(id,value)],
                                        by='id',all=TRUE,suffixes=c('.normal','.tumor'))
        combined.dat[is.na(value.normal),normal.value:=0]
        combined.dat[is.na(value.tumor),value.tumor:=0]
        combined.dat = combined.dat[,.(i,j,id,value=value.normal+value.tumor)]
    }
    if (gm.out){
        return(gM(gr=dt2gr(prepped.data$tiled.target),dat=combined.dat))
    } else{
        return(list(gr=prepped.data$tiled.target,dat=combined.dat))
    }
}

#' Tile target region and prepare data necessary for Hi-C simulation. 
#'
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param target_region a gRange corresponding to the region of interest for the simulation. If left blank, assume the footprint of the walks.
#' @param pix.size resolution of the simulation tiles, default is 100kb. If set to 0, each tile corresponds to half a node in the gGraph.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @return a list with two elements: "tiled.target", a data.table containing the tiles, their coordinates, and compartment identity, and "widthdt", a data.table containing the widths of all nodes from the input gGraph
prep_for_sim = function(walks,target_region=NULL,pix.size=1e5,if.comps=F){
    #first, get nodes from graph
    if (is.null(walks$graph) | is.null(walks$graph$gr$cn)){
        stop('Please supply gGraph with CN in input walks object')
    }
    nodesgr = walks$graph$gr[,c('node.id','cn')] 
    nodesgr = nodesgr[strand(nodesgr)=='+']
    widthdt = gr2dt(nodesgr)[,.(nodewidth=width,node.id)]
    #then, subset to target region or make target region if not provided
    if (is.null(target_region)){
        if (is.null(walks$footprint)){
            allnodes = unique(abs(do.call(c,walks$snode.id)))
            target_region = nodesgr[allnodes]
        } else{
            target_region = walks$footprint
        }
    }
    #now, tile target region, merging with node boundaries where the tumor graph is defined.
    #If target_region extends outside the supplied graph, keep those extra tiles with
    #node.id = NA and tumor cn = 0 so a diploid normal component can still be simulated
    #on the full requested region without rebinning.
    # parts of the target not covered by any node, within the graph's seqlengths (the reference took
    # these from the no-query rows of gr.merge(streduce(nodes %&% target), target, all = TRUE))
    bare = function(g) GenomicRanges::GRanges(as.character(GenomicRanges::seqnames(g)), IRanges::ranges(g))
    target.ns = bare(target_region)
    nodes.ns = bare(nodesgr)
    any.in.target = any(IRanges::overlapsAny(nodes.ns, target.ns))
    gaps.ns = suppressWarnings(GenomicRanges::setdiff(target.ns, nodes.ns))
    sl = GenomeInfoDb::seqlengths(nodesgr)[as.character(GenomicRanges::seqnames(gaps.ns))]
    sl[is.na(sl)] = Inf
    keep = GenomicRanges::start(gaps.ns) <= sl
    gaps.ns = gaps.ns[keep]
    sl = sl[keep]
    over = GenomicRanges::end(gaps.ns) > sl
    if (any(over)) GenomicRanges::end(gaps.ns)[over] = as.integer(sl[over])
    outside.target = streduce(gaps.ns)
    tiled.parts = list()
    if (any.in.target){
        tiled.parts[[length(tiled.parts)+1]] = gr2dt(gr.merge(gr.tile(target_region,pix.size),nodesgr))[,.(start,end,seqnames,width,node.id,cn)]
    }
    if (length(outside.target)){
        outside.dt = gr2dt(gr.tile(outside.target,pix.size))[,.(start,end,seqnames,width)]
        outside.dt[,`:=`(node.id=NA_integer_,cn=0)]
        tiled.parts[[length(tiled.parts)+1]] = outside.dt[,.(start,end,seqnames,width,node.id,cn)]
    }
    tiled.target = rbindlist(tiled.parts,use.names=TRUE,fill=TRUE)
    setorder(tiled.target,seqnames,start,end,node.id)
    tiled.target[,tile.id:=.I]
    tiled.target = tiled.target[,.(start,end,seqnames,tile.id,width,node.id,cn)]
    #check for compartment data
    if (if.comps==FALSE){
        comps.gr = NULL
        tiled.target$compartment = 'A'
    }else{
        stop('Compartments not supported in this version')
        comps.gr = gr.nochr(fastKarOptim::compartment_lookup)
        tiled.target = eval_comps(gr.tile(target_region,pix.size),comps.gr) #need to update eval_comps to dt mode
    }
    return(list(tiled.target=tiled.target,widthdt=widthdt))
}

#' Downstream of prep_for_sim, calculate Hi-C data in reference coordinates given sets of walks on a graph and a tiled target region
#' 
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param tiled.target a data.table produced by prep_for_sim(), containing target tile coordinates, tile.ids, and corresponding node.ids, as well as CN and width
#' @param widthdt a data.table produced by prep_for_sim(), containing width data for all the nodes in the gGraph simulated here.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @param mc.cores accepted for compatibility (defaults to 1). The walk kernel runs single-threaded in C++, where forking over walks costs more than it saves (bench/bench_forward_forks.R); results do not depend on it. Parallelize over candidates instead.
#' @param if.sum boolean deciding whether to sum the simulated Hi-C signal from all walks or whether to return a list of Hi-C maps (defaults to TRUE)
#' @param depth depth of the simulated Hi-C data, in X coverage (defaults to 1)
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @param if.interchr boolean, sets whether to simulate interchromosomal contacts (default is TRUE). If if.sum = FALSE, then this is set to FALSE. 
#' @param gm.out boolean, chooses whether to return a gMatrix or a data.table (defaults to TRUE, returning a gMatrtix)
#' @return either a gMatrix (gm.out = if.sum = T), a list of gMatrices (gm.out = T, if.sum = F), a data.table (gm.out = F, if.sum = T), or a list of data.tables (gm.out = if.sum = F)
simulate_walks <- function(walks,tiled.target,widthdt,if.comps=F,mc.cores=1,if.sum=T,depth=1,model=0,if.interchr=T,gm.out=T){
    # walks is either a gWalk object or a list with elements: snode.id, circular
    # tiled.target is a gr2dt with metadata columns of tile.id, node.id, and cn
    # mc.cores is accepted for compatibility; the C++ kernel is single-threaded and results do not depend on it.
    if (if.comps){
        stop('Compartments not supported in this version')
    }
    snode.id = walks$snode.id
    # subset to walks that hit target region
    hit_target = walks_hit_target(snode.id,tiled.target$node.id)
    if (!any(hit_target)){
        stop('none of the walks hits the target region')
    }
    circular = walks$circular[hit_target]
    #then, check for walk CN and if not provided, assume all 1s
    walk.cn.all = walks$dt$cn
    if(is.null(walk.cn.all)){
        walk.cn = rep(1,sum(hit_target))
    } else{
        walk.cn = walk.cn.all[hit_target]
    }
    # every walk as one run of local rows (tiles in walk order, flipped on the minus strand)
    rows = walk_rows(snode.id[hit_target],tiled.target,widthdt)
    kern = forward_kernel_args(model)
    n.tiles = nrow(tiled.target)
    walk_acc = function(w, cn){
        # local rows of walks w (consecutive), walk CNs cn
        r = seq.int(rows$walk_ptr[w[1]] + 1, length.out = rows$walk_ptr[w[length(w)] + 1] - rows$walk_ptr[w[1]])
        ptr = as.integer(rows$walk_ptr[c(w, w[length(w)] + 1)] - rows$walk_ptr[w[1]])
        forward_walks_cpp(ptr, rows$orig[r], rows$width[r], rows$int_width, as.logical(circular[w]), as.numeric(cn), n.tiles,
                          kern$z, kern$lmin, kern$slope, kern$maxval, kern$maxval_density, kern$model, kern$lr_denom, kern$diagf)
    }
    if(if.sum){
        acc = walk_acc(seq_along(walk.cn), walk.cn)
        interchr = if.interchr & model==0
        total.counts = setDT(forward_assemble_cpp(as.numeric(tiled.target$width), as.numeric(tiled.target$cn), acc, depth,
                                                  fastKarOptim::interchrom_density, interchr))
        if (interchr){
            setattr(total.counts,'sorted','id') # calculate_interchrom() returned a table keyed by id
        }
        if (gm.out){
            return(gM(gr = dt2gr(tiled.target),dat=total.counts))
        } else{
            return(list(gr = tiled.target, dat=total.counts))
        }
    }else{
        target.dat = as.data.table(make_template_dat_cpp(tiled.target))
        setattr(target.dat,'sorted',c('i','j'))
        if (gm.out){
            target.gr = dt2gr(tiled.target)
        }
        scaled_out = lapply(seq_along(walk.cn),function(ii){
            x = data.table::copy(target.dat)
            set(x,j='value',value=depth*(walk.cn[ii]*walk_acc(ii, 1)))
            if (gm.out){
                return(gM(gr=target.gr,dat=x))
            } else{
                return(x)
            }
        })
        return(scaled_out)
    }
}

#' Which walks visit at least one node of the target tiling
#'
#' @param snode.id list of signed node id vectors
#' @param target.nodes node.id column of the tiled target (may contain NA)
#' @return logical vector, one per walk
walks_hit_target = function(snode.id,target.nodes){
    lens = base::lengths(snode.id)
    if (!length(lens)) return(logical(0))
    flat = abs(unlist(snode.id,use.names=FALSE))
    hit = flat %in% target.nodes[!is.na(target.nodes)]
    tabulate(rep.int(seq_along(lens), lens)[hit], nbins = length(lens)) > 0
}

#' Local rows (tiles in walk order) of a set of walks, as simulate_walks() lays them out
#'
#' Each node visit expands into the node's tiles in the target tiling, in tiling order (reversed on
#' the minus strand), or into one untiled row with the node's full width if the node has no tiles.
#' @return list(walk_ptr = 0-based row offsets per walk (length W+1), orig = tile.id or NA,
#'   width = row width (numeric), int_width = whether the reference width column is integer)
walk_rows = function(snode.id,tiled.target,widthdt){
    lens = base::lengths(snode.id)
    flat = unlist(snode.id,use.names=FALSE)
    v = abs(flat)
    tt_node = tiled.target$node.id
    tt_rows = which(!is.na(tt_node))
    tt_rows = tt_rows[order(tt_node[tt_rows], tt_rows)]
    rl = rle(tt_node[tt_rows])
    rl_first = cumsum(c(1L, rl$lengths))[seq_along(rl$values)]
    mt = match(v, rl$values)
    ntile = ifelse(is.na(mt), 0L, rl$lengths[mt])
    nrow_visit = pmax(ntile, 1L)
    vis = rep.int(seq_along(v), nrow_visit)
    k = sequence(nrow_visit)
    minus = flat[vis] < 0
    k[minus] = nrow_visit[vis[minus]] - k[minus] + 1L
    tiled = ntile[vis] > 0
    tile = rep(NA_integer_, length(vis))
    tile[tiled] = tt_rows[rl_first[mt[vis[tiled]]] + k[tiled] - 1L]
    width = tiled.target$width[tile]
    untiled = !tiled
    nodewidth = widthdt$nodewidth[match(v[vis[untiled]], widthdt$node.id)]
    width[untiled] = if (is.integer(width)) as.integer(nodewidth) else nodewidth
    walk_ptr = c(0, cumsum(as.numeric(nrow_visit)))[cumsum(c(0, lens)) + 1]
    list(walk_ptr = walk_ptr,
         orig = as.integer(tiled.target$tile.id[tile]),
         width = as.numeric(width),
         int_width = is.integer(width))
}

#' Model constants for the C++ walk kernel, read from the package data at call time
forward_kernel_args = function(model){
    f = fastKarOptim::splineobj$distance_spline
    z = environment(f)$z
    lmin = z$x[1]
    list(z = z, lmin = lmin, slope = f(lmin, deriv = 1), maxval = 1e8,
         maxval_density = distance_density(1e8), model = model,
         lr_denom = if (model > 0) 10*min(fastKarOptim::small_lookup[d>0]$d)^2 else 1,
         diagf = fastKarOptim::splineobj$diag_spline)
}

#' Pixel id of reference pixel (i, j), i <= j, on an n-tile template: its row in make_template_dat()
#'
#' @param i,j integer tile indices (1-based) with i <= j
#' @param n number of tiles
#' @return integer pixel ids
pixel_id = function(i, j, n){
    i = as.numeric(i); j = as.numeric(j)
    as.integer((i - 1)*n - (i - 1)*(i - 2)/2 + (j - i + 1))
}

#' Hi-C contact density per bp^2 at genomic distance dist (bp), from the Rao 2014 distance spline
#'
#' The spline is trained on distances of 10 kb to 249 Mb. Below 10 kb its cubic (fmm) extrapolation
#' diverges as dist -> 0 (~1e10 at 2 bp), which walks that revisit tiny nodes back to back can reach.
#' Below the first knot we instead extend it as a power law with the slope at that knot (~ -1).
#' @param dist numeric vector of distances in bp (> 0)
#' @return numeric vector of contact densities
distance_density = function(dist){
    ldist = log(dist)
    lmin = environment(fastKarOptim::splineobj$distance_spline)$z$x[1]
    out = fastKarOptim::splineobj$distance_spline(pmax(ldist, lmin))
    below = ldist < lmin
    out[below] = out[below] + fastKarOptim::splineobj$distance_spline(lmin, deriv = 1) * (ldist[below] - lmin)
    exp(out)
}

#' Simulate Hi-C data using spline trained on Rao 2014 data
#' 
#' @param walk.dt data.table prepared by simulate_walks(), each row is a tile to simulate and columns are start,end,midpoints of the tile in local walk coordinates.
#' @param is.circular boolean indicating whether the walk is circular or not (defaults to FALSE)
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @return data.table of long-format gMatrix with columns i,j,id,value
simulate_map<-  function(walk.dt,is.circular=FALSE,model=0){
    dat.new = as.data.table(make_template_dat_cpp(walk.dt))
    if (is.circular){
        totl = walk.dt$end[nrow(walk.dt)]
        dat.new[,dist_temp:= abs(walk.dt$mid[j]-walk.dt$mid[i])]
        dat.new[,dist:=pmin(dist_temp,totl-dist_temp)][,dist_temp:=NULL]
    }else{
        dat.new[,dist:= abs(walk.dt$mid[j]-walk.dt$mid[i])]
    }
    if(model==0){
        maxval = 1e8
        maxval_density = distance_density(maxval) #density at maximum cutoff
        dat.new[dist>0,value:=widthprod*distance_density(dist)]
        dat.new[dist==0,value:=exp(fastKarOptim::splineobj$diag_spline(log(widthprod)))]
        dat.new[dist>maxval,value:=widthprod*maxval_density*maxval/dist] #assume 1/distance scaling for values > 1e8
   }else {
        min.res = min(fastKarOptim::small_lookup[d>0]$d)
        dat.new[,value:=widthprod/(10*min.res^2)*(exp(-dist/model))]
        #some analytical model for short/long reads with readlength as an input parameter
   }
   dat.new[is.na(value),value:=0]
#   return(dat.new[,.(i,j,value,id)])
   return(dat.new[,.(i,j,value,id,dist)])
}

#' Evaluate compartment identity of all tiles
#' 
#' TO BE UPDATED, CURRENTLY DEPRECATED
eval_comps <- function(tiled.walk,comps.gr=NULL){
    nbin = length(tiled.walk)
    if (is.null(comps.gr)){
        afrac = rep(1,nbin)
        bfrac = rep(0,nbin)
    } else if(class(comps.gr)=='GRanges'){
        afrac = tiled.walk %O% (comps.gr %Q% (!is.na(compartment) & compartment=='A'))
        bfrac = tiled.walk %O% (comps.gr %Q% (!is.na(compartment) & compartment=='B'))
    } else if(comps.gr=='rand'){
        afrac = sample(c(0,1),nbin,replace=TRUE)
        bfrac = 1-afrac
    }
    compartment = ifelse(afrac > 0.5 & bfrac < 0.5, 'A', ifelse(afrac < 0.5 & bfrac > 0.5, 'B', NA))
    tiled.walk$compartment = compartment
    return(tiled.walk)
}

#' Calculate inter-chromosomal contacts given walks using one-shot formula
#' 
#' @param intra.dts list of data.tables of intra-chromosomal contacts for each walk, generated by simulate_map
#' @param walk.dts list of data.tables corresponding to each walk, with at least one column orig.id
#' @param tiled.target data.table of target region tiles and their coordinates, given by prep_for_sim()
#' @param walk.cn array giving the CN of each walk in walk.dts
#' @param if.comps boolean, whether or not to consider compartments (defaults to FALSE)
#' @return single data.table corresponding to the total contacts in the target region including both intra- and inter-chromosomal contacts
calculate_interchrom <- function(intra.dts,walk.dts,tiled.target,walk.cn,if.comps=F,mc.cores=1,if.fast=T){
    template = intra.dts[[1]]
    if (!if.comps && same_pixels(intra.dts) && !is.unsorted(template$id, strictly = TRUE)){
        # all maps on one template: plain vector sums, no merge
        raw_density = fastKarOptim::interchrom_density
        cn = tiled.target$cn
        intra = Reduce(`+`, lapply(intra.dts, `[[`, 'value'))
        out = template[,.(i,j,id)]
        setattr(out,'sorted',NULL)
        out[,value:=intra + template$widthprod*raw_density*(cn[i]*cn[j])]
        setkeyv(out,'id')
        return(out)
    }
    target.dt = data.table::copy(tiled.target)[,.(orig.id=tile.id,cn,compartment)]
    template = data.table::copy(intra.dts[[1]])
    if (if.comps){
        lookup_interchr = fastKarOptim::small_lookup[d<0]
        template[,raw_density:=lookup_interchr[interaction==this.interaction]$tot_density/1500,by=this.interaction]
    }else{
        raw_density= fastKarOptim::interchrom_density
    }
    # use one-shot calculation formula: n_ij = rate * widthprod * (CN_i*CN_j)
    mydat = data.table::copy(template)[,titj:=target.dt$cn[i]*target.dt$cn[j]] 
    mydat[,value:=widthprod*raw_density*(titj)] #can subtract sum_sisj if no inter-contacts among same chrom nodes
    intra_sum = sum_matrices(intra.dts)
    combdat = merge.data.table(mydat[,.(inter_value=value,id)],intra_sum[,.(i,j,intra_value=value,id)],by='id',allow.cartesian=T)[,value:=intra_value+inter_value]
    return(combdat[,.(i,j,id,value)])
}

#' TRUE if every table in a list has the same pixel ids, in the same order, with no duplicates
same_pixels = function(matrices){
    id1 = matrices[[1]]$id
    if (anyDuplicated(id1)) return(FALSE)
    all(vapply(matrices, function(m) identical(m$id, id1), logical(1)))
}

#' Sums the contacts coming from many gMatrix data.tables as fast as possible using the 'id' column. 
#' 
#' @param matrices list of data.tables of long-form gMatrices with columns i,j,id,value
#' @return single data.table with same columns, with value column corresponding to the sum of all inputs for a given id
sum_matrices <- function(matrices){
    if (same_pixels(matrices)){
        # all maps on one template: add the value vectors
        outmat = matrices[[1]][,.(i,j,id)]
        setattr(outmat,'sorted',NULL)
        outmat[,value:=Reduce(`+`, lapply(matrices, `[[`, 'value'))]
        return(outmat)
    }
    allmats = do.call(rbind,matrices)
    allmats[,sumval:=sum(value),by='id']
    outmat = unique(allmats[,.(i,j,id,value = sumval)])
    return(outmat)
}

#' Given a set of bins, make the gMatrix data.table corresponding to all pairwise contacts for those bins
#' 
#' @param target.bins either a GRanges, or a data.table. If DT, must have a column "width", giving the width of all the bins
#' @return data.table with columns i,j,id,widthprod corresponding to all pairwise contacts, widthprod is the product of widths of tiles i and j
make_template_dat = function(target.bins,if.comps=FALSE){
	if (inherits(target.bins,'GRanges')){
		return(as.data.table(make_template_dat_cpp(gr2dt(target.bins))))
	}
	else if (inherits(target.bins,'data.table')){
		return(as.data.table(make_template_dat_cpp(target.bins)))
	}
}

# make_template_dat_cpp() is compiled with the package (src/forward_kernel.cpp).

#' Symmetrize a matrix
#' 
#' @param input.mat some input matrix, must be square
#' @return matrix of the same shape as input.mat, with same diagonal, and non-diagonal being sum of matrix and its transpose
#' NOTE: could be sped up with cpp
symmetrize <- function(input.mat){
    output.mat = input.mat + Matrix::t(input.mat)
    Matrix::diag(output.mat) = Matrix::diag(input.mat)
    return(output.mat)
}

#' Directly simulate a diploid reference Hi-C map on an existing tiling
#'
#' This is the fast/reference-specialized equivalent of simulating two linear
#' reference homologs through make_diploid_walks() + simulate_walks(), but it
#' takes the already-prepared tile table from prep_for_sim() and computes the
#' diploid signal directly in reference coordinates.  Using the tumor simulation
#' tiling makes the normal and tumor data.tables directly combinable by pixel id.
#'
#' @param tiled.target data.table produced by prep_for_sim(), containing at least start, end, seqnames, width, and tile.id columns
#' @param if.comps boolean setting whether or not to simulate compartments (currently unsupported)
#' @param depth depth multiplier, matching simulate_walks()
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = read length in bp)
#' @param if.interchr boolean, whether to add the same interchromosomal/background term as calculate_interchrom()
#' @param gm.out boolean, chooses whether to return a gMatrix or a list with gr/dat data.tables
#' @return either a gMatrix or a list with tiled target and contact data.table
simulate_diploid_map = function(tiled.target,if.comps=FALSE,depth=1,model=0,if.interchr=TRUE,gm.out=TRUE){
    if (missing(tiled.target) || is.null(tiled.target)){
        stop('tiled.target must be supplied for direct diploid simulation')
    }
    if (if.comps){
        stop('Compartments not supported in direct diploid simulation')
    }
    tiled.target = data.table::copy(as.data.table(tiled.target))
    required.cols = c('start','end','seqnames','width')
    if (!all(required.cols %in% colnames(tiled.target))){
        stop('tiled.target must contain start, end, seqnames, and width columns')
    }
    if (!'tile.id' %in% colnames(tiled.target)){
        tiled.target[,tile.id:=.I]
    }
    tiled.target[,cn:=2]
    tiled.target[,compartment:='A']

    if (model == 0){
        # C++: same-chromosome Hi-C model x 2, interchromosomal term, NA -> 0, depth
        chrom = match(as.character(tiled.target$seqnames), unique(as.character(tiled.target$seqnames)))
        chrom[is.na(tiled.target$seqnames)] = NA_integer_
        kern = forward_kernel_args(0)
        value = forward_diploid_cpp(as.integer(chrom), as.numeric((tiled.target$start+tiled.target$end)/2), as.numeric(tiled.target$width),
                                    kern$z, kern$lmin, kern$slope, kern$maxval, kern$maxval_density,
                                    if.interchr, fastKarOptim::interchrom_density, depth, kern$diagf)
        out.dat = setDT(forward_assemble_cpp(as.numeric(tiled.target$width), numeric(0), value, 1, 0, FALSE))
    } else{
        # reference code path (errors: same_chrom is looked up as a variable by data.table)
        out.dat = as.data.table(make_template_dat_cpp(tiled.target))
        tiled.target[,mid:=(start+end)/2]
        out.dat[,same_chrom:=tiled.target$seqnames[i]==tiled.target$seqnames[j]]
        out.dat[,dist:=abs(tiled.target$mid[j]-tiled.target$mid[i])]
        min.res = min(fastKarOptim::small_lookup[d>0]$d)
        out.dat[same_chrom,value:=2*widthprod/(10*min.res^2)*(exp(-dist/model))]
        out.dat[is.na(value),value:=0]
        out.dat[,value:=depth*value]
        out.dat = out.dat[,.(i,j,id,value)]
        tiled.target[,mid:=NULL]
    }

    if (gm.out){
        return(gM(gr=dt2gr(tiled.target),dat=out.dat))
    } else{
        return(list(gr=tiled.target,dat=out.dat))
    }
}
