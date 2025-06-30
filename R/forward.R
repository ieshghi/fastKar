# In this file I define some functions used in the prediction of 3D structure at rearrangements.
#' Simulate the Hi-C data of a set of gWalks in reference coordinates using the fastKar framework
#' 
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param target_region a gRange corresponding to the region of interest for the simulation. If left blank, assume the footprint of the walks.
#' @param pix.size resolution of the simulation tiles, default is 100kb. If set to 0, each tile corresponds to half a node in the gGraph.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @param mc.cores number of cores to use in parallel simulation (defaults to 1)
#' @param if.sum boolean deciding whether to sum the simulated Hi-C signal from all walks or whether to return a list of Hi-C maps (defaults to TRUE)
#' @param depth depth of the simulated Hi-C data, in X coverage (defaults to 1)
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @param if.interchr boolean, sets whether to simulate interchromosomal contacts (default is TRUE). If if.sum = FALSE, then this is set to FALSE. 
#' @param gm.out boolean, chooses whether to return a gMatrix or a data.table (defaults to TRUE, returning a gMatrtix)
#' @return either a gMatrix (gm.out = if.sum = T), a list of gMatrices (gm.out = T, if.sum = F), a data.table (gm.out = F, if.sum = T), or a list of data.tables (gm.out = if.sum = F)
forward_simulate <- function(walks,target_region = NULL,pix.size=1e5,if.comps=FALSE,mc.cores=1,if.sum=TRUE,depth=1,model=0,if.interchr=T,gm.out=T){
    prepped.data = prep_for_sim(walks,target_region,pix.size,if.comps)
    return(simulate_walks(walks,prepped.data$tiled.target,prepped.data$widthdt,if.comps,mc.cores,if.sum,depth,model,if.interchr,gm.out))
}

#' Tile target region and prepare data necessary for Hi-C simulation. 
#'
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param target_region a gRange corresponding to the region of interest for the simulation. If left blank, assume the footprint of the walks.
#' @param pix.size resolution of the simulation tiles, default is 100kb. If set to 0, each tile corresponds to half a node in the gGraph.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @return a list with two elements: "tiled.target", a data.table containing the tiles, their coordinates, and compartment identity, and "widthdt", a data.table containing the widths of all nodes from the input gGraph
prep_for_sim <- function(walks,target_region=NULL,pix.size=1e5,if.comps=F){
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

    #now, tile target region, merging with node boundaries
    if (pix.size > 0){ #tiling mode
        tiled.target = gr2dt(gr.merge(gr.tile(target_region,pix.size),nodesgr))
        tiled.target[,tile.id:=.I]
        tiled.target = tiled.target[,.(start,end,seqnames,tile.id,width,node.id,cn)]
    } else{ # split node mode: tiles have variable size, each node is just two tiles
        targetnodes = gr2dt(gr.merge(target_region,nodesgr))
        nodesdt.split = targetnodes[rep(1:.N,each=2)][,side:=ifelse(mod(.I,2)==1,'left','right')][,.(start,end,seqnames,node.id,cn,side)][,width:=ifelse(side=='left',floor((end-start+1)/2),ceil((end-start+1)/2))]
        nodesdt.split[side=='left',end:=start+width-1]
        nodesdt.split[side=='right',start:=end-width+1]
        tiled.target = nodesdt.split[,.(start,end,seqnames,node.id,side,width,cn)][,tile.id:=.I]
    }

    #check for compartment data
    if (if.comps==FALSE){
        comps.gr = NULL
        tiled.target$compartment = 'A'
    }else{
        stop('Compartments not supported in this version')
        comps.gr = gr.nochr(fastKar::compartment_lookup)
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
#' @param mc.cores number of cores to use in parallel simulation (defaults to 1)
#' @param if.sum boolean deciding whether to sum the simulated Hi-C signal from all walks or whether to return a list of Hi-C maps (defaults to TRUE)
#' @param depth depth of the simulated Hi-C data, in X coverage (defaults to 1)
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @param if.interchr boolean, sets whether to simulate interchromosomal contacts (default is TRUE). If if.sum = FALSE, then this is set to FALSE. 
#' @param gm.out boolean, chooses whether to return a gMatrix or a data.table (defaults to TRUE, returning a gMatrtix)
#' @return either a gMatrix (gm.out = if.sum = T), a list of gMatrices (gm.out = T, if.sum = F), a data.table (gm.out = F, if.sum = T), or a list of data.tables (gm.out = if.sum = F)
simulate_walks <- function(walks,tiled.target,widthdt,if.comps=F,mc.cores=1,if.sum=T,depth=1,model=0,if.interchr=T,gm.out=T){
    # walks is either a gWalk object or a list with elements: snode.id, circular
    # tiled.target is a gr2dt with metadata columns of tile.id, node.id, and cn
    # subset to walks that hit target region
    hit_target = unlist(lapply(walks$snode.id,function(w){length(intersect(unique(abs(w)),tiled.target$node.id)) > 0}))
    input_walks = list(graph = walks$graph,snode.id = walks$snode.id[hit_target], circular = walks$circular[hit_target])
    #then, check for walk CN and if not provided, assume all 1s
    if(is.null(walks$dt$cn)){
        walk.cn = rep(1,length(input_walks$snode.id))
    } else{
        walk.cn = walks$dt$cn[hit_target]
    }
    target.dat = as.data.table(make_template_dat_cpp(tiled.target)) #useful to make a data.table for all the gMatrices going forward, this is defined on the tiling of the target so that we know exactly what the outputs will look like. Using a single "ID" column to refer to the pixels accelerates summation further on.
    setkeyv(target.dat,c('i','j'))

    # use walk snode.id to generate appropriate ordered sequences of tile.ids to be simulated, making sure to flip order where strand is negative
    walk.dts = mclapply(input_walks$snode.id,function(walk){
        walk.dt = data.table(node.id=abs(walk),strand=ifelse(sign(walk)==+1,'+','-'))[,order:=.I]
        walk.dt = merge.data.table(walk.dt,tiled.target,by='node.id',all.x=T,sort=F,allow.cartesian = T)
        toflip = walk.dt$strand=='-'
        walk.dt[,step.id:=.I]
        walk.dt[toflip,step.id:=rev(step.id),by=order]
        setkeyv(walk.dt,'step.id')
        walk.dt = merge.data.table(walk.dt,widthdt,by='node.id',all.x=T)
        walk.dt[is.na(width),width:=nodewidth]
        walk.dt[,nodewidth:=NULL]
        # what about when a tile is outside target_region? pull width from nodesgr. OK if compartment is NA, since contacts with that region dont matter
        id.dt = walk.dt[,.(orig.id=tile.id,compartment,width,step.id)][,tile.id:=step.id][,step.id:=NULL]
        setkeyv(id.dt,'tile.id')
        return(id.dt)
    },mc.cores=mc.cores)
    #
    reference.intrachr = mclapply(1:length(walk.dts),function(ii){
        this.dt = copy(walk.dts[[ii]]) 
        this.dt[,end:=cumsum(width)] #end coordinates along the allele are just the sum of widths
        this.dt[,start:=end - width + 1] #start coordinates are the ends minus widths
        this.dt[,mid:=(start+end)/2]
        #
        transfmat.spr = sparseMatrix(i=this.dt[!is.na(orig.id)]$tile.id,j=this.dt[!is.na(orig.id)]$orig.id,x=1,dims=c(nrow(this.dt),nrow(tiled.target))) #this sparse matrix gives us the lift back to reference. Since all bins in liftedwalk are subsets of bins in merwalk, this is just a matrix of ones, but each bin in merwalk can have multiple bins in liftedwalk point to it, so each column can have multiple 1s in it.
        #
        local_sim = simulate_map(this.dt,input_walks$circular[ii],model=model) #run the simulator, get a gMatrix in local coordinates
        if (nrow(local_sim)){
            #
            localmat = sparseMatrix(i=local_sim$i,j=local_sim$j,x=local_sim$value) %>% symmetrize #convert to a matrix and symmetrize
            #
            refmat = (Matrix::t(transfmat.spr) %*% localmat)%*%transfmat.spr #apply coordinate transformation 
            #
            tsparse.ref = as(refmat,'TsparseMatrix')
            dt.ref = data.table(i=tsparse.ref@i+1,j=tsparse.ref@j+1,value=tsparse.ref@x)[i<=j]# change to long format, Tsparsematrix is 0-indexed
            setkeyv(dt.ref,c('i','j'))
            intra.contacts = copy(target.dat)[dt.ref,on=.(i,j),value:=walk.cn[ii]*i.value] #join with the template hi-c map. Also multiply contacts by copy-number of the walk
        }else{
            intra.contacts = data.table::copy(target.dat[,.(i,j,id,value,widthprod)])
        }
    },mc.cores=mc.cores)
    #
    if(if.sum){
        if (if.interchr & model==0){
            total.counts = calculate_interchrom(reference.intrachr,walk.dts,tiled.target,walk.cn,if.comps=if.comps,mc.cores=mc.cores)
        }else{
            total.counts = sum_matrices(reference.intrachr)
        }
        if (gm.out){
            return(gM(gr = dt2gr(tiled.target),dat=total.counts[,value:=depth*value]))
        } else{
            return(list(gr = tiled.target, dat=total.counts[,value:=depth*value]))
        }
    }else{
        scaled_out = lapply(reference.intrachr,function(x){
            if (gm.out){
                return(gM(gr=dt2gr(tiled.target),dat=x[,value:=depth*value]))
            } else{
                return(dat=x[,value:=depth*value])
            }
        })
        return(scaled_out)
    }
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
        maxval_density = exp(fastKar::splineobj$distance_spline(log(maxval))) #density at maximum cutoff
        dat.new[dist>0,value:=widthprod*exp(fastKar::splineobj$distance_spline(log(dist)))]
        dat.new[dist==0,value:=exp(fastKar::splineobj$diag_spline(log(widthprod)))]
        dat.new[dist>maxval,value:=widthprod*maxval_density*maxval/dist] #assume 1/distance scaling for values > 1e8
   }else {
        min.res = min(fastKar::small_lookup[d>0]$d)
        dat.new[,value:=widthprod*1500/(10*min.res^2)*(exp(-this.d/model))]
        #some analytical model for short/long reads with readlength as an input parameter
        #return what the data would be at 1500x depth, adjust accordingly
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
    target.dt = copy(tiled.target)[,.(orig.id=tile.id,cn,compartment)]
    template = data.table::copy(intra.dts[[1]])
    if (if.comps){
        lookup_interchr = fastKar::small_lookup[d<0]
        template[,raw_density:=lookup_interchr[interaction==this.interaction]$tot_density/1500,by=this.interaction]
    }else{
        raw_density= fastKar::interchrom_density
    }
    # use one-shot calculation formula: n_ij = rate * widthprod * (CN_i*CN_j)
    mydat = data.table::copy(template)[,titj:=target.dt$cn[i]*target.dt$cn[j]] 
    mydat[,value:=widthprod*raw_density*(titj)] #can subtract sum_sisj if no inter-contacts among same chrom nodes
    intra_sum = sum_matrices(intra.dts, mc.cores = mc.cores)
    combdat = merge.data.table(mydat[,.(inter_value=value,id)],intra_sum[,.(i,j,intra_value=value,id)],by='id',allow.cartesian=T)[,value:=intra_value+inter_value]
    return(combdat[,.(i,j,id,value)])
}

#' Sums the contacts coming from many gMatrix data.tables as fast as possible using the 'id' column. 
#' 
#' @param matrices list of data.tables of long-form gMatrices with columns i,j,id,value
#' @return single data.table with same columns, with value column corresponding to the sum of all inputs for a given id
#' NOTE: could be sped up with cpp
sum_matrices <- function(matrices, mc.cores = 1){

  batch_attempts = c(20, 5, 2)
    success = FALSE
    error_msg = NULL

    for (batch_size in batch_attempts) {
      tryCatch({
        message("Trying batch size: ", batch_size)
        batches = split(matrices, ceiling(seq_along(matrices) / batch_size))

        batch_sums = mclapply(batches, function(batch) {
          rbindlist(batch)[, .(value = sum(value)), by = .(i, j, id)]
        }, mc.cores = mc.cores)

        outmat = rbindlist(batch_sums)[, .(value = sum(value)), by = .(i, j, id)]
        success = TRUE
        break
      }, error = function(e) {
        message("Failed at batch size ", batch_size, ": ", e$message)
        error_msg <<- e$message
      })
    }
    if (!success) {
      stop("All batch sizes failed. Last error: ", error_msg)
    }

    
    ## sum in batches
    ## batches <- split(matrices, ceiling(seq_along(matrices) / 10))
    ## batch_sums <- mclapply(batches, function(batch) {
    ##   rbindlist(batch)[, .(value = sum(value)), by = .(i, j, id)]
    ## }, mc.cores = 4)
    ## outmat <- rbindlist(batch_sums)[, .(value = sum(value)), by = .(i, j, id)]

    ## original:
    ## allmats = do.call(rbind,matrices)
    ## allmats[,sumval:=sum(value),by='id']
    ## outmat = unique(allmats[,.(i,j,id,value = sumval)])
    return(outmat)
}

#' Given a set of bins, make the gMatrix data.table corresponding to all pairwise contacts for those bins
#' 
#' @param target.bins either a GRanges, or a data.table. If DT, must have a column "width", giving the width of all the bins
#' @return data.table with columns i,j,id,widthprod corresponding to all pairwise contacts, widthprod is the product of widths of tiles i and j
make_template_dat <- function(target.bins,if.comps=FALSE){
    if(inherits(target.bins,'GRanges')){
        widths = width(target.bins) %>% as.numeric
        l = length(target.bins)
    }else if(inherits(target.bins,'data.table')){
        widths = target.bins$width
        l = nrow(target.bins)
    }
    template.dat = data.table(expand.grid(i=1:l,j=1:l))[j>=i][,value:=0]
    setkeyv(template.dat,c('i','j'))
    if(if.comps){
        target.bins[is.na(target.bins$compartment)]$compartment='A'
        template.dat[,this.interaction:=paste0(target.bins[i]$compartment,target.bins[j]$compartment)]
        template.dat[grepl('NA',this.interaction),this.interaction:=NA]
        template.dat[this.interaction=='BA' | this.interaction=='AB',this.interaction:='diff']
        template.dat[this.interaction=='AA' | this.interaction=='BB',this.interaction:='same']
    } else{
        template.dat[,this.interaction:='same']
    }
    template.dat[,widthprod:=widths[i]*widths[j]]
    template.dat[,id:=.I]
    return(template.dat)
}

#' Rcpp version of make_template_dat(), no compartments for now
#' 
#' @param A a data.table corresponding to the target region tiles, with a "width" column giving the size of each tile
#' @return data.table with columns i,j,id,widthprod corresponding to all pairwise contacts, widthprod is the product of widths of tiles i and j
cppFunction('DataFrame make_template_dat_cpp(DataFrame A) {
  NumericVector widths = A["width"];
  int l = widths.size();
  std::vector<int> i_vals, j_vals;
  std::vector<double> widthprod_vals;

  for (int i = 0; i < l; ++i) { // make the combinations where i<=j
    for (int j = i; j < l; ++j) {
      i_vals.push_back(i + 1);
      j_vals.push_back(j + 1);
      widthprod_vals.push_back(widths[i] * widths[j]);
    }
  }
  return DataFrame::create( //return as a dataframe
    _["i"] = i_vals,
    _["j"] = j_vals,
    _["value"] = NumericVector(i_vals.size(), 0),
    _["widthprod"] = widthprod_vals,
    _["id"] = seq(1, i_vals.size())
  );
}')

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