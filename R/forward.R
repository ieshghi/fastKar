# In this file I define some functions used in the prediction of 3D structure at rearrangements.

#this branch will be for me to develop a version of the forward model which is tiling-ambiguous, with helper functions to re-do the tiling as in the old version

#possible further optimizations:
# - one of the bottleneck functions is matrix symmetrization. Can i just write cpp code to do the matmul,
#   assuming a symmetric matrix (basically setting i<-->j when j<i)
forward_simulate <- function(walks,target_region = NULL,pix.size=1e5,if.comps=FALSE,mc.cores=1,if.sum=TRUE,depth=1,model=0,if.interchr=T,gm.out=T){
    #TODO:
    #split into a "prep" function which does the merging of the target region and nodes etc.
    #   and then a "match" function mapping the walks to the tiles and flips them if negative strand

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
    } else{ # split node mode
        targetnodes = gr2dt(gr.merge(target_region,nodesgr))
        nodesdt.split = targetnodes[rep(1:.N,each=2)][,lr:=ifelse(mod(.I,2)==1,'l','r')][,.(start,end,seqnames,node.id,cn,lr)][,width:=ifelse(lr=='l',floor((end-start+1)/2),ceil((end-start+1)/2))]
        nodesdt.split[lr=='l',end:=start+width-1]
        nodesdt.split[lr=='r',start:=end-width+1]
        tiled.target = nodesdt.split[,.(start,end,seqnames,node.id,width,cn)][,tile.id:=.I]
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

    # subset to walks that hit target region
    if (inherits(walks,'gWalk')){
        input_walks = walks %&% target_region
    } else{ # if not gWalk, check which walks hit any of the target node.ids
        hit_target = unlist(lapply(walks$snode.id,function(w){length(intersect(unique(abs(w)),tiled.target$node.id)) > 0}))
        input_walks = list(graph = walks$graph,snode.id = walks$snode.id[hit_target], circular = walks$circular[hit_target])
    }
    #then, check for walk CN and if not provided, assume all 1s
    if(is.null(input_walks$dt$cn)){
        walk.cn = rep(1,length(input_walks$snode.id))
    } else{
        walk.cn = input_walks$dt$cn
    }
    # use walk snode.id to generate appropriate ordered sequences of tile.ids to be simulated, making sure to flip order where strand is negative
    tiled.walks.dt = mclapply(input_walks$snode.id,function(walk){
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
    return(simulate_walks_dt(list(walkdt=tiled.walks.dt,walk.cn=walk.cn,circular=input_walks$circular,target=tiled.target),if.comps=if.comps,mc.cores=mc.cores,if.sum=if.sum,depth=depth,model=model,if.interchr=if.interchr,gm.out=gm.out))
}

simulate_walks_dt <- function(tiling.obj,if.comps=FALSE,mc.cores=1,if.sum=TRUE,depth=1,model=0,if.interchr=T,gm.out = T){
    #tiling obj is the output of the "tile_and_prep" or "split_and_prep" objects
    #each element of walk.dt has four columns:
    # width (width of the tile)
    # tile.id (order tiles appear in the walk),
    # orig.id (id of tile in the target region),
    # and compartment
    tiled.target = tiling.obj$target
    walk.dts = tiling.obj$walkdt
    circulars = tiling.obj$circular
    walk.cn = tiling.obj$walk.cn
    target.dat = as.data.table(make_template_dat_cpp(tiled.target)) #useful to make a data.table for all the gMatrices going forward, this is defined on the tiling of the target so that we know exactly what the outputs will look like. Using a single "ID" column to refer to the pixels accelerates summation further on.
    setkeyv(target.dat,c('i','j'))
    #
    reference.intrachr = mclapply(1:length(walk.dts),function(ii){
        this.dt = copy(walk.dts[[ii]]) 
        this.dt[,end:=cumsum(width)] #end coordinates along the allele are just the sum of widths
        this.dt[,start:=end - width + 1] #start coordinates are the ends minus widths
        this.dt[,mid:=(start+end)/2]
        #
        transfmat.spr = sparseMatrix(i=this.dt[!is.na(orig.id)]$tile.id,j=this.dt[!is.na(orig.id)]$orig.id,x=1,dims=c(nrow(this.dt),nrow(tiled.target))) #this sparse matrix gives us the lift back to reference. Since all bins in liftedwalk are subsets of bins in merwalk, this is just a matrix of ones, but each bin in merwalk can have multiple bins in liftedwalk point to it, so each column can have multiple 1s in it.
        #
        local_sim = simulate_hic_dt(this.dt,circulars[ii],model=model) #run the simulator, get a gMatrix in local coordinates

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
            total.counts = inter_allele_dt(reference.intrachr,walk.dts,tiled.target,walk.cn,if.comps=if.comps,mc.cores=mc.cores)
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

simulate_hic_dt  <-  function(walk.dt,is.circular=FALSE,model=0){
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
   return(dat.new[,.(i,j,value,id)])
}

run_analysis <- function(walks,target_region=NULL,if.comps=FALSE,pix.size=1e5,mc.cores=20,verbose=FALSE,if.sum=TRUE,depth=1,model=0){
    lookup_data = fastKar::medium_lookup
    if (pix.size==1e4){
        lookup_data = fastKar::small_lookup 
    }else if(pix.size==1e6){
        lookup_data = fastKar::large_lookup 
    }
    if (is.null(target_region)){
        target_region = walks$footprint
    }
    #
    if (if.comps==FALSE){
        comps.gr = NULL
    }else if(if.comps=='rand'){
        comps.gr = 'rand'
    }else{
        comps.gr = gr.nochr(fastKar::compartment_lookup)
    }
    input_walks = walks %&% target_region
    walk.cn = input_walks$dt$cn
    if(is.null(walk.cn)){
        walk.cn = rep(1,length(input_walks))
        input_walks$set(cn=walk.cn)
    }
    circulars = input_walks$circular
    nw = length(input_walks)
    #
    tiled.target = eval_comps(gr.tile(target_region,pix.size),comps.gr)
    template.dat = make_template_dat(tiled.target) #useful to make a data.table for all the gMatrices going forward, this is defined on the tiling of the target so that we know exactly what the outputs will look like. Using a single "ID" column to refer to the pixels accelerates summation further on.
    #
    if(verbose){message('Simulating contacts')}
    #
    sim.dats = mclapply(1:nw,function(ii){
        walk = input_walks$grl[[ii]][,c()] #take the ii'th input walk
        toflip = which(as.character(strand(walk))=='-')
        merwalk = gr.merge(walk,tiled.target[,c()]) #merge it with the target bins so edges all match
        strand(merwalk) = strand(walk[merwalk$query.id])
        inds = sort(merwalk$query.id,index.return=TRUE)
        merwalk.st = merwalk[inds$ix]
        tile.walk = gr.tile(merwalk.st,pix.size) #tile with pix.size bins
        tile.walk$orig.id = merwalk.st[tile.walk$query.id]$subject.id #which tile in the target region does this tile map to?
        tile.walk$compartment = tiled.target[tile.walk$orig.id]$compartment

        #the steps that are on negative strand need to be flipped
        tile.walk$step.id = merwalk.st[tile.walk$query.id]$query.id #the ID of the "original" bin inside "merwalk" which this bin was taken from
        id.dt = data.table(step = tile.walk$step.id , tile = tile.walk$tile.id)
        id.dt[,new.id:=tile]
        id.dt[step %in% toflip,new.id:=rev(tile),by=step]
        tile.walk$tile.id = id.dt$new.id
        inds = sort(tile.walk$tile.id,index.return=TRUE)
        tile.walk = tile.walk[inds$ix]

        newends = cumsum(width(tile.walk)) #end coordinates along the allele are just the sum of widths
        newstarts = newends - width(tile.walk) + 1 #start coordinates are the ends minus widths
        liftedwalk = GRanges(seqnames=1,ranges=IRanges(start = newstarts,end=newends))#put this all into a granges
        liftedwalk$compartment = tile.walk$compartment #load compartments
        liftedwalk$orig.id = tile.walk$orig.id
        liftedwalk$tile.id = tile.walk$tile.id #the ID of the bin in its own coordinates
        #
        transfmat.spr = sparseMatrix(i=liftedwalk$tile.id,j=liftedwalk$orig.id,x=1,dims=c(length(liftedwalk),length(tiled.target))) #this sparse matrix gives us the lift back to reference. Since all bins in liftedwalk are subsets of bins in merwalk, this is just a matrix of ones, but each bin in merwalk can have multiple bins in liftedwalk point to it, so each column can have multiple 1s in it.
        #
        local_sim = simulate_hic(liftedwalk,lookup_data,if.comps=if.comps,is.circular=circulars[[ii]],model=model) #run the simulator, get a gMatrix in local coordinates
        if (nrow(local_sim$dat)){
            #
            localmat = local_sim$mat %>% unname %>% as.matrix %>% symmetrize #convert to a matrix and symmetrize
            #
            refmat = (Matrix::t(transfmat.spr) %*% localmat)%*%transfmat.spr #apply coordinate transformation 
            #
            tsparse.ref = as(refmat,'TsparseMatrix')
            dt.ref = data.table(i=tsparse.ref@i+1,j=tsparse.ref@j+1,value=tsparse.ref@x)[i<=j]# change to long format, Tsparsematrix is 0-indexed
            out.dt = merge.data.table(data.table::copy(template.dat),dt.ref,all.x=TRUE,by=c('i','j'))[,.(i,j,id,value.y,value=value.x)] # merge with the template data table so the ID matches the relevant coordinates for future summation
            out.dt[!is.na(value.y),value:=value+value.y]
            out.dt = out.dt[,.(i,j,id,value)]
        }else{
            out.dt = data.table::copy(template.dat[,.(i,j,id,value)])
        }
        return(out.dt[,value:=walk.cn[ii]*value])
        #
    },mc.cores=mc.cores)
    #
    if(if.sum){
       if(verbose){message('Adding inter-allele contacts')}
       filled_counts = inter_allele(sim.dats,tiled.target,input_walks,lookup_data,template.dat,mc.cores=mc.cores,model=model)
       #
       if(verbose){message('Summing contacts and adjusting depth') }
       totalcounts = sum_matrices(filled_counts)
       totalcounts[,value:=abs(depth)*value]
       final_mat = gM(dat = totalcounts,gr = tiled.target)
    }else{
       final_mat = lapply(1:nw,function(i){
           outdat = sim.dats[[i]]
           outdat[,value:=abs(depth)*value]
           return(gM(gr = tiled.target,dat = outdat))
       })
    }
    return(final_mat)
}

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

simulate_hic  <-  function(gr,lookup_data,if.comps=FALSE,is.circular=FALSE,model=0){
    grdt = gr2dt(gr)
    widths = grdt$width
    ends = grdt$end
    mids = gr%>%mid
    nbin =  length(gr)
    dat.new = make_template_dat(gr,if.comps=if.comps)
    if (is.circular){
        totl = ends[nbin]
        dat.new[,dist1:= abs(mids[j]-mids[i])]
        dat.new[,this.d:=pmin(dist1,totl-dist1)][,dist1:=NULL]
    }else{
        dat.new[,this.d := abs(mids[j]-mids[i])]
    }
    if(model==0){
        min.res = min(lookup_data[d>0]$d)
        max.res = 1e8
        inds = (dat.new$this.d < max.res) & (dat.new$this.d >= min.res) & !(is.na(dat.new$this.interaction))
        if (sum(inds) > 0){
            dat.new[inds,density:=interp1(x=lookup_data[d>=0 & interaction==this.interaction]$d,y=lookup_data[d>=0 & interaction==this.interaction]$tot_density,xi=this.d),by=c('this.interaction')] 
        }
        dat.new[this.d < min.res,density:=lookup_data[d==0]$tot_density]
        same.minval = min(lookup_data[interaction=='same' & d>0 & d <=max.res]$tot_density)
        diff.minval = min(lookup_data[interaction=='diff' & d>0 & d <=max.res]$tot_density)
        dat.new[this.d >= max.res & this.interaction=='same',density:=(same.minval*max.res/this.d)]        
        dat.new[this.d >= max.res & this.interaction=='diff',density:=(diff.minval*max.res/this.d)]        
   }else {
        min.res = min(lookup_data[d>0]$d)
        dat.new[,density:=1500/(10*min.res^2)*(exp(-this.d/model))]
        #some analytical model for short/long reads with readlength as an input parameter
        #return what the data would be at 1500x depth, adjust accordingly
   }
   dat.new[,value:=density*widthprod]
   dat.new[is.na(value),value:=0]
   gr$tile.id = 1:length(gr)
   return(gM(gr = gr,dat = dat.new[,.(i,j,value,id)]))
}

inter_allele_dt <- function(intra.dts,walk.dts,tiled.target,walk.cn,if.comps=F,mc.cores=1,if.fast=T){
    target.dt = copy(tiled.target)[,.(orig.id=tile.id,cn,compartment)]
    template = data.table::copy(intra.dts[[1]])
    if (if.comps){
        lookup_interchr = fastKar::small_lookup[d<0]
        template[,raw_density:=lookup_interchr[interaction==this.interaction]$tot_density/1500,by=this.interaction]
    }else{
        lookup_interchr = fastKar::small_lookup[d<0 & interaction=='same']
        raw_density = lookup_interchr$tot_density/1500
    }
    # use one-shot calculation formula: n_ij = rate * widthprod * (total_i*total_j - sum_(k walks) (cn_i * cn_j))
    mydat = data.table::copy(template)[,titj:=target.dt$cn[i]*target.dt$cn[j]] 
    walk_nodecounts = t(do.call(cbind,lapply(1:length(walk.dts),function(i){
        w = walk.dts[[i]]
        wcn = walk.cn[i]
        w[,tilecount:=.N,by=orig.id]
        nodecounts.dt = unique(w[!is.na(orig.id),.(orig.id,tilecount)])
        nodecounts.dt[,tilecount:=tilecount*wcn]
        return((merge.data.table(target.dt,nodecounts.dt,by='orig.id',all.x=T)[is.na(tilecount),tilecount:=0]$tilecount))
    })))
    mydat[,sum_sisj :=colSums(walk_nodecounts[,i, drop = FALSE] * walk_nodecounts[,j, drop = FALSE])]
    mydat[,value:=widthprod*raw_density*(titj-sum_sisj)]
    intra_sum = sum_matrices(intra.dts)
    combdat = merge.data.table(mydat[,.(inter_value=value,id)],intra_sum[,.(i,j,intra_value=value,id)],by='id',allow.cartesian=T)[,value:=intra_value+inter_value]
    return(combdat[,.(i,j,id,value)])
}

inter_allele <- function(rebinned_dats,tiled.target,input_walks,lookup_data,template.dat,mc.cores = 10,model=0){
    if (model==0){
    cn.gr = input_walks$graph$gr[,c('cn')]
    all.gr = gr.val(tiled.target,cn.gr,'cn')
    template = data.table::copy(template.dat)
    lookup_interchr = lookup_data[d<0]
    template[,raw_density:=lookup_interchr[interaction==this.interaction]$tot_density,by=this.interaction]
    #
    filled_matrices = mclapply(1:length(input_walks),function(x){
        this.walk = input_walks[x]
        walk.cn = this.walk$dt$cn
        this.grl = this.walk$grl[[1]][,c()]
        this.grl$my.cn = (this.grl %N% (this.grl - 1))
        this.grl = gr.val(this.grl,cn.gr,'cn')
        this.grl$other.cn = this.grl$cn - (this.grl$my.cn %>% replace(is.na(.),0))
        my.gr = gr.val(all.gr,this.grl-1,'my.cn')
        my.gr = gr.val(my.gr,this.grl-1,'other.cn')
        my.gr$other.cn[is.na(my.gr$other.cn)] = my.gr$cn[is.na(my.gr$other.cn)]
        my.gr$my.cn = my.gr$my.cn %>% replace(is.na(.),0)
        if ((sum(my.gr$my.cn < 0,na.rm=TRUE) + sum(my.gr$other.cn < 0,na.rm=TRUE)) > 0){
            print('CN < 0 is not allowed!')
        }
        mydat = data.table::copy(template)
        mydat[,cnprod1:=(my.gr[i]$my.cn)*(my.gr[j]$other.cn)]
        mydat[,cnprod:= ((cnprod1 %>% replace(is.na(.),0)))]
        mydat[,value:=walk.cn*raw_density*cnprod*widthprod]
        mydat[is.na(value),value:=0]
        # subtract CN of myself from the overall CN from jabba
        # then the number of counts for contact i-j is:
            # cts = (CN of i inside myself)*(CN of j in everyone else)*(p from the lookup table(AA/AB/BB))*(binwidth^2) 
        # add the interchromosomal gmat to my self-contact gmat
        combdat = merge.data.table(mydat[,.(inter_value=value,id)],rebinned_dats[[x]][,.(i,j,intra_value=value,id)],by='id')[,value:=intra_value+inter_value]
        out.dat = combdat[,.(i,j,id,value)]
        return(out.dat)
    },mc.cores=mc.cores)
    }else{
        filled_matrices = rebinned_dats
    }
    return(filled_matrices)
}

sum_matrices <- function(matrices){
    allmats = do.call(rbind,matrices)
    allmats[,sumval:=sum(value),by='id']
    outmat = unique(allmats[,.(i,j,id,value = sumval)])
    return(outmat)
}

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

#Rcpp version of the above, no compartments for now
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

symmetrize <- function(input.mat){
    output.mat = input.mat + Matrix::t(input.mat)
    Matrix::diag(output.mat) = Matrix::diag(input.mat)
    return(output.mat)
}

make_compartment_reference <- function(tiling.res,runname,hic.file,comps.file,write.dir,ifplot=FALSE,if.iqr=FALSE){
    #
    comps.gr = gr.nochr(readRDS(comps.file))
    #
    ref.chroms = as.character(c(2,4))
    list.of.ranges = si2gr(hg_seqlengths(chr=FALSE)) %Q% (seqnames %in% ref.chroms)
    #
    hic.data = GxG::straw(hic.file,list.of.ranges,res=tiling.res %>% as.integer)
    bla = make_template_dat(hic.data$gr)
    hic.contacts = merge.data.table(bla[,.(i,j,id)],hic.data$dat[,.(i,j,value)],all.x=TRUE,by=c('i','j'))
    hic.contacts[is.na(value),value:=0]
    #
    hic.ranges = hic.data$gr
    mids = mid(hic.ranges)
    hic.ranges$tile.id = 1:length(hic.ranges)
    hic.contacts[,tile.ar := as.numeric(width(hic.ranges[i]))*as.numeric(width(hic.ranges[j]))]
    #
    message('Evaluating compartments')
    hic.ranges = eval_comps(hic.ranges,comps.gr)
    hic.ranges$chroms=seqnames(hic.ranges)
    hic.contacts[,interaction:=paste0(hic.ranges[i]$compartment,hic.ranges[j]$compartment)]
    hic.contacts[,intra.chr:=(hic.ranges[i]$chroms == hic.ranges[j]$chroms) %>% as.numeric]
    hic.contacts[,d:=abs(mids[j]-mids[i])]
    hic.contacts[intra.chr==0,d:=-1]
    #
    anon.data = hic.contacts[,.(value,interaction,d,tile.ar,id)]
    anon.data[grepl('NA',interaction),interaction:=NA]
    anon.data = anon.data[!is.na(interaction)]
    anon.data[interaction=='BA' | interaction=='AB',interaction:='diff']
    anon.data[interaction=='AA' | interaction=='BB', interaction:='same']
    anon.data[,this_density:=value/tile.ar]
    message('Averaging densities')
    if(if.iqr){
        iqrmean <- function(data,na.rm=TRUE){
            q1 = quantile(data, 0.25)
            q3 = quantile(data, 0.75)
            iqrdat = data[data >= q1 & data <= q3]
            return(mean(iqrdat,na.rm=na.rm))
        }
            anon.data[,tot_density:=iqrmean(this_density,na.rm=TRUE),by=c('interaction','d')]
    }else{
        anon.data[,tot_density:=mean(this_density,na.rm=TRUE),by=c('interaction','d')]
    }
    #
    lookup.data = unique(anon.data[,.(interaction,d,tot_density)])
    maxval = 1e8
    dat_max=lookup.data[d==maxval]
    lookup.data[d>maxval & interaction=='same',tot_density:=(maxval/d)*dat_max[interaction=='same']$tot_density]
    lookup.data[d>maxval & interaction=='diff',tot_density:=(maxval/d)*dat_max[interaction=='diff']$tot_density]
    setkeyv(lookup.data,c('interaction','d'))
    if(ifplot){
    ppdf(plot(ggplot(lookup.data[d>0 & d<maxval])+geom_point(aes(x=d,y=tot_density,color=interaction)) + geom_line(aes(x=d,y=0.05/d),linetype='dashed') + scale_x_log10()+scale_y_log10()+labs(x='Distance [bp]',y='Count density [1/bp^2]') + scale_color_discrete(name='Interaction',labels=c('different compartment','same compartment'))),width=6,height=4,filename=paste0('contact_junction_background/gm12878',runname,'_lookup_ref_distdecay'))}

    noiseval = anon.data[d<maxval,.(value,tot_density,tile.ar,d,interaction,id)]
    noiseval[,expected:=tot_density*tile.ar]
    indkeep = sample(1:nrow(noiseval),1e6)
    subdata = noiseval[indkeep,]
    
    thetaval= theta.ml(subdata$value,subdata$expected)
    subdata[,nbinom_p := dnbinom(value,size=thetaval,mu=expected)]
    simulated = data.table(value=rnbinom(nrow(subdata),mu=subdata$expected,size=thetaval),expected=subdata$expected)
    simulated[,nbinom_p:=dnbinom(value,size=thetaval,mu=expected)]
    if(ifplot){
    ppdf(qq_pval(subdata$nbinom_p,exp=simulated$nbinom_p,max.y=10,max.x=10),width=5,height=4,filename=paste0('contact_junction_background/gm12878',runname,'_model_qqplot'))}
    lookup.data[,theta:=thetaval]

    if (!is.null(write.dir)){
        saveRDS(lookup.data,file=paste0(write.dir,runname,'_lookupfile.rds'))
#        saveRDS(anon.data,file=paste0(write.dir,runname,'_alldata.rds'))
    }
    return(lookup.data)
}


