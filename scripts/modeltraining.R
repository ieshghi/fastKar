#Script to train (and test?) the model used throughout fastKar from Rao 2014 data

#make very fine tiling of Rao 2014, keep it to chromosome 1 
datafolder = '~/Projects/karyotype-inference/data/'
hic.file = '/gpfs/commons/groups/imielinski_lab/data/PoreC/Rao2014/4DNFI1UEG1HD.hic'
chrom=c(1)
chrom.gr = parse.gr(c('1'))
fine.res = 1e4
finest.tiling = gr.tile(chrom.gr,fine.res)
hic.data.fine = straw(hic.file,gr=chrom.gr,res=fine.res) 
fine.dat = hic.data.fine$dat %>% setkeyv(c('i','j'))
write.csv(hic.data.fine$dat,paste0(datafolder,'10kb_rao_hic_data.csv'))
saveRDS(hic.data.fine$gr,paste0(datafolder,'10kb_tiling_ranges.rds'))
# or just read in this below
fine.dat = read.csv(paste0(datafolder,'10kb_rao_hic_data.csv'))
finest.tiling = readRDS(paste0(datafolder,'10kb_tiling_ranges.rds'))

#coarse-grain at multiple scales to get a sense of width-dependence along with distance-dependence
hic.data.gm = gM(gr=finest.tiling,dat=hic.data.fine)
resolutions = c(1e4,2e4,5e4,1e5,2e5,5e5,1e6,2e6,5e6)
coarsegrained.data = mclapply(resolutions,function(res){
    tiling = gr.tile(parse.gr(paste0(chrom)),res)
    coarse.data = (hic.data.gm$disjoin(tiling))$agg(tiling)
    return(coarse.data)
  },mc.cores=20,mc.preschedule=F)

saveRDS(coarsegrained.data,paste0(datafolder,'coarsegrained_gmatrices.rds'))

coarsegrained.data = readRDS(paste0(datafolder,'coarsegrained_gmatrices.rds'))
dists.and.widths = mclapply(coarsegrained.data,function(gm){
  tiling = gm$gr
  empty.mat = data.table(expand.grid(tiling$tile.id,tiling$tile.id))[,.(i=Var1,j=Var2,value=0)][j>=i] %>% setkeyv(c('i','j'))
  dat = gm$dat
  gr.dt = gr2dt(gm$gr)
  gr.dt[,mid:=(start + end)/2]
  dat[,width.i:=gr.dt[i]$width]
  dat[,width.j:=gr.dt[j]$width]
  dat[,widthprod:=width.i*width.j]
  dat[,dist:=abs((gr.dt[i]$mid)-(gr.dt[j]$mid))]
  empty.mat[,width.i:=gr.dt[i]$width]
  empty.mat[,width.j:=gr.dt[j]$width]
  empty.mat[,widthprod:=width.i*width.j]
  empty.mat[,dist:=abs((gr.dt[i]$mid)-(gr.dt[j]$mid))]
  empty.mat[,num:=.N,by=c('dist','widthprod')]
  dat.out = merge.data.table(dat,empty.mat[,.(i,j,num)],by=c('i','j'))
  dat.out = dat.out[,.(width.i,width.j,widthprod,dist,value,num)]
  return(dat.out)
  },mc.cores=20) %>% rbindlist

saveRDS(dists.and.widths,paste0(datafolder,'dists_and_widths.rds'))

# let's make the model
# first look at non-self contacts
dists.and.widths = readRDS(paste0(datafolder,'dists_and_widths.rds'))
dists.and.widths[,distmean:=sum(value)/num,by=c('dist','widthprod')]
distmeans = unique(dists.and.widths[,.(value=distmean,dist,widthprod)])
saveRDS(distmeans,paste0(datafolder,'distmeans.rds'))

#what is the depth of the Rao sample
hic.file = '/gpfs/commons/groups/imielinski_lab/data/PoreC/Rao2014/4DNFI1UEG1HD.hic'
wholegenome = gr.tile(si2gr(hg_seqlengths()[1:24]),1e7)
depthest.hic = straw(hic.file,res=as.integer(1e7),gr=wholegenome)
rao.haploid.depth = (depthest.hic$value %>% sum)*300/6e9 #this is HAPLOID depth, depth of sequencing per unique base pair (assuming we can tell the difference between alleles)

distmeans = readRDS(paste0(datafolder,'distmeans.rds'))
maxval = 1e8
plot_distmeans = rbind(distmeans[widthprod< 1e10 & dist > 5e6][seq(1, .N, by = 10)],distmeans[widthprod >= 1e10 | dist <= 5e6])
self.hits = distmeans[dist==0]
library(patchwork)
nonorm_plot = ggplot(plot_distmeans[dist>0 & dist < maxval])+geom_point(aes(x=dist,y=value,color=log10(widthprod))) + scale_x_log10()+scale_y_log10() + labs(x='d [bp]',y='counts',color='log(bin area)')+theme(legend.position = "none")
norm_plot = ggplot(plot_distmeans[dist>0 & dist < maxval])+geom_point(aes(x=dist,y=value/widthprod,color=log10(widthprod))) + scale_x_log10()+scale_y_log10() + labs(x='d [bp]',y='counts / area [1/bp^2]',color='log(bin area)') 
diag_plot = ggplot(self.hits) + geom_point(aes(x=widthprod,y=value)) + geom_line(aes(x=widthprod,y=6e-3*(widthprod)^(2/3)),linetype='dashed')+ scale_x_log10() + scale_y_log10() + labs(x='Bin area [bp^2]',y='Diagonal Hi-C counts') 
#ppdf(print(plot_grid(nonorm_plot,norm_plot,diag_plot,nrow=1,align = "hv", axis = "tblr")),width=12,height=3,filename='karyotype_inference/distance_normalization')
y = (nonorm_plot + norm_plot + diag_plot )
ppdf(print(y
  ),width=14,height=3,filename='karyotype_inference/distance_normalization')
#Okay great, we have area dependence

#train a spline function for distance decay
dist.splinedat = distmeans[dist>0 & widthprod==1e8]
dist.splinedat[,x:=log(dist)]
dist.splinedat[,y:=log(value/(2*rao.haploid.depth*widthprod))] #factor of two because GM12878 is in fact diploid
spline_dist = splinefun(dist.splinedat$x,dist.splinedat$y)

#make a plot showing the quality of the fit
plotdat = dist.splinedat[dist < maxval,.(x,y,type='data')]
plotx = linspace(min(plotdat$x),max(plotdat$x))
plotdat2 = data.table(x=plotx,y=spline_dist(plotx),type='fit')
plotdat = rbind(plotdat,plotdat2)
ppdf(plot(ggplot()+geom_point(data=plotdat[type=='data'],aes(x=exp(x),y=exp(y),color=type),size=1) + geom_line(data=plotdat[type=='fit'],aes(x=exp(x),y=exp(y),color=type))+ scale_x_log10() + scale_y_log10() + labs(x='Distance [bp]',y='Count density [bp^-2]')),width=6,height=4,
filename='karyotype_inference/spline_fit')

diag.splinedat = distmeans[dist==0]
diag.splinedat[,x:=log(widthprod)]
diag.splinedat[,y:=log(value/(2*rao.haploid.depth))] #factor of two because GM12878 is in fact diploid
diag_model <- lm(diag.splinedat$y ~ diag.splinedat$x)
coeffs <- coef(diag_model)
spline_diag <- eval(substitute(
  function(x_new) { a + b * x_new },
  list(a = coeffs[1], b = coeffs[2])
))
splineobj=list(diag_spline = spline_diag,distance_spline = spline_dist) 

ppdf(plot(ggplot(diag.splinedat) + geom_line(aes(x=x,y=spline_diag(x))) +geom_point(aes(x=x,y=y))),width=5,height=4,filename='karyotype_inference/diag_linfit')
saveRDS(splineobj,paste0(datafolder,'spline_model.rds'))
save(splineobj,file='~/git/fastKar/data/scalefree_spline.rda') #save spline for the model to use!

# We also need estimates for inter-chromosomal data. We can stick to a bigger tiling for simplicity
resolution = 1e5
interchroms = parse.gr('1,2')
hic.data.coarse = straw(hic.file,gr=interchroms,res=as.integer(resolution)) 

tiling = hic.data.coarse$gr
tiling$tile.id = 1:length(tiling)
dat = hic.data.coarse$dat
gr.dt = gr2dt(tiling)
dat[,interchrom:=!((gr.dt[i]$seqnames)==(gr.dt[j]$seqnames))]
# assuming density is ~constant, we just calculate total interchromosomal area, get total counts and divide by that
interchr.area = as.numeric(width(parse.gr('1')))*width(parse.gr('2')) #just a check
total.interchrom.hits = dat[interchrom==TRUE]$value %>% sum
interchrom_density = total.interchrom.hits/(interchr.area * rao.haploid.depth * 8)

#just checking we're getting the density more or less right
old_data = readRDS('~/Projects/contacts_sv_correction/data/100kb_lookupfile.rds')
interchrom_density_old = old_data[d<0][,tot_density:=tot_density/1500]

jsonlite::write_json(diag.splinedat,path=paste0(datafolder,'diag_splinedat.json'))
jsonlite::write_json(dist.splinedat,path=paste0(datafolder,'dist_splinedat.json'))
jsonlite::write_json(interchrom_density,path=paste0(datafolder,'background_density.json'))

save(interchrom_density,file='~/git/fastKar/data/interchrom_density.rda') #save density for the model to use!
# estimate dispersion

splineobj = readRDS(paste0(datafolder,'spline_model.rds'))
spline_dist = splineobj$distance_spline
spline_diag = splineobj$diag_spline
predict.hic <- function(dist,area,depth=1){ 
  maxval = 1e8
  maxval_density = depth*exp(spline_dist(log(maxval))) #density at maximum cutoff
  out = depth*area*exp(spline_dist(log(dist)))
  out[dist==0] = depth*exp(spline_diag(log(area[dist==0])))
  out[dist>maxval] = area[dist>maxval]*maxval_density*maxval/dist[dist>maxval] #assume 1/distance scaling for values > 1e8
  return(out)
}

dists.and.widths = readRDS(paste0(datafolder,'dists_and_widths.rds'))
dists.and.widths[,mynum:=.N,by=c('widthprod','dist')] 
#the column "num" gives the total number of pixels in Hi-C data with this given resolution and distance. 
#However zeros have been omitted from the dataset so we need to restore them. 
#bringing all of them back will take too much mem, so we sample the data then restore zeros proportionally
noiseval = dists.and.widths[,.(value,tile.ar = widthprod,d = dist,numzeros=num-mynum)]
noiseval[,expected:=predict.hic(d,tile.ar,rao.depth/2)]

mindist = 1e5
maxdist = 1e8
noiseval = noiseval[d > mindist & d < maxdist]

samplesize = 3e6
samplefrac = samplesize/nrow(noiseval)
noiseval[,zerokeep := as.integer(numzeros*samplefrac)]
indkeep = sample(1:nrow(noiseval),samplesize)
subdata_nozeros = noiseval[indkeep,.(value,d,tile.ar,expected,zerokeep)]
zerorows = unique(subdata_nozeros[,.(value=0,d,tile.ar,expected,zerokeep)])[rep(1:.N,zerokeep)]
subdata= rbind(subdata_nozeros,zerorows)

ppdf(plot(ggplot(subdata) + geom_point(aes(x=value,y=expected,color=log(tile.ar)),size=0.5) + geom_line(aes(x=value,y=value),linetype='dashed',color='red')+ scale_x_log10() + scale_y_log10()),width=6,height=5,'karyotype_inference/multiscale_obsexp')

areas = subdata$tile.ar %>% unique
theta = theta.ml(subdata$value,subdata$expected)
subdata[,nbinom_p := dnbinom(value,size=theta,mu=expected)]
simulated = data.table(value=rnbinom(nrow(subdata),mu=subdata$expected,size=theta),expected=subdata$expected)
simulated[,nbinom_p:=dnbinom(value,size=theta,mu=expected)]
ppdf(qq_pval(subdata$nbinom_p,exp=simulated$nbinom_p),width=5,height=4,filename=paste0('karyotype_inference/multiscale_qqplot'))

saveRDS(list(spline_dist,spline_diag,theta),paste0(datafolder,'splines_and_theta.rds'))

subdata2 = copy(subdata)[,density:=value/tile.ar]
subdata2[,predicted_density:=expected/tile.ar]
subdata2[, tranche := cut(log(predicted_density), breaks = 100)]
result <- unique(subdata2[, .(mexp = mean(predicted_density),err = sd(density)/mean(predicted_density),mdat=mean(density)), by = c('tranche','tile.ar')])

ppdf(plot(ggplot(result[tile.ar <= 1e12])+geom_point(aes(x=mexp,y=mdat,color=log(tile.ar))) + scale_x_log10() + scale_y_log10()),width=5,height=3)
ppdf(plot(ggplot(result)+geom_point(aes(x=mexp,y=err*(tile.ar)^(1/4),color=log(tile.ar))) + scale_x_log10() + scale_y_log10()),width=5,height=3)

# test depth
chr12mb = straw(hic.file,gr=parse.gr('1:1-1e8,2:1-1e8'),res=1e6)
chr1walk = makewalk(parse.gr('1:1-1e8'),list(1,1))
chr1_simdat = forward_simulate(chr1walk,target_region=streduce(chr1walk$grl),pix.size=1e6,depth=rao.depth)
noisy_chr1 = make_noisymap(chr1_simdat,theta=3)[[1]]
chr1data = chr12mb$disjoin(chr1_simdat$gr)$agg(chr1_simdat$gr)
combdat = merge.data.table(noisy_chr1$dat[,.(i,j,value)],chr1data$dat[,.(i,j,value)],by=c('i','j'))
#
combdat[,d:=abs(j-i)]
combdat[,xavg:=mean(value.x),by=c('d')]
combdat[,yavg:=mean(value.y),by=c('d')]
plotcomb = unique(combdat[,.(xavg,yavg)])
ppdf(plot(ggplot(plotcomb) +geom_point(aes(x=xavg,y=yavg)) + geom_line(aes(x=xavg,y=xavg),linetype='dashed') + scale_x_log10() + scale_y_log10()),width=5,height=4)
ppdf(plot(c(chr1data$gtrack(cmap.max=5e4),chr1_simdat$gtrack(cmap.max=5e4)),parse.gr('1')),width=5,height=10)
#
chr1ids = gr2dt(chr12mb$gr)[,id:=.I][seqnames==1]$id
chr2ids = gr2dt(chr12mb$gr)[,id:=.I][seqnames==2]$id
interchrom.dat = rbind(chr12mb$dat[(i %in% chr1ids) & (j %in% chr2ids)],chr12mb$dat[(j %in% chr1ids) & (i %in% chr2ids)])
#
interchrwalk = makewalk(c(parse.gr('1:1-1e8'),parse.gr('2:1-1e8')),list(1,1,2,2))
interchrwalk_sim = forward_simulate(interchrwalk,target_region=streduce(interchrwalk$grl),pix.size=1e6,depth=rao.depth)
chr1ids_sim = gr2dt(interchrwalk_sim$gr)[,id:=.I][seqnames==1]$id
chr2ids_sim = gr2dt(interchrwalk_sim$gr)[,id:=.I][seqnames==2]$id
noisy_interchrom = make_noisymap(interchrwalk_sim,theta=3)[[1]]
interchrom.simdat = rbind(noisy_interchrom$dat[(i %in% chr1ids_sim) & (j %in% chr2ids_sim)],noisy_interchrom$dat[(j %in% chr1ids) & (i %in% chr2ids)])

plotdat = rbind(interchrom.dat[,.(value,type='Hi-C')],interchrom.simdat[,.(value,type='Sim')])

ppdf(plot(ggplot(plotdat[value < 1000]) + geom_histogram(aes(x=value,fill=type),bins=50,color='black',alpha=0.3,position='identity') + geom_vline(aes(xintercept = (interchrom.simdat$value %>% mean)),color='chartreuse4',linetype='dashed')),width=5,height=4)
