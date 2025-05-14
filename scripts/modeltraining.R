#Script to train the model used throughout fastKar from Rao 2014 data

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

# first look at non-self contacts
dists.and.widths = readRDS(paste0(datafolder,'dists_and_widths.rds'))
dists.and.widths[,distmean:=sum(value)/num,by=c('dist','widthprod')]
distmeans = unique(dists.and.widths[,.(value=distmean,dist,widthprod)])
saveRDS(distmeans,paste0(datafolder,'distmeans.rds'))

distmeans = readRDS(paste0(datafolder,'distmeans.rds'))
maxval = 1e8
library(gridExtra)
nonorm_plot = plot(ggplot(distmeans[dist>0 & dist < maxval])+geom_point(aes(x=1/dist,y=value,color=log10(widthprod))) + scale_x_log10()+scale_y_log10() + labs(x='1/distance [1/bp]',y='counts',color='log(bin area)')+    theme(legend.position = "none"))
norm_plot = plot(ggplot(distmeans[dist>0 & dist < maxval])+geom_point(aes(x=1/dist,y=value/widthprod,color=log10(widthprod))) + scale_x_log10()+scale_y_log10() + labs(x='1/distance [1/bp]',y='counts / bin area [1/bp^2]',color='log(bin area)'))
ppdf(grid.arrange(nonorm_plot,norm_plot,nrow=1),width=10,height=4,filename='karyotype_inference/distance_normalization')

#Okay great, we have area dependence

# what about self-contacts
self.hits = distmeans[dist==0]
ppdf(plot(ggplot(self.hits) + geom_point(aes(x=widthprod,y=value)) + geom_line(aes(x=widthprod,y=1e-2*(widthprod)^(2/3)),linetype='dashed')+ scale_x_log10() + scale_y_log10()
 + labs(x='Bin area [bp^2]',y='Diagonal Hi-C counts')), 
width=5,height=3,filename='karyotype_inference/diagonal_areadep')

#train a spline function for distance decay
dist.splinedat = distmeans[dist>0]
dist.splinedat[,x:=log(dist)]
dist.splinedat[,y:=log(2*value/(1500*widthprod))] #Factor of 2 is important!! 1500x is the depth of the Rao sequencing, but it is "diploid" depth, while we are trying to simulate depth of single alleles
spline_dist = splinefun(dist.splinedat$x,dist.splinedat$y)

#make a plot showing the quality of the fit
plotdat = dist.splinedat[dist < maxval,.(x,y,type='data')]
plotx = linspace(min(plotdat$x),max(plotdat$x))
plotdat2 = data.table(x=plotx,y=spline_dist(plotx),type='fit')
plotdat = rbind(plotdat,plotdat2)
ppdf(plot(ggplot()+geom_point(data=plotdat[type=='data'],aes(x=exp(x),y=exp(y),color=type),size=1) + geom_line(data=plotdat[type=='fit'],aes(x=exp(x),y=exp(y),color=type))+ scale_x_log10() + scale_y_log10() + labs(x='Distance [bp]',y='Count density [bp^-2]')),width=6,height=4,
filename='karyotype_inference/spline_fit')

#train a spline function for area dependence near diagonal
diag.splinedat = distmeans[dist==0]
diag.splinedat[,x:=log(widthprod)]
diag.splinedat[,y:=log(2*value/(1500))] #Again, we have a factor of 2 from the DIPLOID coverage
diag_model <- lm(diag.splinedat$y ~ diag.splinedat$x)
coeffs <- coef(diag_model)
spline_diag = function(x) {coeffs[1] + coeffs[2] * x}
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
interchrom_density = total.interchrom.hits/(interchr.area * 750 * 4)
# each allele has 750x coverage, then interchromosomal hits at position i,j go like density* CN(i) * CN(j) * area

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
noiseval[,expected:=predict.hic(d,tile.ar,750)]

mindist = 1e5
maxdist = 1e8
noiseval = noiseval[d > mindist & d < maxdist]

#let's evenly sample from all tiling sizes
data_per_ar = split(noiseval,by='tile.ar')
unique_areas = noiseval$tile.ar %>% unique
samplesize = 3e5
samples_per_ar = samplesize / length(unique_areas)
subdata = mclapply(data_per_ar,function(x){
  if (nrow(x)<samples_per_ar){
    return(NULL)
  }
  x[,zerokeep := as.integer(numzeros*samples_per_ar/nrow(x))]
  indkeep = sample(1:nrow(x),samples_per_ar)
  subdata_nozeros = x[indkeep,.(value,d,tile.ar,expected,zerokeep)]
  zerorows = unique(subdata_nozeros[,.(value=0,d,tile.ar,expected,zerokeep)])[rep(1:.N,zerokeep)]
  return(rbind(subdata_nozeros,zerorows))
},mc.cores=5) %>% rbindlist

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
