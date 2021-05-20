#!/usr/bin/R
#devtools::install_github('rsh249/cRacle')
library(cRacle)
library(dplyr)
library(reshape2)
library(rgbif)
library(terra)
library(rgbif)

nclus = 8
#Read fossil encoding
cb = read.table('data_foss/fossil_coding.txt', sep='\t', header=T)
cb.melt = melt(cb, "fossil_taxon", c('genus1', 'genus2', 'genus3'), 'grp', na.rm=T)
cb.melt = cb.melt[order(cb.melt$fossil_taxon),]
cb.melt = na.omit(cb.melt)

#get clim
if(file.exists('wc2.1_2.5m_bio.zip')) {
  
} else {
  download.file('https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip', 'wc2.1_2.5m_bio.zip')
  system('mkdir bio2.5')
  system('unzip wc2.1_2.5m_bio.zip -d bio2.5')
  
}

# ENVIREM
if(file.exists('envirem_global_2.5.zip')){
  
} else {
  download.file('https://deepblue.lib.umich.edu/data/downloads/ms35t870p', 'envirem_global_2.5.zip')
  system('mkdir envirem2.5')
  system('unzip envirem_global_2.5.zip -d envirem2.5')
}


f1 = list.files('envirem2.5', pattern=".tif", full.names = T)
f2 = list.files('bio2.5', pattern=".tif", full.names = T)
cl1 = raster::stack(f1)
cl2 = raster::stack(f2)
ext = raster::extent(cl1)
cl2 = raster::crop(cl2, ext)
cl = raster::stack(cl1, cl2)


# get data for modern genera
#rgbif with citations

keys <- sapply(cb.melt$value, function(x) name_backbone(name=x, kingdom="Viridiplantae")$genusKey, USE.NAMES=F)
nombre = name_lookup(higherTaxonKey=keys, limit = 3000, status = "accepted")
print(unique(nombre$data$genus))

## Get GBIF Citation
# #occ_search(taxonKey=keys, limit=5, hasCoordinate=TRUE)
gbif_auth = read.csv('R/gbif_auth', stringsAsFactors = F)
x = occ_download(
  pred_in("taxonKey", unlist(keys)),
  user=gbif_auth$gbif_user[1],
  pwd=gbif_auth$gbif_pwd[1],
  email=gbif_auth$gbif_email[1]
)
citeGBIF = occ_download_meta(x) %>% gbif_citation()
#### End GBIF Citation




#cRacle
tax_tab = nombre$data$species #? nombre$data$scientificName
if (file.exists('data/dist.tab')) {
  dat = data.table::fread('data/dist.tab', sep='\t')
} else {
  # dat = getextr(
  #   tax_tab,
  #   cl,
  #   maxrec = 100000,
  #   schema = 'raw',
  #   rm.outlier = F,
  #   nmin = 5,
  #   parallel = T,
  #   nclus = nclus
  # )
  
  # do one download for every 100 names
  dat = gbif_dl(tax_tab,
                maxrec= 500000,
                gbif_user=gbif_auth$gbif_user[1],
                gbif_pw=gbif_auth$gbif_pwd[1],
                gbif_email=gbif_auth$gbif_email[1])
  
  # assemble dat from list of data frames
  data.table::fwrite(dat,
                     file = 'data/dist.tab',
                     sep = '\t')

}

dat.ext = extraction(data = dat, clim = cl, schema = 'flat', factor=8, nmin=5)
data.table::fwrite(dat.ext,
                   file = 'data/dist_ext.tab',
                   sep = '\t')

# Assemble CRACLE compatible distribution table 
#regroup by genus
dat_grp = dat.ext %>% 
  rename(species=tax) %>% 
  left_join(nombre$data, by="species") %>% 
  select(ind_id, genus, lat, lon, scientificName, names(cl)) %>% 
  rename(value=genus) %>%
  rename(name=scientificName) %>%
  left_join(cb.melt, by = 'value') %>%
  select(ind_id, fossil_taxon, lat, lon, names(cl)) %>%
  rename(tax = fossil_taxon) %>%
  group_by(tax) %>%
  sample_n(500)

dat_grp %>% group_by(tax) %>% summarize(count=n())
dat.ext %>% group_by(tax) %>% summarize(count=n())


#densities
dens_list = dens_obj(
  ex = dat_grp %>% na.omit,
  clim = cl,
  manip = 'condi',
  n = 1024,
  parallel = T,
  nclus = nclus,
  bg.n = 250, 
  clip='range'
) 

# Rename groups with equal weight / fossil


# get cracle estimates by site
and = and_fun(dens_list)
opt = get_optim(and)

# summarize and visualize
dir.create('figs')
for(zz in 1:length(names(cl))) {
  pdf(paste('figs/', names(cl[[zz]]), ".pdf", sep=''), height = 6, width = 9)
  multiplot(
    dens_list,
    names(cl[[zz]]),
    l.pos = 'topleft',
    l.cex = 0.7,
    type = '.kde',
    col = wesanderson::wes_palette('FantasticFox1')
  )
  addplot(and, names(cl[[zz]]), col = 'black')
  abline(v = opt$conintkde[[zz]])
  text(9,
       0.2,
       cex = 0.7,
       labels = round(opt$conintkde[[zz]][1], digits = 1))
  text(22,
       0.2,
       cex = 0.7,
       labels = round(opt$conintkde[[zz]][2], digits = 1))
  abline(
    v = median(opt$origk[[zz]]),
    lwd = 2,
    col = 'black'
  )
  text(28,
       0.3,
       cex = 0.7,
       labels = paste('Optimum = ', round(median(opt$origk[[zz]]), digits = 1)))
  #addplot(and, names(cl[[1]]), type='.kde', col='black')
  dev.off()
}



est = as.data.frame(opt$conintkde)
hold = cl
for(i in 1:raster::nlayers(cl)){
  hold[[i]] <- ((cl[[i]] > est[1,i]) + (cl[[i]] < est[2,i])) == 2
}

shold = sum(hold)

pdf('figs/map_allvars.pdf')
plot(shold, col=viridis::viridis(99))
dev.off()

pdf('figs/map_allvars_maj_ruls.pdf')
plot(shold>(0.5*35), col=viridis::viridis(99))
dev.off()


# 
# 
# #save.image('may13.RData')
# 
# 
# ### Try building pdfs for all taxa and then merge (by or_fun()) to describe each fossil
# #densities
# dens2 = dens_obj(
#   ex = dat.ext,
#   clim = cl,
#   manip = 'condi',
#   n = 256,
#   parallel = T,
#   nclus = nclus,
#   bg.n = 250
# ) 
# 
# #define names map object to map modern taxonomy to fossil taxa
# 
# name.map = dat.ext %>% 
#   rename(scientificName=tax) %>% 
#   left_join(nombre$data, by="scientificName") %>% 
#   select(ind_id, genus, lat, lon, scientificName, names(cl)) %>% 
#   rename(value=genus) %>%
#   rename(name=scientificName) %>%
#   left_join(cb.melt, by = 'value') %>%
#   select( fossil_taxon, name) %>% unique() %>%
#   arrange(fossil_taxon)
# 
# dens_list2 = list()
# dens2.names = sapply(dens2, function(x) {return(x$name)})
# for(i in 1:length(unique(name.map$fossil_taxon))){
#   foss = unique(name.map$fossil_taxon)[i]
#   print(foss)
#   map = name.map[which(name.map$fossil_taxon == foss),]
#   targets = which(dens2.names %in% map$name)
#   dens_list2[[i]] = or_fun(dens2[targets])
#   dens_list2[[i]]$name = foss
# }
# 
# and2 = and_fun(dens_list2)
# opt2 = get_optim(and2)
# 
# # summarize and visualize
# pdf('equal_weight.pdf', height = 6, width = 9)
# multiplot(
#   dens_list2,
#   names(cl[[4]]),
#   l.pos = 'topleft',
#   l.cex = 0.7,
#   type = '.kde',
#   col=wesanderson::wes_palette('FantasticFox1')
# )
# abline(v = opt2$conintkde[[1]])
# text(round(opt2$conintkde[[1]][1], digits = 1)-0.5, 0.2, cex=0.7, labels = round(opt2$conintkde[[1]][1], digits = 1))
# text(round(opt2$conintkde[[1]][2], digits = 1)+0.5, 0.2, cex=0.7, labels = round(opt2$conintkde[[1]][2], digits = 1))
# abline(
#   v = median(opt2$origk[[1]]),
#   lwd = 2,
#   col = wesanderson::wes_palette('Royal1')[1]
# )
# text(28, 0.3, cex=0.7, labels = paste('Optimum = ', round(median(opt2$origk[[1]]), digits = 1)))
# #addplot(and, names(cl[[1]]), type='.kde', col='black')
# dev.off()
# 
# 
# #write results tables
# 
# 
# ## map figure
library(ggplot2)
# 
wc_df = as.data.frame(cl[[1]], xy=TRUE)
ggplot() +
  geom_raster(data = wc_df, aes(x = x, y = y)) +
  geom_point(data=dat_grp, aes(x=lon, y=lat), col=as.integer(dat_grp$tax)) +
  coord_quickmap() +
  theme_bw() +
  scale_fill_gradientn(colours=c('navy', 'white', 'darkred'),
                       na.value = "black") +
  theme(legend.position = 'bottom')
# 
# 
# 
