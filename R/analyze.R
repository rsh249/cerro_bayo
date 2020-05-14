#!/usr/bin/R
#devtools::install_github('rsh249/cRacle')
library(cRacle)
library(dplyr)
library(reshape2)
library(rgbif)
nclus = 24
#Read fossil encoding
cb = read.table('data/fossil_coding.txt', sep='\t', header=T)
cb.melt = melt(cb, "fossil_taxon", c('genus1', 'genus2', 'genus3'), 'grp', na.rm=T)
cb.melt = cb.melt[order(cb.melt$fossil_taxon),]
cb.melt = na.omit(cb.melt)

#get clim
system('wget https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip')
system('mkdir 2.5m')
system('unzip wc2.1_2.5m_bio.zip -d bio2.5')
cl = raster::stack(list.files('bio2.5/', pattern=".tif", full.names = T))

# get data for modern genera
#rgbif with citations

keys <- sapply(cb.melt$value, function(x) name_backbone(name=x, kingdom="Plantae")$genusKey, USE.NAMES=F)
nombre = name_lookup(higherTaxonKey=keys, limit = 3000, status = "accepted")
print(unique(nombre$data$genus))

#cRacle
tax_tab = nombre$data$scientificName
if (file.exists('data/dist.tab')) {
  dat = data.table::fread('data/dist.tab', sep='\t')
  dat.ext = data.table::fread('data/dist_ext.tab', sep='\t')
  
} else {
  dat = getextr(
    tax_tab,
    cl,
    maxrec = 100000,
    schema = 'raw',
    rm.outlier = F,
    nmin = 5,
    parallel = T,
    nclus = nclus
  )
  data.table::fwrite(dat,
                     file = 'data/dist.tab',
                     sep = '\t')
  dat.ext = extraction(dat[, 1:4], cl, schema = 'flat', factor = 16, nmin=15)
  data.table::fwrite(dat.ext,
                     file = 'data/dist_ext.tab',
                     sep = '\t')
}


# Assemble CRACLE compatible distribution table 
#regroup by genus
dat_grp = dat.ext %>% 
  rename(scientificName=tax) %>% 
  left_join(nombre$data, by="scientificName") %>% 
  select(ind_id, genus, lat, lon, scientificName, names(cl)) %>% 
  rename(value=genus) %>%
  rename(name=scientificName) %>%
  left_join(cb.melt, by = 'value') %>%
  select(ind_id, fossil_taxon, lat, lon, names(cl)) %>%
  rename(tax = fossil_taxon) 
dat_grp %>% group_by(tax) %>% summarize(count=n())
dat.ext %>% group_by(tax) %>% summarize(count=n())


#densities
dens_list = dens_obj(
  ex = dat_grp,
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
pdf('first_pass.pdf', height = 6, width = 9)
multiplot(
  dens_list,
  names(cl[[1]]),
  l.pos = 'topleft',
  l.cex = 0.7,
  type = '.kde',
  col=wesanderson::wes_palette('FantasticFox1')
)
abline(v = opt$conintkde[[1]])
text(9, 0.2, cex=0.7, labels = round(opt$conintkde[[1]][1], digits = 1))
text(22, 0.2, cex=0.7, labels = round(opt$conintkde[[1]][2], digits = 1))
abline(
  v = median(opt$origk[[1]]),
  lwd = 2,
  col = wesanderson::wes_palette('Royal1')[1]
)
text(28, 0.3, cex=0.7, labels = paste('Optimum = ', round(median(opt$origk[[1]]), digits =
                                                            1)))
#addplot(and, names(cl[[1]]), type='.kde', col='black')
dev.off()




#save.image('may13.RData')


### Try building pdfs for all taxa and then merge (by or_fun()) to describe each fossil
#densities
dens2 = dens_obj(
  ex = dat.ext,
  clim = cl,
  manip = 'condi',
  n = 256,
  parallel = T,
  nclus = nclus,
  bg.n = 250
) 

#define names map object to map modern taxonomy to fossil taxa

name.map = dat.ext %>% 
  rename(scientificName=tax) %>% 
  left_join(nombre$data, by="scientificName") %>% 
  select(ind_id, genus, lat, lon, scientificName, names(cl)) %>% 
  rename(value=genus) %>%
  rename(name=scientificName) %>%
  left_join(cb.melt, by = 'value') %>%
  select( fossil_taxon, name) %>% unique() %>%
  arrange(fossil_taxon)

dens_list2 = list()
dens2.names = sapply(dens2, function(x) {return(x$name)})
for(i in 1:length(unique(name.map$fossil_taxon))){
  foss = unique(name.map$fossil_taxon)[i]
  print(foss)
  map = name.map[which(name.map$fossil_taxon == foss),]
  targets = which(dens2.names %in% map$name)
  dens_list2[[i]] = or_fun(dens2[targets])
  dens_list2[[i]]$name = foss
}

and2 = and_fun(dens_list2)
opt2 = get_optim(and2)

# summarize and visualize
pdf('equal_weight.pdf', height = 6, width = 9)
multiplot(
  dens_list2,
  names(cl[[4]]),
  l.pos = 'topleft',
  l.cex = 0.7,
  type = '.kde',
  col=wesanderson::wes_palette('FantasticFox1')
)
abline(v = opt2$conintkde[[1]])
text(round(opt2$conintkde[[1]][1], digits = 1)-0.5, 0.2, cex=0.7, labels = round(opt2$conintkde[[1]][1], digits = 1))
text(round(opt2$conintkde[[1]][2], digits = 1)+0.5, 0.2, cex=0.7, labels = round(opt2$conintkde[[1]][2], digits = 1))
abline(
  v = median(opt2$origk[[1]]),
  lwd = 2,
  col = wesanderson::wes_palette('Royal1')[1]
)
text(28, 0.3, cex=0.7, labels = paste('Optimum = ', round(median(opt2$origk[[1]]), digits = 1)))
#addplot(and, names(cl[[1]]), type='.kde', col='black')
dev.off()


#write results tables


## map figure
library(ggplot2)

wc_df = as.data.frame(cl, xy=TRUE)
ggplot() +
  geom_raster(data = wc_df, aes(x = x, y = y, fill = wc2.1_2.5m_bio_1)) +
  geom_point(data=dat_grp, aes(x=lon, y=lat), col=as.integer(dat_grp$tax)) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=c('navy', 'white', 'darkred'),
                       na.value = "black")



