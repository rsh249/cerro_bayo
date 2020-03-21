#!/usr/bin/R
#devtools::install_github('rsh249/cRacle')
library(cRacle)
library(dplyr)
library(reshape2)
library(rgbif)

#Read fossil encoding
cb = read.table('data/fossil_coding.txt', sep='\t', header=T)
cb.melt = melt(cb, "fossil_taxon", c('genus1', 'genus2', 'genus3'), 'grp', na.rm=T)
cb.melt = cb.melt[order(cb.melt$fossil_taxon),]
cb.melt = na.omit(cb.melt)

#get clim
cl = raster::stack(list.files('/usr/share/data/wc2.0/bio2.5/', pattern=".tif", full.names = T))

# get data for modern genera
#rgbif with citations

keys <- sapply(cb.melt$value, function(x) name_backbone(name=x, kingdom="Plantae")$genusKey, USE.NAMES=F)
nombre = name_lookup(higherTaxonKey=keys, limit = 3000, status = "accepted")
print(unique(nombre$data$genus))

#cRacle
tax_tab = nombre$data$scientificName
if (file.exists('data/dist_ext.tab')) {
  dat = read.table('data/dist_ext.tab')
} else {
  dat = getextr(
    tax_tab,
    cl,
    maxrec = 100000,
    schema = 'raw',
    rm.outlier = F,
    nmin = 5,
    parallel = T,
    nclus = 40
  )
  write.table(dat,
              file = 'data/dist.tab',
              sep = '\t')
  dat.ext = extraction(dat[, 1:4], cl, schema = 'flat', factor = 16)
  write.table(dat.ext,
              file = 'data/dist_ext.tab',
              sep = '\t')
}

#read with vroom

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


#densities
dens_list = dens_obj(
  ex = dat_grp,
  clim = cl,
  manip = 'condi',
  n = 256,
  parallel = T,
  nclus = 12,
  bg.n = 250
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

save.image('friday.RData')
