#devtools::install_github('rsh249/cRacle')

library(cRacle)
library(rgbif)
#download or load climate data
if(file.exists('wc2.1_2.5m_bio.zip')) {
  
} else {
  download.file('https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip', 'wc2.1_2.5m_bio.zip')
  system('mkdir bio2.5')
  system('unzip wc2.1_2.5m_bio.zip -d bio2.5')
  
}

# try get_extr()
dat = getextr(
  "Pinus strobus",
  cl,
  maxrec = 1000,
  schema = 'raw',
  rm.outlier = F,
  nmin = 5,
  parallel = F
)

keys <- sapply(cb.melt$value, function(x) name_backbone(name=x, kingdom="Viridiplantae")$genusKey, USE.NAMES=F)
nombre = name_lookup(higherTaxonKey=keys, limit = 3000, status = "accepted")
print(unique(nombre$data$genus))

## Get GBIF Citation

# #occ_search(taxonKey=keys, limit=5, hasCoordinate=TRUE)
gbif_auth = read.csv('R/gbif_auth', stringsAsFactors = F)
x = occ_download(
  pred_in("taxonKey", unlist(keys)), ### change code above to get object 'keys' for output of getextr()
  user=gbif_auth$gbif_user[1],
  pwd=gbif_auth$gbif_pwd[1],
  email=gbif_auth$gbif_email[1]
)
citeGBIF = occ_download_meta(x) %>% gbif_citation()
