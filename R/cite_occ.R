library(rgbif)
# #occ_search(taxonKey=keys, limit=5, hasCoordinate=TRUE)

gbif_auth = read.csv('R/gbif_auth', stringsAsFactors = F)
x = occ_download(
  pred_in("taxonKey", unlist(keys)),
  user=gbif_auth$gbif_user[1],
  pwd=gbif_auth$gbif_pwd[1],
  email=gbif_auth$gbif_email[1]
)
occ_download_meta(x) %>% gbif_citation()
# 


d = occ_download_get(x, overwrite=T)
cit = gbif_citation(d)

system(paste('unzip', d[[1]]))
occ = vroom::vroom('occurrence.txt', delim='\t') # fails! columns are not matched with header

