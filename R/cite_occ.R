library(rgbif)
# #occ_search(taxonKey=keys, limit=5, hasCoordinate=TRUE)
# x = occ_download(
#   pred_in("taxonKey", unlist(keys)),
#   user='rsh249',
#   pwd='Ghost4',
#   email='rsh249@cornell.edu'
# )
# 
# d = occ_download_get(x, overwrite=T)
# cit = gbif_citation(d)
# 
# system(paste('unzip', d[[1]]))
# occ = vroom::vroom('occurrence.txt', delim='\t') # fails! columns are not matched with header
# 
