#-----------------------
# Function to process shapeinfo
# It also contains a bunch of functions 
source("functions/namesfunctions.R")

process_shapeinfo <- function(x){
  regs <- translate_regions(x$tx_rgn_descr_nl)
  provs <- translate_provinces(x$tx_prov_descr_nl,
                               x$tx_prov_descr_fr,
                               regs)
  distr <- translate_districts(x$tx_adm_dstr_descr_nl,
                               x$tx_adm_dstr_descr_fr,
                               regs)
  x %>%
    select(cd_munty_refnis,
           names_munty = tx_munty_descr_nl,
           cd_dstr_refnis,
           cd_prov_refnis,
           cd_rgn_refnis) %>%
    mutate(names_dstr = distr,
           names_prov = provs,
           names_rgn = regs)
}
