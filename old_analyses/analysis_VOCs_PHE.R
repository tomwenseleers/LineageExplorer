# ANALYSIS OF GROWTH ADVANTAGE OF B.1.617 OVER B.1.1.7 BASED ON DATA FROM INCOMING TRAVELLERS FROM INDIA ####
# DATA FROM PHE Technical briefing 9, SARS-CoV-2 variants of concern and variants under investigation in England, 22 April 2021

# last update 8 MAY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-08")
today_num = as.numeric(today)
today # "2021-05-08"
plotdir = "PHE"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# ANALYSIS PHE REPORT NR. 9 ####

phe = read.csv(".//data/PHE/PHE_sequences_from_travellers_Fig3B.csv")
phe$week_commencing = as.Date(phe$week_commencing)
phe$lineage = factor(phe$lineage, levels=c("B.1.1.7","B.1.617","other"))
phe$DATE_NUM = as.numeric(phe$week_commencing)
head(phe)

# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_phe_multi = nnet::multinom(lineage ~ scale(DATE_NUM), data=phe, weights=count, maxit=1000)
fit2_phe_multi = nnet::multinom(lineage ~ ns(DATE_NUM, df=2), data=phe, weights=count, maxit=1000)
BIC(fit1_phe_multi, fit2_phe_multi) 
# fit1_phe_multi fits best (lowest BIC)

# multinomial fits taking into account overdispersion
fit1_phe_mblogit = mblogit(lineage ~ scale(DATE_NUM), data=phe, weights=count, dispersion=TRUE, from.table=TRUE, control=mclogit.control(maxit=100))
fit2_phe_mblogit = nnet::multinom(lineage ~ ns(DATE_NUM, df=2), data=phe, weights=count, dispersion=TRUE, from.table=TRUE, control=mclogit.control(maxit=100))
BIC(fit1_phe_mblogit) # best BIC
BIC(fit2_phe_mblogit)

# estimated growth advantage of B.1.617 vs B.1.1.7 based on multinomial fit
max(as.Date(phe$DATE_NUM)) # 2021-04-05
delta_r_phe = data.frame(confint(emtrends(fit1_phe_multi, trt.vs.ctrl ~ lineage,  
                                             var="DATE_NUM",  mode="latent",
                                             at=list(DATE_NUM=max(phe$DATE_NUM))), 
                                    adjust="none", df=NA)$contrasts)
delta_r_phe
#            contrast     estimate         SE df   asymp.LCL asymp.UCL
# 1 B.1.617 - B.1.1.7  0.157379208 0.03018054 NA  0.09822644 0.2165320
# 2   other - B.1.1.7 -0.006248478 0.01407530 NA -0.03383556 0.0213386


# pairwise tests for differences in growth rates among lineages
contr_phe = data.frame(emtrends(fit1_phe_multi, pairwise ~ lineage|1, 
                                   var="DATE_NUM",  mode="latent",
                                   at=list(DATE_NUM=max(sanger$DATE_NUM)),
                                   adjust="none", df=NA)$contrasts)
contr_phe

# tests for differences in growth rate with UK variant B.1.1.7
contr_phe_againstB117 = data.frame(emtrends(fit1_phe_multi, trt.vs.ctrl ~ lineage|1, 
                                               var="DATE_NUM",  mode="latent",
                                               at=list(DATE_NUM=max(sanger$DATE_NUM)),
                                               adjust="none", df=NA)$contrasts)
contr_phe_againstB117
# contrast     estimate         SE df    z.ratio      p.value
# 1 B.1.617 - B.1.1.7  0.157379208 0.03018054 NA  5.2145925 1.842218e-07
# 2   other - B.1.1.7 -0.006248478 0.01407530 NA -0.4439322 6.570916e-01


# ANALYSIS PHE REPORT NR. 10 ####
# GROWTH ADVANTAGE OF B.1.617.+ VS. B.1.1.7 AMONG NON-TRAVEL RELATED CASES

phe = read.csv(".//data/PHE/PHE_REPORT10_B1617_TRAVELHISTORY.csv")
phe$DATE = as.Date(phe$DATE)
phe$LINEAGE = factor(phe$LINEAGE, levels=c("B.1.1.7","B.1.617.1","B.1.617.2","B.1.617.3"))
phe$DATE_NUM = as.numeric(phe$DATE)
phe=phe[phe$TRAVEL_INDICATOR=="Not travel-associated",]
head(phe)

# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_phe_multi = nnet::multinom(LINEAGE ~ DATE_NUM, data=phe, weights=COUNT, maxit=1000)
fit2_phe_multi = nnet::multinom(LINEAGE ~ ns(DATE_NUM, df=2), data=phe, weights=COUNT, maxit=1000)
BIC(fit1_phe_multi, fit2_phe_multi) 
# fit1_phe_multi fits best (lowest BIC)

# multinomial fits taking into account overdispersion
fit1_phe_mblogit = mblogit(LINEAGE ~ DATE_NUM, data=phe, weights=COUNT, dispersion=TRUE, from.table=TRUE, control=mclogit.control(maxit=100))
fit2_phe_mblogit = nnet::multinom(LINEAGE ~ ns(DATE_NUM, df=2), data=phe, weights=COUNT, dispersion=TRUE, from.table=TRUE, control=mclogit.control(maxit=100))
BIC(fit1_phe_mblogit) # best BIC
BIC(fit2_phe_mblogit)

# estimated growth advantage of B.1.617 vs B.1.1.7 based on multinomial fit
max(as.Date(phe$DATE)) # 2021-04-26
delta_r_phe = data.frame(confint(emtrends(fit1_phe_multi, trt.vs.ctrl ~ LINEAGE,  
                                          var="DATE_NUM",  mode="latent",
                                          at=list(DATE_NUM=max(phe$DATE_NUM))), 
                                 adjust="none", df=NA)$contrasts)
delta_r_phe


# pairwise tests for differences in growth rates among lineages
contr_phe = data.frame(emtrends(fit1_phe_multi, pairwise ~ lineage|1, 
                                var="DATE_NUM",  mode="latent",
                                at=list(DATE_NUM=max(sanger$DATE_NUM)),
                                adjust="none", df=NA)$contrasts)
contr_phe

# tests for differences in growth rate with UK variant B.1.1.7
contr_phe_againstB117 = data.frame(emtrends(fit1_phe_multi, trt.vs.ctrl ~ lineage|1, 
                                            var="DATE_NUM",  mode="latent",
                                            at=list(DATE_NUM=max(sanger$DATE_NUM)),
                                            adjust="none", df=NA)$contrasts)
contr_phe_againstB117
