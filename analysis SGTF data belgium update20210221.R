# ANALYSIS OF S-GENE TARGET FAILURE (S DROPOUT) DATA FROM BELGIUM TO INFER CONTAGIOUSNESS OF NEW VARIANT OF CONCERN B.1.1.7 / 501Y.V1 ####
# PLUS INTERNATIONAL COMPARISON (USING DATA FROM THE UK, DENMARK, SWITZERLAND & THE US)
# AND FIRST ASSESSMENT OF GROWTH ADVANTAGE OF THE SOUTH AFRICAN VOC 501Y.V2 & BRAZILIAN VOC 501Y.V3 BASED ON SEQUENCING DATA
# Tom Wenseleers & Niel Hens
# All Belgian data provided by Emmanuel André

# Data provided by Emmanuel André (BE), COG-UK, PHE & N. Davies (UK), 
# Statens Serum Institut & Danish Covid-19 Genome Consortium (DK, https://www.covid19genomics.dk/statistics), 
# Christian Althaus, Swiss Viollier Sequencing Consortium, Institute of Medical Virology, University of Zurich, 
# Swiss National Covid-19 Science Task Force (Switzerland, https://ispmbern.github.io/covid-19/variants/, 
# https://ispmbern.github.io/covid-19/variants/data & https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv)
# and Helix, San Mateo, CA, Karthik Gangavarapu & Kristian G. Andersen (US, https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data, https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1)

# last update 22 FEBR. 2021

library(lme4)
library(splines)
library(purrr)
library(readxl)
library(effects)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(quantreg)
library(gamm4)
# install from https://github.com/tomwenseleers/export
# library(devtools)
# devtools::install_github("tomwenseleers/export")
library(export) 
library(afex)
library(dfoptim)
library(optimx)
library(lubridate)
library(zoo)
library(gridExtra)
library(sf)
library(broom)
# unloadNamespace("emmeans") # install latest development version of emmeans to add support for mblogit models & to fix bug in v1.5.4 with multinom models
library(devtools)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(broom)
library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)


dat="2021_02_21" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-02-22")
today_num = as.numeric(today)
today # "2021-02-22"

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOC 501Y.V1,VOC 501Y.V2&VOC 501Y.V3 IN BELGIUM BASED ON BASELINE SEQUENCING DATA ####

# data taken from latest Sciensano weekly report "COVID 19 WEKELIJKS EPIDEMIOLOGISCH BULLETIN (19 FEBRUARI 2021), p. 24"
# https://covid-19.sciensano.be/nl/covid-19-epidemiologische-situatie
# (baseline sequencing results, i.e. randomly sampled)

be_seqdata = read.csv(".\\data\\be_latest\\sequencing_501YV1_501YV2_501YV3.csv")
# data is split up in baseline surveillance (randomly sampled) and active surveillance (from travellers, known outbreaks &
# S dropout sequencing), below I will use the randomly sampled baseline surveillance part
be_seqdata$collection_date = as.Date(be_seqdata$collection_date)
be_seqdata$baselinesurv_n_wild_type = be_seqdata$baselinesurv_total_sequenced-be_seqdata$baselinesurv_n_501Y.V1-be_seqdata$baselinesurv_n_501Y.V2-be_seqdata$baselinesurv_n_501Y.V3
be_seqdata$baselinesurv_prop501YV1 = be_seqdata$baselinesurv_n_501Y.V1 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_prop501YV2 = be_seqdata$baselinesurv_n_501Y.V2 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_prop501YV3 = be_seqdata$baselinesurv_n_501Y.V3 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_n_501Y.V1plusV2plusV3 = be_seqdata$baselinesurv_n_501Y.V1+be_seqdata$baselinesurv_n_501Y.V2+be_seqdata$baselinesurv_n_501Y.V3
be_seqdata$baselinesurv_n_501Y.propV1V2V3 = be_seqdata$baselinesurv_n_501Y.V1plusV2plusV3 / be_seqdata$baselinesurv_total_sequenced

be_seqdata$activesurv_n_wild_type = be_seqdata$activesurv_total_sequenced-be_seqdata$activesurv_n_501Y.V1-be_seqdata$activesurv_n_501Y.V2-be_seqdata$activesurv_n_501Y.V3
be_seqdata$activesurv_prop501YV1 = be_seqdata$activesurv_n_501Y.V1 / be_seqdata$activesurv_total_sequenced
be_seqdata$activesurv_prop501YV2 = be_seqdata$activesurv_n_501Y.V2 / be_seqdata$activesurv_total_sequenced
be_seqdata$activesurv_prop501YV3 = be_seqdata$activesurv_n_501Y.V3 / be_seqdata$activesurv_total_sequenced
be_seqdata$activesurv_n_501Y.V1plusV2plusV3 = be_seqdata$activesurv_n_501Y.V1+be_seqdata$activesurv_n_501Y.V2+be_seqdata$activesurv_n_501Y.V3
be_seqdata$activesurv_n_501Y.propV1V2V3 = be_seqdata$activesurv_n_501Y.V1plusV2plusV3 / be_seqdata$activesurv_total_sequenced


be_seqdata$basplusactivesurv_n_501Y.V1 = be_seqdata$baselinesurv_n_501Y.V1+be_seqdata$activesurv_n_501Y.V1
be_seqdata$basplusactivesurv_n_501Y.V2 = be_seqdata$baselinesurv_n_501Y.V2+be_seqdata$activesurv_n_501Y.V2
be_seqdata$basplusactivesurv_n_501Y.V3 = be_seqdata$baselinesurv_n_501Y.V3+be_seqdata$activesurv_n_501Y.V3
be_seqdata$basplusactivesurv_total_sequenced = be_seqdata$baselinesurv_total_sequenced+be_seqdata$activesurv_total_sequenced 
be_seqdata$basplusactivesurv_n_wild_type = be_seqdata$basplusactivesurv_total_sequenced-be_seqdata$basplusactivesurv_n_501Y.V1-be_seqdata$basplusactivesurv_n_501Y.V2-be_seqdata$basplusactivesurv_n_501Y.V3
be_seqdata$basplusactivesurv_prop501YV1 = be_seqdata$basplusactivesurv_n_501Y.V1 / be_seqdata$basplusactivesurv_total_sequenced
be_seqdata$basplusactivesurv_prop501YV2 = be_seqdata$basplusactivesurv_n_501Y.V2 / be_seqdata$basplusactivesurv_total_sequenced
be_seqdata$basplusactivesurv_prop501YV3 = be_seqdata$basplusactivesurv_n_501Y.V3 / be_seqdata$basplusactivesurv_total_sequenced
be_seqdata$basplusactivesurv_n_501Y.V1plusV2plusV3 = be_seqdata$basplusactivesurv_n_501Y.V1+be_seqdata$basplusactivesurv_n_501Y.V2+be_seqdata$basplusactivesurv_n_501Y.V3
be_seqdata$basplusactivesurv_n_501Y.propV1V2V3 = be_seqdata$basplusactivesurv_n_501Y.V1plusV2plusV3 / be_seqdata$basplusactivesurv_total_sequenced


head(be_seqdata)


# BASELINE SURVEILLANCE DATA ####

be_basseqdata_long = gather(be_seqdata[,c("collection_date",
                                          "baselinesurv_n_wild_type",
                                          "baselinesurv_n_501Y.V1",
                                          "baselinesurv_n_501Y.V2",
                                          "baselinesurv_n_501Y.V3",
                                          "baselinesurv_n_501Y.V1plusV2plusV3",
                                          "baselinesurv_total_sequenced")], 
                            variant, count, c("baselinesurv_n_wild_type",
                                              "baselinesurv_n_501Y.V1",
                                              "baselinesurv_n_501Y.V2",
                                              "baselinesurv_n_501Y.V3",
                                              "baselinesurv_n_501Y.V1plusV2plusV3"), factor_key=TRUE)
be_basseqdata_long$variant = factor(be_basseqdata_long$variant, 
                                    levels=c("baselinesurv_n_wild_type","baselinesurv_n_501Y.V1","baselinesurv_n_501Y.V2","baselinesurv_n_501Y.V3","baselinesurv_n_501Y.V1plusV2plusV3"), 
                                    labels=c("wild type", "501Y.V1", "501Y.V2", "501Y.V3", "501Y.V1+V2+V3"))
be_basseqdata_long$collection_date_num = as.numeric(be_basseqdata_long$collection_date)
be_basseqdata_long$prop = be_basseqdata_long$count / be_basseqdata_long$baselinesurv_total_sequenced

baseline_sequencing = ggplot(data=be_basseqdata_long[be_basseqdata_long$variant!="501Y.V1+V2+V3",], 
                             aes(x=collection_date, 
                                 y=count, fill=variant, group=variant)) +
  # facet_wrap(~LABORATORY) +
  geom_area(aes(fill=variant), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("variant", values=c("darkgrey","blue","red","green3"), 
                    labels=c("wild type","501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-03-01")), expand=c(0,0)) +
  ylab("Share among newly diagnosed infections") +
  xlab("Collection date") +
  # ggtitle("Test outcomes") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right")
baseline_sequencing

saveRDS(baseline_sequencing, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.pdf"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
be_seq_mfit0 = nnet::multinom(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                              subset=be_basseqdata_long$variant!="501Y.V1+V2+V3", maxit=1000)
be_seq_mfit1 = nnet::multinom(variant ~ scale(ns(collection_date_num, df=2), center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                              subset=be_basseqdata_long$variant!="501Y.V1+V2+V3", maxit=1000) 
be_seq_mfit2 = nnet::multinom(variant ~ scale(ns(collection_date_num, df=3), center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                              subset=be_basseqdata_long$variant!="501Y.V1+V2+V3", maxit=1000) 
be_seq_mfitallVOC = nnet::multinom(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                                   subset=be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="501Y.V1+V2+V3", maxit=1000) 
BIC(be_seq_mfit0, be_seq_mfit1, be_seq_mfit2) # mfit0 fits best
#              df      BIC
# be_seq_mfit0  6 2405.575
# be_seq_mfit1  9 2420.704
# be_seq_mfit2 12 2443.147
summary(be_seq_mfit0)

# growth rate advantage compared to wild type
delta_r_501V1_501YV2_501YV3 = data.frame(confint(emtrends(be_seq_mfit0, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_501V1_501YV2_501YV3) = delta_r_501V1_501YV2_501YV3[,"contrast"]
delta_r_501V1_501YV2_501YV3 = delta_r_501V1_501YV2_501YV3[,-1]
delta_r_501V1_501YV2_501YV3
#                       estimate   asymp.LCL  asymp.UCL
# 501Y.V1 - wild type 0.08294642 0.06919210 0.09670073
# 501Y.V2 - wild type 0.07684561 0.04674708 0.10694414
# 501Y.V3 - wild type 0.18033281 0.03023209 0.33043353

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_501V1_501YV2_501YV3*4.7) 
#                     estimate asymp.LCL asymp.UCL
# 501Y.V1 - wild type 1.476757  1.384311  1.575375
# 501Y.V2 - wild type 1.435014  1.245717  1.653075
# 501Y.V3 - wild type 2.333955  1.152681  4.725803

# for all 3 variants together
# growth rate advantage compared to wild type
delta_r_allVOCs = data.frame(confint(emtrends(be_seq_mfitallVOC, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_allVOCs) = delta_r_allVOCs[,"contrast"]
delta_r_allVOCs = delta_r_allVOCs[,-1]
delta_r_allVOCs
#                               estimate  asymp.LCL  asymp.UCL
# (501Y.V1+V2+V3) - wild type 0.08273544 0.06978487 0.09568601

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_allVOCs*4.7) 
#                             estimate asymp.LCL asymp.UCL
# (501Y.V1+V2+V3) - wild type   1.475293  1.388174   1.56788




# PS: mblogit fit would also be possible & would take into account overdispersion
be_basseqdata_long$obs = factor(1:nrow(be_basseqdata_long))
be_seq_mblogitfit = mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
                            # random = ~ 1|obs,
                            weights = count, data=be_basseqdata_long, 
                            subset=be_basseqdata_long$variant!="501Y.V1+V2+V3",
                            dispersion = FALSE)
be_seq_mblogitfit1 = mblogit(variant ~ scale(ns(collection_date_num, df=2), center=TRUE, scale=FALSE),
                             # random = ~ 1|obs,
                             weights = count, data=be_basseqdata_long, 
                             subset=be_basseqdata_long$variant!="501Y.V1+V2+V3",
                             dispersion = FALSE)
BIC(be_seq_mblogitfit, be_seq_mblogitfit1)
#                    df      BIC
# be_seq_mblogitfit   6 2405.575
# be_seq_mblogitfit1  9 2420.665
summary(be_seq_mblogitfit)

dispersion(mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
                   # random = ~ 1|obs,
                   weights = count, data=be_basseqdata_long,
                   subset = be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="501Y.V1+V2+V3",
                   dispersion = TRUE), method="Afroz") # dispersion coefficient = 3.2

# growth rate advantage compared to wild type
delta_r_501V1_501YV2_501YV3 = data.frame(confint(emtrends(be_seq_mblogitfit, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent", df=NA), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_501V1_501YV2_501YV3) = delta_r_501V1_501YV2_501YV3[,"contrast"]
delta_r_501V1_501YV2_501YV3 = delta_r_501V1_501YV2_501YV3[,-1]
delta_r_501V1_501YV2_501YV3
#                       estimate  asymp.LCL  asymp.UCL
# 501Y.V1 - wild type 0.08294652 0.06919221 0.09670084
# 501Y.V2 - wild type 0.07684610 0.04674750 0.10694470
# 501Y.V3 - wild type 0.18032114 0.03022225 0.33042003




# plot multinomial model fit

# library(effects)
# plot(Effect("collection_date_num",be_seq_mfit0), style="stacked")

extrapolate = 90
date.from = as.numeric(as.Date("2020-11-01")) # min(be_basseqdata_long$collection_date_num)
date.to = max(be_basseqdata_long$collection_date_num)+extrapolate

be_seq_mfit0_preds = data.frame(emmeans(be_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfit0_preds$collection_date = as.Date(be_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
be_seq_mfit0_preds$variant = factor(be_seq_mfit0_preds$variant, levels=c("wild type","501Y.V1","501Y.V2","501Y.V3"),
                                    labels=c("wild type","501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)"))

be_seq_mfitallVOCs_preds = data.frame(emmeans(be_seq_mfitallVOC, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfitallVOCs_preds$collection_date = as.Date(be_seq_mfitallVOCs_preds$collection_date_num, origin="1970-01-01")
be_seq_mfitallVOCs_preds$variant = factor(be_seq_mfitallVOCs_preds$variant, levels=c("wild type","501Y.V1+V2+V3"),
                                          labels=c("wild type","501Y.V1+V2+V3"))

be_basseqdata_long2 = be_basseqdata_long[be_basseqdata_long$variant!="wild type",]
be_basseqdata_long2$variant = droplevels(be_basseqdata_long2$variant) 
be_basseqdata_long2$variant = factor(be_basseqdata_long2$variant, levels=c("501Y.V1","501Y.V2","501Y.V3","501Y.V1+V2+V3"),
                                     labels=c("501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)","501Y.V1+V2+V3"))

muller_be_seq_mfit0 = ggplot(data=be_seq_mfit0_preds, 
                             aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(be_basseqdata_long$collection_date)+1, 
           xmax=as.Date("2021-03-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("variant", values=c("darkgrey","blue","red","green3")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-03-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Relative abundance") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")
muller_be_seq_mfit0

saveRDS(muller_be_seq_mfit0, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.pdf"), width=7, height=5)

multinom_501YV1_501YV2_501YV3 = ggarrange(baseline_sequencing+ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")+xlab("")+ylab("Share"), 
                                          muller_be_seq_mfit0+ggtitle("Multinomial fit plus extrapolation")+ylab("Share"), ncol=1, 
                                          common.legend = TRUE, legend="bottom")
multinom_501YV1_501YV2_501YV3
saveRDS(multinom_501YV1_501YV2_501YV3, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.pptx"), width=7, height=10)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.png"), width=7, height=10)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.pdf"), width=7, height=10)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS
be_seq_mfit0_preds2 = rbind(be_seq_mfit0_preds[be_seq_mfit0_preds$variant!="wild type",],
                            be_seq_mfitallVOCs_preds[be_seq_mfitallVOCs_preds$variant!="wild type",])
be_seq_mfit0_preds2$variant = droplevels(be_seq_mfit0_preds2$variant)


# on response scale:
plot_multinom_501YV1_501YV2_501YV3_response = qplot(data=be_seq_mfit0_preds2, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share of diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(be_basseqdata_long2$collection_date), as.Date("2021-03-01")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("blue","red","green3","black")) +
  scale_colour_manual("variant", values=c("blue","red","green3","black")) +
  geom_point(data=be_basseqdata_long2,
             aes(x=collection_date, y=100*prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("number\nsequenced", trans="sqrt",
                        range=c(1, 4), limits=c(1,10^3), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_multinom_501YV1_501YV2_501YV3_response

saveRDS(plot_multinom_501YV1_501YV2_501YV3_response, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2_501YV3_multinomial fit_model preds_response.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pdf"), width=8, height=6)


# on logit scale:

be_seq_mfit0_preds3 = be_seq_mfit0_preds2
ymin = 0.001
ymax = 0.900001
be_seq_mfit0_preds3$asymp.LCL[be_seq_mfit0_preds3$asymp.LCL<ymin] = ymin
be_seq_mfit0_preds3$asymp.UCL[be_seq_mfit0_preds3$asymp.UCL<ymin] = ymin
be_seq_mfit0_preds3$asymp.UCL[be_seq_mfit0_preds3$asymp.UCL>ymax] = ymax
be_seq_mfit0_preds3$prob[be_seq_mfit0_preds3$prob<ymin] = ymin

plot_multinom_501YV1_501YV2_501YV3 = qplot(data=be_seq_mfit0_preds3, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share of diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("blue","red","green3","black")) +
  scale_colour_manual("variant", values=c("blue","red","green3","black")) +
  geom_point(data=be_basseqdata_long2,
             aes(x=collection_date, y=prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(1, 4), limits=c(10,10^3), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-03-01")), ylim=c(0.001, 0.9000001), expand=c(0,0))
plot_multinom_501YV1_501YV2_501YV3


saveRDS(plot_multinom_501YV1_501YV2_501YV3, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.pdf"), width=8, height=6)


# estimated proportion of 501Y.V1 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18680 0.6153158 0.06310555 NA 0.4916312 0.7390004      2021-02-22

# estimated proportion of 501Y.V1 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18687 0.6551634 0.1519726 NA 0.3573025 0.9530242      2021-03-01

# estimated proportion of 501Y.V2 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18680 0.0840299 0.02993692 NA 0.02535461 0.1427052      2021-02-22

# estimated proportion of 501Y.V2 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18687 0.08573114 0.04243251 NA 0.002564944 0.1688973      2021-03-01

# estimated proportion of 501Y.V3 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18680 0.0529844 0.07281237 NA -0.08972521  0.195694      2021-02-22

# estimated proportion of 501Y.V3 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18687 0.1115476 0.1950751 NA -0.2707925 0.4938877      2021-03-01



# estimated proportion of one of the three VOCs among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18680 0.7429313 0.03203194 NA 0.6801498 0.8057128      2021-02-22

# estimated proportion of one of the three VOCs among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18687 0.8375905 0.02882499 NA 0.7810946 0.8940864      2021-03-01




# the time at which new lab diagnoses would be by more than 50%, 75% 90% by one of the three VOCs :
be_seq_mfit0_preds2_subs = be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date>as.Date("2021-01-01")&
                                                 be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
# >50% by 10th of February [8th Feb - 12th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]

# >75% by 23d of February [19th Feb - 27th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]

# >90% by 8th of March [2nd March - 14th of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]


# the time at which new infections would be by more than 50%, 75% 90% by one of the three VOCs
# (counting 7 days between infection & diagnosis) :
# >50% by 3d of Feb [1st Feb - 5th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]-7

# >75% by 16th of February [12th Feb - 20th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]-7

# >90% by 1st of March [23 Febr - 7th of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]-7



# 2. ESTIMATE PROPORTION OF S DROPOUT SAMPLES THAT ARE B.1.1.7 / 501Y.V1 IN BELGIUM IN FUNCTION OF TIME BASED ON SEQUENCING DATA ####

datBE_b117 = read.csv(paste0(".//data//", dat, "//sequencing_Sdropouts.csv"), check.names=F) # n_b117/n_sgtf_seq = prop of S dropout samples that were B.1.1.7
datBE_b117$country = "Belgium"
datBE_b117$collection_date = as.Date(datBE_b117$collection_date)
datBE_b117$collection_date_num = as.numeric(datBE_b117$collection_date)
datBE_b117$propB117 = datBE_b117$n_b117/datBE_b117$n_sgtf_seq
datBE_b117$obs = factor(1:nrow(datBE_b117))
datBE_b117

write.csv(datBE_b117, file=".\\data\\be_latest\\sequencing_Sdropouts.csv", row.names=FALSE)

fit_seq = glmer(cbind(n_b117,n_sgtf_seq-n_b117) ~ (1|obs)+scale(collection_date_num), family=binomial(logit), data=datBE_b117)
summary(fit_seq)

# implied growth rate advantage of 501Y.V1 over other earlier strains showing S dropout:
as.data.frame(emtrends(fit_seq, ~ 1, var="collection_date_num"))[,c(2,5,6)]
#   collection_date_num.trend  asymp.LCL asymp.UCL
# 1             0.1107948 0.07844251  0.143147

# with a generation time of 4.7 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of
exp(4.7*as.data.frame(emtrends(fit_seq, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1              1.683265  1.445825  1.959699

plot(fit_seq)

# PLOT MODEL FIT
extrapolate = 90 # nr of days to extrapolate fit into the future
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                   at=list(collection_date_num=seq(as.numeric(min(datBE_b117$collection_date)),
                                                               as.numeric(max(datBE_b117$collection_date))+extrapolate)), 
                                   type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")

# prop of S dropout samples among newly diagnosed infections that are now estimated to be B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$collection_date==today,]
#    collection_date_num      prob         SE  df asymp.LCL asymp.UCL collection_date
# 62           18680 0.997214 0.002161191 Inf 0.9873354 0.9993923      2021-02-22

# prop of S dropout samples among new infections that are now estimated to be B.1.1.7 / 501Y.V1 (using 7 days for time from infection to diagnosis)
fitseq_preds[fitseq_preds$collection_date==(today+7),]
#    collection_date_num     prob          SE  df asymp.LCL asymp.UCL collection_date
# 69           18687 0.9987151 0.001143097 Inf  0.992683 0.9997756      2021-03-01

# from 13th of Jan 2021 >80% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.80,"collection_date"][1]

# from 20th of Jan 2021 >90% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.90,"collection_date"][1]

# from 11th of Feb 2021 >99% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.99,"collection_date"][1]

# on logit scale:
plot_fitseq = qplot(data=fitseq_preds, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~laboratory) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), fill=I("#b0c4de"), alpha=I(1)) +
  geom_line(aes(y=prob), colour=I("steelblue"), alpha=I(1)) +
  ylab("S dropout samples that are 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("REPRESENTATION OF 501Y.V1 AMONG S DROPOUT SAMPLES") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-12-01"),today), 
    ylim=c(0.01,0.999002), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=datBE_b117, 
             aes(x=collection_date, y=propB117, size=n_sgtf_seq), colour=I("steelblue"), alpha=I(1)) +
  scale_size_continuous("number of S dropout\nsamples sequenced", trans="sqrt", 
                        range=c(1, 6), limits=c(1,
                                                   10^(round(log10(max(datBE_b117$n_sgtf_seq)),0)+1) ), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right") + xlab("Collection date")
plot_fitseq

# saveRDS(plot_fitseq, file = paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_link scale.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_link scale.pptx"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_link scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_link scale.pdf"), width=8, height=6)


# same on response scale:
plot_fitseq_response = qplot(data=fitseq_preds, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~laboratory) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL), fill=I("#b0c4de"), alpha=I(1)) +
  geom_line(aes(y=100*prob), colour=I("steelblue"), alpha=I(1)) +
  ylab("S dropout samples that are 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("REPRESENTATION OF 501Y.V1 AMONG S DROPOUT SAMPLES") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-12-01"),today), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=datBE_b117, 
             aes(x=collection_date, y=100*propB117, size=n_sgtf_seq), colour=I("steelblue"), alpha=I(1)) +
  scale_size_continuous("number of S dropout\nsamples sequenced", trans="sqrt", 
                        range=c(1, 6), limits=c(1,10^(round(log10(max(datBE_b117$n_sgtf_seq)),0)+1) ), 
                        breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right") + xlab("Collection date")
plot_fitseq_response

saveRDS(plot_fitseq, file = paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.pdf"), width=8, height=6)




# 3. ANALYSIS OF Ct VALUES OF S-DROPOUT & NON-S-DROPOUT SAMPLES IN BELGIUM ####

# Read in Ct data of all valid PCRs
file_dec1 = paste0(".//data//", dat, "//PCR December 2020 1 to 20.xlsx")
file_dec2 = paste0(".//data//", dat, "//PCR December 2020 21 to 31.xlsx")
file_jan = paste0(".//data//", dat, "//PCR January 2021 complete.xlsx")
file_feb1 = paste0(".//data//", dat, "//PCR February 2021 1 to 9.xlsx")
file_feb2 = paste0(".//data//", dat, "//PCR February 2021 10 to 21 Feb.xlsx")
sheets = excel_sheets(file_jan)
ctdata_dec1 = map_df(sheets, ~ read_excel(file_dec1, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_dec1$Analysis_date), origin="1899-12-30")) # "2020-12-01" "2020-12-20"
ctdata_dec2 = map_df(sheets, ~ read_excel(file_dec2, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_dec2$Analysis_date), origin="1899-12-30")) # "2020-12-21" "2020-12-31"
ctdata_jan = map_df(sheets, ~ read_excel(file_jan, sheet = .x, skip = 1, 
                                       col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                       col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_jan$Analysis_date), origin="1899-12-30")) # "2021-01-01" "2021-01-31"
ctdata_feb1 = map_df(sheets, ~ read_excel(file_feb1, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_feb1$Analysis_date), origin="1899-12-30")) # "2021-01-01" "2021-01-08"
ctdata_feb2 = map_df(sheets, ~ read_excel(file_feb2, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_feb2$Analysis_date), origin="1899-12-30")) # "2021-01-30" "2021-02-21"
# PS there is still something wrong with file ctdata_feb2 (not the right date range & negatives missing) & ctdata_feb1 misses data 9th of Febr
ctdata = bind_rows(ctdata_dec1, ctdata_dec2, ctdata_jan, ctdata_feb1, ctdata_feb2)
range(as.Date(as.numeric(ctdata$Analysis_date), origin="1899-12-30")) # "2020-12-01" "2021-02-21"
ctdata$Laboratory[ctdata$Laboratory=="ULG - FF 3.x"] = "ULG"
unique(ctdata$Laboratory) 
# "ULG"     "Namur"            "UMons - Jolimont" "UZ leuven"        "UZA"              "UZ Gent"          "ULB"             
# "Saint LUC - UCL"
unique(ctdata$Outcome) 
unique(ctdata_dec1$Outcome) 
unique(ctdata_dec2$Outcome) 
unique(ctdata_jan$Outcome) 
unique(ctdata_feb1$Outcome) 
unique(ctdata_feb2$Outcome) # all labs now use "Detected" for "Detected" or "Positive"
ctdata$Outcome[ctdata$Outcome=="Detected"] = "Positive"
ctdata$Outcome[ctdata$Outcome=="Not detected"] = "Negative"
unique(ctdata$Outcome) # "Positive" "Negative"
ctdata$Analysis_date = as.Date(as.numeric(ctdata$Analysis_date), origin="1899-12-30")
sum(is.na(ctdata$Analysis_date)) # 0
range(ctdata$Analysis_date) # "2020-12-01" - "2021-02-21"
ctdata$collection_date = ctdata$Analysis_date-1 # collection date = analysis date-1 
sum(is.na(ctdata$collection_date)) # 0
ctdata$collection_date_num = as.numeric(ctdata$collection_date)
range(ctdata$collection_date) # "2020-11-30" "2021-02-20"
ctdata$group = interaction(ctdata$Outcome, ctdata$S_dropout)
ctdata$group = droplevels(ctdata$group)
unique(ctdata$group) # Positive.0 Positive.1 Negative.0
ctdata$Outcome = factor(ctdata$Outcome)
unique(ctdata$Outcome) # Positive Negative
ctdata$Laboratory = factor(ctdata$Laboratory)
ctdata$S_dropout = factor(ctdata$S_dropout)
head(ctdata)
str(ctdata)
nrow(ctdata) # 716677

ctdata_onlypos = ctdata[ctdata$Outcome=="Positive",] # subset with only the positive samples
ctdata_onlypos = bind_rows(ctdata_onlypos[ctdata_onlypos$S_dropout=="0",], ctdata_onlypos[ctdata_onlypos$S_dropout=="1",])

# ANALYSIS OF Ct VALUES OF S DROPOUT & NON-S DROPOUT SAMPLES

# plot & analysis of Ct values of all labs for dates from 13th of Jan onward when >80% of all S dropouts were B.1.1.7 / 501Y.V1

(fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1]) # "2021-01-13", from 13th of Jan >80% of all S dropouts are B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections
# this is the same criterion that was used for the SGTF analysis in the UK (N. Davies, pers. comm.)
subs = (ctdata_onlypos$collection_date > (fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1])) & 
       (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30) 
ctdata_onlypos_subs = ctdata_onlypos[subs,]
ctdata_onlypos_subs = ctdata_onlypos_subs[!(is.na(ctdata_onlypos_subs$S_dropout)|
                                              is.na(ctdata_onlypos_subs$N_cq)|
                                              is.na(ctdata_onlypos_subs$ORF1_cq)|
                                              (ctdata_onlypos_subs$ORF1_cq==0)),]

cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.76

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
# Laboratory correlation_Ct_N_ORF1ab
# Namur                       Namur             0.947866198
# Saint LUC - UCL   Saint LUC - UCL             0.981080423
# ULB                           ULB             0.981429368
# ULG                           ULG             0.955847687
# UMons - Jolimont UMons - Jolimont            -0.282023995
# UZ Gent                   UZ Gent            -0.003990932
# UZ leuven               UZ leuven             0.976912170
# UZA                           UZA             0.881584309

ctcorplot_all_labs = qplot(data=ctdata_onlypos_subs, x=ORF1_cq, y=N_cq, group=Laboratory, fill=S_dropout, colour=S_dropout, size=I(3), shape=I(16)) +
  facet_wrap(~Laboratory) + 
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Ct value ORF1ab gene") + ylab("Ct value N gene") +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctcorplot_all_labs
# PS weird results for UMons & UZ Gent, and to some extent UZA, not sure of the cause
saveRDS(ctcorplot_all_labs, file = paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.pptx"), width=7, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.png"), width=7, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.pdf"), width=7, height=6)


ctplot_rawdataN_all_labs = qplot(data=ctdata_onlypos_subs, x=collection_date, y=N_cq, group=S_dropout, 
                        colour=S_dropout, fill=S_dropout, geom="point", size=I(1), shape=I(16)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "N gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataN_all_labs
# PS UZA has suspect Ct values before the 21st of Jan & UMons & UGhent also have very different ranges for the Ct values
ctdata_onlypos_subs[ctdata_onlypos_subs$Laboratory=="UZA"&ctdata_onlypos_subs$N_cq<=15,"collection_date"][1,] # "2021-01-21"
ctplot_rawdataORF1_all_labs = qplot(data=ctdata_onlypos_subs, x=collection_date, y=ORF1_cq, group=S_dropout, 
                                 colour=S_dropout, fill=S_dropout, geom="point", size=I(1), shape=I(16)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "ORF1ab gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataORF1_all_labs

ctplots_rawdata_all_labs = ggarrange(ctplot_rawdataN_all_labs+xlab("")+theme(axis.text.x = element_blank()), 
                                     ctplot_rawdataORF1_all_labs,
                                     ncol=1, common.legend=TRUE, legend="right")
ctplots_rawdata_all_labs

saveRDS(ctplots_rawdata_all_labs, file = paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.pptx"), width=7, height=9)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.png"), width=7, height=9)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.pdf"), width=7, height=9)



# plot & analysis of Ct values of all labs except Mons (low S dropout counts & weird Ct value distribution), 
# ULG (low sample size), UZ Gent (heavily involved in active surveillance, so not representative, plus weird Ct value distribution) &
# UZA (heavily involved in active surveillance, so not representative, plus weird Ct value distribution)
# for dates from 13th of Jan onward when >80% of all S dropouts were B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections

# PS we could still consider including the data from UZA since the 21st of Jan, as Ct values 
# from then onwards are in the same range as the other labs

labs_to_remove = c("UMons - Jolimont", "ULG", "UZ Gent", "UZA")
sel_labs = setdiff(unique(ctdata_onlypos_subs$Laboratory), labs_to_remove) 
sel_labs # "UZ leuven"       "Saint LUC - UCL" "ULB"             "Namur" 
# we use data from these 4 labs as the data distribution was comparable for these
# they also had large sample size & were not heavily involved in active surveillance
# sel_labs = unique(ctdata_onlypos$Laboratory) # to select data from all the labs, but distribution not comparable for all
# we use the subset of timepoints (from 13th Jan 2021 onwards) where >80% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
(fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1]) # "2021-01-13", from 13th of Jan >80% of all S dropouts are B.1.1.7 / 501Y.V1
# we also just use the positive samples with relatively strong signal, (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30)
# to not include pos samples with very low viral titers (indicative of old infections etc)
# this is the same criterion that was used for the SGTF analysis in the UK (N. Davies, pers. comm.)
subs = (ctdata_onlypos$collection_date > (fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1])) &
  (ctdata_onlypos$Laboratory %in% sel_labs) & (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30)
ctdata_onlypos_subs = ctdata_onlypos[subs,]
ctdata_onlypos_subs = ctdata_onlypos_subs[!(is.na(ctdata_onlypos_subs$S_dropout)|
                                              is.na(ctdata_onlypos_subs$N_cq)|
                                              is.na(ctdata_onlypos_subs$ORF1_cq)|
                                              (ctdata_onlypos_subs$ORF1_cq==0)),]
ctdata_onlypos_subs$Laboratory = droplevels(ctdata_onlypos_subs$Laboratory)
  
# make joint dataset for integrated analysis of both genes to estimate average effect across both sets of genes
ctdata_onlypos_subs_bothgenes = rbind(data.frame(ctdata_onlypos_subs, Gene="N gene", Ct=ctdata_onlypos_subs$N_cq), 
                                      data.frame(ctdata_onlypos_subs, Gene="ORF1ab gene", Ct=ctdata_onlypos_subs$ORF1_cq))
# we define a high viral load as one where the Ct value was 1.25x lower than in the non-S dropout sample group
# which was a Ct value < 15.099 for the N gene and < 16.001 for the ORF1ab gene
thresh_N = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","N_cq"]))/1.25
thresh_N # 15.099
thresh_ORF1 = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","ORF1_cq"]))/1.25
thresh_ORF1 # 16.001
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="N gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="N gene"]<thresh_N
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"]<thresh_ORF1


# check correlation between Ct values for N & ORF1ab gene
cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.97, t=549.68, p<2E-16

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                      Laboratory correlation_Ct_N_ORF1ab
# Namur                     Namur               0.9478662
# Saint LUC - UCL Saint LUC - UCL               0.9810804
# ULB                         ULB               0.9814294
# UZ leuven             UZ leuven               0.9769122


ctcorplot_sellabs = qplot(data=ctdata_onlypos_subs, x=ORF1_cq, y=N_cq, group=Laboratory, fill=S_dropout, colour=S_dropout, size=I(3), shape=I(16)) +
  facet_wrap(~Laboratory) + 
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Ct value ORF1ab gene") + ylab("Ct value N gene") +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept=thresh_N, colour=alpha("black", 1), lwd=I(0.3), lty=I(2)) +
  geom_vline(xintercept=thresh_ORF1, colour=alpha("black", 1), lwd=I(0.3), lty=I(2))
ctcorplot_sellabs
saveRDS(ctcorplot_sellabs, file = paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.pptx"), width=7, height=4.5)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.png"), width=7, height=4.5)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.pdf"), width=7, height=4.5)


ctplot_rawdataN = qplot(data=ctdata_onlypos_subs, x=collection_date, y=N_cq, group=S_dropout, 
                        colour=S_dropout, fill=S_dropout, geom="point", shape=I(16), size=I(2)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "N gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataN # there is no obvious temporal patterns

ctplot_rawdataORF1 = qplot(data=ctdata_onlypos_subs, x=collection_date, y=ORF1_cq, group=S_dropout, 
                           colour=S_dropout, fill=S_dropout, geom="point", shape=I(16), size=I(2)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "ORF1ab gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataORF1 # there is no obvious temporal patterns

ctplots_rawdata = ggarrange(ctplot_rawdataN+xlab("")+theme(axis.text.x = element_blank()), 
                            ctplot_rawdataORF1,
                            ncol=1, common.legend=TRUE, legend="right")
ctplots_rawdata # there is no obvious temporal patterns

saveRDS(ctplots_rawdata, file = paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.pptx"), width=6, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.png"), width=6, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.pdf"), width=6, height=8)




# quantile/median regression to compare median Ct values of both genes across S dropout & non-S dropout samples in the different labs
qr_bothgenes0 = rq(Ct ~ Gene + S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes1 = rq(Ct ~ Gene * S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes2 = rq(Ct ~ Gene + S_dropout * Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes3 = rq(Ct ~ Gene * Laboratory + S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes4 = rq(Ct ~ (Gene + Laboratory + S_dropout)^2, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes5 = rq(Ct ~ Gene * Laboratory * S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
AIC(qr_bothgenes0, k=-1) # 214755.8
AIC(qr_bothgenes1, k=-1) # 214734
AIC(qr_bothgenes2, k=-1) # 214572.4 
AIC(qr_bothgenes3, k=-1) # 214767.4
AIC(qr_bothgenes4, k=-1) # 214556.7 # fits data best based on BIC criterion (lowest value, PS: here AIC with k<0 returns BIC)
AIC(qr_bothgenes5, k=-1) # 214585.4

summary(qr_bothgenes4)
# tau: [1] 0.5
# 
# Coefficients:
#   Value     Std. Error t value   Pr(>|t|) 
# (Intercept)             18.30352   0.05320  344.02131   0.00000
# Gene1                   -0.78890   0.05221  -15.11005   0.00000
# Laboratory1             -0.50165   0.09395   -5.33940   0.00000
# Laboratory2             -0.61455   0.09444   -6.50718   0.00000
# Laboratory3             -0.64705   0.08850   -7.31106   0.00000
# S_dropout1               1.22808   0.05277   23.27324   0.00000
# Gene1:Laboratory1        0.16578   0.09228    1.79652   0.07242
# Gene1:Laboratory2        0.07128   0.09594    0.74288   0.45756
# Gene1:Laboratory3       -0.18283   0.08814   -2.07426   0.03806
# Gene1:S_dropout1         0.23072   0.05343    4.31848   0.00002
# Laboratory1:S_dropout1  -0.03735   0.09378   -0.39828   0.69043
# Laboratory2:S_dropout1  -0.15855   0.09451   -1.67757   0.09344
# Laboratory3:S_dropout1  -0.72295   0.08858   -8.16187   0.00000

qr_emmeans_bylab99 = data.frame(emmeans(qr_bothgenes4, ~ Laboratory + Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans_bylab99
qr_emmeans99 = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans99
qr_emmeans = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.95)) # median Ct values + 95% CLs
qr_emmeans
#          Gene S_dropout   emmean        SE    df lower.CL upper.CL
# 1      N gene         0 18.97342 0.0896920 32327 18.79763 19.14922
# 2 ORF1ab gene         0 20.08977 0.1005441 32327 19.89270 20.28685
# 3      N gene         1 16.05582 0.1206785 32327 15.81929 16.29236
# 4 ORF1ab gene         1 18.09507 0.1098416 32327 17.87978 18.31037

# mean difference in median Ct value of 2.46, which is highly significant across both genes: p<0.0001
contrast(emmeans(qr_bothgenes4, ~ S_dropout, level=0.95), method="pairwise") 
# contrast estimate    SE    df t.ratio p.value
# 0 - 1        2.46 0.106 32327 23.273  <.0001
# mean difference in median Ct value of 2.92 of rN gene and 1.99 for ORF1ab gene
confint(contrast(emmeans(qr_bothgenes4, ~ S_dropout|Gene, level=0.95), method="pairwise"))
# Gene = N gene:
#   contrast estimate   SE    df lower.CL upper.CL
# 0 - 1        2.92 0.15 32327     2.62     3.21
# 
# Gene = ORF1ab gene:
#   contrast estimate   SE    df lower.CL upper.CL
# 0 - 1        1.99 0.15 32327     1.70     2.29
# 
# Results are averaged over the levels of: Laboratory 
# Confidence level used: 0.95 


# violin plots by gene & lab & S dropout with expected marginal means+99% CLs of best fitting median regression model
ctviolinplots_bylab = ggplot(data=ctdata_onlypos_subs_bothgenes, aes(x=factor(S_dropout), y=Ct, fill=factor(S_dropout))) +
  geom_violin(alpha=1, colour=NA, trim=TRUE, draw_quantiles=TRUE, adjust=2, scale="width") +
  geom_crossbar(data=qr_emmeans_bylab99, aes(x=factor(S_dropout), y=emmean, ymin=lower.CL, ymax=upper.CL, group=Gene)) +
  # stat_summary(fun.data=data_summary,  
  #             geom="pointrange", aes(color=factor(S_dropout))) +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  # geom_point(aes(colour=factor(S_dropout))) +
  facet_wrap(~ Gene + Laboratory,ncol=4) +
  scale_colour_manual("", values=alpha(c("steelblue","lightcoral"), 1), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("steelblue","lightcoral"), 1), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("") + ylab("Ct value") + 
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("0","1"), labels=c("S pos","S dropout")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctviolinplots_bylab

saveRDS(ctviolinplots_bylab, file = paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.pptx"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.png"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.pdf"), width=6, height=6)

# violin plots by gene & S dropout with expected marginal means+99% CLs of best fitting median regression model
ctviolinplots = ggplot(data=ctdata_onlypos_subs_bothgenes, aes(x=factor(S_dropout), y=Ct, fill=factor(S_dropout))) +
  geom_violin(alpha=1, colour=NA, trim=TRUE, draw_quantiles=TRUE, adjust=2, scale="width") +
  geom_crossbar(data=qr_emmeans99, aes(x=factor(S_dropout), y=emmean, ymin=lower.CL, ymax=upper.CL, group=Gene, lwd=I(0.1))) +
  # stat_summary(fun.data=data_summary,  
  #             geom="pointrange", aes(color=factor(S_dropout))) +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  # geom_point(aes(colour=factor(S_dropout))) +
  facet_wrap(~Gene,ncol=4) +
  scale_colour_manual("", values=muted(c("steelblue","lightcoral"), l=55), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=muted(c("steelblue","lightcoral"), l=55), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("") + ylab("Ct value") + 
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("0","1"), labels=c("S pos","S dropout")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctviolinplots

saveRDS(ctviolinplots, file = paste0(".\\plots\\",dat,"\\Ct values_violin plots.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.pptx"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.png"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.pdf"), width=6, height=6)




# binomial GLMMs to test for differences in prop with high viral load (Ct values a factor of 1.25 lower than median Ct in non-S dropout samples) :

fitct_highvirload_0A = glmer(high_viral_load ~ (1|Laboratory) + Gene + S_dropout, family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_1A = glmer(high_viral_load ~ (1|Laboratory) + Gene + S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_2A = glmer(high_viral_load ~ (1|Laboratory) + Gene + S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_3A = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene + S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_4A = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene + S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_0B = glmer(high_viral_load ~ (1|Laboratory) + Gene * S_dropout, family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_1B = glmer(high_viral_load ~ (1|Laboratory) + Gene * S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_2B = glmer(high_viral_load ~ (1|Laboratory) + Gene * S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_3B = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene * S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_4B = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene * S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)

BIC(fitct_highvirload_0A, fitct_highvirload_1A, fitct_highvirload_2A, fitct_highvirload_3A, fitct_highvirload_4A,
    fitct_highvirload_0B, fitct_highvirload_1B, fitct_highvirload_2B, fitct_highvirload_3B, fitct_highvirload_4B)
# df      BIC
# fitct_highvirload_0A  4 41190.92
# fitct_highvirload_1A  5 41192.97
# fitct_highvirload_2A  6 41202.62
# fitct_highvirload_3A  6 41203.28
# fitct_highvirload_4A  7 41212.94
# fitct_highvirload_0B  5 41167.42
# fitct_highvirload_1B  6 41169.44
# fitct_highvirload_2B  9 41199.54
# fitct_highvirload_3B  7 41179.76
# fitct_highvirload_4B 10 41209.86


# fitct_highvirload_0B best model
summary(fitct_highvirload_0B) # S dropout samples more frequently have high viral load based on N gene Ct values
# Random effects:
#   Groups     Name        Variance Std.Dev.
# Laboratory (Intercept) 0.03284  0.1812  
# Number of obs: 32340, groups:  Laboratory, 4
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.59899    0.09152  -6.545 5.95e-11 ***
#   Gene1             0.12286    0.01242   9.892  < 2e-16 ***
#   S_dropout1       -0.16234    0.01249 -12.993  < 2e-16 ***
#   Gene1:S_dropout1 -0.07224    0.01242  -5.816 6.02e-09 ***

plot(allEffects(fitct_highvirload_0B))


# odds to encounter high viral load samples based on Ct values of both genes (high vir load = Ct values >1.25x lower than median in non-S dropout samples)
# 1.38x [1.32-1.45x] 95% CLs increased among S dropout samples
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout, type="response"), method="revpairwise", type="response"))
# contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.38 0.0346 Inf      1.32      1.45

# odds ratio = 1.6 [1.49-1.71] for N gene & 1.2 [1.12-1.28] for ORF1ab gene 
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout|Gene, type="response"), method="revpairwise", type="response"))
# Gene = N gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.6 0.0554 Inf      1.49      1.71
# 
# Gene = ORF1ab gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.2 0.0429 Inf      1.12      1.28

fitct_highvirload_emmeans = as.data.frame(emmeans(fitct_highvirload_0B, ~ S_dropout+Gene, type="response"))
fitct_highvirload_emmeans$S_dropout = factor(fitct_highvirload_emmeans$S_dropout)
fitct_highvirload_emmeans$Gene = factor(fitct_highvirload_emmeans$Gene)
fitct_highvirload_emmeans
# S_dropout        Gene      prob         SE  df asymp.LCL asymp.UCL
# 1         0      N gene 0.3294412 0.02052853 Inf 0.2905225 0.3708482
# 2         1      N gene 0.4399057 0.02338181 Inf 0.3947135 0.4861166
# 3         0 ORF1ab gene 0.3074735 0.01980499 Inf 0.2700726 0.3475874
# 4         1 ORF1ab gene 0.3471130 0.02158259 Inf 0.3061012 0.3905270

ct_highvirload_binGLMM = ggplot(data=fitct_highvirload_emmeans, 
                                               aes(x=S_dropout, y=prob*100, fill=S_dropout, group=S_dropout)) +
  facet_wrap(~ Gene) +
  geom_col(colour=NA, position=position_dodge2(width=0.8, padding=0.5)) +
  geom_linerange(aes(ymin=asymp.LCL*100, ymax=asymp.UCL*100), position=position_dodge2(width=0.8, padding=0.5)) +
  scale_fill_manual("", values=muted(c("steelblue","lightcoral"), l=55), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_y_continuous(breaks=seq(0,100,by=10), expand=c(0,0)) +
  ylab("High viral load samples (%)") + xlab("Gene") + coord_cartesian(ylim=c(0,50)) +
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("0","1"), labels=c("S pos","S dropout"), expand=c(0.3,0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ct_highvirload_binGLMM

ctplots_multipanel = ggarrange(ctviolinplots+xlab(""), 
                               ct_highvirload_binGLMM,
                                ncol=1, legend=NULL, common.legend=FALSE)
ctplots_multipanel

saveRDS(ctplots_all_rawdata, file = paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.pdf"), width=8, height=6)




# 4. ESTIMATE GROWTH RATE AND TRANSMISSION ADVANTAGE OF B.1.1.7 / 501Y.V1 IN BELGIUM BASED ON S-GENE TARGET FAILURE DATA ####

# we remove ULG - FF 3.x due to low sample size & also just use positive samples with Ct values for N & ORF1ab < 30 to focus on active, recent infections
# we also remove UZ Gent data because this lab was very heavily involved in active surveillance, and so 
# its results could bias the overall inferred growth rate advantage (plus this lab also had weird Ct patterns)
sel_labs = setdiff(unique(ctdata$Laboratory), c("ULG - FF 3.x", "ULG", "UZ Gent"))  # unique(ctdata$Laboratory) 
# setdiff(unique(ctdata$Laboratory), c("UMons - Jolimont", "ULG - FF 3.x", "UZ Gent", "UZA")) 
# setdiff(unique(ctdata$Laboratory), c("UMons - Jolimont", "UZ Gent","UZA","ULG - FF 3.x")) 
sel_labs  
pos_ctbelow30 = (ctdata$Laboratory %in% sel_labs) & (((ctdata$Outcome=="Positive")&(ctdata$N_cq<30)&(ctdata$ORF1_cq<30)))
pos_ctbelow30 = pos_ctbelow30[!is.na(pos_ctbelow30)]
pos_ctabove30 = (ctdata$Laboratory %in% sel_labs) & (((ctdata$Outcome=="Positive")&(ctdata$N_cq>=30)|(ctdata$ORF1_cq>=30)))
pos_ctabove30 = pos_ctabove30[!is.na(pos_ctabove30)]
100*sum(pos_ctabove30) / (sum(pos_ctabove30)+sum(pos_ctbelow30)) # 25% of positives have Ct > 30

subs = which((ctdata$Laboratory %in% sel_labs) & (((ctdata$Outcome=="Positive")&(ctdata$N_cq<30)&(ctdata$ORF1_cq<30))|
                                            (ctdata$Outcome=="Negative")))
ctdata_subs = ctdata
ctdata_subs = ctdata_subs[subs,]
nrow(ctdata) # 716677
nrow(ctdata_subs) # 556431

ctdata_subs$group = factor(ctdata_subs$group, 
                           levels=c("Negative.0","Positive.0","Positive.1"), 
                           labels=c("negative","S_pos","S_dropout"))
ctdata_subs$Laboratory = droplevels(ctdata_subs$Laboratory)

# aggregated counts by date (sample date) and Laboratory
data_ag = as.data.frame(table(ctdata_subs$collection_date, ctdata_subs$Laboratory, ctdata_subs$group), check.names=F)
colnames(data_ag) = c("collection_date", "LABORATORY", "GROUP", "COUNT")
data_ag_wide = spread(data_ag, GROUP, COUNT)
colnames(data_ag_wide)[colnames(data_ag_wide) %in% c("negative","S_pos","S_dropout")] = c("n_neg","n_spos","n_sgtf")
data_ag_wide$n_pos = data_ag_wide$n_spos+data_ag_wide$n_sgtf
data_ag_wide$total = data_ag_wide$n_neg + data_ag_wide$n_pos
data_ag_wide$collection_date = as.Date(data_ag_wide$collection_date)
data_ag_wide$collection_date_num = as.numeric(data_ag_wide$collection_date)
# calculate prop of S dropout that is actually B.1.1.7 / 501Y.V1 estimated from binomial GLMM:
# (using expected marginal mean calculated using emmeans, taking into account random effects)
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                     at=list(collection_date_num=seq(min(data_ag_wide$collection_date_num),
                                                                     max(data_ag_wide$collection_date_num))),
                                     type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")
# prob that S dropout was B.1.1.7 / 501Y.V1
data_ag_wide$TRUEPOS = fitseq_preds$prob[match(data_ag_wide$collection_date, fitseq_preds$collection_date)] 
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate
data_ag_wide$est_n_B117 = data_ag_wide$n_sgtf * data_ag_wide$TRUEPOS 
## estimates props 501Y.V2 & 501Y.V3 from multinomial fit, we use this to also estimate the nr of wild type pos samples (excluding either of the 3 VOCs)
#data_ag_wide$prop501Y.V2 = be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V2 (South African)","prob"][match(data_ag_wide$collection_date,be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V2 (South African)","collection_date"])]
#data_ag_wide$prop501Y.V3 = be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V3 (Brazilian)","prob"][match(data_ag_wide$collection_date,be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V3 (Brazilian)","collection_date"])]
#data_ag_wide$est_npos_wildtype = (data_ag_wide$n_pos - data_ag_wide$est_n_B117)*(1-data_ag_wide$prop501Y.V2)
#data_ag_wide$est_npos_wildtype[data_ag_wide$est_n_B117>data_ag_wide$est_npos_wildtype] = data_ag_wide$est_n_B117[data_ag_wide$est_n_B117>data_ag_wide$est_npos_wildtype]
data_ag_wide$propB117 = data_ag_wide$est_n_B117 / data_ag_wide$n_pos
#data_ag_wide$propB117amongwildtype = data_ag_wide$est_n_B117 / data_ag_wide$est_npos_wildtype
#data_ag_wide$propB117amongwildtype[data_ag_wide$propB117amongwildtype>1] = 1
data_ag_wide$obs = factor(1:nrow(data_ag_wide))
data_ag_wide = data_ag_wide[data_ag_wide$total != 0, ]
head(data_ag_wide)
tail(data_ag_wide, 70)

write.csv(data_ag_wide, file=".\\data\\be_latest\\be_B117_by lab.csv", row.names=FALSE)


# aggregated counts by date over all Laboratories
data_ag_byday = as.data.frame(table(ctdata_subs$collection_date, ctdata_subs$group), check.names=F)
colnames(data_ag_byday) = c("collection_date", "GROUP", "COUNT")
data_ag_byday_wide = spread(data_ag_byday, GROUP, COUNT)
colnames(data_ag_byday_wide)[colnames(data_ag_byday_wide) %in% c("negative","S_pos","S_dropout")] = c("n_neg","n_spos","n_sgtf")
data_ag_byday_wide$n_pos = data_ag_byday_wide$n_spos+data_ag_byday_wide$n_sgtf
data_ag_byday_wide$total = data_ag_byday_wide$n_neg + data_ag_byday_wide$n_pos
data_ag_byday_wide$collection_date = as.Date(data_ag_byday_wide$collection_date)
data_ag_byday_wide$collection_date_num = as.numeric(data_ag_byday_wide$collection_date)
# calculate prop of S dropout that is actually B.1.1.7 / 501Y.V1 estimated from binomial GLMM:
# (using expected marginal mean calculated using emmeans, taking into account random effects)
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                     at=list(collection_date_num=seq(min(data_ag_byday_wide$collection_date_num),
                                                                     max(data_ag_byday_wide$collection_date_num))),
                                     type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")
# prob that S dropout was B.1.1.7 / 501Y.V1
data_ag_byday_wide$TRUEPOS = fitseq_preds$prob[match(data_ag_byday_wide$collection_date, fitseq_preds$collection_date)] 
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate
data_ag_byday_wide$est_n_B117 = data_ag_byday_wide$n_sgtf * data_ag_byday_wide$TRUEPOS
## est prop 501Y.V2 & V3 among positive tests from multinomial fit
#data_ag_byday_wide$prop501Y.V2 = be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V2 (South African)","prob"][match(data_ag_byday_wide$collection_date,be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V2 (South African)","collection_date"])]
#data_ag_byday_wide$prop501Y.V3 = be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V3 (Brazilian)","prob"][match(data_ag_byday_wide$collection_date,be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="501Y.V3 (Brazilian)","collection_date"])]
#data_ag_byday_wide$est_npos_wildtype = (data_ag_byday_wide$n_pos - data_ag_byday_wide$est_n_B117)*(1-data_ag_byday_wide$prop501Y.V2)
#data_ag_byday_wide$est_npos_wildtype[data_ag_byday_wide$est_n_B117>data_ag_byday_wide$est_npos_wildtype] = data_ag_byday_wide$est_n_B117[data_ag_byday_wide$est_n_B117>data_ag_byday_wide$est_npos_wildtype]
data_ag_byday_wide$propB117 = data_ag_byday_wide$est_n_B117 / data_ag_byday_wide$n_pos
#data_ag_byday_wide$propB117amongwildtype = data_ag_byday_wide$est_n_B117 / data_ag_byday_wide$est_npos_wildtype
data_ag_byday_wide$obs = factor(1:nrow(data_ag_byday_wide))
data_ag_byday_wide = data_ag_byday_wide[data_ag_byday_wide$total != 0, ]
head(data_ag_byday_wide)

write.csv(data_ag_byday_wide, file=".\\data\\be_latest\\be_B117_total.csv", row.names=FALSE)



# 4.1 ESTIMATE GROWTH RATE & TRANSMISSION ADVANTAGE OF 501Y.V1 USING BINOMIAL GLMM (LOGISTIC FIT) ####

# fit common-slope and separate-slopes binomial GLM
set_sum_contrasts()
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit1)
fit1 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)+LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings2)  # common slope model, with lab coded as fixed factor

fit2 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)*LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings) # separate slopes model, with lab coded as fixed factor
BIC(fit1,fit2) 
#      df      BIC
# fit1  8 2176.638
# fit2 13 2198.915


# common-slope model fit1 fits best, i.e. rate at which 501Y.V1 is displacing other strains constant across regions/labs

summary(fit1)
# Random effects:
#   Groups Name        Variance Std.Dev.
# obs    (Intercept) 0.3144   0.5607  
# Number of obs: 495, groups:  obs, 495
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                -2.99346    0.05936 -50.430  < 2e-16 ***
#   scale(collection_date_num)  2.14210    0.05919  36.188  < 2e-16 ***
#   LABORATORY1                -0.18886    0.08902  -2.122   0.0339 *  
#   LABORATORY2                 0.36197    0.08059   4.491 7.08e-06 ***
#   LABORATORY3                 0.48638    0.08080   6.020 1.75e-09 ***
#   LABORATORY4                -0.89337    0.09285  -9.621  < 2e-16 ***
#   LABORATORY5                -0.15269    0.08343  -1.830   0.0672 .  

# growth rate advantage (differences in growth rate between 501Y.V1 and old strains):
# results common-slope model
fit1_emtrends = as.data.frame(emtrends(fit1, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$emtrends)
fit1_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL  asymp.UCL
# 1                0.08935008 0.08451078 0.09418937

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*4.7) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.521878  1.487654   1.55689

# with a generation time of 5.5 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*5.5) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.634645  1.591711  1.678737


# tests for differences in date of introduction
# UCL, ULB & UZA (& Ghent if included in data) earlier than avg, Namur, Mons & UZ Leuven later than avg
emmeans(fit1,eff~LABORATORY)$contrasts 
# contrast                  estimate     SE  df z.ratio p.value
# Namur effect                -0.189 0.0890 Inf -2.122  0.0407 
# (Saint LUC - UCL) effect     0.362 0.0806 Inf  4.491  <.0001 
# ULB effect                   0.486 0.0808 Inf  6.020  <.0001 
# (UMons - Jolimont) effect   -0.893 0.0929 Inf -9.621  <.0001 
# UZ leuven effect            -0.153 0.0834 Inf -1.830  0.0672 
# UZA effect                   0.387 0.0802 Inf  4.818  <.0001 
# # 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: fdr method for 8 tests

# results of growth rate advantage of separate-slopes model fit2 by lab/region:                         
fit2_emtrends = emtrends(fit2, revpairwise ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$emtrends
fit2_emtrends
# LABORATORY       collection_date_num.trend      SE  df asymp.LCL asymp.UCL
# Namur                               0.0952 0.00662 Inf    0.0822    0.1081
# Saint LUC - UCL                     0.0813 0.00499 Inf    0.0715    0.0911
# ULB                                 0.0909 0.00530 Inf    0.0805    0.1013
# UMons - Jolimont                    0.1039 0.00763 Inf    0.0890    0.1189
# UZ leuven                           0.0895 0.00560 Inf    0.0786    0.1005
# UZA                                 0.0832 0.00500 Inf    0.0734    0.0930

# no lab/region displays above-average growth rate of B.1.1.7
fit2_contrasts = emtrends(fit2, eff ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$contrasts
fit2_contrasts
# contrast                   estimate      SE  df z.ratio p.value
# Namur effect               0.004492 0.00587 Inf  0.766  0.9704 
# (Saint LUC - UCL) effect  -0.009385 0.00467 Inf -2.008  0.2396 
# ULB effect                 0.000215 0.00490 Inf  0.044  1.0000 
# (UMons - Jolimont) effect  0.013247 0.00669 Inf  1.980  0.2544 
# UZ leuven effect          -0.001127 0.00509 Inf -0.222  1.0000 
# UZA effect                -0.007442 0.00468 Inf -1.590  0.5092 
# 
# P value adjustment: sidak method for 7 tests 



# PLOT MODEL FIT

# for best fitting common slope model fit1
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-05-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit1_preds = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                                         # by="LABORATORY", 
                                         at=list(collection_date_num=seq(date.from,
                                                                         date.to)), 
                                         type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_preds$collection_date = as.Date(fit1_preds$collection_date_num, origin="1970-01-01")


total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1))$sdcor, function (x) x^2))) 
fit1_preds_bylab = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                                   by="LABORATORY", 
                                   at=list(collection_date_num=seq(date.from,
                                                               date.to)), 
                                    type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_preds_bylab$collection_date = as.Date(fit1_preds_bylab$collection_date_num, origin="1970-01-01")
# order labs by estimated date of introduction (intercepts)
dfemmeanslabs = as.data.frame(emmeans(fit1,~LABORATORY))
levels_BE = as.character(dfemmeanslabs$LABORATORY[order(dfemmeanslabs$emmean,decreasing=T)])
fit1_preds_bylab$LABORATORY = factor(fit1_preds_bylab$LABORATORY, 
                                     levels=levels_BE)


# estimated share of 501Y.V1 among currently diagnosed infections based on fit1
fit1_preds[fit1_preds$collection_date==today,]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 27              18680 0.6846541 0.01451628 Inf 0.6556019  0.712453      2021-02-22

# estimated share of 501Y.V1 among new infections (assuming time between infection & diagnosis of 7 days)
fit1_preds[fit1_preds$collection_date==(today+7),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 34               18687 0.7968834 0.01345049 Inf 0.7693272 0.8220421      2021-03-01


sum(tail(data_ag_byday_wide$est_n_B117, 14))/sum(tail(data_ag_byday_wide$n_pos,14)) 
# 47% of the samples of last 2 weeks in the dataset were estimated to be by British variant
# PS: with data 31/1 this was 15.4%  
# note: this is not the same as the estimated prop of the new infections or new diagnoses today that are of the British
# variant, which are much higher, see above)


# implied Re of wild type and 501Y.V1 given this predicted share of 501Y.V1 among all infections ####
# under a particular fitted transmission advantage
# based on the fact that the overall Re is a weighted average of the Re of the individual variants
# functions to calculate Re of wild type and of 501Y.V1 based on overall Re value and prop of positives that is 501Y.V1 propB117
# and transmission advantage of 501Y.V1 M
M_fitted = exp(fit1_emtrends[,c(2,5,6)]*4.7)[1,2] # 1.487654
Re_wild_type = function (Re, propB117, M=M_fitted) {
  Re / (1-propB117+M*propB117)
}
Re_B117 = function (Re, propB117, M=M_fitted) {
  M*Re / (1-propB117+M*propB117)
}

Re_cases = read.csv(".//Re_fits//Re_cases.csv") 
# Re values calculated from instantaneous growth rate in nr of new cases 
# with instant growth rate = derivative/emtrends in function of sample date of GAM fit on new cases 
# gam(cbind(NEWCASES, totpop-NEWCASES) ~ s(DATE_NUM, bs="cs", k=23, fx=F) + WEEKDAY + s(log(TESTS_ALL), bs="cs", k=5, fx=F), family=binomial(cloglog), data=cases_tot) 
# and with Re calculated using R.from.r with gamma_mean=4.7, gamma_sd=2.9
# Re here is at time of diagnosis, we use a shift of 7 days to recalculate to Re at day of infection
Re_cases$DATE = as.Date(Re_cases$DATE_NUM, origin="1970-01-01")
Re_cases$collection_date_num = Re_cases$DATE_NUM
Re_cases$propB117 = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                      # by="LABORATORY", 
                      at=list(collection_date_num=seq(min(Re_cases$collection_date_num),
                                                      max(Re_cases$collection_date_num))), 
                      type="response"), bias.adjust = TRUE, sigma = total.SD)$prob

head(Re_cases)
Re_cases$Re_WT = Re_wild_type(Re=Re_cases$Re, propB117=Re_cases$propB117)
Re_cases$Re_WT_LOWER = Re_wild_type(Re=Re_cases$Re_LOWER, propB117=Re_cases$propB117)
Re_cases$Re_WT_UPPER = Re_wild_type(Re=Re_cases$Re_UPPER, propB117=Re_cases$propB117)
Re_cases$Re_B117 = Re_B117(Re=Re_cases$Re, propB117=Re_cases$propB117)
Re_cases$Re_B117_LOWER = Re_B117(Re=Re_cases$Re_LOWER, propB117=Re_cases$propB117)
Re_cases$Re_B117_UPPER = Re_B117(Re=Re_cases$Re_UPPER, propB117=Re_cases$propB117)

# d = as.Date(max(Re_cases$DATE))
# tag = paste("@TWenseleers\n",dat)

# PS we plot Re at time of infection here by shifting plot by 7 days
qplot(data=Re_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("grey")) +
  geom_line() + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous(limits=c(1/3,3), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re OF 501Y.V1 (red) AND WILD TYPE (blue) IN BELGIUM\nBASED ON NEW CONFIRMED CASES AND\nESTIMATED 501Y.V1 TRANSMISSION ADVANTAGE") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  geom_ribbon(aes(y=Re_WT, ymin=Re_WT_LOWER, ymax=Re_WT_UPPER), fill=I("blue"), alpha=I(0.5)) +
  geom_line(aes(y=Re_WT), fill=I("blue"), alpha=I(0.8)) +
  geom_ribbon(data=Re_cases[Re_cases$DATE>=as.Date("2021-01-01"),], 
              aes(y=Re_B117, ymin=Re_B117_LOWER, ymax=Re_B117_UPPER), fill=I("red"), alpha=I(0.5)) +
  geom_line(data=Re_cases[Re_cases$DATE>=as.Date("2021-01-01"),],
            aes(y=Re_B117), fill=I("blue"), alpha=I(0.8)) +
  coord_cartesian(xlim=c(as.Date("2020-08-20"),today),
                  ylim=c(0.6,1.6)) +
  ylab("Re at time of infection")
ggsave(file=paste0(".//plots//",dat,"//Re_cases_Re_501YV1_Re_wildtype.png"), width=7, height=5)
Re_cases[Re_cases$DATE==max(Re_cases$DATE),]
# DATE_NUM          r          SE      df     r_LOWER    r_UPPER       DATE      Re Re_LOWER Re_UPPER collection_date_num  propB117
# 365    18686 0.01226322 0.001377151 320.037 0.009553808 0.01497263 2021-02-28 1.05867 1.045529 1.071913               18686 0.7827269
# Re_WT Re_WT_LOWER Re_WT_UPPER  Re_B117 Re_B117_LOWER Re_B117_UPPER
# 365 0.7662085   0.7566977    0.775793 1.139854      1.125705      1.154112




# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50%, 75% 90% by 501Y.V1 :
fit1_preds$collection_date[fit1_preds[,"prob"]>=0.5][1]-7 # >50% by 6th of February [5th Feb - 7 Feb] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.5][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.5][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.75][1]-7 # >75% by 19th of February [18 Feb - 21 Feb] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.75][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.75][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.9][1]-7 # >90% by 4th of March [2 March - 7 March] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.9][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.9][1]-7



# PLOT MODEL FIT common-slope model fit1

# plot for the whole of Belgium on response scale:
plot_fit1_response = qplot(data=fit1_preds, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~LABORATORY) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, 
                  # fill=LABORATORY
  ), 
  fill=I("#b0c4de"), 
  alpha=I(1)) +
  geom_line(aes(y=100*prob, 
                # colour=LABORATORY
  ), 
  colour=I("steelblue"), 
  alpha=I(1)) +
  ylab("Share of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC 202012/01 BY NHS REGION") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(data_ag_byday_wide$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_byday_wide, 
             aes(x=collection_date, y=100*propB117, size=n_pos,
                 # colour=LABORATORY
             ), 
             colour=I("steelblue"), 
             alpha=I(1)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(10,10^round(log10(max(data_ag_byday_wide$n_pos)+1),0)), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date")
plot_fit1_response


saveRDS(plot_fit1_response, file = paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_response scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_response scale.pdf"), width=8, height=6)


# plot for the whole of Belgium on logit scale:
plot_fit1_link = qplot(data=fit1_preds, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~LABORATORY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  # fill=LABORATORY
  ), 
  fill=I("#b0c4de"), 
  alpha=I(1)) +
  geom_line(aes(y=prob, 
                # colour=LABORATORY
  ), 
  colour=I("steelblue"), 
  alpha=I(1)) +
  ylab("Share of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC 202012/01 BY NHS REGION") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(data_ag_byday_wide$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.0001,0.999), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_byday_wide, 
             aes(x=collection_date, y=propB117, size=n_pos,
                 # colour=LABORATORY
             ), 
             colour=I("steelblue"), 
             alpha=I(1)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(10,10^round(log10(max(data_ag_byday_wide$n_pos)+1),0)), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date")
plot_fit1_link


saveRDS(plot_fit1_response, file = paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_logit scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_logit scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_501YV1_Belgium_logit scale.pdf"), width=8, height=6)


# plot per lab on logit scale:
plot_fit1 = qplot(data=fit1_preds_bylab, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~LABORATORY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=LABORATORY
                  ), 
              # fill=I("steelblue"), 
              alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=LABORATORY
                ), 
            # colour=I("steelblue"), 
            alpha=I(0.8)) +
  ylab("Share of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(data_ag_byday_wide$collection_date), as.Date("2021-04-01")-1), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.01,0.99), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_wide,  
             aes(x=collection_date, y=propB117, size=n_pos,
                 colour=LABORATORY
                 ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="identity", 
                        range=c(1, 3), limits=c(1,max(data_ag_wide$n_pos)), breaks=c(1,10,100)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
plot_fit1

saveRDS(plot_fit1, file = paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_501YV1_Belgium by lab.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_501YV1_Belgium by lab.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_501YV1_Belgium by lab.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_501YV1_Belgium by lab.pdf"), width=8, height=6)





# 4.2 ESTIMATE GROWTH RATE & R VALUE OF 501Y.V1 & WILD TYPE STRAINS SEPARATELY USING MULTINOMIAL MODEL ####

# Function to recalculate Malthusian growth rate to effective reproduction number Rt
# (assuming generation time is gamma distributed)
# Ref: Park et al. 2020, https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
R.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) {
  k = (gamma_sd / gamma_mean)^2
  R = (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

data_ag_wide2 = data_ag_wide 
data_ag_wide2$n_b117 = round(data_ag_wide2$est_n_B117, 0)
data_ag_wide2$n_spos = data_ag_wide2$n_pos-data_ag_wide2$n_b117
(data_ag_wide2$n_neg+data_ag_wide2$n_spos+data_ag_wide2$n_b117)==data_ag_wide2$total # check
data_ag_wide2 = data_ag_wide2[,c("collection_date","LABORATORY","n_neg","n_spos","n_b117")]
head(data_ag_wide2)
data_ag_long = gather(data_ag_wide2, outcome, count, n_neg:n_b117, factor_key=TRUE)
data_ag_long$outcome = factor(data_ag_long$outcome, levels=c("n_neg","n_spos","n_b117"))
data_ag_long$collection_date_num = as.numeric(data_ag_long$collection_date)

test_outcomes = ggplot(data=data_ag_long, 
      aes(x=collection_date, 
      y=count, fill=outcome, group=outcome)) +
  facet_wrap(~LABORATORY) +
  geom_area(aes(fill=outcome), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("test outcome", values=c("darkgrey","steelblue","lightcoral"), labels=c("negative","S positive","S dropout")) +
  ylab("Share") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  # ggtitle("Test outcomes") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right") +
  coord_cartesian(xlim=c(min(data_ag_long$collection_date), as.Date("2021-02-10"))) 
# PS - I still need to fix a bug here - the negative samples from the last file from february (since febr 10) are missing....
test_outcomes

saveRDS(test_outcomes, file = paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.pdf"), width=7, height=5)

test_outcomes_pos = ggplot(data=data_ag_long[data_ag_long$outcome!="n_neg",], 
                       aes(x=collection_date, 
                           y=count, fill=outcome, group=outcome)) +
  facet_wrap(~LABORATORY) +
  geom_area(aes(fill=outcome), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("test outcome", values=c("steelblue","lightcoral"), labels=c("S positive","S dropout")) +
  ylab("Share") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  # ggtitle("Test outcomes") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right")
test_outcomes_pos

saveRDS(test_outcomes, file = paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.pdf"), width=7, height=5)


# multinomial spline fit on test outcome data (negative / positive wild type / positive 501Y.V1
# to be able to estimate growth rate and Rt of 501Y.V1 and wild type separately

set.seed(1)
# we use data from the 14th of Jan onwards, as data has been approx randomly sampled from then on
sel_labs = unique(data_ag_long$LABORATORY)
date.from = as.Date("2021-01-14")
date.to = as.Date("2021-02-10") 
# PS now just using data up till Feb 10, because I have a bug in the code that causes negative samples to be left out since then
# STILL TO FIX
data_ag_long_subs = data_ag_long[(data_ag_long$LABORATORY %in% sel_labs)&(data_ag_long$collection_date>=date.from)&(data_ag_long$collection_date<=date.to),]
data_ag_long_subs$LABORATORY = droplevels(data_ag_long_subs$LABORATORY)

mfit0 = nnet::multinom(outcome ~ scale(collection_date_num) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mfit1 = nnet::multinom(outcome ~ scale(collection_date_num) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mfit2 = nnet::multinom(outcome ~ ns(collection_date_num, df=2) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit3 = nnet::multinom(outcome ~ ns(collection_date_num, df=3) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit4 = nnet::multinom(outcome ~ ns(collection_date_num, df=4) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit5 = nnet::multinom(outcome ~ ns(collection_date_num, df=5) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
BIC(mfit0, mfit1, mfit2, mfit3, mfit4, mfit5) # mfit3 fits best, df splines tuned based on BIC
#       df      BIC
# mfit0 16 150913.8
# mfit1 28 150894.7
# mfit2 18 149513.9
# mfit3 56 149434.8
summary(mfit4)

plot(Effect("collection_date_num",mfit4), style="stacked")
plot(Effect("collection_date_num",mfit4), confint=list(style="bands"), rug=FALSE)

# average growth rates of S-positive/wild type & S dropout/501Y.V1 cases evaluated today based on best fit multinomial model mfit4
emtrends(mfit4, ~outcome|1, var="collection_date_num",  at=list(collection_date_num=today_num), mode="latent")
# emmeans_1.5.3 output:
# outcome collection_date_num.trend      SE df lower.CL upper.CL
# n_neg                    -0.05648 0.00646 60  -0.0694  -0.0436
# n_spos                   -0.00276 0.00882 60  -0.0204   0.0149
# n_b117                    0.05924 0.01103 60   0.0372   0.0813
R.from.r(-0.00276) # Rt S pos / wild type = 0.99
R.from.r(-0.0204) # Rt S pos / wild type LCL = 0.91
R.from.r(0.0149) # Rt S pos / wild type UCL = 1.07

R.from.r(0.05924) # Rt of S dropout = 1.30 [1.18-1.43]
R.from.r(0.0372) # Rt of S dropout = 0.78
R.from.r(0.0813) # Rt of S dropout = 0.95


R.from.r(0.05924)/R.from.r(-0.00276) # Rt of 501Y.V1 = 1.32x times higher than of wild type


# implied growth rate & transmission advantage
delta_r = data.frame(confint(contrast(emtrends(mfit4, ~outcome|1, var="collection_date_num",  
                                               at=list(outcome=c("n_spos","n_b117"),
                                                       collection_date_num=today_num), mode="latent"), method="revpairwise")))[,c(2,5,6)]
delta_r # growth advantage
#       estimate lower.CL  upper.CL
# 1 0.06199392 0.02419276 0.09979507
exp(delta_r*4.7) # transmission advantage
#    estimate lower.CL upper.CL
# 1 1.338262 1.120423 1.598454


# growth rates and Re values of the 501Y.V1 variant and the wild type calculated over time
extrapolate = 0
r_and_Re_B117_wildtype = data.frame(emtrends(mfit3, ~outcome|1, var="collection_date_num", by=c("collection_date_num"), # by=c("collection_date_num","LABORATORY"),  
                                  at=list(collection_date_num=seq(min(data_ag_long_subs$collection_date_num),
                                                                  today_num+extrapolate),
                                          outcome=c("n_spos","n_b117")), mode="latent"))
r_and_Re_B117_wildtype$collection_date = as.Date(r_and_Re_B117_wildtype$collection_date_num, origin="1970-01-01")
colnames(r_and_Re_B117_wildtype)[colnames(r_and_Re_B117_wildtype) %in% c("collection_date_num.trend","lower.CL","upper.CL")] = c("r","r.LCL","r.UCL")
r_and_Re_B117_wildtype$Re = R.from.r(r_and_Re_B117_wildtype$r)
r_and_Re_B117_wildtype$Re.LCL = R.from.r(r_and_Re_B117_wildtype$r.LCL)
r_and_Re_B117_wildtype$Re.UCL = R.from.r(r_and_Re_B117_wildtype$r.UCL)

plot_Re_B117_WT = qplot(data=r_and_Re_B117_wildtype, x=collection_date, y=Re, ymin=Re.LCL, ymax=Re.UCL, geom="ribbon", alpha=I(0.5), 
      fill=outcome , colour=NULL, group=outcome ) +
  # facet_wrap(~LABORATORY) +
  geom_line(aes(colour=outcome )) + theme_hc() + xlab("") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  # scale_y_continuous(trans="log2", breaks=c(1/seq(3,1),seq(1,4)),
  #                   labels=round(c(1/seq(3,1),seq(1,4)),2)) +
  coord_cartesian(xlim=c(min(data_ag_long_subs$collection_date),
                         max(r_and_Re_B117_wildtype$collection_date)), 
                  ylim=c(0.6,1.8), 
                  expand=c(0,0)) +
  geom_hline(yintercept=1, colour=alpha(I("black"),0.2)) +
  # theme(legend.position = "none") +
  ggtitle("Effective reproduction nr. Re of 501Y.V1 and wild type") +
  guides(colour=FALSE) +
  scale_colour_manual("", values=c("steelblue","lightcoral")) +
  scale_fill_manual("", values=c("steelblue","lightcoral"), labels=c("S positive (wild type)","S dropout (501Y.V1)")) +
  theme(legend.position = "bottom")
  # labs(tag = tag) +
  # theme(plot.tag.position = "bottomright",
  #      plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_Re_B117_WT

# implied growth rate advantage + 95% CLs over time
extrapolate = 0
growthadvantage = data.frame(confint(contrast(emtrends(mfit3, ~outcome|1, var="collection_date_num", by="collection_date_num", 
                                               at=list(outcome=c("n_spos","n_b117"),
                                                       collection_date_num=seq(min(data_ag_long$collection_date_num),
                                                                               today_num+extrapolate)), 
                                               mode="latent"), method="revpairwise")))
colnames(growthadvantage)[colnames(growthadvantage) %in% c("estimate","lower.CL","upper.CL")] = c("delta_r","delta_r.LCL","delta_r.UCL")
growthadvantage$transmadv = (exp(4.7*growthadvantage$delta_r)-1)*100
growthadvantage$transmadv.LCL = (exp(4.7*growthadvantage$delta_r.LCL)-1)*100
growthadvantage$transmadv.UCL = (exp(4.7*growthadvantage$delta_r.UCL)-1)*100
growthadvantage$collection_date = as.Date(growthadvantage$collection_date_num, origin="1970-01-01")

plot_growthadvB117 = qplot(data=growthadvantage, x=collection_date, y=transmadv, 
                            ymin=transmadv.LCL, ymax=transmadv.UCL, geom="ribbon", alpha=I(0.5), 
      fill=I("lightcoral") , colour=NULL ) +
  geom_line(aes(colour=I("lightcoral"))) + theme_hc() + xlab("") + ylab("Transmission advantage (%)") +
  geom_hline(yintercept=0, colour=alpha(I("black"),0.2)) +
  ggtitle("Transmission advantage of 501Y.V1 over wild type") +
  coord_cartesian(xlim=c(min(data_ag_long_subs$collection_date),
                         max(r_and_Re_B117_wildtype$collection_date)), 
                  expand=c(0,0))
# labs(tag = tag) +
# theme(plot.tag.position = "bottomright",
#      plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_growthadvB117

multipanel_Re_growthadv_multinom = ggarrange(plot_Re_B117_WT, plot_growthadvB117, ncol=1)
multipanel_Re_growthadv_multinom

saveRDS(multipanel_Re_growthadv_multinom, file = paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.pptx"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.png"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.pdf"), width=7, height=8)








# 5. SOME INTERNATIONAL COMPARISONS ####

# SEE SECTION 5. IN https://github.com/nicholasdavies/newcovid/blob/master/multinomial_logistic_fits/multinomial%20logistic%20fits_FINAL.R
