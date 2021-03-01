# ANALYSIS OF S-GENE TARGET FAILURE (S DROPOUT) DATA FROM BELGIUM TO INFER CONTAGIOUSNESS OF NEW VARIANT OF CONCERN B.1.1.7 / 501Y.V1 ####
# PLUS INTERNATIONAL COMPARISON (USING DATA FROM THE UK, DENMARK, SWITZERLAND & THE US)
# AND FIRST ASSESSMENT OF GROWTH ADVANTAGE OF THE SOUTH AFRICAN VOC 501Y.V2 & BRAZILIAN VOC 501Y.V3 BASED ON SEQUENCING DATA
# Tom Wenseleers & Niel Hens
# All Belgian data provided by Emmanuel André

# Associated report: 
# https://github.com/tomwenseleers/newcovid_belgium/blob/main/reports/Genomic%20surveillance%20update_23%20feb%202021.pdf

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


dat="2021_02_23" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-02-23")
today_num = as.numeric(today)
today # "2021-02-23"

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOC 501Y.V1,VOC 501Y.V2&VOC 501Y.V3 IN BELGIUM BASED ON BASELINE SEQUENCING DATA ####
# (baseline sequencing results, i.e. randomly sampled)

be_seqdata = read.csv(paste0(".\\data\\",dat,"\\sequencing_501YV1_501YV2_501YV3.csv"))
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
range(be_seqdata$collection_date)


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
be_seq_mfitallVOC = nnet::multinom(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                                   subset=be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="501Y.V1+V2+V3", maxit=1000) 
summary(be_seq_mfit0)

# growth rate advantage compared to wild type
delta_r_501V1_501YV2_501YV3 = data.frame(confint(emtrends(be_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_501V1_501YV2_501YV3) = delta_r_501V1_501YV2_501YV3[,"contrast"]
delta_r_501V1_501YV2_501YV3 = delta_r_501V1_501YV2_501YV3[,-1]
delta_r_501V1_501YV2_501YV3
#                       estimate  asymp.LCL  asymp.UCL
# 501Y.V1 - wild type 0.06966334 0.06119481 0.07813188
# 501Y.V2 - wild type 0.05340663 0.03519376 0.07161950
# 501Y.V3 - wild type 0.12353360 0.05273034 0.19433686

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(be_seq_mfit0, pairwise ~ variant|1, 
                   var="collection_date_num",  mode="latent",
                   at=list(collection_date_num=today_num), 
          df=NA)$contrasts
# contrast            estimate      SE df z.ratio p.value
# wild type - 501Y.V1  -0.0697 0.00432 NA -16.123 <.0001 
# wild type - 501Y.V2  -0.0534 0.00929 NA  -5.747 <.0001 
# wild type - 501Y.V3  -0.1235 0.03612 NA  -3.420 0.0035 
# 501Y.V1 - 501Y.V2     0.0163 0.00971 NA   1.675 0.3371 
# 501Y.V1 - 501Y.V3    -0.0539 0.03614 NA  -1.490 0.4433 
# 501Y.V2 - 501Y.V3    -0.0701 0.03711 NA  -1.890 0.2325 
# 
# Degrees-of-freedom method: user-specified 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# pairwise contrasts in growth rate (here without Tukey correction)
emtrends(be_seq_mfit0, pairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none")$contrasts
# contrast            estimate      SE df z.ratio p.value
# wild type - 501Y.V1  -0.0697 0.00432 NA -16.123 <.0001 
# wild type - 501Y.V2  -0.0534 0.00929 NA  -5.747 <.0001 
# wild type - 501Y.V3  -0.1235 0.03612 NA  -3.420 0.0006 
# 501Y.V1 - 501Y.V2     0.0163 0.00971 NA   1.675 0.0940 
# 501Y.V1 - 501Y.V3    -0.0539 0.03614 NA  -1.490 0.1361 
# 501Y.V2 - 501Y.V3    -0.0701 0.03711 NA  -1.890 0.0588 
# 
# Degrees-of-freedom method: user-specified 


# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_501V1_501YV2_501YV3*4.7) 
#                     estimate asymp.LCL asymp.UCL
# 501Y.V1 - wild type 1.387381  1.333245  1.443715
# 501Y.V2 - wild type 1.285324  1.179878  1.400195
# 501Y.V3 - wild type 1.787124  1.281245  2.492742

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




# # PS: mblogit fit would also be possible & would take into account overdispersion
# be_basseqdata_long$obs = factor(1:nrow(be_basseqdata_long))
# be_seq_mblogitfit = mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                             # random = ~ 1|obs,
#                             weights = count, data=be_basseqdata_long, 
#                             subset=be_basseqdata_long$variant!="501Y.V1+V2+V3",
#                             dispersion = FALSE)
# dispersion(mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                    # random = ~ 1|obs,
#                    weights = count, data=be_basseqdata_long,
#                    subset = be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="501Y.V1+V2+V3",
#                    dispersion = TRUE), method="Afroz") # dispersion coefficient = 3.2


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
           xmax=as.Date("2021-04-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("variant", values=c("darkgrey","blue","red","green3")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-04-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share among newly diagnosed infections") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")
muller_be_seq_mfit0

saveRDS(muller_be_seq_mfit0, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.pdf"), width=7, height=5)

multinom_501YV1_501YV2_501YV3 = ggarrange(baseline_sequencing+ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")+xlab("")+ylab("Share")+
                                            coord_cartesian(xlim=c(as.Date("2020-12-01"),as.Date("2021-04-01"))), 
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
  ylab("Share among newly diagnosed infections (%)") +
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
  ylab("Share among newly diagnosed infections (%)") +
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
# 501Y.V1 (British)               18681 0.615664 0.02353748 NA 0.5695314 0.6617966      2021-02-23

# estimated proportion of 501Y.V1 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18688 0.6958476 0.03143403 NA 0.6342381 0.7574572      2021-03-02

# estimated proportion of 501Y.V2 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18681 0.0586159 0.01160612 NA 0.03586833 0.08136348      2021-02-23

# estimated proportion of 501Y.V2 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18688 0.05912409 0.01512296 NA 0.02948362 0.08876455      2021-03-02

# estimated proportion of 501Y.V3 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18681 0.01987965 0.01105983 NA -0.001797221 0.04155653      2021-02-23

# estimated proportion of 501Y.V3 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18688 0.03276032 0.02512898 NA -0.01649158 0.08201221      2021-03-02



# estimated proportion of one of the three VOCs among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18681 0.6910092 0.02021192 NA 0.6513945 0.7306238      2021-02-23

# estimated proportion of one of the three VOCs among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18688 0.7827108 0.02055526 NA 0.7424233 0.8229984      2021-03-02




# the time at which new lab diagnoses would be by more than 50%, 75% 90% by one of the three VOCs :
be_seq_mfit0_preds2_subs = be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date>as.Date("2021-01-01")&
                                                 be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
# >50% by 12th of February [10th Feb - 13th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]

# >75% by 28th of February [25th Feb - 3d March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]

# >90% by 16th of March [11th March - 21st of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]


# the time at which new infections would be by more than 50%, 75% 90% by one of the three VOCs
# (counting 7 days between infection & diagnosis) :
# >50% by 5th of Feb [3d Feb - 6th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]-7

# >75% by 21st of February [18th Feb - 24th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]-7

# >90% by 9th of March [4th March - 14th of March] 95% CLs
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
# 62           18681 0.9975054 0.001975727 Inf 0.9882887 0.9994728      2021-02-23

# prop of S dropout samples among new infections that are now estimated to be B.1.1.7 / 501Y.V1 (using 7 days for time from infection to diagnosis)
fitseq_preds[fitseq_preds$collection_date==(today+7),]
#    collection_date_num     prob          SE  df asymp.LCL asymp.UCL collection_date
# 69           18688 0.9988497 0.001042089 Inf 0.9932354 0.9998054      2021-03-02

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
file_feb2 = paste0(".//data//", dat, "//PCR February 2021 10 to 23.xlsx")
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
range(as.Date(as.numeric(ctdata_feb2$Analysis_date), origin="1899-12-30")) # "2021-02-10" "2021-02-23"
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
unique(ctdata_feb2$Outcome)
ctdata$Outcome[ctdata$Outcome=="Detected"] = "Positive"
ctdata$Outcome[ctdata$Outcome=="Not detected"] = "Negative"
unique(ctdata$Outcome) # "Positive" "Negative"
ctdata$Analysis_date = as.Date(as.numeric(ctdata$Analysis_date), origin="1899-12-30")
sum(is.na(ctdata$Analysis_date)) # 0
range(ctdata$Analysis_date) # "2020-12-01" - "2021-02-23"
ctdata$collection_date = ctdata$Analysis_date-1 # collection date = analysis date-1 
sum(is.na(ctdata$collection_date)) # 0
ctdata$collection_date_num = as.numeric(ctdata$collection_date)
range(ctdata$collection_date) # "2020-11-30" "2021-02-22"
ctdata$group = interaction(ctdata$Outcome, ctdata$S_dropout)
ctdata$group = droplevels(ctdata$group)
unique(ctdata$group) # Positive.0 Positive.1 Negative.0
ctdata$Outcome = factor(ctdata$Outcome)
unique(ctdata$Outcome) # Positive Negative
ctdata$Laboratory = factor(ctdata$Laboratory)
ctdata$S_dropout = factor(ctdata$S_dropout)
head(ctdata)
str(ctdata)
nrow(ctdata) # 651273

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

cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.75

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                        Laboratory correlation_Ct_N_ORF1ab
# Namur                       Namur              0.94723044
# Saint LUC - UCL   Saint LUC - UCL              0.98145030
# ULB                           ULB              0.98106595
# ULG                           ULG              0.96102935
# UMons - Jolimont UMons - Jolimont             -0.28112232
# UZ Gent                   UZ Gent             -0.01491232
# UZ leuven               UZ leuven              0.97663099
# UZA                           UZA              0.86651032

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



# plot & analysis of Ct values for 4 labs with most comparable Ct value distributions ("UZ leuven","Saint LUC - UCL","ULB","Namur") 
# for dates from 13th of Jan onward when >80% of all S dropouts were B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections

# labs_to_remove = c("UMons - Jolimont", "ULG", "UZ Gent", "UZA")
# sel_labs = setdiff(unique(ctdata_onlypos_subs$Laboratory), labs_to_remove) 
sel_labs = c("UZ leuven","Saint LUC - UCL","ULB","Namur")  
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
# which was a Ct value < 15.04 for the N gene and < 15.98 for the ORF1ab gene
thresh_N = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","N_cq"]))/1.25
thresh_N # 15.04
thresh_ORF1 = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","ORF1_cq"]))/1.25
thresh_ORF1 # 15.98
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="N gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="N gene"]<thresh_N
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"]<thresh_ORF1


# check correlation between Ct values for N & ORF1ab gene
cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.97, t=504.75, p<2E-16

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                      Laboratory correlation_Ct_N_ORF1ab
# Namur                     Namur               0.9472304
# Saint LUC - UCL Saint LUC - UCL               0.9814503
# ULB                         ULB               0.9810659
# UZ leuven             UZ leuven               0.9766310


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
AIC(qr_bothgenes0, k=-1) # 182174.7
AIC(qr_bothgenes1, k=-1) # 182157.5
AIC(qr_bothgenes2, k=-1) # 182070.4
AIC(qr_bothgenes3, k=-1) # 182191.7
AIC(qr_bothgenes4, k=-1) # 182064.8 # fits data best based on BIC criterion (lowest value, PS: here AIC with k<0 returns BIC)
AIC(qr_bothgenes5, k=-1) # 182094.3

summary(qr_bothgenes4)
# tau: [1] 0.5
# 
# Coefficients:
#   Value     Std. Error t value   Pr(>|t|) 
# (Intercept)             18.26135   0.06120  298.38200   0.00000
# Gene1                   -0.75436   0.05983  -12.60889   0.00000
# Laboratory1             -0.48108   0.10872   -4.42474   0.00001
# Laboratory2             -0.50353   0.10757   -4.68105   0.00000
# Laboratory3             -0.70263   0.09900   -7.09712   0.00000
# S_dropout1               1.19384   0.06094   19.58885   0.00000
# Gene1:Laboratory1        0.17739   0.10578    1.67689   0.09358
# Gene1:Laboratory2        0.08599   0.10879    0.79040   0.42930
# Gene1:Laboratory3       -0.21286   0.09807   -2.17046   0.02998
# Gene1:S_dropout1         0.26298   0.06108    4.30551   0.00002
# Laboratory1:S_dropout1  -0.14241   0.10853   -1.31224   0.18945
# Laboratory2:S_dropout1  -0.00246   0.10750   -0.02291   0.98173
# Laboratory3:S_dropout1  -0.68401   0.09858   -6.93856   0.00000

qr_emmeans_bylab99 = data.frame(emmeans(qr_bothgenes4, ~ Laboratory + Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans_bylab99
qr_emmeans99 = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans99
qr_emmeans = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.95)) # median Ct values + 95% CLs
qr_emmeans
# Gene S_dropout   emmean        SE    df lower.CL upper.CL
# 1      N gene         0 18.96380 0.1056879 27411 18.75665 19.17095
# 2 ORF1ab gene         0 19.94657 0.1123086 27411 19.72644 20.16671
# 3      N gene         1 16.05017 0.1349678 27411 15.78563 16.31472
# 4 ORF1ab gene         1 18.08485 0.1306900 27411 17.82869 18.34101

# mean difference in median Ct value of 2.39, which is highly significant across both genes: p<0.0001
contrast(emmeans(qr_bothgenes4, ~ S_dropout, level=0.95), method="pairwise") 
# contrast estimate    SE    df t.ratio p.value
# 0 - 1        2.39 0.122 27411 19.589  <.0001 
# mean difference in median Ct value of 2.91 for N gene and 1.86 for ORF1ab gene
confint(contrast(emmeans(qr_bothgenes4, ~ S_dropout|Gene, level=0.95), method="pairwise"))
# Gene = N gene:
#   contrast estimate    SE    df lower.CL upper.CL
# 0 - 1        2.91 0.172 27411     2.58     3.25
# 
# Gene = ORF1ab gene:
#   contrast estimate    SE    df lower.CL upper.CL
# 0 - 1        1.86 0.174 27411     1.52     2.20
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
# fitct_highvirload_0A  4 34873.68
# fitct_highvirload_1A  5 34870.51
# fitct_highvirload_2A  6 34878.14
# fitct_highvirload_3A  6 34880.66
# fitct_highvirload_4A  7 34888.29
# fitct_highvirload_0B  5 34855.42
# fitct_highvirload_1B  6 34852.21
# fitct_highvirload_2B  9 34879.89
# fitct_highvirload_3B  7 34862.36
# fitct_highvirload_4B 10 34890.07


# fitct_highvirload_0B almost the best model (fitct_highvirload_1B only very slightly better)
summary(fitct_highvirload_0B) # S dropout samples more frequently have high viral load based on N gene Ct values
# Random effects:
#   Groups     Name        Variance Std.Dev.
# Laboratory (Intercept) 0.0305   0.1746  
# Number of obs: 27424, groups:  Laboratory, 4
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.60361    0.08834  -6.832 8.35e-12 ***
#   Gene1             0.11979    0.01350   8.870  < 2e-16 ***
#   S_dropout1       -0.16823    0.01357 -12.399  < 2e-16 ***
#   Gene1:S_dropout1 -0.07201    0.01350  -5.333 9.69e-08 ***

plot(allEffects(fitct_highvirload_0B))


# odds to encounter high viral load samples based on Ct values of both genes (high vir load = Ct values >1.25x lower than median in non-S dropout samples)
# 1.40x [1.33-1.48x] 95% CLs increased among S dropout samples
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout, type="response"), method="revpairwise", type="response"))
# contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.4 0.038 Inf      1.33      1.48

# odds ratio = 1.62 [1.50-1.74] for N gene & 1.21 [1.12-1.31] for ORF1ab gene 
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout|Gene, type="response"), method="revpairwise", type="response"))
# Gene = N gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.62 0.0609 Inf      1.50      1.74
# 
# Gene = ORF1ab gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.21 0.0471 Inf      1.12      1.31


fitct_highvirload_emmeans = as.data.frame(emmeans(fitct_highvirload_0B, ~ S_dropout+Gene, type="response"))
fitct_highvirload_emmeans$S_dropout = factor(fitct_highvirload_emmeans$S_dropout)
fitct_highvirload_emmeans$Gene = factor(fitct_highvirload_emmeans$Gene)
fitct_highvirload_emmeans
# S_dropout        Gene      prob         SE  df asymp.LCL asymp.UCL
# 1         0      N gene 0.3264970 0.01980731 Inf 0.2889238 0.3664383
# 2         1      N gene 0.4394043 0.02277381 Inf 0.3953723 0.4844112
# 3         0 ORF1ab gene 0.3058449 0.01914409 Inf 0.2696590 0.3445955
# 4         1 ORF1ab gene 0.3481491 0.02107303 Inf 0.3080651 0.3905052

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
sel_labs = setdiff(unique(ctdata$Laboratory), c("ULG - FF 3.x", "ULG")) # , "UZ Gent"))  # unique(ctdata$Laboratory) # 
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
nrow(ctdata) # 651273
nrow(ctdata_subs) # 503320

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
fit1_22jan = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)+LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0&data_ag_wide$collection_date>=as.Date("2021-01-01")&data_ag_wide$collection_date<=as.Date("2021-01-22"), control=glmersettings2)  # common slope model, with lab coded as fixed factor
fit1 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)+LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings2)  # common slope model, with lab coded as fixed factor

fit2 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)*LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings) # separate slopes model, with lab coded as fixed factor
BIC(fit1,fit2) 
#      df      BIC
# fit1  8 2129.427
# fit2 13 2150.914


# common-slope model fit1 fits best, i.e. rate at which 501Y.V1 is displacing other strains constant across regions/labs

summary(fit1)
# Random effects:
#   Groups Name        Variance Std.Dev.
# obs    (Intercept) 0.3134   0.5598  
# Number of obs: 499, groups:  obs, 499
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                -2.94002    0.05825 -50.469  < 2e-16 ***
#   scale(collection_date_num)  2.15662    0.05825  37.024  < 2e-16 ***
#   LABORATORY1                -0.14878    0.08866  -1.678   0.0933 .  
# LABORATORY2                 0.35760    0.08104   4.413 1.02e-05 ***
#   LABORATORY3                 0.47501    0.08079   5.879 4.12e-09 ***
#   LABORATORY4                -0.89283    0.09430  -9.468  < 2e-16 ***
#   LABORATORY5                -0.17346    0.08395  -2.066   0.0388 *  

# growth rate advantage (differences in growth rate between 501Y.V1 and old strains):
# results common-slope model
fit1_emtrends = as.data.frame(emtrends(fit1, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$emtrends)
fit1_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL  asymp.UCL
# 1                0.08835838 0.08368093 0.09303582

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*4.7) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.514801  1.481863  1.548472

# with a generation time of 5.5 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*5.5) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.625753  1.584462   1.668121

# results original first first report of 28 jan (using data from jan 1 tot jan 22)
fit1_22jan_emtrends = as.data.frame(emtrends(fit1_22jan, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$emtrends)
fit1_22jan_emtrends[,c(2,5,6)]
# collection_date_num.trend  asymp.LCL asymp.UCL
# 1                 0.1111769 0.08615571 0.1361981

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_22jan_emtrends[,c(2,5,6)]*4.7) 
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.686291    1.4992   1.89673



# tests for differences in date of introduction
# UCL, ULB, UZA (& Ghent) earlier than avg, Namur, Mons & UZ Leuven later than avg
emmeans(fit1,eff~LABORATORY)$contrasts 
# contrast                  estimate     SE  df z.ratio p.value
# Namur effect                -0.149 0.0887 Inf -1.678  0.0933 
# (Saint LUC - UCL) effect     0.358 0.0810 Inf  4.413  <.0001 
# ULB effect                   0.475 0.0808 Inf  5.879  <.0001 
# (UMons - Jolimont) effect   -0.893 0.0943 Inf -9.468  <.0001 
# UZ leuven effect            -0.173 0.0840 Inf -2.066  0.0466 
# UZA effect                   0.382 0.0801 Inf  4.774  <.0001  
# # 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: fdr method for 8 tests

# results of growth rate advantage of separate-slopes model fit2 by lab/region:                         
fit2_emtrends = emtrends(fit2, revpairwise ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$emtrends
fit2_emtrends
# LABORATORY       collection_date_num.trend      SE  df asymp.LCL asymp.UCL
# Namur                               0.0957 0.00630 Inf    0.0834    0.1080
# Saint LUC - UCL                     0.0804 0.00488 Inf    0.0708    0.0899
# ULB                                 0.0891 0.00507 Inf    0.0792    0.0991
# UMons - Jolimont                    0.1031 0.00743 Inf    0.0886    0.1177
# UZ leuven                           0.0869 0.00539 Inf    0.0764    0.0975
# UZA                                 0.0828 0.00481 Inf    0.0734    0.0922

# no lab/region displays above-average growth rate of B.1.1.7
fit2_contrasts = emtrends(fit2, eff ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$contrasts
fit2_contrasts
# contrast                   estimate      SE  df z.ratio p.value
# Namur effect               0.006014 0.00559 Inf  1.075  0.8633 
# (Saint LUC - UCL) effect  -0.009313 0.00455 Inf -2.046  0.2209 
# ULB effect                -0.000548 0.00469 Inf -0.117  1.0000 
# (UMons - Jolimont) effect  0.013457 0.00651 Inf  2.069  0.2103 
# UZ leuven effect          -0.002737 0.00490 Inf -0.559  0.9942 
# UZA effect                -0.006874 0.00450 Inf -1.526  0.5573 
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

# original fit of original report of jan 28 using data from jan 1 to jan 22
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1_22jan))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit1_22jan_preds = as.data.frame(emmeans(fit1_22jan, ~ collection_date_num, 
                                   # by="LABORATORY", 
                                   at=list(collection_date_num=seq(date.from,
                                                                   date.to)), 
                                   type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_22jan_preds$collection_date = as.Date(fit1_22jan_preds$collection_date_num, origin="1970-01-01")


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
# 27              18681 0.6957358 0.0141653 Inf 0.6673482 0.7228281      2021-02-23

# estimated share of 501Y.V1 among new infections (assuming time between infection & diagnosis of 7 days)
fit1_preds[fit1_preds$collection_date==(today+7),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 34               18688 0.8042867 0.01294345 Inf 0.7777467 0.8284802      2021-03-02


# estimated share of 501Y.V1 among currently diagnosed infections predicted based on jan 22 data
fit1_22jan_preds[fit1_22jan_preds$collection_date==today,]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 27              18681 0.8663892 0.0578731 Inf 0.7127089 0.9459225      2021-02-23

# estimated share of 501Y.V1 among new infections (assuming time between infection & diagnosis of 7 days) predicted based on jan 22 data
fit1_22jan_preds[fit1_22jan_preds$collection_date==(today+7),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 34               18688 0.9326176 0.03746546 Inf 0.8139173 0.9782164      2021-03-02


sum(tail(data_ag_byday_wide$est_n_B117, 14))/sum(tail(data_ag_byday_wide$n_pos,14)) 
# 49.95% of the samples of last 2 weeks in the dataset were estimated to be by British variant
# PS: with data 31/1 this was 15.4%  
# note: this is not the same as the estimated prop of the new infections or new diagnoses today that are of the British
# variant, which are much higher, see above)


# implied Re of wild type and 501Y.V1 given this predicted share of 501Y.V1 among all infections ####
# under a particular fitted transmission advantage
# based on the fact that the overall Re is a weighted average of the Re of the individual variants
# functions to calculate Re of wild type and of 501Y.V1 based on overall Re value and prop of positives that is 501Y.V1 propB117
# and transmission advantage of 501Y.V1 M
M_fitted = exp(fit1_emtrends[,c(2,5,6)]*4.7)[1,1] 
M_fitted # 1.52
Re_wild_type = function (Re, propB117, M=M_fitted) {
  Re / (1-propB117+M*propB117)
}
Re_B117 = function (Re, propB117, M=M_fitted) {
  M*Re / (1-propB117+M*propB117)
}

Re_cases = read.csv(paste0(".//data//",dat,"//Re_cases.csv")) 
Re_cases = read.csv(paste0(".//data//be_latest//Re_cases.csv")) 
# Re values calculated from instantaneous growth rate in nr of new cases 
# with instantaneous growth rate calculated as the first derivative (calculated using emtrends) to the GAM fit on new cases 
# gam(cbind(NEWCASES, totpop-NEWCASES) ~ s(DATE_NUM, bs="cs", k=32, fx=F) + 
#                                         WEEKDAY + s(log(TESTS_ALL), bs="cs", k=5, fx=F), family=binomial(cloglog), data=cases_tot) 
# and with Re calculated using R.from.r with gamma_mean=4.7, gamma_sd=2.9
# note that Re here is calculated at time of diagnosis
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

qplot(data=Re_cases, x=DATE, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.7), fill=I("darkgrey")) +
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
  ylab("Re at time of diagnosis")
ggsave(file=paste0(".//plots//",dat,"//Re_cases_Re_501YV1_Re_wildtype.png"), width=7, height=5)
Re_cases[Re_cases$DATE==max(Re_cases$DATE),]
# DATE_NUM          r          SE       df    r_LOWER    r_UPPER       DATE       Re Re_LOWER Re_UPPER collection_date_num  propB117
# 373    18694 0.02624037 0.001789314 321.1952 0.02272011 0.02976063 2021-03-08 1.128085 1.110345    1.146               18694 0.8798404
# Re_WT Re_WT_LOWER Re_WT_UPPER  Re_B117 Re_B117_LOWER Re_B117_UPPER
# 373 0.7709675   0.7588433   0.7832109 1.176857       1.15835      1.195546

Re_cases[Re_cases$DATE==today,]
# DATE_NUM          r          SE       df    r_LOWER    r_UPPER       DATE       Re Re_LOWER Re_UPPER collection_date_num  propB117
# 360    18681 0.02549623 0.001687455 321.1952 0.02217637 0.02881609 2021-02-23 1.124321 1.107621 1.141176               18681 0.7072377
# Re_WT Re_WT_LOWER Re_WT_UPPER  Re_B117 Re_B117_LOWER Re_B117_UPPER
# 360 0.8192741   0.8071049   0.8315563 1.250595       1.23202      1.269344




# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50%, 75% 90% by 501Y.V1 :
fit1_preds$collection_date[fit1_preds[,"prob"]>=0.5][1]-7 # >50% by 6th of February [5th Feb - 8 Feb] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.5][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.5][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.75][1]-7 # >75% by 20th of February [18 Feb - 22 Feb] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.75][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.75][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.9][1]-7 # >90% by 5th of March [2 March - 7 March] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.9][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.9][1]-7



# PLOT MODEL FIT common-slope model fit1

# plot for the whole of Belgium on response scale:

plot_fit1_response = qplot(data=fit1_preds, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~LABORATORY) +
  # geom_ribbon(data=fit1_22jan_preds, aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, 
  #                                        # fill=LABORATORY
  # ), 
  # fill=I("lightgrey"), 
  # alpha=I(1)) +
  # geom_line(data=fit1_22jan_preds, aes(y=100*prob, 
  #                                      # colour=LABORATORY
  # ), 
  # colour=I("grey"), 
  # alpha=I(1)) +
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

plot_fit1_responseB = qplot(data=fit1_preds, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~LABORATORY) +
  geom_ribbon(data=fit1_22jan_preds, aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, 
                  # fill=LABORATORY
  ), 
  fill=I("lightgrey"), 
  alpha=I(1)) +
  geom_line(data=fit1_22jan_preds, aes(y=100*prob, 
                # colour=LABORATORY
  ), 
  colour=I("grey"), 
  alpha=I(1)) +
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
plot_fit1_responseB


saveRDS(plot_fit1_responseB, file = paste0(".\\plots\\",dat,"\\Fig4B_fit1_binomGLMM_501YV1_Belgium_response scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig4B_fit1_binomGLMM_501YV1_Belgium_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4B_fit1_binomGLMM_501YV1_Belgium_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4B_fit1_binomGLMM_501YV1_Belgium_response scale.pdf"), width=8, height=6)


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
  scale_fill_manual("test outcome", values=c("darkgrey","grey40","blue"), labels=c("negative","positive (other strains)","S dropout (501Y.V1)")) +
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
  scale_fill_manual("test outcome", values=c("darkgrey","blue"), labels=c("positive (other strains)","S dropout (501Y.V1)")) +
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
date.to = today+30 
data_ag_long_subs = data_ag_long[(data_ag_long$LABORATORY %in% sel_labs)&(data_ag_long$collection_date>=date.from)&(data_ag_long$collection_date<=date.to),]
data_ag_long_subs$LABORATORY = droplevels(data_ag_long_subs$LABORATORY)

mfit0 = nnet::multinom(outcome ~ scale(collection_date_num) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mfit1 = nnet::multinom(outcome ~ scale(collection_date_num) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mfit2 = nnet::multinom(outcome ~ ns(collection_date_num, df=2) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit3 = nnet::multinom(outcome ~ ns(collection_date_num, df=3) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit4 = nnet::multinom(outcome ~ ns(collection_date_num, df=4) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit5 = nnet::multinom(outcome ~ ns(collection_date_num, df=5) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit6 = nnet::multinom(outcome ~ ns(collection_date_num, df=6) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
BIC(mfit0, mfit1, mfit2, mfit3, mfit4, mfit5, mfit6) # mfit3 fits best, df splines tuned based on BIC
#       df      BIC
# mfit0 16 150913.8
# mfit1 28 150894.7
# mfit2 18 149513.9
# mfit3 56 149434.8
summary(mfit5)

plot(Effect("collection_date_num",mfit5), style="stacked")
plot(Effect("collection_date_num",mfit5), confint=list(style="bands"), rug=FALSE)

# average growth rates of S-positive/wild type & S dropout/501Y.V1 cases evaluated today based on best fit multinomial model mfit4
emtrends(mfit5, ~outcome|1, var="collection_date_num",  at=list(collection_date_num=today_num), mode="latent")
# emmeans_1.5.3 output:
# outcome collection_date_num.trend      SE df lower.CL upper.CL
# n_neg                      0.0336 0.00636 72  0.02088   0.0462
# n_spos                    -0.0613 0.00966 72 -0.08053  -0.0420
# n_b117                     0.0277 0.00980 72  0.00818   0.047213
R.from.r(-0.0613) # Rt S pos / wild type = 0.74 [0.66-0.81]
R.from.r(-0.08053) # Rt S pos / wild type LCL = 0.66
R.from.r(-0.0420) # Rt S pos / wild type UCL = 0.81

R.from.r(0.0277) # Rt of S dropout = 1.14 [1.04-1.24]
R.from.r(0.00818) # Rt of S dropout = 1.04
R.from.r(0.047213) # Rt of S dropout = 1.24


R.from.r(0.0277)/R.from.r(-0.0613) # Rt of 501Y.V1 = 1.54x times higher than of wild type


# implied growth rate & transmission advantage
delta_r = data.frame(confint(contrast(emtrends(mfit4, ~outcome|1, var="collection_date_num",  
                                               at=list(outcome=c("n_spos","n_b117"),
                                                       collection_date_num=today_num), mode="latent"), method="revpairwise")))[,c(2,5,6)]
delta_r # growth advantage
#       estimate lower.CL  upper.CL
# 1 0.04360504 0.01757831 0.06963176
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

# saveRDS(multipanel_Re_growthadv_multinom, file = paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.pptx"), width=7, height=8)
# ggsave(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.png"), width=7, height=8)
# ggsave(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.pdf"), width=7, height=8)








# 5. SOME INTERNATIONAL COMPARISONS ####

# adapted from section 5. in https://github.com/nicholasdavies/newcovid/blob/master/multinomial_logistic_fits/multinomial%20logistic%20fits_FINAL.R
# cf. associated paper https://cmmid.github.io/topics/covid19/uk-novel-variant.html
# https://cmmid.github.io/topics/covid19/reports/uk-novel-variant/2021_02_06_Transmissibility_and_severity_of_VOC_202012_01_in_England_v2.pdf


# 5. INTERNATIONAL COMPARISONS: COMPETITIVE ADVANTAGE OF 501Y.V1 IN THE UK, DENMARK, SWITZERLAND & THE USA ####

# GIVEN THAT THE EFFECTIVE REPRODUCTION NUMBER R = (1 + k * r * g)^(1 / k) 
# (Park et al. 2020) WHEN GENERATION TIME IS GAMMA DISTRIBUTED (with mean g and k=(SD/g)^2), 
# WHICH IS APPROX EQUAL TO exp(r*g) WITH r=MALTHUSIAN GROWTH RATE, IT FOLLOWS THAT
# THE EXPECTED MULTIPLICATIVE DIFFERENCE IN THE R VALUE OF TWO COMPETING VARIANTS,
# ASSUMING IDENTICAL GENERATION TIMES, EQUALS exp((r_new-r_old)*g) = exp(delta_r*g),
# WHERE THE DIFFERENCE IN MALTHUSIAN GROWTH RATE delta_r IS SOMETIMES REFERRED TO AS
# THE SELECTION RATE (TRAVISANO & LENSKI 1996) AND delta_r*g IS THE DIMENSIONLESS
# SELECTION COEFFICIENT sT OF CHEVIN (2011).

M.from.delta_r = function (delta_r, g=4.7) { 
  delta_R = exp(delta_r*g)
  return( delta_R ) 
}

M.from.delta_r(0.088, 4.7)

# This function calculates the expected multiplicative effect on R M for 
# gamma distributed generation time gamma(mean=4.7d) (Nishiura et al. 2020)
# It works on an input dataframe df with delta_r values and returns the original
# data frame plus the estimate of M as a dataframe with extra columns
# with column names coln
M.from.delta_r_df = function (df, g=4.7, 
                              coln=c("M","M.LCL","M.UCL")) { 
  df_num = df[,which(unlist(lapply(df, is.numeric))), drop=F]
  df_nonnum = df[,which(!unlist(lapply(df, is.numeric))), drop=F]
  df_out1 = apply(df_num, 2, function (delta_r) M.from.delta_r(delta_r, g))
  if (class(df_out1)[1]=="numeric") df_out1=as.data.frame(t(df_out1), check.names=F)
  df_out = data.frame(df_out1, check.names=F)
  if (!is.null(coln)) colnames(df_out) = coln
  return( data.frame(df_nonnum, df_num, df_out, check.names=F) )
}


# 5.1. DATA UK : ANALYSIS OF PILLAR 2 S-GENE TARGET FAILURE DATA ####

levels_UKregions = c("South East","London","East of England",
                     "South West","Midlands","North East and Yorkshire",
                     "Scotland","North West","Wales")

# Pillar 2 S gene targeted failure data (SGTF) (S dropout)
sgtfdata_uk = read.csv("https://github.com/nicholasdavies/newcovid/raw/master/fitting_data/sgtf-2021-01-18.csv") 
sgtfdata_uk$other = sgtfdata_uk$other+sgtfdata_uk$sgtf
colnames(sgtfdata_uk) = c("collection_date","REGION","SGTF","TOTAL")
# modelled proportion of S dropout that was actually the VOC
sgtfdata_uk_truepos = read.csv("https://github.com/nicholasdavies/newcovid/raw/master/data/sgtfvoc.csv") 
sgtfdata_uk$TRUEPOS = sgtfdata_uk_truepos$sgtfv[match(interaction(sgtfdata_uk$REGION, sgtfdata_uk$collection_date),
                                                      interaction(sgtfdata_uk_truepos$nhs_name, sgtfdata_uk_truepos$date))] # modelled proportion of S dropout samples that were actually the VOC
sgtfdata_uk$est_n_B117 = sgtfdata_uk$SGTF * sgtfdata_uk$TRUEPOS
sgtfdata_uk$COUNTRY = "UK"
sgtfdata_uk = sgtfdata_uk[,c("collection_date","COUNTRY","REGION","est_n_B117","TOTAL")]
colnames(sgtfdata_uk)[which(colnames(sgtfdata_uk)=="TOTAL")] = "n_pos"
range(sgtfdata_uk$collection_date) # "2020-10-01" "2021-01-17"
sgtfdata_uk$collection_date = as.Date(sgtfdata_uk$collection_date)
sgtfdata_uk$collection_date_num = as.numeric(sgtfdata_uk$collection_date)
sgtfdata_uk$REGION = factor(sgtfdata_uk$REGION, levels=levels_UKregions)
sgtfdata_uk$REGION = droplevels(sgtfdata_uk$REGION)
sgtfdata_uk$obs = factor(1:nrow(sgtfdata_uk))
sgtfdata_uk$propB117 = sgtfdata_uk$est_n_B117 / sgtfdata_uk$n_pos
head(sgtfdata_uk)

set_sum_contrasts()
glmersettings = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e5)) # bobyqa, PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B"))
glmersettings3 = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb"))
glmersettings4 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
fit_ukSGTF_1 = glmer(cbind(est_n_B117, n_pos-est_n_B117 ) ~ (1|obs)+scale(collection_date_num)+REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings)  # common slope model
fit_ukSGTF_2 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)*REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings3) # heter slope model
fit_ukSGTF_3 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(ns(collection_date_num,df=3))+REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings3) # with additive spline term
fit_ukSGTF_4 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+ns(collection_date_num,df=3)*REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings3) # with spline term in interaction with region
BIC(fit_ukSGTF_1, fit_ukSGTF_2, fit_ukSGTF_3, fit_ukSGTF_4) 
# separate-slopes 3 df spline model fit_be_uk2_4 best
# df      BIC
# fit_ukSGTF_1  9 4902.696
# fit_ukSGTF_2 15 4769.405
# fit_ukSGTF_3 11 4905.592
# fit_ukSGTF_4 29 4428.474


# model fit_ukSGTF_4 best

summary(fit_ukSGTF_4)

# GROWTH RATE AND TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model fit_ukSGTF_4, we get
fit_ukSGTF_4_growthrates_avg_model2h = as.data.frame(emtrends(fit_ukSGTF_4, ~ 1, var="collection_date_num",
                                                              at=list(sample_date_num=as.numeric(seq(as.Date("2020-11-01"),
                                                                                                     max(sgtfdata_uk$collection_date), by=1)))))[,-c(3,4)] 
colnames(fit_ukSGTF_4_growthrates_avg_model2h)[2] = "logistic_growth_rate"
fit_ukSGTF_4_growthrates_avg_model2h = M.from.delta_r_df(fit_ukSGTF_4_growthrates_avg_model2h)
fit_ukSGTF_4_growthrates_avg_model2h
# 1         logistic_growth_rate asymp.LCL asymp.UCL        M    M.LCL    M.UCL
# 1 overall            0.1093807 0.1074623  0.111299 1.672115 1.657106 1.687259

# growth rates per region for model fit_ukSGTF_4
fit_ukSGTF_4_growthrates_region_model2h = as.data.frame(emtrends(fit_ukSGTF_4, ~ REGION, var="collection_date_num",
                                                                 at=list(sample_date_num=as.numeric(seq(as.Date("2020-11-01"),
                                                                                                        max(sgtfdata_uk$collection_date), by=1)))))[,-c(3,4)] 
colnames(fit_ukSGTF_4_growthrates_region_model2h)[2] = "logistic_growth_rate"
fit_ukSGTF_4_growthrates_region_model2h = M.from.delta_r_df(fit_ukSGTF_4_growthrates_region_model2h)
fit_ukSGTF_4_growthrates_region_model2h
#                     REGION logistic_growth_rate  asymp.LCL asymp.UCL        M    M.LCL    M.UCL
# 1               South East           0.09898888 0.09592910 0.1020487 1.592409 1.569672 1.615474
# 2                   London           0.11321265 0.11007799 0.1163473 1.702503 1.677604 1.727771
# 3          East of England           0.11927245 0.11526389 0.1232810 1.751689 1.718996 1.785004
# 4               South West           0.11793937 0.10916388 0.1267149 1.740748 1.670412 1.814046
# 5                 Midlands           0.09964142 0.09541678 0.1038661 1.597300 1.565897 1.629333
# 6 North East and Yorkshire           0.10880077 0.10423148 0.1133700 1.667564 1.632133 1.703763
# 7               North West           0.10780908 0.10264109 0.1129771 1.659809 1.619979 1.700619


# PLOT MODEL FIT

# spline model fit_ukSGTF_4
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-04-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_ukSGTF_4))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_ukSGTF_4_preds = as.data.frame(emmeans(fit_ukSGTF_4, ~ collection_date_num, 
                                           by=c("REGION"), 
                                           at=list(collection_date_num=seq(date.from,
                                                                           date.to)), 
                                           type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_ukSGTF_4_preds$collection_date = as.Date(fit_ukSGTF_4_preds$collection_date_num, origin="1970-01-01")

n = length(levels(fit_ukSGTF_4_preds$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
# reg_cols[2:n] = rev(reg_cols[2:n])

fit_ukSGTF_4_preds$REGION = factor(fit_ukSGTF_4_preds$REGION, levels=unique(fit_ukSGTF_4_preds$REGION))
sgtfdata_uk$REGION = factor(sgtfdata_uk$REGION, levels=levels(fit_ukSGTF_4_preds$REGION))

# PLOT MODEL FIT (logit scale):
plot_UK_SGTF = qplot(data=fit_ukSGTF_4_preds, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~COUNTRY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("UK") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0.001,99.9), expand=c(0,0)) 
plot_UK_SGTF


# PLOT MODEL FIT (response scale):
plot_UK_SGTF_response = qplot(data=fit_ukSGTF_4_preds, x=collection_date, y=prob*100, geom="blank") +
  # facet_wrap(~COUNTRY) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117*100, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("UK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_UK_SGTF_response


# plot model fit (response scale) (data UK + Belgium combined)
fit1_preds$REGION = "Belgium"
preds_UK_plus_BE = rbind(fit_ukSGTF_4_preds, fit1_preds)

plot_UK_SGTF_BE_response = qplot(data=preds_UK_plus_BE, x=collection_date, y=prob*100, geom="blank") +
  # facet_wrap(~COUNTRY) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=c(reg_cols,"steelblue")) +
  scale_fill_manual("", values=c(reg_cols, "steelblue")) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117*100, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  geom_point(data=data_ag_byday_wide, 
             aes(x=collection_date, y=propB117*100, size=n_pos
             ), 
             colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("SPREAD OF VARIANT 501Y.V1 IN THE\nUK & BELGIUM BASED ON S DROPOUT DATA") +
  theme(plot.title = element_text(hjust = 0.5))
plot_UK_SGTF_BE_response
saveRDS(plot_UK_SGTF_BE_response, file = paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.pdf"), width=8, height=6)


# 5.2. DATA DENMARK: SEQUENCING DATA ####

# Data source: Danish Covid-19 Genome Consortium & the Statens Serum Institut, https://www.covid19genomics.dk/statistics

data_denmark = read.csv(".//data//dk//data_denmark_20210224.csv", sep=";", dec=",")
data_denmark$percent = NULL
data_denmark$Region = gsub("SjÃ¦lland","Sjælland",data_denmark$Region)
data_denmark$WEEK = sapply(data_denmark$Week, function(s) as.numeric(strsplit(s, "W")[[1]][[2]]))
data_denmark$date = as.Date(NA)
data_denmark$date[data_denmark$WEEK>=42] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK>=42] - 1 ) + 1
data_denmark$date[data_denmark$WEEK<42] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK<42] - 1 ) + 6 
data_denmark$date_num = as.numeric(data_denmark$date)
data_denmark$obs = factor(1:nrow(data_denmark))
colnames(data_denmark)[colnames(data_denmark) %in% c("yes")] = "n_B117"
data_denmark$propB117 = data_denmark$n_B117 / data_denmark$total

data_denmark_whole = data_denmark[data_denmark$Region=="Whole Denmark",]
data_denmark = data_denmark[data_denmark$Region!="Whole Denmark",]
levels_DK = c("Syddanmark","Sjælland","Nordjylland","Hovedstaden","Midtjylland")
data_denmark$Region = factor(data_denmark$Region, levels=levels_DK)
range(data_denmark$date) # "2020-11-5" "2021-02-28"

fit_denmark1 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + Region + scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark2 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + Region * scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark3 = glmer(cbind(n_B117,total-n_B117) ~ (1|Region/obs) + scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark4 = glmer(cbind(n_B117,total-n_B117) ~ (date_num||Region/obs) + scale(date_num), family=binomial(logit), data=data_denmark)
BIC(fit_denmark1, fit_denmark2, fit_denmark3, fit_denmark4)
#              df      BIC
# fit_denmark1  7 541.6427
# fit_denmark2 11 505.6540
# fit_denmark3  4 533.2144
# fit_denmark4  6 541.9784

summary(fit_denmark2)

# common-slope model fit_denmark2 with nested random intercepts fits best

#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model fit_denmark2, we get
dk_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_denmark2, ~ 1, var="date_num"))[,-c(3,4)] 
colnames(dk_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
dk_growthrates_avg_B117vsallother = M.from.delta_r_df(dk_growthrates_avg_B117vsallother)
dk_growthrates_avg_B117vsallother
# 1 logistic_growth_rate  asymp.LCL  asymp.UCL        M    M.LCL    M.UCL
# 1 overall            0.0776125 0.07254026 0.08268473 1.440195 1.406268 1.474941


# PLOT MODEL FIT
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-04-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_denmark3))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_denmark_preds = as.data.frame(emmeans(fit_denmark2, ~ date_num, 
                                          at=list(date_num=seq(date.from,
                                                               date.to)), 
                                          type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_denmark_preds$date = as.Date(fit_denmark_preds$date_num, origin="1970-01-01")

# n = length(levels(fit_denmark_preds$Region))
# reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]

# PLOT MODEL FIT (response scale)
plot_denmark = qplot(data=fit_denmark_preds, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL # , 
                  # fill=Region
  ), 
  fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob# , 
                # colour=Region
  ), 
  colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-10-01"),as.Date("2021-04-01")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_denmark_whole, 
             aes(x=date, y=propB117, size=total,
                 # colour=Region
             ), 
             colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\nsequences", trans="sqrt", 
                        range=c(1, 4), limits=c(1,max(data_denmark_whole$total)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("DENMARK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_denmark


# PLOT MODEL FIT (response scale)
plot_denmark_response = qplot(data=fit_denmark_preds, x=date, y=prob*100, geom="blank") +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL # , 
                  # fill=Region
  ), 
  fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100 # , 
                # colour=Region
  ), 
  colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-10-01"),as.Date("2021-04-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  # scale_color_manual("", values=reg_cols) +
  # scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_denmark_whole, 
             aes(x=date, y=propB117*100, size=total # ,
                 # colour=Region
             ), 
             colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,max(data_denmark_whole$total)), 
                        breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("DENMARK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_denmark_response



# 5.3. DATA SWITZERLAND : SEQUENCING & RT-PCR RE-SCREENING DATA ####

# Data source: https://ispmbern.github.io/covid-19/variants (contact: Christian Althaus) 
# & https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv (https://ibz-shiny.ethz.ch/covidDashboard/variant-plot/index.html, contact: Tanja Stadler)

data_geneva = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_GE.csv")
data_geneva$date = as.Date(data_geneva$date)
data_geneva$lab = "Geneva"
colnames(data_geneva)[colnames(data_geneva) %in% c("N501Y")] = c("n_B117")
head(data_geneva)
data_zurich = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_ZH.csv")
data_zurich$date = as.Date(data_zurich$date)
data_zurich$lab = "Zürich"
colnames(data_zurich)[colnames(data_zurich) %in% c("N501Y")] = c("n_B117")
head(data_zurich)
data_bern = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_BE.csv")
data_bern$date = as.Date(data_bern$date)
data_bern$lab = "Bern"
colnames(data_bern)[colnames(data_bern) %in% c("N501Y")] = c("n_B117")
head(data_bern)

data_viollier_risch = read.csv("https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv")
data_viollier_risch[is.na(data_viollier_risch)] = 0
data_viollier_risch$date = as.Date(NA)
data_viollier_risch$date[data_viollier_risch$week>=51] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_viollier_risch$week[data_viollier_risch$week>=51] - 1 ) + 1
data_viollier_risch$date[data_viollier_risch$week<51] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_viollier_risch$week[data_viollier_risch$week<51] - 1 ) + 6 # PS dates were made to match the ones given in https://ispmbern.github.io/covid-19/variants/data/variants_CH.csv
colnames(data_viollier_risch)[colnames(data_viollier_risch) %in% c("n","b117")] = c("total","n_B117")
data_viollier_risch = data_viollier_risch[,c("date","total","n_B117","lab")]

data_switzerland = rbind(data_geneva, data_zurich, data_bern, data_viollier_risch)[,c("date","lab","n_B117","total")]
# write_csv(data_switzerland, file=".//data//ch//data_switzerland_20210224.csv")

# Details data:
# Viollier data = sequencing of a random subset of all positive cases by ETH/Tanja Stadler (covers large parts of Switzerland, though with a bias towards German speaking Switzerland) - as it is sequencing data is 1-2 weeks later than N501Y screening
# Risch - Taqpath + N501Y re-screening = faster (covers primarily German speaking Switzerland)
# Samples are provided and screened by Labor Risch. Genomic characterization is performed by Labor Risch, the University Hospital Basel (Clinical Mircobiology) and the University Hospitals of Geneva (Group Eckerle and Group Kaiser). 
# Geneva - centre de reference pour infections virales emergentes / university hospital Geneva - N501Y and WGS currently:
# Samples that were sent to the Geneva University Hospitals for primary diagnosis of SARS-CoV-2. All positives were re-screened for 501Y using RT-PCR (mostly B.1.1.7). To cover the period of November and December 2020, we use sequence data from randomly chosen samples from Geneva that were submitted to GISAID by the Swiss Viollier Sequencing Consortium from ETH Zurich.
# Bern: Samples from SARS-CoV-2-positive cases that were re-screened for 501Y using RT-PCR at the Institute for Infectious Diseases, University of Bern.
# Zurich: Samples from SARS-CoV-2-positive cases from the University Hospital Zurich and test centers at Limmattal Hospital in Schlieren (ZH) and Spital Männedorf that were re-screened for 501Y using RT-PCR at the Institute of Medical Virology, University of Zurich. In addition, we use SARS-CoV-2-positive samples from Kantonsspital Winterthur and its walk-in test center that were re-screened for 501Y using RT-PCR.

# data_switzerland = read_csv(file=".//multinomial_logistic_fits//data//ch//data_switzerland_20210216.csv", col_names=TRUE) 
data_switzerland = data.frame(data_switzerland)
data_switzerland$date = as.Date(data_switzerland$date)
data_switzerland$lab = factor(data_switzerland$lab, levels=c("Geneva","Zürich","Bern","Viollier","Risch"),
                              labels=c("Geneva","Zürich","Bern","Switzerland","Switzerland"))
data_switzerland$date_num = as.numeric(data_switzerland$date)
data_switzerland$obs = factor(1:nrow(data_switzerland))
data_switzerland$propB117 = data_switzerland$n_B117 / data_switzerland$total
range(data_switzerland$date) # "2020-11-02" "2021-02-18"

fit_switerland1 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + lab + scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland2 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + lab * scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland3 = glmer(cbind(n_B117,total-n_B117) ~ (1|lab/obs) + scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland4 = glmer(cbind(n_B117,total-n_B117) ~ (date_num||lab/obs) + scale(date_num), family=binomial(logit), data=data_switzerland)
BIC(fit_switerland1, fit_switerland2, fit_switerland3, fit_switerland4)
#                 df      BIC
# fit_switerland1  6 646.8101
# fit_switerland2  9 659.2752
# fit_switerland3  4 656.9695
# fit_switerland4  6 666.7944

# fit fit_switerland1 best

summary(fit_switerland1)


#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model fit_switerland1, we get
ch_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_switerland1, ~ 1, var="date_num"))[,-c(3,4)] 
colnames(ch_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
ch_growthrates_avg_B117vsallother = M.from.delta_r_df(ch_growthrates_avg_B117vsallother)
ch_growthrates_avg_B117vsallother
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL     M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall            0.09292921 0.08611969 0.09973872 1.547696 1.498947 1.598031



# PLOT MODEL FIT
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-04-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_switerland1))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_switzerland_preds = as.data.frame(emmeans(fit_switerland1, ~ date_num, 
                                              by=c("lab"), 
                                              at=list(date_num=seq(date.from,
                                                                   date.to)), 
                                              type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_switzerland_preds$date = as.Date(fit_switzerland_preds$date_num, origin="1970-01-01")

n = length(levels(fit_switzerland_preds$lab))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
# reg_cols[2:n] = rev(reg_cols[2:n])

# PLOT MODEL FIT (logit scale):
plot_switzerland = qplot(data=fit_switzerland_preds, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=lab
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=lab
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-11-01"),as.Date("2021-04-01")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_switzerland, 
             aes(x=date, y=propB117, size=total,
                 colour=lab
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("SWITZERLAND") +
  theme(plot.title = element_text(hjust = 0.5))
plot_switzerland


# PLOT MODEL FIT (response scale):
plot_switzerland_response = qplot(data=fit_switzerland_preds, x=date, y=prob*100, geom="blank") +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=lab
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=lab
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-11-01"),as.Date("2021-03-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_switzerland, 
             aes(x=date, y=propB117*100, size=total,
                 colour=lab
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(500,1000,2000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("SWITZERLAND") +
  theme(plot.title = element_text(hjust = 0.5))
plot_switzerland_response




# 5.3. DATA USA : S-GENE TARGET FAILURE DATA ####

# Data source: Helix® COVID-19 Surveillance, https://github.com/myhelix/helix-covid19db
# see preprint https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1 & https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data

us_data = read.csv("https://github.com/myhelix/helix-covid19db/raw/master/counts_by_state.csv")
# write.csv(us_data, file=".//data//us//data_us_20210224.csv", row.names=F)

us_data$collection_date = as.Date(us_data$collection_date)
us_data$collection_date_num = as.numeric(us_data$collection_date)
us_data$obs = factor(1:nrow(us_data))
# us_data = us_data[us_data$state %in% sel_states,]
us_data$state = factor(us_data$state)

range(us_data$collection_date) # "2020-09-05" "2021-02-21"

fit_us_propB117amongSGTF = glmer(cbind(B117, sequenced_SGTF-B117) ~ (1|state)+scale(collection_date_num), 
                                 family=binomial(logit), data=us_data)

# implied growth rate advantage of B.1.1.7 over other earlier strains showing S dropout:
as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)]
#   collection_date_num.trend  asymp.LCL asymp.UCL
# 1                0.09847144 0.08327195 0.1136709

# with a generation time of 4.7 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of B.1.1.7 over other strains showing S dropout of
exp(4.7*as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.588541  1.479018  1.706174


# FIT FOR WHOLE US + PLOT

fitted_truepos = predict(fit_us_propB117amongSGTF, newdat=us_data, type="response") 
# fitted true positive rate, ie prop of S dropout samples that are B.1.1.7 for dates & states in helix_sgtf

us_data$est_n_B117 = us_data$all_SGTF*fitted_truepos # estimated nr of B.1.1.7 samples
us_data$propB117 = us_data$est_n_B117/us_data$positive
fit_us1 = glmer(cbind(est_n_B117, positive-est_n_B117) ~ (1|state/obs)+scale(collection_date_num), 
                family=binomial(logit), data=us_data) # random intercepts by state
fit_us2 = glmer(cbind(est_n_B117, positive-est_n_B117) ~ (collection_date_num||state/obs)+scale(collection_date_num), 
                family=binomial(logit), data=us_data) # random intercepts+slopes by state, with uncorrelated intercepts & slopes
BIC(fit_us1, fit_us2) # random intercept model fit_us1 is best
# df      BIC
# fit_us1  4 2253.782
# fit_us2  6 2269.764
summary(fit_us1)

#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all states, using the most parsimonious model fit_us1, we get
us_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_us1, ~ 1, var="collection_date_num"))[,-c(3,4)] 
colnames(us_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
us_growthrates_avg_B117vsallother = M.from.delta_r_df(us_growthrates_avg_B117vsallother)
us_growthrates_avg_B117vsallother
# 1 logistic_growth_rate asymp.LCL  asymp.UCL        M  M.LCL    M.UCL
# 1 overall           0.08399167 0.0805767 0.08740664 1.484029 1.4604 1.508041


# plot model fit fit_us

date.to = as.numeric(as.Date("2021-06-01"))
# sel_states = intersect(rownames(ranef(fit_us)$state)[order(ranef(fit_us1)$state[,1], decreasing=T)],states_gt_500)[1:16] # unique(helix_sgtf$state[helix_sgtf$propB117>0.03])
# rem_states = c("NY","NJ","MN","IL","AL","OH","MI") # states with too few data points we don't want to show on plot
# sel_states = setdiff(sel_states,rem_states)


# sel_states = unique(us_data$state)
# we fitted our model on all the available data from all states, but below we will plot just
# the 9 states with the most data
# sel_states=c("FL","NY","CA","NJ","GA","TX","OH","PA","LA","IL","MI","MA","NC","IN","AZ")
# sel_states=c("FL","CA","GA","TX","PA","LA","IL","MI","MA","NC","IN","AZ")
sel_states=c("FL","CA","GA","TX","PA","MA","NC","IN","AZ")
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_us1))$sdcor, function (x) x^2))) 
fit_us_preds = as.data.frame(emmeans(fit_us1, ~ collection_date_num, 
                                     # by="state", 
                                     at=list(collection_date_num=seq(min(us_data$collection_date_num),
                                                                     date.to)), 
                                     type="link"), bias.adjust = TRUE, sigma = total.SD)
fit_us_preds$collection_date = as.Date(fit_us_preds$collection_date_num, origin="1970-01-01")
fit_us_preds2 = do.call(rbind,lapply(unique(us_data$state), function(st) { ranintercs = ranef(fit_us1)$state
raninterc = ranintercs[rownames(ranintercs)==st,]
data.frame(state=st, fit_us_preds, raninterc=raninterc)}))
fit_us_preds2$prob = plogis(fit_us_preds2$emmean+fit_us_preds2$raninterc)
fit_us_preds2$prob.asymp.LCL = plogis(fit_us_preds2$asymp.LCL+fit_us_preds2$raninterc)
fit_us_preds2$prob.asymp.UCL = plogis(fit_us_preds2$asymp.UCL+fit_us_preds2$raninterc)
fit_us_preds2 = fit_us_preds2[as.character(fit_us_preds2$state) %in% sel_states,]
fit_us_preds2$state = droplevels(fit_us_preds2$state)
fit_us_preds2$state = factor(fit_us_preds2$state, # we order states by random intercept, ie date of introduction
                             levels=intersect(rownames(ranef(fit_us1)$state)[order(ranef(fit_us1)$state[,1], decreasing=T)],
                                              sel_states))

# PLOT MODEL FIT (logit scale)
plot_us = qplot(data=fit_us_preds2, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~state, nrow=3) +
  geom_ribbon(aes(y=prob, ymin=prob.asymp.LCL, ymax=prob.asymp.UCL, colour=NULL, 
                  fill=state
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=state
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_us_preds$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.001,0.9990001), expand=c(0,0)) +
  scale_color_discrete("state", h=c(0, 240), c=180, l=55) +
  scale_fill_discrete("state", h=c(0, 240), c=180, l=55) +
  geom_point(data=us_data[us_data$state %in% sel_states,],  
             aes(x=collection_date, y=propB117, size=positive,
                 colour=state
             ), pch=I(16),
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="sqrt", 
                        range=c(1, 2), limits=c(1,max(us_data$positive)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("") # +
# theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
# ggtitle("US") +
# theme(plot.title = element_text(hjust = 0.5))
plot_us


# PLOT MODEL FIT (response scale)
plot_us_response = qplot(data=fit_us_preds2, x=collection_date, y=prob*100, geom="blank") +
  facet_wrap(~state) +
  geom_ribbon(aes(y=prob*100, ymin=prob.asymp.LCL*100, ymax=prob.asymp.UCL*100, colour=NULL, 
                  fill=state
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=state
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_us_preds2$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("state", h=c(0, 240), c=180, l=55) +
  scale_fill_discrete("state", h=c(0, 240), c=180, l=55) +
  geom_point(data=us_data[us_data$state %in% sel_states,],  
             aes(x=collection_date, y=propB117*100, size=positive,
                 colour=state
             ), pch=I(16),
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 2), limits=c(1,max(us_data$positive)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("") # +
# theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
# ggtitle("US") +
# theme(plot.title = element_text(hjust = 0.5))

plot_us_response






# 5.4. MULTIPANEL PLOT INTERNATIONAL COMPARISONS ####

fit_uk_preds2 = fit_ukSGTF_4_preds
fit_uk_preds2$country = "UK"
colnames(fit_uk_preds2)[2] = "REGION"
colnames(fit_uk_preds2)[1] = "date_num"
colnames(fit_uk_preds2)[8] = "date"
fit_switzerland_preds2 = fit_switzerland_preds
fit_switzerland_preds2$country = "Switzerland"
colnames(fit_switzerland_preds2)[2] = "REGION"
colnames(fit_switzerland_preds2)[1] = "date_num"
colnames(fit_switzerland_preds2)[8] = "date"
fit_denmark_preds2 = fit_denmark_preds
fit_denmark_preds2$country = "Denmark"
fit_denmark_preds2$REGION = "Denmark"
colnames(fit_denmark_preds2)[1] = "date_num"
colnames(fit_denmark_preds2)[7] = "date"
fit_us_preds3 = fit_us_preds2
fit_us_preds3$country = "USA"
fit_us_preds3 = fit_us_preds3[,-which(colnames(fit_us_preds3) %in% c("asymp.LCL","asymp.UCL"))]
colnames(fit_us_preds3)[1] = "REGION"
colnames(fit_us_preds3)[2] = "date_num"
colnames(fit_us_preds3)[6] = "date"
colnames(fit_us_preds3)[9] = "asymp.LCL"
colnames(fit_us_preds3)[10] = "asymp.UCL"
fit_us_preds3 = fit_us_preds3[fit_us_preds3$REGION %in% c("FL","CA"),]
fit_us_preds3$REGION = factor(fit_us_preds3$REGION, levels=c("FL","CA"), labels=c("Florida","California"))
fit_us_preds3 = fit_us_preds3[,c("date_num","REGION","prob","SE","df","asymp.LCL","asymp.UCL","date","country")]
fit_be_preds = fit1_preds
fit_be_preds$country = "Belgium"
colnames(fit_be_preds)[1] = "date_num"
colnames(fit_be_preds)[7] = "date"

fits_international = rbind(fit_uk_preds2,fit_switzerland_preds2,fit_us_preds3,
                           fit_denmark_preds2,fit_be_preds)
fits_international$country = factor(fits_international$country, levels=c("UK","Switzerland","USA","Denmark","Belgium"))

sgtfdata_uk2 = sgtfdata_uk
sgtfdata_uk2$country = "UK"
colnames(sgtfdata_uk2)[colnames(sgtfdata_uk2) %in% c("collection_date","n_pos")] = c("date","total")
sgtfdata_uk2 = sgtfdata_uk2[,c("date","country","REGION","propB117","total")]

data_switzerland2 = data_switzerland
data_switzerland2$country = "Switzerland"
colnames(data_switzerland2)[colnames(data_switzerland2) %in% c("lab")] = c("REGION")
data_switzerland2 = data_switzerland2[,c("date","country","REGION","propB117","total")]

data_denmark2 = data_denmark_whole
data_denmark2$country = "Denmark"
data_denmark2$REGION = "Denmark"
data_denmark2 = data_denmark2[,c("date","country","REGION","propB117","total")]

data_us2 = data.frame(us_data)
data_us2$country = "USA"
colnames(data_us2)[1] = "REGION"
colnames(data_us2)[2] = "date"
colnames(data_us2)[3] = "total"
data_us2 = data_us2[,c("date","country","REGION","propB117","total")]
data_us2 = data_us2[data_us2$REGION %in% c("FL","CA"),]
data_us2$REGION = factor(data_us2$REGION, levels=c("FL","CA"), labels=c("Florida","California"))

data_belgium = data_ag_byday_wide
data_belgium$country = "Belgium"
data_belgium$REGION = "Belgium"
colnames(data_belgium)[1] = "date"
data_belgium = data_belgium[,c("date","country","REGION","propB117","total")]

data_international = rbind(sgtfdata_uk2, data_switzerland2, data_us2,
                           data_denmark2, data_belgium)
data_international$country = factor(data_international$country, levels=c("UK","Switzerland","USA","Denmark","Belgium"))

# n1 = length(levels(fit_uk_preds2$REGION))
# n2 = length(levels(fit_switzerland_preds2$REGION))
# n3 = length(levels(fit_denmark_preds2$REGION))
# reg_cols = c(hcl(h = seq(290, 0, length = n1), l = 50, c = 255),
#              muted(hcl(h = seq(290, 0, length = n2+n3), l = 50, c = 255), c=200, l=40))

# ymin = 0.001
ymax = 0.999
data_international$propB117[data_international$propB117>ymax] = ymax
fits_international$prob[fits_international$prob>ymax] = ymax
fits_international$asymp.LCL[fits_international$asymp.LCL>ymax] = ymax
fits_international$asymp.UCL[fits_international$asymp.UCL>ymax] = ymax

fits_international$REGION = factor(fits_international$REGION, levels=levels(fits_international$REGION))
data_international$REGION = factor(data_international$REGION, levels=levels(fits_international$REGION))

# PLOT MODEL FITS (response scale)
plot_international = qplot(data=fits_international, x=date, y=prob, geom="blank") +
  facet_wrap(~country, nrow=2, scales="fixed") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9") # ,
                      # limits = c(ymin,ymax+1E-7)
  ) +
  # scale_color_manual("", values=reg_cols) +
  # scale_fill_manual("", values=reg_cols) +
  scale_color_discrete("region", h=c(0, 290), c=180, l=55) +
  scale_fill_discrete("region", h=c(0, 290), c=180, l=55) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_international, 
             aes(x=date, y=propB117, size=total, # shape=country,
                 colour=REGION, fill=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 2), limits=c(1,max(data_international$total)), breaks=c(100,1000,10000)) +
  # scale_shape_manual(values=21:25) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("") +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2),
    fill = guide_legend(order = 2),
    size = guide_legend(order = 3)
  ) + 
  coord_cartesian( 
    xlim=c(as.Date("2020-09-01"),as.Date("2021-03-31")),
    ylim=c(ymin,ymax+1E-7), 
    expand=FALSE) 
# ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
# theme(plot.title = element_text(hjust = 0.5))
plot_international

saveRDS(plot_international, file = paste0(".\\plots\\",dat,"\\Fig7_international_data_UK_CH_USA_DK_BE.rds"))
graph2ppt(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE.pptx"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE.png"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE.pdf"), width=9, height=7)




# PLOT MODEL FITS (response scale)
plot_international_response = qplot(data=fits_international, x=date, y=prob*100, geom="blank") +
  facet_wrap(~country, nrow=2) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + 
  xlab("") + 
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_color_discrete("region", h=c(0, 290), c=180, l=55) +
  scale_fill_discrete("region", h=c(0, 290), c=180, l=55) +
  #   scale_color_manual("", values=reg_cols) +
  #  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_international, 
             aes(x=date, y=propB117*100, size=total, # shape=country,
                 colour=REGION, fill=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="identity", 
                        range=c(1, 2), limits=c(1,max(data_international$total)), breaks=c(100,1000,10000)) +
  # scale_shape_manual(values=21:25) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2),
    fill = guide_legend(order = 2),
    size = guide_legend(order = 3)
  ) +
  coord_cartesian( 
    xlim=c(as.Date("2020-09-01"),as.Date("2021-04-01")-1),
    ylim=c(0,100), expand=c(0,0))
# +
# ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
# theme(plot.title = element_text(hjust = 0.5))
plot_international_response

saveRDS(plot_international_response, file = paste0(".\\plots\\",dat,"\\Fig7_international_data_UK_CH_USA_DK_BE_response.rds"))
graph2ppt(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE_response.pptx"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE_response.png"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE_response.pdf"), width=9, height=7)


plot_us2 = plot_us + coord_cartesian(xlim=c(as.Date("2020-11-01"), as.Date("2021-03-31")),
                                     ylim=c(0.001,99.9), expand=c(0,0)) # + ggtitle("SPREAD OF VARIANT B.1.1.7 IN THE US")
plot_us2

saveRDS(plot_us2, file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.rds"))
graph2ppt(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.pptx"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.png"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.pdf"), width=9, height=7)

