# ANALYSIS OF S-GENE TARGET FAILURE (S DROPOUT) DATA FROM BELGIUM TO INFER CONTAGIOUSNESS OF NEW VARIANT OF CONCERN B.1.1.7 / 501Y.V1 ####
# PLUS INTERNATIONAL COMPARISON (USING DATA FROM THE UK, DENMARK, SWITZERLAND & THE US)
# AND ASSESSMENT OF GROWTH ADVANTAGE OF THE SOUTH AFRICAN VOC 501Y.V2 & BRAZILIAN VOC 501Y.V3 BASED ON BASELINE SURVEILLANCE SEQUENCING DATA
# Tom Wenseleers
# All Belgian data provided by Emmanuel André

# Data provided by Emmanuel André (BE), COG-UK, PHE & N. Davies (UK), 
# Statens Serum Institut & Danish Covid-19 Genome Consortium (DK, https://www.covid19genomics.dk/statistics), 
# Christian Althaus, Swiss Viollier Sequencing Consortium, Institute of Medical Virology, University of Zurich, 
# Swiss National Covid-19 Science Task Force (Switzerland, https://ispmbern.github.io/covid-19/variants/, 
# https://ispmbern.github.io/covid-19/variants/data & https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv)
# and Helix, San Mateo, CA, Karthik Gangavarapu & Kristian G. Andersen (US, https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data, https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1)

# last update 30 APRIL 2021

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


dat="2021_04_30" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-04-30")
today_num = as.numeric(today)
today # "2021-04-30"

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOC 501Y.V1,VOC 501Y.V2&VOC 501Y.V3 IN BELGIUM BASED ON BASELINE SURVEILLANCE SEQUENCING DATA ####
# (baseline surveillance sequencing results, i.e. randomly sampled)
# data from weekly Sciensano report of 9/4/2021

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
range(be_seqdata$collection_date) # "2020-12-03" "2021-04-22"


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

baseline_surveillance = ggplot(data=be_basseqdata_long[be_basseqdata_long$variant!="501Y.V1+V2+V3",], 
                             aes(x=collection_date, 
                                 y=count, fill=variant, group=variant)) +
  # facet_wrap(~LABORATORY) +
  geom_area(aes(fill=variant), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("variant", values=c("grey75","red","blue","green3"), 
                    labels=c("wild type","501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-04-30")), expand=c(0,0)) +
  ylab("Share among newly diagnosed infections") +
  xlab("Collection date") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")
baseline_surveillance

# saveRDS(baseline_surveillance, file = paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3.png"), width=7, height=5)
# ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3.pdf"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
be_seq_mfit0 = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=be_basseqdata_long, 
                              subset=be_basseqdata_long$variant!="501Y.V1+V2+V3", maxit=1000)
be_seq_mfitallVOC = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=be_basseqdata_long, 
                                   subset=be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="501Y.V1+V2+V3", maxit=1000) 
BIC(be_seq_mfit0, be_seq_mfitallVOC)
summary(be_seq_mfit0)

# growth rate advantage compared to wild type
delta_r_501V1_501YV2_501YV3 = data.frame(confint(emtrends(be_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_501V1_501YV2_501YV3) = delta_r_501V1_501YV2_501YV3[,"contrast"]
delta_r_501V1_501YV2_501YV3 = delta_r_501V1_501YV2_501YV3[,-1]
delta_r_501V1_501YV2_501YV3
# estimate   asymp.LCL  asymp.UCL
# 501Y.V1 - wild type  0.037730104  0.03239347 0.04306674
# 501Y.V2 - wild type -0.008211584 -0.01880695 0.00238378
# 501Y.V3 - wild type  0.058270608  0.04725384 0.06928738

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
                   var="collection_date_num",  mode="latent",
                   at=list(collection_date_num=today_num), 
          df=NA)$contrasts
# contrast            estimate      SE df z.ratio p.value
# 501Y.V1 - wild type  0.03773 0.00272 NA 13.857  <.0001 
# 501Y.V2 - wild type -0.00821 0.00541 NA -1.519  0.4260 
# 501Y.V2 - 501Y.V1   -0.04594 0.00503 NA -9.125  <.0001 
# 501Y.V3 - wild type  0.05827 0.00562 NA 10.367  <.0001 
# 501Y.V3 - 501Y.V1    0.02054 0.00510 NA  4.030  0.0003 
# 501Y.V3 - 501Y.V2    0.06648 0.00702 NA  9.473  <.0001 
# 
# Degrees-of-freedom method: user-specified 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_501V1_501YV2_501YV3*4.7) 
# estimate asymp.LCL asymp.UCL
# 501Y.V1 - wild type 1.1940268 1.1644505  1.224354
# 501Y.V2 - wild type 0.9621408 0.9154014  1.011267
# 501Y.V3 - wild type 1.3150463 1.2486875  1.384931

# for all 3 variants together
# growth rate advantage compared to wild type
delta_r_allVOCs = data.frame(confint(emtrends(be_seq_mfitallVOC, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_allVOCs) = delta_r_allVOCs[,"contrast"]
delta_r_allVOCs = delta_r_allVOCs[,-1]
delta_r_allVOCs
#                               estimate  asymp.LCL  asymp.UCL
# (501Y.V1+V2+V3) - wild type 0.06314745 0.05977179 0.06652311

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_allVOCs*4.7) 
#                             estimate asymp.LCL asymp.UCL
# (501Y.V1+V2+V3) - wild type   1.345537  1.324357  1.367055




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

extrapolate = 30*6
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
           xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("variant", values=c("grey75","red","blue","green3")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-05-31")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share among newly diagnosed infections") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")
muller_be_seq_mfit0

# saveRDS(muller_be_seq_mfit0, file = paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit.png"), width=7, height=5)
# ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit.pdf"), width=7, height=5)


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
  coord_cartesian(xlim=c(min(be_basseqdata_long2$collection_date), as.Date("2021-05-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("red","blue","green3","black")) +
  scale_colour_manual("variant", values=c("red","blue","green3","black")) +
  geom_point(data=be_basseqdata_long2,
             aes(x=collection_date, y=100*prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("number\nsequenced", trans="sqrt",
                        range=c(1, 6), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_multinom_501YV1_501YV2_501YV3_response

# saveRDS(plot_multinom_501YV1_501YV2_501YV3_response, file = paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2_501YV3_multinomial fit_model preds_response.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds_response.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pdf"), width=8, height=6)


# on logit scale:

be_seq_mfit0_preds3 = be_seq_mfit0_preds2
ymin = 0.001
ymax = 0.990001
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
  scale_fill_manual("variant", values=c("red","blue","green3","black")) +
  scale_colour_manual("variant", values=c("red","blue","green3","black")) +
  geom_point(data=be_basseqdata_long2,
             aes(x=collection_date, y=prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(1, 6), limits=c(10,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom_501YV1_501YV2_501YV3


# saveRDS(plot_multinom_501YV1_501YV2_501YV3, file = paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.pdf"), width=8, height=6)


# estimated proportion of 501Y.V1 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18747 0.8624182 0.0145765 NA 0.8338487 0.8909876      2021-04-30

# estimated proportion of 501Y.V1 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18754 0.8565886 0.01968549 NA 0.8180058 0.8951715      2021-05-07

# estimated proportion of 501Y.V2 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18747 0.009103845 0.001811677 NA 0.005553024 0.01265467      2021-04-30

# estimated proportion of 501Y.V2 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18754 0.006555619 0.001531077 NA 0.003554763 0.009556476      2021-05-07

# estimated proportion of 501Y.V3 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18747 0.1011785 0.01453796 NA 0.07268463 0.1296724      2021-04-30

# estimated proportion of 501Y.V3 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)                18754 0.1160345 0.01987099 NA 0.07708807 0.1549809      2021-05-07



# estimated proportion of one of the three VOCs among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18747 0.9705554 0.003289346 NA 0.9641084 0.9770024      2021-04-30

# estimated proportion of one of the three VOCs among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18754 0.9769296 0.00300129 NA 0.9710472 0.9828121      2021-05-07




# the time at which new lab diagnoses would be by more than 50%, 75% 90% by one of the three VOCs :
be_seq_mfit0_preds2_subs = be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date>as.Date("2021-01-01")&
                                                 be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
# >50% by 10th of February [10th Feb - 11th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]

# >75% by 2nd of March [1st Mar - 3d March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]

# >90% by 26th of March [25th March - 28th of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]


# the time at which new infections would be by more than 50%, 75% 90% by one of the three VOCs
# (counting 7 days between infection & diagnosis) :
# >50% by 3d of Feb [3d Feb - 4th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]-7

# >75% by 23d of February [22nd Feb - 24th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]-7

# >90% by 19th of March [18th March - 21st of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]-7



# BASELINE PLUS ACTIVE SURVEILLANCE DATA ####

be_basplusactseqdata_long = gather(be_seqdata[,c("collection_date",
                                                 "basplusactivesurv_n_wild_type",
                                                 "basplusactivesurv_n_501Y.V1",
                                                 "basplusactivesurv_n_501Y.V2",
                                                 "basplusactivesurv_n_501Y.V3",
                                                 "basplusactivesurv_n_501Y.V1plusV2plusV3",
                                                 "basplusactivesurv_total_sequenced")], 
                                   variant, count, c("basplusactivesurv_n_wild_type",
                                                     "basplusactivesurv_n_501Y.V1",
                                                     "basplusactivesurv_n_501Y.V2",
                                                     "basplusactivesurv_n_501Y.V3",
                                                     "basplusactivesurv_n_501Y.V1plusV2plusV3"), factor_key=TRUE)
be_basplusactseqdata_long$variant = factor(be_basplusactseqdata_long$variant, 
                                           levels=c("basplusactivesurv_n_wild_type","basplusactivesurv_n_501Y.V1","basplusactivesurv_n_501Y.V2","basplusactivesurv_n_501Y.V3","basplusactivesurv_n_501Y.V1plusV2plusV3"), 
                                           labels=c("wild type", "501Y.V1", "501Y.V2", "501Y.V3", "501Y.V1+V2+V3"))
be_basplusactseqdata_long$collection_date_num = as.numeric(be_basplusactseqdata_long$collection_date)
be_basplusactseqdata_long$prop = be_basplusactseqdata_long$count / be_basplusactseqdata_long$basplusactivesurv_total_sequenced

basplusact_sequencing = ggplot(data=be_basplusactseqdata_long[be_basplusactseqdata_long$variant!="501Y.V1+V2+V3",], 
                               aes(x=collection_date, 
                                   y=count, fill=variant, group=variant)) +
  # facet_wrap(~LABORATORY) +
  geom_area(aes(fill=variant), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("variant", values=c("darkgrey","red","blue","green3"), 
                    labels=c("wild type","501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-03-01")), expand=c(0,0)) +
  ylab("Share among newly diagnosed infections") +
  xlab("Collection date") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline plus active surveillance)") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right")
basplusact_sequencing

# saveRDS(basplusact_sequencing, file = paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3.png"), width=7, height=5)
# ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3.pdf"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
be_basplusact_seq_mfit0 = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=be_basplusactseqdata_long, 
                                         subset=be_basplusactseqdata_long$variant!="501Y.V1+V2+V3"&
                                           be_basplusactseqdata_long$collection_date>=as.Date("2021-02-01"), 
                                         maxit=1000)
be_basplusact_seq_mfitallVOC = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=be_basplusactseqdata_long, 
                                              subset=be_basplusactseqdata_long$variant=="wild type"|be_basplusactseqdata_long$variant=="501Y.V1+V2+V3"&
                                                be_basplusactseqdata_long$collection_date>=as.Date("2021-02-01"), 
                                              maxit=1000) 
summary(be_basplusact_seq_mfit0)
BIC(be_basplusact_seq_mfit0)

# growth rate advantage compared to wild type
delta_r_501V1_501YV2_501YV3_basplusact = data.frame(confint(emtrends(be_basplusact_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                                     var="collection_date_num",  mode="latent",
                                                                     at=list(collection_date_num=today_num)), 
                                                            adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_501V1_501YV2_501YV3_basplusact) = delta_r_501V1_501YV2_501YV3_basplusact[,"contrast"]
delta_r_501V1_501YV2_501YV3_basplusact = delta_r_501V1_501YV2_501YV3_basplusact[,-1]
delta_r_501V1_501YV2_501YV3_basplusact

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(be_basplusact_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA)$contrasts

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_501V1_501YV2_501YV3_basplusact*4.7) 

# for all 3 variants together
# growth rate advantage compared to wild type
delta_r_allVOCs_basplusact = data.frame(confint(emtrends(be_basplusact_seq_mfitallVOC, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_allVOCs_basplusact) = delta_r_allVOCs_basplusact[,"contrast"]
delta_r_allVOCs_basplusact = delta_r_allVOCs_basplusact[,-1]
delta_r_allVOCs_basplusact

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_allVOCs_basplusact*4.7) 

# # PS: mblogit fit would also be possible & would take into account overdispersion

extrapolate = 90
date.from = as.numeric(as.Date("2020-11-01")) # min(be_basplusactseqdata_long$collection_date_num)
date.to = max(be_basplusactseqdata_long$collection_date_num)+extrapolate

be_basplusact_seq_mfit0_preds = data.frame(emmeans(be_basplusact_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_basplusact_seq_mfit0_preds$collection_date = as.Date(be_basplusact_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
be_basplusact_seq_mfit0_preds$variant = factor(be_basplusact_seq_mfit0_preds$variant, levels=c("wild type","501Y.V1","501Y.V2","501Y.V3"),
                                               labels=c("wild type","501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)"))

be_basplusact_seq_mfitallVOCs_preds = data.frame(emmeans(be_basplusact_seq_mfitallVOC, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_basplusact_seq_mfitallVOCs_preds$collection_date = as.Date(be_basplusact_seq_mfitallVOCs_preds$collection_date_num, origin="1970-01-01")
be_basplusact_seq_mfitallVOCs_preds$variant = factor(be_basplusact_seq_mfitallVOCs_preds$variant, levels=c("wild type","501Y.V1+V2+V3"),
                                                     labels=c("wild type","501Y.V1+V2+V3"))

be_basplusactseqdata_long2 = be_basplusactseqdata_long[be_basplusactseqdata_long$variant!="wild type",]
be_basplusactseqdata_long2$variant = droplevels(be_basplusactseqdata_long2$variant) 
be_basplusactseqdata_long2$variant = factor(be_basplusactseqdata_long2$variant, levels=c("501Y.V1","501Y.V2","501Y.V3","501Y.V1+V2+V3"),
                                            labels=c("501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)","501Y.V1+V2+V3"))

muller_be_basplusact_seq_mfit0 = ggplot(data=be_basplusact_seq_mfit0_preds, 
                                        aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(be_basplusactseqdata_long$collection_date)+1, 
           xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("variant", values=c("darkgrey","red","blue","green3")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-05-31")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share among newly diagnosed infections") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline plus active surveillance)")
muller_be_basplusact_seq_mfit0

# saveRDS(muller_be_basplusact_seq_mfit0, file = paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit.png"), width=7, height=5)
# ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit.pdf"), width=7, height=5)




# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS
be_basplusact_seq_mfit0_preds2 = rbind(be_basplusact_seq_mfit0_preds[be_basplusact_seq_mfit0_preds$variant!="wild type",],
                                       be_basplusact_seq_mfitallVOCs_preds[be_basplusact_seq_mfitallVOCs_preds$variant!="wild type",])
be_basplusact_seq_mfit0_preds2$variant = droplevels(be_basplusact_seq_mfit0_preds2$variant)


# on response scale:
plot_multinom_501YV1_501YV2_501YV3_basplusact_response = qplot(data=be_basplusact_seq_mfit0_preds2, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline plus active surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(be_basplusactseqdata_long2$collection_date), as.Date("2021-05-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("blue","red","green3","black")) +
  scale_colour_manual("variant", values=c("blue","red","green3","black")) +
  geom_point(data=be_basplusactseqdata_long2[be_basplusactseqdata_long2$collection_date>=as.Date("2021-02-01"),],
             aes(x=collection_date, y=100*prop, size=basplusactivesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("number\nsequenced", trans="sqrt",
                        range=c(1, 6), limits=c(0,10^5), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_multinom_501YV1_501YV2_501YV3_basplusact_response

# saveRDS(plot_multinom_501YV1_501YV2_501YV3_basplusact_response, file = paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2_501YV3_multinomial fit_model preds_response.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pdf"), width=8, height=6)


# on logit scale:

be_basplusact_seq_mfit0_preds3 = be_basplusact_seq_mfit0_preds2
ymin = 0.001
ymax = 0.990001
be_basplusact_seq_mfit0_preds3$asymp.LCL[be_basplusact_seq_mfit0_preds3$asymp.LCL<ymin] = ymin
be_basplusact_seq_mfit0_preds3$asymp.UCL[be_basplusact_seq_mfit0_preds3$asymp.UCL<ymin] = ymin
be_basplusact_seq_mfit0_preds3$asymp.UCL[be_basplusact_seq_mfit0_preds3$asymp.UCL>ymax] = ymax
be_basplusact_seq_mfit0_preds3$prob[be_basplusact_seq_mfit0_preds3$prob<ymin] = ymin

plot_multinom_501YV1_501YV2_501YV3_basplusact = qplot(data=be_basplusact_seq_mfit0_preds3, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline plus active surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("blue","red","green3","black")) +
  scale_colour_manual("variant", values=c("blue","red","green3","black")) +
  geom_point(data=be_basplusactseqdata_long2[be_basplusactseqdata_long2$collection_date>=as.Date("2021-02-01"),],
             aes(x=collection_date, y=prop, size=basplusactivesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(1, 4), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom_501YV1_501YV2_501YV3_basplusact


# saveRDS(plot_multinom_501YV1_501YV2_501YV3, file = paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\basplusact_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.pdf"), width=8, height=6)




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
# 62           18718 0.9999585 5.789292e-05 Inf 0.9993608 0.9999973      2021-04-01

# prop of S dropout samples among new infections that are now estimated to be B.1.1.7 / 501Y.V1 (using 7 days for time from infection to diagnosis)
fitseq_preds[fitseq_preds$collection_date==(today+7),]
#    collection_date_num     prob          SE  df asymp.LCL asymp.UCL collection_date
# 69           18725 0.9999809 2.884737e-05 Inf 0.9996314  0.999999      2021-04-08

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
ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_link scale.png"), width=8, height=6)
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

# saveRDS(plot_fitseq, file = paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_501YV1_binomial GLMM_response scale.pdf"), width=8, height=6)




# 3. ANALYSIS OF Ct VALUES OF S-DROPOUT & NON-S-DROPOUT SAMPLES IN BELGIUM ####

# Read in Ct data of all valid PCRs
file_dec1 = paste0(".//data//", dat, "//1-10 december 2020.xlsx")
file_dec2 = paste0(".//data//", dat, "//11-20 december 2020.xlsx")
file_dec3 = paste0(".//data//", dat, "//21-31 december 2020.xlsx")
file_jan1 = paste0(".//data//", dat, "//1-9 january 2021.xlsx")
file_jan2 = paste0(".//data//", dat, "//10-19 january 2021.xlsx")
file_jan3 = paste0(".//data//", dat, "//20-31 january 2021.xlsx")
file_feb1 = paste0(".//data//", dat, "//1-14 february 2021.xlsx")
file_feb2 = paste0(".//data//", dat, "//15-28 february 2021.xlsx")
file_mar1 = paste0(".//data//", dat, "//1-15 march 2021.xlsx")
sheets = excel_sheets(file_jan1)
ctdata_dec1 = map_df(sheets, ~ read_excel(file_dec1, sheet = .x, skip = 1,
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"),
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric")))
range(as.Date(as.numeric(ctdata_dec1$Analysis_date), origin="1899-12-30")) # "2020-12-01" "2020-12-10"
ctdata_dec2 = map_df(sheets, ~ read_excel(file_dec2, sheet = .x, skip = 1,
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"),
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric")))
range(as.Date(as.numeric(ctdata_dec2$Analysis_date), origin="1899-12-30")) # "2020-12-11" "2020-12-20"
ctdata_dec3 = map_df(sheets, ~ read_excel(file_dec3, sheet = .x, skip = 1,
                                          col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"),
                                          col_types=c("text","text","text","numeric","numeric","numeric","numeric")))
range(as.Date(as.numeric(ctdata_dec3$Analysis_date), origin="1899-12-30")) # "2020-12-21" "2020-12-31"
ctdata_jan1 = map_df(sheets, ~ read_excel(file_jan1, sheet = .x, skip = 1, 
                                       col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                       col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_jan1$Analysis_date), origin="1899-12-30")) # "2021-01-01" "2021-01-09"
ctdata_jan2 = map_df(sheets, ~ read_excel(file_jan2, sheet = .x, skip = 1, 
                                          col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                          col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_jan2$Analysis_date), origin="1899-12-30")) # "2021-01-10" "2021-01-19"
ctdata_jan3 = map_df(sheets, ~ read_excel(file_jan3, sheet = .x, skip = 1, 
                                          col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                          col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_jan3$Analysis_date), origin="1899-12-30")) # "2021-01-20" "2021-01-31"
ctdata_feb1 = map_df(sheets, ~ read_excel(file_feb1, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_feb1$Analysis_date), origin="1899-12-30")) # "2021-02-01" "2021-01-14"
ctdata_feb2 = map_df(sheets, ~ read_excel(file_feb2, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_feb2$Analysis_date), origin="1899-12-30")) # "2021-02-15" "2021-02-28"
ctdata_mar1 = map_df(sheets, ~ read_excel(file_mar1, sheet = .x, skip = 1, 
                                          col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                          col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
range(as.Date(as.numeric(ctdata_mar1$Analysis_date), origin="1899-12-30")) # "2021-03-01" "2021-03-15"
ctdata = bind_rows(ctdata_dec1, ctdata_dec2, ctdata_dec3, ctdata_jan1, ctdata_jan2, ctdata_jan3, ctdata_feb1, ctdata_feb2, ctdata_mar1)
range(as.Date(as.numeric(ctdata$Analysis_date), origin="1899-12-30")) # "2020-12-01" "2021-03-15"
ctdata$Laboratory[ctdata$Laboratory=="ULG - FF 3.x"] = "ULG"
unique(ctdata$Laboratory) 
# "ULG"              "Namur"            "UMons - Jolimont" "UZ leuven"        "UZA"              "UZ Gent"          "ULB"             
# "Saint LUC - UCL" 
unique(ctdata$Outcome) 
# unique(ctdata_dec1$Outcome) 
# unique(ctdata_dec2$Outcome) 
unique(ctdata_dec1$Outcome) 
unique(ctdata_dec2$Outcome) 
unique(ctdata_dec3$Outcome) 
unique(ctdata_jan1$Outcome) 
unique(ctdata_jan2$Outcome) 
unique(ctdata_jan3$Outcome) 
unique(ctdata_feb1$Outcome) 
unique(ctdata_feb2$Outcome)
unique(ctdata_mar1$Outcome)
ctdata$Outcome[ctdata$Outcome=="Detected"] = "Positive"
ctdata$Outcome[ctdata$Outcome=="Not detected"] = "Negative"
unique(ctdata$Outcome) # "Positive" "Negative"
ctdata$Analysis_date = as.Date(as.numeric(ctdata$Analysis_date), origin="1899-12-30")
sum(is.na(ctdata$Analysis_date)) # 0
range(ctdata$Analysis_date) # "2020-12-01" - "2021-03-15"
ctdata$collection_date = ctdata$Analysis_date-1 # collection date = analysis date-1 
sum(is.na(ctdata$collection_date)) # 0
ctdata$collection_date_num = as.numeric(ctdata$collection_date)
range(ctdata$collection_date) # "2020-11-30" "2021-03-14"
ctdata$group = interaction(ctdata$Outcome, ctdata$S_dropout)
ctdata$group = droplevels(ctdata$group)
unique(ctdata$group) # Positive.0 Positive.1 Negative.0
ctdata$Outcome = factor(ctdata$Outcome)
unique(ctdata$Outcome) # Positive Negative
ctdata$Laboratory = factor(ctdata$Laboratory)
ctdata$S_dropout = factor(ctdata$S_dropout)
head(ctdata)
str(ctdata)
nrow(ctdata) # 811436

ctdata_onlypos = ctdata[ctdata$Outcome=="Positive",] # subset with only the positive samples
ctdata_onlypos = bind_rows(ctdata_onlypos[ctdata_onlypos$S_dropout=="0",], ctdata_onlypos[ctdata_onlypos$S_dropout=="1",])

# ANALYSIS OF Ct VALUES OF S DROPOUT & NON-S DROPOUT SAMPLES

# plot & analysis of Ct values of all labs for dates from 13th of Jan onward when >80% of all S dropouts were B.1.1.7 / 501Y.V1

(fitseq_preds[fitseq_preds$prob>0.9,"collection_date"][1]) # "2021-01-20", from 13th of Jan >90% of all S dropouts are B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections
# this is the same criterion that was used for the SGTF analysis in the UK (N. Davies, pers. comm.)
subs = (ctdata_onlypos$collection_date > (fitseq_preds[fitseq_preds$prob>0.9,"collection_date"][1])) & 
       (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30) 
ctdata_onlypos_subs = ctdata_onlypos[subs,]
ctdata_onlypos_subs = ctdata_onlypos_subs[!(is.na(ctdata_onlypos_subs$S_dropout)|
                                              is.na(ctdata_onlypos_subs$N_cq)|
                                              is.na(ctdata_onlypos_subs$ORF1_cq)|
                                              (ctdata_onlypos_subs$ORF1_cq==0)),]

cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.897

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                        Laboratory correlation_Ct_N_ORF1ab
# Namur                       Namur               0.9520116
# Saint LUC - UCL   Saint LUC - UCL               0.9814627
# ULB                           ULB               0.9823499
# ULG                           ULG               0.9589549
# UMons - Jolimont UMons - Jolimont               0.8903198
# UZ Gent                   UZ Gent               0.2733935
# UZ leuven               UZ leuven               0.9769981
# UZA                           UZA               0.9596583

ctcorplot_all_labs = qplot(data=ctdata_onlypos_subs, x=ORF1_cq, y=N_cq, group=Laboratory, fill=S_dropout, colour=S_dropout, size=I(3), shape=I(16)) +
  facet_wrap(~Laboratory) + 
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Ct value ORF1ab gene") + ylab("Ct value N gene") +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctcorplot_all_labs
# PS there is problem in exported Ct values for UZ Gent, this is being fixed
# saveRDS(ctcorplot_all_labs, file = paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.pptx"), width=7, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.png"), width=7, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.pdf"), width=7, height=6)


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
# PS UZA has suspect Ct values before the 21st of Jan (software export problem), UMons after 10th of March 
# U Gent also still has confirmed software problem with exported Ct values, this is being fixed
ctdata_onlypos_subs[ctdata_onlypos_subs$Laboratory=="UZA"&ctdata_onlypos_subs$N_cq<=15,"collection_date"][1,] # "2021-01-21"
tail(ctdata_onlypos_subs[ctdata_onlypos_subs$Laboratory=="UMons - Jolimont"&ctdata_onlypos_subs$N_cq<=15,"collection_date"],1) # "2021-03-10"
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

# saveRDS(ctplots_rawdata_all_labs, file = paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.pptx"), width=7, height=9)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.png"), width=7, height=9)
# ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.pdf"), width=7, height=9)



# plot & analysis of Ct values for labs with most comparable Ct value distributions (all except "UZ Gent", still a software problem there)
# ("UZ leuven","Saint LUC - UCL","ULB","Namur","ULG","Umons - Jolimont") 
# for dates from 20th of Jan onward when >90% of all S dropouts were B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections

labs_to_remove = c("UZ Gent")
sel_labs = setdiff(unique(ctdata_onlypos$Laboratory), labs_to_remove) 
# sel_labs = c("UZ leuven","Saint LUC - UCL","ULB","Namur","ULG","Umons - Jolimont","UZA"))  
# we use data from these 4 labs as the data distribution was comparable for these
# they also had large sample size & were not heavily involved in active surveillance
# sel_labs = unique(ctdata_onlypos$Laboratory) # to select data from all the labs, but distribution not comparable for all
# we use the subset of timepoints (from 13th Jan 2021 onwards) where >80% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
(fitseq_preds[fitseq_preds$prob>0.9,"collection_date"][1]) # "2021-01-20", from 20th of Jan >90% of all S dropouts are B.1.1.7 / 501Y.V1
# we also just use the positive samples with relatively strong signal, (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30)
# to not include pos samples with very low viral titers (indicative of old infections etc)
# this is the same criterion that was used for the SGTF analysis in the UK (N. Davies, pers. comm.)
subs = (ctdata_onlypos$collection_date > (fitseq_preds[fitseq_preds$prob>0.9,"collection_date"][1])) &
  (ctdata_onlypos$Laboratory %in% sel_labs) & (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30) 
ctdata_onlypos_subs = ctdata_onlypos[subs,]
ctdata_onlypos_subs = ctdata_onlypos_subs[!is.na(ctdata_onlypos_subs$collection_date)|
                                            !(is.na(ctdata_onlypos_subs$S_dropout)|
                                              is.na(ctdata_onlypos_subs$N_cq)|
                                              is.na(ctdata_onlypos_subs$ORF1_cq)|
                                              (ctdata_onlypos_subs$ORF1_cq==0)),]
delrows = (ctdata_onlypos_subs$Laboratory=="UMons - Jolimont") & (ctdata_onlypos_subs$collection_date>as.Date("2021-03-10")) # temporary fix for UMons data
ctdata_onlypos_subs = ctdata_onlypos_subs[!delrows,]
nrow(ctdata_onlypos_subs) # 24863
ctdata_onlypos_subs$Laboratory = droplevels(ctdata_onlypos_subs$Laboratory)

  
# make joint dataset for integrated analysis of both genes to estimate average effect across both sets of genes
ctdata_onlypos_subs_bothgenes = rbind(data.frame(ctdata_onlypos_subs, Gene="N gene", Ct=ctdata_onlypos_subs$N_cq), 
                                      data.frame(ctdata_onlypos_subs, Gene="ORF1ab gene", Ct=ctdata_onlypos_subs$ORF1_cq))
# we define a high viral load as one where the Ct value was 1.25x lower than in the non-S dropout sample group
# which was a Ct value < 15.02 for the N gene and < 16.08 for the ORF1ab gene
thresh_N = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","N_cq"]))/1.25
thresh_N # 15.02
thresh_ORF1 = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","ORF1_cq"]))/1.25
thresh_ORF1 # 16.08
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="N gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="N gene"]<thresh_N
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"]<thresh_ORF1


# check correlation between Ct values for N & ORF1ab gene
cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.97, t=593.77, p<2E-16

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                        Laboratory correlation_Ct_N_ORF1ab
# Namur                       Namur               0.9520116
# Saint LUC - UCL   Saint LUC - UCL               0.9814627
# ULB                           ULB               0.9823499
# ULG                           ULG               0.9589549
# UMons - Jolimont UMons - Jolimont               0.9640394
# UZ leuven               UZ leuven               0.9769981
# UZA                           UZA               0.9596583


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
# saveRDS(ctcorplot_sellabs, file = paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_selected labs.pptx"), width=7, height=4.5)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_selected labs.png"), width=7, height=4.5)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_selected labs.pdf"), width=7, height=4.5)


ctplot_rawdataN = qplot(data=ctdata_onlypos_subs, x=collection_date, y=N_cq, group=S_dropout, 
                        colour=S_dropout, fill=S_dropout, geom="point", shape=I(16), size=I(1)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "N gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataN # there is no obvious temporal patterns

ctplot_rawdataORF1 = qplot(data=ctdata_onlypos_subs, x=collection_date, y=ORF1_cq, group=S_dropout, 
                           colour=S_dropout, fill=S_dropout, geom="point", shape=I(16), size=I(1)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "ORF1ab gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataORF1 # there is no obvious temporal patterns

ctplots_rawdata = ggarrange(ctplot_rawdataN+xlab("")+theme(axis.text.x = element_blank()), 
                            ctplot_rawdataORF1,
                            ncol=1, common.legend=TRUE, legend="right")
ctplots_rawdata # there is no obvious temporal patterns

# saveRDS(ctplots_rawdata, file = paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_selected labs.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_selected labs.pptx"), width=7, height=9)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_selected labs.png"), width=7, height=9)
# ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_selected labs.pdf"), width=7, height=9)




# quantile/median regression to compare median Ct values of both genes across S dropout & non-S dropout samples in the different labs
qr_bothgenes0 = rq(Ct ~ Gene + S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes1 = rq(Ct ~ Gene * S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes2 = rq(Ct ~ Gene + S_dropout * Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes3 = rq(Ct ~ Gene * Laboratory + S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes4 = rq(Ct ~ (Gene + Laboratory + S_dropout)^2, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes5 = rq(Ct ~ Gene * Laboratory * S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
AIC(qr_bothgenes0, k=-1) # 324433.7
AIC(qr_bothgenes1, k=-1) # 324421.8
AIC(qr_bothgenes2, k=-1) # 324162.2 # fits data best (PS: here AIC with k<0 returns BIC)
AIC(qr_bothgenes3, k=-1) # 324460.7
AIC(qr_bothgenes4, k=-1) # 324165.3 # 2nd best (PS: here AIC with k<0 returns BIC)
AIC(qr_bothgenes5, k=-1) # 324222.6

summary(qr_bothgenes4) # I will continue with qr_bothgenes4 since it is more general
# tau: [1] 0.5
# 
# Coefficients:
#   Value     Std. Error t value   Pr(>|t|) 
# (Intercept)             18.36546   0.04195  437.81624   0.00000
# Gene1                   -0.79366   0.04038  -19.65526   0.00000
# Laboratory1             -0.31046   0.09113   -3.40670   0.00066
# Laboratory2             -0.08966   0.12354   -0.72577   0.46798
# Laboratory3             -0.89461   0.10564   -8.46865   0.00000
# Laboratory4              1.00597   0.12218    8.23371   0.00000
# Laboratory5             -0.01596   0.08307   -0.19215   0.84762
# Laboratory6              0.90304   0.09064    9.96285   0.00000
# S_dropout1               1.04311   0.04171   25.00875   0.00000
# Gene1:Laboratory1        0.06176   0.08834    0.69913   0.48447
# Gene1:Laboratory2        0.04126   0.10728    0.38459   0.70054
# Gene1:Laboratory3       -0.14644   0.10057   -1.45614   0.14536
# Gene1:Laboratory4        0.61564   0.11659    5.28031   0.00000
# Gene1:Laboratory5       -0.16804   0.08291   -2.02665   0.04270
# Gene1:Laboratory6       -0.37334   0.08979   -4.15809   0.00003
# Gene1:S_dropout1         0.17895   0.04108    4.35574   0.00001
# Laboratory1:S_dropout1   0.07104   0.09117    0.77916   0.43589
# Laboratory2:S_dropout1   0.50464   0.12360    4.08270   0.00004
# Laboratory3:S_dropout1  -0.33311   0.10550   -3.15738   0.00159
# Laboratory4:S_dropout1  -0.48933   0.12258   -3.99204   0.00007
# Laboratory5:S_dropout1  -0.61896   0.08337   -7.42432   0.00000
# Laboratory6:S_dropout1   0.99234   0.08868   11.19036   0.00000

qr_emmeans_bylab99 = data.frame(emmeans(qr_bothgenes4, ~ Laboratory + Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans_bylab99
qr_emmeans99 = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans99
qr_emmeans = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.95)) # median Ct values + 95% CLs
qr_emmeans
# Gene S_dropout   emmean         SE    df lower.CL upper.CL
# 1      N gene         0 18.79386 0.09847595 49704 18.60085 18.98688
# 2 ORF1ab gene         0 20.02328 0.08823135 49704 19.85035 20.19622
# 3      N gene         1 16.34974 0.06568341 49704 16.22100 16.47848
# 4 ORF1ab gene         1 18.29496 0.07398727 49704 18.14995 18.43998

# mean difference in median Ct value of 2.09, which is highly significant across both genes: p<0.0001
contrast(emmeans(qr_bothgenes4, ~ S_dropout, level=0.95), method="pairwise") 
# contrast estimate    SE    df t.ratio p.value
# 0 - 1           2.09 0.0834 49704 25.009  <.0001
# mean difference in median Ct value of 2.39 for N gene and 1.60 for ORF1ab gene
confint(contrast(emmeans(qr_bothgenes4, ~ S_dropout|Gene, level=0.95), method="pairwise"))
# Gene = N gene:
#   contrast estimate    SE    df lower.CL upper.CL
# 0 - 1        2.44 0.120 49704     2.21     2.68
# 
# Gene = ORF1ab gene:
#   contrast estimate    SE    df lower.CL upper.CL
# 0 - 1        1.73 0.114 49704     1.50     1.95
# 
# Results are averaged over the levels of: Laboratory 
# Confidence level used: 0.95 

# violin plots by gene & lab & S dropout with expected marginal means+99% CLs of best fitting median regression model
ctviolinplots_bylab = ggplot(data=ctdata_onlypos_subs_bothgenes, aes(x=factor(S_dropout), y=Ct, fill=factor(S_dropout))) +
  geom_violin(alpha=1, colour=NA, trim=TRUE, draw_quantiles=TRUE, adjust=2, scale="width") +
  geom_crossbar(data=qr_emmeans_bylab99, aes(x=factor(S_dropout), y=emmean, ymin=lower.CL, ymax=upper.CL, group=Gene, lwd=I(0.1))) +
  # stat_summary(fun.data=data_summary,  
  #             geom="pointrange", aes(color=factor(S_dropout))) +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  # geom_point(aes(colour=factor(S_dropout))) +
  facet_wrap(~ Gene + Laboratory,ncol=length(levels(ctdata_onlypos_subs_bothgenes$Laboratory))) +
  scale_colour_manual("", values=alpha(c("steelblue","lightcoral"), 1), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("steelblue","lightcoral"), 1), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("") + ylab("Ct value") + 
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("0","1"), labels=c("S pos","S dropout")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctviolinplots_bylab

# saveRDS(ctviolinplots_bylab, file = paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.pptx"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.pdf"), width=6, height=6)

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

# saveRDS(ctviolinplots, file = paste0(".\\plots\\",dat,"\\Ct values_violin plots.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.pptx"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.png"), width=6, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.pdf"), width=6, height=6)




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
# fitct_highvirload_0A  4 63605.16
# fitct_highvirload_1A  5 63611.29
# fitct_highvirload_2A  6 63610.33
# fitct_highvirload_3A  6 63622.01
# fitct_highvirload_4A  7 63621.06
# fitct_highvirload_0B  5 63585.95
# fitct_highvirload_1B  6 63592.06
# fitct_highvirload_2B  9 63610.56
# fitct_highvirload_3B  7 63602.79
# fitct_highvirload_4B 10 63621.30


# fitct_highvirload_0B the best model
summary(fitct_highvirload_0B) # S dropout samples more frequently have high viral load based on N gene Ct values
# Random effects:
#   Groups     Name        Variance Std.Dev.
# Laboratory (Intercept) 0.04825  0.2197  
# Number of obs: 49726, groups:  Laboratory, 7
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.698582   0.084025  -8.314  < 2e-16 ***
#   Gene1             0.119089   0.009509  12.524  < 2e-16 ***
#   S_dropout1       -0.177787   0.009562 -18.593  < 2e-16 ***
#   Gene1:S_dropout1 -0.052105   0.009509  -5.480 4.26e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(allEffects(fitct_highvirload_0B))


# odds to encounter high viral load samples based on Ct values of both genes (high vir load = Ct values >1.25x lower than median in non-S dropout samples)
# 1.43x [1.37-1.48x] 95% CLs increased among S dropout samples
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout, type="response"), method="revpairwise", type="response"))
# contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.43 0.0273 Inf      1.37      1.48

# odds ratio = 1.58 [1.50-1.67] for N gene & 1.29 [1.22-1.36] for ORF1ab gene 
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout|Gene, type="response"), method="revpairwise", type="response"))
# Gene = N gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.58 0.0421 Inf      1.50      1.67
# 
# Gene = ORF1ab gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.29 0.0352 Inf      1.22      1.36
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log odds ratio scale


fitct_highvirload_emmeans = as.data.frame(emmeans(fitct_highvirload_0B, ~ S_dropout+Gene, type="response"))
fitct_highvirload_emmeans$S_dropout = factor(fitct_highvirload_emmeans$S_dropout)
fitct_highvirload_emmeans$Gene = factor(fitct_highvirload_emmeans$Gene)
fitct_highvirload_emmeans
# S_dropout        Gene      prob         SE  df asymp.LCL asymp.UCL
# 1         0      N gene 0.3080217 0.01826386 Inf 0.2734220 0.3449212
# 2         1      N gene 0.4134791 0.02071798 Inf 0.3735444 0.4545852
# 3         0 ORF1ab gene 0.2802236 0.01730975 Inf 0.2475807 0.3153663
# 4         1 ORF1ab gene 0.3335908 0.01902555 Inf 0.2973965 0.3718587

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

# saveRDS(ctplots_all_rawdata, file = paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.pdf"), width=8, height=6)




# 4. ESTIMATE GROWTH RATE AND TRANSMISSION ADVANTAGE OF B.1.1.7 / 501Y.V1 IN BELGIUM BASED ON S-GENE TARGET FAILURE DATA ####

# we use data from all labs here except U Gent (still Ugentec software problem with exported Ct values) & ULG (low sample size)
# and a Ct cutoff of 30 for the N & ORF1 gene to be able to focus on recent infections
sel_labs = setdiff(unique(ctdata$Laboratory), c("UZ Gent", "ULG")) 
# sel_labs = unique(ctdata$Laboratory)
# setdiff(unique(ctdata$Laboratory), c("UMons - Jolimont", "ULG", "UZ Gent", "UZA")) 
# setdiff(unique(ctdata$Laboratory), c("UMons - Jolimont", "UZ Gent","UZA","ULG")) 
sel_labs  

# ctcutoff = max(c(ctdata$N_cq[!is.na(ctdata$N_cq)], ctdata$ORF1_cq[!is.na(ctdata$ORF1_cq)])) + 1 # i.e. no Ct cutoff
ctcutoff = 30
  
subs = which((ctdata$Laboratory %in% sel_labs) & (((ctdata$Outcome=="Positive")&(ctdata$N_cq<ctcutoff)&(ctdata$ORF1_cq<ctcutoff))|
                                            (ctdata$Outcome=="Negative")))
ctdata_subs = ctdata
ctdata_subs = ctdata_subs[subs,]
nrow(ctdata) # 811436
nrow(ctdata_subs) # 647611

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
tail(data_ag_byday_wide)

write.csv(data_ag_byday_wide, file=".\\data\\be_latest\\be_B117_total.csv", row.names=FALSE)



# 4.1 ESTIMATE GROWTH RATE & TRANSMISSION ADVANTAGE OF 501Y.V1 USING BINOMIAL GLMM (LOGISTIC FIT) ####

# fit common-slope and separate-slopes binomial GLM with or without natural cubic spline terms in function of time
set_sum_contrasts()
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit1)
# fit1_22jan is fit as used in original report of jan 28 using data from jan 1 to jan 22
fit1_22jan = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)+LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0&data_ag_wide$collection_date>=as.Date("2021-01-01")&data_ag_wide$collection_date<=as.Date("2021-01-22"), control=glmersettings2)  # common slope model, with lab coded as fixed factor
fit1 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)+LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings2)  # common slope model, with lab coded as fixed factor
fit2 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)*LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings) # separate slopes model, with lab coded as fixed factor
fit3 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+ns(collection_date_num,df=2)+LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings2)  # with 2 df natural cubic spline term ifo date + laboratory
fit4 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+ns(collection_date_num,df=2)*LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings) # with 2 df natural cubic spline term ifo date * laboratory
BIC(fit1,fit2,fit3,fit4) 
#      df      BIC
# fit1  8 3077.410
# fit2 13 3032.537
# fit3  9 2895.935
# fit4 19 2801.960


# separate-slope spline model fit4 fits best, i.e. rate at which 501Y.V1 is displacing other strains differs across regions/labs & has nonlinearities

summary(fit4)
# Random effects:
#   Groups Name        Variance Std.Dev.
# obs    (Intercept) 0.1414   0.3761  
# Number of obs: 626, groups:  obs, 626
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  -8.95196    0.30379 -29.467  < 2e-16 ***
#   ns(collection_date_num, df = 2)1             15.02083    0.56773  26.458  < 2e-16 ***
#   ns(collection_date_num, df = 2)2              6.26252    0.13508  46.360  < 2e-16 ***
#   LABORATORY1                                  -0.81835    0.79734  -1.026   0.3047    
# LABORATORY2                                   0.79891    0.54856   1.456   0.1453    
# LABORATORY3                                   0.06999    0.60899   0.115   0.9085    
# LABORATORY4                                  -0.77304    0.85072  -0.909   0.3635    
# LABORATORY5                                  -0.78477    0.65790  -1.193   0.2329    
# ns(collection_date_num, df = 2)1:LABORATORY1  1.39387    1.48494   0.939   0.3479    
# ns(collection_date_num, df = 2)2:LABORATORY1  0.37805    0.35019   1.080   0.2803    
# ns(collection_date_num, df = 2)1:LABORATORY2 -1.06974    1.04028  -1.028   0.3038    
# ns(collection_date_num, df = 2)2:LABORATORY2 -1.34958    0.24246  -5.566 2.60e-08 ***
#   ns(collection_date_num, df = 2)1:LABORATORY3  0.67804    1.14839   0.590   0.5549    
# ns(collection_date_num, df = 2)2:LABORATORY3 -0.49155    0.26523  -1.853   0.0638 .  
# ns(collection_date_num, df = 2)1:LABORATORY4 -0.04171    1.56841  -0.027   0.9788    
# ns(collection_date_num, df = 2)2:LABORATORY4  2.24299    0.38283   5.859 4.66e-09 ***
#   ns(collection_date_num, df = 2)1:LABORATORY5  1.17749    1.23453   0.954   0.3402    
# ns(collection_date_num, df = 2)2:LABORATORY5 -0.30125    0.28974  -1.040   0.2985 

# growth rate advantage (differences in growth rate between 501Y.V1 and old strains):

# avg growth advantage of 501Y.V1 based on fit4 today
fit_emtrends = as.data.frame(emtrends(fit4, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$emtrends)
fit_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL  asymp.UCL
# 1               0.03298381 0.02776099 0.03820663

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_emtrends[,c(2,5,6)]*4.7) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.167686  1.139371  1.196704

# avg growth advantage of 501Y.V1 on the 1st of January 2021
fit_emtrends_1jan = as.data.frame(emtrends(fit4, revpairwise ~ 1, var="collection_date_num", 
                                      at=list(collection_date_num=as.numeric(as.Date("2021-01-01"))), mode="link", adjust="Tukey")$emtrends)
fit_emtrends_1jan[,c(2,5,6)]
# collection_date_num.trend asymp.LCL  asymp.UCL
# 1                 0.1352384 0.1244043 0.1460726

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_emtrends_1jan[,c(2,5,6)]*4.7) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1.888194  1.794453  1.986832

# results original first report of 28 jan (using data from jan 1 tot jan 22 and using fit1)
fit1_22jan_emtrends = as.data.frame(emtrends(fit1_22jan, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$emtrends)
fit1_22jan_emtrends[,c(2,5,6)]
# collection_date_num.trend  asymp.LCL asymp.UCL
# 1                 0.1124202 0.08875747  0.136083

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_22jan_emtrends[,c(2,5,6)]*4.7) 
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.696175  1.517646  1.895705



# tests for differences in date of introduction based on fit1
# UCL, ULB, UZA earlier than avg, Namur, Mons & UZ Leuven later than avg
emmeans(fit1,eff~LABORATORY)$contrasts 
# contrast                  estimate     SE  df z.ratio p.value
# Namur effect               -0.0815 0.0777 Inf -1.048  0.2946 
# (Saint LUC - UCL) effect    0.0911 0.0720 Inf  1.266  0.2466 
# ULB effect                  0.3505 0.0726 Inf  4.826  <.0001 
# (UMons - Jolimont) effect  -0.4522 0.0767 Inf -5.894  <.0001 
# UZ leuven effect           -0.2669 0.0743 Inf -3.595  0.0005 
# UZA effect                  0.3590 0.0713 Inf  5.035  <.0001 
# 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: fdr method for 8 tests


# results of growth rate advantage of separate-slopes spline model fit4 by lab/region evaluated today:  
# large heterogeneity in growth rate advantage of 501Y.V1 - large in Mons, Namur & Antwerp
fit4_emtrends = emtrends(fit4, revpairwise ~ LABORATORY, var="collection_date_num", 
                         at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$emtrends
fit4_emtrends
# LABORATORY       collection_date_num.trend      SE  df asymp.LCL asymp.UCL
# Namur                               0.0310 0.00713 Inf   0.01707    0.0450
# Saint LUC - UCL                     0.0084 0.00590 Inf  -0.00317    0.0200
# ULB                                 0.0154 0.00643 Inf   0.00282    0.0280
# UMons - Jolimont                    0.0885 0.00727 Inf   0.07426    0.1028
# UZ leuven                           0.0161 0.00651 Inf   0.00332    0.0288
# UZA                                 0.0384 0.00578 Inf   0.02711    0.0498

# tests for above- or below-average growth rate of B.1.1.7 evaluated today:
# above-avg rapid spread in Mons, below-avg rapid spread in Brussels (Saint Luc & ULB) & Leuven
fit4_contrasts = emtrends(fit4, eff ~ LABORATORY, var="collection_date_num", 
                          at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$contrasts
fit4_contrasts
# contrast                  estimate      SE  df z.ratio p.value
# Namur effect              -0.00194 0.00641 Inf -0.302  0.9998 
# (Saint LUC - UCL) effect  -0.02458 0.00551 Inf -4.464  <.0001 
# ULB effect                -0.01756 0.00589 Inf -2.983  0.0170 
# (UMons - Jolimont) effect  0.05553 0.00651 Inf  8.531  <.0001 
# UZ leuven effect          -0.01691 0.00594 Inf -2.845  0.0264 
# UZA effect                 0.00546 0.00542 Inf  1.007  0.8957 
# 
# P value adjustment: sidak method for 7 tests 


# results of growth rate advantage of separate-slopes spline model fit4 by lab/region evaluated on 1st of Jan:  
# large heterogeneity in growth rate advantage of 501Y.V1 - large in Mons, Namur & Antwerp
fit4_emtrends_1jan = emtrends(fit4, revpairwise ~ LABORATORY, var="collection_date_num", 
                         at=list(collection_date_num=as.numeric(as.Date("2021-01-01"))), mode="link", adjust="Tukey")$emtrends
fit4_emtrends_1jan
# LABORATORY       collection_date_num.trend      SE  df asymp.LCL asymp.UCL
# Namur                                0.148 0.01632 Inf    0.1165     0.180
# Saint LUC - UCL                      0.129 0.01054 Inf    0.1079     0.149
# ULB                                  0.144 0.01202 Inf    0.1203     0.167
# UMons - Jolimont                     0.128 0.01733 Inf    0.0936     0.161
# UZ leuven                            0.148 0.01315 Inf    0.1227     0.174
# UZA                                  0.115 0.00982 Inf    0.0954     0.134

# tests for above- or below-average growth rate of B.1.1.7 evaluated on 1st of Jan:
# no sign differences
fit4_contrasts_1jan = emtrends(fit4, eff ~ LABORATORY, var="collection_date_num", 
                          at=list(collection_date_num=as.numeric(as.Date("2021-01-01"))), mode="link", adjust="Tukey")$contrasts
fit4_contrasts_1jan
# contrast                  estimate      SE  df z.ratio p.value
# Namur effect               0.01321 0.01442 Inf  0.916  0.9310 
# (Saint LUC - UCL) effect  -0.00669 0.01019 Inf -0.656  0.9864 
# ULB effect                 0.00862 0.01124 Inf  0.767  0.9701 
# (UMons - Jolimont) effect -0.00773 0.01518 Inf -0.509  0.9965 
# UZ leuven effect           0.01318 0.01204 Inf  1.094  0.8533 
# UZA effect                -0.02059 0.00971 Inf -2.121  0.1870 
# 
# P value adjustment: sidak method for 7 tests 



# PLOT MODEL FIT

# best fitting separate-slope spline fit fit4, overall avg model prediction over all labs
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-05-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit4))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_preds = as.data.frame(emmeans(fit4, ~ collection_date_num, 
                                         # by="LABORATORY", 
                                         at=list(collection_date_num=seq(date.from,
                                                                         date.to)), 
                                         type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_preds$collection_date = as.Date(fit_preds$collection_date_num, origin="1970-01-01")

# original fit of original report of jan 28 using data from jan 1 to jan 22
fit1_22jan.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1_22jan))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit1_22jan_preds = as.data.frame(emmeans(fit1_22jan, ~ collection_date_num, 
                                   # by="LABORATORY", 
                                   at=list(collection_date_num=seq(date.from,
                                                                   date.to)), 
                                   type="response"), bias.adjust = TRUE, sigma = fit1_22jan.SD)
fit1_22jan_preds$collection_date = as.Date(fit1_22jan_preds$collection_date_num, origin="1970-01-01")

# best fitting separate-slope spline fit fit4, by lab
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit4))$sdcor, function (x) x^2))) 
fit_preds_bylab = as.data.frame(emmeans(fit4, ~ collection_date_num, 
                                   by="LABORATORY", 
                                   at=list(collection_date_num=seq(date.from,
                                                               date.to)), 
                                    type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_preds_bylab$collection_date = as.Date(fit_preds_bylab$collection_date_num, origin="1970-01-01")
# order labs by estimated date of introduction (intercepts)
dfemmeanslabs = as.data.frame(emmeans(fit1,~LABORATORY))
levels_BE = as.character(dfemmeanslabs$LABORATORY[order(dfemmeanslabs$emmean,decreasing=T)])
fit_preds_bylab$LABORATORY = factor(fit_preds_bylab$LABORATORY, 
                                     levels=levels_BE)


# estimated share of 501Y.V1 among currently diagnosed infections based on best fitting model
fit_preds[fit_preds$collection_date==today,]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 27              18718 0.8305255 0.01398165 Inf 0.8013776 0.8562271      2021-04-01

# estimated share of 501Y.V1 among new infections (assuming time between infection & diagnosis of 7 days)
fit_preds[fit_preds$collection_date==(today+7),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 34               18725 0.8601258 0.01411272 Inf 0.8301598 0.8855993      2021-04-08


# estimated share of 501Y.V1 among currently diagnosed infections predicted based on jan 22 data used in first report
fit1_22jan_preds[fit1_22jan_preds$collection_date==today,]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 27              18718 0.9978622 0.001998574 Inf 0.9867597 0.9996587      2021-04-01

# estimated share of 501Y.V1 among new infections (assuming time between infection & diagnosis of 7 days) predicted based on jan 22 data used in first report
fit1_22jan_preds[fit1_22jan_preds$collection_date==(today+7),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 34               18725 0.9990254 0.0009944078 Inf 0.9928372 0.9998683      2021-04-08


# estimated share of 501Y.V1 among currently diagnosed infections in different regions based on best fitting model
fit_preds_bylab[fit_preds_bylab$collection_date==today,]
# collection_date_num       LABORATORY      prob         SE  df asymp.LCL asymp.UCL collection_date
# 213                18718            Namur 0.8188628 0.038155725 Inf 0.7322974 0.8823128      2021-04-01
# 456                18718  Saint LUC - UCL 0.6375110 0.050885199 Inf 0.5333954 0.7305289      2021-04-01
# 699                18718              ULB 0.7695901 0.043186029 Inf 0.6748710 0.8435647      2021-04-01
# 942                18718 UMons - Jolimont 0.9710271 0.007510052 Inf 0.9520894 0.9826316      2021-04-01
# 1185               18718        UZ leuven 0.6672784 0.053579783 Inf 0.5557769 0.7632224      2021-04-01
# 1428               18718              UZA 0.8892086 0.022257140 Inf 0.8376718 0.9259734      2021-04-01

# estimated share of 501Y.V1 among new infections (assuming time between infection & diagnosis of 7 days) in different regions
fit_preds_bylab[fit_preds_bylab$collection_date==(today+7),]
# collection_date_num       LABORATORY      prob         SE  df asymp.LCL asymp.UCL collection_date
# 220                18725            Namur 0.8483952 0.039257961 Inf 0.7551563 0.9107133      2021-04-08
# 463                18725  Saint LUC - UCL 0.6505617 0.058828367 Inf 0.5288308 0.7559119      2021-04-08
# 706                18725              ULB 0.7877582 0.047878212 Inf 0.6797222 0.8670279      2021-04-08
# 949                18725 UMons - Jolimont 0.9841714 0.004929228 Inf 0.9709778 0.9914264      2021-04-08
# 1192               18725        UZ leuven 0.6910551 0.060673937 Inf 0.5621256 0.7964552      2021-04-08
# 1435               18725              UZA 0.9128367 0.021088844 Inf 0.8618354 0.9463080      2021-04-08


sum(tail(data_ag_byday_wide$est_n_B117, 14))/sum(tail(data_ag_byday_wide$n_pos,14)) 
# 68.28% of the samples of last 2 weeks in the dataset were estimated to be by British variant
# PS: with data 31/1 this was 15.4%  
# note: this is not the same as the estimated prop of the new infections or new diagnoses today that are of the British
# variant, which are much higher, see above)


# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50%, 75% 90% by 501Y.V1 :
fit_preds$collection_date[fit_preds[,"prob"]>=0.5][1]-7 # >50% by 9th of February [7th Feb - 10 Feb] 95% CLs
fit_preds$collection_date[fit_preds[,"asymp.UCL"]>=0.5][1]-7
fit_preds$collection_date[fit_preds[,"asymp.LCL"]>=0.5][1]-7

fit_preds$collection_date[fit_preds[,"prob"]>=0.75][1]-7 # >75% by 10th of March [7 March - 15 March] 95% CLs
fit_preds$collection_date[fit_preds[,"asymp.UCL"]>=0.75][1]-7
fit_preds$collection_date[fit_preds[,"asymp.LCL"]>=0.75][1]-7

fit_preds$collection_date[fit_preds[,"prob"]>=0.9][1]-7 # >90% by 13th of April [6 April - 24 April] 95% CLs
fit_preds$collection_date[fit_preds[,"asymp.UCL"]>=0.9][1]-7
fit_preds$collection_date[fit_preds[,"asymp.LCL"]>=0.9][1]-7



# PLOT MODEL FIT best model fit4

# plot for the whole of Belgium on response scale:

plot_fit_response = qplot(data=fit_preds, x=collection_date, y=100*prob, geom="blank") +
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
plot_fit_response


# saveRDS(plot_fit_response, file = paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_response scale.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_response scale.pdf"), width=8, height=6)

plot_fit_responseB = qplot(data=fit_preds, x=collection_date, y=100*prob, geom="blank") +
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
plot_fit_responseB


# saveRDS(plot_fit_responseB, file = paste0(".\\plots\\",dat,"\\Fig4B_fit4_binomGLMM_501YV1_Belgium_response scale.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig4B_fit4_binomGLMM_501YV1_Belgium_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4B_fit4_binomGLMM_501YV1_Belgium_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig4B_fit4_binomGLMM_501YV1_Belgium_response scale.pdf"), width=8, height=6)


# plot for the whole of Belgium on logit scale:
plot_fit_link = qplot(data=fit_preds, x=collection_date, y=prob, geom="blank") +
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
plot_fit_link


# saveRDS(plot_fit_response, file = paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_logit scale.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_logit scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit4_binomGLMM_501YV1_Belgium_logit scale.pdf"), width=8, height=6)


# plot per lab on logit scale:
n = length(levels(fit_preds_bylab$LABORATORY))
reg_cols_BE = hcl(h = seq(0, 280, length = n), l = 55, c = 200)[1:n]

plot_fit_bylab = qplot(data=fit_preds_bylab, x=collection_date, y=prob, geom="blank") +
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
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_color_manual("", values=reg_cols_BE) +
  scale_fill_manual("", values=reg_cols_BE) +
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
plot_fit_bylab

# saveRDS(plot_fit_bylab, file = paste0(".\\plots\\",dat,"\\Fig5_fit4_binomGLMM_501YV1_Belgium by lab.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig5_fit4_binomGLMM_501YV1_Belgium by lab.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig5_fit4_binomGLMM_501YV1_Belgium by lab.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig5_fit4_binomGLMM_501YV1_Belgium by lab.pdf"), width=8, height=6)




# 4.2. PLOT OF Re VALUE of WILD TYPE AND 501Y.V1, CALCULATED FROM PREDICTED SHARE OF 501Y.V1 AMONG ALL INFECTIONS & FITTED TRANSMISSION ADVANTAGE ####
# based on the fact that the overall Re is a weighted average of the Re of the individual variants

# fitted Re values based on Sciensano new confirmed case & testing data
Re_cases = read.csv(paste0(".//data//",dat,"//Re_cases.csv")) 
# Re_cases = read.csv(paste0(".//data//be_latest//Re_cases.csv")) 
# Re values calculated from instantaneous growth rate in nr of new cases 
# with instantaneous growth rate calculated as the first derivative (calculated using emtrends) to the GAM fit on new cases 
# gam(cbind(NEWCASES, totpop-NEWCASES) ~ s(DATE_NUM, bs="cs", k=32, fx=F) + 
#                                         WEEKDAY + s(log(TESTS_ALL), bs="cs", k=5, fx=F), family=binomial(cloglog), data=cases_tot) 
# and with Re calculated using R.from.r with gamma_mean=4.7, gamma_sd=2.9
# note that Re here is calculated at time of diagnosis
Re_cases$DATE = as.Date(Re_cases$DATE_NUM, origin="1970-01-01")
Re_cases$collection_date_num = Re_cases$DATE_NUM

date.from.num = min(Re_cases$collection_date_num)
date.to.num = max(Re_cases$collection_date_num)

# estimated transmission advantage of 501Y.V1 (potentially time-varying when using a binomial spline GLMM)
M_fitted = exp(as.data.frame(emtrends(fit4, ~ 1, var="collection_date_num", by="collection_date_num",
                                                 at=list(collection_date_num=seq(date.from.num,
                                                                                 date.to.num)), 
                                  mode="link"))$collection_date_num.trend*4.7)
Re_cases$M = M_fitted

# functions to calculate Re of wild type and of 501Y.V1 based on overall Re value and prop of positives that is 501Y.V1 propB117
# and transmission advantage of 501Y.V1 M

Re_wild_type = function (Re, propB117, M) {
  Re / (1-propB117+M*propB117)
}
Re_B117 = function (Re, propB117, M) {
  M*Re / (1-propB117+M*propB117)
}

Re_cases$propB117 = as.data.frame(emmeans(fit4, ~ collection_date_num, 
                                          # by="LABORATORY", 
                                          at=list(collection_date_num=seq(min(Re_cases$collection_date_num),
                                                                          max(Re_cases$collection_date_num))), 
                                          type="response"), bias.adjust = TRUE, sigma = total.SD)$prob

head(Re_cases)
Re_cases$Re_WT = Re_wild_type(Re=Re_cases$Re, propB117=Re_cases$propB117, M=Re_cases$M)
Re_cases$Re_WT_LOWER = Re_wild_type(Re=Re_cases$Re_LOWER, propB117=Re_cases$propB117, M=Re_cases$M)
Re_cases$Re_WT_UPPER = Re_wild_type(Re=Re_cases$Re_UPPER, propB117=Re_cases$propB117, M=Re_cases$M)
Re_cases$Re_B117 = Re_B117(Re=Re_cases$Re, propB117=Re_cases$propB117, M=Re_cases$M)
Re_cases$Re_B117_LOWER = Re_B117(Re=Re_cases$Re_LOWER, propB117=Re_cases$propB117, M=Re_cases$M)
Re_cases$Re_B117_UPPER = Re_B117(Re=Re_cases$Re_UPPER, propB117=Re_cases$propB117, M=Re_cases$M)

# d = as.Date(max(Re_cases$DATE))
# tag = paste("@TWenseleers\n",dat)

qplot(data=Re_cases, x=DATE, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.7), fill=I("darkgrey")) +
  geom_line() + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous(limits=c(1/3,3), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("R WAARDE VAN DE BRITSE VARIANT (rood)\nEN HET OUDE TYPE (blauw)") +
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
                  ylim=c(0.6,1.8)) +
  scale_y_continuous(breaks=seq(0.6,1.8,by=0.1), labels=seq(0.6,1.8,by=0.1)) +
  ylab("R waarde")
ggsave(file=paste0(".//plots//",dat,"//Re_cases_Re_501YV1_Re_wildtype.png"), width=7, height=5)
ggsave(file=paste0(".//plots//",dat,"//Re_cases_Re_501YV1_Re_wildtype.pdf"), width=7, height=5)
Re_cases[Re_cases$DATE==max(Re_cases$DATE),]
# DATE_NUM          r          SE       df    r_LOWER    r_UPPER       DATE      Re Re_LOWER Re_UPPER collection_date_num        M  propB117     Re_WT
# 237    18711 0.02632812 0.001408553 204.0343 0.02355093 0.02910531 2021-03-25 1.12853 1.114516 1.142652               18711 1.167686 0.7962589 0.9955965
# Re_WT_LOWER Re_WT_UPPER  Re_B117 Re_B117_LOWER Re_B117_UPPER
# 237   0.9832337    1.008055 1.162544      1.148108      1.177092


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
# separate-slopes 3 df spline model fit_ukSGTF_4 best
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
# fit_preds_bylab2 = fit_preds_bylab
# colnames(fit_preds_bylab2)[2] = "REGION"
fit_preds2 = fit_preds
fit_preds2$REGION = "Belgium"
preds_UK_plus_BE = rbind(fit_ukSGTF_4_preds, fit_preds2) # fit_preds_bylab2
n = length(levels(fit_ukSGTF_4_preds$REGION))
reg_cols = hcl(h = seq(300, 20, length = n), l = 50, c = 200)[1:n]
reg_cols_UK_BE = c(reg_cols, "steelblue") # muted(reg_cols_BE, l=65, c=200)
# data_ag_wide2 = data_ag_wide
# colnames(data_ag_wide2)[2] = "REGION"
data_ag_byday_wide2 = data_ag_byday_wide
# colnames(data_ag_wide2)[2] = "REGION"
data_ag_byday_wide2$REGION = "Belgium"

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
    xlim=c(as.Date("2020-10-01"),as.Date("2021-04-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols_UK_BE) +
  scale_fill_manual("", values=reg_cols_UK_BE) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117*100, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  geom_point(data=data_ag_byday_wide2, 
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
  ggtitle("SPREAD OF VARIANT 501Y.V1 IN THE\nUK & BELGIUM BASED ON S DROPOUT DATA") +
  theme(plot.title = element_text(hjust = 0.5))
plot_UK_SGTF_BE_response
# saveRDS(plot_UK_SGTF_BE_response, file = paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig6_UK plus BE S dropout data_response scale.pdf"), width=8, height=6)


# 5.2. DATA DENMARK: SEQUENCING DATA ####

# Data source: Danish Covid-19 Genome Consortium & the Statens Serum Institut, https://www.covid19genomics.dk/statistics

data_denmark = read.csv(".//data//dk//data_denmark_20210319.csv", sep=";", dec=",")
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
data_denmark = data_denmark[data_denmark$date>="2020-11-19",] # we use data from 19th of Nov onwards

data_denmark_whole = data_denmark[data_denmark$Region=="Whole Denmark",]
data_denmark = data_denmark[data_denmark$Region!="Whole Denmark",]
levels_DK = c("Syddanmark","Sjælland","Nordjylland","Hovedstaden","Midtjylland")
data_denmark$Region = factor(data_denmark$Region, levels=levels_DK)
range(data_denmark$date) # "2020-11-19" "2021-03-11"

fit_denmark1 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + Region + scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark2 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + Region * scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark3 = glmer(cbind(n_B117,total-n_B117) ~ (1|Region/obs) + scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark4 = glmer(cbind(n_B117,total-n_B117) ~ (date_num||Region/obs) + scale(date_num), family=binomial(logit), data=data_denmark)

BIC(fit_denmark1, fit_denmark2, fit_denmark3, fit_denmark4)
# df      BIC
# fit_denmark1  7 681.1196
# fit_denmark2 11 648.3259
# fit_denmark3  4 675.8928
# fit_denmark4  6 684.6456

summary(fit_denmark2)

# common-slope model fit_denmark2 with nested random intercepts fits best

#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model fit_denmark2, we get
dk_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_denmark2, ~ 1, var="date_num"))[,-c(3,4)] 
colnames(dk_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
dk_growthrates_avg_B117vsallother = M.from.delta_r_df(dk_growthrates_avg_B117vsallother)
dk_growthrates_avg_B117vsallother
# 1 logistic_growth_rate  asymp.LCL  asymp.UCL        M    M.LCL    M.UCL
# 1 overall            0.08215316 0.07792146 0.08638487 1.471261 1.442288 1.500816


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
# write_csv(data_switzerland, file=".//data//ch//data_switzerland_20210307.csv")

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
range(data_switzerland$date) # "2020-11-02" "2021-03-04"

fit_switerland1 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + lab + scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland2 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + lab * scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland3 = glmer(cbind(n_B117,total-n_B117) ~ (1|lab/obs) + scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland4 = glmer(cbind(n_B117,total-n_B117) ~ (date_num||lab/obs) + scale(date_num), family=binomial(logit), data=data_switzerland)
BIC(fit_switerland1, fit_switerland2, fit_switerland3, fit_switerland4)
#                 df      BIC
# fit_switerland1  6 861.0248
# fit_switerland2  9 865.8782
# fit_switerland3  4 870.3925
# fit_switerland4  6 880.7673

# fit fit_switerland1 best

summary(fit_switerland1)


#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model fit_switerland1, we get
ch_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_switerland1, ~ 1, var="date_num"))[,-c(3,4)] 
colnames(ch_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
ch_growthrates_avg_B117vsallother = M.from.delta_r_df(ch_growthrates_avg_B117vsallother)
ch_growthrates_avg_B117vsallother
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL     M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall            0.0837439 0.0788791 0.0886087 1.482302 1.448794 1.516585



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




# 5.4. DATA USA : S-GENE TARGET FAILURE DATA ####

# Data source: Helix® COVID-19 Surveillance, https://github.com/myhelix/helix-covid19db
# see preprint https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1 & https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data

us_data = read.csv("https://github.com/myhelix/helix-covid19db/raw/master/counts_by_state.csv")
# write.csv(us_data, file=".//data//us//data_us_20210224.csv", row.names=F)

us_data$collection_date = as.Date(us_data$collection_date)
us_data$collection_date_num = as.numeric(us_data$collection_date)
us_data$obs = factor(1:nrow(us_data))
# us_data = us_data[us_data$state %in% sel_states,]
us_data$state = factor(us_data$state)

range(us_data$collection_date) # "2020-09-05" "2021-04-22"

us_data2 = us_data
us_data2 = us_data2[us_data2$sequenced_SGTF!=0,]
us_data2$propB117 = us_data2$B117/us_data2$sequenced_SGTF
us_data2 = us_data2[us_data2$propB117<=1,]
range(us_data2$propB117)

range(us_data2$collection_date) # "2020-12-19" "2021-04-06"

fit_us_propB117amongSGTF = glmer(cbind(B117, sequenced_SGTF-B117) ~ (1|state)+scale(collection_date_num), 
                                 family=binomial(logit), data=us_data2)
BIC(fit_us_propB117amongSGTF) 

# implied growth rate advantage of B.1.1.7 over other earlier strains showing S dropout:
as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)]
#   collection_date_num.trend  asymp.LCL asymp.UCL
# 1                 0.06099678 0.05686935 0.06512421

# with a generation time of 4.7 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of B.1.1.7 over other strains showing S dropout of
exp(4.7*as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.332004  1.306414  1.358096

# with a generation time of 5.5 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of B.1.1.7 over other strains showing S dropout of
exp(5.5*as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.398615  1.367223  1.430728


# FIT FOR WHOLE US + PLOT

fitted_truepos = predict(fit_us_propB117amongSGTF, newdat=us_data, type="response", re.form=NA) 
# fitted true positive rate, ie prop of S dropout samples that are B.1.1.7 for dates & states in helix_sgtf

us_data$est_n_B117 = us_data$all_SGTF*fitted_truepos # estimated nr of B.1.1.7 samples
us_data$propB117 = us_data$est_n_B117/us_data$positive
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit1)
fit_us1 = glmer(cbind(est_n_B117, positive-est_n_B117) ~ (1|state/obs)+scale(collection_date_num), 
                family=binomial(logit), data=us_data) # random intercepts by state
fit_us2 = glmer(cbind(est_n_B117, positive-est_n_B117) ~ (collection_date_num||state/obs)+scale(collection_date_num), 
                family=binomial(logit), data=us_data, control=glmersettings2) # random intercepts+slopes by state, with uncorrelated intercepts & slopes
BIC(fit_us1, fit_us2) # random intercept model fit_us1 is best
# df      BIC
# fit_us1  4 8309.006
# fit_us2  6 8961.308
summary(fit_us1)

#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all states, using the most parsimonious model fit_us1, we get
us_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_us1, ~ 1, var="collection_date_num"))[,-c(3,4)] 
colnames(us_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
us_growthrates_avg_B117vsallother = M.from.delta_r_df(us_growthrates_avg_B117vsallother)
us_growthrates_avg_B117vsallother
# 1 logistic_growth_rate asymp.LCL  asymp.UCL        M  M.LCL    M.UCL
# 1 overall           0.06815698 0.06698542 0.06932853 1.377593 1.370028 1.385199

# with a generation time of 4.7 days this would translate to a multiplicative effect on Rt of
exp(4.7*as.data.frame(emtrends(fit_us1, ~ 1, var="collection_date_num"))[c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.377593  1.370028  1.385199 

# with a generation time of 5.5 days this would translate to a multiplicative effect on Rt of
exp(5.5*as.data.frame(emtrends(fit_us1, ~ 1, var="collection_date_num"))[c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.454793  1.445449  1.464197



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
# sel_states=c("MI","MN","NY","TX","GA","FL","CA","MA","PA","IN","IL","AZ")
sel_states=c("CA","CO","CT","FL","IL","MD","MA","MI","MN","NJ","NY","PA")
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
fit_us_preds2$state = factor(fit_us_preds2$state, levels=
                               fit_us_preds2[fit_us_preds2$collection_date==today,"state"][order(fit_us_preds2[fit_us_preds2$collection_date==today,"prob"], decreasing=T)])

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
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(as.Date("2020-12-01"), as.Date("2021-04-30")), 
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
  xlab("") +
  labs(tag = "@TWenseleers\ndata Helix") +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
# ggtitle("US") +
# theme(plot.title = element_text(hjust = 0.5))
plot_us
ggsave(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.png"), width=9, height=7)


# PLOT MODEL FIT (response scale)
plot_us_response = qplot(data=fit_us_preds2, x=collection_date, y=prob*100, geom="blank") +
  facet_wrap(~state, nrow=3) +
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
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(as.Date("2020-12-01"), as.Date("2021-04-30")), 
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
  xlab("") +
  labs(tag = "@TWenseleers\ndata Helix") +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# +
# theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
# ggtitle("US") +
# theme(plot.title = element_text(hjust = 0.5))

plot_us_response
ggsave(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state_response.png"), width=9, height=7)






# 5.5. MULTIPANEL PLOT INTERNATIONAL COMPARISONS ####

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
fit_us_preds3 = fit_us_preds3[fit_us_preds3$REGION %in% c("FL","CA","TX","GA","MI"),]
fit_us_preds3$REGION = factor(fit_us_preds3$REGION, levels=c("FL","CA","TX","GA","MI"), labels=c("Florida","California","Texas","Georgia","Michigan"))
fit_us_preds3 = fit_us_preds3[,c("date_num","REGION","prob","SE","df","asymp.LCL","asymp.UCL","date","country")]
fit_be_preds = fit_preds
fit_be_preds$REGION = "Belgium"
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
data_us2 = data_us2[data_us2$REGION %in% c("FL","CA","TX","GA","MI"),]
data_us2$REGION = factor(data_us2$REGION, levels=c("FL","CA","TX","GA","MI"), labels=c("Florida","California","Texas","Georgia","Michigan"))

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
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
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
    xlim=c(as.Date("2020-09-01"),as.Date("2021-04-30")),
    ylim=c(ymin,ymax+1E-7), 
    expand=FALSE) 
# ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
# theme(plot.title = element_text(hjust = 0.5))
plot_international

# saveRDS(plot_international, file = paste0(".\\plots\\",dat,"\\Fig7_international_data_UK_CH_USA_DK_BE.rds"))
# graph2ppt(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE.pptx"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE.png"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE.pdf"), width=9, height=7)
# ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE.pdf"), width=9, height=7)




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
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
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
    xlim=c(as.Date("2020-09-01"),as.Date("2021-04-30")),
    ylim=c(0,100), expand=c(0,0))
# +
# ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
# theme(plot.title = element_text(hjust = 0.5))
plot_international_response

# saveRDS(plot_international_response, file = paste0(".\\plots\\",dat,"\\Fig7_international_data_UK_CH_USA_DK_BE_response.rds"))
# graph2ppt(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE_response.pptx"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE_response.png"), width=9, height=7)
ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE_response.pdf"), width=9, height=7)
# ggsave(file = paste0(".\\plots\\",dat,"\\Fig7_internat_data_UK_CH_USA_DK_BE_response.pdf"), width=9, height=7)


plot_us2 = plot_us + coord_cartesian(xlim=c(as.Date("2020-11-01"), as.Date("2021-03-31")),
                                     ylim=c(0.001,99.9), expand=c(0,0)) # + ggtitle("SPREAD OF VARIANT B.1.1.7 IN THE US")
plot_us2

# saveRDS(plot_us2, file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.rds"))
# graph2ppt(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.pptx"), width=9, height=7)
# ggsave(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.png"), width=9, height=7)
# ggsave(file = paste0(".\\plots\\",dat,"\\Fig8_US_data_by state.pdf"), width=9, height=7)



# 5.6. PLOTS OF ESTIMATED NR OF NEW INFECTIONS BY DIFFERENT VOCs VS. WILD TYPE BASED ON LOGISTIC FIT TO S DROPOUT DATA & MULTINOMIAL FIT TO BASELINE SURVEILLANCE SEQ DATA ####

# fit POISSON GAM TO TOTAL SCIENSANO CASE DATA, correcting for weekday & testing intensity ####

source("scripts/downloadData.R") # download latest data with new confirmed cases per day from Sciensano website, code adapted from https://github.com/JoFAM/covidBE_analysis by Joris Meys
range(cases_tot$DATE) # "2020-03-01" "2021-04-28"

cases_tot = cases_tot[cases_tot$DATE>=as.Date("2020-08-01"),]

# smooth out weekday effects in case nrs using GAM & correct for unequal testing intensity
fit_cases_BE = gam(CASES ~ s(DATE_NUM, bs="cs", k=20, fx=F) + 
                     WEEKDAY + BANKHOLIDAY +
                     s(TESTS_ALL, bs="cs", k=8, fx=F),
                   family=poisson(log), data=cases_tot,
) 
BIC(fit_cases_BE)
# expected marginal mean cases for average weekdata at mean testing effort of last week
df = as.data.frame(emmeans(fit_cases_BE, ~ DATE_NUM, 
                           at = list(DATE_NUM=cases_tot$DATE_NUM,
                                     BANKHOLIDAY = factor("no"),
                                     TESTS_ALL=max(cases_tot$TESTS_ALL)
                           ),
                           type="response"))
colnames(df)[2] = "CASES_SMOOTH"
colnames(df)[5] = "CASES_LCL"
colnames(df)[6] = "CASES_UCL"
cases_tot$CASES_SMOOTH = df$CASES_SMOOTH # Poisson GAM smoothed data
cases_tot$CASES_LCL = df$CASES_LCL
cases_tot$CASES_UCL = df$CASES_UCL

# fit Poisson GAM to cases by districts where national testing labs are based ####
# TO DO
rawmunicipalities$cases[rawmunicipalities$cases=="<5"] = "2.5"
rawmunicipalities$cases = as.numeric(rawmunicipalities$cases)
sel_districts = c("Antwerpen","Brussel-Hoofdstad","Leuven","Gent","Arrondissement Luik","Arrondissement Namen","Arrondissement Bergen")
cases_distr = rawmunicipalities[rawmunicipalities$DISTRICT %in% sel_districts,]
cases_distr$DISTRICT = factor(cases_distr$DISTRICT, levels=sel_districts, labels=c("Antwerpen","Brussel","Leuven","Gent","Luik","Namen","Bergen"))
cases_distr$DATE = as.Date(cases_distr$DATE)
cases_distr$DATE_NUM = as.numeric(cases_distr$DATE)
cases_distr = cases_distr[with(cases_distr, order(DISTRICT,DATE)),]
cases_distr$newcases = unlist(lapply(unique(cases_distr$DISTRICT), function(distr) c(0,diff(cases_distr$cases[cases_distr$DISTRICT==distr]))))
cases_distr$newcases[cases_distr$newcases<0] = 0
range(rawmunicipalities$DATE)
cases_distr$DISTRDATE = interaction(cases_distr$DISTRICT, cases_distr$DATE)
cases_distr$WEEKDAY = as.factor(weekdays(cases_distr$DATE))
cases_distr = cases_distr[complete.cases(cases_distr),]
fit_cases_distr = gam(newcases ~ # s(DATE_NUM, bs="cs", k=30, fx=T) + 
                       s(DATE_NUM, bs="cs", k=30, by=DISTRICT, fx=T) + 
                        DISTRICT + WEEKDAY, # evt  + s(log(TESTS_ALL+1), k=5, fx=T)
                     family=poisson, data=cases_distr) 
BIC(fit_cases_distr)
cases_distr$newcases_fitted = predict(fit_cases_distr, type="response")

emmeans_cases_distr = as.data.frame(emmeans(fit_cases_distr, ~ DATE_NUM+DISTRICT, 
                                                 at = list(DATE_NUM=unique(cases_distr$DATE_NUM))
                                                 ),
                                                 type="response")
emmeans_cases_distr$DATE = as.Date(emmeans_cases_distr$DATE_NUM, origin="1970-01-01")
qplot(data=cases_distr, x=DATE, y=newcases, group=DISTRICT, colour=DISTRICT, fill=DISTRICT, geom="col")
qplot(data=cases_distr, x=DATE, y=newcases_fitted, group=DISTRICT, colour=DISTRICT, fill=DISTRICT, geom="col")
qplot(data=emmeans_cases_distr, x=DATE, y=rate, group=DISTRICT, colour=DISTRICT, fill=DISTRICT, geom="area")

head(data_ag)
data_ag_wide


# fit Poisson GAM to cases by province where national testing labs are based ####
# TO DO
cases_prov = rawcases[rawcases$PROVINCE!="All"&rawcases$REGION!="Belgium"&rawcases$SEX=="All"&rawcases$AGEGROUP!="unknown"&rawcases$AGEGROUP=="All"&rawcases$PROVINCE!="unknown",]
sel_provs = c("Antwerpen","Brussels","VlaamsBrabant","OostVlaanderen","Liège","Namur","Hainaut")
cases_prov = cases_prov[cases_prov$PROVINCE %in% sel_provs,] 
cases_prov$DATE = as.Date(cases_prov$DATE)
cases_prov$DATE_NUM = as.numeric(cases_prov$DATE)
cases_prov$PROVDATE = interaction(cases_prov$PROVINCE,cases_prov$DATE)
cases_prov$TESTS_ALL = tests_prov$TESTS_ALL[match(cases_prov$PROVDATE,tests_prov$PROVDATE)]
cases_prov$PROVINCE = factor(cases_prov$PROVINCE, levels=c("Antwerpen","Brussels","VlaamsBrabant","OostVlaanderen","Liège","Namur","Hainaut"),
                             labels=c("Antwerpen","Brussel","Vlaams Brabant","Oost Vlaanderen","Luik","Namen","Henegouwen"))
cases_prov$WEEKDAY = as.factor(weekdays(cases_prov$DATE))
cases_prov = cases_prov[complete.cases(cases_prov),]
fit_cases_prov = gam(CASES ~ # s(DATE_NUM, bs="cs", k=30, fx=T) + 
                       s(DATE_NUM, bs="cs", k=30, by=PROVINCE, fx=T) + 
                       PROVINCE + WEEKDAY, # + s(log(TESTS_ALL+1), k=5, fx=T), 
                     family=poisson, data=cases_prov) 
BIC(fit_cases_prov) # 46345.07
cases_prov$CASES_fitted = predict(fit_cases_prov, type="response")

emmeans_cases_prov = as.data.frame(emmeans(fit_cases_prov, ~ DATE_NUM+PROVINCE, 
                                            at = list(DATE_NUM=unique(cases_distr$DATE_NUM))
                                          ), type="response")
emmeans_cases_prov$DATE = as.Date(emmeans_cases_prov$DATE_NUM, origin="1970-01-01")
qplot(data=cases_prov, x=DATE, y=CASES, group=PROVINCE, colour=PROVINCE, fill=PROVINCE, geom="col")
qplot(data=cases_prov, x=DATE, y=CASES_fitted, group=PROVINCE, colour=PROVINCE, fill=PROVINCE, geom="col")
qplot(data=emmeans_cases_prov, x=DATE, y=rate, group=PROVINCE, colour=PROVINCE, fill=PROVINCE, geom="area")

# head(data_ag)
# data_ag_wide

# fit of prop of 501Y.V1 by district based on S dropout data
# TO DO
data_ag_wide$DISTRICT = factor(data_ag_wide$LABORATORY, levels=c("Namur","Saint LUC - UCL","ULB","UMons - Jolimont","UZ leuven","UZA"),
                               labels=c("Namen","Brussel","Brussel","Bergen","Leuven","Antwerpen"))
fit4_distr = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+ns(collection_date_num,df=2)*DISTRICT, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings) # with 2 df natural cubic spline term ifo date * DISTRICT
prop501YV1_distr = data.frame(emmeans(fit4_distr, ~ collection_date_num+DISTRICT, at=list(collection_date_num=as.numeric(cases_tot$DATE)), 
                                                    type="response", df=NA)) # proportions of 3 VOCs as fit based on mulinomial fit


# hospitalisations & ICU admissions by region/city
# TO DO


# plot for Belgium (BASED ON MULTINOMIAL FIT TO BASELINE SURVEILLANCE SEQUENCING DATA & POISSON GAM FIT TO SCIENSANO TOTAL CASE DATA) ####

props = data.frame(emmeans(be_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=as.numeric(cases_tot$DATE)), 
                           mode="prob", df=NA)) # proportions of 3 VOCs as fit based on mulinomial fit
cases_tot1 = cases_tot
cases_tot1$prop501YV1 =  props[props$variant=="501Y.V1","prob"]
cases_tot1$prop501YV1[cases_tot1$DATE<as.Date("2020-12-01")] =  0
cases_tot1$prop501YV1_LCL =  props[props$variant=="501Y.V1","asymp.LCL"] 
cases_tot1$prop501YV1_LCL[cases_tot1$DATE<as.Date("2020-12-01")] =  0
cases_tot1$prop501YV1_UCL =  props[props$variant=="501Y.V1","asymp.UCL"] 
cases_tot1$prop501YV1_UCL[cases_tot1$DATE<as.Date("2020-12-01")] =  0

cases_tot1$prop501YV2 =  props[props$variant=="501Y.V2","prob"]
cases_tot1$prop501YV2[cases_tot1$DATE<as.Date("2020-12-01")] =  0
cases_tot1$prop501YV2_LCL =  props[props$variant=="501Y.V2","asymp.LCL"] 
cases_tot1$prop501YV2_LCL[cases_tot1$DATE<as.Date("2020-12-01")] =  0
cases_tot1$prop501YV2_UCL =  props[props$variant=="501Y.V2","asymp.UCL"] 
cases_tot1$prop501YV2_UCL[cases_tot1$DATE<as.Date("2020-12-01")] =  0

cases_tot1$prop501YV3 =  props[props$variant=="501Y.V3","prob"]
cases_tot1$prop501YV3[cases_tot1$DATE<as.Date("2020-12-01")] =  0
cases_tot1$prop501YV3_LCL =  props[props$variant=="501Y.V3","asymp.LCL"] 
cases_tot1$prop501YV3_LCL[cases_tot1$DATE<as.Date("2020-12-01")] =  0
cases_tot1$prop501YV3_UCL =  props[props$variant=="501Y.V3","asymp.UCL"] 
cases_tot1$prop501YV3_UCL[cases_tot1$DATE<as.Date("2020-12-01")] =  0

data_cases_BE_501YV1 = data.frame(cases_tot1, variant="501Y.V1")
data_cases_BE_501YV1$CASES = cases_tot1$CASES*cases_tot1$prop501YV1
data_cases_BE_501YV1$CASES_SMOOTH = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV1
data_cases_BE_501YV1$CASES_LCL = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV1_LCL
data_cases_BE_501YV1$CASES_UCL = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV1_UCL

data_cases_BE_501YV2 = data.frame(cases_tot1, variant="501Y.V2")
data_cases_BE_501YV2$CASES = cases_tot1$CASES*cases_tot1$prop501YV2
data_cases_BE_501YV2$CASES_SMOOTH = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV2
data_cases_BE_501YV2$CASES_LCL = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV2_LCL
data_cases_BE_501YV2$CASES_UCL = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV2_UCL

data_cases_BE_501YV3 = data.frame(cases_tot1, variant="501Y.V3")
data_cases_BE_501YV3$CASES = cases_tot1$CASES*cases_tot1$prop501YV3
data_cases_BE_501YV3$CASES_SMOOTH = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV3
data_cases_BE_501YV3$CASES_LCL = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV3_LCL
data_cases_BE_501YV3$CASES_UCL = cases_tot1$CASES_SMOOTH*cases_tot1$prop501YV3_UCL

data_cases_BE_wildtype = data.frame(cases_tot1, variant="wild type")
data_cases_BE_wildtype$CASES = cases_tot1$CASES*(1-(cases_tot1$prop501YV1+cases_tot1$prop501YV2+cases_tot1$prop501YV3))
data_cases_BE_wildtype$CASES_SMOOTH = cases_tot1$CASES_SMOOTH*(1-(cases_tot1$prop501YV1+cases_tot1$prop501YV2+cases_tot1$prop501YV3))
data_cases_BE_wildtype$CASES_LCL = cases_tot1$CASES_SMOOTH*(1-(cases_tot1$prop501YV1_LCL+cases_tot1$prop501YV2_LCL+cases_tot1$prop501YV3_LCL))
data_cases_BE_wildtype$CASES_UCL = cases_tot1$CASES_SMOOTH*(1-(cases_tot1$prop501YV1_UCL+cases_tot1$prop501YV2_UCL+cases_tot1$prop501YV3_UCL))

data_cases_BE_total = data.frame(cases_tot1, variant="total")

data_cases_BE = rbind(data_cases_BE_501YV1, data_cases_BE_501YV2, data_cases_BE_501YV3, data_cases_BE_wildtype, data_cases_BE_total)
data_cases_BE$variant = factor(data_cases_BE$variant, levels=c("wild type", "501Y.V1", "501Y.V2", "501Y.V3", "total"), 
                               labels=c("wild type","501Y.V1 (UK)", "501Y.V2 (South Africa)", "501Y.V3 (Brazil)", "total"))
data_cases_BE_NL = data_cases_BE
data_cases_BE_NL$variant = factor(data_cases_BE$variant, levels=c("wild type","501Y.V1 (UK)", "501Y.V2 (South Africa)", "501Y.V3 (Brazil)", "total"), 
                               labels=c("oude varianten","Britse variant", "Zuid-Afrikaanse variant", "Braziliaanse variant", "totaal"))

d = as.Date(max(data_cases_BE$DATE))
tag = paste("@TWenseleers\ndata Sciensano & Emmanuel André\n",d)

data_cases_BE2 = data_cases_BE
data_cases_BE2[data_cases_BE2$CASES<=0.1] = 0.101
data_cases_BE2[data_cases_BE2$CASES_LCL<=0.1] = 0.101

qplot(data=data_cases_BE2[data_cases_BE2$DATE>=as.Date("2020-09-01"),], # data_cases_BE$variant!="total"&
      x=DATE, y=CASES_SMOOTH, group=variant, colour=variant, fill=variant, geom="blank") +
  # geom_area(position="stack") +
  geom_line(lwd=I(0.75)) +
  geom_ribbon(aes(ymin=CASES_LCL, ymax=CASES_UCL), alpha=I(0.3), colour=NA) + # TO DO: fix bug in confidence intervals
  # geom_smooth(lwd=I(0.75), se=FALSE, method = "gam", formula = y ~ s(x, bs = "cs", k=15)) +
  scale_colour_manual("", values=c("deepskyblue2","red","blue","green3", "grey40")) +
  scale_fill_manual("", values=c("deepskyblue2","red","blue","green3", "grey40")) +
  ylab("Estimated new infections per day") + xlab("") + ggtitle("Estimated new infections by UK, South African & Brazilian variants in Belgium") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  scale_y_continuous(breaks=10^seq(0,4), minor_breaks=0, expand=FALSE, limits=c(0.1,3E4)) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), expand=c(0,0)) +
  coord_trans(y="log10", ylim=c(10,max(data_cases_BE$CASES_SMOOTH)), expand=FALSE) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_501YV2_501YV3_wildtype_baseline surveillance_log10 scale_smoothed.png"), width=7, height=5)
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_501YV2_501YV3_wildtype_baseline surveillance_log10 scale_smoothed.pdf"), width=7, height=5)

qplot(data=data_cases_BE[data_cases_BE$variant!="total",], 
                                    x=DATE, y=CASES, group=variant, colour=variant, fill=variant, geom="blank") +
  geom_area(position="stack") +
  scale_colour_manual("", values=c("grey75","red","blue","green3")) +
  scale_fill_manual("", values=c("grey75","red","blue","green3")) +
  ylab("New infections per day") + xlab("") + ggtitle("Estimated new infections by UK, South African & Brazilian variants in Belgium") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), expand=c(0,0)) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_501YV2_501YV3_wildtype_baseline surveillance_stacked_raw.png"), width=7, height=5)
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_501YV2_501YV3_wildtype_baseline surveillance_stacked_raw.pdf"), width=7, height=5)



# plot for Belgium (BASED ON LOGISTIC FIT TO S DROPOUT DATA & POISSON GAM FIT TO SCIENSANO TOTAL CASE DATA) ####

cases_tot2 = cases_tot

df_fit4_emmeans = data.frame(emmeans(fit4, ~ collection_date_num, at=list(collection_date_num=as.numeric(cases_tot2$DATE)), 
                                  type="response", df=NA)) # proportion of 501Y.V1 based on S dropout data
cases_tot2$prop501YV1 =  df_fit4_emmeans$prob 
cases_tot2$prop501YV1_LCL =  df_fit4_emmeans$asymp.LCL 
cases_tot2$prop501YV1_UCL =  df_fit4_emmeans$asymp.UCL 

data_cases_BE_501YV1 = data.frame(cases_tot2, variant="501Y.V1")
data_cases_BE_501YV1$CASES = cases_tot2$CASES*cases_tot2$prop501YV1
data_cases_BE_501YV1$CASES_SMOOTH = cases_tot2$CASES_SMOOTH*cases_tot2$prop501YV1
data_cases_BE_501YV1$CASES_LCL = cases_tot2$CASES_SMOOTH*cases_tot2$prop501YV1_LCL
data_cases_BE_501YV1$CASES_UCL = cases_tot2$CASES_SMOOTH*cases_tot2$prop501YV1_UCL

data_cases_BE_wildtype = data.frame(cases_tot2, variant="wild type")
data_cases_BE_wildtype$CASES = cases_tot2$CASES*(1-cases_tot2$prop501YV1)
data_cases_BE_wildtype$CASES_SMOOTH = cases_tot2$CASES_SMOOTH*(1-cases_tot2$prop501YV1)
data_cases_BE_wildtype$CASES_LCL = cases_tot2$CASES_SMOOTH*(1-cases_tot2$prop501YV1_LCL)
data_cases_BE_wildtype$CASES_UCL = cases_tot2$CASES_SMOOTH*(1-cases_tot2$prop501YV1_UCL)

data_cases_BE_total = data.frame(cases_tot2, variant="total")

data_cases_BE = rbind(data_cases_BE_501YV1, data_cases_BE_wildtype, data_cases_BE_total)
data_cases_BE$variant = factor(data_cases_BE$variant, levels=c("wild type", "501Y.V1", "total"), 
                               labels=c("wild type & other variants","501Y.V1", "total"))

d = as.Date(max(data_cases_BE$DATE))
tag = paste("@TWenseleers\ndata Sciensano & Emmanuel André\n",d)

qplot(data=data_cases_BE[data_cases_BE$DATE>=as.Date("2020-09-01"),], # data_cases_BE$variant!="total"&
      x=DATE, y=CASES_SMOOTH, group=variant, colour=variant, fill=variant, geom="blank") +
  # geom_area(position="stack") +
  geom_line(lwd=I(0.75)) +
  geom_ribbon(aes(ymin=CASES_LCL, ymax=CASES_UCL), alpha=I(0.3), colour=NA) +
  # geom_smooth(lwd=I(0.75), se=FALSE, method = "gam", formula = y ~ s(x, bs = "cs", k=15)) +
  scale_colour_manual("", values=c("deepskyblue2","red","grey40")) +
  scale_fill_manual("", values=c("deepskyblue2","red","grey40")) +
  ylab("Estimated new confirmed cases per day") + xlab("") + ggtitle("Estimated infections by UK SARS-CoV2 variant 501Y.V1\nin Belgium (S dropout data)") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) +
  coord_trans(y="log10", ylim=c(100,2E4), expand=FALSE) +
  # coord_trans(y="log10", expand=FALSE, clip="off", ylim=c(101,3E4)) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  scale_y_continuous(breaks=10^seq(0,4), minor_breaks = 0, expand=c(0,0)) 
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_wildtype_S dropout_BE_log10 scale.png"), width=7, height=5)

qplot(data=data_cases_BE[data_cases_BE$variant!="total",], 
      x=DATE, y=CASES, group=variant, colour=variant, fill=variant, geom="blank") +
  geom_area(position="stack") +
  scale_colour_manual("", values=c("grey75","red2")) +
  scale_fill_manual("", values=c("grey75","red2")) +
  ylab("Estimated new confirmed cases per day") + xlab("") + ggtitle("Estimated infections by UK SARS-CoV2 variant 501Y.V1\nin Belgium (S dropout data)") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), expand=c(0,0)) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_wildtype_S dropout_BE_stacked_raw.png"), width=7, height=5)

qplot(data=data_cases_BE[data_cases_BE$variant!="total",], 
      x=DATE, y=CASES_SMOOTH, group=variant, colour=variant, fill=variant, geom="blank") +
  geom_area(position="stack") +
  scale_colour_manual("", values=c("grey75","red2")) +
  scale_fill_manual("", values=c("grey75","red2")) +
  ylab("Estimated new confirmed cases per day") + xlab("") + ggtitle("Estimated infections by UK SARS-CoV2 variant 501Y.V1\nin Belgium (S dropout data)") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), expand=c(0,0)) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_wildtype_S dropout_BE_stacked.png"), width=7, height=5)



# plot for Florida ####

us_data_by_state = read.csv("https://github.com/nytimes/covid-19-data/raw/master/us-states.csv")
us_data_by_state$date = as.Date(us_data_by_state$date)
us_data_by_state$state = factor(us_data_by_state$state, 
                                levels=c("Washington","Illinois","California",
                                         "Arizona","Massachusetts","Wisconsin",
                                         "Texas","Nebraska","Utah","Oregon",
                                         "Florida","New York","Rhode Island",
                                         "Georgia","New Hampshire","North Carolina",
                                         "New Jersey","Colorado","Maryland","Nevada",
                                         "Tennessee","Hawaii","Indiana","Kentucky","Minnesota",
                                         "Oklahoma","Pennsylvania","South Carolina","District of Columbia",
                                         "Kansas","Missouri","Vermont","Virginia","Connecticut",
                                         "Iowa","Louisiana","Ohio","Michigan","South Dakota",
                                         "Arkansas","Delaware","Mississippi","New Mexico","North Dakota",
                                         "Wyoming","Alaska","Maine","Alabama","Idaho","Montana",
                                         "Puerto Rico","Virgin Islands","Guam","West Virginia","Northern Mariana Islands"))
data_florida = us_data_by_state[us_data_by_state$state=="Florida",]
data_florida$newcases = c(0,diff(data_florida$cases))
data_florida$newcases[data_florida$newcases<0] = 0
# plot(data_florida$newcases)
data_florida$propB117 =  data_us2$propB117[match(interaction(data_florida$date, data_florida$state),
                                                    interaction(data_us2$date,data_us2$REGION))]
data_florida = data_florida[!is.na(data_florida$propB117),]

data_florida_501YV1 = data.frame(data_florida, variant="501Y.V1")
data_florida_501YV1$newcases = data_florida_501YV1$newcases*data_florida_501YV1$propB117
data_florida_wildtype = data.frame(data_florida, variant="wild type")
data_florida_wildtype$newcases = data_florida_wildtype$newcases*(1-data_florida_wildtype$propB117)
data_florida = rbind(data_florida_501YV1, data_florida_wildtype)
data_florida$variant = factor(data_florida$variant, levels=c("wild type", "501Y.V1"), labels=c("wild type","501Y.V1 (British)"))

d = as.Date(max(data_florida$date))
tag = paste("@TWenseleers\n                           data Helix & NYT\n",d)

plotcases501YV1_Sdropout_Florida = qplot(data=data_florida[data_florida$variant!="total",], x=date, y=newcases, group=variant, colour=variant, fill=variant, geom="blank") +
  geom_area(position="stack") +
  scale_colour_manual("", values=c("grey75","red2")) +
  scale_fill_manual("", values=c("grey75","red2")) +
  ylab("New confirmed cases") + xlab("") + ggtitle("Estimated infections by UK SARS-CoV2 variant 501Y.V1\nin Florida (S dropout data)") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-11-01"),NA), expand=c(0,0))
plotcases501YV1_Sdropout_Florida
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_vs_wildtype_S dropout_Florida.png"), width=7, height=5)




# plot for Michigan ####

us_data_by_state = read.csv("https://github.com/nytimes/covid-19-data/raw/master/us-states.csv")
us_data_by_state$date = as.Date(us_data_by_state$date)
us_data_by_state$state = factor(us_data_by_state$state, 
                                levels=c("Washington","Illinois","California",
                                         "Arizona","Massachusetts","Wisconsin",
                                         "Texas","Nebraska","Utah","Oregon",
                                         "Florida","New York","Rhode Island",
                                         "Georgia","New Hampshire","North Carolina",
                                         "New Jersey","Colorado","Maryland","Nevada",
                                         "Tennessee","Hawaii","Indiana","Kentucky","Minnesota",
                                         "Oklahoma","Pennsylvania","South Carolina","District of Columbia",
                                         "Kansas","Missouri","Vermont","Virginia","Connecticut",
                                         "Iowa","Louisiana","Ohio","Michigan","South Dakota",
                                         "Arkansas","Delaware","Mississippi","New Mexico","North Dakota",
                                         "Wyoming","Alaska","Maine","Alabama","Idaho","Montana",
                                         "Puerto Rico","Virgin Islands","Guam","West Virginia","Northern Mariana Islands"))
data_michigan = us_data_by_state[us_data_by_state$state=="Michigan",]
data_michigan$newcases = c(0,diff(data_michigan$cases))
data_michigan$newcases[data_michigan$newcases<0] = 0
# plot(data_michigan$newcases)
data_michigan$propB117 =  data_us2$propB117[match(interaction(data_michigan$date, data_michigan$state),
                                                 interaction(data_us2$date,data_us2$REGION))]
data_michigan = data_michigan[!is.na(data_michigan$propB117),]

data_michigan_501YV1 = data.frame(data_michigan, variant="501Y.V1")
data_michigan_501YV1$newcases = data_michigan_501YV1$newcases*data_michigan_501YV1$propB117
data_michigan_wildtype = data.frame(data_michigan, variant="wild type")
data_michigan_wildtype$newcases = data_michigan_wildtype$newcases*(1-data_michigan_wildtype$propB117)
data_michigan = rbind(data_michigan_501YV1, data_michigan_wildtype)
data_michigan$variant = factor(data_michigan$variant, levels=c("wild type", "501Y.V1"), labels=c("wild type","501Y.V1 (British)"))

d = as.Date(max(data_michigan$date))
tag = paste("@TWenseleers\n                           data Helix & NYT\n",d)

plotcases501YV1_Sdropout_Michigan = qplot(data=data_michigan[data_michigan$variant!="total",], x=date, y=newcases, group=variant, colour=variant, fill=variant, geom="blank") +
  geom_area(position="stack") +
  scale_colour_manual("", values=c("grey75","red2")) +
  scale_fill_manual("", values=c("grey75","red2")) +
  ylab("New confirmed cases") + xlab("") + ggtitle("Estimated infections by UK SARS-CoV2 variant 501Y.V1\nin Michigan (S dropout data)") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-11-01"),NA), expand=c(0,0))
plotcases501YV1_Sdropout_Michigan
ggsave(file=paste0(".//plots//",dat,"//confirmed_cases_by_501YV1_vs_wildtype_S dropout_Michigan.png"), width=7, height=5)
