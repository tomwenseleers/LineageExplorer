# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN BELGIUM ####
# Tom Wenseleers

# Data: weekly Sciensano report from the 21st of May 2021

# Tom Wenseleers, last update 21 MAY 2021

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


dat="2021_05_21" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-21")
today_num = as.numeric(today)
today # "2021-05-21"

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOC 501Y.V1,VOC 501Y.V2&VOC 501Y.V3 IN BELGIUM BASED ON BASELINE SURVEILLANCE SEQUENCING DATA ####
# (baseline surveillance sequencing results, i.e. randomly sampled)
# data from weekly Sciensano report of 8/5/2021

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
range(be_seqdata$collection_date) # "2020-12-03" "2021-05-13"


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
#                        estimate   asymp.LCL    asymp.UCL
# 501Y.V1 - wild type  0.01918685  0.01495988  0.02341383
# 501Y.V2 - wild type -0.03186883 -0.04110012 -0.02263754
# 501Y.V3 - wild type  0.02227665  0.01438616  0.03016713

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA)$contrasts
# contrast            estimate      SE df z.ratio p.value
# 501Y.V1 - wild type  0.01919 0.00216 NA   8.897 <.0001 
# 501Y.V2 - wild type -0.03187 0.00471 NA  -6.766 <.0001 
# 501Y.V2 - 501Y.V1   -0.05106 0.00444 NA -11.507 <.0001 
# 501Y.V3 - wild type  0.02228 0.00403 NA   5.533 <.0001 
# 501Y.V3 - 501Y.V1    0.00309 0.00351 NA   0.880 0.8154 
# 501Y.V3 - 501Y.V2    0.05415 0.00558 NA   9.704 <.0001 
# 
# Degrees-of-freedom method: user-specified 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# confidence intervals:
confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA))
# contrast            estimate      SE df asymp.LCL asymp.UCL
# 501Y.V1 - wild type  0.01919 0.00216 NA   0.01365    0.0247
# 501Y.V2 - wild type -0.03187 0.00471 NA  -0.04397   -0.0198
# 501Y.V2 - 501Y.V1   -0.05106 0.00444 NA  -0.06245   -0.0397
# 501Y.V3 - wild type  0.02228 0.00403 NA   0.01193    0.0326
# 501Y.V3 - 501Y.V1    0.00309 0.00351 NA  -0.00593    0.0121
# 501Y.V3 - 501Y.V2    0.05415 0.00558 NA   0.03981    0.0685
# 
# Degrees-of-freedom method: user-specified 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 4 estimates 


# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_501V1_501YV2_501YV3*4.7) 
#                     estimate asymp.LCL asymp.UCL
# 501Y.V1 - wild type 1.0943693 1.0728422 1.1163283
# 501Y.V2 - wild type 0.8608943 0.8243414 0.8990682
# 501Y.V3 - wild type 1.1103777 1.0699533 1.1523295

# with confidence intervals (in % increase or decrease):
exp(data.frame(confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
                                var="collection_date_num",  mode="latent",
                                at=list(collection_date_num=today_num), 
                                df=NA))$contrasts)[,c(2,5,6)]*4.7)
#   estimate asymp.LCL asymp.UCL
# 1 1.0943693 1.0662392 1.1232415
# 2 0.8608943 0.8133016 0.9112721
# 3 0.7866580 0.7456246 0.8299496
# 4 1.1103777 1.0576935 1.1656862
# 5 1.0146280 0.9724997 1.0585813
# 6 1.2897956 1.2057596 1.3796885


# for all 3 variants together
# growth rate advantage compared to wild type
delta_r_allVOCs = data.frame(confint(emtrends(be_seq_mfitallVOC, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_allVOCs) = delta_r_allVOCs[,"contrast"]
delta_r_allVOCs = delta_r_allVOCs[,-1]
delta_r_allVOCs
#                               estimate  asymp.LCL  asymp.UCL
# (501Y.V1+V2+V3) - wild type 0.05702273 0.05479781 0.05924765

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_allVOCs*4.7) 
#                             estimate asymp.LCL asymp.UCL
# (501Y.V1+V2+V3) - wild type   1.307356  1.293756  1.321099




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
                        range=c(0.01, 6), limits=c(1,10^4), breaks=c(10,100,1000)) +
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
                        range=c(0.01, 6), limits=c(10,10^4), breaks=c(10,100,1000)) +
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
# 501Y.V1 (British)               18768 0.8901471 0.009323046 NA 0.8718743   0.90842      2021-05-21

# estimated proportion of 501Y.V1 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18775 0.8931314 0.01106152 NA 0.8714512 0.9148116      2021-05-28

# estimated proportion of 501Y.V2 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18768 0.0031991 0.0007005142 NA 0.001826118 0.004572083      2021-05-21

# estimated proportion of 501Y.V2 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18775 0.002245272 0.0005595896 NA 0.001148497 0.003342048      2021-05-28

# estimated proportion of 501Y.V3 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18768 0.07471791 0.008925923 NA 0.05722342 0.0922124      2021-05-21

# estimated proportion of 501Y.V3 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)                18775 0.07660753 0.0107839 NA 0.05547147 0.09774359      2021-05-2



# estimated proportion of one of the three VOCs among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18768 0.9666985 0.003415944 NA 0.9600034 0.9733936      2021-05-21

# estimated proportion of one of the three VOCs among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18775 0.9704742 0.00344916 NA  0.963714 0.9772344      2021-05-28




# the time at which new lab diagnoses would be by more than 50%, 75% 90% by one of the three VOCs :
be_seq_mfit0_preds2_subs = be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date>as.Date("2021-01-01")&
                                                 be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
# >50% by 10th of February [10th Feb - 11th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]

# >75% by 2nd of March [28 Feb - 3d March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]

# >90% by 27th of March [25th March - 29 of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]


# the time at which new infections would be by more than 50%, 75% 90% by one of the three VOCs
# (counting 7 days between infection & diagnosis) :
# >50% by 3d of Feb [3d Feb - 4th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]-7

# >75% by 23d of February [21 Feb - 24th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]-7

# >90% by 20th of March [18th March - 22nd of March] 95% CLs
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
                        range=c(0.01, 6), limits=c(0,10^5), breaks=c(10,100,1000)) +
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