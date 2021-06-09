# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN BELGIUM ####
# Tom Wenseleers

# Data: weekly Sciensano report from the 28th of May 2021 with baseline surveillance data on B.1.617.2 & data from later week manually added from GISAID (attempt to use only baseline surveillance)

# Tom Wenseleers, last update 9 JUNE 2021

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


dat="2021_06_09" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
# today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-09")
today_num = as.numeric(today)

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOCs B.1.1.7, B.1.351, P.1 and B.1.617.2 IN BELGIUM BASED ON BASELINE SURVEILLANCE SEQUENCING DATA ####
# (baseline surveillance sequencing results, i.e. randomly sampled)
# data from weekly Sciensano report of 21/5/2021 + baseline surveillance data on B.1.617.2 manually added from GISAID

be_seqdata = read.csv(paste0(".\\data\\",dat,"\\baseline_surveillance_4VOCs_withB16172.csv"))
be_seqdata$collection_date = as.Date(be_seqdata$collection_date)
be_seqdata$baselinesurv_n_wild_type = be_seqdata$baselinesurv_total_sequenced-be_seqdata$baselinesurv_n_B.1.1.7-be_seqdata$baselinesurv_n_B.1.351-be_seqdata$baselinesurv_n_P.1-be_seqdata$baselinesurv_n_B.1.617.2
be_seqdata$baselinesurv_propB117 = be_seqdata$baselinesurv_n_B.1.1.7 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_propB1351 = be_seqdata$baselinesurv_n_B.1.351 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_propP1 = be_seqdata$baselinesurv_n_P.1 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_propB16172 = be_seqdata$baselinesurv_n_B.1.617.2 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_n_allVOCs = be_seqdata$baselinesurv_n_B.1.1.7+be_seqdata$baselinesurv_n_B.1.351+be_seqdata$baselinesurv_n_P.1+be_seqdata$baselinesurv_n_B.1.617.2
be_seqdata$baselinesurv_propVOCs = be_seqdata$baselinesurv_n_allVOCs / be_seqdata$baselinesurv_total_sequenced

head(be_seqdata)
range(be_seqdata$collection_date) # "2020-12-03" "2021-06-03"


# BASELINE SURVEILLANCE DATA ####

be_basseqdata_long = gather(be_seqdata[,c("collection_date",
                                          "baselinesurv_n_wild_type",
                                          "baselinesurv_n_B.1.1.7",
                                          "baselinesurv_n_B.1.351",
                                          "baselinesurv_n_P.1",
                                          "baselinesurv_n_B.1.617.2",
                                          "baselinesurv_n_allVOCs",
                                          "baselinesurv_total_sequenced")], 
                            variant, count, c("baselinesurv_n_wild_type",
                                              "baselinesurv_n_B.1.1.7",
                                              "baselinesurv_n_B.1.351",
                                              "baselinesurv_n_P.1",
                                              "baselinesurv_n_B.1.617.2",
                                              "baselinesurv_n_allVOCs"), factor_key=TRUE)
be_basseqdata_long$variant = factor(be_basseqdata_long$variant, 
                                    levels=c("baselinesurv_n_wild_type","baselinesurv_n_B.1.1.7","baselinesurv_n_B.1.351","baselinesurv_n_P.1","baselinesurv_n_B.1.617.2","baselinesurv_n_allVOCs"), 
                                    labels=c("wild type", "B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "all VOCs"))
be_basseqdata_long$collection_date_num = as.numeric(be_basseqdata_long$collection_date)
be_basseqdata_long$prop = be_basseqdata_long$count / be_basseqdata_long$baselinesurv_total_sequenced

baseline_surveillance = ggplot(data=be_basseqdata_long[be_basseqdata_long$variant!="all VOCs",], 
                               aes(x=collection_date, 
                                   y=count, fill=variant, group=variant)) +
  # facet_wrap(~LABORATORY) +
  geom_area(aes(fill=variant), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("variant", values=c("grey75","red","blue","green3","magenta"), 
                    labels=c("wild type","B.1.1.7 (UK)","B.1.351 (South Africa)","P.1 (Brazil)","B.1.617.2 (India)")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01",NA)), expand=c(0,0)) +
  ylab("Share") +
  xlab("Collection date") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance, sequencing+VOC PCR)")
baseline_surveillance

ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_4VOCs.png"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
be_basseqdata_long2 = be_basseqdata_long # in fits we recode B.1.1.7 as reference strain
be_basseqdata_long2$variant = relevel(be_basseqdata_long2$variant, ref="B.1.1.7")
be_seq_mfit0 = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=be_basseqdata_long2, 
                              subset=be_basseqdata_long2$variant!="all VOCs", maxit=1000)
be_seq_mfitallVOC = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=be_basseqdata_long2, 
                                   subset=be_basseqdata_long2$variant=="wild type"|be_basseqdata_long2$variant=="all VOCs", maxit=1000) 
BIC(be_seq_mfit0, be_seq_mfitallVOC)
summary(be_seq_mfit0)

# growth rate advantage compared to UK type B.1.1.7
delta_r_4VOCs = data.frame(confint(emtrends(be_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_4VOCs) = delta_r_4VOCs[,"contrast"]
delta_r_4VOCs = delta_r_4VOCs[,-1]
delta_r_4VOCs
#                        estimate   asymp.LCL    asymp.UCL
# wild type - B.1.1.7 -0.013223553 -0.01741276 -0.009034349
# B.1.351 - B.1.1.7   -0.048768687 -0.05759128 -0.039946097
# P.1 - B.1.1.7        0.007791472  0.00249279  0.013090154
# B.1.617.2 - B.1.1.7  0.097419995  0.06834051  0.126499477

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA)$contrasts
# contrast            estimate      SE df z.ratio p.value
# wild type - B.1.1.7   -0.01322 0.00214 NA  -6.187 <.0001 
# B.1.351 - B.1.1.7     -0.04877 0.00450 NA -10.834 <.0001 
# B.1.351 - wild type   -0.03555 0.00477 NA  -7.444 <.0001 
# P.1 - B.1.1.7          0.00779 0.00270 NA   2.882 0.0323 
# P.1 - wild type        0.02102 0.00336 NA   6.248 <.0001 
# P.1 - B.1.351          0.05656 0.00520 NA  10.885 <.0001 
# B.1.617.2 - B.1.1.7    0.09742 0.01484 NA   6.566 <.0001 
# B.1.617.2 - wild type  0.11064 0.01498 NA   7.384 <.0001 
# B.1.617.2 - B.1.351    0.14619 0.01551 NA   9.427 <.0001 
# B.1.617.2 - P.1        0.08963 0.01503 NA   5.963 <.0001
# 
# Degrees-of-freedom method: user-specified 
# P value adjustment: tukey method for comparing a family of 5 estimates 

# confidence intervals:
confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA))
# contrast              estimate      SE df asymp.LCL asymp.UCL
# wild type - B.1.1.7   -0.01322 0.00214 NA -0.019054  -0.00739
# B.1.351 - B.1.1.7     -0.04877 0.00450 NA -0.061047  -0.03649
# B.1.351 - wild type   -0.03555 0.00477 NA -0.048570  -0.02252
# P.1 - B.1.1.7          0.00779 0.00270 NA  0.000417   0.01517
# P.1 - wild type        0.02102 0.00336 NA  0.011841   0.03019
# P.1 - B.1.351          0.05656 0.00520 NA  0.042387   0.07073
# B.1.617.2 - B.1.1.7    0.09742 0.01484 NA  0.056949   0.13789
# B.1.617.2 - wild type  0.11064 0.01498 NA  0.069771   0.15152
# B.1.617.2 - B.1.351    0.14619 0.01551 NA  0.103890   0.18849
# B.1.617.2 - P.1        0.08963 0.01503 NA  0.048631   0.13063


# implied transmission advantage (assuming no immune evasion advantage of B.1.351, if there is such an advantage, transm advantage would be less)
exp(delta_r_4VOCs*4.7) 
#                       estimate asymp.LCL asymp.UCL
# wild type - B.1.1.7 0.9397413 0.9214194 0.9584274
# B.1.351 - B.1.1.7   0.7951593 0.7628613 0.8288247
# P.1 - B.1.1.7       1.0372987 1.0117850 1.0634557
# B.1.617.2 - B.1.1.7 1.5807098 1.3787819 1.8122108

# with confidence intervals (in % increase or decrease):
exp(data.frame(confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
                                var="collection_date_num",  mode="latent",
                                at=list(collection_date_num=today_num), 
                                df=NA))$contrasts)[,c(2,5,6)]*4.7)
#   estimate asymp.LCL asymp.UCL
# 1  0.9397413 0.9143397 0.9658486
# 2  0.7951593 0.7505693 0.8423982
# 3  0.8461470 0.7959027 0.8995633
# 4  1.0372987 1.0019620 1.0738816
# 5  1.1038131 1.0572283 1.1524506
# 6  1.3045169 1.2204463 1.3943786
# 7  1.5807098 1.3069013 1.9118838
# 8  1.6820692 1.3880838 2.0383184
# 9  1.9879159 1.6295147 2.4251451
# 10 1.5238714 1.2567952 1.8477028


# # PS: mblogit fit would also be possible & would take into account overdispersion
# be_basseqdata_long$obs = factor(1:nrow(be_basseqdata_long))
# be_seq_mblogitfit = mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                             # random = ~ 1|obs,
#                             weights = count, data=be_basseqdata_long, 
#                             subset=be_basseqdata_long$variant!="all VOCs",
#                             dispersion = FALSE)
# dispersion(mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                    # random = ~ 1|obs,
#                    weights = count, data=be_basseqdata_long,
#                    subset = be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="all VOCs",
#                    dispersion = TRUE), method="Afroz") # dispersion coefficient = 3.2


# plot multinomial model fit

# library(effects)
# plot(Effect("collection_date_num",be_seq_mfit0), style="stacked")

# extrapolate = 30*6
date.from = as.numeric(as.Date("2020-02-01")) # min(be_basseqdata_long$collection_date_num)
date.to = as.numeric(as.Date("2021-07-31")) # max(be_basseqdata_long$collection_date_num)+extrapolate

be_seq_mfit0_preds = data.frame(emmeans(be_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfit0_preds$collection_date = as.Date(be_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
be_seq_mfit0_preds$variant = factor(be_seq_mfit0_preds$variant, levels=c("wild type","B.1.1.7","B.1.351","P.1","B.1.617.2"),
                                    labels=c("wild type","B.1.1.7 (UK)","B.1.351 (South Africa)","P.1 (Brazil)","B.1.617.2 (India)"))
be_seq_mfit0_preds_allVOCs = data.frame(variant="all VOCs", 
                                        collection_date_num=seq(date.from, date.to), 
                                        collection_date=as.Date(seq(date.from, date.to), origin="1970-01-01"),
                                        SE=NA,
                                        df=NA)
be_seq_mfit0_preds_allVOCs$prob = rowSums(reshape(be_seq_mfit0_preds[be_seq_mfit0_preds$variant!="wild type",c("collection_date","variant","prob")], 
                                                  idvar="collection_date", timevar="variant", direction="wide")[,-1])
be_seq_mfit0_preds_allVOCs$asymp.LCL = rowSums(reshape(be_seq_mfit0_preds[be_seq_mfit0_preds$variant!="wild type",c("collection_date","variant","asymp.LCL")], 
                                                  idvar="collection_date", timevar="variant", direction="wide")[,-1])
be_seq_mfit0_preds_allVOCs$asymp.UCL = rowSums(reshape(be_seq_mfit0_preds[be_seq_mfit0_preds$variant!="wild type",c("collection_date","variant","asymp.UCL")], 
                                                  idvar="collection_date", timevar="variant", direction="wide")[,-1])


be_seq_mfit0_preds2 = rbind(be_seq_mfit0_preds[be_seq_mfit0_preds$variant!="wild type",], be_seq_mfit0_preds_allVOCs)
be_seq_mfit0_preds2$variant = droplevels(be_basseqdata_long2$variant) 

be_basseqdata_long2 = be_basseqdata_long[be_basseqdata_long$variant!="wild type",]
be_basseqdata_long2$variant = droplevels(be_basseqdata_long2$variant) 
be_basseqdata_long2$variant = factor(be_basseqdata_long2$variant,
                                     levels=c("B.1.1.7","B.1.351","P.1","B.1.617.2","all VOCs"),
                                     labels=c("B.1.1.7 (UK)","B.1.351 (South Africa)","P.1 (Brazil)","B.1.617.2 (India)","all VOCs"))

muller_be_seq_mfit0 = ggplot(data=be_seq_mfit0_preds, 
                             aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(be_basseqdata_long$collection_date)+1, 
           xmax=as.Date("2021-07-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("variant", values=c("grey75","red","blue","green3","magenta")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-07-31")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share among newly diagnosed infections") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance, sequencing+VOC PCR)")
muller_be_seq_mfit0

ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_4VOCs_multinomial fit.png"), width=7, height=5)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS


# on response scale:
plot_multinom_4VOCs_response = qplot(data=be_seq_mfit0_preds2[be_seq_mfit0_preds2$variant!="all VOCs",], x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance, sequencing+VOC PCRs)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-07-31")), expand=c(0,0)) +
  coord_cartesian(xlim=c(min(be_basseqdata_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("red","blue","green3","magenta","black")) +
  scale_colour_manual("variant", values=c("red","blue","green3","magenta","black")) +
  geom_point(data=be_basseqdata_long2[be_basseqdata_long2$variant!="all VOCs",],
             aes(x=collection_date, y=100*prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_multinom_4VOCs_response

ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_4VOCs_multinomial fit_model preds_response.png"), width=8, height=6)


# on logit scale:

be_seq_mfit0_preds3 = be_seq_mfit0_preds2
ymin = 0.001
ymax = 0.990001
be_seq_mfit0_preds3$asymp.LCL[be_seq_mfit0_preds3$asymp.LCL<ymin] = ymin
be_seq_mfit0_preds3$asymp.UCL[be_seq_mfit0_preds3$asymp.UCL<ymin] = ymin
be_seq_mfit0_preds3$asymp.UCL[be_seq_mfit0_preds3$asymp.UCL>ymax] = ymax
be_seq_mfit0_preds3$prob[be_seq_mfit0_preds3$prob<ymin] = ymin

plot_multinom_4VOCs = qplot(data=be_seq_mfit0_preds3[be_seq_mfit0_preds3$variant!="all VOCs",], x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance, sequencing+VOC PCRs)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-07-31")), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("red","blue","green3","magenta","black")) +
  scale_colour_manual("variant", values=c("red","blue","green3","magenta","black")) +
  geom_point(data=be_basseqdata_long2[be_basseqdata_long2$variant!="all VOCs",],
             aes(x=collection_date, y=prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(10,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date") +
  coord_cartesian(xlim=c(min(be_basseqdata_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom_4VOCs


ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_4VOCs_multinomial fit_model preds.png"), width=8, height=6)


# estimated share of different variants of concern among lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today,]
#                     variant collection_date_num        prob          SE df    asymp.LCL   asymp.UCL collection_date
# 2471           B.1.1.7 (UK)               18787 0.778425705 0.0280347831 NA 0.723478540 0.833372870      2021-06-09
# 2473 B.1.351 (South Africa)               18787 0.001495316 0.0004006385 NA 0.000710079 0.002280553      2021-06-09
# 2474           P.1 (Brazil)               18787 0.093999825 0.0103104489 NA 0.073791716 0.114207933      2021-06-09
# 2475      B.1.617.2 (India)               18787 0.104045189 0.0301796102 NA 0.044894240 0.163196138      2021-06-09
# 9901               all VOCs               18787 0.974868041 0.0113782437 NA 0.952567093 0.997168989      2021-06-09

# estimated share of different variants of concern among new infections today (assuming 1 week between infection & diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7),]
#                      variant collection_date_num        prob          SE df    asymp.LCL   asymp.UCL collection_date
# 2506            B.1.1.7 (UK)               18794 0.7047033262 0.0560044086 NA 0.5949367023 0.814469950      2021-06-16
# 2508  B.1.351 (South Africa)               18794 0.0009621933 0.0002952326 NA 0.0003835481 0.001540839      2021-06-16
# 2509            P.1 (Brazil)               18794 0.0898675184 0.0128603565 NA 0.0646616828 0.115073354      2021-06-16
# 2510       B.1.617.2 (India)               18794 0.1862832980 0.0634814912 NA 0.0618618615 0.310704734      2021-06-16
# 10041               all VOCs               18794 0.9770194857 0.0114410641 NA 0.9545954122 0.999443559      2021-06-16

