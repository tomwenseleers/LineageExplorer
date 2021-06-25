# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN BELGIUM ####
# Tom Wenseleers

# Data: baseline surveillance whole genome sequencing results (weeks 49-53 2020, 1-24 2021) + 
# VOC PCR baseline surveillance results (weeks 23-24 2021) (Emmanuel Andre, federal test platform & National Reference Lab) 
# cf weekly Sciensano reports, section 3.4.1,
# https://covid-19.sciensano.be/nl/covid-19-epidemiologische-situatie (federal test platform) & "Genomic surveillance of SARS-CoV-2 in Belgium", 
# reports, https://www.uzleuven.be/nl/laboratoriumgeneeskunde/genomic-surveillance-sars-cov-2-belgium

# Tom Wenseleers, last update 22 JUNE 2021

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


dat="2021_06_22" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
# today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-22")
today_num = as.numeric(today)

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOCs B.1.1.7 (alpha), B.1.351 (beta), P.1 (gamma), B.1.617.1 (kappa) & B.1.617.2 (delta)
# IN BELGIUM BASED ON BASELINE SURVEILLANCE SEQUENCING & VOC PCR DATA ####
# (baseline surveillance, i.e. randomly sampled, excluding travellers & surge testing)

be_baseline = read.csv(paste0(".\\data\\",dat,"\\VOC_baseline_surveillance_belgium_22 june 2021.csv"))
be_baseline$week_startdate = as.Date(be_baseline$week_startdate)
be_baseline$collection_date = be_baseline$week_startdate+3.5 # we use the week midpoint
# except for last VOC PCR datapoint which is from 14 & 15/6 only, so we take the midpoint of that
be_baseline$collection_date[be_baseline$week_startdate==max(be_baseline$week_startdate)&be_baseline$data=="VOC_PCR"] = be_baseline$week_startdate[be_baseline$week_startdate==max(be_baseline$week_startdate)&be_baseline$data=="VOC_PCR"]+0.5 
be_baseline$total = be_baseline$alpha+be_baseline$beta+be_baseline$gamma+be_baseline$kappa+be_baseline$delta+be_baseline$other
be_baseline$prop_alpha = be_baseline$alpha / be_baseline$total
be_baseline$prop_beta = be_baseline$beta / be_baseline$total
be_baseline$prop_gamma = be_baseline$gamma / be_baseline$total
be_baseline$prop_kappa = be_baseline$kappa / be_baseline$total
be_baseline$prop_delta = be_baseline$delta / be_baseline$total
be_baseline$n_allVOCs = be_baseline$alpha+be_baseline$beta+be_baseline$gamma+be_baseline$kappa+be_baseline$delta
be_baseline$prop_allVOCs = be_baseline$n_allVOCs / be_baseline$total

head(be_baseline)
range(be_baseline$collection_date) # "2020-12-03" "2021-06-17" 


# BASELINE SURVEILLANCE DATA ####

be_baseline_long = gather(be_baseline[,c("data",
                                         "collection_date",
                                          "other",
                                          "alpha",
                                          "beta",
                                          "gamma",
                                          "kappa",
                                          "delta",
                                          "n_allVOCs",
                                          "total")], 
                            variant, count, c("other",
                                              "alpha",
                                              "beta",
                                              "gamma",
                                              "kappa",
                                              "delta",
                                              "n_allVOCs"), factor_key=TRUE)
be_baseline_long$variant = factor(be_baseline_long$variant, 
                                    levels=c("other","alpha","beta","gamma","kappa","delta","n_allVOCs"), 
                                    labels=c("other", "B.1.1.7 (alpha)", "B.1.351 (beta)", "P.1 (gamma)", "B.1.617.1 (kappa)", "B.1.617.2 (delta)", "all VOCs"))
levels_VARIANTS = c("B.1.1.7 (alpha)","B.1.351 (beta)","P.1 (gamma)","B.1.617.1 (kappa)","B.1.617.2 (delta)","other")
colours_VARIANTS = c("#0085FF","#9A9D00","cyan3",muted("magenta"),"magenta","grey70")

be_baseline_long$collection_date_num = as.numeric(be_baseline_long$collection_date)
be_baseline_long$prop = be_baseline_long$count / be_baseline_long$total

# aggregated WGS + VOC PCR data by week
be_baseline_long2 = be_baseline_long[,c("collection_date","variant","count")]
be_baseline_long2 = be_baseline_long2[be_baseline_long2$variant!="all VOCs",]
be_baseline_long2$variant = droplevels(be_baseline_long2$variant)
be_baseline_long2 = be_baseline_long2[rep(seq_len(nrow(be_baseline_long2)), be_baseline_long2$count),] # convert to long format
be_baseline_long2$count = NULL
nrow(be_baseline_long2) # n=25829

data_agbyweek1 = as.data.frame(table(be_baseline_long2[,c("collection_date", "variant")]))
colnames(data_agbyweek1) = c("collection_date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ collection_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$collection_date, data_agbyweek1_sum$collection_date)]
sum(data_agbyweek1[data_agbyweek1$variant=="B.1.617.2 (delta)","total"]) == nrow(be_baseline_long2) # TRUE
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$collection_date))
data_agbyweek1$variant = factor(data_agbyweek1$variant, levels=levels_VARIANTS)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total


# Muller plot raw data
muller_be_raw = ggplot(data=data_agbyweek1, 
                               aes(x=collection_date, 
                                   y=count, fill=variant, group=variant)) +
  # facet_wrap(~LABORATORY) +
  geom_area(aes(fill=variant), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("", values=colours_VARIANTS) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01",NA)), expand=c(0,0)) +
  ylab("Share") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance, whole genome sequencing+VOC PCR)")
muller_be_raw

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_belgium_raw data.png"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
data_agbyweek2 = data_agbyweek1 # in fits we recode B.1.1.7 as reference strain
data_agbyweek2$variant = relevel(data_agbyweek2$variant, ref="B.1.1.7 (alpha)")
be_seq_mfit0 = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=data_agbyweek2, maxit=1000)
BIC(be_seq_mfit0)
summary(be_seq_mfit0)

# growth rate advantage per day compared to UK type B.1.1.7
delta_r = data.frame(confint(emtrends(be_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r) = delta_r[,"contrast"]
delta_r = delta_r[,-1]
delta_r
#                                         estimate    asymp.LCL     asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.031874766 -0.038800402 -0.0249491308
# P.1 (gamma) - B.1.1.7 (alpha)        0.005833074  0.001892746  0.0097734019
# B.1.617.1 (kappa) - B.1.1.7 (alpha) -0.011912821 -0.101827009  0.0780013660
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.091767774  0.078593463  0.1049420846
# other - B.1.1.7 (alpha)             -0.003832043 -0.007272370 -0.0003917152

# i.e. delta has a 9.2% [7.9-10.5%] growth rate advantage per day over alpha

# pairwise contrasts in growth rate today (no Tukey correction applied)
emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none")$contrasts
# contrast                              estimate      SE df z.ratio p.value
# B.1.351 (beta) - B.1.1.7 (alpha)      -0.03187 0.00353 NA  -9.021 <.0001 
# P.1 (gamma) - B.1.1.7 (alpha)          0.00583 0.00201 NA   2.901 0.0037 
# P.1 (gamma) - B.1.351 (beta)           0.03771 0.00402 NA   9.382 <.0001 
# B.1.617.1 (kappa) - B.1.1.7 (alpha)   -0.01191 0.04588 NA  -0.260 0.7951 
# B.1.617.1 (kappa) - B.1.351 (beta)     0.01996 0.04601 NA   0.434 0.6644 
# B.1.617.1 (kappa) - P.1 (gamma)       -0.01775 0.04591 NA  -0.387 0.6991 
# B.1.617.2 (delta) - B.1.1.7 (alpha)    0.09177 0.00672 NA  13.652 <.0001 
# B.1.617.2 (delta) - B.1.351 (beta)     0.12364 0.00759 NA  16.288 <.0001 
# B.1.617.2 (delta) - P.1 (gamma)        0.08593 0.00695 NA  12.363 <.0001 
# B.1.617.2 (delta) - B.1.617.1 (kappa)  0.10368 0.04635 NA   2.237 0.0253 
# other - B.1.1.7 (alpha)               -0.00383 0.00176 NA  -2.183 0.0290 
# other - B.1.351 (beta)                 0.02804 0.00378 NA   7.428 <.0001 
# other - P.1 (gamma)                   -0.00967 0.00260 NA  -3.721 0.0002 
# other - B.1.617.1 (kappa)              0.00808 0.04591 NA   0.176 0.8603 
# other - B.1.617.2 (delta)             -0.09560 0.00693 NA -13.792 <.0001 

# pairwise contrasts in growth rate today with confidence intervals:
confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none"))$contrasts
# contrast                              estimate      SE df asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)      -0.03187 0.00353 NA  -0.03880 -0.024949
# P.1 (gamma) - B.1.1.7 (alpha)          0.00583 0.00201 NA   0.00189  0.009773
# P.1 (gamma) - B.1.351 (beta)           0.03771 0.00402 NA   0.02983  0.045585
# B.1.617.1 (kappa) - B.1.1.7 (alpha)   -0.01191 0.04588 NA  -0.10183  0.078001
# B.1.617.1 (kappa) - B.1.351 (beta)     0.01996 0.04601 NA  -0.07022  0.110141
# B.1.617.1 (kappa) - P.1 (gamma)       -0.01775 0.04591 NA  -0.10773  0.072238
# B.1.617.2 (delta) - B.1.1.7 (alpha)    0.09177 0.00672 NA   0.07859  0.104942
# B.1.617.2 (delta) - B.1.351 (beta)     0.12364 0.00759 NA   0.10876  0.138521
# B.1.617.2 (delta) - P.1 (gamma)        0.08593 0.00695 NA   0.07231  0.099559
# B.1.617.2 (delta) - B.1.617.1 (kappa)  0.10368 0.04635 NA   0.01284  0.194524
# other - B.1.1.7 (alpha)               -0.00383 0.00176 NA  -0.00727 -0.000392
# other - B.1.351 (beta)                 0.02804 0.00378 NA   0.02064  0.035442
# other - P.1 (gamma)                   -0.00967 0.00260 NA  -0.01476 -0.004575
# other - B.1.617.1 (kappa)              0.00808 0.04591 NA  -0.08189  0.098056
# other - B.1.617.2 (delta)             -0.09560 0.00693 NA  -0.10919 -0.082014


# implied increase in infectiousness (due to combination of increased transmissibility and/or immune escape)
# assuming generation time of 4.7 days (Nishiura et al. 2020)
# delta has a 54% [45-64%] increased infectiousness compared to alpha
exp(delta_r*4.7) 
#                                      estimate asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    0.8608703 0.8332997 0.8893531
# P.1 (gamma) - B.1.1.7 (alpha)       1.0277947 1.0089356 1.0470063
# B.1.617.1 (kappa) - B.1.1.7 (alpha) 0.9455483 0.6196584 1.4428299
# B.1.617.2 (delta) - B.1.1.7 (alpha) 1.5392704 1.4468507 1.6375935
# other - B.1.1.7 (alpha)             0.9821506 0.9663974 0.9981606


# # PS: mblogit fit would also be possible & would take into account overdispersion
# be_baseline_long$obs = factor(1:nrow(be_baseline_long))
# be_seq_mblogitfit = mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                             # random = ~ 1|obs,
#                             weights = count, data=be_baseline_long, 
#                             subset=be_baseline_long$variant!="all VOCs",
#                             dispersion = FALSE)
# dispersion(mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                    # random = ~ 1|obs,
#                    weights = count, data=be_baseline_long,
#                    subset = be_baseline_long$variant=="wild type"|be_baseline_long$variant=="all VOCs",
#                    dispersion = TRUE), method="Afroz") # dispersion coefficient = 3.2


# plot multinomial model fit ####

# library(effects)
# plot(Effect("collection_date_num",be_seq_mfit0), style="stacked")

date.from = min(be_baseline_long$collection_date_num) 
date.to = as.numeric(as.Date("2021-07-31")) 

be_seq_mfit0_preds = data.frame(emmeans(be_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfit0_preds$collection_date = as.Date(be_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
be_seq_mfit0_preds$variant = factor(be_seq_mfit0_preds$variant, levels=levels_VARIANTS)

muller_be_seq_mfit0 = ggplot(data=be_seq_mfit0_preds, 
                             aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(be_baseline_long$collection_date)+1, 
           xmax=as.Date("2021-07-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("", values=colours_VARIANTS) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-07-31")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance, sequencing+VOC PCR)")
muller_be_seq_mfit0

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_belgium_multinomial fit.png"), width=7, height=5)

library(ggpubr)
ggarrange(muller_be_raw+coord_cartesian(xlim=c(as.Date("2020-12-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_be_seq_mfit0+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_belgium_raw data plus multinomial fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_belgium_raw data plus multinomial fit multipanel.pdf"), width=8, height=6)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# on response scale:
plot_multinom_response = qplot(data=be_seq_mfit0_preds, 
                                     x=collection_date, y=100*prob, geom="blank") +
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
  coord_cartesian(xlim=c(min(be_baseline_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  scale_colour_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  geom_point(data=be_baseline_long[be_baseline_long$variant!="all VOCs",],
             aes(x=collection_date, y=100*prop, size=total,
                 colour=variant, shape=data
             ),
             alpha=I(1)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_multinom_response

ggsave(file=paste0(".\\plots\\",dat,"\\belgium_baseline_surveillance_multinomial fits_response scale.png"), width=8, height=6)


# on logit scale:

be_seq_mfit0_preds2 = be_seq_mfit0_preds
ymin = 0.001
ymax = 0.990001
be_seq_mfit0_preds2$asymp.LCL[be_seq_mfit0_preds2$asymp.LCL<ymin] = ymin
be_seq_mfit0_preds2$asymp.UCL[be_seq_mfit0_preds2$asymp.UCL<ymin] = ymin
be_seq_mfit0_preds2$asymp.UCL[be_seq_mfit0_preds2$asymp.UCL>ymax] = ymax
be_seq_mfit0_preds2$prob[be_seq_mfit0_preds2$prob<ymin] = ymin

plot_multinom = qplot(data=be_seq_mfit0_preds2[be_seq_mfit0_preds2$variant!="all VOCs",], x=collection_date, y=prob, geom="blank") +
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
  scale_fill_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  scale_colour_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  geom_point(data=be_baseline_long[be_baseline_long$variant!="all VOCs",],
             aes(x=collection_date, y=prop, size=total,
                 colour=variant, shape=data
             ),
             alpha=I(1)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(10,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date") +
  coord_cartesian(xlim=c(min(be_baseline_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom


ggsave(file=paste0(".\\plots\\",dat,"\\belgium_baseline_surveillance_multinomial fits_link scale.png"), width=8, height=6)


# estimated share of different variants of concern among lab diagnoses today
be_seq_mfit0_preds[as.character(be_seq_mfit0_preds$collection_date)==as.character(today),]
#                variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1207   B.1.1.7 (alpha)             18800.5 0.5865506133 0.0273196594 NA  0.5330050648 0.6400961618      2021-06-22
# 1208    B.1.351 (beta)             18800.5 0.0021248942 0.0004788357 NA  0.0011863935 0.0030633948      2021-06-22
# 1209       P.1 (gamma)             18800.5 0.0785017693 0.0071403801 NA  0.0645068815 0.0924966570      2021-06-22
# 1210 B.1.617.1 (kappa)             18800.5 0.0002101933 0.0003539612 NA -0.0004835579 0.0009039444      2021-06-22
# 1211 B.1.617.2 (delta)             18800.5 0.3100620704 0.0312537806 NA  0.2488057859 0.3713183548      2021-06-22
# 1212             other             18800.5 0.0225504596 0.0025506211 NA  0.0175513341 0.0275495852      2021-06-22
  
# estimated share of different variants of concern among new infections today (assuming 1 week between infection & diagnosis)
be_seq_mfit0_preds[as.character(be_seq_mfit0_preds$collection_date)==as.character(today+7),]
#                variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1249   B.1.1.7 (alpha)             18807.5 0.4576721277 0.0402486772 NA  0.3787861699 0.5365580855      2021-06-29
# 1250    B.1.351 (beta)             18807.5 0.0013264321 0.0003445511 NA  0.0006511244 0.0020017399      2021-06-29
# 1251       P.1 (gamma)             18807.5 0.0638059705 0.0080048857 NA  0.0481166827 0.0794952582      2021-06-29
# 1252 B.1.617.1 (kappa)             18807.5 0.0001508871 0.0003002436 NA -0.0004375795 0.0007393536      2021-06-29
# 1253 B.1.617.2 (delta)             18807.5 0.4599146867 0.0470136723 NA  0.3677695822 0.5520597911      2021-06-29
# 1254             other             18807.5 0.0171298959 0.0024716702 NA  0.0122855114 0.0219742805      2021-06-29

  

# estimated date that B.1.617.2 would make out >50% of all lab diagnoses: "2021-07-01" ["2021-06-28"-"2021-07-06"] 95% CLs (7 days earlier for infections)
be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","prob"] >= 0.5)[1]]
be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.LCL"] >= 0.5)[1]]
be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.UCL"] >= 0.5)[1]]



# CALCULATION OF EFFECTIVE REPRODUCTION NUMBER BASED ON CASE DATA & Re OF VARIANTS THROUGH TIME ####

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

source("scripts/downloadData.R") # download latest data with new confirmed cases per day from Sciensano website, code adapted from https://github.com/JoFAM/covidBE_analysis by Joris Meys
range(cases_tot$DATE) # "2020-03-01" "2021-06-23"
cases_tot = cases_tot[cases_tot$DATE>=as.Date("2020-08-01"),]

# smooth out weekday effects in case nrs using GAM & correct for unequal testing intensity
fit_cases = gam(CASES ~ s(DATE_NUM, bs="cs", k=25, fx=F) + 
                     WEEKDAY + BANKHOLIDAY +
                     s(TESTS_ALL, bs="cs", k=8, fx=F),
                   family=poisson(log), data=cases_tot,
) 
BIC(fit_cases)

# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to),
                                             BANKHOLIDAY="no"
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[2] = "r"
colnames(avg_r_cases)[5] = "r_LOWER"
colnames(avg_r_cases)[6] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN BELGIUM AT MOMENT OF INFECTION BASED ON NEW CASES\n(data Sciensano)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
                                              function (dat) { 
                                                wt = as.data.frame(emmeans(be_seq_mfit0, ~ variant , at=list(collection_date_num=dat), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(be_seq_mfit0,  ~ variant , by=c("collection_date_num"),
                                                               var="collection_date_num", mode="latent",
                                                               at=list(collection_date_num=dat))
                                                out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                # sum(out$estimate*wt) # should sum to zero
                                                return(out) } ))

above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_VARIANTS), 
                                       labels=levels(data_agbyweek2$variant))
above_avg_r_variants$variant = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$collection_date_num, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2020-12-03" "2021-07-30"
above_avg_r_variants$avg_r = avg_r_cases$r[match(above_avg_r_variants$collection_date,
                                                 avg_r_cases$DATE)]  # average growth rate of all lineages calculated from case nrs
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                collection_date_num=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                # REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE,
                variant="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$variant = factor(above_avg_r_variants$variant, levels=c(levels_VARIANTS,"avg"))
above_avg_r_variants$prob = be_seq_mfit0_preds$prob[match(interaction(above_avg_r_variants$collection_date_num,
                                                                                      above_avg_r_variants$variant),
                                                                          interaction(be_seq_mfit0_preds$collection_date_num,
                                                                                      be_seq_mfit0_preds$variant))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 2
ymin = 1/ymax
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$variant %in% c("other"))|above_avg_r_variants2$collection_date>max(cases_tot$DATE)),], 
      x=collection_date, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
      group=variant, linetype=I(0)) +
  # facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN BELGIUM\n(based on case data Sciensano & multinomial fit to\nbaseline surveillance lineage frequencies federal test platform)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$DATE))) +
  scale_fill_manual("variant", values=c(head(colours_VARIANTS,-1),"black")) +
  scale_colour_manual("variant", values=c(head(colours_VARIANTS,-1),"black")) +
  theme(legend.position="right",  
        axis.title.x=element_blank()) 

ggsave(file=paste0(".\\plots\\",dat,"\\belgium_Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\belgium_Re values per variant_avgRe_from_cases_with clipping.pdf"), width=8, height=6)
