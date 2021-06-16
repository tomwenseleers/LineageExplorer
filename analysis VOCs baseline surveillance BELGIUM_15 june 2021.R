# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN BELGIUM ####
# Tom Wenseleers

# Data: baseline surveillance whole genome sequencing results (weeks 49-53 2020, 1-23 2021) + 
# VOC PCR baseline surveillance results (weeks 23-24 2021) (Emmanuel Andre, federal test platform & National Reference Lab) 
# cf weekly Sciensano reports, section 3.4.1,
# https://covid-19.sciensano.be/nl/covid-19-epidemiologische-situatie (federal test platform) & "Genomic surveillance of SARS-CoV-2 in Belgium", 
# reports, https://www.uzleuven.be/nl/laboratoriumgeneeskunde/genomic-surveillance-sars-cov-2-belgium

# Tom Wenseleers, last update 15 JUNE 2021

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


dat="2021_06_15" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
# today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-15")
today_num = as.numeric(today)

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOCs B.1.1.7, B.1.351, P.1 and B.1.617.2 IN BELGIUM BASED ON BASELINE SURVEILLANCE SEQUENCING DATA ####
# (baseline surveillance, i.e. randomly sampled, excluding travellers & surge testing)

be_baseline = read.csv(paste0(".\\data\\",dat,"\\VOC_baseline_surveillance_belgium_15 june 2021.csv"))
be_baseline$week_startdate = as.Date(be_baseline$week_startdate)
be_baseline$collection_date = be_baseline$week_startdate+3.5 # we use the week midpoint
be_baseline$collection_date[be_baseline$week_startdate==max(be_baseline$week_startdate)] = be_baseline$week_startdate[be_baseline$week_startdate==max(be_baseline$week_startdate)]+0.5 # except for last datapoint which is from 14 & 15/6
be_baseline$prop_alpha = be_baseline$alpha / be_baseline$total
be_baseline$prop_beta = be_baseline$beta / be_baseline$total
be_baseline$prop_gamma = be_baseline$gamma / be_baseline$total
be_baseline$prop_delta = be_baseline$delta / be_baseline$total
be_baseline$n_allVOCs = be_baseline$alpha+be_baseline$beta+be_baseline$gamma+be_baseline$delta
be_baseline$prop_allVOCs = be_baseline$n_allVOCs / be_baseline$total

head(be_baseline)
range(be_baseline$collection_date) # "2020-12-03" "2021-06-14" 


# BASELINE SURVEILLANCE DATA ####

be_baseline_long = gather(be_baseline[,c("data",
                                         "collection_date",
                                          "other",
                                          "alpha",
                                          "beta",
                                          "gamma",
                                          "delta",
                                          "n_allVOCs",
                                          "total")], 
                            variant, count, c("other",
                                              "alpha",
                                              "beta",
                                              "gamma",
                                              "delta",
                                              "n_allVOCs"), factor_key=TRUE)
be_baseline_long$variant = factor(be_baseline_long$variant, 
                                    levels=c("other","alpha","beta","gamma","delta","n_allVOCs"), 
                                    labels=c("other", "B.1.1.7 (alpha)", "B.1.351 (beta)", "P.1 (gamma)", "B.1.617.2 (delta)", "all VOCs"))
levels_VARIANTS = c("B.1.1.7 (alpha)","B.1.351 (beta)","P.1 (gamma)","B.1.617.2 (delta)","other")
colours_VARIANTS = c("#0085FF","#9A9D00","cyan3","magenta","grey70")

be_baseline_long$collection_date_num = as.numeric(be_baseline_long$collection_date)
be_baseline_long$prop = be_baseline_long$count / be_baseline_long$total

# aggregated WGS + VOC PCR data by week
be_baseline_long2 = be_baseline_long[,c("collection_date","variant","count")]
be_baseline_long2 = be_baseline_long2[be_baseline_long2$variant!="all VOCs",]
be_baseline_long2$variant = droplevels(be_baseline_long2$variant)
be_baseline_long2 = be_baseline_long2[rep(seq_len(nrow(be_baseline_long2)), be_baseline_long2$count),] # convert to long format
be_baseline_long2$count = NULL
nrow(be_baseline_long2) # 24955

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
delta_r_4VOCs = data.frame(confint(emtrends(be_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_4VOCs) = delta_r_4VOCs[,"contrast"]
delta_r_4VOCs = delta_r_4VOCs[,-1]
delta_r_4VOCs
#                                         estimate   asymp.LCL    asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.035489822 -0.042786441 -0.0281932037
# P.1 (gamma) - B.1.1.7 (alpha)        0.006053728  0.001738389  0.0103690668
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.086126188  0.068457745  0.1037946310
# other - B.1.1.7 (alpha)             -0.004315329 -0.007882082 -0.0007485765

# pairwise contrasts in growth rate today (no Tukey correction applied)
emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none")$contrasts
# contrast                            estimate      SE df z.ratio p.value
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.03549 0.00372 NA -9.533  <.0001 
# P.1 (gamma) - B.1.1.7 (alpha)        0.00605 0.00220 NA  2.750  0.0060 
# P.1 (gamma) - B.1.351 (beta)         0.04154 0.00428 NA  9.718  <.0001 
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.08613 0.00901 NA  9.554  <.0001 
# B.1.617.2 (delta) - B.1.351 (beta)   0.12162 0.00975 NA 12.472  <.0001 
# B.1.617.2 (delta) - P.1 (gamma)      0.08007 0.00922 NA  8.682  <.0001 
# other - B.1.1.7 (alpha)             -0.00432 0.00182 NA -2.371  0.0177 
# other - B.1.351 (beta)               0.03117 0.00396 NA  7.865  <.0001 
# other - P.1 (gamma)                 -0.01037 0.00278 NA -3.731  0.0002 
# other - B.1.617.2 (delta)           -0.09044 0.00918 NA -9.848  <.0001 
# 
# Degrees-of-freedom method: user-specified

# pairwise contrasts in growth rate today with confidence intervals:
confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none"))$contrasts
# contrast                            estimate      SE df asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.03549 0.00372 NA  -0.04279 -0.028193
# P.1 (gamma) - B.1.1.7 (alpha)        0.00605 0.00220 NA   0.00174  0.010369
# P.1 (gamma) - B.1.351 (beta)         0.04154 0.00428 NA   0.03316  0.049922
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.08613 0.00901 NA   0.06846  0.103795
# B.1.617.2 (delta) - B.1.351 (beta)   0.12162 0.00975 NA   0.10250  0.140729
# B.1.617.2 (delta) - P.1 (gamma)      0.08007 0.00922 NA   0.06200  0.098149
# other - B.1.1.7 (alpha)             -0.00432 0.00182 NA  -0.00788 -0.000749
# other - B.1.351 (beta)               0.03117 0.00396 NA   0.02341  0.038943
# other - P.1 (gamma)                 -0.01037 0.00278 NA  -0.01582 -0.004922
# other - B.1.617.2 (delta)           -0.09044 0.00918 NA  -0.10844 -0.072442
# 
# Degrees-of-freedom method: user-specified 
# Confidence level used: 0.95 

# pairwise contrasts in growth rate evaluated on the 21st of January, when B.1.1.7 made up 16% of all diagnosed infections, with confidence intervals:
confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
                 var="collection_date_num",  mode="latent",
                 at=list(collection_date_num=as.numeric(as.Date("2021-01-21"))), 
                 df=NA, adjust="none"))$contrasts
# contrast                            estimate      SE df asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)     -0.0172 0.00392 NA   -0.0249  -0.00955
# P.1 (gamma) - B.1.1.7 (alpha)         0.0262 0.00545 NA    0.0155   0.03693
# P.1 (gamma) - B.1.351 (beta)          0.0435 0.00660 NA    0.0305   0.05641
# B.1.617.2 (delta) - B.1.1.7 (alpha)   0.0641 0.08623 NA   -0.1049   0.23310
# B.1.617.2 (delta) - B.1.351 (beta)    0.0813 0.08632 NA   -0.0878   0.25052
# B.1.617.2 (delta) - P.1 (gamma)       0.0379 0.08636 NA   -0.1314   0.20713
# other - B.1.1.7 (alpha)              -0.0774 0.00201 NA   -0.0814  -0.07351
# other - B.1.351 (beta)               -0.0602 0.00396 NA   -0.0680  -0.05245
# other - P.1 (gamma)                  -0.1037 0.00569 NA   -0.1148  -0.09253
# other - B.1.617.2 (delta)            -0.1415 0.08625 NA   -0.3106   0.02750


# implied increase in infectiousness (due to combination of increased transmissibility and/or immune escape)
# assuming generation time of 4.7 days (Nishiura et al. 2020)
exp(delta_r_4VOCs*4.7) 
#                                      estimate asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    0.8463670 0.8178337 0.8758959
# P.1 (gamma) - B.1.1.7 (alpha)       1.0288612 1.0082039 1.0499417
# B.1.617.2 (delta) - B.1.1.7 (alpha) 1.4989923 1.3795418 1.6287857
# other - B.1.1.7 (alpha)             0.9799222 0.9636320 0.9964879


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
plot_multinom_4VOCs_response = qplot(data=be_seq_mfit0_preds, 
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
plot_multinom_4VOCs_response

ggsave(file=paste0(".\\plots\\",dat,"\\belgium_baseline_surveillance_multinomial fits_response scale.png"), width=8, height=6)


# on logit scale:

be_seq_mfit0_preds2 = be_seq_mfit0_preds
ymin = 0.001
ymax = 0.990001
be_seq_mfit0_preds2$asymp.LCL[be_seq_mfit0_preds2$asymp.LCL<ymin] = ymin
be_seq_mfit0_preds2$asymp.UCL[be_seq_mfit0_preds2$asymp.UCL<ymin] = ymin
be_seq_mfit0_preds2$asymp.UCL[be_seq_mfit0_preds2$asymp.UCL>ymax] = ymax
be_seq_mfit0_preds2$prob[be_seq_mfit0_preds2$prob<ymin] = ymin

plot_multinom_4VOCs = qplot(data=be_seq_mfit0_preds2[be_seq_mfit0_preds2$variant!="all VOCs",], x=collection_date, y=prob, geom="blank") +
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
plot_multinom_4VOCs


ggsave(file=paste0(".\\plots\\",dat,"\\belgium_baseline_surveillance_multinomial fits_link scale.png"), width=8, height=6)


# estimated share of different variants of concern among lab diagnoses today
be_seq_mfit0_preds[as.character(be_seq_mfit0_preds$collection_date)==as.character(today),]
#              variant collection_date_num        prob          SE df    asymp.LCL   asymp.UCL collection_date
# 971   B.1.1.7 (alpha)             18793.5 0.71612256 0.0222310180 NA 0.672550566 0.759694556      2021-06-15
# 972    B.1.351 (beta)             18793.5 0.00258820 0.0005626088 NA 0.001485507 0.003690893      2021-06-15
# 973       P.1 (gamma)             18793.5 0.09236308 0.0077862360 NA 0.077102339 0.107623823      2021-06-15
# 974 B.1.617.2 (delta)             18793.5 0.15979264 0.0245195428 NA 0.111735223 0.207850065      2021-06-15
# 975             other             18793.5 0.02913351 0.0030060993 NA 0.023241668 0.035025361      2021-06-15

# estimated share of different variants of concern among new infections today (assuming 1 week between infection & diagnosis)
be_seq_mfit0_preds[as.character(be_seq_mfit0_preds$collection_date)==as.character(today+7),]
#                variant collection_date_num        prob          SE df    asymp.LCL   asymp.UCL collection_date
# 1006   B.1.1.7 (alpha)             18800.5 0.631072240 0.0397607605 NA 0.5531425811 0.709001898      2021-06-22
# 1007    B.1.351 (beta)             18800.5 0.001779092 0.0004420257 NA 0.0009127372 0.002645446      2021-06-22
# 1008       P.1 (gamma)             18800.5 0.084916842 0.0093772513 NA 0.0665377674 0.103295917      2021-06-22
# 1009 B.1.617.2 (delta)             18800.5 0.257322285 0.0458682023 NA 0.1674222605 0.347222310      2021-06-22
# 1010             other             18800.5 0.024909541 0.0031596918 NA 0.0187166591 0.031102423      2021-06-22

# estimated date that B.1.617.2 would make out >50% of all lab diagnoses: "2021-07-05" ["2021-06-29"-"2021-07-14"] 95% CLs (7 days earlier for infections)
be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","prob"] >= 0.5)[1]]
be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.LCL"] >= 0.5)[1]]
be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(be_seq_mfit0_preds[be_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.UCL"] >= 0.5)[1]]


