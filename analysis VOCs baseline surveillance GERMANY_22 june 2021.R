# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN GERMANY ####
# Tom Wenseleers

# Data: baseline surveillance whole genome sequencing results (weeks 1-22 2021+week 23 manually added) 
# Robert Kock Insitute, https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Berichte-VOC-tab.html
# and https://github.com/KITmetricslab/covid19-forecast-hub-de/blob/master/data-truth/RKI/variants/variants_of_concern_sample.csv
# and https://github.com/KITmetricslab/covid19-forecast-hub-de/blob/master/data-truth/RKI/variants/variants_of_interest_sample.csv

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


dat="DE" # (path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
# filedate = as.Date(gsub("_","-",dat)) # file date
# filedate_num = as.numeric(filedate)
# today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-22")
today_num = as.numeric(today)

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOCs B.1.1.7 (alpha), B.1.351 (beta), P.1 (gamma), B.1.617.1 (kappa) & B.1.617.2 (delta)
# IN GERMANY BASED ON BASELINE SURVEILLANCE SEQUENCING & VOC PCR DATA ####
# (baseline surveillance, i.e. randomly sampled, excluding travellers & surge testing)

VOC_baseline = read.csv("https://raw.githubusercontent.com/KITmetricslab/covid19-forecast-hub-de/master/data-truth/RKI/variants/variants_of_concern_sample.csv")
VOI_baseline = read.csv("https://raw.githubusercontent.com/KITmetricslab/covid19-forecast-hub-de/master/data-truth/RKI/variants/variants_of_interest_sample.csv")
de_baseline = cbind(VOC_baseline, VOI_baseline[,-1])
de_baseline = de_baseline[,-which(grepl("prop",colnames(de_baseline)))]
colnames(de_baseline) = gsub("_count","",colnames(de_baseline) )
colnames(de_baseline) = c("week","alpha","beta","gamma","delta","A.23.1","A.27","B.1.1.318","B.1.427","B.1.429",
                          "B.1.525","B.1.526","kappa","B.1.620","P.2","P.3")
n=2000
de_baseline = rbind(de_baseline, data.frame(week=23, # from latest RKI report, total sample size not mentioned though, set at 2000
                                            alpha=n*74.1/100,
                                            beta=n*0.7/100,
                                            gamma=n*0.7/100,
                                            delta=n*15.1/100,
                                            A.23.1=0,
                                            A.27=0,
                                            B.1.1.318=n*1.7/100,
                                            B.1.427=0,
                                            B.1.429=0,
                                            B.1.525=n*1.0/100,
                                            B.1.526=0,
                                            kappa=n*0.5/100,
                                            B.1.620=0,
                                            P.2=0,
                                            P.3=0))
de_baseline$other = de_baseline$A.23.1+de_baseline$A.27+de_baseline$B.1.1.318+de_baseline$B.1.427+de_baseline$B.1.429+
  de_baseline$B.1.525+de_baseline$B.1.526+de_baseline$B.1.620+de_baseline$P.2+de_baseline$P.3

de_baseline$week_startdate = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( de_baseline$week - 1 )
de_baseline$collection_date = de_baseline$week_startdate+3.5 # we use the week midpoint

de_baseline$total = as.numeric(de_baseline$alpha+de_baseline$beta+de_baseline$gamma+de_baseline$kappa+de_baseline$delta+de_baseline$other)
de_baseline$prop_alpha = de_baseline$alpha / de_baseline$total
de_baseline$prop_beta = de_baseline$beta / de_baseline$total
de_baseline$prop_gamma = de_baseline$gamma / de_baseline$total
de_baseline$prop_kappa = de_baseline$kappa / de_baseline$total
de_baseline$prop_delta = de_baseline$delta / de_baseline$total
de_baseline$n_allVOCs = as.numeric(de_baseline$alpha+de_baseline$beta+de_baseline$gamma+de_baseline$kappa+de_baseline$delta)
de_baseline$prop_allVOCs = de_baseline$n_allVOCs / de_baseline$total

head(de_baseline)
range(de_baseline$collection_date) # "2021-01-04" "2021-06-07"


# BASELINE SURVEILLANCE DATA ####

de_baseline_long = gather(de_baseline[,c("collection_date",
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
de_baseline_long$variant = factor(de_baseline_long$variant, 
                                    levels=c("other","alpha","beta","gamma","kappa","delta","n_allVOCs"), 
                                    labels=c("other", "B.1.1.7 (alpha)", "B.1.351 (beta)", "P.1 (gamma)", "B.1.617.1 (kappa)", "B.1.617.2 (delta)", "all VOCs"))
levels_VARIANTS = c("B.1.1.7 (alpha)","B.1.351 (beta)","P.1 (gamma)","B.1.617.1 (kappa)","B.1.617.2 (delta)","other")
colours_VARIANTS = c("#0085FF","#9A9D00","cyan3",muted("magenta"),"magenta","grey70")

de_baseline_long$collection_date_num = as.numeric(de_baseline_long$collection_date)
de_baseline_long$prop = de_baseline_long$count / de_baseline_long$total

# aggregated WGS + VOC PCR data by week
de_baseline_long2 = de_baseline_long[,c("collection_date","variant","count")]
de_baseline_long2 = de_baseline_long2[de_baseline_long2$variant!="all VOCs",]
de_baseline_long2$variant = droplevels(de_baseline_long2$variant)
de_baseline_long2 = de_baseline_long2[rep(seq_len(nrow(de_baseline_long2)), de_baseline_long2$count),] # convert to long format
de_baseline_long2$count = NULL
nrow(de_baseline_long2) # n=52821

data_agbyweek1 = as.data.frame(table(de_baseline_long2[,c("collection_date", "variant")]))
colnames(data_agbyweek1) = c("collection_date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ collection_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$collection_date, data_agbyweek1_sum$collection_date)]
sum(data_agbyweek1[data_agbyweek1$variant=="B.1.617.2 (delta)","total"]) == nrow(de_baseline_long2) # TRUE
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$collection_date))
data_agbyweek1$variant = factor(data_agbyweek1$variant, levels=levels_VARIANTS)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total


# Muller plot raw data
muller_de_raw = ggplot(data=data_agbyweek1, 
                               aes(x=collection_date, 
                                   y=count, fill=variant, group=variant)) +
  # facet_wrap(~LABORATORY) +
  geom_area(aes(fill=variant), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("", values=colours_VARIANTS) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  ylab("Share") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Germany\n(baseline surveillance, data RKI)")
muller_de_raw

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_germany_raw data.png"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
data_agbyweek2 = data_agbyweek1 # in fits we recode B.1.1.7 as reference strain
data_agbyweek2$variant = relevel(data_agbyweek2$variant, ref="B.1.1.7 (alpha)")
de_seq_mfit0 = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=data_agbyweek2, maxit=1000)
BIC(de_seq_mfit0)
summary(de_seq_mfit0)

# growth rate advantage per day compared to UK type B.1.1.7
delta_r = data.frame(confint(emtrends(de_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r) = delta_r[,"contrast"]
delta_r = delta_r[,-1]
delta_r
#                                        estimate   asymp.LCL    asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.007736407 -0.01321147 -0.0022613388
# P.1 (gamma) - B.1.1.7 (alpha)        0.041441655  0.02647645  0.0564068584
# B.1.617.1 (kappa) - B.1.1.7 (alpha)  0.002846667 -0.01258396  0.0182772902
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.061536995  0.05439065  0.0686833383
# other - B.1.1.7 (alpha)             -0.004765093 -0.00883810 -0.0006920857

# pairwise contrasts in growth rate today (no Tukey correction applied)
emtrends(de_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none")$contrasts
# contrast                              estimate      SE df z.ratio p.value
# B.1.351 (beta) - B.1.1.7 (alpha)      -0.00774 0.00279 NA  -2.769 0.0056 
# P.1 (gamma) - B.1.1.7 (alpha)          0.04144 0.00764 NA   5.428 <.0001 
# P.1 (gamma) - B.1.351 (beta)           0.04918 0.00812 NA   6.058 <.0001 
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    0.00285 0.00787 NA   0.362 0.7177 
# B.1.617.1 (kappa) - B.1.351 (beta)     0.01058 0.00834 NA   1.269 0.2046 
# B.1.617.1 (kappa) - P.1 (gamma)       -0.03859 0.01095 NA  -3.524 0.0004 
# B.1.617.2 (delta) - B.1.1.7 (alpha)    0.06154 0.00365 NA  16.877 <.0001 
# B.1.617.2 (delta) - B.1.351 (beta)     0.06927 0.00458 NA  15.132 <.0001 
# B.1.617.2 (delta) - P.1 (gamma)        0.02010 0.00842 NA   2.387 0.0170 
# B.1.617.2 (delta) - B.1.617.1 (kappa)  0.05869 0.00865 NA   6.786 <.0001 
# other - B.1.1.7 (alpha)               -0.00477 0.00208 NA  -2.293 0.0218 
# other - B.1.351 (beta)                 0.00297 0.00345 NA   0.862 0.3886 
# other - P.1 (gamma)                   -0.04621 0.00790 NA  -5.850 <.0001 
# other - B.1.617.1 (kappa)             -0.00761 0.00813 NA  -0.936 0.3491 
# other - B.1.617.2 (delta)             -0.06630 0.00417 NA -15.889 <.0001 
# 
# Degrees-of-freedom method: user-specified 

# pairwise contrasts in growth rate today with confidence intervals:
confint(emtrends(de_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none"))$contrasts
# contrast                              estimate      SE df asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)      -0.00774 0.00279 NA  -0.01321 -0.002261
# P.1 (gamma) - B.1.1.7 (alpha)          0.04144 0.00764 NA   0.02648  0.056407
# P.1 (gamma) - B.1.351 (beta)           0.04918 0.00812 NA   0.03327  0.065090
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    0.00285 0.00787 NA  -0.01258  0.018277
# B.1.617.1 (kappa) - B.1.351 (beta)     0.01058 0.00834 NA  -0.00577  0.026935
# B.1.617.1 (kappa) - P.1 (gamma)       -0.03859 0.01095 NA  -0.06006 -0.017131
# B.1.617.2 (delta) - B.1.1.7 (alpha)    0.06154 0.00365 NA   0.05439  0.068683
# B.1.617.2 (delta) - B.1.351 (beta)     0.06927 0.00458 NA   0.06030  0.078246
# B.1.617.2 (delta) - P.1 (gamma)        0.02010 0.00842 NA   0.00360  0.036595
# B.1.617.2 (delta) - B.1.617.1 (kappa)  0.05869 0.00865 NA   0.04174  0.075640
# other - B.1.1.7 (alpha)               -0.00477 0.00208 NA  -0.00884 -0.000692
# other - B.1.351 (beta)                 0.00297 0.00345 NA  -0.00378  0.009725
# other - P.1 (gamma)                   -0.04621 0.00790 NA  -0.06169 -0.030725
# other - B.1.617.1 (kappa)             -0.00761 0.00813 NA  -0.02354  0.008320
# other - B.1.617.2 (delta)             -0.06630 0.00417 NA  -0.07448 -0.058123


# implied increase in infectiousness (due to combination of increased transmissibility and/or immune escape)
# assuming generation time of 4.7 days (Nishiura et al. 2020)
exp(delta_r*4.7) 
#                                      estimate asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    0.964292 0.9397946 0.9894280
# P.1 (gamma) - B.1.1.7 (alpha)       1.215039 1.1325133 1.3035773
# B.1.617.1 (kappa) - B.1.1.7 (alpha) 1.013469 0.9425705 1.0897009
# B.1.617.2 (delta) - B.1.1.7 (alpha) 1.335391 1.2912827 1.3810053
# other - B.1.1.7 (alpha)             0.977853 0.9593119 0.9967525


# # PS: mblogit fit would also be possible & would take into account overdispersion
# de_baseline_long$obs = factor(1:nrow(de_baseline_long))
# de_seq_mblogitfit = mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                             # random = ~ 1|obs,
#                             weights = count, data=de_baseline_long, 
#                             subset=de_baseline_long$variant!="all VOCs",
#                             dispersion = FALSE)
# dispersion(mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                    # random = ~ 1|obs,
#                    weights = count, data=de_baseline_long,
#                    subset = de_baseline_long$variant=="wild type"|de_baseline_long$variant=="all VOCs",
#                    dispersion = TRUE), method="Afroz") # dispersion coefficient = 3.2


# plot multinomial model fit ####

# library(effects)
# plot(Effect("collection_date_num",de_seq_mfit0), style="stacked")

date.from = min(de_baseline_long$collection_date_num) 
date.to = as.numeric(as.Date("2021-07-31")) 

de_seq_mfit0_preds = data.frame(emmeans(de_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
de_seq_mfit0_preds$collection_date = as.Date(de_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
de_seq_mfit0_preds$variant = factor(de_seq_mfit0_preds$variant, levels=levels_VARIANTS)

muller_de_seq_mfit0 = ggplot(data=de_seq_mfit0_preds, 
                             aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(de_baseline_long$collection_date)+1, 
           xmax=as.Date("2021-07-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("", values=colours_VARIANTS) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-07-31")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Germany\n(baseline surveillance, data RKI)")
muller_de_seq_mfit0

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_germany_multinomial fit.png"), width=7, height=5)

library(ggpubr)
ggarrange(muller_de_raw+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_de_seq_mfit0+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_germany_raw data plus multinomial fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_germany_raw data plus multinomial fit multipanel.pdf"), width=8, height=6)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# on response scale:
plot_multinom_response = qplot(data=de_seq_mfit0_preds, 
                                     x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Germany\n(baseline surveillance, data RKI)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-07-31")), expand=c(0,0)) +
  coord_cartesian(xlim=c(min(de_baseline_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  scale_colour_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  geom_point(data=de_baseline_long[de_baseline_long$variant!="all VOCs",],
             aes(x=collection_date, y=100*prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_multinom_response

ggsave(file=paste0(".\\plots\\",dat,"\\germany_baseline_surveillance_multinomial fits_response scale.png"), width=8, height=6)


# on logit scale:

de_seq_mfit0_preds2 = de_seq_mfit0_preds
ymin = 0.001
ymax = 0.990001
de_seq_mfit0_preds2$asymp.LCL[de_seq_mfit0_preds2$asymp.LCL<ymin] = ymin
de_seq_mfit0_preds2$asymp.UCL[de_seq_mfit0_preds2$asymp.UCL<ymin] = ymin
de_seq_mfit0_preds2$asymp.UCL[de_seq_mfit0_preds2$asymp.UCL>ymax] = ymax
de_seq_mfit0_preds2$prob[de_seq_mfit0_preds2$prob<ymin] = ymin

plot_multinom = qplot(data=de_seq_mfit0_preds2[de_seq_mfit0_preds2$variant!="all VOCs",], x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Germany\n(baseline surveillance, data RKI)") +
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
  geom_point(data=de_baseline_long[de_baseline_long$variant!="all VOCs",],
             aes(x=collection_date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(10,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date") +
  coord_cartesian(xlim=c(min(de_baseline_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom


ggsave(file=paste0(".\\plots\\",dat,"\\germany_baseline_surveillance_multinomial fits_link scale.png"), width=8, height=6)


# estimated share of different variants of concern among lab diagnoses today
de_seq_mfit0_preds[as.character(de_seq_mfit0_preds$collection_date)==as.character(today),]
#                variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1015   B.1.1.7 (alpha)             18800.5 0.685307026 0.0208590613 NA 0.6444240170 0.726190035      2021-06-22
# 1016    B.1.351 (beta)             18800.5 0.006133378 0.0009862228 NA 0.0042004164 0.008066339      2021-06-22
# 1017       P.1 (gamma)             18800.5 0.009729608 0.0030323823 NA 0.0037862474 0.015672968      2021-06-22
# 1018 B.1.617.1 (kappa)             18800.5 0.002901282 0.0010294536 NA 0.0008835904 0.004918974      2021-06-22
# 1019 B.1.617.2 (delta)             18800.5 0.280155623 0.0216289839 NA 0.2377635932 0.322547652      2021-06-22
# 1020             other             18800.5 0.015773084 0.0017912123 NA 0.0122623720 0.019283795      2021-06-22
  
# estimated share of different variants of concern among new infections today (assuming 1 week between infection & diagnosis)
de_seq_mfit0_preds[as.character(de_seq_mfit0_preds$collection_date)==as.character(today+7),]
#                variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1057   B.1.1.7 (alpha)             18807.5 0.594196070 0.0293973460 NA 0.5365783307 0.651813810      2021-06-29
# 1058    B.1.351 (beta)             18807.5 0.005037617 0.0009235754 NA 0.0032274425 0.006847791      2021-06-29
# 1059       P.1 (gamma)             18807.5 0.011275222 0.0040803865 NA 0.0032778110 0.019272632      2021-06-29
# 1060 B.1.617.1 (kappa)             18807.5 0.002566189 0.0010493915 NA 0.0005094194 0.004622959      2021-06-29
# 1061 B.1.617.2 (delta)             18807.5 0.373697486 0.0307565987 NA 0.3134156602 0.433979312      2021-06-29
# 1062             other             18807.5 0.013227416 0.0017544109 NA 0.0097888342 0.016665998      2021-06-29

  

# estimated date that B.1.617.2 would make out >50% of all lab diagnoses: "2021-07-08" ["2021-07-03"-"2021-07-14"] 95% CLs (7 days earlier for infections)
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","prob"] >= 0.5)[1]]
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.LCL"] >= 0.5)[1]]
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.UCL"] >= 0.5)[1]]



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
range(cases_tot$DATE) # "2020-03-01" "2021-06-21"
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
  ggtitle("Re IN GERMANY AT MOMENT OF INFECTION BASED ON NEW CASES\n(data Sciensano)") +
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
                                                wt = as.data.frame(emmeans(de_seq_mfit0, ~ variant , at=list(collection_date_num=dat), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(de_seq_mfit0,  ~ variant , by=c("collection_date_num"),
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
above_avg_r_variants$prob = de_seq_mfit0_preds$prob[match(interaction(above_avg_r_variants$collection_date_num,
                                                                                      above_avg_r_variants$variant),
                                                                          interaction(de_seq_mfit0_preds$collection_date_num,
                                                                                      de_seq_mfit0_preds$variant))]
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
qplot(data=above_avg_r_variants2[!(above_avg_r_variants2$variant %in% c("other")),], 
      x=collection_date, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
      group=variant, linetype=I(0)) +
  # facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN GERMANY\n(based on case data Sciensano & multinomial fit to\nbaseline surveillance lineage frequencies federal test platform)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$DATE))) +
  scale_fill_manual("variant", values=c(head(colours_VARIANTS,-1),"black")) +
  scale_colour_manual("variant", values=c(head(colours_VARIANTS,-1),"black")) +
  theme(legend.position="right",  
        axis.title.x=element_blank())

ggsave(file=paste0(".\\plots\\",dat,"\\germany_Re values per variant_with clipping.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\germany_Re values per variant_with clipping.pdf"), width=8, height=6)
