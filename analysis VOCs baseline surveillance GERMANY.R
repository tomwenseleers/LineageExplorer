# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN GERMANY ####
# Tom Wenseleers

# Data: baseline surveillance whole genome sequencing results Robert Kock Insitute, 
# https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/VOC_VOI_Tabelle.html;jsessionid=63E82DB77A4BE622FA4E641582DA0E25.internet091?nn=13490888
# https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/VOC_VOI_Tabelle.xlsx?__blob=publicationFile
# also see
# https://github.com/KITmetricslab/covid19-forecast-hub-de/blob/master/data-truth/RKI/variants/variants_of_concern_sample.csv
# and https://github.com/KITmetricslab/covid19-forecast-hub-de/blob/master/data-truth/RKI/variants/variants_of_interest_sample.csv

# Tom Wenseleers, last update 16 JULY 2021

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


plotdir="DE" # (path in //data)
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))
# today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-07-16")
today_num = as.numeric(today)

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOCs B.1.1.7 (alpha), B.1.351 (beta), P.1 (gamma), B.1.617.1 (kappa) & B.1.617.2 (delta)
# IN GERMANY BASED ON BASELINE SURVEILLANCE SEQUENCING & VOC PCR DATA ####
# (baseline surveillance, i.e. randomly sampled, excluding travellers & surge testing)

library(rio)
VOC_baseline = rio::import(file="https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/VOC_VOI_Tabelle.xlsx?__blob=publicationFile", which=1)
VOI_baseline = rio::import(file="https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/VOC_VOI_Tabelle.xlsx?__blob=publicationFile", which=2)
# PS these files should update every Wednesday  
# VOC_baseline = read.csv("https://raw.githubusercontent.com/KITmetricslab/covid19-forecast-hub-de/master/data-truth/RKI/variants/variants_of_concern_sample.csv")
# VOI_baseline = read.csv("https://raw.githubusercontent.com/KITmetricslab/covid19-forecast-hub-de/master/data-truth/RKI/variants/variants_of_interest_sample.csv")
de_baseline = cbind(VOC_baseline, VOI_baseline[,-1])
de_baseline = de_baseline[-which(grepl("-",de_baseline[,1])),]
propVOCandVOIs = rowSums(de_baseline[,grepl("Anteil",colnames(de_baseline))])/100 # total share that all VOCs and VOIs made up out of all samples
nVOCandVOIs = rowSums(de_baseline[,grepl("Anzahl",colnames(de_baseline))]) # total number of all VOC and VOI samples
de_baseline = de_baseline[,-which(grepl("Anteil",colnames(de_baseline)))]
colnames(de_baseline) = gsub("_Anzahl","",colnames(de_baseline) )
de_baseline = de_baseline[,-which(grepl("Gesamt",colnames(de_baseline)))]
colnms = colnames(de_baseline)
colnms[colnms=="KW"] = "week"
colnms[colnms=="B.1.1.7"] = "alpha"
colnms[colnms=="B.1.351"] = "beta"
colnms[colnms=="B.1.617.2"] = "delta"
colnms[colnms=="P.1"] = "gamma"
colnms[colnms=="B.1.617.1"] = "kappa"
colnames(de_baseline) = colnms
de_baseline$total = round(nVOCandVOIs/propVOCandVOIs,0) # total number of baseline surveillance samples sequenced
de_baseline$other = de_baseline$total-(de_baseline$alpha+de_baseline$beta+de_baseline$gamma+de_baseline$delta+de_baseline$kappa)
de_baseline$week = as.numeric(gsub("KW","",de_baseline$week))
de_baseline$week_startdate = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( de_baseline$week - 1 )
de_baseline$collection_date = de_baseline$week_startdate+3.5 # we use the week midpoint

de_baseline$prop_alpha = de_baseline$alpha / de_baseline$total
de_baseline$prop_beta = de_baseline$beta / de_baseline$total
de_baseline$prop_gamma = de_baseline$gamma / de_baseline$total
de_baseline$prop_kappa = de_baseline$kappa / de_baseline$total
de_baseline$prop_delta = de_baseline$delta / de_baseline$total
de_baseline$n_allVOCs = as.numeric(de_baseline$alpha+de_baseline$beta+de_baseline$gamma+de_baseline$kappa+de_baseline$delta)
de_baseline$prop_allVOCs = de_baseline$n_allVOCs / de_baseline$total

head(de_baseline)
range(de_baseline$collection_date) # "2021-01-04" "2021-06-28"


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
nrow(de_baseline_long2) # n=78132

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

ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_germany_raw data.png"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
data_agbyweek2 = data_agbyweek1 # in fits we recode B.1.1.7 as reference strain
data_agbyweek2$variant = relevel(data_agbyweek2$variant, ref="B.1.1.7 (alpha)")
de_seq_mfit0 = nnet::multinom(variant ~ ns(collection_date_num, df=4), weights=count, data=data_agbyweek2, maxit=1000)
BIC(de_seq_mfit0) # df=4 gave best BIC
summary(de_seq_mfit0)

# growth rate advantage per day compared to UK type B.1.1.7
delta_r = data.frame(confint(emtrends(de_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r) = delta_r[,"contrast"]
delta_r = delta_r[,-1]
delta_r
# estimate   asymp.LCL    asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.02618750 -0.05083319 -0.001541808
# P.1 (gamma) - B.1.1.7 (alpha)        0.03685398  0.01775720  0.055950757
# B.1.617.1 (kappa) - B.1.1.7 (alpha) -0.04904879 -0.12737200  0.029274417
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.14700884  0.13926183  0.154755841
# other - B.1.1.7 (alpha)              0.02647736  0.01901636  0.033938353

# pairwise contrasts in growth rate today (no Tukey correction applied)
emtrends(de_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none")$contrasts
# contrast                              estimate      SE df z.ratio p.value
# B.1.351 (beta) - B.1.1.7 (alpha)       -0.0262 0.01257 NA  -2.083 0.0373 
# P.1 (gamma) - B.1.1.7 (alpha)           0.0369 0.00974 NA   3.782 0.0002 
# P.1 (gamma) - B.1.351 (beta)            0.0630 0.01584 NA   3.981 0.0001 
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    -0.0490 0.03996 NA  -1.227 0.2197 
# B.1.617.1 (kappa) - B.1.351 (beta)     -0.0229 0.04186 NA  -0.546 0.5850 
# B.1.617.1 (kappa) - P.1 (gamma)        -0.0859 0.04110 NA  -2.090 0.0366 
# B.1.617.2 (delta) - B.1.1.7 (alpha)     0.1470 0.00395 NA  37.193 <.0001 
# B.1.617.2 (delta) - B.1.351 (beta)      0.1732 0.01313 NA  13.194 <.0001 
# B.1.617.2 (delta) - P.1 (gamma)         0.1102 0.01023 NA  10.763 <.0001 
# B.1.617.2 (delta) - B.1.617.1 (kappa)   0.1961 0.04014 NA   4.884 <.0001 
# other - B.1.1.7 (alpha)                 0.0265 0.00381 NA   6.955 <.0001 
# other - B.1.351 (beta)                  0.0527 0.01302 NA   4.046 0.0001 
# other - P.1 (gamma)                    -0.0104 0.01032 NA  -1.006 0.3146 
# other - B.1.617.1 (kappa)               0.0755 0.04011 NA   1.883 0.0597 
# other - B.1.617.2 (delta)              -0.1205 0.00515 NA -23.400 <.0001 
# 
# Degrees-of-freedom method: user-specified 

# pairwise contrasts in growth rate today with confidence intervals:
confint(emtrends(de_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none"))$contrasts
# contrast                              estimate      SE df asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)       -0.0262 0.01257 NA  -0.05083  -0.00154
# P.1 (gamma) - B.1.1.7 (alpha)           0.0369 0.00974 NA   0.01776   0.05595
# P.1 (gamma) - B.1.351 (beta)            0.0630 0.01584 NA   0.03200   0.09408
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    -0.0490 0.03996 NA  -0.12737   0.02927
# B.1.617.1 (kappa) - B.1.351 (beta)     -0.0229 0.04186 NA  -0.10490   0.05918
# B.1.617.1 (kappa) - P.1 (gamma)        -0.0859 0.04110 NA  -0.16646  -0.00534
# B.1.617.2 (delta) - B.1.1.7 (alpha)     0.1470 0.00395 NA   0.13926   0.15476
# B.1.617.2 (delta) - B.1.351 (beta)      0.1732 0.01313 NA   0.14747   0.19892
# B.1.617.2 (delta) - P.1 (gamma)         0.1102 0.01023 NA   0.09010   0.13021
# B.1.617.2 (delta) - B.1.617.1 (kappa)   0.1961 0.04014 NA   0.11738   0.27474
# other - B.1.1.7 (alpha)                 0.0265 0.00381 NA   0.01902   0.03394
# other - B.1.351 (beta)                  0.0527 0.01302 NA   0.02715   0.07818
# other - P.1 (gamma)                    -0.0104 0.01032 NA  -0.03060   0.00985
# other - B.1.617.1 (kappa)               0.0755 0.04011 NA  -0.00309   0.15415
# other - B.1.617.2 (delta)              -0.1205 0.00515 NA  -0.13063  -0.11044
# 
# Degrees-of-freedom method: user-specified 
# Confidence level used: 0.95  


# implied increase in infectiousness (due to combination of increased transmissibility and/or immune escape)
# assuming generation time of 4.7 days (Nishiura et al. 2020)
exp(delta_r*4.7) 
# estimate asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    0.8841918 0.7874810 0.9927797
# P.1 (gamma) - B.1.1.7 (alpha)       1.1891202 1.0870405 1.3007858
# B.1.617.1 (kappa) - B.1.1.7 (alpha) 0.7941132 0.5495539 1.1475047
# B.1.617.2 (delta) - B.1.1.7 (alpha) 1.9955936 1.9242391 2.0695940
# other - B.1.1.7 (alpha)             1.1325181 1.0934927 1.1729363


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

ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_germany_multinomial fit.png"), width=7, height=5)

library(ggpubr)
ggarrange(muller_de_raw+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_de_seq_mfit0+ggtitle("Multinomial spline fit (4 df)"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_germany_raw data plus multinomial fit multipanel.png"), width=8, height=8)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_germany_raw data plus multinomial fit multipanel.pdf"), width=8, height=6)


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

ggsave(file=paste0(".\\plots\\",plotdir,"\\germany_baseline_surveillance_multinomial fits_response scale.png"), width=8, height=6)


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


ggsave(file=paste0(".\\plots\\",plotdir,"\\germany_baseline_surveillance_multinomial fits_link scale.png"), width=8, height=6)


# estimated share of different variants of concern among lab diagnoses today
de_seq_mfit0_preds[as.character(de_seq_mfit0_preds$collection_date)==as.character(today),]
# variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1159   B.1.1.7 (alpha)             18824.5 1.445045e-02 1.956974e-03 NA  1.061485e-02 1.828605e-02      2021-07-16
# 1160    B.1.351 (beta)             18824.5 3.922837e-05 2.519093e-05 NA -1.014496e-05 8.860169e-05      2021-07-16
# 1161       P.1 (gamma)             18824.5 1.155857e-03 4.976082e-04 NA  1.805630e-04 2.131151e-03      2021-07-16
# 1162 B.1.617.1 (kappa)             18824.5 2.114425e-06 4.142213e-06 NA -6.004164e-06 1.023301e-05      2021-07-16
# 1163 B.1.617.2 (delta)             18824.5 9.808596e-01 2.568697e-03 NA  9.758250e-01 9.858941e-01      2021-07-16
# 1164             other             18824.5 3.492773e-03 7.278999e-04 NA  2.066115e-03 4.919430e-03      2021-07-16
  
# estimated share of different variants of concern among new infections today (assuming 1 week between infection & diagnosis)
de_seq_mfit0_preds[as.character(de_seq_mfit0_preds$collection_date)==as.character(today+7),]
# variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1201   B.1.1.7 (alpha)             18831.5 5.226087e-03 8.556809e-04 NA  3.548984e-03 6.903191e-03      2021-07-23
# 1202    B.1.351 (beta)             18831.5 1.181093e-05 8.657433e-06 NA -5.157331e-06 2.877918e-05      2021-07-23
# 1203       P.1 (gamma)             18831.5 5.410506e-04 2.710240e-04 NA  9.853265e-06 1.072248e-03      2021-07-23
# 1204 B.1.617.1 (kappa)             18831.5 5.424708e-07 1.213876e-06 NA -1.836682e-06 2.921624e-06      2021-07-23
# 1205 B.1.617.2 (delta)             18831.5 9.927001e-01 1.195519e-03 NA  9.903569e-01 9.950433e-01      2021-07-23
# 1206             other             18831.5 1.520402e-03 3.712686e-04 NA  7.927290e-04 2.248075e-03      2021-07-23

  

# estimated date that B.1.617.2 would make out >50% of all lab diagnoses: "2021-06-19" ["2021-06-18"-"2021-06-19"] 95% CLs (7 days earlier for infections)
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","prob"] >= 0.5)[1]]
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.LCL"] >= 0.5)[1]]
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.UCL"] >= 0.5)[1]]




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot = as.data.frame(get_national_data(countries = "Germany"))
cases_tot = cases_tot[cases_tot$date>=as.Date("2020-08-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
# cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
# cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
k=25
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F) + 
                  WEEKDAY, # BANKHOLIDAY +
                  # s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=poisson(log), data=cases_tot,
                method = "REML",
                knots = list(DATE_NUM = c(min(cases_tot$DATE_NUM)-14,
                                          seq(min(cases_tot$DATE_NUM)+0.5*diff(range(cases_tot$DATE_NUM))/(k-2), 
                                              max(cases_tot$DATE_NUM)-0.5*diff(range(cases_tot$DATE_NUM))/(k-2), length.out=k-2),
                                          max(cases_tot$DATE_NUM)+14))
) 
BIC(fit_cases)



# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

de_seq_mfit0_preds$totcases = cases_tot$cases_new[match(round(de_seq_mfit0_preds$collection_date_num),cases_tot$DATE_NUM)]
de_seq_mfit0_preds$cases = de_seq_mfit0_preds$totcases * de_seq_mfit0_preds$prob
de_seq_mfit0_preds$cases[de_seq_mfit0_preds$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM, at=list(DATE_NUM=seq(date.from, date.to, by=0.5), BANHOLIDAY="no"), type="response"))
de_seq_mfit0_preds$smoothed_totcases = cases_emmeans$rate[match(de_seq_mfit0_preds$collection_date_num,cases_emmeans$DATE_NUM)]
de_seq_mfit0_preds$smoothed_cases = de_seq_mfit0_preds$smoothed_totcases * de_seq_mfit0_preds$prob
de_seq_mfit0_preds$smoothed_cases[de_seq_mfit0_preds$smoothed_cases<=0.001] = NA

ggplot(data=de_seq_mfit0_preds, 
       aes(x=collection_date, y=cases, group=variant)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),max(cases_tot$date)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN GERMANY\n(case data & multinomial fit to baseline surveillance data RKI)") +
  scale_fill_manual("variant", values=colours_VARIANTS) +
  scale_colour_manual("variant", values=colours_VARIANTS) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=de_seq_mfit0_preds, 
       aes(x=collection_date-7, y=smoothed_cases, group=variant)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),today), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN GERMANY\n(case data & multinomial fit to baseline surveillance data RKI)") +
  scale_fill_manual("variant", values=colours_VARIANTS) +
  scale_colour_manual("variant", values=colours_VARIANTS) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data.png"), width=8, height=6)


# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}


# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to)#,
                                             # BANKHOLIDAY="no"
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
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN GERMANY AT MOMENT OF INFECTION BASED ON NEW CASES\n(data RKI)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from,
                                                 date.to), 
                                             function (d) { 
                                               wt = as.data.frame(emmeans(de_seq_mfit0, ~ variant , at=list(collection_date_num=d), type="response"))$prob   # important: these should sum to 1
                                               # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                               cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                               names(cons) = seq_along(cons)
                                               EMT = emtrends(de_seq_mfit0,  ~ variant , by=c("collection_date_num"),
                                                              var="collection_date_num", mode="latent",
                                                              at=list(collection_date_num=d))
                                               out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                               # sum(out$estimate*wt) # should sum to zero
                                               return(out) } ))
above_avg_r_variants = above_avg_r_variants0
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_VARIANTS), 
                                       labels=levels_VARIANTS)
above_avg_r_variants$variant = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$collection_date_num, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2021-01-04" "2021-07-30"
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
ymax = 2.5
ymin = 1/3
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
      x=collection_date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
      group=variant, linetype=I(0)) +
  # facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN GERMANY\nAT MOMENT OF INFECTION\n(based on case data & multinomial fit to\nbaseline surveillance lineage frequencies RKI)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(head(colours_VARIANTS,-1),"black")) +
  scale_colour_manual("variant", values=c(head(colours_VARIANTS,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)

