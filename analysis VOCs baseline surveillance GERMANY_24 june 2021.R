# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN GERMANY ####
# Tom Wenseleers

# Data: baseline surveillance whole genome sequencing results (weeks 1-22 2021+week 23 manually added) 
# Robert Kock Insitute, 
# https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/VOC_VOI_Tabelle.html;jsessionid=63E82DB77A4BE622FA4E641582DA0E25.internet091?nn=13490888
# https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/VOC_VOI_Tabelle.xlsx?__blob=publicationFile
# also see
# https://github.com/KITmetricslab/covid19-forecast-hub-de/blob/master/data-truth/RKI/variants/variants_of_concern_sample.csv
# and https://github.com/KITmetricslab/covid19-forecast-hub-de/blob/master/data-truth/RKI/variants/variants_of_interest_sample.csv

# Tom Wenseleers, last update 24 JUNE 2021

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
# today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-24")
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
colnames(de_baseline) = c("week","alpha","beta","delta","gamma","A.23.1","A.27","B.1.1.318","B.1.324.1","B.1.427","B.1.429",
                          "B.1.525","B.1.526","kappa","B.1.620","C.36.3","C.37","P.2","P.3")
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
nrow(de_baseline_long2) # n=74039

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
#                                        estimate    asymp.LCL    asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.03832583 -0.056168678 -0.020482974
# P.1 (gamma) - B.1.1.7 (alpha)       -0.02692698 -0.051958992 -0.001894972
# B.1.617.1 (kappa) - B.1.1.7 (alpha) -0.07877264 -0.124909270 -0.032636016
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.07563355  0.063110535  0.088156566
# other - B.1.1.7 (alpha)              0.01277389  0.006642653  0.018905120

# pairwise contrasts in growth rate today (no Tukey correction applied)
emtrends(de_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none")$contrasts
# contrast                              estimate      SE df z.ratio p.value
# B.1.351 (beta) - B.1.1.7 (alpha)       -0.0383 0.00910 NA -4.210  <.0001 
# P.1 (gamma) - B.1.1.7 (alpha)          -0.0269 0.01277 NA -2.108  0.0350 
# P.1 (gamma) - B.1.351 (beta)            0.0114 0.01563 NA  0.729  0.4659 
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    -0.0788 0.02354 NA -3.346  0.0008 
# B.1.617.1 (kappa) - B.1.351 (beta)     -0.0404 0.02521 NA -1.605  0.1086 
# B.1.617.1 (kappa) - P.1 (gamma)        -0.0518 0.02674 NA -1.939  0.0526 
# B.1.617.2 (delta) - B.1.1.7 (alpha)     0.0756 0.00639 NA 11.837  <.0001 
# B.1.617.2 (delta) - B.1.351 (beta)      0.1140 0.01108 NA 10.282  <.0001 
# B.1.617.2 (delta) - P.1 (gamma)         0.1026 0.01422 NA  7.214  <.0001 
# B.1.617.2 (delta) - B.1.617.1 (kappa)   0.1544 0.02437 NA  6.336  <.0001 
# other - B.1.1.7 (alpha)                 0.0128 0.00313 NA  4.083  <.0001 
# other - B.1.351 (beta)                  0.0511 0.00952 NA  5.368  <.0001 
# other - P.1 (gamma)                     0.0397 0.01309 NA  3.033  0.0024 
# other - B.1.617.1 (kappa)               0.0915 0.02372 NA  3.860  0.0001 
# other - B.1.617.2 (delta)              -0.0629 0.00700 NA -8.975  <.0001 
# 
# Degrees-of-freedom method: user-specified 

# pairwise contrasts in growth rate today with confidence intervals:
confint(emtrends(de_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none"))$contrasts
# contrast                              estimate      SE df asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)       -0.0383 0.00910 NA  -0.05617 -0.020483
# P.1 (gamma) - B.1.1.7 (alpha)          -0.0269 0.01277 NA  -0.05196 -0.001895
# P.1 (gamma) - B.1.351 (beta)            0.0114 0.01563 NA  -0.01924  0.042039
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    -0.0788 0.02354 NA  -0.12491 -0.032636
# B.1.617.1 (kappa) - B.1.351 (beta)     -0.0404 0.02521 NA  -0.08985  0.008954
# B.1.617.1 (kappa) - P.1 (gamma)        -0.0518 0.02674 NA  -0.10426  0.000573
# B.1.617.2 (delta) - B.1.1.7 (alpha)     0.0756 0.00639 NA   0.06311  0.088157
# B.1.617.2 (delta) - B.1.351 (beta)      0.1140 0.01108 NA   0.09224  0.135683
# B.1.617.2 (delta) - P.1 (gamma)         0.1026 0.01422 NA   0.07470  0.130423
# B.1.617.2 (delta) - B.1.617.1 (kappa)   0.1544 0.02437 NA   0.10664  0.202167
# other - B.1.1.7 (alpha)                 0.0128 0.00313 NA   0.00664  0.018905
# other - B.1.351 (beta)                  0.0511 0.00952 NA   0.03244  0.069758
# other - P.1 (gamma)                     0.0397 0.01309 NA   0.01405  0.065356
# other - B.1.617.1 (kappa)               0.0915 0.02372 NA   0.04506  0.138034
# other - B.1.617.2 (delta)              -0.0629 0.00700 NA  -0.07659 -0.049132
# 
# Degrees-of-freedom method: user-specified 
# Confidence level used: 0.95 


# implied increase in infectiousness (due to combination of increased transmissibility and/or immune escape)
# assuming generation time of 4.7 days (Nishiura et al. 2020)
exp(delta_r*4.7) 
#                                      estimate asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    0.8351605 0.7679791 0.9082188
# P.1 (gamma) - B.1.1.7 (alpha)       0.8811241 0.7833253 0.9911332
# B.1.617.1 (kappa) - B.1.1.7 (alpha) 0.6905745 0.5559519 0.8577957
# B.1.617.2 (delta) - B.1.1.7 (alpha) 1.4268621 1.3453033 1.5133653
# other - B.1.1.7 (alpha)             1.0618761 1.0317129 1.0929211


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
          muller_de_seq_mfit0+ggtitle("Multinomial spline fit (4 df)"), ncol=1)

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
# 1027   B.1.1.7 (alpha)             18802.5 0.5819896569 0.0385345768 NA  5.064633e-01 0.6575160395      2021-06-24
# 1028    B.1.351 (beta)             18802.5 0.0015778321 0.0006041294 NA  3.937603e-04 0.0027619040      2021-06-24
# 1029       P.1 (gamma)             18802.5 0.0019162643 0.0009309424 NA  9.165078e-05 0.0037408778      2021-06-24
# 1030 B.1.617.1 (kappa)             18802.5 0.0001239659 0.0001142523 NA -9.996444e-05 0.0003478962      2021-06-24
# 1031 B.1.617.2 (delta)             18802.5 0.3505364869 0.0423267259 NA  2.675776e-01 0.4334953452      2021-06-24
# 1032             other             18802.5 0.0638557939 0.0083718940 NA  4.744718e-02 0.0802644047      2021-06-24
  
# estimated share of different variants of concern among new infections today (assuming 1 week between infection & diagnosis)
de_seq_mfit0_preds[as.character(de_seq_mfit0_preds$collection_date)==as.character(today+7),]
#                variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1069   B.1.1.7 (alpha)             18809.5 4.656347e-01 5.129599e-02 NA  3.650964e-01 0.5661729808      2021-07-01
# 1070    B.1.351 (beta)             18809.5 9.653349e-04 4.374601e-04 NA  1.079289e-04 0.0018227409      2021-07-01
# 1071       P.1 (gamma)             18809.5 1.269772e-03 7.354532e-04 NA -1.716897e-04 0.0027112339      2021-07-01
# 1072 B.1.617.1 (kappa)             18809.5 5.714241e-05 6.213738e-05 NA -6.464462e-05 0.0001789294      2021-07-01
# 1073 B.1.617.2 (delta)             18809.5 4.762050e-01 5.724575e-02 NA  3.640054e-01 0.5884045830      2021-07-01
# 1074             other             18809.5 5.586808e-02 9.611452e-03 NA  3.702998e-02 0.0747061836      2021-07-01

  

# estimated date that B.1.617.2 would make out >50% of all lab diagnoses: "2021-07-03" ["2021-06-28"-"2021-07-10"] 95% CLs (7 days earlier for infections)
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","prob"] >= 0.5)[1]]
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.LCL"] >= 0.5)[1]]
de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(de_seq_mfit0_preds[de_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.UCL"] >= 0.5)[1]]



