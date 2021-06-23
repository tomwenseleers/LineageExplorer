# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN NETHERLANDS ####
# Tom Wenseleers

# Data: baseline surveillance whole genome sequencing, https://www.rivm.nl/coronavirus-covid-19/virus/varianten
# downloaded on 23th of June

# Tom Wenseleers, last update 23 JUNE 2021

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


dat="nl" # desired path in //data
suppressWarnings(dir.create(paste0(".//plots//",dat)))
# filedate = as.Date(gsub("_","-",dat)) # file date
# filedate_num = as.numeric(filedate)
# today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-23")
today_num = as.numeric(today)

set_sum_contrasts() # we use effect coding for all models

# 1. ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOCs B.1.1.7 (alpha), B.1.351 (beta), P.1 (gamma), B.1.617.1 (kappa) & B.1.617.2 (delta)
# IN NETHERLANDS BASED ON BASELINE SURVEILLANCE SEQUENCING & VOC PCR DATA ####
# (baseline surveillance, i.e. randomly sampled, excluding travellers & surge testing)

nl_baseline = read.csv(paste0(".\\data\\",dat,"\\baseline_surveillance_netherlands_RIVM_23 june 2021.csv"))
nl_baseline$week_startdate = as.Date(nl_baseline$week_startdate)
nl_baseline$collection_date = nl_baseline$week_startdate+3.5 # we use the week midpoint
nl_baseline$total = nl_baseline$alpha+nl_baseline$beta+nl_baseline$gamma+nl_baseline$kappa+nl_baseline$delta+nl_baseline$other
nl_baseline$prop_alpha = nl_baseline$alpha / nl_baseline$total
nl_baseline$prop_beta = nl_baseline$beta / nl_baseline$total
nl_baseline$prop_gamma = nl_baseline$gamma / nl_baseline$total
nl_baseline$prop_kappa = nl_baseline$kappa / nl_baseline$total
nl_baseline$prop_delta = nl_baseline$delta / nl_baseline$total
nl_baseline$n_allVOCs = nl_baseline$alpha+nl_baseline$beta+nl_baseline$gamma+nl_baseline$kappa+nl_baseline$delta
nl_baseline$prop_allVOCs = nl_baseline$n_allVOCs / nl_baseline$total

head(nl_baseline)
range(nl_baseline$collection_date) # "2020-12-03" "2021-06-03" 


# BASELINE SURVEILLANCE DATA ####

nl_baseline_long = gather(nl_baseline[,c("data",
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
nl_baseline_long$variant = factor(nl_baseline_long$variant, 
                                    levels=c("other","alpha","beta","gamma","kappa","delta","n_allVOCs"), 
                                    labels=c("other", "B.1.1.7 (alpha)", "B.1.351 (beta)", "P.1 (gamma)", "B.1.617.1 (kappa)", "B.1.617.2 (delta)", "all VOCs"))
levels_VARIANTS = c("B.1.1.7 (alpha)","B.1.351 (beta)","P.1 (gamma)","B.1.617.1 (kappa)","B.1.617.2 (delta)","other")
colours_VARIANTS = c("#0085FF","#9A9D00","cyan3",muted("magenta"),"magenta","grey70")

nl_baseline_long$collection_date_num = as.numeric(nl_baseline_long$collection_date)
nl_baseline_long$prop = nl_baseline_long$count / nl_baseline_long$total

# aggregated WGS + VOC PCR data by week
nl_baseline_long2 = nl_baseline_long[,c("collection_date","variant","count")]
nl_baseline_long2 = nl_baseline_long2[nl_baseline_long2$variant!="all VOCs",]
nl_baseline_long2$variant = droplevels(nl_baseline_long2$variant)
nl_baseline_long2 = nl_baseline_long2[rep(seq_len(nrow(nl_baseline_long2)), nl_baseline_long2$count),] # convert to long format
nl_baseline_long2$count = NULL
nrow(nl_baseline_long2) # n=27772

data_agbyweek1 = as.data.frame(table(nl_baseline_long2[,c("collection_date", "variant")]))
colnames(data_agbyweek1) = c("collection_date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ collection_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$collection_date, data_agbyweek1_sum$collection_date)]
sum(data_agbyweek1[data_agbyweek1$variant=="B.1.617.2 (delta)","total"]) == nrow(nl_baseline_long2) # TRUE
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$collection_date))
data_agbyweek1$variant = factor(data_agbyweek1$variant, levels=levels_VARIANTS)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total


# Muller plot raw data
muller_nl_raw = ggplot(data=data_agbyweek1, 
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
  ggtitle("Spread of SARS-CoV2 variants of concern in Netherlands\n(baseline surveillance, whole genome sequencing+VOC PCR)")
muller_nl_raw

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_netherlands_raw data.png"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
data_agbyweek2 = data_agbyweek1 # in fits we recode B.1.1.7 as reference strain
data_agbyweek2$variant = relevel(data_agbyweek2$variant, ref="B.1.1.7 (alpha)")
nl_seq_mfit0 = nnet::multinom(variant ~ ns(collection_date_num, df=2), weights=count, data=data_agbyweek2, maxit=1000)
BIC(nl_seq_mfit0)
summary(nl_seq_mfit0)

# growth rate advantage per day compared to UK type B.1.1.7
delta_r = data.frame(confint(emtrends(nl_seq_mfit0, trt.vs.ctrl ~ variant|1, 
                                                          var="collection_date_num",  mode="latent",
                                                          at=list(collection_date_num=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r) = delta_r[,"contrast"]
delta_r = delta_r[,-1]
delta_r
#                                          estimate    asymp.LCL    asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    -0.0422601122 -0.051416318 -0.033103906
# P.1 (gamma) - B.1.1.7 (alpha)       -0.0002997852 -0.009123517  0.008523946
# B.1.617.1 (kappa) - B.1.1.7 (alpha) -0.2694433281 -0.698286088  0.159399432
# B.1.617.2 (delta) - B.1.1.7 (alpha)  0.0891788549  0.058584113  0.119773597
# other - B.1.1.7 (alpha)             -0.0944917866 -0.104807428 -0.084176145

# pairwise contrasts in growth rate today (no Tukey correction applied)
emtrends(nl_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none")$contrasts
# contrast                              estimate      SE df z.ratio p.value
# B.1.351 (beta) - B.1.1.7 (alpha)       -0.0423 0.00467 NA  -9.046 <.0001 
# P.1 (gamma) - B.1.1.7 (alpha)          -0.0003 0.00450 NA  -0.067 0.9469 
# P.1 (gamma) - B.1.351 (beta)            0.0420 0.00648 NA   6.479 <.0001 
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    -0.2694 0.21880 NA  -1.231 0.2182 
# B.1.617.1 (kappa) - B.1.351 (beta)     -0.2272 0.21885 NA  -1.038 0.2992 
# B.1.617.1 (kappa) - P.1 (gamma)        -0.2691 0.21885 NA  -1.230 0.2188 
# B.1.617.2 (delta) - B.1.1.7 (alpha)     0.0892 0.01561 NA   5.713 <.0001 
# B.1.617.2 (delta) - B.1.351 (beta)      0.1314 0.01630 NA   8.066 <.0001 
# B.1.617.2 (delta) - P.1 (gamma)         0.0895 0.01622 NA   5.516 <.0001 
# B.1.617.2 (delta) - B.1.617.1 (kappa)   0.3586 0.21936 NA   1.635 0.1021 
# other - B.1.1.7 (alpha)                -0.0945 0.00526 NA -17.953 <.0001 
# other - B.1.351 (beta)                 -0.0522 0.00692 NA  -7.547 <.0001 
# other - P.1 (gamma)                    -0.0942 0.00696 NA -13.542 <.0001 
# other - B.1.617.1 (kappa)               0.1750 0.21887 NA   0.799 0.4241 
# other - B.1.617.2 (delta)              -0.1837 0.01648 NA -11.148 <.0001 

# pairwise contrasts in growth rate today with confidence intervals:
confint(emtrends(nl_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA, adjust="none"))$contrasts
# contrast                              estimate      SE df asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)       -0.0423 0.00467 NA  -0.05142  -0.03310
# P.1 (gamma) - B.1.1.7 (alpha)          -0.0003 0.00450 NA  -0.00912   0.00852
# P.1 (gamma) - B.1.351 (beta)            0.0420 0.00648 NA   0.02927   0.05465
# B.1.617.1 (kappa) - B.1.1.7 (alpha)    -0.2694 0.21880 NA  -0.69829   0.15940
# B.1.617.1 (kappa) - B.1.351 (beta)     -0.2272 0.21885 NA  -0.65612   0.20176
# B.1.617.1 (kappa) - P.1 (gamma)        -0.2691 0.21885 NA  -0.69807   0.15979
# B.1.617.2 (delta) - B.1.1.7 (alpha)     0.0892 0.01561 NA   0.05858   0.11977
# B.1.617.2 (delta) - B.1.351 (beta)      0.1314 0.01630 NA   0.09950   0.16338
# B.1.617.2 (delta) - P.1 (gamma)         0.0895 0.01622 NA   0.05768   0.12127
# B.1.617.2 (delta) - B.1.617.1 (kappa)   0.3586 0.21936 NA  -0.07131   0.78856
# other - B.1.1.7 (alpha)                -0.0945 0.00526 NA  -0.10481  -0.08418
# other - B.1.351 (beta)                 -0.0522 0.00692 NA  -0.06580  -0.03867
# other - P.1 (gamma)                    -0.0942 0.00696 NA  -0.10782  -0.08056
# other - B.1.617.1 (kappa)               0.1750 0.21887 NA  -0.25402   0.60392
# other - B.1.617.2 (delta)              -0.1837 0.01648 NA  -0.21596  -0.15138


# implied increase in infectiousness (due to combination of increased transmissibility and/or immune escape)
# assuming generation time of 4.7 days (Nishiura et al. 2020)
# delta has a 54% [45-64%] increased infectiousness compared to alpha
exp(delta_r*4.7) 
#                                      estimate  asymp.LCL asymp.UCL
# B.1.351 (beta) - B.1.1.7 (alpha)    0.8198593 0.78532574 0.8559114
# P.1 (gamma) - B.1.1.7 (alpha)       0.9985920 0.95802584 1.0408759
# B.1.617.1 (kappa) - B.1.1.7 (alpha) 0.2818490 0.03755516 2.1152591
# B.1.617.2 (delta) - B.1.1.7 (alpha) 1.5206542 1.31698539 1.7558199
# other - B.1.1.7 (alpha)             0.6413940 0.61103874 0.6732572


# # PS: mblogit fit would also be possible & would take into account overdispersion
# nl_baseline_long$obs = factor(1:nrow(nl_baseline_long))
# nl_seq_mblogitfit = mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                             # random = ~ 1|obs,
#                             weights = count, data=nl_baseline_long, 
#                             subset=nl_baseline_long$variant!="all VOCs",
#                             dispersion = FALSE)
# dispersion(mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
#                    # random = ~ 1|obs,
#                    weights = count, data=nl_baseline_long,
#                    subset = nl_baseline_long$variant=="wild type"|nl_baseline_long$variant=="all VOCs",
#                    dispersion = TRUE), method="Afroz") # dispersion coefficient = 3.2


# plot multinomial model fit ####

# library(effects)
# plot(Effect("collection_date_num",nl_seq_mfit0), style="stacked")

date.from = min(nl_baseline_long$collection_date_num) 
date.to = as.numeric(as.Date("2021-07-31")) 

nl_seq_mfit0_preds = data.frame(emmeans(nl_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
nl_seq_mfit0_preds$collection_date = as.Date(nl_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
nl_seq_mfit0_preds$variant = factor(nl_seq_mfit0_preds$variant, levels=levels_VARIANTS)

muller_nl_seq_mfit0 = ggplot(data=nl_seq_mfit0_preds, 
                             aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(nl_baseline_long$collection_date)+1, 
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
  ggtitle("Spread of SARS-CoV2 variants of concern in Netherlands\n(baseline surveillance)")
muller_nl_seq_mfit0

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_netherlands_multinomial fit.png"), width=7, height=5)

library(ggpubr)
ggarrange(muller_nl_raw+coord_cartesian(xlim=c(as.Date("2020-12-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_nl_seq_mfit0+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_netherlands_raw data plus multinomial fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\muller plot_netherlands_raw data plus multinomial fit multipanel.pdf"), width=8, height=6)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# on response scale:
plot_multinom_response = qplot(data=nl_seq_mfit0_preds, 
                                     x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Netherlands\n(baseline surveillance, sequencing+VOC PCRs)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-07-31")), expand=c(0,0)) +
  coord_cartesian(xlim=c(min(nl_baseline_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  scale_colour_manual("variant", values=colours_VARIANTS) + # c("red","blue","green3","magenta","black")
  geom_point(data=nl_baseline_long[nl_baseline_long$variant!="all VOCs",],
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

ggsave(file=paste0(".\\plots\\",dat,"\\netherlands_baseline_surveillance_multinomial fits_response scale.png"), width=8, height=6)


# on logit scale:

nl_seq_mfit0_preds2 = nl_seq_mfit0_preds
ymin = 0.001
ymax = 0.990001
nl_seq_mfit0_preds2$asymp.LCL[nl_seq_mfit0_preds2$asymp.LCL<ymin] = ymin
nl_seq_mfit0_preds2$asymp.UCL[nl_seq_mfit0_preds2$asymp.UCL<ymin] = ymin
nl_seq_mfit0_preds2$asymp.UCL[nl_seq_mfit0_preds2$asymp.UCL>ymax] = ymax
nl_seq_mfit0_preds2$prob[nl_seq_mfit0_preds2$prob<ymin] = ymin

plot_multinom = qplot(data=nl_seq_mfit0_preds2[nl_seq_mfit0_preds2$variant!="all VOCs",], x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Netherlands\n(baseline surveillance)") +
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
  geom_point(data=nl_baseline_long[nl_baseline_long$variant!="all VOCs",],
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
  coord_cartesian(xlim=c(min(nl_baseline_long2$collection_date), as.Date("2021-07-31")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom


ggsave(file=paste0(".\\plots\\",dat,"\\netherlands_baseline_surveillance_multinomial fits_link scale.png"), width=8, height=6)


# estimated share of different variants of concern among lab diagnoses today
nl_seq_mfit0_preds[as.character(nl_seq_mfit0_preds$collection_date)==as.character(today),]
# 1213   B.1.1.7 (alpha)             18801.5 8.492118e-01 5.677340e-02 NA  7.379380e-01 9.604856e-01      2021-06-23
# 1214    B.1.351 (beta)             18801.5 6.513709e-04 2.215123e-04 NA  2.172148e-04 1.085527e-03      2021-06-23
# 1215       P.1 (gamma)             18801.5 1.672786e-02 3.888757e-03 NA  9.106034e-03 2.434968e-02      2021-06-23
# 1216 B.1.617.1 (kappa)             18801.5 8.956972e-11 1.102425e-09 NA -2.071143e-09 2.250282e-09      2021-06-23
# 1217 B.1.617.2 (delta)             18801.5 1.334025e-01 5.781059e-02 NA  2.009585e-02 2.467092e-01      2021-06-23
# 1218             other             18801.5 6.438863e-06 2.985200e-06 NA  5.879792e-07 1.228975e-05      2021-06-23
  
# estimated share of different variants of concern among new infections today (assuming 1 week between infection & diagnosis)
nl_seq_mfit0_preds[as.character(nl_seq_mfit0_preds$collection_date)==as.character(today+7),]
# variant collection_date_num         prob           SE df     asymp.LCL    asymp.UCL collection_date
# 1255   B.1.1.7 (alpha)             18808.5 7.613280e-01 1.029830e-01 NA  5.594850e-01 9.631710e-01      2021-06-30
# 1256    B.1.351 (beta)             18808.5 4.344210e-04 1.692666e-04 NA  1.026645e-04 7.661775e-04      2021-06-30
# 1257       P.1 (gamma)             18808.5 1.496528e-02 4.283658e-03 NA  6.569462e-03 2.336109e-02      2021-06-30
# 1258 B.1.617.1 (kappa)             18808.5 1.217847e-11 1.684683e-10 NA -3.180133e-10 3.423703e-10      2021-06-30
# 1259 B.1.617.2 (delta)             18808.5 2.232693e-01 1.049958e-01 NA  1.748135e-02 4.290573e-01      2021-06-30
# 1260             other             18808.5 2.979231e-06 1.530049e-06 NA -1.961114e-08 5.978073e-06      2021-06-30

  

# estimated date that B.1.617.2 would make out >50% of all lab diagnoses: "2021-07-14" ["2021-07-03"->31/7/2021] 95% CLs (7 days earlier for infections)
nl_seq_mfit0_preds[nl_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(nl_seq_mfit0_preds[nl_seq_mfit0_preds$variant=="B.1.617.2 (delta)","prob"] >= 0.5)[1]]
nl_seq_mfit0_preds[nl_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(nl_seq_mfit0_preds[nl_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.LCL"] >= 0.5)[1]]
nl_seq_mfit0_preds[nl_seq_mfit0_preds$variant=="B.1.617.2 (delta)","collection_date"][which(nl_seq_mfit0_preds[nl_seq_mfit0_preds$variant=="B.1.617.2 (delta)","asymp.UCL"] >= 0.5)[1]]


