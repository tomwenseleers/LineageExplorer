# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN BELGIUM ####
# Tom Wenseleers

# Data: weekly Sciensano report from the 21st of May 2021 with baseline surveillance data on B.1.617.2 manually added from GISAID

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
range(be_seqdata$collection_date) # "2020-12-03" "2021-05-06"


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
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-04-30")), expand=c(0,0)) +
  ylab("Share") +
  xlab("Collection date") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right") +
  ggtitle("Spread of the SARS-CoV2 variants of concern in Belgium\n(baseline surveillance)")
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
# wild type - B.1.1.7 -0.023674900 -0.027655907 -0.01969389
# B.1.351 - B.1.1.7   -0.048702341 -0.056676301 -0.04072838
# P.1 - B.1.1.7        0.004365824 -0.002000211  0.01073186
# B.1.617.2 - B.1.1.7  0.100101646  0.025601231  0.17460206

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA)$contrasts
# contrast            estimate      SE df z.ratio p.value
# wild type - B.1.1.7   -0.02367 0.00203 NA -11.656 <.0001 
# B.1.351 - B.1.1.7     -0.04870 0.00407 NA -11.971 <.0001 
# B.1.351 - wild type   -0.02503 0.00434 NA  -5.765 <.0001 
# P.1 - B.1.1.7          0.00437 0.00325 NA   1.344 0.6635 
# P.1 - wild type        0.02804 0.00374 NA   7.491 <.0001 
# P.1 - B.1.351          0.05307 0.00513 NA  10.343 <.0001 
# B.1.617.2 - B.1.1.7    0.10010 0.03801 NA   2.633 0.0644 
# B.1.617.2 - wild type  0.12378 0.03806 NA   3.252 0.0101 
# B.1.617.2 - B.1.351    0.14880 0.03823 NA   3.893 0.0009 
# B.1.617.2 - P.1        0.09574 0.03813 NA   2.511 0.0882 
# 
# Degrees-of-freedom method: user-specified 
# P value adjustment: tukey method for comparing a family of 5 estimates 

# confidence intervals:
confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA))
# contrast              estimate      SE df asymp.LCL asymp.UCL
# wild type - B.1.1.7   -0.02367 0.00203 NA  -0.02922   -0.0181
# B.1.351 - B.1.1.7     -0.04870 0.00407 NA  -0.05980   -0.0376
# B.1.351 - wild type   -0.02503 0.00434 NA  -0.03687   -0.0132
# P.1 - B.1.1.7          0.00437 0.00325 NA  -0.00449    0.0132
# P.1 - wild type        0.02804 0.00374 NA   0.01783    0.0383
# P.1 - B.1.351          0.05307 0.00513 NA   0.03907    0.0671
# B.1.617.2 - B.1.1.7    0.10010 0.03801 NA  -0.00358    0.2038
# B.1.617.2 - wild type  0.12378 0.03806 NA   0.01995    0.2276
# B.1.617.2 - B.1.351    0.14880 0.03823 NA   0.04453    0.2531
# B.1.617.2 - P.1        0.09574 0.03813 NA  -0.00827    0.1997


# implied transmission advantage (assuming no immune evasion advantage of B.1.351, if there is such an advantage, transm advantage would be less)
exp(delta_r_4VOCs*4.7) 
#                       estimate asymp.LCL asymp.UCL
# wild type - B.1.1.7 0.8946953 0.8781106 0.9115933
# B.1.351 - B.1.1.7   0.7954073 0.7661490 0.8257829
# P.1 - B.1.1.7       1.0207313 0.9906431 1.0517335
# B.1.617.2 - B.1.1.7 1.6007587 1.1278642 2.2719300

# with confidence intervals (in % increase or decrease):
exp(data.frame(confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
                                var="collection_date_num",  mode="latent",
                                at=list(collection_date_num=today_num), 
                                df=NA))$contrasts)[,c(2,5,6)]*4.7)
#   estimate asymp.LCL asymp.UCL
# 1  0.8946953 0.8716977 0.9182997
# 2  0.7954073 0.7549828 0.8379963
# 3  0.8890258 0.8408989 0.9399073
# 4  1.0207313 0.9790993 1.0641336
# 5  1.1408703 1.0874128 1.1969558
# 6  1.2832814 1.2015813 1.3705365
# 7  1.6007587 0.9832958 2.6059590
# 8  1.7891663 1.0982976 2.9146164
# 9  2.0125020 1.2327790 3.2853936
# 10 1.5682469 0.9618705 2.5568912


# for all 4 variants together compared to wild type
# growth rate advantage compared to wild type
delta_r_allVOCs = data.frame(confint(emtrends(be_seq_mfitallVOC, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_allVOCs) = delta_r_allVOCs[,"contrast"]
delta_r_allVOCs = delta_r_allVOCs[,-1]
delta_r_allVOCs
#                        estimate  asymp.LCL  asymp.UCL
# (all VOCs) - wild type 0.05969318 0.05715543 0.06223092

# implied transmission advantage (assuming no immune evasion advantage of B.1.351, if there is such an advantage, transm advantage would be less)
exp(delta_r_allVOCs*4.7) 
#                          estimate asymp.LCL asymp.UCL
# (all VOCs) - wild type   1.323868  1.308172  1.339753




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
date.from = as.numeric(as.Date("2020-11-01")) # min(be_basseqdata_long$collection_date_num)
date.to = as.numeric(as.Date("2021-06-14")) # max(be_basseqdata_long$collection_date_num)+extrapolate

be_seq_mfit0_preds = data.frame(emmeans(be_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfit0_preds$collection_date = as.Date(be_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
be_seq_mfit0_preds$variant = factor(be_seq_mfit0_preds$variant, levels=c("wild type","B.1.1.7","B.1.351","P.1","B.1.617.2"),
                                    labels=c("wild type","B.1.1.7 (UK)","B.1.351 (South Africa)","P.1 (Brazil)","B.1.617.2 (India)"))

be_seq_mfitallVOCs_preds = data.frame(emmeans(be_seq_mfitallVOC, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfitallVOCs_preds$collection_date = as.Date(be_seq_mfitallVOCs_preds$collection_date_num, origin="1970-01-01")
be_seq_mfitallVOCs_preds$variant = factor(be_seq_mfitallVOCs_preds$variant, levels=c("wild type","all VOCs"),
                                          labels=c("wild type","all VOCs"))

be_basseqdata_long2 = be_basseqdata_long[be_basseqdata_long$variant!="wild type",]
be_basseqdata_long2$variant = droplevels(be_basseqdata_long2$variant) 
be_basseqdata_long2$variant = factor(be_basseqdata_long2$variant, levels=c("B.1.1.7","B.1.351","P.1","B.1.617.2","all VOCs"),
                                     labels=c("B.1.1.7 (UK)","B.1.351 (South Africa)","P.1 (Brazil)","B.1.617.2 (India)","all VOCs"))

muller_be_seq_mfit0 = ggplot(data=be_seq_mfit0_preds, 
                             aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(be_basseqdata_long$collection_date)+1, 
           xmax=as.Date("2021-06-14"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("variant", values=c("grey75","red","blue","green3","magenta")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-06-14")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share among newly diagnosed infections") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance)")
muller_be_seq_mfit0

ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_4VOCs_multinomial fit.png"), width=7, height=5)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS
be_seq_mfit0_preds2 = rbind(be_seq_mfit0_preds[be_seq_mfit0_preds$variant!="wild type",],
                            be_seq_mfitallVOCs_preds[be_seq_mfitallVOCs_preds$variant!="wild type",])
be_seq_mfit0_preds2$variant = droplevels(be_seq_mfit0_preds2$variant)


# on response scale:
plot_multinom_4VOCs_response = qplot(data=be_seq_mfit0_preds2, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(be_basseqdata_long2$collection_date), as.Date("2021-06-14")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("red","blue","green3","magenta","black")) +
  scale_colour_manual("variant", values=c("red","blue","green3","magenta","black")) +
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

plot_multinom_4VOCs = qplot(data=be_seq_mfit0_preds3, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of SARS-CoV2 variants of concern in Belgium\n(baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("red","blue","green3","magenta","black")) +
  scale_colour_manual("variant", values=c("red","blue","green3","magenta","black")) +
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-14")), ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom_4VOCs


ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_4VOCs_multinomial fit_model preds.png"), width=8, height=6)


# estimated share of different variants of concern among lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today,]
#                     variant collection_date_num        prob          SE df    asymp.LCL   asymp.UCL collection_date
# 1006           B.1.1.7 (UK)               18768 0.856663378 0.030619336 NA  0.796650582 0.916676174      2021-05-21
# 1008 B.1.351 (South Africa)               18768 0.003290052 0.000701096 NA  0.001915929 0.004664175      2021-05-21
# 1009           P.1 (Brazil)               18768 0.074080279 0.009002785 NA  0.056435144 0.091725413      2021-05-21
# 1010      B.1.617.2 (India)               18768 0.040189945 0.032833985 NA -0.024163484 0.104543374      2021-05-21
# 4041               all VOCs               18768 0.971982204 0.002876337 NA  0.966344687 0.977619721      2021-05-21

# estimated share of different variants of concern among new infections today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7),]
# variant collection_date_num        prob          SE df    asymp.LCL   asymp.UCL collection_date
# 1041           B.1.1.7 (UK)               18775 0.825133249 0.072058831 NA  0.683900535 0.966365962      2021-05-28
# 1043 B.1.351 (South Africa)               18775 0.002253505 0.000570838 NA  0.001134683 0.003372327      2021-05-28
# 1044           P.1 (Brazil)               18775 0.073567986 0.011869240 NA  0.050304703 0.096831268      2021-05-28
# 1045      B.1.617.2 (India)               18775 0.078009312 0.079721752 NA -0.078242451 0.234261076      2021-05-28
# 4181               all VOCs               18775 0.975931498 0.002796955 NA  0.970449567 0.981413430      2021-05-28

