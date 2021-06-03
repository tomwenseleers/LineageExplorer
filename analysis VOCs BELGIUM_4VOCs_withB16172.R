# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT SARS-CoV2 VARIANTS OF CONCERN IN BELGIUM ####
# Tom Wenseleers

# Data: weekly Sciensano report from the 28th of May 2021 with baseline surveillance data on B.1.617.2 & data from later week manually added from GISAID (attempt to use only baseline surveillance)

# Tom Wenseleers, last update 3 JUNE 2021

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


dat="2021_06_03" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-03")
today_num = as.numeric(today)
today # "2021-06-03"

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
range(be_seqdata$collection_date) # "2020-12-03" "2021-05-27"


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
# wild type - B.1.1.7 -0.018197206 -0.022659064 -0.013735349
# B.1.351 - B.1.1.7   -0.049185901 -0.058141981 -0.040229821
# P.1 - B.1.1.7        0.002114385 -0.004069589  0.008298359
# B.1.617.2 - B.1.1.7  0.075281080  0.039830001  0.110732160

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA)$contrasts
# contrast            estimate      SE df z.ratio p.value
# wild type - B.1.1.7   -0.01820 0.00228 NA  -7.994 <.0001 
# B.1.351 - B.1.1.7     -0.04919 0.00457 NA -10.764 <.0001 
# B.1.351 - wild type   -0.03099 0.00488 NA  -6.356 <.0001 
# P.1 - B.1.1.7          0.00211 0.00316 NA   0.670 0.9628 
# P.1 - wild type        0.02031 0.00381 NA   5.334 <.0001 
# P.1 - B.1.351          0.05130 0.00549 NA   9.347 <.0001 
# B.1.617.2 - B.1.1.7    0.07528 0.01809 NA   4.162 0.0003 
# B.1.617.2 - wild type  0.09348 0.01823 NA   5.129 <.0001 
# B.1.617.2 - B.1.351    0.12447 0.01866 NA   6.671 <.0001 
# B.1.617.2 - P.1        0.07317 0.01832 NA   3.995 0.0006 
# 
# Degrees-of-freedom method: user-specified 
# P value adjustment: tukey method for comparing a family of 5 estimates 

# confidence intervals:
confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
         var="collection_date_num",  mode="latent",
         at=list(collection_date_num=today_num), 
         df=NA))
# contrast              estimate      SE df asymp.LCL asymp.UCL
# wild type - B.1.1.7   -0.01820 0.00228 NA  -0.02441   -0.0120
# B.1.351 - B.1.1.7     -0.04919 0.00457 NA  -0.06165   -0.0367
# B.1.351 - wild type   -0.03099 0.00488 NA  -0.04429   -0.0177
# P.1 - B.1.1.7          0.00211 0.00316 NA  -0.00649    0.0107
# P.1 - wild type        0.02031 0.00381 NA   0.00992    0.0307
# P.1 - B.1.351          0.05130 0.00549 NA   0.03633    0.0663
# B.1.617.2 - B.1.1.7    0.07528 0.01809 NA   0.02594    0.1246
# B.1.617.2 - wild type  0.09348 0.01823 NA   0.04376    0.1432
# B.1.617.2 - B.1.351    0.12447 0.01866 NA   0.07357    0.1754
# B.1.617.2 - P.1        0.07317 0.01832 NA   0.02320    0.1231


# implied transmission advantage (assuming no immune evasion advantage of B.1.351, if there is such an advantage, transm advantage would be less)
exp(delta_r_4VOCs*4.7) 
#                       estimate asymp.LCL asymp.UCL
# wild type - B.1.1.7 0.9180285 0.8989772 0.9374835
# B.1.351 - B.1.1.7   0.7936016 0.7608894 0.8277202
# P.1 - B.1.1.7       1.0099872 0.9810547 1.0397729
# B.1.617.2 - B.1.1.7 1.4245003 1.2058696 1.6827698

# with confidence intervals (in % increase or decrease):
exp(data.frame(confint(emtrends(be_seq_mfit0, revpairwise ~ variant|1, 
                                var="collection_date_num",  mode="latent",
                                at=list(collection_date_num=today_num), 
                                df=NA))$contrasts)[,c(2,5,6)]*4.7)
#   estimate asymp.LCL asymp.UCL
# 1  0.9180285 0.8916222 0.9452168
# 2  0.7936016 0.7484452 0.8414825
# 3  0.8644629 0.8120818 0.9202228
# 4  1.0099872 0.9699477 1.0516794
# 5  1.1001697 1.0477456 1.1552170
# 6  1.2726627 1.1861879 1.3654417
# 7  1.4245003 1.1296728 1.7962732
# 8  1.5516951 1.2283499 1.9601562
# 9  1.7949817 1.4130845 2.2800895
# 10 1.4104143 1.1152230 1.7837404


# for all 4 variants together compared to wild type
# growth rate advantage compared to wild type
delta_r_allVOCs = data.frame(confint(emtrends(be_seq_mfitallVOC, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_allVOCs) = delta_r_allVOCs[,"contrast"]
delta_r_allVOCs = delta_r_allVOCs[,-1]
delta_r_allVOCs
#                        estimate  asymp.LCL  asymp.UCL
# (all VOCs) - wild type 0.0521646 0.05042665 0.05390256

# implied transmission advantage (assuming no immune evasion advantage of B.1.351, if there is such an advantage, transm advantage would be less)
exp(delta_r_allVOCs*4.7) 
#                          estimate asymp.LCL asymp.UCL
# (all VOCs) - wild type   1.277843  1.267448  1.288324




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
# 2441           B.1.1.7 (UK)               18781 0.831649993 0.0247327379 NA 0.7831747171 0.880125268      2021-06-03
# 2443 B.1.351 (South Africa)               18781 0.001966892 0.0005082554 NA 0.0009707294 0.002963054      2021-06-03
# 2444           P.1 (Brazil)               18781 0.073279632 0.0092191202 NA 0.0552104887 0.091348776      2021-06-03
# 2445      B.1.617.2 (India)               18781 0.072348310 0.0256006520 NA 0.0221719536 0.122524666      2021-06-03
# 9781               all VOCs               18781 0.977274802 0.0028387363 NA 0.9717109814 0.982838623      2021-06-03

# estimated share of different variants of concern among new infections today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7),]
# variant collection_date_num        prob          SE df    asymp.LCL   asymp.UCL collection_date
# 2476           B.1.1.7 (UK)               18788 0.793383680 0.0472500049 NA 0.7007753717 0.885991988      2021-06-10
# 2478 B.1.351 (South Africa)               18788 0.001329826 0.0003910474 NA 0.0005633877 0.002096265      2021-06-10
# 2479           P.1 (Brazil)               18788 0.070950231 0.0108833771 NA 0.0496192042 0.092281259      2021-06-10
# 2480      B.1.617.2 (India)               18788 0.116904214 0.0513524609 NA 0.0162552403 0.217553188      2021-06-10
# 992                all VOCs               18788 0.979853791 0.0028227106 NA 0.9743213803 0.985386202      2021-06-10

