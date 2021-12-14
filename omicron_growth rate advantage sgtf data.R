# ANALYSIS OF GROWTH RATE ADVANTAGE OF OMICRON (B.1.1.529) IN South Africa, England, Denmark & Belgium BASED ON SGTF (S gene target failure / S dropout) DATA
# data: 
# Belgium: Federal Test Platform, provided by Emmanuel André
# Denmark: Statens Serum Institut
# England: UK Health Security Agency (traced off graph by Alex Selby)
# South Africa: Lesley Scott & NHLS team (traced off graph by Alex Selby)

# T. Wenseleers
# last update 13 DECEMBER 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)
library(stringr)
library(lubridate)
library(dplyr)  
library(splines)
library(tidyr)
library(tidyselect)
library(effects)
library(MASS)
library(nlme)
library(lme4)

today = as.Date(Sys.time()) 
# today = as.Date("2021-12-10")
today_num = as.numeric(today)
plotdir = "omicron_sgtf"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

extrapolate = 30
date.from = as.Date("2021-09-14")
date.to = today+extrapolate
dateseq = as.numeric(seq(date.from, date.to, by=1))

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

# import GISAID data to check prop of sgtf samples that are Omicron in function of time (GISAID_sgtf_isomicron)
# source(".//parse_GISAID.R") # load GISAID data
GISAID_sgtf_isomicron = readRDS(".//data//GISAID//GISAID_sgtf_isomicron.rds")
GISAID_sgtf_isomicron$prop_omicron = GISAID_sgtf_isomicron$Omicron/(GISAID_sgtf_isomicron$Omicron+GISAID_sgtf_isomicron$Other)
selcountries = c("South Africa", "Scotland", "Belgium", "England")
GISAID_sgtf_isomicron_subs = GISAID_sgtf_isomicron[GISAID_sgtf_isomicron$country %in% selcountries &
                                                     GISAID_sgtf_isomicron$date>as.Date("2021-07-01"),]

# logistic fit of proportion of sgtf sequences (with Spike_H69del/Spike_V70del) that are Omicron
fit_sgtf_isomicron = glm(cbind(Omicron,Other) ~ country*DATE_NUM, family=binomial, data=GISAID_sgtf_isomicron_subs)
# plot(allEffects(fit_sgtf_isomicron))

fit_sgtf_isomicron_preds = as.data.frame(emmeans(fit_sgtf_isomicron, ~ DATE_NUM*country, at=list(DATE_NUM=as.numeric(seq(as.Date("2021-07-01"),
                                                                                                                         date.to, by=1))), type="response"))
fit_sgtf_isomicron_preds$date = as.Date(fit_sgtf_isomicron_preds$DATE_NUM, origin="1970-01-01")
fit_sgtf_isomicron_preds$country = factor(fit_sgtf_isomicron_preds$country, levels=selcountries)

qplot(data=fit_sgtf_isomicron_preds, x=date, y=prob, ymin=asymp.LCL, ymax=asymp.UCL, geom="blank", fill=country, colour=country) +
  geom_ribbon(aes(colour=NULL), alpha=I(0.2)) +
  geom_line() +
  geom_point(data=GISAID_sgtf_isomicron_subs, aes(y=prop_omicron, ymin=NULL, ymax=NULL)) +
  facet_wrap(~ country) +
  theme_hc() +
  ylab("Proportion of SGTF sequences that are Omicron") +
  ggtitle("Proportion of SGTF sequences\n(with Spike_H69del/Spike_V70del)\nthat are Omicron", "(GISAID data, logistic fit)") +
  xaxis

# import variant-specific PCR data or SGTF data (now proxy for B.1.1.529 / Omicron)
varpcr_dk = read.csv(".//data//omicron_sgtf//variantpcr_denmark.csv") # variant-specfc PCR data for DK
# sgtf_dk = read.csv(".//data//omicron_sgtf//sgtf_denmark.csv") 
sgtf_sa = read.csv(".//data//omicron_sgtf//sgtf_south africa.csv") # data traced from graph (data not open at the moment)
sgtf_sa$date = as.Date(sgtf_sa$date)
sgtf_sa = sgtf_sa[sgtf_sa$date>=as.Date("2021-11-05"),] # data limited to >= Nov 1 when from sequencing data nearly all S dropout was Omicron
# sgtf_eng = read.csv(".//data//omicron_sgtf//sgtf_england.csv") # data traced from Twitter graph
sgtf_eng = read.csv(".//data//omicron_sgtf//sgtf_england_official.csv") # official data (but closed at the moment, so not on github for now)
# sgtf_eng = sgtf_eng[sgtf_eng$date>=as.Date("2021-12-01"),] # data limited to >= Dec 1 when from sequencing data nearly all S dropout was Omicron
sgtf_scot = read.csv(".//data//omicron_sgtf//sgtf_scotland.csv") 
sgtf_be = read.csv(".//data//omicron_sgtf//sgtf_belgium.csv") 

sgtf = rbind(sgtf_sa, sgtf_eng, sgtf_scot, sgtf_be)
sgtf$DATE_NUM = as.numeric(sgtf$date)
# Omicron count is sgtf count * prop sgtf samples inferred to be Omicron (from logistic fit on GISAID data)
prop_sgtf_omicron = predict(fit_sgtf_isomicron, newdata=sgtf, type="response")
sgtf$omicron = sgtf$sgtf * prop_sgtf_omicron
sgtf$prop_omicron = sgtf$omicron / sgtf$pos_tests
sgtf$sgtf = NULL
sgtf$prop_sgtf = NULL
sgtf$DATE_NUM = NULL
sgtf = rbind(sgtf, varpcr_dk)
sgtf$date = as.Date(sgtf$date)
sgtf$date_num = as.numeric(sgtf$date)
levels_country = c("South Africa", "England", "Scotland", "Denmark", "Belgium")
sgtf$country = factor(sgtf$country, levels=levels_country)
sgtf$prop <- sgtf$prop_omicron <- sgtf$omicron/sgtf$pos_tests
sgtf$non_omicron = sgtf$pos_tests-sgtf$omicron
sgtf = sgtf[sgtf$date >= as.Date("2021-11-01"),]
sgtf$obs = as.factor(1:nrow(sgtf)) # for observation-level random effect to take into account overdispersion in binomial GLMM
sgtf$random = factor(1) # fake random effect to be able to run glmmPQL
names(sgtf)
head(sgtf)[,c("country","date","omicron","prop_omicron")]


# ANALYSIS OF SGTF DATA USING LOGISTIC REGRESSION ####

# fit binomial GLM to SGTF data

fit_sgtf0 = glm(cbind(omicron, non_omicron) ~ date_num + country, family=binomial(logit), data=sgtf) 
fit_sgtf1 = glm(cbind(omicron, non_omicron) ~ date_num * country, family=binomial(logit), data=sgtf) 
AIC(fit_sgtf0, fit_sgtf1)
#            df      AIC
# fit_sgtf0  6 885.3442
# fit_sgtf1 10 797.6278 # fits best

summary(fit_sgtf1)

# fit_sgtf0 = glmmPQL(cbind(sgtf, non_sgtf) ~ date_num * country, family=binomial(logit), correlation=corAR1(), random=~1|obs, data=sgtf) # taking into account lag-1 autocorrelation in residuals
# fit_sgtf1 = glmmPQL(cbind(sgtf, non_sgtf) ~ date_num * country, family=binomial(logit), correlation=corAR1(), random=~1|obs, data=sgtf) # taking into account lag-1 autocorrelation in residuals
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf)

plot(allEffects(fit_sgtf1, residuals=T))

# growth rate advantage delta_r of Omicron over Delta (difference in growth rate per day)
deltar_sgtf_bycountry = as.data.frame(emtrends(fit_sgtf1, ~ date_num, by="country", var="date_num", at=list(date_num=max(dateseq))))[,c(2,3,6,7)] # growth rate advantage of Omicron over Delta
rownames(deltar_sgtf_bycountry) = deltar_sgtf_bycountry$country
deltar_sgtf_bycountry$country = NULL
colnames(deltar_sgtf_bycountry) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_bycountry
#                delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# South Africa 0.3429570         0.2555922         0.4303217
# England      0.3932711         0.3801003         0.4064418
# Scotland     0.3257174         0.2983505         0.3530843
# Denmark      0.3299872         0.3153420         0.3446324
# Belgium      0.2439724         0.2130805         0.2748644

# mean growth rate advantage of Omicron over Delta across countries (mean difference in growth rate per day)
deltar_sgtf = as.data.frame(emtrends(fit_sgtf1, ~ date_num, var="date_num", at=list(date_num=max(dateseq))))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
colnames(deltar_sgtf) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf
#     delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# 1 0.327181         0.3074591         0.3469029
deltar_sgtf_char = sapply(deltar_sgtf, function (x) sprintf(as.character(round(x,2)), 2) )
deltar_sgtf_char = paste0(deltar_sgtf_char[1], " [", deltar_sgtf_char[2], "-", deltar_sgtf_char[3],"] 95% CLs")

# transmission advantage of Omicron over Delta (using generation time of 4.7 days)
# i.e. how much higher the effective reproduction number Re value of Omicron is than that of Delta
# and how many more people the virus infects over the course of a single generation
# (due to higher immune escape, i.e. more frequent reinfection of people with immunity built up
# as a result of vaccination or prior exposure)
transmadv_sgtf_by_country = exp(deltar_sgtf_bycountry*4.7) 
colnames(transmadv_sgtf_by_country) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf_by_country
#              transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# South Africa  5.012314            3.324381            7.557284
# England       6.349487            5.968357            6.754957
# Scotland      4.622204            4.064323            5.256661
# Denmark       4.715901            4.402214            5.051940
# Belgium       3.147695            2.722303            3.639559

# mean transmission advantage of Omicron over Delta across countries (using generation time of 4.7 days)
# i.e. mean fold difference in the effective reproduction number Re of Omicron compared to that of Delta
transmadv_sgtf = exp(deltar_sgtf*4.7)
colnames(transmadv_sgtf) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf
#   transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# 1  4.65411            4.242098            5.106139
transmadv_sgtf_char = sapply(transmadv_sgtf, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_char = paste0(transmadv_sgtf_char[1], " [", transmadv_sgtf_char[2], "-", transmadv_sgtf_char[3],"] 95% CLs")

# plot model predictions
emmeans_sgtf = as.data.frame(emmeans(fit_sgtf1, ~ date_num+country, at=list(date_num=dateseq, country=levels_country)), type="response")
colnames(emmeans_sgtf) = c("date_num", "country", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
emmeans_sgtf$date = as.Date(emmeans_sgtf$date_num, origin="1970-01-01")

# plot of share of Omicron among confirmed infections (on logit scale)
qplot(data=sgtf, x=date, y=prop, geom="point", colour=country, fill=country, size=I(1.8)) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf, aes(y=prob), alpha=I(0.5), size=1) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  xlim(date.from, date.to) +
  coord_cartesian(xlim=c(as.Date("2021-10-14"),NA), ylim=c(0.001,0.99)) +
  # scale_size_continuous(range=c()) +
  ggtitle("Spread of Omicron variant inferred from\nS gene target failure (SGTF) data (SA, UK, BE)\n& variant PCR data (DK)",
          "(SGTF counts adjusted by proportion that was estimated to be Omicron\nbased on GISAID data, data Lesley Scott & NHLS team, UKHSA,\nPublic Health Scotland, Statens Serum Institut, Federal Test Platform Belgium)") +
  ylab("Share of Omicron among confirmed cases (%)") +
  annotate("text", x = c(as.Date("2021-10-16")), y = c(0.94),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char,
                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + xaxis +
  theme(plot.subtitle=element_text(size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data.png"), width=8, height=6)

# plot of share of Omicron among confirmed infections (on linear scale)
qplot(data=sgtf, x=date, y=100*prop, geom="point", colour=country, fill=country, size=I(1.8)) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob*100, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf, aes(y=prob*100), alpha=I(0.5), size=1) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  xlim(date.from, date.to) +
  coord_cartesian(xlim=c(as.Date("2021-10-14"),NA), ylim=c(0,100)) +
  # scale_size_continuous(range=c()) +
  ggtitle("Spread of Omicron variant inferred from\nS gene target failure (SGTF) data (SA, UK, BE)\n& variant PCR data (DK)",
          "(SGTF counts adjusted by proportion that was estimated to be Omicron\nbased on GISAID data, data Lesley Scott & NHLS team, UKHSA,\nPublic Health Scotland, Statens Serum Institut, Federal Test Platform Belgium)") +
  ylab("Share of Omicron among confirmed cases (%)") +
  annotate("text", x = c(as.Date("2021-10-16")), y = c(80),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char,
                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + xaxis +
  theme(plot.subtitle=element_text(size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data_linear scale.png"), width=8, height=6)



