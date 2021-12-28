# ANALYSIS OF GROWTH RATE ADVANTAGE OF OMICRON IN South Africa, England, Denmark & Belgium BASED ON SGTF (S gene target failure / S dropout) DATA
# data: 
# Belgium: Federal Test Platform, provided by Emmanuel André
# Denmark: Statens Serum Institut (https://covid19.ssi.dk/virusvarianter/delta-pcr)
# England: UK Health Security Agency (https://www.gov.uk/government/publications/covid-19-omicron-daily-overview)
# Scotland: https://www.gov.scot/publications/coronavirus-covid-19-additional-data-and-information/
# South Africa: Lesley Scott & NHLS team (traced off graph by Alex Selby)
# GISAID data to infer the proportion of S dropout (with spike Spike_H69del/Spike_V70del) that are Omicron (https://www.gisaid.org/)

# T. Wenseleers
# last update 28 DECEMBER 2021

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
library(scam)
library(mgcv)

today = as.Date(Sys.time()) 
# today = as.Date("2021-12-10")
today_num = as.numeric(today)
plotdir = "omicron_sgtf"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

extrapolate = 30
date.from = as.Date("2021-09-01")
date.to = today+extrapolate
dateseq = as.numeric(seq(date.from, date.to, by=1))

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

# 1. Import GISAID data to check prop of sgtf samples that are Omicron in function of time (GISAID_sgtf_isomicron)
# source(".//parse_GISAID.R") # load GISAID data
GISAID_sgtf_isomicron = readRDS(".//data//GISAID//GISAID_sgtf_isomicron.rds")
GISAID_sgtf_isomicron$prop_omicron = GISAID_sgtf_isomicron$Omicron/(GISAID_sgtf_isomicron$Omicron+GISAID_sgtf_isomicron$Other)
GISAID_sgtf_isomicron$total_seqs_sgtf = GISAID_sgtf_isomicron$Omicron+GISAID_sgtf_isomicron$Other
selcountries = c("South Africa", "Scotland", "Belgium", "England")
GISAID_sgtf_isomicron_subs = GISAID_sgtf_isomicron[GISAID_sgtf_isomicron$country %in% selcountries &
                                                     GISAID_sgtf_isomicron$date>as.Date("2021-07-01"),]

# 2. Logistic fit of proportion of sgtf sequences (with Spike_H69del/Spike_V70del) that are Omicron
fit_sgtf_isomicron = glm(cbind(Omicron,Other) ~ country*DATE_NUM, family=binomial, data=GISAID_sgtf_isomicron_subs)
# plot(allEffects(fit_sgtf_isomicron))

fit_sgtf_isomicron_preds = as.data.frame(emmeans(fit_sgtf_isomicron, ~ DATE_NUM*country, at=list(DATE_NUM=as.numeric(seq(as.Date("2021-07-01"),
                                                                                                                         date.to, by=1))), type="response"))
fit_sgtf_isomicron_preds$date = as.Date(fit_sgtf_isomicron_preds$DATE_NUM, origin="1970-01-01")
fit_sgtf_isomicron_preds$country = factor(fit_sgtf_isomicron_preds$country, levels=selcountries)

qplot(data=fit_sgtf_isomicron_preds, x=date, y=prob, ymin=asymp.LCL, ymax=asymp.UCL, geom="blank", fill=country, colour=country) +
  geom_ribbon(aes(colour=NULL), alpha=I(0.4)) +
  geom_line() +
  geom_point(data=GISAID_sgtf_isomicron_subs, aes(y=prop_omicron, ymin=NULL, ymax=NULL, size=total_seqs_sgtf)) +
  scale_size_continuous("nr. of GISAID sequences\nwith SGTF", range=c(1,4), breaks=c(10, 100, 1000, 10000)) +
  # facet_wrap(~ country) +
  theme_hc() +
  ylab("Proportion of SGTF sequences that are Omicron") +
  ggtitle("Proportion of GISAID sequences with SGTF\n(with Spike_H69del/Spike_V70del)\nthat are Omicron", "(GISAID data, logistic fit)") +
  xaxis + theme(legend.position="right") +
  coord_cartesian(xlim=c(as.Date("2021-09-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\estimated prop sgtf seqs that are Omicron from GISAID data.png"), width=8, height=6)

# 3. Import variant-specific PCR data or SGTF data (now proxy for Omicron)
varpcr_dk = read.csv(".//data//omicron_sgtf//variantpcr_denmark.csv") # variant-specfc PCR data for DK from https://covid19.ssi.dk/virusvarianter/delta-pcr
sgtf_sa = read.csv(".//data//omicron_sgtf//sgtf_south africa.csv") # data traced from graph (data not open at the moment)
sgtf_sa$date = as.Date(sgtf_sa$date)
sgtf_sa = sgtf_sa[sgtf_sa$date>=as.Date("2021-11-05"),] # data limited to >= Nov 5 when from sequencing data nearly all S dropout was Omicron

# get latest SGTF data for England from https://www.gov.uk/government/publications/covid-19-omicron-daily-overview
html = paste(readLines("https://www.gov.uk/government/publications/covid-19-omicron-daily-overview"), collapse="\n")
library(stringr)
matched <- str_match_all(html, "<a class=\"govuk-link\" href=\"(.*?)\"")
url_pattern <- "http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\\(\\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+"
links = str_match_all(html, url_pattern)[[1]]
link = links[grepl("sgtf_regionepicurve", links)]
link_byLTLA = links[grepl("sgtf_totalepicurve", links)]
sgtf_eng = read.csv(link)
sgtf_eng$sgtf = factor(sgtf_eng$sgtf, levels=c("Cases with confirmed S-gene","Cases with confirmed SGTF"), labels=c("no_sgtf","sgtf"))
# data for England by LTLA also at 
# see https://www.gov.uk/government/publications/covid-19-omicron-daily-overview
sgtf_eng = spread(sgtf_eng, sgtf, n)
sgtf_eng[is.na(sgtf_eng)] = 0
sgtf_eng$date = dmy(sgtf_eng$specimen_date) # as.Date(sgtf_eng$specimen_date)

sgtf_eng_byregion = sgtf_eng
sgtf_eng_byregion = sgtf_eng_byregion[,c("UKHSA_region","date","no_sgtf","sgtf")]
sgtf_eng_byregion = sgtf_eng_byregion %>%  # sum counts per date & UKHSA_region
  group_by(UKHSA_region, date) %>% 
  summarise(across(where(is.numeric), list(sum = sum))) 
sgtf_eng_byregion = as.data.frame(sgtf_eng_byregion)
colnames(sgtf_eng_byregion) = c("UKHSA_region", "date","no_sgtf","sgtf")
sgtf_eng_byregion$country = "England"
sgtf_eng_byregion$pos_tests = sgtf_eng_byregion$sgtf + sgtf_eng_byregion$no_sgtf
sgtf_eng_byregion$prop_sgtf = sgtf_eng_byregion$sgtf/sgtf_eng_byregion$pos_tests

sgtf_eng = sgtf_eng %>%  # sum counts over ONS regions
  group_by(date) %>% 
  summarise(across(c(no_sgtf,sgtf), list(sum = sum)))
sgtf_eng = as.data.frame(sgtf_eng)
colnames(sgtf_eng) = c("date","no_sgtf","sgtf")
sgtf_eng$country = "England"
sgtf_eng$pos_tests = sgtf_eng$sgtf + sgtf_eng$no_sgtf
sgtf_eng$no_sgtf = NULL
sgtf_eng$prop_sgtf = sgtf_eng$sgtf/sgtf_eng$pos_tests
sgtf_eng$comment = ""
sgtf_eng$source = "UKHSA, https://www.gov.uk/government/publications/covid-19-omicron-daily-overview"
write.csv(sgtf_eng,".//data//omicron_sgtf//sgtf_england.csv",row.names=F)

sgtf_scot = read.csv(".//data//omicron_sgtf//sgtf_scotland.csv") # SGTF data Scotland from https://www.gov.scot/publications/coronavirus-covid-19-additional-data-and-information/
sgtf_be = read.csv(".//data//omicron_sgtf//sgtf_belgium.csv") # SGTF data Belgium from Federal Test Platform (contact: Emmanuel André)
# sgtf_be = sgtf_be[sgtf_be$date>=as.Date("2021-12-11"),] # we include data from 11 dec onwards

# TO DO still add SGTF data from San Diego (Andersen Laboratory @ Scripps Research): 
# https://raw.githubusercontent.com/andersen-lab/SARS-CoV-2_SGTF_San-Diego/main/SGTF_San_Diego.csv
# https://github.com/andersen-lab/SARS-CoV-2_SGTF_San-Diego

# add data from Geneva
# https://www.hug.ch/laboratoire-virologie

# add data from Amsterdam
# https://twitter.com/ARGOSamsterdam/status/1473047414235967495
# data reported to RIVM: https://www.rivm.nl/coronavirus-covid-19/virus/varianten/omikronvariant?s=09

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
sgtf$prop = sgtf$prop_omicron <- sgtf$omicron/sgtf$pos_tests
sgtf$non_omicron = sgtf$pos_tests-sgtf$omicron
sgtf = sgtf[sgtf$date >= as.Date("2021-11-01"),]
sgtf$obs = as.factor(1:nrow(sgtf)) # for observation-level random effect to take into account overdispersion in binomial GLMM (this would be a quasibinomial GLM)
# sgtf$random = factor(1) # fake random effect to be able to run glmmPQL
names(sgtf)
head(sgtf)[,c("country","date","omicron","prop_omicron")]
# sgtf = sgtf[sgtf$prop_omicron>=0.01,] # only include data where omicron made up >1% of all cases
write.csv(sgtf, ".//data//omicron_sgtf//raw_data_share_omicron_SA_ENG_SCOT_DK_BE.csv", row.names=F)

# 4. ANALYSIS OF SGTF DATA USING LOGISTIC REGRESSION ####

# fit using regular binomial GLM (not taking into account temporal autocorrelation or overdispersion)
# fit_sgtf0 = glm(cbind(omicron, non_omicron) ~ scale(date_num) + country, family=binomial(logit), data=sgtf) 
# fit_sgtf1 = glm(cbind(omicron, non_omicron) ~ scale(date_num) * country, family=binomial(logit), data=sgtf) 
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf1)

# fit using a separate-slope logistic regression that allows for overdispersion via the inclusion of an observation-level random effect
fit_sgtf1 = glmmPQL(cbind(omicron, non_omicron) ~ scale(date_num) * country, family=binomial(logit), random=~1|obs, data=sgtf)
fit_sgtf2 = glmmPQL(cbind(omicron, non_omicron) ~ ns(date_num, df=4, Boundary.knots = c(min(date_num),
                                                                                        max(date_num)+1)) * country, family=binomial(logit), 
                                                  random=~1|obs, 
                                                  data=sgtf[sgtf$country!="South Africa",])
# AIC(fit_sgtf0, fit_sgtf1)
summary(fit_sgtf2)

# # fit using a binomial GAM 
# fit_sgtf2 = gam(cbind(omicron, non_omicron) ~ # s(date_num, k=5, bs="cs", fx=T) + 
#                                               s(date_num, k=5, bs="ps", by=country, fx=T), 
#                 family=quasibinomial(logit), 
#                 data=sgtf[sgtf$country!="South Africa",]) # using quasibinomial to take into account overdispersion
# fit_sgtf2 = gamm(cbind(omicron, non_omicron) ~ s(date_num, bs="cs") + s(date_num, bs="cs", by=country), family=binomial(logit), 
#                  random=list(obs=~1), # observation level random effect to take into account overdispersion 
#                  data=sgtf[sgtf$country!="South Africa",]) 
# library(gamm4)
# fit_sgtf2 = gamm4(cbind(omicron, non_omicron) ~ s(date_num, bs="cs") + s(date_num, bs="cs", by=country), family=binomial(logit), 
#                 random = ~ (1 | obs), # observation level random effect to take into account overdispersion 
#                 data=sgtf[sgtf$country!="South Africa",]) 
# AIC(fit_sgtf2)
# summary(fit_sgtf2)

# # fit using a shape-constrained strictly increasing scam logistic spline
# fit_sgtf2 = scam(cbind(omicron, non_omicron) ~ s(date_num, k=4, bs="mpi", by=country), family=binomial(logit), data=sgtf[sgtf$country!="South Africa",])  # or bs="micx"
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf2)

# note allowing lag-1 temporally autocorrelated residuals would be possible using
# fit_sgtf0 = glmmPQL(cbind(omicron, non_omicron) ~ scale(date_num) + country, family=binomial(logit), random=~1|obs, correlation=corAR1(), data=sgtf) 
# fit_sgtf1 = glmmPQL(cbind(omicron, non_omicron) ~ scale(date_num) * country, family=binomial(logit), random=~1|obs, correlation=corAR1(), data=sgtf) 
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf1)

plot(allEffects(fit_sgtf1, residuals=T))
plot(allEffects(fit_sgtf2, residuals=T))

# growth rate advantage delta_r of Omicron over Delta (difference in growth rate per day)
deltar_sgtf_bycountry = as.data.frame(emtrends(fit_sgtf2, ~ date_num, by="country", var="date_num", at=list(date_num=max(dateseq)) ))[,c(2,3,6,7)] # growth rate advantage of Omicron over Delta
rownames(deltar_sgtf_bycountry) = deltar_sgtf_bycountry$country
deltar_sgtf_bycountry$country = NULL
colnames(deltar_sgtf_bycountry) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_bycountry
#            delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# England  0.1286577        0.07556503         0.1817505
# Scotland 0.1230953        0.05492756         0.1912631
# Denmark  0.1501815        0.07438797         0.2259751
# Belgium  0.3194380        0.24187857         0.3969975

# mean growth rate advantage of Omicron over Delta across countries (mean difference in growth rate per day)
deltar_sgtf = as.data.frame(emtrends(fit_sgtf2, ~ date_num, var="date_num", at=list(date_num=max(dateseq))))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
colnames(deltar_sgtf) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf
#     delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# 1 0.1803432         0.1456789         0.2150075
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
#          transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# England   1.830687            1.426403            2.349558
# Scotland  1.783447            1.294545            2.456988
# Denmark   2.025574            1.418533            2.892389
# Belgium   4.487783            3.116870            6.461673

# mean transmission advantage of Omicron over Delta across countries (using generation time of 4.7 days)
# i.e. mean fold difference in the effective reproduction number Re of Omicron compared to that of Delta
transmadv_sgtf = exp(deltar_sgtf*4.7)
colnames(transmadv_sgtf) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf
#   transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# 1  2.334068            1.983158             2.74707
transmadv_sgtf_char = sapply(transmadv_sgtf, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_char = paste0(transmadv_sgtf_char[1], " [", transmadv_sgtf_char[2], "-", transmadv_sgtf_char[3],"] 95% CLs")

# plot model predictions
emmeans_sgtf = as.data.frame(emmeans(fit_sgtf2, ~ date_num+country, at=list(date_num=dateseq, country=levels_country[!levels_country=="South Africa"]),
                                     data=sgtf[sgtf$country!="South Africa",] ), type="response")
colnames(emmeans_sgtf) = c("date_num", "country", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
emmeans_sgtf$date = as.Date(emmeans_sgtf$date_num, origin="1970-01-01")

write.csv(emmeans_sgtf, ".//data//omicron_sgtf//fit_logistic_mixed_model_share_omicron_SA_ENG_SCOT_DK_BE.csv", row.names=F)

# plot of share of Omicron among confirmed infections (on logit scale)
qplot(data=sgtf[sgtf$country!="South Africa",], x=date, y=prop, geom="point", colour=country, fill=country, size=pos_tests, alpha=(1/pos_tests), shape=I(16)) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, size=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf, aes(y=prob), alpha=I(0.5), size=1) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_size_continuous("nr. of positive tests", range=c(2,7), breaks=c(100, 1000, 10000)) +
  scale_alpha_continuous("", range=c(0.5,1), guide=FALSE) +
  xlim(date.from, date.to) +
  # scale_size_continuous(range=c()) +
  ggtitle("Spread of Omicron variant inferred from\nS gene target failure (SGTF) data (UK, BE)\n& variant PCR data (DK)",
          "(SGTF counts adjusted by proportion that was estimated to be Omicron\nbased on GISAID data, data UKHSA, Public Health Scotland,\nStatens Serum Institut, Federal Test Platform Belgium)") + # Lesley Scott & NHLS team, 
  ylab("Share of Omicron among confirmed cases (%)") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(0.94),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char,
                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),NA), ylim=c(0.001,0.99)) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data.png"), width=8, height=6)

# plot of share of Omicron among confirmed infections (on linear scale)
qplot(data=sgtf[sgtf$country!="South Africa",], x=date, y=100*prop, geom="point", colour=country, fill=country, size=pos_tests, alpha=(1/pos_tests), shape=I(16)) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob*100, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, size=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf, aes(y=prob*100), alpha=I(0.5), size=1) +
  scale_size_continuous("nr. of positive tests", range=c(2,7), breaks=c(100, 1000, 10000)) +
  scale_alpha_continuous("", range=c(0.5,1), guide=FALSE) +
  xlim(date.from, date.to) +
# scale_size_continuous(range=c()) +
  ggtitle("Spread of Omicron variant inferred from\nS gene target failure (SGTF) data (UK, BE)\n& variant PCR data (DK)",
          "(SGTF counts adjusted by proportion that was estimated to be Omicron\nbased on GISAID data, data UKHSA, Public Health Scotland,\nStatens Serum Institut, Federal Test Platform Belgium)") + # Lesley Scott & NHLS team, 
  ylab("Share of Omicron among confirmed cases (%)") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(80),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char,
                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),NA), ylim=c(0,100)) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data_linear scale.png"), width=8, height=6)


emmeans_sgtf[emmeans_sgtf$date==today,] # estimated prop Omicron among confirmed cases today
# date_num  country      prob         SE  df asymp.LCL asymp.UCL       date
# 119    18989  England 0.9375624 0.01291325 114 0.9065420 0.9587550 2021-12-28
# 268    18989 Scotland 0.8246276 0.03662615 114 0.7400652 0.8859206 2021-12-28
# 417    18989  Denmark 0.8668608 0.03707203 114 0.7750703 0.9248260 2021-12-28
# 566    18989  Belgium 0.7754987 0.05492457 114 0.6490051 0.8658300 2021-12-28

emmeans_sgtf[emmeans_sgtf$date==(today+7),] # estimated prop Omicron among new infections today (1 week before diagnosis)
# date_num  country      prob         SE  df asymp.LCL asymp.UCL       date
# 126    18996  England 0.9736535 0.01035461 114 0.9432186 0.9879831 2022-01-04
# 275    18996 Scotland 0.9175630 0.03691390 114 0.8089089 0.9669600 2022-01-04
# 424    18996  Denmark 0.9490566 0.02828403 114 0.8539419 0.9834332 2022-01-04
# 573    18996  Belgium 0.9699881 0.01703548 114 0.9102294 0.9903867 2022-01-04

emmeans_sgtf[emmeans_sgtf$date==as.Date("2022-01-02"),]
emmeans_sgtf[emmeans_sgtf$date==as.Date("2021-12-24"),]





# 5. ANALYSIS OF ENGLISH DATA AT THE LEVEL OF ONS REGIONS ####

# fit using a separate-slope logistic regression that allows for overdispersion via the inclusion of an observation-level random effect

sgtf_eng_byregion$obs = factor(1:nrow(sgtf_eng_byregion))
sgtf_eng_byregion$date_num = as.numeric(sgtf_eng_byregion$date)
fit_sgtf_eng_byregion1 = glmmPQL(cbind(sgtf, no_sgtf) ~ ns(date_num, df=2) * UKHSA_region, family=binomial(logit), random=~1|obs, 
                                 data=sgtf_eng_byregion[sgtf_eng_byregion$date>=as.Date("2021-12-01"),]) 
# AIC(fit_sgtf_eng_byregion1)
summary(fit_sgtf_eng_byregion1)

plot(allEffects(fit_sgtf_eng_byregion1, residuals=T))

# growth rate advantage delta_r of Omicron over Delta (difference in growth rate per day)
deltar_sgtf_eng_byregion = as.data.frame(emtrends(fit_sgtf_eng_byregion1, ~ date_num, by="UKHSA_region", var="date_num", at=list(date_num=max(dateseq))))[,c(2,3,6,7)] # growth rate advantage of Omicron over Delta
rownames(deltar_sgtf_eng_byregion) = deltar_sgtf_eng_byregion$UKHSA_region
deltar_sgtf_eng_byregion$UKHSA_region = NULL
colnames(deltar_sgtf_eng_byregion) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_eng_byregion
#                         delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# East Midlands        0.17940010       0.090311103         0.2684891
# East of England      0.13702996       0.046373789         0.2276861
# London               0.08341477      -0.006912061         0.1737416
# North East           0.25822436       0.157210443         0.3592383
# North West           0.25393983       0.167987928         0.3398917
# South East           0.15935888       0.070448951         0.2482688
# South West           0.20755803       0.107842160         0.3072739
# West Midlands        0.21990413       0.129395599         0.3104127
# Yorkshire and Humber 0.19620127       0.100645121         0.2917574

# mean growth rate advantage of Omicron over Delta across ONS regions (mean difference in growth rate per day)
deltar_sgtf_eng_avgoverregions = as.data.frame(emtrends(fit_sgtf_eng_byregion1, ~ date_num, var="date_num", at=list(date_num=max(dateseq))))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
colnames(deltar_sgtf_eng_avgoverregions) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_eng_avgoverregions
#     delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# 1 0.1883368         0.1574893         0.2191843
deltar_sgtf_eng_avgoverregions_char = sapply(deltar_sgtf_eng_avgoverregions, function (x) sprintf(as.character(round(x,2)), 2) )
deltar_sgtf_eng_avgoverregions_char = paste0(deltar_sgtf_eng_avgoverregions_char[1], " [", deltar_sgtf_eng_avgoverregions_char[2], "-", deltar_sgtf_eng_avgoverregions_char[3],"] 95% CLs")

# transmission advantage of Omicron over Delta (using generation time of 4.7 days)
# i.e. how much higher the effective reproduction number Re value of Omicron is than that of Delta
# and how many more people the virus infects over the course of a single generation
# (due to higher immune escape, i.e. more frequent reinfection of people with immunity built up
# as a result of vaccination or prior exposure)
transmadv_sgtf_eng_byregion = exp(deltar_sgtf_eng_byregion*4.7) 
colnames(transmadv_sgtf_eng_byregion) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf_eng_byregion
#                      transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# East Midlands         2.323746           1.5287680            3.532122
# East of England       1.904160           1.2435334            2.915743
# London                1.480011           0.9680353            2.262760
# North East            3.365762           2.0936084            5.410924
# North West            3.298663           2.2023903            4.940622
# South East            2.114856           1.3925131            3.211902
# South West            2.652553           1.6600673            4.238406
# West Midlands         2.811026           1.8370469            4.301396
# Yorkshire and Humber  2.514681           1.6048528            3.940312

# mean transmission advantage of Omicron over Delta across countries (using generation time of 4.7 days)
# i.e. mean fold difference in the effective reproduction number Re of Omicron compared to that of Delta
transmadv_sgtf_eng_avgoverregions = exp(deltar_sgtf_eng_avgoverregions*4.7)
colnames(transmadv_sgtf_eng_avgoverregions) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf_eng_avgoverregions
#   transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# 1  2.423428            2.096355            2.801531
transmadv_sgtf_eng_avgoverregions_char = sapply(transmadv_sgtf_eng_avgoverregions, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_eng_avgoverregions_char = paste0(transmadv_sgtf_eng_avgoverregions_char[1], " [", transmadv_sgtf_eng_avgoverregions_char[2], "-", transmadv_sgtf_eng_avgoverregions_char[3],"] 95% CLs")

# plot model predictions
emmeans_sgtf_eng_byregion = as.data.frame(emmeans(fit_sgtf_eng_byregion1, ~ date_num+UKHSA_region, at=list(date_num=dateseq, UKHSA_region=unique(sgtf_eng_byregion$UKHSA_region) )), type="response")
colnames(emmeans_sgtf_eng_byregion) = c("date_num", "UKHSA_region", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
emmeans_sgtf_eng_byregion$date = as.Date(emmeans_sgtf_eng_byregion$date_num, origin="1970-01-01")
propomicron_byregion = emmeans_sgtf_eng_byregion[emmeans_sgtf_eng_byregion$date==as.Date("2021-12-07"),] # estimated prop Omicron among confirmed cases on Dec 10
levels_regions = levels(propomicron_today_byregion$UKHSA_region)[order(propomicron_byregion$prob, decreasing=T)] # sort regions by order of introduction of Omicron
emmeans_sgtf_eng_byregion$UKHSA_region = factor(emmeans_sgtf_eng_byregion$UKHSA_region,
                                                levels=levels_regions)
sgtf_eng_byregion$UKHSA_region = factor(sgtf_eng_byregion$UKHSA_region,
                                        levels=levels_regions)
write.csv(emmeans_sgtf_eng_byregion, ".//data//omicron_sgtf//fit_logistic_mixed_model_share_omicron_ENG_byregion.csv", row.names=F)

# plot of share of Omicron among confirmed infections (on logit scale)
qplot(data=sgtf_eng_byregion[sgtf_eng_byregion$date>=as.Date("2021-12-01"),], 
      x=date, y=prop_sgtf, geom="point", colour=UKHSA_region, fill=UKHSA_region, size=pos_tests, alpha=(1/pos_tests), shape=I(16)) + 
  # facet_wrap(~ UKHSA_region) +
  geom_ribbon(data=emmeans_sgtf_eng_byregion, aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, size=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf_eng_byregion, aes(y=prob), alpha=I(0.5), size=1) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_size_continuous("nr. of positive tests", range=c(2,7), breaks=c(100, 1000, 10000)) +
  scale_alpha_continuous("", range=c(0.5,1), guide=FALSE) +
  xlim(date.from, date.to) +
  # scale_size_continuous(range=c()) +
  ggtitle("Spread of the Omicron variant in England,\nusing S gene target failure (SGTF) data as a proxy",
          "(data UKHSA, mixed logistic 2 df spline fit with observation-level random effect to take into account overdispersion)") +  
  ylab("Share of tests with S dropout (% Omicron)") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(0.94),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_eng_avgoverregions_char,
                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_eng_avgoverregions_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),as.Date("2022-01-01")), ylim=c(0.001,0.99))  +
  scale_colour_hue(h=c(10, 310), c=100) +
  scale_fill_hue(h=c(10, 310), c=100)

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data england by region.png"), width=8, height=6)

# plot of share of Omicron among confirmed infections (on linear scale)
qplot(data=sgtf_eng_byregion[sgtf_eng_byregion$date>=as.Date("2021-12-01"),], x=date, y=100*prop_sgtf, 
      geom="point", colour=UKHSA_region, fill=UKHSA_region, size=pos_tests, alpha=(1/pos_tests), shape=I(16)) + 
  geom_ribbon(data=emmeans_sgtf_eng_byregion, aes(y=prob*100, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, size=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf_eng_byregion, aes(y=prob*100), alpha=I(0.5), size=1) +
  scale_size_continuous("nr. of positive tests", range=c(2,7), breaks=c(100, 1000, 10000)) +
  scale_alpha_continuous("", range=c(0.5,1), guide=FALSE) +
  xlim(date.from, date.to) +
  # scale_size_continuous(range=c()) +
  ggtitle("Spread of the Omicron variant in England,\nusing S gene target failure (SGTF) data as a proxy",
          "(data UKHSA, mixed logistic 2 df spline fit with observation-level random effect to take into account overdispersion)") +  
  ylab("Share of tests with S dropout (% Omicron)") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(80),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_eng_avgoverregions_char,
                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_eng_avgoverregions_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),as.Date("2022-01-01")), ylim=c(0,100)) +
  scale_colour_hue(h=c(10, 310), c=100) +
  scale_fill_hue(h=c(10, 310), c=100)


ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data_linear scale.png"), width=8, height=6)







# PART BELOW NOT FINISHED YET

# 6. PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
install.packages("covidregionaldata",
                 repos = "https://epiforecasts.r-universe.dev"
)
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales) 

# cases_tot = as.data.frame(get_national_data(countries = levels_country, source="WHO"))
# cases_tot = cases_tot[,c("date","country","cases_new","hosp_new","deaths_new","tested_new")]
# get_available_datasets("national")
cases_tot_level2 = as.data.frame(get_national_data(countries = c("United Kingdom"), level="2", source="Google")) # get data England & Scotland
cases_tot_eng_scot = cases_tot_level2[cases_tot_level2$subregion %in% c("England","Scotland"),c("date","subregion","cases_new","hosp_new","deaths_new","tested_new")]
colnames(cases_tot_eng_scot)[2] = "country"

cases_tot_sa = read.csv("https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=cumulative")
cases_tot_sa = cases_tot_sa[cases_tot_sa$date!="",]
cases_tot_sa$tests_daily = abs(cases_tot_sa$tests_daily)
cases_tot_sa$country="South Africa"
cases_tot_sa$cases_new = cases_tot_sa$cases_daily
cases_tot_sa$hosp_new = NA
cases_tot_sa$deaths_new = cases_tot_sa$deaths_daily
cases_tot_sa$tested_new = cases_tot_sa$new_public_tests+cases_tot_sa$new_private_tests
cases_tot_sa = cases_tot_sa[,c("date","country","cases_new","hosp_new","deaths_new","tested_new")]
cases_tot_sa$date = dmy(cases_tot_sa$date)

cases_tot = as.data.frame(get_national_data(countries = levels_country[!levels_country %in% c("England","Scotland","South Africa")], source="WHO"))
cases_tot = cases_tot[,c("date","country","cases_new","hosp_new","deaths_new","tested_new")]

cases_tot = rbind(cases_tot_eng_scot, cases_tot_sa, cases_tot)

cases_tot = cases_tot[cases_tot$date>=as.Date("2021-09-01"),]
cases_tot$date_num = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
# cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
cases_tot$country = factor(cases_tot$country, levels=levels_country)

# smooth out weekday effects in case nrs using GAM (& potentially correct for unequal testing intensity)
library(mgcv)
k=7
fit_cases = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F, by=country) + country +
                  WEEKDAY, # + 
                # s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=poisson(log), data=cases_tot,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+0.5*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-0.5*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 
BIC(fit_cases)


# STACKED AREA CHART OF NEW CASES BY VARIANT (LOGISTIC OR MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####
emmeans_sgtf2 = emmeans_sgtf
emmeans_sgtf2$prob = 1-emmeans_sgtf2$prob
emmeans_sgtf2$asymp.LCL = 1-emmeans_sgtf2$asymp.LCL
emmeans_sgtf2$asymp.UCL = 1-emmeans_sgtf2$asymp.UCL

emmeans_lineages = rbind(data.frame(LINEAGE="Omicron", emmeans_sgtf),
                         data.frame(LINEAGE="Other", emmeans_sgtf2))
emmeans_lineages$totcases = cases_tot$cases_new[match(interaction(emmeans_lineages$date_num,emmeans_lineages$country),
                                                            interaction(cases_tot$date_num,cases_tot$country))]
emmeans_lineages$cases = emmeans_lineages$totcases * emmeans_lineages$prob
emmeans_lineages$cases[emmeans_lineages$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ date_num|country, at=list(date_num=seq(as.numeric(date.from), as.numeric(date.to), by=1)
), type="response"))
emmeans_lineages$smoothed_totcases = cases_emmeans$rate[match(interaction(emmeans_lineages$date_num,emmeans_lineages$country),
                                                                    interaction(cases_emmeans$date_num,cases_emmeans$country))]
emmeans_lineages$smoothed_cases = emmeans_lineages$smoothed_totcases * emmeans_lineages$prob
emmeans_lineages$smoothed_cases[emmeans_lineages$smoothed_cases<=0.001] = NA
lineage_cols = c("red2","grey50")

ggplot(data=emmeans_lineages, 
       aes(x=date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ country, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-09-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day by variant_stacked area chart_raw case data.png"), width=8, height=6)

ggplot(data=emmeans_lineages, 
       aes(x=date, y=cases+1, group=LINEAGE)) + 
  facet_wrap(~ country) +
  geom_col(data=emmeans_lineages[emmeans_lineages$LINEAGE=="Other",], aes(lwd=I(1.2), colour=NULL), fill=lineage_cols[2], position="identity", alpha=I(0.4), width=I(2)) +
  geom_col(data=emmeans_lineages[emmeans_lineages$LINEAGE=="Omicron",], aes(lwd=I(1.2), colour=NULL), fill=lineage_cols[1], position="identity", alpha=I(0.4), width=I(2)) +
  xaxis +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-09-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day by variant_stacked area chart_raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=emmeans_lineages,
       aes(x=date-7, y=smoothed_cases, group=LINEAGE)) +
  facet_wrap(~ country, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-09-01"),today)) # as.Date("2022-01-01")

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day by variant_stacked area chart_smoothed.png"), width=8, height=6)

ggplot(data=emmeans_lineages, 
       aes(x=date, y=smoothed_cases+1, group=LINEAGE)) + 
  facet_wrap(~ country, scale="free") +
  geom_area(data=emmeans_lineages[emmeans_lineages$LINEAGE=="Other",], aes(lwd=I(1.2), colour=NULL), fill=lineage_cols[2], position="identity", alpha=I(0.4), width=I(2)) +
  geom_area(data=emmeans_lineages[emmeans_lineages$LINEAGE=="Omicron",], aes(lwd=I(1.2), colour=NULL), fill=lineage_cols[1], position="identity", alpha=I(0.4), width=I(2)) +
  xaxis +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-09-01"),today)) # as.Date("2022-01-01")

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=emmeans_lineages, 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE, group=LINEAGE)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10()

# ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=emmeans_lineages,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE)) +
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE, group=LINEAGE)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),today), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit smoothed case data.png"), width=8, height=6)






# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####
# TO DO : need to finish this part

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

# calculate r from Re
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L127
r.from.Re <- function(Re, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  r <- (Re^k - 1) / (k * gamma_mean)
  return(r)
}



# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, by="country", var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to, by=3)
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[3] = "r"
colnames(avg_r_cases)[6] = "r_LOWER"
colnames(avg_r_cases)[7] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  facet_wrap(~ country) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES IN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants1 = do.call(rbind, lapply(levels_countries, function(country) { do.call(rbind, 
                                                                                            lapply(seq(date.from,
                                                                                                       date.to, by=3), 
                                                                                                   function (d) { 
                                                                                                     wt = as.data.frame(emmeans(fit1_africa_multi, ~ LINEAGE , by="country", 
                                                                                                                                at=list(DATE_NUM=d, country=country), type="response"))$prob   # important: these should sum to 1
                                                                                                     # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) 
                                                                                                     # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
                                                                                                     cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                                                                     names(cons) = seq_along(cons)
                                                                                                     EMT = emtrends(fit1_africa_multi,  ~ LINEAGE , by=c("DATE_NUM", "country"),
                                                                                                                    var="DATE_NUM", mode="latent",
                                                                                                                    at=list(DATE_NUM=d, country=country))
                                                                                                     out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                                                                     # sum(out$estimate*wt) # should sum to zero
                                                                                                     return(out) } )) } ))
above_avg_r_variants = above_avg_r_variants1
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE), 
                                       labels=levels_LINEAGE)
above_avg_r_variants$LINEAGE = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2021-01-01" "2021-07-30"
# average growth rate of all lineages calculated from case nrs
above_avg_r_variants$avg_r = avg_r_cases$r[match(interaction(above_avg_r_variants$collection_date,above_avg_r_variants$country),
                                                 interaction(avg_r_cases$DATE,avg_r_cases$country))]  
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                country=avg_r_cases$country,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE,
                LINEAGE="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$LINEAGE = factor(above_avg_r_variants$LINEAGE, levels=c(levels_LINEAGE,"avg"))
above_avg_r_variants$prob = emmeans_lineages$prob[match(interaction(round(above_avg_r_variants$DATE_NUM),
                                                                          as.character(above_avg_r_variants$LINEAGE),
                                                                          as.character(above_avg_r_variants$country)),
                                                              interaction(round(emmeans_lineages$DATE_NUM),
                                                                          as.character(emmeans_lineages$LINEAGE),
                                                                          as.character(emmeans_lineages$country)))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 2
ymin = 1/2
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$LINEAGE %in% c("other"))|(above_avg_r_variants2$collection_date>=max(cases_tot$date))),], # |above_avg_r_variants2$collection_date>max(cases_tot$DATE)
      x=collection_date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE, fill=LINEAGE, alpha=I(0.5),
      group=LINEAGE, linetype=I(0)) +
  facet_wrap(~ country) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES BY VARIANT IN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-11-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  scale_colour_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)



