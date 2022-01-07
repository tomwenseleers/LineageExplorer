# ANALYSIS OF GROWTH RATE ADVANTAGE OF OMICRON IN South Africa, England, Denmark & Belgium BASED ON SGTF (S gene target failure / S dropout) DATA
# data: 
# Belgium: Federal Test Platform, provided by Emmanuel André
# Denmark: Statens Serum Institut (https://covid19.ssi.dk/virusvarianter/delta-pcr)
# England: UK Health Security Agency (https://www.gov.uk/government/publications/covid-19-omicron-daily-overview)
# Scotland: https://www.gov.scot/publications/coronavirus-covid-19-additional-data-and-information/
# South Africa: Lesley Scott & NHLS team (traced off graph by Alex Selby)
# GISAID data to infer the proportion of S dropout (with spike Spike_H69del/Spike_V70del) that are Omicron (https://www.gisaid.org/)

# T. Wenseleers
# last update 6 JANUARY 2021

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

# FUNCTIONS TO CALCULATE EFFECTIVE REPRODUCTION NUMBER R FROM
# THE INSTANTANEOUS GROWTH RATE r (per day), ASSUMING GAMMA DISTRIBUTED GENERATION TIME
# from epiforecasts package growth_to_R
# https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# see https://www.sciencedirect.com/science/article/pii/S1755436518300847?via%3Dihub and
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) {
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

# for Alpha variant: mean & SD of generation time taken as the mean
# of the intrinsic generation time (mean 5.5 days, SD 3.8 days, k=(3.8/5.5)^2=0.48) & 
# the household generation time (mean 4.5 days, SD 3.3 days, k=(3.3/4.5)^2=0.54), i.e.
# mean GT = mean(c(5.5, 4.5)) = 5.0 days and SD = mean(c(3.8, 3.3)) = 3.55 days (k=0.50)
# given in suppl Table S3 of Hart et al. 2022, https://www.medrxiv.org/content/10.1101/2021.10.21.21265216v1)

# for Delta variant: mean & SD of generation time taken as the mean
# of the intrinsic generation time (mean 4.6 days, SD 3.1 days, k=(3.1/4.6)^2=0.45) & 
# the household generation time (mean 3.2 days, SD 2.4 days, k=(2.4/3.2)^2=0.56), i.e.
# mean GT = mean(c(4.6, 3.2)) = 3.9 days and SD = mean(c(3.1, 2.4)) = 2.75 days (k=0.50)
# given in suppl Table S3 of Hart et al. 2022, https://www.medrxiv.org/content/10.1101/2021.10.21.21265216v1)

# for Omicron variant: mean & SD of generation time taken as given in
# Kim et al. 2022, https://www.medrxiv.org/content/10.1101/2021.12.25.21268301v1
# also see https://www.jkms.org/DOIx.php?id=10.3346/jkms.2021.36.e346
# Given that another S Korean study reported a serial interval of 3.26 days 
# (https://academic.oup.com/jid/advance-article/doi/10.1093/infdis/jiab586/6448309?login=true),
# one could also opt for a more modest 32% reduction in the mean and SD of the generation
# time from Delta to Omicron, i.e. mean GT = 3.13 days and SD = 2.22 days, as in
# https://github.com/fvalka/covid19-austria-omicron/blob/master/reports/omicron_austria.pdf

Re.from.r.variant <- function(r, variant="Delta") {
  gamma_mean = 5.0*(variant=="Alpha") + 3.9*(variant=="Delta") + 2.22*(variant=="Omicron")
  gamma_sd = 3.55*(variant=="Alpha") + 2.75*(variant=="Delta") + 1.62*(variant=="Omicron")
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}



today = as.Date(Sys.time()) 
# today = as.Date("2021-12-10")
today_num = as.numeric(today)
plotdir = "omicron_sgtf"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

tag = paste("@TWenseleers\n",today)

extrapolate = 30
date.from = as.Date("2021-08-01")
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
varpcr_dk$region = varpcr_dk$country
sgtf_sa = read.csv(".//data//omicron_sgtf//sgtf_south africa.csv") # data traced from graph (data not open at the moment)
sgtf_sa$date = as.Date(sgtf_sa$date)
sgtf_sa = sgtf_sa[sgtf_sa$date>=as.Date("2021-11-05"),] # data limited to >= Nov 5 when from sequencing data nearly all S dropout was Omicron
sgtf_sa$region = sgtf_sa$country

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
colnames(sgtf_eng_byregion) = c("region", "date","no_sgtf","sgtf")
sgtf_eng_byregion$country = "England"
sgtf_eng_byregion$pos_tests = sgtf_eng_byregion$sgtf + sgtf_eng_byregion$no_sgtf
sgtf_eng_byregion$prop_sgtf = sgtf_eng_byregion$sgtf/sgtf_eng_byregion$pos_tests
sgtf_eng_byregion$no_sgtf = NULL
sgtf_eng_byregion$comment = NA
sgtf_eng_byregion$source = "UKHSA, https://www.gov.uk/government/publications/covid-19-omicron-daily-overview"

sgtf_eng = sgtf_eng %>%  # sum counts over ONS regions
  group_by(date) %>% 
  summarise(across(c(no_sgtf,sgtf), list(sum = sum)))
sgtf_eng = as.data.frame(sgtf_eng)
colnames(sgtf_eng) = c("date","no_sgtf","sgtf")
sgtf_eng$country = "England"
sgtf_eng$region = sgtf_eng$country
sgtf_eng$pos_tests = sgtf_eng$sgtf + sgtf_eng$no_sgtf
sgtf_eng$no_sgtf = NULL
sgtf_eng$prop_sgtf = sgtf_eng$sgtf/sgtf_eng$pos_tests
sgtf_eng$comment = ""
sgtf_eng$source = "UKHSA, https://www.gov.uk/government/publications/covid-19-omicron-daily-overview"
write.csv(sgtf_eng,".//data//omicron_sgtf//sgtf_england.csv",row.names=F)

# sgtf_eng_byltla = read.csv(link_byLTLA)
# sgtf_eng_byltla$sgtf = factor(sgtf_eng_byltla$sgtf, levels=c("Cases with detectable S-gene","Cases with SGTF"), labels=c("no_sgtf","sgtf"))
# sgtf_eng_byltla = spread(sgtf_eng_byltla, sgtf, n)
# sgtf_eng_byltla[is.na(sgtf_eng_byltla)] = 0
# sgtf_eng_byltla$date = dmy(sgtf_eng_byltla$specimen_date) # as.Date(sgtf_eng$specimen_date)

sgtf_scot = read.csv(".//data//omicron_sgtf//sgtf_scotland.csv") # SGTF data Scotland from https://www.gov.scot/publications/coronavirus-covid-19-additional-data-and-information/
sgtf_scot$region = sgtf_scot$country
sgtf_scot = sgtf_scot[complete.cases(sgtf_scot),]
sgtf_be = read.csv(".//data//omicron_sgtf//sgtf_belgium.csv") # SGTF data Belgium from Federal Test Platform (contact: Emmanuel André)
sgtf_be$region = sgtf_be$country
# sgtf_be = sgtf_be[sgtf_be$date>=as.Date("2021-12-11"),] # we include data from 11 dec onwards

# TO DO still add SGTF data from San Diego (Andersen Laboratory @ Scripps Research): 
# https://raw.githubusercontent.com/andersen-lab/SARS-CoV-2_SGTF_San-Diego/main/SGTF_San_Diego.csv
# https://github.com/andersen-lab/SARS-CoV-2_SGTF_San-Diego

# add data from Geneva
# https://www.hug.ch/laboratoire-virologie

# add data from Amsterdam
# https://twitter.com/ARGOSamsterdam/status/1473047414235967495
# data reported to RIVM: https://www.rivm.nl/coronavirus-covid-19/virus/varianten/omikronvariant?s=09

# sgtf = rbind(sgtf_sa, sgtf_eng_byregion[sgtf_eng_byregion$region=="London",], sgtf_scot, sgtf_be) # just using data for London here
sgtf = rbind(sgtf_sa, sgtf_eng, sgtf_scot, sgtf_be)
sgtf$DATE_NUM = as.numeric(sgtf$date)
# Omicron count is sgtf count * prop sgtf samples inferred to be Omicron (from logistic fit on GISAID data)
prop_sgtf_omicron = predict(fit_sgtf_isomicron, newdata=sgtf, type="response")
# sgtf$country[sgtf$country=="England"] = "England (London)"
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
head(sgtf)[,c("country","date","omicron","non_omicron","prop_omicron")]
# sgtf = sgtf[sgtf$prop_omicron>=0.01,] # only include data where omicron made up >1% of all cases

sgtf2=sgtf
sgtf2$region = NULL
write.csv(sgtf2, ".//data//omicron_sgtf//raw_data_share_omicron_SA_ENG_SCOT_DK_BE.csv", row.names=F)

# 4. ANALYSIS OF SGTF DATA USING LOGISTIC REGRESSION ####

# fit using regular binomial GLM (not taking into account temporal autocorrelation or overdispersion)
# fit_sgtf0 = glm(cbind(omicron, non_omicron) ~ scale(date_num) + country, family=binomial(logit), data=sgtf) 
# fit_sgtf1 = glm(cbind(omicron, non_omicron) ~ scale(date_num) * country, family=binomial(logit), data=sgtf) 
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf1)

# fit using a separate-slope logistic regression that allows for overdispersion via the inclusion of an observation-level random effect
fit_sgtf2 = glmmPQL(cbind(omicron, non_omicron) ~ ns(date_num, df=5, Boundary.knots = c(min(date_num)-0.75,
                                                                                        max(date_num)+0.75)) * country, family=binomial(logit), 
                                                  random=~1|obs, 
                                                  data=sgtf[sgtf$country!="South Africa",])
# AIC(fit_sgtf0, fit_sgtf1)
summary(fit_sgtf2)

# # # fit using a binomial GAM 
# fit_sgtf2 = gam(cbind(omicron, non_omicron) ~ 
#                                               s(date_num, k=5, bs="cs") +
#                                               # s(date_num, k=4, bs="cs", by=country, fx=F)+
#                                               country+
#                                               country*date_num,
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

# plot(allEffects(fit_sgtf1, residuals=T))
plot(allEffects(fit_sgtf2, residuals=T))

# growth rate advantage delta_r / s of Omicron over Delta (difference in growth rate per day)
deltar_sgtf_bycountry = as.data.frame(emtrends(fit_sgtf2, ~ date_num, by="country", var="date_num", at=list(date_num=today_num) ))[,c(2,3,6,7)] # growth rate advantage of Omicron over Delta
rownames(deltar_sgtf_bycountry) = deltar_sgtf_bycountry$country
deltar_sgtf_bycountry$country = NULL
colnames(deltar_sgtf_bycountry) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_bycountry
#            delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# England  0.2367685        0.07567385         0.3978632
# Scotland 0.1716302        0.06932510         0.2739352
# Denmark  0.1572643        0.08042112         0.2341074
# Belgium  0.1899732        0.13390528         0.2460411

# mean growth rate advantage of Omicron over Delta across countries (mean difference in growth rate per day)
deltar_sgtf = as.data.frame(emtrends(fit_sgtf2, ~ date_num, var="date_num", at=list(date_num=today_num)))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
colnames(deltar_sgtf) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf
#     delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# 1 0.188909          0.135602         0.2422161
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
# England   3.042903            1.427132            6.488017
# Scotland  2.240416            1.385177            3.623700
# Denmark   2.094138            1.459333            3.005082
# Belgium   2.442138            1.876400            3.178449

# mean transmission advantage of Omicron over Delta across countries (using generation time of 4.7 days)
# i.e. mean fold difference in the effective reproduction number Re of Omicron compared to that of Delta
transmadv_sgtf = exp(deltar_sgtf*4.7)
colnames(transmadv_sgtf) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf
#   transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# 1 2.429954            1.891423            3.121818
transmadv_sgtf_char = sapply(transmadv_sgtf, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_char = paste0(transmadv_sgtf_char[1], " [", transmadv_sgtf_char[2], "-", transmadv_sgtf_char[3],"] 95% CLs")

# plot model predictions
emmeans_sgtf = as.data.frame(emmeans(fit_sgtf2, ~ date_num+country, at=list(date_num=dateseq, country=levels_country[!levels_country=="South Africa"]),
                                     data=sgtf[sgtf$country!="South Africa",] ), type="response")
colnames(emmeans_sgtf) = c("date_num", "country", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
emmeans_sgtf$date = as.Date(emmeans_sgtf$date_num, origin="1970-01-01")

write.csv(emmeans_sgtf, ".//data//omicron_sgtf//fit_logistic_mixed_model_share_omicron_SA_ENG_SCOT_DK_BE.csv", row.names=F)

# growth rate advantage s of Omicron over Delta (=logistic growth rate/slope=difference in growth rate Omicron & Delta)
emtrends_sgtf = as.data.frame(emtrends(fit_sgtf2, ~ date_num, var="date_num", by="country",
                                       at=list(date_num=dateseq, country=levels_country[!levels_country=="South Africa"])), type="link")
colnames(emtrends_sgtf) = c("date_num", "country", "s", "SE", "df", "s.LCL", "s.UCL")
emtrends_sgtf$date = as.Date(emtrends_sgtf$date_num, origin="1970-01-01")
emtrends_sgtf$prop_omicron = emmeans_sgtf$prob[match(interaction(emtrends_sgtf$country, emtrends_sgtf$date),
                                                     interaction(emmeans_sgtf$country, emmeans_sgtf$date))]
emtrends_sgtf$s[emtrends_sgtf$prop_omicron<0.01] = NA
emtrends_sgtf$s.LCL[emtrends_sgtf$prop_omicron<0.01] = NA
emtrends_sgtf$s.UCL[emtrends_sgtf$prop_omicron<0.01] = NA


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
  xlab("Date of diagnosis") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(0.98),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char, " per day"), #,
#                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  # xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),NA), ylim=c(0.001,0.99)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

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
  annotate("text", x = c(as.Date("2021-12-01")), y = c(90),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char," per day"),
  #                         " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
          color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + 
  # xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),NA), ylim=c(0,100)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data_linear scale.png"), width=8, height=6)


emmeans_sgtf[emmeans_sgtf$date==today,] # estimated prop Omicron among confirmed cases today
# date_num  country      prob          SE  df asymp.LCL asymp.UCL       date
# 159    18998  England 0.9943191 0.004755222 136 0.9706926 0.9989200 2022-01-06
# 348    18998 Scotland 0.9592776 0.018398255 136 0.9027375 0.9835491 2022-01-06
# 537    18998  Denmark 0.9674316 0.010667004 136 0.9383005 0.9830569 2022-01-06
# 726    18998  Belgium 0.9359844 0.013598502 136 0.9032316 0.9581650 2022-01-06

emmeans_sgtf[emmeans_sgtf$date==(today+7),] # estimated prop Omicron among new infections today (1 week before diagnosis)
# date_num  country      prob          SE  df asymp.LCL asymp.UCL       date
# 166    19005  England 0.9989120 0.001531203 136 0.9826417 0.9999328 2022-01-13
# 355    19005 Scotland 0.9873930 0.010258482 136 0.9388373 0.9975039 2022-01-13
# 544    19005  Denmark 0.9889274 0.006627752 136 0.9642631 0.9966288 2022-01-13
# 733    19005  Belgium 0.9822295 0.007306245 136 0.9602472 0.9921554 2022-01-13

emmeans_sgtf[emmeans_sgtf$date==as.Date("2022-01-03"),]
emmeans_sgtf[emmeans_sgtf$date==as.Date("2021-12-25"),]


# plot of growth rate advantage of Omicron over Delta (difference in growth rate per day)
qplot(data=emtrends_sgtf[emtrends_sgtf$date>=as.Date("2021-12-01"),], x=date, y=s, ymin=s.LCL, ymax=s.UCL, geom="line", colour=country, fill=country) + 
  geom_ribbon(aes(colour=NULL), alpha=I(0.5)) +
  geom_line(alpha=I(0.5), size=1) +
  scale_size_continuous("nr. of positive tests", range=c(2,7), breaks=c(100, 1000, 10000)) +
  scale_alpha_continuous("", range=c(0.5,1), guide=FALSE) +
  xlim(date.from, date.to) +
  # scale_size_continuous(range=c()) +
  ggtitle("Growth rate advantage of Omicron over Delta",
          "inferred from slope of logistic spline fit to S gene target\nfailure (SGTF) data (UK, BE) & variant PCR data (DK);\nSGTF counts adjusted by proportion that was estimated to be Omicron\nbased on GISAID data; data UKHSA, Public Health Scotland,\nStatens Serum Institut, Federal Test Platform Belgium") + # Lesley Scott & NHLS team, 
  ylab("Growth rate advantage of Omicron over Delta\n(difference in growth rate per day)") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(0.98),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char, " per day"), #,
           #                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("Date of diagnosis") + 
  # xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),max(sgtf$date)), ylim=c(0,0.5)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\growth rate advantage omicron logistic fit sgtf data.png"), width=8, height=6)



# 5. ANALYSIS OF ENGLISH DATA AT THE LEVEL OF ONS REGIONS ####

# fit using a separate-slope logistic regression that allows for overdispersion via the inclusion of an observation-level random effect

sgtf_eng_byregion$no_sgtf = sgtf_eng_byregion$pos_tests-sgtf_eng_byregion$sgtf
sgtf_eng_byregion$obs = factor(1:nrow(sgtf_eng_byregion))
sgtf_eng_byregion$date_num = as.numeric(sgtf_eng_byregion$date)
fit_sgtf_eng_byregion1 = glmmPQL(cbind(sgtf, no_sgtf) ~ ns(date_num, df=2, Boundary.knots = c(min(date_num)-0.75,
                                                                                              max(date_num)+0.75)) * region, family=binomial(logit), random=~1|obs, 
                                 data=sgtf_eng_byregion[sgtf_eng_byregion$date>=as.Date("2021-12-01"),]) 
# AIC(fit_sgtf_eng_byregion1)
summary(fit_sgtf_eng_byregion1)

plot(allEffects(fit_sgtf_eng_byregion1, residuals=T))

# growth rate advantage delta_r / s of Omicron over Delta (difference in growth rate per day)
deltar_sgtf_eng_byregion = as.data.frame(emtrends(fit_sgtf_eng_byregion1, ~ date_num, by="region", var="date_num", at=list(date_num=today_num)))[,c(2,3,6,7)] # growth rate advantage of Omicron over Delta
rownames(deltar_sgtf_eng_byregion) = deltar_sgtf_eng_byregion$region
deltar_sgtf_eng_byregion$region = NULL
colnames(deltar_sgtf_eng_byregion) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_eng_byregion
#                         delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# London               0.06859292         0.0444971        0.09268874
# East of England      0.12729800         0.1019584        0.15263760
# South East           0.13560509         0.1103240        0.16088617
# North West           0.18515796         0.1622892        0.20802673
# East Midlands        0.16516529         0.1380261        0.19230445
# South West           0.16510925         0.1300734        0.20014509
# West Midlands        0.16724048         0.1390742        0.19540675
# Yorkshire and Humber 0.16742928         0.1409524        0.19390618
# North East           0.20842761         0.1758435        0.24101176


# mean growth rate advantage of Omicron over Delta across ONS regions (mean difference in growth rate per day)
deltar_sgtf_eng_avgoverregions = as.data.frame(emtrends(fit_sgtf_eng_byregion1, ~ date_num, var="date_num", at=list(date_num=today_num)))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
colnames(deltar_sgtf_eng_avgoverregions) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_eng_avgoverregions
#     delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# 1 0.1544473         0.1452144         0.1636802
deltar_sgtf_eng_avgoverregions_char = sapply(deltar_sgtf_eng_avgoverregions, function (x) sprintf(as.character(round(x,2)), 2) )
deltar_sgtf_eng_avgoverregions_char = paste0(deltar_sgtf_eng_avgoverregions_char[1], " [", deltar_sgtf_eng_avgoverregions_char[2], "-", deltar_sgtf_eng_avgoverregions_char[3],"] 95% CLs")

# transmission advantage of Omicron over Delta (using fixed generation time of 4.7 days)
# i.e. how much higher the effective reproduction number Re value of Omicron is than that of Delta
# and how many more people the virus infects over the course of a single generation
# (due to higher immune escape, i.e. more frequent reinfection of people with immunity built up
# as a result of vaccination or prior exposure)
transmadv_sgtf_eng_byregion = exp(deltar_sgtf_eng_byregion*4.7) 
colnames(transmadv_sgtf_eng_byregion) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf_eng_byregion
#                      transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# London                1.380419            1.232613            1.545948
# East of England       1.819025            1.614789            2.049092
# South East            1.891450            1.679545            2.130092
# North West            2.387490            2.144184            2.658403
# East Midlands         2.173365            1.913096            2.469044
# South West            2.172793            1.842908            2.561728
# West Midlands         2.194667            1.922543            2.505308
# Yorkshire and Humber  2.196615            1.939589            2.487701
# North East            2.663417            2.285224            3.104198



# mean transmission advantage of Omicron over Delta across countries (using generation time of 4.7 days)
# i.e. mean fold difference in the effective reproduction number Re of Omicron compared to that of Delta
transmadv_sgtf_eng_avgoverregions = exp(deltar_sgtf_eng_avgoverregions*4.7)
colnames(transmadv_sgtf_eng_avgoverregions) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf_eng_avgoverregions
#   transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# 1  2.066595            1.978834            2.158248
transmadv_sgtf_eng_avgoverregions_char = sapply(transmadv_sgtf_eng_avgoverregions, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_eng_avgoverregions_char = paste0(transmadv_sgtf_eng_avgoverregions_char[1], " [", transmadv_sgtf_eng_avgoverregions_char[2], "-", transmadv_sgtf_eng_avgoverregions_char[3],"] 95% CLs")

# plot model predictions
emmeans_sgtf_eng_byregion = as.data.frame(emmeans(fit_sgtf_eng_byregion1, ~ date_num+region, at=list(date_num=dateseq, 
                                                                                                     region=unique(sgtf_eng_byregion$region) )), type="response")
colnames(emmeans_sgtf_eng_byregion) = c("date_num", "region", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
emmeans_sgtf_eng_byregion$date = as.Date(emmeans_sgtf_eng_byregion$date_num, origin="1970-01-01")
propomicron_byregion = emmeans_sgtf_eng_byregion[emmeans_sgtf_eng_byregion$date==as.Date("2021-12-07"),] # estimated prop Omicron among confirmed cases on Dec 10
levels_regions = levels(propomicron_byregion$region)[order(propomicron_byregion$prob, decreasing=T)] # sort regions by order of introduction of Omicron
emmeans_sgtf_eng_byregion$region = factor(emmeans_sgtf_eng_byregion$region,
                                                levels=levels_regions)
sgtf_eng_byregion$region = factor(sgtf_eng_byregion$region,
                                        levels=levels_regions)
write.csv(emmeans_sgtf_eng_byregion, ".//data//omicron_sgtf//fit_logistic_mixed_model_share_omicron_ENG_byregion.csv", row.names=F)

# growth rate advantage s of Omicron over Delta (=logistic growth rate/slope=difference in growth rate Omicron & Delta)
emtrends_sgtf_eng_byregion = as.data.frame(emtrends(fit_sgtf_eng_byregion1, ~ date_num, var="date_num", by="region",
                                       at=list(date_num=dateseq, region=unique(sgtf_eng_byregion$region))), type="link")
colnames(emtrends_sgtf_eng_byregion) = c("date_num", "region", "s", "SE", "df", "s.LCL", "s.UCL")
emtrends_sgtf_eng_byregion$date = as.Date(emtrends_sgtf_eng_byregion$date_num, origin="1970-01-01")
emtrends_sgtf_eng_byregion$prop_omicron = emmeans_sgtf_eng_byregion$prob[match(interaction(emtrends_sgtf_eng_byregion$region, emtrends_sgtf_eng_byregion$date),
                                                                               interaction(emmeans_sgtf_eng_byregion$region, emmeans_sgtf_eng_byregion$date))]
emtrends_sgtf_eng_byregion$s[emtrends_sgtf_eng_byregion$prop_omicron<0.01] = NA
emtrends_sgtf_eng_byregion$s.LCL[emtrends_sgtf_eng_byregion$prop_omicron<0.01] = NA
emtrends_sgtf_eng_byregion$s.UCL[emtrends_sgtf_eng_byregion$prop_omicron<0.01] = NA


# plot of share of Omicron among confirmed infections (on logit scale)
qplot(data=sgtf_eng_byregion[sgtf_eng_byregion$date>=as.Date("2021-12-01"),], 
      x=date, y=prop_sgtf, geom="point", colour=region, fill=region, size=pos_tests, alpha=(1/pos_tests), shape=I(16)) + 
  # facet_wrap(~ region) +
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
  annotate("text", x = c(as.Date("2021-12-01")), y = c(0.98),
            label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_eng_avgoverregions_char, " per day"),
  #                         " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_eng_avgoverregions_char),
            color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + 
  # xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),as.Date("2022-01-01")), ylim=c(0.001,0.99))  +
  scale_colour_hue(h=c(10, 310), c=100) +
  scale_fill_hue(h=c(10, 310), c=100) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data england by region.png"), width=8, height=6)

# plot of share of Omicron among confirmed infections (on linear scale)
qplot(data=sgtf_eng_byregion[sgtf_eng_byregion$date>=as.Date("2021-12-01"),], x=date, y=100*prop_sgtf, 
      geom="point", colour=region, fill=region, size=pos_tests, alpha=(1/pos_tests), shape=I(16)) + 
  geom_ribbon(data=emmeans_sgtf_eng_byregion, aes(y=prob*100, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, size=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf_eng_byregion, aes(y=prob*100), alpha=I(0.5), size=1) +
  scale_size_continuous("nr. of positive tests", range=c(2,7), breaks=c(100, 1000, 10000)) +
  scale_alpha_continuous("", range=c(0.5,1), guide=FALSE) +
  xlim(date.from, date.to) +
  # scale_size_continuous(range=c()) +
  ggtitle("Spread of the Omicron variant in England,\nusing S gene target failure (SGTF) data as a proxy",
          "(data UKHSA, mixed logistic 2 df spline fit with observation-level random effect to take into account overdispersion)") +  
  ylab("Share of tests with S dropout (% Omicron)") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(90),
            label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_eng_avgoverregions_char," per day"),
  #                         " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_eng_avgoverregions_char),
            color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") + 
  # xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),as.Date("2022-01-01")), ylim=c(0,100)) +
  scale_colour_hue(h=c(10, 310), c=100) +
  scale_fill_hue(h=c(10, 310), c=100) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))


ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data england by region_linear scale.png"), width=8, height=6)



# plot of growth rate advantage of Omicron over Delta (difference in growth rate per day)
qplot(data=emtrends_sgtf_eng_byregion[emtrends_sgtf_eng_byregion$date>=as.Date("2021-12-01"),], x=date, y=s, ymin=s.LCL, ymax=s.UCL, geom="line", colour=region, fill=region) + 
  geom_ribbon(aes(colour=NULL), alpha=I(0.5)) +
  geom_line(alpha=I(0.5), size=1) +
  scale_size_continuous("nr. of positive tests", range=c(2,7), breaks=c(100, 1000, 10000)) +
  scale_alpha_continuous("", range=c(0.5,1), guide=FALSE) +
  xlim(date.from, date.to) +
  # scale_size_continuous(range=c()) +
  ggtitle("Growth rate advantage of Omicron over Delta in England",
          "inferred from slope of mixed logistic 2 df spline fit to S gene target\nfailure (SGTF) data UKHSA") +
  ylab("Growth rate advantage of Omicron over Delta\n(difference in growth rate per day)") +
  annotate("text", x = c(as.Date("2021-12-01")), y = c(0.98),
           label = paste0("Avg. growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char, " per day"), #,
           #                          " per day\n\nAvg. transmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("Date of diagnosis") + 
  # xaxis +
  theme(plot.subtitle=element_text(size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-12-01"),max(sgtf$date)), ylim=c(0,0.5)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\growth rate advantage omicron logistic fit sgtf data england by region.png"), width=8, height=6)

# TO DO - correlate growth rate advantage with some measured from the UK corona dashboard, e.g. cumulative incidence over the past month etc...
qplot(data=emtrends_sgtf_eng_byregion, x=prop_omicron, y=s, group=region, colour=region, geom="line") +
  facet_wrap(~region) +
  scale_x_continuous(trans="logit") +
  scale_y_continuous(trans="exp") +
  coord_cartesian(xlim=c(0.05, 0.95))

# some plots on the situation in the UK
# daily new cases by ONS region, with estimated share of Omicron inferred from SGTF data :
cases_uk_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDate&format=csv")
cases_uk_region$date = as.Date(cases_uk_region$date)
cases_uk_region$DATE_NUM = as.numeric(cases_uk_region$date)
cases_uk_region$REGION = factor(cases_uk_region$areaName, levels=levels_REGION)
cases_uk_region$areaName = NULL
cases_uk_region = cases_uk_region[cases_uk_region$date<=(max(cases_uk_region$date)-3),] # cut off data from last 3 days (incomplete)
emmeans_sgtf_eng_byregion2 = emmeans_sgtf_eng_byregion
emmeans_sgtf_eng_byregion2$region = factor(emmeans_sgtf_eng_byregion2$region,
                                           levels=c("London","East of England","South East","North West","East Midlands","South West","West Midlands",       
                                                    "Yorkshire and Humber","North East"),
                                           labels=c("London","East of England","South East","North West","East Midlands","South West","West Midlands",       
                                                    "Yorkshire and The Humber","North East"))
cases_uk_region$prop_omicron = emmeans_sgtf_eng_byregion2$prob[match(interaction(cases_uk_region$date, cases_uk_region$REGION),
                                                                    interaction(emmeans_sgtf_eng_byregion2$date, emmeans_sgtf_eng_byregion2$region))] # proportion that was Omicron as inferred from SGTF data
cases_uk_region$prop_omicron[is.na(cases_uk_region$prop_omicron)] = 0
cases_uk_region$newCasesBySpecimenDate_Other = cases_uk_region$newCasesBySpecimenDate*(1-cases_uk_region$prop_omicron)
cases_uk_region$newCasesBySpecimenDate_Omicron = cases_uk_region$newCasesBySpecimenDate*cases_uk_region$prop_omicron
cases_uk_region_long = gather(cases_uk_region, variant, newCasesBySpecimenDate, newCasesBySpecimenDate_Other:newCasesBySpecimenDate_Omicron, factor_key=T)
cases_uk_region_long$variant = factor(cases_uk_region_long$variant,
                                      levels=c("newCasesBySpecimenDate_Omicron","newCasesBySpecimenDate_Other"),
                                      labels=c("Omicron","Other"))

ggplot(data=cases_uk_region_long, 
       aes(x=date, y=newCasesBySpecimenDate+0.01, group=variant, colour=variant, fill=variant)) +
  facet_wrap(~REGION, scale="free_y") +
  xaxis +
  # geom_point(cex=I(0.2)) +
  # geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  # geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 60, fx=T)) + 
  # scale_y_log10() + 
  stat_smooth(se=FALSE, geom="area", lwd=I(1.1), n=1000,
              method = 'glm', formula = y ~ ns(x, df=70), method.args=list(family="poisson"), # span=0.03,
              alpha=1, aes(fill=variant, colour=NULL), position="stack") +
  ylab("New cases (per day)") + xlab("specimen date") +
  ggtitle("DAILY DIAGNOSED SARS-CoV2 CASES IN ENGLAND", "case data gov.uk, https://coronavirus.data.gov.uk/details/download\nwith Omicron share estimated using mixed logistic spline fit to UKHSA SGTF data,\nhttps://www.gov.uk/government/publications/covid-19-omicron-daily-overview\n(note: low case ascertainment during 1st wave)") +
  theme_hc() +
  theme(legend.position="top") +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  scale_fill_manual("", values=c("red3","grey65")) # +
  # scale_y_continuous(trans=log1p_trans()) +
  # coord_trans(y=expm1_trans())

# ggplot(data=cases_uk_region, 
#        aes(x=date, y=newCasesBySpecimenDate, group=REGION, colour=REGION, fill=REGION)) +
#   facet_wrap(~REGION, scale="free_y") +
#   xaxis +
#   # geom_point(cex=I(0.2)) +
#   geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
#   geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 60, fx=T)) + 
#   # scale_y_log10() + 
#   ylab("New cases (per day)") + 
#   ggtitle("DAILY DIAGNOSED SARS-CoV2 CASES IN ENGLAND", "data gov.uk, https://coronavirus.data.gov.uk/details/download\n(note: low case ascertainment during 1st wave)") +
#   theme_hc() +
#   theme(legend.position="none") +
#   labs(tag=tag) +
#   theme(plot.tag.position = "bottomright",
#         plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases England by ONS region by variant.png"), width=10, height=6)




# daily  new cases by ONS region & age group :
cases_uk_age_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateAgeDemographics&format=csv")
cases_uk_age_region$date = as.Date(cases_uk_age_region$date)
cases_uk_age_region$DATE_NUM = as.numeric(cases_uk_age_region$date)
cases_uk_age_region$REGION = factor(cases_uk_age_region$areaName, levels=levels_REGION)
cases_uk_age_region$areaName = NULL
cases_uk_age_region=cases_uk_age_region[cases_uk_age_region$age!="unassigned",]


ggplot(data=cases_uk_age_region[!cases_uk_age_region$age %in% c("00_59","60+"),], 
       aes(x=date, y=cases, group=age, colour=age, fill=age)) +
  xaxis +
  facet_wrap(~REGION, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  # geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 50, fx=T), method.args=list(family="poisson"), n=1000, lwd=I(0.5)) + 
  scale_y_log10() + 
  ylab("diagnosed cases (per day)") + 
  ggtitle("DAILY DIAGNOSED SARS-CoV2 CASES BY AGE IN ENGLAND", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="right") +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases England by age and ONS region.png"), width=10, height=6)


# daily hospital admissions by NHS region
hosps_uk_nhsregion = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=newAdmissions&metric=hospitalCases&format=csv") # &metric=newCasesBySpecimenDate gives NAs
hosps_uk_nhsregion$date = as.Date(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$DATE_NUM = as.numeric(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$REGION = factor(hosps_uk_nhsregion$areaName, levels=levels_NHSREGION)

ggplot(data=hosps_uk_nhsregion, 
       aes(x=date, y=newAdmissions, group=REGION, colour=REGION, fill=REGION)) +
  xaxis +
  facet_wrap(~REGION, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 70, fx=T), method.args=list(family="poisson"), n=1000) + 
  # scale_y_log10() + 
  ylab("New hospital admissions (per day)") + 
  ggtitle("DAILY COVID HOSPITAL ADMISSIONS IN ENGLAND", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="none") +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\hospital admissions England by NHS region.png"), width=10, height=6)

ggplot(data=hosps_uk_nhsregion, 
       aes(x=date, y=hospitalCases, group=REGION, colour=REGION, fill=REGION)) +
  xaxis +
  facet_wrap(~REGION, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 70, fx=T), method.args=list(family="poisson"), n=1000) + 
  # scale_y_log10() + 
  ylab("Hospitalised Covid patients") + 
  ggtitle("HOSPITALISED COVID PATIENTS IN ENGLAND", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="none") +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\hospitalised patients in England by NHS region.png"), width=10, height=6)

# hospitalisations by age :
hosps_eng_age = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=cumAdmissionsByAge&format=csv")
hosps_eng_age$date = as.Date(hosps_eng_age$date)
hosps_eng_age$DATE_NUM = as.numeric(hosps_eng_age$date)

hosps_eng_age = data.frame(hosps_eng_age %>% group_by(age) %>% mutate(Diff = value - lag(value,1,order_by=date))) # see https://stackoverflow.com/questions/26291988/how-to-create-a-lag-variable-within-each-group/26292059
hosps_eng_age$age = factor(hosps_eng_age$age, levels=c("0_to_5", "6_to_17", "18_to_64", "65_to_84", "85+"))
hosps_eng_age$week = cut(hosps_eng_age$date, breaks="week")
hosps_eng_age_byweek = hosps_eng_age %>% group_by(week, age) %>% summarise(hosps=sum(Diff))  
hosps_eng_age_byweek$hosps[is.na(hosps_eng_age_byweek$hosps)] = 0

ggplot(data=hosps_eng_age, 
       aes(x=date, y=Diff, group=age, colour=age, fill=age)) +
  xaxis +
  facet_wrap(~ age, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL), position="identity") +
  # geom_line() +
  # geom_area(position="stack") +
  stat_smooth(se=FALSE, geom="line", lwd=I(1.1),
              method = 'gam', formula = y ~ s(x, k = 60, fx=T), method.args=list(family="poisson"), n=1000,
              alpha=1, aes(fill=age, colour=age), position="identity") +
  # geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 40, fx=T)) + 
  # scale_y_log10() + 
  ylab("new hospitalisations (per day)") + 
  ggtitle("DAILY COVID HOSPITAL ADMISSIONS BY AGE IN ENGLAND", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="none") +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\hospital admissions in England by age.png"), width=10, height=6)





# PLOTS OF NEW CASES PER DAY & Re VALUES BY VARIANT FOR BELGIUM (for TaqPath tests only + also for all Sciensano case data) ####

# for TaqPath tests only ####

source(".//scripts//downloadData.R") # download Belgian case data
head(sgtf)
sgtf_long = gather(sgtf, variant, count, all_of(c("omicron","non_omicron")), factor_key=TRUE)
sgtf_long$country
sgtf_long$date_num = as.numeric(sgtf_long$date)
sgtf_long$WEEKDAY = weekdays(sgtf_long$date)
sgtf_long$BANKHOLIDAY = bankholiday(sgtf_long$date)
sgtf_long$variant = factor(sgtf_long$variant, levels=c("omicron", "non_omicron"), labels=c("Omicron","Delta"))


# smooth out weekday effects in variant case nrs using negative binomial GAM
library(mgcv)
k=4
fit_cases = gam(count ~ s(date_num, bs="cs", k=k, m=c(2), fx=F, by=variant) + variant + WEEKDAY + BANKHOLIDAY,
                family=nb(), 
                data=sgtf_long[sgtf_long$country=="Belgium",],
                method = "REML",
                knots = list(date_num = c(min(sgtf_long$date_num)-14,
                                          seq(min(sgtf_long$date_num)+0.5*diff(range(sgtf_long$date_num))/(k-2), 
                                              max(sgtf_long$date_num)-0.5*diff(range(sgtf_long$date_num))/(k-2), length.out=k-2),
                                          max(sgtf_long$date_num)+14))
) 
BIC(fit_cases)

emmeans_casesbyvariant = as.data.frame(emmeans(fit_cases, ~ date_num+variant, at=list(date_num=as.numeric(seq(as.Date("2021-11-27"),
                                                                                                              as.Date("2022-01-18"),
                                                                                                   by=1)),
                                                                                      BANKHOLIDAY=factor("no")), type="response"))
emmeans_casesbyvariant$date = as.Date(emmeans_casesbyvariant$date_num, origin="1970-01-01")
emmeans_casesbyvariant$count= emmeans_casesbyvariant$response
lineage_cols = c("red2","grey50")

# calculate growth rates & Re values per variant, using variant-specific generation time
growth_rates_variants = as.data.frame(emtrends(fit_cases, ~ date_num, var="date_num", by="variant",
                              at=list(date_num=seq(as.numeric(as.Date("2021-11-27")),
                                                   today_num, # as.numeric(as.Date(max(sgtf_be$date)))
                                                   by=1),
                                      BANKHOLIDAY=factor("no")
                                      # TESTS_ALL=max(cases_tot$TESTS_ALL)
                                      ),
                              type="link"))
colnames(growth_rates_variants)[3] = "r"
colnames(growth_rates_variants)[6] = "r_LOWER"
colnames(growth_rates_variants)[7] = "r_UPPER"
growth_rates_variants$date = as.Date(growth_rates_variants$date_num, origin="1970-01-01") - 7 # -7 to get time of infection
growth_rates_variants$Re = Re.from.r.variant(growth_rates_variants$r, growth_rates_variants$variant)
growth_rates_variants$Re_LOWER = Re.from.r.variant(growth_rates_variants$r_LOWER, growth_rates_variants$variant)
growth_rates_variants$Re_UPPER = Re.from.r.variant(growth_rates_variants$r_UPPER, growth_rates_variants$variant)

growth_rates_variants_wide = spread(growth_rates_variants[,c("date_num","variant","Re")], variant, Re)
colnames(growth_rates_variants_wide) = c("date_num","Re_Omicron","Re_Delta")
growth_rates_variants_wide_lower = spread(growth_rates_variants[,c("date_num","variant","Re_LOWER")], variant, Re_LOWER)
colnames(growth_rates_variants_wide_lower) = c("date_num","Re_Omicron_LOWER","Re_Delta_LOWER")
growth_rates_variants_wide_lower$date_num = NULL
growth_rates_variants_wide_upper = spread(growth_rates_variants[,c("date_num","variant","Re_UPPER")], variant, Re_UPPER)
colnames(growth_rates_variants_wide_upper) = c("date_num","Re_Omicron_UPPER","Re_Delta_UPPER")
growth_rates_variants_wide_upper$date_num = NULL
growth_rates_variants_wide = cbind(growth_rates_variants_wide, growth_rates_variants_wide_lower, growth_rates_variants_wide_upper)
growth_rates_variants_wide$transmadv_Omicron = growth_rates_variants_wide$Re_Omicron/growth_rates_variants_wide$Re_Delta
growth_rates_variants_wide$date = as.Date(growth_rates_variants_wide$date_num, origin="1970-01-01") - 7 # -7 to get time of infection

# plot of Re values per variant
ggplot(data=growth_rates_variants, 
       aes(x=date, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, colour=variant, fill=variant, group=variant)) + 
  geom_ribbon(aes(colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2))) +
  # xaxis + 
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("Rt value") + xlab("Estimated date of infection\n(7 days before diagnosis)") +
  ggtitle("Rt VALUES OF OMICRON & DELTA IN BELGIUM","data Federal Test Platform Belgium,\nbased on negative binomial GAM fit to S dropout\n& S positive counts, adjusted for proportion\nof S dropouts that were Omicron (GISAID data)\nusing gamma distributed generation time with mean and SD of\n2.2 and 1.6 days for Omicron (Kim et al. 2022) &\n3.9 and 2.8 days for Delta (Hart et al. 2022) ") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  geom_hline(yintercept=1, colour=I("red")) +
  coord_cartesian(ylim=c(0.6, 1.8)) +
  scale_y_continuous(breaks=seq(0.4,1.8,by=0.2), trans="log2") +
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  # theme(axis.title.x=element_blank()) +
  theme(plot.subtitle=element_text(size=10)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\Rt by variant_belgium.png"), width=8, height=6)

# plot transmission advantage (Re Omicron / Re Delta) over time
ggplot(data=growth_rates_variants_wide[growth_rates_variants_wide$date>=as.Date("2021-12-06"),], 
       aes(x=date, y=transmadv_Omicron, 
           # ymin=Re_LOWER, ymax=Re_UPPER
           )) + 
  # geom_ribbon(aes(colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2))) +
  # xaxis + 
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("Transmission advantage (Rt Omicron / Rt Delta)") + xlab("Estimated date of infection") +
  ggtitle("TRANSMISSION ADVANTAGE OF OMICRON IN BELGIUM","data Federal Test Platform Belgium,\nbased on negative binomial GAM fit to S dropout\n& S positive counts, adjusted for proportion\nof S dropouts that were Omicron (GISAID data)\nusing gamma distributed generation time with mean and SD of\n2.2 and 1.6 days for Omicron (Kim et al. 2022) &\n3.9 and 2.8 days for Delta (Hart et al. 2022) ") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  # geom_hline(yintercept=1, colour=I("red")) +
  coord_cartesian(ylim=c(1, 2.5)) +
  scale_y_continuous(breaks=seq(0.4,2.5,by=0.2)) +
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  # theme(axis.title.x=element_blank()) +
  theme(plot.subtitle=element_text(size=10)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\transmission advantage omicron_belgium.png"), width=8, height=6)



# more plots

ggplot(data=sgtf_long[sgtf_long$country=="Belgium",], 
       aes(x=date, y=count, colour=variant, group=variant)) + 
  geom_point() +
  geom_ribbon(data=emmeans_casesbyvariant, aes(ymin=lower.CL, ymax=upper.CL, fill=variant, colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  geom_line(data=emmeans_casesbyvariant, aes(lwd=I(1.2))) +
  # xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("") +
  ggtitle("RISE OF OMICRON AMONG TaqPath DIAGNOSED\nPOSITIVE SAMPLES IN BELGIUM","data Federal Test Platform Belgium,\nbased on negative binomial GAM fit to S dropout\n& S positive counts, adjusted for proportion\nof S dropouts that were Omicron (GISAID data)") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  annotate("rect", xmin=max(sgtf_long[sgtf_long$country=="Belgium","date"])+1, 
           xmax=I(as.Date("2022-01-18")+2), ymin=0, ymax=80000, alpha=0.4, fill="white") + # extrapolated part
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10() +
  theme(axis.title.x=element_blank()) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\Taqpath cases per day by variant_belgium_log10 y scale.png"), width=8, height=6)

ggplot(data=sgtf_long[sgtf_long$country=="Belgium",], 
       aes(x=date, y=count, colour=variant, group=variant)) + 
  geom_point() +
  geom_ribbon(data=emmeans_casesbyvariant, aes(ymin=lower.CL, ymax=upper.CL, fill=variant, colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  geom_line(data=emmeans_casesbyvariant, aes(lwd=I(1.2))) +
  # xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("") +
  ggtitle("RISE OF OMICRON AMONG TaqPath DIAGNOSED\nPOSITIVE SAMPLES IN BELGIUM","data Federal Test Platform Belgium,\nbased on negative binomial GAM fit to S dropout\n& S positive counts, adjusted for proportion\nof S dropouts that were Omicron (GISAID data)") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  annotate("rect", xmin=max(sgtf_long[sgtf_long$country=="Belgium","date"])+1, 
           xmax=I(as.Date("2022-01-18")+2), ymin=0, ymax=80000, alpha=0.4, fill="white") + # extrapolated part
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  theme(axis.title.x=element_blank()) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\Taqpath cases per day by variant_belgium_linear y scale.png"), width=8, height=6)

ggplot(data=sgtf_long[sgtf_long$country=="Belgium",], 
       aes(x=date, y=count, colour=variant, group=variant)) + 
  # geom_point() +
  geom_area(data=emmeans_casesbyvariant, aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  # geom_ribbon(data=emmeans_casesbyvariant, aes(ymin=lower.CL, ymax=upper.CL, fill=variant, colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  # geom_line(data=emmeans_casesbyvariant, aes(lwd=I(1.2))) +
  # xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("") +
  ggtitle("RISE OF OMICRON AMONG TaqPath DIAGNOSED\nPOSITIVE SAMPLES IN BELGIUM","data Federal Test Platform Belgium,\nbased on negative binomial GAM fit to S dropout\n& S positive counts, adjusted for proportion\nof S dropouts that were Omicron (GISAID data)") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  annotate("rect", xmin=max(sgtf_long[sgtf_long$country=="Belgium","date"])+1, 
           xmax=I(as.Date("2022-01-18")), ymin=0, ymax=60000, alpha=0.4, fill="white") + # extrapolated part
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  theme(axis.title.x=element_blank()) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\Taqpath cases per day by variant_belgium_stacked area chart.png"), width=8, height=6)

ggplot(data=sgtf_long[sgtf_long$country=="Belgium",], 
       aes(x=date, y=count+1, colour=variant, group=variant)) + 
  # geom_point() +
  geom_area(data=emmeans_casesbyvariant[emmeans_casesbyvariant$variant=="Omicron",], aes(lwd=I(1.2), colour=NULL, fill=lineage_cols[2], group=variant, alpha=I(0.3)), position="stack") +
  geom_area(data=emmeans_casesbyvariant[emmeans_casesbyvariant$variant=="Delta",], aes(lwd=I(1.2), colour=NULL, fill=lineage_cols[1], group=variant, alpha=I(0.3)), position="stack") +
  # geom_ribbon(data=emmeans_casesbyvariant, aes(ymin=lower.CL, ymax=upper.CL, fill=variant, colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  # geom_line(data=emmeans_casesbyvariant, aes(lwd=I(1.2))) +
  # xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("") +
  ggtitle("RISE OF OMICRON AMONG TaqPath DIAGNOSED\nPOSITIVE SAMPLES IN BELGIUM","data Federal Test Platform Belgium,\nbased on negative binomial GAM fit to S dropout\n& S positive counts, adjusted for proportion\nof S dropouts that were Omicron (GISAID data)") +
  scale_fill_manual("", values=lineage_cols, labels=c("Omicron", "Delta")) +
  scale_colour_manual("", values=lineage_cols, labels=c("Omicron", "Delta")) +
  annotate("rect", xmin=max(sgtf_long[sgtf_long$country=="Belgium","date"])+1, 
           xmax=I(as.Date("2022-01-18")), ymin=0, ymax=60000, alpha=0.4, fill="white") + # extrapolated part
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  theme(axis.title.x=element_blank()) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\Taqpath cases per day by variant_belgium_overlaid area chart.png"), width=8, height=6)








# for all cases (data Sciensano) ####

cases_tot$prop_omicron = emmeans_sgtf$prob[match(interaction("Belgium", cases_tot$DATE_NUM),
                                                 interaction(emmeans_sgtf$country, emmeans_sgtf$date_num))] # prop omicron, from logistic spline fit to SGTF data
cases_tot$prop_omicron[is.na(cases_tot$prop_omicron)] = 0
cases_tot$CASES_Omicron = cases_tot$CASES*cases_tot$prop_omicron
cases_tot$CASES_Delta = cases_tot$CASES*(1-cases_tot$prop_omicron)
cases_tot2 = cases_tot[!is.na(cases_tot$CASES_Omicron),]

cases_tot_long = gather(cases_tot2[,c("DATE","CASES_Omicron","CASES_Delta")], variant, count, all_of(c("CASES_Omicron","CASES_Delta")), factor_key=TRUE)
cases_tot_long$variant = factor(cases_tot_long$variant, levels=c("CASES_Omicron", "CASES_Delta"), labels=c("Omicron","Delta"))
cases_tot_long$date_num = as.numeric(cases_tot_long$DATE)
cases_tot_long$WEEKDAY = weekdays(cases_tot_long$DATE)
cases_tot_long$BANKHOLIDAY = bankholiday(cases_tot_long$DATE)
cases_tot_long$TESTS_ALL = cases_tot$TESTS_ALL[match(cases_tot_long$DATE, cases_tot$DATE)]

# growth rate advantage of Omicron over Delta s over time in Belgium
s = as.data.frame(emtrends(fit_sgtf2, ~ date_num, by="country", var="date_num", at=list(date_num=seq(min(cases_tot_long$date_num),
                                                                                                     max(cases_tot_long$date_num),
                                                                                                     by=1),
                                                                                        country="Belgium") )) # growth rate advantage s of Omicron over Delta
s$DATE = as.Date(s$date_num, origin="1970-01-01") 
cases_tot_long$s = s$date_num.trend[match(cases_tot_long$DATE, s$DATE)]
range(cases_tot_long$DATE) # "2020-03-01" "2022-01-04"

# negative binomial GAM fit to Sciensano data (estimated count of new cases by each variant)
k=45
fit_cases_all = gam(count ~ s(date_num, bs = "cs", k=k, m = c(2), fx=F, by=variant) + variant +
                  WEEKDAY + BANKHOLIDAY +
                  s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=nb(link="log"), data=cases_tot_long,
                method = "REML",
                knots = list(date_num = c(min(cases_tot_long$date_num)-10,
                                          seq(min(cases_tot_long$date_num)+
                                                0.75*diff(range(cases_tot_long$date_num))/(k-2), max(cases_tot_long$date_num)-
                                                0.75*diff(range(cases_tot_long$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot_long$date_num)+10)) # see https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/
)
BIC(fit_cases_all) # 2358.354

emmeans_casesbyvariant_all = as.data.frame(emmeans(fit_cases_all, ~ date_num+variant, at=list(date_num=as.numeric(seq(min(cases_tot_long$date_num),
                                                                                                                      max(cases_tot_long$date_num)+14,
                                                                                                              by=1)),
                                                                                      BANKHOLIDAY=factor("no"),
                                                                                      TESTS_ALL=max(cases_tot$TESTS_ALL)), type="response"))
emmeans_casesbyvariant_all$date = as.Date(emmeans_casesbyvariant_all$date_num, origin="1970-01-01")
emmeans_casesbyvariant_all$count= emmeans_casesbyvariant_all$response
lineage_cols = c("red2","grey60")

# calculate growth rates & Re values per variant, using variant-specific generation time
growth_rates_variants_all = as.data.frame(emtrends(fit_cases_all, ~ date_num, var="date_num", by="variant",
                                               at=list(date_num=seq(min(cases_tot_long$date_num),
                                                                    max(cases_tot_long$date_num)+7,
                                                                    by=1),
                                                       BANKHOLIDAY=factor("no"),
                                                       TESTS_ALL=max(cases_tot$TESTS_ALL)
                                               ),
                                               type="link"))
colnames(growth_rates_variants_all)[3] = "r"
colnames(growth_rates_variants_all)[6] = "r_LOWER"
colnames(growth_rates_variants_all)[7] = "r_UPPER"
growth_rates_variants_all$date = as.Date(growth_rates_variants_all$date_num, origin="1970-01-01") - 7 # -7 to get time of infection
growth_rates_variants_all$Re = Re.from.r.variant(growth_rates_variants_all$r, growth_rates_variants_all$variant)
growth_rates_variants_all$Re_LOWER = Re.from.r.variant(growth_rates_variants_all$r_LOWER, growth_rates_variants_all$variant)
growth_rates_variants_all$Re_UPPER = Re.from.r.variant(growth_rates_variants_all$r_UPPER, growth_rates_variants_all$variant)
growth_rates_variants_all$prop_omicron = emmeans_sgtf$prob[match(interaction("Belgium", growth_rates_variants_all$date_num),
                                                                 interaction(emmeans_sgtf$country, emmeans_sgtf$date_num))] # prop omicron, from logistic spline fit to SGTF data
growth_rates_variants_all$Re[growth_rates_variants_all$variant=="Omicron"][growth_rates_variants_all$prop_omicron[growth_rates_variants_all$variant=="Omicron"]<0.01] = NA
growth_rates_variants_all$Re_LOWER[growth_rates_variants_all$variant=="Omicron"][growth_rates_variants_all$prop_omicron[growth_rates_variants_all$variant=="Omicron"]<0.01] = NA
growth_rates_variants_all$Re_UPPER[growth_rates_variants_all$variant=="Omicron"][growth_rates_variants_all$prop_omicron[growth_rates_variants_all$variant=="Omicron"]<0.01] = NA

# TO DO: add avg Re

growth_rates_variants_all_wide = spread(growth_rates_variants_all[,c("date_num","variant","Re")], variant, Re)
colnames(growth_rates_variants_all_wide) = c("date_num","Re_Omicron","Re_Delta")
growth_rates_variants_all_wide_lower = spread(growth_rates_variants_all[,c("date_num","variant","Re_LOWER")], variant, Re_LOWER)
colnames(growth_rates_variants_all_wide_lower) = c("date_num","Re_Omicron_LOWER","Re_Delta_LOWER")
growth_rates_variants_all_wide_lower$date_num = NULL
growth_rates_variants_all_wide_upper = spread(growth_rates_variants_all[,c("date_num","variant","Re_UPPER")], variant, Re_UPPER)
colnames(growth_rates_variants_all_wide_upper) = c("date_num","Re_Omicron_UPPER","Re_Delta_UPPER")
growth_rates_variants_all_wide_upper$date_num = NULL
growth_rates_variants_all_wide = cbind(growth_rates_variants_all_wide, growth_rates_variants_all_wide_lower, growth_rates_variants_all_wide_upper)
growth_rates_variants_all_wide$transmadv_Omicron = growth_rates_variants_all_wide$Re_Omicron/growth_rates_variants_all_wide$Re_Delta
growth_rates_variants_all_wide$date = as.Date(growth_rates_variants_all_wide$date_num, origin="1970-01-01") - 7 # -7 to get time of infection

# plot of Re values per variant
ggplot(data=growth_rates_variants_all, 
       aes(x=date, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, colour=variant, fill=variant, group=variant)) + 
  geom_ribbon(aes(colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2))) +
  # xaxis + 
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("Rt value") + xlab("Estimated date of infection\n(7 days before diagnosis)") +
  ggtitle("Rt VALUES OF OMICRON & DELTA IN BELGIUM","based on negative binomial GAM fit to estimated Omicron & Delta cases,\nusing Sciensano case data and share of Omicron estimated from SGTF\ndata Federal Test Platform, adjusted for proportion of S dropouts that\nwere Omicron (GISAID data), using gamma distributed generation time\nwith mean and SD of 2.2 and 1.6 days for Omicron (Kim et al. 2022) &\n3.9 and 2.8 days for Delta (Hart et al. 2022)") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  geom_hline(yintercept=1, colour=I("red")) +
  coord_cartesian(ylim=c(0.6, 1.8), xlim=c(as.Date("2021-08-01"),as.Date(max(sgtf_be$date))-7)) +
  scale_y_continuous(breaks=seq(0.4,2,by=0.2), trans="log2") +
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  # theme(axis.title.x=element_blank()) +
  theme(plot.subtitle=element_text(size=10)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\Rt by variant_all cases_belgium.png"), width=8, height=6)

# plot transmission advantage (Re Omicron / Re Delta) over time
ggplot(data=growth_rates_variants_all_wide, 
       aes(x=date, y=transmadv_Omicron, 
           # ymin=Re_LOWER, ymax=Re_UPPER
       )) + 
  # geom_ribbon(aes(colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2))) +
  # xaxis + 
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("Transmission advantage (Rt Omicron / Rt Delta)") + xlab("Estimated date of infection\n(7 days before diagnosis)") +
  ggtitle("TRANSMISSION ADVANTAGE OF OMICRON IN BELGIUM","based on negative binomial GAM fit to estimated Omicron & Delta cases,\nusing Sciensano case data and share of Omicron estimated from SGTF\ndata Federal Test Platform, adjusted for proportion of S dropouts that\nwere Omicron (GISAID data), using gamma distributed generation time\nwith mean and SD of 2.2 and 1.6 days for Omicron (Kim et al. 2022) &\n3.9 and 2.8 days for Delta (Hart et al. 2022)") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  # geom_hline(yintercept=1, colour=I("red")) +
  coord_cartesian(ylim=c(1, 2.5), xlim=c(as.Date("2021-12-01"),as.Date(max(sgtf_be$date))-7)) +
  scale_y_continuous(breaks=seq(0.4,2.5,by=0.2)) +
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  # theme(axis.title.x=element_blank()) +
  theme(plot.subtitle=element_text(size=10)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\transmission advantage omicron_all cases_belgium.png"), width=8, height=6)



# more plots

ggplot(data=cases_tot_long, 
       aes(x=DATE, y=count, colour=variant, group=variant)) + 
  geom_point() +
  geom_ribbon(data=emmeans_casesbyvariant_all, aes(x=date, ymin=lower.CL, ymax=upper.CL, fill=variant, colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  geom_line(data=emmeans_casesbyvariant_all, aes(x=date, lwd=I(1.2))) +
  # xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("") +
  ggtitle("DAILY DIAGNOSED SARS-CoV2 CASES IN BELGIUM", "case data Sciensano with Omicron share estimated using mixed logistic spline fit\nto SGTF data Federal Test Platform, with correction for proportion\nof S dropouts that were Omicron based on GISAID data\nand cases by each variant extrapolated using a negative binomial GAM fit\nwith correction for weekday & bank holiday effects and variable testing intensity") +
  scale_fill_manual("", values=lineage_cols ) +
  scale_colour_manual("", values=lineage_cols ) +
  annotate("rect", xmin=max(sgtf_long[sgtf_long$country=="Belgium","date"])+1, 
           xmax=I(as.Date("2022-01-18")), ymin=0, ymax=120000, alpha=0.2, fill="white") + # extrapolated part
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10() +
  theme(axis.title.x=element_blank()) +
  coord_cartesian(xlim=c(as.Date("2021-11-27"),NA), ylim=c(100,120000)) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day by variant_all cases_belgium_log10 y scale.png"), width=8, height=6)

# ggplot(data=cases_tot_long, 
#        aes(x=DATE, y=count, colour=variant, group=variant)) + 
#   geom_point() +
#   geom_ribbon(data=emmeans_casesbyvariant_all, aes(x=date, ymin=lower.CL, ymax=upper.CL, fill=variant, colour=NULL), alpha=I(0.4)) +
#   # facet_wrap(~ country, scale="free") +
#   geom_line(data=emmeans_casesbyvariant_all, aes(x=date, lwd=I(1.2))) +
#   # xaxis +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right") + 
#   ylab("New confirmed cases per day") + xlab("") +
#   ggtitle("RISE OF OMICRON AMONG TaqPath DIAGNOSED\nPOSITIVE SAMPLES IN BELGIUM","data Federal Test Platform Belgium,\nbased on negative binomial GAM fit to S dropout\n& S positive counts, adjusted for proportion\nof S dropouts that were Omicron (GISAID data)") +
#   scale_fill_manual("", values=lineage_cols ) +
#   scale_colour_manual("", values=lineage_cols ) +
#   annotate("rect", xmin=max(sgtf_long[sgtf_long$country=="Belgium","date"])+1, 
#            xmax=I(as.Date("2022-01-18")+4), ymin=0, ymax=120000, alpha=0.4, fill="white") + # extrapolated part
#   # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
#   # scale_y_log10() +
#   theme(axis.title.x=element_blank()) +
#   labs(tag=tag) +
#   theme(plot.tag.position = "bottomright",
#         plot.tag = element_text(vjust = 1, hjust = 1, size=8))
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day by variant_all cases_belgium_linear y scale.png"), width=8, height=6)

emmeans_casesbyvariant_all$variant2 = factor(emmeans_casesbyvariant_all$variant, levels=c("Omicron", "Delta"), labels=c("Omicron", "Other"))
ggplot(data=emmeans_casesbyvariant_all, 
       aes(x=date, y=count, colour=variant2, group=variant2)) + 
  # geom_point() +
  geom_area(data=emmeans_casesbyvariant_all, aes(x=date, lwd=I(1.2), colour=NULL, fill=variant2, group=variant2), position="stack") +
  # geom_ribbon(data=emmeans_casesbyvariant, aes(ymin=lower.CL, ymax=upper.CL, fill=variant, colour=NULL), alpha=I(0.4)) +
  # facet_wrap(~ country, scale="free") +
  # geom_line(data=emmeans_casesbyvariant, aes(lwd=I(1.2))) +
  # xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("") +
  ggtitle("DAILY DIAGNOSED SARS-CoV2 CASES IN BELGIUM", "case data Sciensano with Omicron share estimated using mixed logistic spline fit\nto SGTF data Federal Test Platform, with correction for proportion\nof S dropouts that were Omicron based on GISAID data\nand cases by each variant extrapolated using a negative binomial GAM fit\nwith correction for weekday & bank holiday effects and variable testing intensity") +
  scale_fill_manual("", values=c("red2","grey65") ) +
  scale_colour_manual("", values=c("red2","grey65") ) +
  annotate("rect", xmin=max(sgtf_long[sgtf_long$country=="Belgium","date"])+1, 
           xmax=I(as.Date("2022-01-18")), ymin=0, ymax=80000, alpha=0.4, fill="white") + # extrapolated part
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  # scale_y_log10() +
  theme(axis.title.x=element_blank()) +
  labs(tag=tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day by variant_all cases_belgium_stacked area chart.png"), width=8, height=6)























# CODE BELOW NOT FINISHED YET - DISCARD

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



