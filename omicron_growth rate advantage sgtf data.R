# ANALYSIS OF GROWTH RATE ADVANTAGE OF OMICRON IN South Africa, England, Denmark & Belgium BASED ON SGTF (S gene target failure / S dropout) DATA
# data: 
# Belgium: Federal Test Platform, provided by Emmanuel André
# Denmark: Statens Serum Institut (https://covid19.ssi.dk/virusvarianter/delta-pcr)
# England: UK Health Security Agency (https://www.gov.uk/government/publications/covid-19-omicron-daily-overview)
# South Africa: Lesley Scott & NHLS team (traced off graph by Alex Selby)
# GISAID data to infer the proportion of S dropout (with spike Spike_H69del/Spike_V70del) that are Omicron (https://www.gisaid.org/)

# T. Wenseleers
# last update 19 DECEMBER 2021

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
  geom_ribbon(aes(colour=NULL), alpha=I(0.2)) +
  geom_line() +
  geom_point(data=GISAID_sgtf_isomicron_subs, aes(y=prop_omicron, ymin=NULL, ymax=NULL, size=total_seqs_sgtf)) +
  scale_size_continuous("nr. of GISAID sequences\nwith SGTF", range=c(1,4), breaks=c(1, 10, 50)) +
  facet_wrap(~ country) +
  theme_hc() +
  ylab("Proportion of SGTF sequences that are Omicron") +
  ggtitle("Proportion of GISAID sequences with SGTF\n(with Spike_H69del/Spike_V70del)\nthat are Omicron", "(GISAID data, logistic fit)") +
  xaxis + theme(legend.position="right")

ggsave(file=paste0(".\\plots\\",plotdir,"\\estimated prop sgtf seqs that are Omicron from GISAID data.png"), width=8, height=6)

# 3. Import variant-specific PCR data or SGTF data (now proxy for Omicron)
varpcr_dk = read.csv(".//data//omicron_sgtf//variantpcr_denmark.csv") # variant-specfc PCR data for DK
# sgtf_dk = read.csv(".//data//omicron_sgtf//sgtf_denmark.csv") 
sgtf_sa = read.csv(".//data//omicron_sgtf//sgtf_south africa.csv") # data traced from graph (data not open at the moment)
sgtf_sa$date = as.Date(sgtf_sa$date)
sgtf_sa = sgtf_sa[sgtf_sa$date>=as.Date("2021-11-05"),] # data limited to >= Nov 1 when from sequencing data nearly all S dropout was Omicron
sgtf_eng = read.csv("https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1042224/sgtf_regionepicurve_2021-12-17.csv")
sgtf_eng$sgtf = factor(sgtf_eng$sgtf, levels=c("Cases with confirmed S-gene","Cases with confirmed SGTF"), labels=c("no_sgtf","sgtf"))
# data for England by LTLA are here 
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1041668/sgtf_LAdaily_2021-12-15__1_.csv
# see https://www.gov.uk/government/publications/covid-19-omicron-daily-overview
sgtf_eng = spread(sgtf_eng, sgtf, n)
sgtf_eng[is.na(sgtf_eng)] = 0
sgtf_eng$date = dmy(sgtf_eng$specimen_date)
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
sgtf_scot = read.csv(".//data//omicron_sgtf//sgtf_scotland.csv") 
sgtf_be = read.csv(".//data//omicron_sgtf//sgtf_belgium.csv") 

# TO DO still add SGTF data from San Diego (Andersen Laboratory @ Scripps Research): 
# https://raw.githubusercontent.com/andersen-lab/SARS-CoV-2_SGTF_San-Diego/main/SGTF_San_Diego.csv
# https://github.com/andersen-lab/SARS-CoV-2_SGTF_San-Diego

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

write.csv(sgtf, ".//data//omicron_sgtf//raw_data_share_omicron_SA_ENG_SCOT_DK_BE.csv", row.names=F)

# 4. ANALYSIS OF SGTF DATA USING LOGISTIC REGRESSION ####

# fit using regular binomial GLM (not taking into account temporal autocorrelation or overdispersion)
# 
# fit_sgtf0 = glm(cbind(omicron, non_omicron) ~ scale(date_num) + country, family=binomial(logit), data=sgtf) 
# fit_sgtf1 = glm(cbind(omicron, non_omicron) ~ scale(date_num) * country, family=binomial(logit), data=sgtf) 
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf1)

# fit using a separate-slope logistic regression that allows for overdispersion via the inclusion of an observation-level random effect

fit_sgtf1 = glmmPQL(cbind(omicron, non_omicron) ~ scale(date_num) * country, family=binomial(logit), random=~1|obs, data=sgtf) 
# AIC(fit_sgtf0, fit_sgtf1)
summary(fit_sgtf1)

# note allowing lag-1 temporally autocorrelated residuals would be possible using
# fit_sgtf0 = glmmPQL(cbind(omicron, non_omicron) ~ scale(date_num) + country, family=binomial(logit), random=~1|obs, correlation=corAR1(), data=sgtf) 
# fit_sgtf1 = glmmPQL(cbind(omicron, non_omicron) ~ scale(date_num) * country, family=binomial(logit), random=~1|obs, correlation=corAR1(), data=sgtf) 
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf1)

plot(allEffects(fit_sgtf1, residuals=T))

# growth rate advantage delta_r of Omicron over Delta (difference in growth rate per day)
deltar_sgtf_bycountry = as.data.frame(emtrends(fit_sgtf1, ~ date_num, by="country", var="date_num", at=list(date_num=max(dateseq))))[,c(2,3,6,7)] # growth rate advantage of Omicron over Delta
rownames(deltar_sgtf_bycountry) = deltar_sgtf_bycountry$country
deltar_sgtf_bycountry$country = NULL
colnames(deltar_sgtf_bycountry) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf_bycountry
#                delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# South Africa 0.2546775         0.2000379         0.3093170
# England      0.3605907         0.3450973         0.3760841
# Scotland     0.2838093         0.2580729         0.3095457
# Denmark      0.3710232         0.3372933         0.4047530
# Belgium      0.2758321         0.2442708         0.3073934

# mean growth rate advantage of Omicron over Delta across countries (mean difference in growth rate per day)
deltar_sgtf = as.data.frame(emtrends(fit_sgtf1, ~ date_num, var="date_num", at=list(date_num=max(dateseq))))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
colnames(deltar_sgtf) = c("delta_r","delta_r_asymp.LCL","delta_r_asymp.UCL")
deltar_sgtf
#     delta_r delta_r_asymp.LCL delta_r_asymp.UCL
# 1 0.3091865         0.2936666         0.3247065
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
# South Africa  3.310119            2.560438            4.279302
# England       5.445428            5.062991            5.856752
# Scotland      3.795832            3.363367            4.283903
# Denmark       5.719085            4.880652            6.701551
# Belgium       3.656151            3.152113            4.240787

# mean transmission advantage of Omicron over Delta across countries (using generation time of 4.7 days)
# i.e. mean fold difference in the effective reproduction number Re of Omicron compared to that of Delta
transmadv_sgtf = exp(deltar_sgtf*4.7)
colnames(transmadv_sgtf) = c("transmadv","transmadv_asymp.LCL","transmadv_asymp.UCL")
transmadv_sgtf
#   transmadv transmadv_asymp.LCL transmadv_asymp.UCL
# 1  4.276679            3.975828            4.600295
transmadv_sgtf_char = sapply(transmadv_sgtf, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_char = paste0(transmadv_sgtf_char[1], " [", transmadv_sgtf_char[2], "-", transmadv_sgtf_char[3],"] 95% CLs")

# plot model predictions
emmeans_sgtf = as.data.frame(emmeans(fit_sgtf1, ~ date_num+country, at=list(date_num=dateseq, country=levels_country)), type="response")
colnames(emmeans_sgtf) = c("date_num", "country", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
emmeans_sgtf$date = as.Date(emmeans_sgtf$date_num, origin="1970-01-01")

write.csv(emmeans_sgtf, ".//data//omicron_sgtf//fit_logistic_mixed_model_share_omicron_SA_ENG_SCOT_DK_BE.csv", row.names=F)

# plot of share of Omicron among confirmed infections (on logit scale)
qplot(data=sgtf, x=date, y=prop, geom="point", colour=country, fill=country, size=pos_tests) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, size=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf, aes(y=prob), alpha=I(0.5), size=1) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_size_continuous("nr. of positive tests", range=c(1,4), breaks=c(100, 1000, 10000)) +
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
qplot(data=sgtf, x=date, y=100*prop, geom="point", colour=country, fill=country, size=pos_tests) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob*100, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, size=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf, aes(y=prob*100), alpha=I(0.5), size=1) +
  scale_size_continuous("nr. of positive tests", range=c(1,4), breaks=c(100, 1000, 10000)) +
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


# 5.  

