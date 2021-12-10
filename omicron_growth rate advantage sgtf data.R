# ANALYSIS OF GROWTH RATE ADVANTAGE OF OMICRON (B.1.1.529) IN South Africa, England, Denmark & Belgium BASED ON SGTF (S gene target failure / S dropout) DATA
# data: 
# Belgium: Federal Test Platform, provided by Emmanuel André
# Denmark: Statens Serum Institut
# England: UK Health Security Agency (traced off graph by Alex Selby)
# South Africa: Lesley Scott & NHLS team (traced off graph by Alex Selby)

# T. Wenseleers
# last update 10 DECEMBER 2021

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

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

# import SGTF data (now proxy for B.1.1.529 / Omicron)
sgtf_sa = read.csv(".//data//omicron_sgtf//sgtf_south africa.csv")
sgtf_sa$date = as.Date(sgtf_sa$date)
sgtf_sa = sgtf_sa[sgtf_sa$date>=as.Date("2021-11-01"),] 
sgtf_eng = read.csv(".//data//omicron_sgtf//sgtf_england.csv") 
sgtf_dk = read.csv(".//data//omicron_sgtf//sgtf_denmark.csv") 
sgtf_be = read.csv(".//data//omicron_sgtf//sgtf_belgium.csv") 

sgtf = rbind(sgtf_sa, sgtf_eng, sgtf_dk, sgtf_be)
sgtf$date = as.Date(sgtf$date)
sgtf$date_num = as.numeric(sgtf$date)
levels_country = c("South Africa", "England", "Denmark", "Belgium")
sgtf$country = factor(sgtf$country, levels=levels_country)
sgtf$prop = sgtf$sgtf/sgtf$pos_tests
sgtf$non_sgtf = sgtf$pos_tests-sgtf$sgtf
# sgtf$baseline = 0.02 # 2% baseline subtracted
sgtf$obs = as.factor(1:nrow(sgtf)) # for observation-level random effect to take into account overdispersion
sgtf$random = factor(1) # fake random effect to be able to run glmmPQL
names(sgtf)
head(sgtf)[,c(1:4)]


# ANALYSIS OF SGTF DATA USING LOGISTIC REGRESSION ####

# fit (quasi)binomial GLM to SGTF data

fit_sgtf0 = glm(cbind(sgtf, non_sgtf) ~ date_num + country, family=binomial(logit), data=sgtf) # taking into account overdispersion via quasibinomial family
fit_sgtf1 = glm(cbind(sgtf, non_sgtf) ~ date_num * country, family=binomial(logit), data=sgtf) # taking into account overdispersion via quasibinomial family
AIC(fit_sgtf0, fit_sgtf1)
#            df      AIC
# fit_sgtf0  5 728.1557
# fit_sgtf1  8 668.5033 # fits best

summary(fit_sgtf1)

# fit_sgtf0 = glmmPQL(cbind(sgtf, non_sgtf) ~ date_num * country, family=binomial(logit), correlation=corAR1(), random=~1|obs, data=sgtf) # taking into account lag-1 autocorrelation in residuals
# fit_sgtf1 = glmmPQL(cbind(sgtf, non_sgtf) ~ date_num * country, family=binomial(logit), correlation=corAR1(), random=~1|obs, data=sgtf) # taking into account lag-1 autocorrelation in residuals
# AIC(fit_sgtf0, fit_sgtf1)
# summary(fit_sgtf)

plot(allEffects(fit_sgtf1))

# growth rate advantage of Omicron over Delta (difference in growth rate per day)
deltar_sgtf_bycountry = as.data.frame(emtrends(fit_sgtf1, ~ date_num, by="country", var="date_num", at=list(date_num=max(dateseq))))[,c(2,3,6,7)] # growth rate advantage of Omicron over Delta
rownames(deltar_sgtf_bycountry) = deltar_sgtf_bycountry$country
deltar_sgtf_bycountry$country = NULL
deltar_sgtf_bycountry
#              date_num.trend asymp.LCL asymp.UCL
# South Africa      0.2884629 0.2185518 0.3583740
# England           0.3986284 0.3829438 0.4143129
# Denmark           0.3234410 0.2966248 0.3502572
# Belgium           0.3082760 0.2490408 0.3675112

# avg growth rate advantage of Omicron over Delta across countries (mean difference in growth rate per day)
deltar_sgtf = as.data.frame(emtrends(fit_sgtf1, ~ date_num, var="date_num", at=list(date_num=max(dateseq))))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
deltar_sgtf
# date_num.trend asymp.LCL asymp.UCL
# 1      0.3297021 0.3055134 0.3538908
deltar_sgtf_char = sapply(deltar_sgtf, function (x) sprintf(as.character(round(x,2)), 2) )
deltar_sgtf_char = paste0(deltar_sgtf_char[1], " [", deltar_sgtf_char[2], "-", deltar_sgtf_char[3],"] 95% CLs")

# transmission advantage of Omicron over Delta (using generation time of 4.7 days)
transmadv_sgtf_by_country = exp(deltar_sgtf_bycountry*4.7) 
transmadv_sgtf_by_country
#              date_num.trend asymp.LCL asymp.UCL
# South Africa       3.879770  2.793216  5.388990
# England            6.511393  6.048656  7.009530
# Denmark            4.573015  4.031492  5.187276
# Belgium            4.258416  3.223577  5.625460

# avg transmission advantage of Omicron over Delta across countries (using generation time of 4.7 days)
transmadv_sgtf = exp(deltar_sgtf*4.7)
transmadv_sgtf_char = sapply(transmadv_sgtf, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_char = paste0(transmadv_sgtf_char[1], " [", transmadv_sgtf_char[2], "-", transmadv_sgtf_char[3],"] 95% CLs")

# plot model predictions
extrapolate = 30
date.from = as.Date("2021-10-14")
date.to = today+extrapolate

dateseq = as.numeric(seq(date.from, date.to, by=1))
emmeans_sgtf = as.data.frame(emmeans(fit_sgtf1, ~ date_num+country, at=list(date_num=dateseq, country=levels_country)), type="response")
colnames(emmeans_sgtf) = c("date_num", "country", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
emmeans_sgtf$date = as.Date(emmeans_sgtf$date_num, origin="1970-01-01")

qplot(data=sgtf, x=date, y=prop, geom="point", colour=country, fill=country, size=I(2)) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), alpha=I(0.5)) +
  geom_line(data=emmeans_sgtf, aes(y=prob), alpha=I(0.5), size=1) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  xlim(date.from, date.to) +
  coord_cartesian(ylim=c(0.001,0.99)) +
  ggtitle("Spread of Omicron variant inferred from SGTF data") +
  ylab("Share of Omicron among confirmed cases (%)") +
  annotate("text", x = c(as.Date("2021-10-14")), y = c(0.94),
           label = paste0("Growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char,
                          " per day\n\nTransmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  theme(legend.position="right") +
  xlab("") # + xaxis

ggsave(file=paste0(".\\plots\\",plotdir,"\\spread omicron logistic fit sgtf data.png"), width=8, height=6)


