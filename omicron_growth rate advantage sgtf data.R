# ANALYSIS OF GROWTH RATE ADVANTAGE OF OMICRON (B.1.1.529) IN BELGIUM BASED ON SGTF (S gene target failure / S dropout) DATA
# original data: Federal Test Platform, provided by Emmanuel André

# T. Wenseleers
# last update 9 DECEMBER 2021

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
# today = as.Date("2021-12-09")
today_num = as.numeric(today)
plotdir = "omicron_sgtf"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

# import SGTF data (now proxy for B.1.1.529 / Omicron)
sgtf = read.csv(".//data//omicron_sgtf//sgtf_belgium.csv") 
sgtf$date = as.Date(sgtf$date)
sgtf$date_num = as.numeric(sgtf$date)
sgtf$prop = sgtf$sgtf/sgtf$pos_tests
sgtf$non_sgtf = sgtf$pos_tests-sgtf$sgtf
# sgtf$baseline = 0.02 # 2% baseline subtracted
sgtf$obs = as.factor(1:nrow(sgtf)) # for observation-level random effect to take into account overdispersion
sgtf$random = factor(1) # fake random effect to be able to run glmmPQL
names(sgtf)


# ANALYSIS OF SGTF DATA USING LOGISTIC REGRESSION ####

# fit (quasi)binomial GLM to SGTF data

# fit_sgtf = glm(cbind(sgtf, non_sgtf) ~ date_num, family=binomial(logit), data=sgtf)
# summary(fit_sgtf)

fit_sgtf = glm(cbind(sgtf, non_sgtf) ~ date_num, family=quasibinomial(logit), data=sgtf) # taking into account overdispersion via quasibinomial family
summary(fit_sgtf)

# fit_sgtf = glmmPQL(cbind(sgtf, non_sgtf) ~ date_num, family=binomial(logit), correlation=corAR1(), random=~1|obs, data=sgtf) # taking into account lag-1 autocorrelation in residuals
# summary(fit_sgtf)

plot(allEffects(fit_sgtf))

extrapolate = 30
date.from = as.Date("2021-10-14")
date.to = today+extrapolate

dateseq = as.numeric(seq(date.from, date.to, by=1))
emmeans_sgtf = as.data.frame(emmeans(fit_sgtf, ~ date_num, at=list(date_num=dateseq)), type="response")
colnames(emmeans_sgtf) = c("date_num","prob","SE","df","asymp.LCL","asymp.UCL")
emmeans_sgtf$date = as.Date(emmeans_sgtf$date_num, origin="1970-01-01")

deltar_sgtf = as.data.frame(emtrends(fit_sgtf, ~ date_num, var="date_num", at=list(date_num=max(dateseq))))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
deltar_sgtf_char = sapply(deltar_sgtf, function (x) sprintf(as.character(round(x,2)), 2) )
deltar_sgtf_char = paste0(deltar_sgtf_char[1], " [", deltar_sgtf_char[2], "-", deltar_sgtf_char[3],"] 95% CLs")
transmadv_sgtf = exp(deltar_sgtf*4.7) # transmission advantage with gen time of 4.7 days
transmadv_sgtf_char = sapply(transmadv_sgtf, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_char = paste0(transmadv_sgtf_char[1], " [", transmadv_sgtf_char[2], "-", transmadv_sgtf_char[3],"] 95% CLs")

qplot(data=sgtf, x=date, y=prop, geom="point", colour=I("steelblue"), size=I(2)) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), alpha=I(0.5), fill=I("steelblue")) +
  geom_line(data=emmeans_sgtf, aes(y=prob), colour=I("steelblue"), alpha=I(0.5), size=1) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  xlim(date.from, date.to) +
  coord_cartesian(ylim=c(0.001,0.99)) +
  ggtitle("Spread of Omicron in Belgium inferred from SGTF data", "Data Federal Test Platform") +
  ylab("Share of Omicron among confirmed cases (%)") +
  annotate("text", x = c(as.Date("2021-10-14")), y = c(0.94), 
           label = paste0("Growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char,
                          " per day\n\nTransmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
           color="black", hjust=0, size=3) +
  theme_hc() +
  xlab("") # + xaxis

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_spread omicron logistic fit sgtf data.png"), width=8, height=6)


