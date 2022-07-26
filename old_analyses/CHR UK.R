# CALCULATION OF CASE HOSPITALISATION RATE IN THE UK
# data gov.uk

library(tidyr)
library(Hmisc)
library(ggplot2)
library(scales)
library(ggthemes)
library(mgcv)
library(lubridate)
library(emmeans)

plotdir = "UK_SANGER"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

cases_hosps_uk = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newAdmissions&metric=newDeaths28DaysByDeathDate&metric=newTestsByPublishDate&metric=newCasesByPublishDate&metric=cumPeopleVaccinatedCompleteByPublishDate&format=csv")
cases_hosps_uk$date = as.Date(cases_hosps_uk$date)
cases_hosps_uk = cases_hosps_uk[order(cases_hosps_uk$date,decreasing=F),]
lag=7
cases_hosps_uk$newCasesByPublishDate = Lag(cases_hosps_uk$newCasesByPublishDate,lag)
cases_hosps_uk$newTestsByPublishDate = Lag(cases_hosps_uk$newTestsByPublishDate,lag)
cases_hosps_uk$CHR = 100*cases_hosps_uk$newAdmissions / cases_hosps_uk$newCasesByPublishDate
cases_hosps_uk$CHR[cases_hosps_uk$CHR>=100] = 100
cases_hosps_uk$CFR = 100*cases_hosps_uk$newDeaths28DaysByDeathDate / cases_hosps_uk$newCasesByPublishDate
cases_hosps_uk$CFR[cases_hosps_uk$CFR>=100] = 100
cases_hosps_uk$perc_fully_vaccinated = 100*cases_hosps_uk$cumPeopleVaccinatedCompleteByPublishDate/67100000
cases_hosps_uk$perc_fully_vaccinated[is.na(cases_hosps_uk$perc_fully_vaccinated)] = 0
cases_hosps_uk$posrate = 100*cases_hosps_uk$newCasesByPublishDate/cases_hosps_uk$newTestsByPublishDate
cases_hosps_uk$reciprposrate = 1/cases_hosps_uk$posrate # reciprocal of the positivity rate
cases_hosps_uk$DATE_NUM = as.numeric(cases_hosps_uk$date)
cases_hosps_uk$weekday = factor(weekdays(cases_hosps_uk$date))
cases_hosps_uk$lognewTestsByPublishDate = log(cases_hosps_uk$newTestsByPublishDate)
# cases_hosps_uk = cases_hosps_uk[complete.cases(cases_hosps_uk),]

fit_CHR = gam(cbind(newAdmissions, newCasesByPublishDate) ~ s(DATE_NUM, bs="cs", k=30) + weekday # offset(log(newTestsByPublishDate)),
                ,# lognewTestsByPublishDate, 
              family=binomial(cloglog), data=cases_hosps_uk)
BIC(fit_CHR)

CHR_emmeans = as.data.frame(emmeans(fit_CHR, ~ DATE_NUM, 
                                    at=list(DATE_NUM=seq(min(cases_hosps_uk$DATE_NUM), max(cases_hosps_uk$DATE_NUM)) #, 
                                                                     # newTestsByPublishDate=max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)
                                                                     # lognewTestsByPublishDate=log(10*max(cases_hosps_uk$newTestsByPublishDate, na.rm=T))
                                                                     # reciprposrate=0.001
), # offset = log(max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)/1.8),
   type="response"))
cases_hosps_uk$CHR_prob = CHR_emmeans$prob[match(cases_hosps_uk$DATE_NUM,CHR_emmeans$DATE_NUM)]
cases_hosps_uk$CHR_prob_LCL = CHR_emmeans$lower.CL[match(cases_hosps_uk$DATE_NUM,CHR_emmeans$DATE_NUM)]
cases_hosps_uk$CHR_prob_UCL = CHR_emmeans$upper.CL[match(cases_hosps_uk$DATE_NUM,CHR_emmeans$DATE_NUM)]


fit_CFR = gam(cbind(newDeaths28DaysByDeathDate, newCasesByPublishDate) ~ s(DATE_NUM, bs="cs", k=30) + weekday # offset(log(newTestsByPublishDate)),
              ,# lognewTestsByPublishDate, 
              family=binomial(cloglog), data=cases_hosps_uk)
BIC(fit_CFR)

CFR_emmeans = as.data.frame(emmeans(fit_CFR, ~ DATE_NUM, 
                                    at=list(DATE_NUM=seq(min(cases_hosps_uk$DATE_NUM), max(cases_hosps_uk$DATE_NUM)) #, 
                                            # newTestsByPublishDate=max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)
                                            # lognewTestsByPublishDate=log(10*max(cases_hosps_uk$newTestsByPublishDate, na.rm=T))
                                            # reciprposrate=0.001
                                    ), # offset = log(max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)/1.8),
                                    type="response"))
cases_hosps_uk$CFR_prob = CFR_emmeans$prob[match(cases_hosps_uk$DATE_NUM,CFR_emmeans$DATE_NUM)]
cases_hosps_uk$CFR_prob_LCL = CFR_emmeans$lower.CL[match(cases_hosps_uk$DATE_NUM,CFR_emmeans$DATE_NUM)]
cases_hosps_uk$CFR_prob_UCL = CFR_emmeans$upper.CL[match(cases_hosps_uk$DATE_NUM,CFR_emmeans$DATE_NUM)]



fit_cases = gam(newCasesByPublishDate ~ s(DATE_NUM, bs="cs", k=30) + weekday + # offset(log(newTestsByPublishDate)),
                  lognewTestsByPublishDate, #  s(reciprposrate, bs="cs", k=6),
                            family=poisson(log), data=cases_hosps_uk)
BIC(fit_cases)
# marginal mean new cases for average weekday if testing would have been done at intensity that would have been 10 x maximum testing intensity on any dat
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM, at=list(DATE_NUM=seq(min(cases_hosps_uk$DATE_NUM), max(cases_hosps_uk$DATE_NUM)), 
                                                                     # newTestsByPublishDate=max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)
                                                                     lognewTestsByPublishDate=log(10*max(cases_hosps_uk$newTestsByPublishDate, na.rm=T))
                                                                     # reciprposrate=0.001
                                                                     ), 
                                      # offset = log(max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)/1.8),
                                      type="response"))
cases_hosps_uk$newcases_smoothed = cases_emmeans$rate[match(cases_hosps_uk$DATE_NUM,cases_emmeans$DATE_NUM)]
cases_hosps_uk$newcases_smoothed_LCL = cases_emmeans$lower.CL[match(cases_hosps_uk$DATE_NUM,cases_emmeans$DATE_NUM)]
cases_hosps_uk$newcases_smoothed_UCL = cases_emmeans$upper.CL[match(cases_hosps_uk$DATE_NUM,cases_emmeans$DATE_NUM)]

qplot(data=cases_hosps_uk, x=date, y=newCasesByPublishDate, geom=c("line")) +
  geom_line(aes(y=newcases_smoothed), colour=I("red"))

qplot(data=cases_hosps_uk, x=date, y=100*newCasesByPublishDate/newcases_smoothed, geom=c("line")) + coord_cartesian(ylim=c(0,100))

cases_hosps_uk$IHR = 100*cases_hosps_uk$newAdmissions / cases_hosps_uk$newcases_smoothed

# INFECTION HOSPITALISATION RATE
qplot(data=cases_hosps_uk, x=date, y=IHR, geom=c("blank")) + 
  geom_point(pch=I(16), aes(colour=perc_fully_vaccinated)) + # I("steelblue")
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 20), fill=I("black"), colour=I("black"), alpha=I(0.5)) + 
  scale_y_log10() +
  ylab("IHR (%)") +
  ggtitle("INFECTION HOSPITALISATION RATE IN THE UK\nIN FUNCTION OF VACCINATION (data gov.uk)") +
  coord_cartesian(ylim=c(1,100)) +
  scale_colour_gradient("Fully vaccinated (%)", low="blue",high="red") + 
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) 
ggsave(file=paste0(".\\plots\\",plotdir,"\\infection hospitalisation rate vaccination UK.png"), width=8, height=6)

# CASE HOSPITALISATION RATE
# qplot(data=cases_hosps_uk, x=date, y=CHR, geom=c("blank")) + 
#   geom_point(pch=I(16), aes(colour=perc_fully_vaccinated)) + # I("steelblue")
#   geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 20), fill=I("black"), colour=I("black"), alpha=I(0.5)) + 
#   scale_y_log10() +
#   ylab("Case (-7 days) hospitalisation rate (%)") +
#   ggtitle("CASE (-7 DAYS) HOSPITALISATION RATE IN THE UK\nIN FUNCTION OF VACCINATION (data gov.uk)") +
#   coord_cartesian(ylim=c(1,20)) +
#   scale_colour_gradient("Fully vaccinated (%)", low="blue",high="red") +
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
#                      limits=as.Date(c("2020-08-01",NA)), expand=c(0,0)) 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\case hospitalisation rate vaccination UK.png"), width=8, height=6)
# 
qplot(data=cases_hosps_uk, x=date, y=100*CHR_prob, ymin=100*CHR_prob_LCL, ymax=100*CHR_prob_UCL, geom=c("blank")) + 
  geom_point(pch=I(16), aes(colour=perc_fully_vaccinated)) + # I("steelblue")
  geom_point(aes(y=CHR, colour=perc_fully_vaccinated), pch=I(16)) + # I("steelblue")
  geom_line() +
  geom_ribbon(alpha=I(0.5)) +
  # geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 20), fill=I("black"), colour=I("black"), alpha=I(0.5)) + 
  scale_y_log10() +
  ylab("Case (-7 days) hospitalisation rate (%)") +
  ggtitle("CASE (-7 DAYS) HOSPITALISATION RATE IN THE UK\nIN FUNCTION OF VACCINATION (data gov.uk)") +
  coord_cartesian(ylim=c(1,20)) +
  scale_colour_gradient("Fully vaccinated (%)", low="blue",high="red") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-08-01",NA)), expand=c(0,0)) 
ggsave(file=paste0(".\\plots\\",plotdir,"\\case hospitalisation rate vaccination UK_binomial fit.png"), width=8, height=6)


# CASE FATALITY RATE
# qplot(data=cases_hosps_uk, x=date, y=CFR, geom=c("blank")) + 
#   geom_point(pch=I(16), aes(colour=perc_fully_vaccinated)) + # I("steelblue")
#   geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 20), fill=I("black"), colour=I("black"), alpha=I(0.5)) + 
#   scale_y_log10() +
#   ylab("Case fatality rate (%)") +
#   ggtitle("CASE FATALITY RATE IN THE UK\nIN FUNCTION OF VACCINATION (data gov.uk)") +
#   coord_cartesian(ylim=c(0.1,10)) +
#   scale_colour_gradient("Fully vaccinated (%)", low="blue",high="red") +
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
#                      limits=as.Date(c("2020-08-01",NA)), expand=c(0,0)) 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\case fatality rate vaccination UK.png"), width=8, height=6)
# 
qplot(data=cases_hosps_uk, x=date, y=100*CFR_prob, ymin=100*CFR_prob_LCL, ymax=100*CFR_prob_UCL, geom=c("blank")) + 
  geom_point(pch=I(16), aes(colour=perc_fully_vaccinated)) + # I("steelblue")
  geom_point(aes(y=CFR, colour=perc_fully_vaccinated), pch=I(16)) + # I("steelblue")
  geom_line() +
  geom_ribbon(alpha=I(0.5)) +
  # geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 20), fill=I("black"), colour=I("black"), alpha=I(0.5)) + 
  scale_y_log10() +
  ylab("Case fatality rate (%)") +
  ggtitle("CASE FATALITY RATE IN THE UK\nIN FUNCTION OF VACCINATION (data gov.uk)") +
  coord_cartesian(ylim=c(0.1,10)) +
  scale_colour_gradient("Fully vaccinated (%)", low="blue",high="red") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-08-01",NA)), expand=c(0,0)) 
ggsave(file=paste0(".\\plots\\",plotdir,"\\case fatality rate vaccination UK_binomial fit.png"), width=8, height=6)


# CALCULATION OF AGE-SPECIFIC CASE FATALITY RATE ####

cases_byage = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDateAgeDemographics&format=csv")
admissions  = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newAdmissions&format=csv")
deaths_byage = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv")
cases_deaths_byage = cases_byage
cases_deaths_byage$rollingSum = NULL
cases_deaths_byage$rollingRate = NULL
cases_deaths_byage$deaths = deaths_byage$deaths[match(interaction(cases_byage$date,
                                                                  cases_byage$age),
                                                      interaction(deaths_byage$date,
                                                                  deaths_byage$age))]
cases_deaths_byage$CFR = 100*cases_deaths_byage$deaths/cases_deaths_byage$cases
cases_deaths_byage$date = as.Date(cases_deaths_byage$date)
cases_deaths_byage$DATE_NUM = as.numeric(cases_deaths_byage$date)
cases_deaths_byage$weekday = factor(weekdays(cases_deaths_byage$date))
cases_deaths_byage = cases_deaths_byage[cases_deaths_byage$age!="unassigned",]
cases_deaths_byage = cases_deaths_byage[cases_deaths_byage$age!="60+",]
cases_deaths_byage = cases_deaths_byage[cases_deaths_byage$age!="00_59",]
cases_deaths_byage$age = factor(cases_deaths_byage$age)
cases_deaths_byage$deaths = as.numeric(cases_deaths_byage$deaths)
cases_deaths_byage$cases = as.numeric(cases_deaths_byage$cases)

fit_CFR_byage = gam(cbind(deaths, cases) ~ s(DATE_NUM, bs="cs", k=30, by=age) + age + weekday # offset(log(newTestsByPublishDate)),
              ,# lognewTestsByPublishDate, 
              family=binomial(cloglog), data=cases_deaths_byage)
BIC(fit_CFR_byage)

CFR_emmeans_byage = as.data.frame(emmeans(fit_CFR_byage, ~ DATE_NUM, by="age",
                                    at=list(DATE_NUM=seq(min(cases_hosps_uk$DATE_NUM), max(cases_hosps_uk$DATE_NUM)) #, 
                                            # newTestsByPublishDate=max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)
                                            # lognewTestsByPublishDate=log(10*max(cases_hosps_uk$newTestsByPublishDate, na.rm=T))
                                            # reciprposrate=0.001
                                    ), # offset = log(max(cases_hosps_uk$newTestsByPublishDate, na.rm=T)/1.8),
                                    type="response"))
cases_deaths_byage$CFR_prob = CFR_emmeans_byage$prob[match(interaction(cases_deaths_byage$DATE_NUM,
                                                                       cases_deaths_byage$age),
                                                           interaction(CFR_emmeans_byage$DATE_NUM,
                                                                       CFR_emmeans_byage$age))]
cases_deaths_byage$CFR_prob_LCL = CFR_emmeans_byage$lower.CL[match(interaction(cases_deaths_byage$DATE_NUM,
                                                                               cases_deaths_byage$age),
                                                                   interaction(CFR_emmeans_byage$DATE_NUM,
                                                                               CFR_emmeans_byage$age))]
cases_deaths_byage$CFR_prob_UCL = CFR_emmeans_byage$upper.CL[match(interaction(cases_deaths_byage$DATE_NUM,
                                                                               cases_deaths_byage$age),
                                                                   interaction(CFR_emmeans_byage$DATE_NUM,
                                                                               CFR_emmeans_byage$age))]



# CASE FATALITY RATE BY AGE
qplot(data=cases_deaths_byage, x=date, y=CFR, geom=c("point"), colour=age) + 
  facet_wrap(~ age, scale="free") +
  # geom_point(pch=I(16), aes(colour=perc_fully_vaccinated)) + # I("steelblue")
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 20), fill=I("black"), colour=I("black"), alpha=I(0.5)) + 
  scale_y_log10() +
  ylab("CFR (%)") +
  ggtitle("CASE FATALITY RATE BY AGE IN ENGLAND (data gov.uk)") +
  # coord_cartesian(ylim=c(0.01,100)) +
  # scale_colour_gradient("Fully vaccinated (%)", low="blue",high="red") + 
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-08-01",NA)), expand=c(0,0)) + 
  theme(legend.position="none")
ggsave(file=paste0(".\\plots\\",plotdir,"\\case fatality rate by age england.png"), width=10, height=8)

qplot(data=cases_deaths_byage, x=date, y=100*CFR_prob, ymin=100*CFR_prob_LCL, ymax=100*CFR_prob_UCL, geom=c("blank"), colour=age) + 
  facet_wrap(~ age, scale="free") +
  geom_point(pch=I(16), aes(colour=age)) + # I("steelblue")
  geom_point(aes(y=CFR, colour=age), pch=I(16)) + # I("steelblue")
  geom_ribbon(alpha=I(0.5), colour=NA) +
  geom_line(lwd=I(1.2), colour=I("black")) +
  # geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 20), fill=I("black"), colour=I("black"), alpha=I(0.5)) + 
  scale_y_log10() +
  ylab("Case fatality rate (%)") +
  ggtitle("CASE FATALITY RATE BY AGE IN ENGLAND (data gov.uk)") +
  coord_cartesian(ylim=c(0.001,50)) +
  # scale_colour_gradient("Fully vaccinated (%)", low="blue",high="red") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-08-01",NA)), expand=c(0,0)) 
ggsave(file=paste0(".\\plots\\",plotdir,"\\case fatality rate by age england_binomial fit.png"), width=8, height=6)





vaccinations_byage = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=vaccinationsAgeDemographics&format=csv")
tests = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newTestsByPublishDate&format=csv")

