# CALCULATION OF GLOBAL DEATH TOLL CAUSED BY DIFFERENT VARIANTS 
# CALCULATED FROM GISAID LINEAGE FREQUENCIES & ANALYSIS OF TOTAL DEATH TOLL BY THE ECONOMIST

# PLUS CALCULATION OF CUMULATIVE TOTAL OF PREVIOUSLY INFECTED PEOPLE

# T. Wenseleers, 9 Nov. 2021

fit_economist = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/output-for-interactive/second_map.csv")
# check "export_regions.csv", and "export_country.csv"
# https://github.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/tree/main/output-data
# for extimated excess mortality through time by country & continent
economist_excessdeaths_bycountry = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/export_country.csv")
economist_excessdeaths_bycountry$estimated_daily_excess_deaths_per_million = economist_excessdeaths_bycountry$estimated_daily_excess_deaths*1E6 / economist_excessdeaths_bycountry$population
economist_cumexcessdeaths_bycountry = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/output-for-interactive/by_location_per_100k_cumulative.csv")
economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million = 10*economist_cumexcessdeaths_bycountry$estimate 
economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million_march = economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million[match(interaction(economist_cumexcessdeaths_bycountry$location,"2021-03-01"),
                                                                                                                                                                              interaction(economist_cumexcessdeaths_bycountry$location,economist_cumexcessdeaths_bycountry$date))]
economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million_sincemarch = economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million - economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million_march
library(countrycode)
economist_cumexcessdeaths_bycountry$iso3c = countrycode(economist_cumexcessdeaths_bycountry$location, "country.name", "iso3c")

OWID = read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv")
OWID$estimated_daily_excess_deaths_per_million = economist_excessdeaths_bycountry$estimated_daily_excess_deaths_per_million[match(interaction(OWID$date,
                                                                                                                                              OWID$iso_code),
                                                                                                                                  interaction(economist_excessdeaths_bycountry$date,
                                                                                                                                              economist_excessdeaths_bycountry$iso3c))]

# OWID$estimated_cumulative_excess_deaths_per_million = 10*fit_economist$cumulative_estimated_daily_excess_deaths_per_100k[match(OWID$iso_code,fit_economist$iso3c)]
OWID$estimated_cumulative_excess_deaths_per_million = economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million[
                                                                                     match(interaction(OWID$iso_code,OWID$date),
                                                                                          interaction(economist_cumexcessdeaths_bycountry$iso3c,
                                                                                                      economist_cumexcessdeaths_bycountry$date))]
OWID$estimated_cumulative_excess_deaths_per_million_sincemarch = economist_cumexcessdeaths_bycountry$estimated_cumulative_excess_deaths_per_million_sincemarch[
  match(interaction(OWID$iso_code,OWID$date),
        interaction(economist_cumexcessdeaths_bycountry$iso3c,
                    economist_cumexcessdeaths_bycountry$date))]

names(OWID)
unique(OWID$continent)
library(ggplot2)
library(ggthemes)
(as.Date(max(OWID$date))-4) # "2021-11-04"
(as.Date(max(fit_economist$date))) # 2021-11-08

ggplot(OWID[(OWID$continent=="Europe")&(OWID$date==(as.Date(max(OWID$date))-7)),], aes(x=people_vaccinated_per_hundred, y=new_cases_per_million))+geom_point()+geom_text(aes(label=location, vjust=-0.7))  
ggplot(OWID[(OWID$continent=="Europe")&(OWID$date==(as.Date(max(OWID$date))-7)),], aes(x=people_vaccinated_per_hundred, y=hosp_patients_per_million))+geom_point()+geom_text(aes(label=location, vjust=-0.7))
ggplot(OWID[(OWID$continent=="Europe")&(OWID$date==(as.Date(max(OWID$date))-7)),], aes(x=people_vaccinated_per_hundred, y=icu_patients_per_million))+geom_point()+geom_text(aes(label=location, vjust=-0.7))
ggplot(OWID[(OWID$continent=="Europe")&(OWID$date==(as.Date(max(OWID$date))-7)),], aes(x=people_vaccinated_per_hundred, y=new_deaths_smoothed_per_million))+geom_point()+geom_text(aes(label=location, vjust=-0.7))
ggplot(OWID[(OWID$continent=="Europe")&(OWID$date==(as.Date(max(OWID$date))-7)),], aes(x=people_vaccinated_per_hundred, y=positive_rate))+geom_point()+geom_text(aes(label=location, vjust=-0.7))
ggplot(OWID[(OWID$continent=="Europe")&(OWID$date==(as.Date(max(OWID$date))-7)),], aes(x=total_boosters_per_hundred, y=new_deaths_smoothed_per_million))+geom_point()+geom_text(aes(label=location, vjust=-0.7))


tag = paste("@TWenseleers\n2021-10-25\ndata Our World in Data")
qplot(data=OWID[OWID$continent=="Europe"&OWID$date==(as.Date(max(OWID$date))-4),], # "2021-10-14"
      x=people_fully_vaccinated_per_hundred, # 
      y=new_deaths_smoothed_per_million,
      size=population,
      colour=I("steelblue"),
      alpha=I(0.7),
      pch=I(16),
      geom="point") + 
  geom_text(aes(label=location, vjust=I(-0.8)), size=I(3)) + # or excess_mortality or new_deaths_smoothed_per_million
  scale_size_continuous(range=c(2,6), trans="sqrt") +
  theme_few() +
  theme(legend.position="none") +
  ylab("Confirmed Covid deaths per million") +
  xlab("People fully vaccinated (%)") +
  ggtitle("Daily excess deaths per million vs.\npeople fully vaccinated against Covid in Europe") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
library(export)
graph2png(file="confirmed covid deaths vs vaccination coverage Europe.png", width=7, height=5)

tag = paste("@TWenseleers\n2021-10-25\ndata Our World in Data & The Economist")
qplot(data=OWID[OWID$continent=="Europe"&OWID$date==(as.Date(max(OWID$date))-3),], # "2021-10-14"
      x=people_fully_vaccinated_per_hundred, # 
      y=estimated_daily_excess_deaths_per_million,
      size=population,
      colour=I("steelblue"),
      alpha=I(0.7),
      pch=I(16),
      geom="point") + 
  geom_text(aes(label=location, vjust=I(-0.8)), size=I(3)) + # or excess_mortality or new_deaths_smoothed_per_million
  scale_size_continuous(range=c(2,6), trans="sqrt") +
  theme_few() +
  theme(legend.position="none") +
  ylab("Estimated daily excess deaths per million") +
  xlab("People fully vaccinated (%)") +
  ggtitle("Daily excess deaths per million vs.\npeople fully vaccinated against Covid in Europe") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
graph2png(file="estimated excess deaths vs vaccination coverage Europe.png", width=7, height=5)

tag = paste("@TWenseleers\n2021-10-25\ndata Our World in Data & The Economist")
qplot(data=OWID[OWID$continent=="Europe"&OWID$date==(as.Date(max(OWID$date))-3),], # "2021-10-14"
      x=people_fully_vaccinated_per_hundred, # 
      y=estimated_cumulative_excess_deaths_per_million,
      size=population,
      colour=I("steelblue"),
      alpha=I(0.7),
      pch=I(16),
      geom="point") + 
  geom_text(aes(label=location, vjust=I(-0.8)), size=I(3)) + # or excess_mortality or new_deaths_smoothed_per_million
  scale_size_continuous(range=c(2,6), trans="sqrt") +
  theme_few() +
  theme(legend.position="none") +
  ylab("Cumulative excess deaths per million since 1/1/2020") +
  xlab("People fully vaccinated (%)") +
  ggtitle("Cumulative excess deaths per million vs.\npeople fully vaccinated against Covid in Europe") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
graph2png(file="estimated cumulative excess deaths vs vaccination coverage Europe.png", width=7, height=5)

tag = paste("@TWenseleers\n2021-10-25\ndata Our World in Data & The Economist")
qplot(data=OWID[OWID$continent=="Europe"&OWID$date==(as.Date(max(OWID$date))-3),], # "2021-10-14"
      x=people_fully_vaccinated_per_hundred, # 
      y=estimated_cumulative_excess_deaths_per_million_sincemarch,
      size=population,
      colour=I("steelblue"),
      alpha=I(0.7),
      pch=I(16),
      geom="point") + 
  geom_text(aes(label=location, vjust=I(-0.8)), size=I(3)) + # or excess_mortality or new_deaths_smoothed_per_million
  scale_size_continuous(range=c(2,6), trans="sqrt") +
  theme_few() +
  theme(legend.position="none") +
  ylab("Cumulative excess deaths per million") +
  xlab("People fully vaccinated (%)") +
  ggtitle("Cumulative excess deaths per million since 1/3/2021 vs.\npeople fully vaccinated against Covid in Europe") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
graph2png(file="estimated cumulative excess deaths since March vs vaccination coverage Europe.png", width=7, height=5)




dat = OWID[OWID$continent=="Europe"&OWID$date==(as.Date(max(OWID$date))-3),c("new_deaths_smoothed_per_million","estimated_daily_excess_deaths_per_million","estimated_cumulative_excess_deaths_per_million","estimated_cumulative_excess_deaths_per_million_sincemarch",
                                                                             "people_fully_vaccinated_per_hundred","median_age","gdp_per_capita","location")]
dat = dat[complete.cases(dat),]
dat$estimated_cumulative_excess_deaths_per_million_beforemarch = dat$estimated_cumulative_excess_deaths_per_million-dat$estimated_cumulative_excess_deaths_per_million_sincemarch

fit_newdeath = lm(new_deaths_smoothed_per_million ~ people_fully_vaccinated_per_hundred+
                       median_age+gdp_per_capita, data=dat)
summary(fit_newdeath)
AIC(fit_newdeath)

library(effects)
plot(effect(mod=fit_newdeath, term="people_fully_vaccinated_per_hundred", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Daily confirmed Covid deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Daily confirmed Covid deaths per million vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale
graph2png(file="confirmed covid deaths vs vaccination coverage Europe adjusted for age GDP per capita.png", width=7, height=5)

fit_excessdeath = lm(estimated_daily_excess_deaths_per_million ~ people_fully_vaccinated_per_hundred+
                       median_age+gdp_per_capita, data=dat)
summary(fit_excessdeath)
AIC(fit_excessdeath)

library(effects)
plot(effect(mod=fit_excessdeath, term="people_fully_vaccinated_per_hundred", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Estimated daily excess deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Daily excess deaths per million vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale
graph2png(file="estimated excess deaths vs vaccination coverage Europe adjusted for age GDP per capita.png", width=7, height=5)

fit_cumexcessdeath = lm(estimated_cumulative_excess_deaths_per_million ~ people_fully_vaccinated_per_hundred+
                       median_age+gdp_per_capita, data=dat)
summary(fit_cumexcessdeath)
AIC(fit_cumexcessdeath)

library(effects)
plot(effect(mod=fit_cumexcessdeath, term="people_fully_vaccinated_per_hundred", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Cumulative excess deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Cumulative excess deaths per million vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale
graph2png(file="estimated cumulative excess deaths vs vaccination coverage Europe adjusted for age GDP per capita.png", width=7, height=5)


fit_cumexcessdeathsincemarch = lm(estimated_cumulative_excess_deaths_per_million_sincemarch ~ people_fully_vaccinated_per_hundred+
                          median_age+gdp_per_capita+estimated_cumulative_excess_deaths_per_million_beforemarch, data=dat)
summary(fit_cumexcessdeathsincemarch)
AIC(fit_cumexcessdeathsincemarch)

library(effects)
plot(effect(mod=fit_cumexcessdeathsincemarch, term="people_fully_vaccinated_per_hundred", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Cumulative excess deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Cumulative excess deaths per million since 1/3/2021 vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale
graph2png(file="estimated cumulative excess deaths since march vs vaccination coverage Europe adjusted for age GDP per capita.png", width=7, height=5)








dat = OWID[OWID$date==(as.Date(max(OWID$date))-3),c("new_deaths_smoothed_per_million","estimated_daily_excess_deaths_per_million","estimated_cumulative_excess_deaths_per_million",
                                                                             "people_fully_vaccinated_per_hundred","median_age","gdp_per_capita","continent",
                                                    "location")]
dat = dat[complete.cases(dat),]

fit_newdeath = lm(new_deaths_smoothed_per_million ~ people_fully_vaccinated_per_hundred+
                    median_age+gdp_per_capita, data=dat)
summary(fit_newdeath)
AIC(fit_newdeath)

library(effects)
plot(effect(mod=fit_newdeath, term="people_fully_vaccinated_per_hundred", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Daily confirmed Covid deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Daily confirmed Covid deaths per million vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale
graph2png(file="confirmed covid deaths vs vaccination coverage world adjusted for age GDP per capita.png", width=7, height=5)

fit_excessdeath = lm(estimated_daily_excess_deaths_per_million ~ people_fully_vaccinated_per_hundred+
                       median_age+gdp_per_capita, data=dat)
summary(fit_excessdeath)
AIC(fit_excessdeath)

library(effects)
plot(effect(mod=fit_excessdeath, term="people_fully_vaccinated_per_hundred", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Estimated daily excess deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Daily excess deaths per million vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale
graph2png(file="estimated excess deaths vs vaccination coverage world adjusted for age GDP per capita.png", width=7, height=5)


fit_cumexcessdeath = lm(estimated_cumulative_excess_deaths_per_million ~ people_fully_vaccinated_per_hundred*continent+
                       log(median_age)+gdp_per_capita, data=dat)
summary(fit_cumexcessdeath)
AIC(fit_cumexcessdeath)

library(effects)
plot(effect(mod=fit_cumexcessdeath, term="people_fully_vaccinated_per_hundred", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Estimated cumulative excess deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Cumulative excess deaths per million vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale
graph2png(file="estimated cumulative excess deaths vs vaccination coverage world adjusted for age GDP per capita.png", width=7, height=5)
plot(Effect(c("people_fully_vaccinated_per_hundred","continent"), fit_cumexcessdeath, residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="People fully vaccinated (%)", ylab="Estimated cumulative excess deaths per million",
     # axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), 
     type="rescale",
     main="Cumulative excess deaths per million vs.\npeople fully vaccinated against Covid in Europe\n(adjusted for differences in median age & GDP per capita)",
     id=list(n=nrow(dat), 
             labels=dat$location)) # on semilog scale








qplot(data=OWID[OWID$date==(as.Date(max(OWID$date))-7),], # "2021-10-14"
      x=people_fully_vaccinated_per_hundred, # 
      y=estimated_daily_excess_deaths_per_million,
      size=population,
      colour=I("steelblue"),
      alpha=I(0.7),
      pch=I(16),
      geom="point") + 
  geom_text(aes(label=location, vjust=I(-0.8))) + # or excess_mortality or new_deaths_smoothed_per_million
  scale_size_continuous(range=c(2,6), trans="sqrt")



  
# CUMULATIVE % INFECTED

fit_economist$population = OWID$population[match(fit_economist$iso3c, OWID$iso_code)]
fit_economist$population[fit_economist$iso3c=="PRK"] = 25670000 # North Korea, from World Bank, 2019
# fit_economist[is.na(fit_economist$population),]
fit_economist$continent = OWID$continent[match(fit_economist$iso3c, OWID$iso_code)]
fit_economist$continent[fit_economist$iso3c=="PRK"] = "Asia" # North Korea

range(fit_economist$implied_infections_per_100_persons, na.rm=T)
range(fit_economist$implied_infections_per_100_persons_top_95, na.rm=T)
range(fit_economist$implied_infections_per_100_persons_bot_95, na.rm=T)

fit_economist$implied_infections_per_100_persons_clipped = fit_economist$implied_infections_per_100_persons 
fit_economist$implied_infections_per_100_persons_clipped[fit_economist$implied_infections_per_100_persons_clipped<0] = 0
fit_economist$implied_infections_per_100_persons_clipped[fit_economist$implied_infections_per_100_persons_clipped>100] = 100

fit_economist$implied_infections_per_100_persons_top_95_clipped = fit_economist$implied_infections_per_100_persons_top_95
fit_economist$implied_infections_per_100_persons_top_95_clipped[fit_economist$implied_infections_per_100_persons_top_95_clipped<0] = 0
fit_economist$implied_infections_per_100_persons_top_95_clipped[fit_economist$implied_infections_per_100_persons_top_95_clipped>100] = 100

fit_economist$implied_infections_per_100_persons_bot_95_clipped = fit_economist$implied_infections_per_100_persons_bot_95
fit_economist$implied_infections_per_100_persons_bot_95_clipped[fit_economist$implied_infections_per_100_persons_bot_95_clipped<0] = 0
fit_economist$implied_infections_per_100_persons_bot_95_clipped[fit_economist$implied_infections_per_100_persons_bot_95_clipped>100] = 100

# estimated global % that was already infected by the virus on 8 Nov 2021
max(fit_economist$date) # "2021-11-08"
# 61% [28-70%]
weighted.mean(fit_economist$implied_infections_per_100_persons_clipped,
              w=fit_economist$population, na.rm=T)
weighted.mean(fit_economist$implied_infections_per_100_persons_bot_95_clipped,
              w=fit_economist$population, na.rm=T)
weighted.mean(fit_economist$implied_infections_per_100_persons_top_95_clipped,
              w=fit_economist$population, na.rm=T)

# estimated % already infected by the virus on 20 Sept 2021
# 
do.call(rbind, lapply(unique(fit_economist$continent), function (cont) { 
    list(continent = cont,
         implied_infections_per_100_persons = weighted.mean(fit_economist$implied_infections_per_100_persons_clipped[fit_economist$continent==cont],
                       w=fit_economist$population[fit_economist$continent==cont], na.rm=T),
         implied_infections_per_100_persons_bot_95 = weighted.mean(fit_economist$implied_infections_per_100_persons_bot_95_clipped[fit_economist$continent==cont],
                       w=fit_economist$population[fit_economist$continent==cont], na.rm=T),
         implied_infections_per_100_persons_top_95 = weighted.mean(fit_economist$implied_infections_per_100_persons_top_95_clipped[fit_economist$continent==cont],
                       w=fit_economist$population[fit_economist$continent==cont], na.rm=T)
    )
    } ))
# continent            implied_infections_per_100_persons implied_infections_per_100_persons_bot_95 implied_infections_per_100_persons_top_95
# [1,] "North America" 51.04899                           45.60529                                  54.17133                                 
# [2,] "Asia"          58.13057                           20.29919                                  65.45027                                 
# [3,] "Africa"        71.19838                           17.76086                                  93.77909                                 
# [4,] "Europe"        32.08875                           30.51862                                  33.78389                                 
# [5,] "South America" 70.04475                           64.50574                                  78.76023                                 
# [6,] "Oceania"       2.370656                           0.1146887                                 17.89174      



# 

