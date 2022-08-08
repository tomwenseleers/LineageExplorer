# EXCESS MORTALITY IN BELGIUM, STABEL, 
# https://statbel.fgov.be/en/about-statbel/what-we-do/visualisations/mortality

library(lubridate)
library(ggplot2)
library(ggthemes)
library(quantreg)
library(splines)
library(tidyr)
library(emmeans)
library(effects)
library(splines2)
library(Hmisc)
# install.packages("covidregionaldata",
#                  repos = "https://epiforecasts.r-universe.dev"
# )
library(covidregionaldata)

today = as.Date(Sys.time())
today_num = as.numeric(today)
today # "2021-08-05"
plotdir = "excess_mortality_belgium"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))
tag = paste("@TWenseleers\n",today)

mort = read.csv("https://www.dropbox.com/s/z4t12xdd2w6gj01/belstat_excess_mortality_belgium.csv?dl=1") # read.csv(".//data//excess_mortality//belstat_excess_mortality_belgium.csv")
mort_long = gather(mort, YEAR, nr_deaths, X2017:X2022, factor_key=TRUE)

cases_deaths_BE = get_national_data(countries="Belgium", source="who", level=1) # new cases & deaths

mort_long$DATE = as.Date(paste0(gsub("X","",mort_long$YEAR), "-", mort_long$DAY))
mort_long$DATE_NUM = as.numeric(mort_long$DATE)
mort_long$DAY_SINCE_JAN_2020 = mort_long$DATE_NUM-mort_long[mort_long$DATE==as.Date("2020-01-01"),"DATE_NUM"]+1
mort_long$DAY_OF_YEAR = yday(mort_long$DATE)
mort_long = mort_long[!is.na(mort_long$nr_deaths),]
mort_long$YEAR_FACTOR = factor(gsub("X","",mort_long$YEAR))  
mort_long$YEAR = as.numeric(as.character(mort_long$YEAR_FACTOR))
mort_long = mort_long[order(mort_long$DATE),]
mort_long$cases_new = cases_deaths_BE$cases_new[match(mort_long$DATE, cases_deaths_BE$date)]
mort_long$deaths_new = cases_deaths_BE$deaths_new[match(mort_long$DATE, cases_deaths_BE$date)]
mort_long$days_since_jan2017 = 1:nrow(mort_long)

# cyclical mSpline weighted median regression with temporal trend to fit baseline
# note that I use 1/(nr_deaths+0.1) weights as approx. 1/variance weights (appropriate for Poisson counts)
fit_cycl = rq(nr_deaths ~ mSpline(DAY_OF_YEAR, df=3, Boundary.knots=c(1,365), periodic=TRUE)+YEAR,
                     weights = 1/(nr_deaths+0.1), # approx 1/variance weights for Poisson counts
                     tau=0.5, data=mort_long[mort_long$YEAR<=2019,])
AIC(fit_cycl)

mort_long$fitted = predict(fit_cycl, newdata=data.frame(DAY_OF_YEAR=mort_long$DAY_OF_YEAR,
                                                        YEAR=mort_long$YEAR))
mort_long$excess_deaths = mort_long$nr_deaths - mort_long$fitted
mort_long$excess_deaths_since_start_pandemic = mort_long$excess_deaths
mort_long$excess_deaths_since_start_pandemic[mort_long$DATE<as.Date("2020-03-01")] = 0
mort_long$cum_excess_deaths_since_start_pandemic = cumsum(mort_long$excess_deaths_since_start_pandemic)
tail(mort_long$cum_excess_deaths_since_start_pandemic,1)
# 25926.7
sum(mort_long$excess_deaths_since_start_pandemic)
# 25926.7
sum(mort_long$excess_deaths[mort_long$YEAR=="2017"]) # 2566.04
sum(mort_long$excess_deaths[mort_long$YEAR=="2018"]) # 3261.277
sum(mort_long$excess_deaths[mort_long$YEAR=="2019"]) # 1040.513
sum(c(sum(mort_long$excess_deaths[mort_long$YEAR=="2017"]), sum(mort_long$excess_deaths[mort_long$YEAR=="2018"], sum(mort_long$excess_deaths[mort_long$YEAR=="2019"]))))/3
# 2289.277
sum(mort_long$excess_deaths[mort_long$YEAR=="2020"]) # 18490.53
sum(mort_long$excess_deaths[mort_long$YEAR=="2021"]) # 3942.986
sum(mort_long$excess_deaths[mort_long$YEAR=="2022"]) # 3146.508
sum(c(sum(mort_long$excess_deaths[mort_long$YEAR=="2020"]), sum(mort_long$excess_deaths[mort_long$YEAR=="2021"]), sum(mort_long$excess_deaths[mort_long$YEAR=="2022"])))
# 25580.02
sum(mort_long[mort_long$DAY_OF_YEAR>215&mort_long$DAY_OF_YEAR<238&mort_long$YEAR==2020,"excess_deaths"]) # 1643 - heatwave related
sum(mort_long$excess_deaths_since_start_pandemic)
# 25926.7
sum(mort_long$excess_deaths_since_start_pandemic)-sum(mort_long[mort_long$DAY_OF_YEAR>215&mort_long$DAY_OF_YEAR<238&mort_long$YEAR==2020,"excess_deaths"])
# 24283.9

# daily number of deaths+fitted baseline
qplot(data=mort_long, x=DAY_OF_YEAR, y=nr_deaths, group=YEAR_FACTOR, colour=YEAR_FACTOR, group=YEAR_FACTOR, geom="line") +
  facet_wrap(~YEAR_FACTOR, ncol=1) +
  scale_colour_hue("YEAR") +
  xlab("day since start of year") +
  ylab("daily number of deaths (all causes)") +
  geom_line(aes(y=fitted, lwd=I(1.5), alpha=I(0.5))) +
  theme_hc() +
  theme(legend.position="none") +
  labs(title="All-cause mortality in Belgium (period 2017-2022)",
       subtitle="Data Statbel, statbel.fgov.be/en/about-statbel/what-we-do/visualisations/mortality
baseline mortality fitted using weighted cyclical mSpline median regression
using model rq(nr_deaths ~ mSpline(DAY_OF_YEAR, df=3, Boundary.knots=c(1,365), periodic=T)+YEAR, 
weights = 1/(nr_deaths+0.1), tau=0.5, data=data[data$YEAR<=2019,])") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".\\plots\\",plotdir,"\\mortality Belgium with baseline.png"), width=9, height=13)

# this was effect of a heatwave
sum(mort_long[mort_long$DAY_OF_YEAR>215&mort_long$DAY_OF_YEAR<238&mort_long$YEAR==2020,"excess_deaths"]) # 1643
qplot(data=mort_long[mort_long$DAY_OF_YEAR>215&mort_long$DAY_OF_YEAR<238&mort_long$YEAR==2020,], 
      x=DAY_OF_YEAR, y=nr_deaths, group=YEAR_FACTOR, colour=YEAR_FACTOR, group=YEAR_FACTOR, geom="line") +
  # facet_wrap(~YEAR_FACTOR, ncol=1) +
  scale_colour_hue("YEAR") +
  xlab("day since start of year") +
  ylab("daily number of deaths (all causes)") +
  geom_line(aes(y=fitted, lwd=I(1.5), alpha=I(0.5))) +
  theme_hc() +
  theme(legend.position="none")

# inferred excess mortality
qplot(data=mort_long[mort_long$DATE>=as.Date("2020-03-01"),], 
      x=DATE, y=cum_excess_deaths_since_start_pandemic, geom="area", fill=I("grey40")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  scale_x_date(date_breaks = "months" , date_labels = "%b-%y") +
  theme_hc() +
  xlab("") +
  ylab("cumulative excess mortality since start pandemic") +
  coord_cartesian(expand=0) +
  labs(title="Cumulative excess mortality in Belgium (period 2017-2022)",
       subtitle="Data Statbel, statbel.fgov.be/en/about-statbel/what-we-do/visualisations/mortality
baseline mortality fitted using weighted cyclical mSpline median regression
using model rq(nr_deaths ~ mSpline(days_since_jan2017, df=3, Boundary.knots=c(1,365), periodic=T)+YEAR, 
weights = 1/(nr_deaths+0.1), tau=0.5, data=data[data$YEAR<=2019,])") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cumulative excess mortality Belgium.png"), width=9, height=8)

qplot(data=mort_long[mort_long$DATE>=as.Date("2020-03-01"),],
      x=Lag(deaths_new,-3),
      y=excess_deaths, geom="line", colour=YEAR_FACTOR) 

mort_long = mort_long %>% 
  mutate(excess_deaths_7dmovingavg = rollmean(excess_deaths, 7, na.pad = T),
         deaths_new_7dmovingavg = rollmean(deaths_new, 7, na.pad = T))
qplot(data=mort_long[mort_long$DATE>=as.Date("2020-03-01"),],
      x=Lag(deaths_new_7dmovingavg,-3),
      y=excess_deaths_7dmovingavg, geom="line", colour=YEAR_FACTOR) 



