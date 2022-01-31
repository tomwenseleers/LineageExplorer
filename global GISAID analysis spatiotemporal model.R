# GLOBAL ANALYSIS OF GISAID DATA USING SPATIOTEMPORAL FIT TO GISAID LINEAGE FREQUENCIES ####
# T. Wenseleers
# 24 JANUARY 2022

library(mgcv)
library(nnet)
library(splines)
library(splines2)
library(ggplot2)
library(ggthemes)
library(scales)
library(tidyr)
library(dplyr)
library(ggpubr)
library(countrycode)

# 1. PARSE GISAID DATA & SETUP SOME USER PARAMETERS, COLOURS ETC ####

# TO PARSE LATEST GISAID METADATA DOWNLOAD LATEST METADATA FROM GISAID AND EXECUTE
# source("parse_GISAID.R")
# (it would also be possible to use the JSON stream, using source("download_GISAID_JSON.R"), 
# but I haven't finished the code for the parsing, see parse_GISAID_JSON.R - this would
# allow some more filtering, e.g. removing imported travel-related cases etc)

today = as.Date(Sys.time()) 
today_num = as.numeric(today)
plotdir = "GISAID_spatiotemporal_model"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

tag = paste("@TWenseleers\n",today)

VARIANT_levels = c("Other", "Alpha", "Beta", "Gamma", "Delta", "Omicron")
VARIANT_colours = c("grey70","#0085FF","#9A9D00","cyan3","magenta","blue")

GISAID_ftable = read.csv(".//data//GISAID//GISAID_freq_VARIANT_byweek_country.csv")
GISAID_ftable$VARIANT = factor(GISAID_ftable$VARIANT, levels=VARIANT_levels)
countries_levels = sort(unique(GISAID_ftable$country))
continents_levels = sort(unique(GISAID_ftable$continent))
GISAID_ftable$country = factor(GISAID_ftable$country, levels=countries_levels)
GISAID_ftable$continent = factor(GISAID_ftable$continent, levels=continents_levels)
GISAID_ftable$week_midpoint = as.Date(GISAID_ftable$week_midpoint)
GISAID_ftable$date = GISAID_ftable$week_midpoint
countries_long_lat = read.csv(".//data//GISAID//countries_long_lat_centroids_curated.csv")
GISAID_ftable$iso3c = countries_long_lat$iso3c[match(GISAID_ftable$country, countries_long_lat$country)]
unique(GISAID_ftable$country[is.na(countries_long_lat$iso3c[match(GISAID_ftable$country, countries_long_lat$country)])])
# Crimea         Kosovo         Sint Eustatius   no iso3c code
# TO DO: fix this in parse_GISAID and recode these as "Ukraine", "Serbia" and "Caribbean Netherland"


# 2. MULTINOMIAL SPATIOTEMPORAL TENSOR SPLINE FITS ####
set.seed(1)
# note: read in the saved fits below to skip this
# fitA = nnet::multinom(VARIANT ~ country + continent*ns(date_num, df=2), 
#                      weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fitB = nnet::multinom(VARIANT ~ continent*ns(date_num, df=2), 
#                       weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit0 = nnet::multinom(VARIANT ~ ns(date_num, df=2) + 
#                              mSpline(latitude, df=5, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                              mSpline(longitude, df=5, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                            weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit1 = nnet::multinom(VARIANT ~ continent + ns(date_num, df=2) + 
#                              mSpline(latitude, df=5, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                              mSpline(longitude, df=5, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                            weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit2 = nnet::multinom(VARIANT ~ continent * ns(date_num, df=2) + 
#                              mSpline(latitude, df=5, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                              mSpline(longitude, df=5, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                            weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit2B = nnet::multinom(VARIANT ~ continent * ns(date_num, df=2) + 
#                         mSpline(latitude, df=10, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                         mSpline(longitude, df=10, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                       weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
fit2C = nnet::multinom(VARIANT ~ continent * ns(date_num, df=2) + # BEST PER-CONTINENT SPATIAL FIT
                         mSpline(latitude, df=12, Boundary.knots=c(-90, 90), periodic=FALSE) :
                         mSpline(longitude, df=12, Boundary.knots=c(-180, 180), periodic=TRUE),
                       weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit3 = nnet::multinom(VARIANT ~ country + continent * ns(date_num, df=2) + 
#                         mSpline(latitude, df=5, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                         mSpline(longitude, df=5, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                       weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
fit3B = nnet::multinom(VARIANT ~ country + continent * ns(date_num, df=2) + # BEST PER-COUNTRY SPATIAL FIT
                        mSpline(latitude, df=10, Boundary.knots=c(-90, 90), periodic=FALSE) : 
                        mSpline(longitude, df=10, Boundary.knots=c(-180, 180), periodic=TRUE), 
                      weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit3C = nnet::multinom(VARIANT ~ country + continent * ns(date_num, df=2) + 
#                          mSpline(latitude, df=12, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                          mSpline(longitude, df=12, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                        weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit4 = nnet::multinom(VARIANT ~ country + continent + ns(date_num, df=2) + 
#                              mSpline(latitude, df=5, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                              mSpline(longitude, df=5, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                            weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit5 = nnet::multinom(VARIANT ~ country + ns(date_num, df=2) + 
#                         mSpline(latitude, df=5, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                         mSpline(longitude, df=5, Boundary.knots=c(-180, 180), periodic=TRUE), 
#                       weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)
# fit6 = nnet::multinom(VARIANT ~ country + continent + # with spline ifo lat x long and time
#                         (mSpline(latitude, df=5, Boundary.knots=c(-90, 90), periodic=FALSE) : 
#                         mSpline(longitude, df=5, Boundary.knots=c(-180, 180), periodic=TRUE)) : 
#                         ns(date_num, df=2), 
#                       weights=Freq, data=GISAID_ftable, maxit=10000, MaxNWts=10000)

BIC(fitA, fitB, fit0, fit1, fit2, fit2B, fit2C, fit3, fit3B, fit3C, fit4, fit5, fit6) 
# df     BIC
# fitA  1085 3675885
# fitB    90 4406583
# fit0   140 3990355
# fit1   165 3962192
# fit2   215 3902727
# fit2B  575 3735257
# fit2C  775 3694524
# fit3  1085 3675889
# fit3B 1085 3675876
# fit3C 1085 3675959
# fit4  1035 3715640
# fit5  1035 3705240
# fit6  1275 3877216
# best per-country fit: fit3B (best especially for countries with a decent amount of data)
# best per-continent fit: fit2C (best for countries with limited data or countries not in GISAID)

# PS1 an mclogit::mblogit fit would also be possible, and would also allow overdispersion to be taken into account
# PS2 a glmnet fit with some LASSO or elastic net regularisation might perform better for fit3C & would be other way to prevent overfitting
# PS3 I also tried an mgcv gam fit, but gam doesn't support multinomial models with frequency weights unfortunately, and counts are
# too large & dataset too large to make it feasible to enter the data in long format; gam supports "splines on a sphere" though, which could be better
# fit3 = gam(formula = list(VARIANT ~ country + continent + te(longitude, latitude, k=c(5, 5), bs=c("cc", "cs")) + s(date_num, k=2, bs="cs", by=continent),
#                                   ~ country + continent + te(longitude, latitude, k=c(5, 5), bs=c("cc", "cs")) + s(date_num, k=2, bs="cs", by=continent),
#                                   ~ country + continent + te(longitude, latitude, k=c(5, 5), bs=c("cc", "cs")) + s(date_num, k=2, bs="cs", by=continent),
#                                   ~ country + continent + te(longitude, latitude, k=c(5, 5), bs=c("cc", "cs")) + s(date_num, k=2, bs="cs", by=continent),
#                                   ~ country + continent + te(longitude, latitude, k=c(5, 5), bs=c("cc", "cs")) + s(date_num, k=2, bs="cs", by=continent),
#                                   ~ country + continent + te(longitude, latitude, k=c(5, 5), bs=c("cc", "cs")) + s(date_num, k=2, bs="cs", by=continent),
#            data=GISAID_ftable, family=multinom(K=6), weight=Freq, method="REML") # PS multinom doesn't work with Frequency weights, so no good...


# saveRDS(fitA, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fitA.rds")
# saveRDS(fitB, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fitB.rds")
# saveRDS(fit0, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit0.rds")
# saveRDS(fit1, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit1.rds")
# saveRDS(fit2, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit2.rds")
# saveRDS(fit2B, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit2B.rds")
# saveRDS(fit2C, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit2C.rds")
# saveRDS(fit3, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit3.rds")
# saveRDS(fit3B, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit3B.rds")
# saveRDS(fit3C, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit3C.rds")
# saveRDS(fit4, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit4.rds")
# saveRDS(fit5, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit5.rds")
# saveRDS(fit6, file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit6.rds")
# PS using higher df for the splines resulted in unstable fits; I now used more or less the highest df that would still return a stable fit

fit2C = readRDS(file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit2C.rds")
fit3B = readRDS(file=".//fits//global_GISAID_analysis_spatiotemporal_model/fit3B.rds")



# 3. PLOT MODEL PREDICTIONS ####

# countries with sufficient data for which per-country spatial fit fit3B looks best; 
# for the remaining countries & any other countries not in GISAID we will use the per-continent spatial fit fit2B:
# TO DO: define good automatic threshold to identify these countries?
countries_sufficient_data = c("Botswana", "Gambia", "Ghana", "Kenya", "Malawi",
                              "Mauritius", "Mayotte", "Morocco", "Mozambique", "Republic of the Congo",
                              "Reunion", "Senegal", "South Africa", "Uganda",
                              "Bangladesh", "Cambodia", "China", "Georgia", "Hong Kong",
                              "India", "Indonesia", "Iran", "Iraq", "Israel",              
                              "Japan", "Jordan", "Kazakhstan", "Malaysia",            
                              "Nepal", "Oman", "Pakistan", "Singapore", "South Korea",
                              "Sri Lanka", "Thailand", "Austria", "Belgium", "Bosnia and Herzegovina",
                              "Bulgaria", "Croatia", "Czech Republic", "Denmark",
                              "Estonia", "Finland","France", "Germany", "Gibraltar", "Greece",
                              "Iceland", "Ireland", "Italy", "Lithuania", "Luxembourg",
                              "Malta", "Netherlands", "North Macedonia", "Norway", "Poland",
                              "Portugal", "Romania", "Russia", "Serbia", "Slovakia",
                              "Slovenia", "Spain", "Sweden", "Switzerland","Turkey", 
                              "United Kingdom", "Canada", "Costa Rica", "Guadeloupe",                       
                              "Martinique", "Mexico", "Puerto Rico", "Sint Maarten", "USA", 
                              "Australia", "New Zealand", "Argentina", "Aruba",
                              "Bonaire", "Brazil", "Chile", "Colombia", "Curacao",           
                              "Ecuador", "French Guiana", "Paraguay", "Peru",
                              "Suriname", "Trinidad and Tobago")
iso3c_sufficient_data = countries_long_lat$iso3c[match(countries_sufficient_data, countries_long_lat$country)] # countrycode(countries_sufficient_data, "country.name", "iso3c")

# check which countries were in the analysis of the global death toll of Covid of The Economist
fit_economist = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/output-for-interactive/second_map.csv")

# remaining countries to get predicted lineage frequencies for based on per-continent spatial fit fit2C
iso3c_remaining = unique(fit_economist$iso3c)[!(unique(fit_economist$iso3c) %in% iso3c_sufficient_data)]
countries_remaining = countrycode(iso3c_remaining, "iso3c", "country.name")

# # check "export_regions.csv", and "export_country.csv"
# # https://github.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/tree/main/output-data
# # for extimated excess mortality through time by country & continent
# economist_excessdeaths_bycountry = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/export_country.csv")
# economist_excessdeaths_bycountry$estimated_daily_excess_deaths_per_million = economist_excessdeaths_bycountry$estimated_daily_excess_deaths*1E6 / economist_excessdeaths_bycountry$population
# OWID = read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv")
# OWID$estimated_daily_excess_deaths_per_million = economist_excessdeaths_bycountry$estimated_daily_excess_deaths_per_million[match(interaction(OWID$date,
#                                                                                                                                               OWID$iso_code),
#                                                                                                                                   interaction(economist_excessdeaths_bycountry$date,
#                                                                                                                                               economist_excessdeaths_bycountry$iso3c))]

# mode predictions of spatial fits by country and by continent

fit_percountry = fit3B # per-country spatial fit to use for countries with sufficient data
fit_percontinent = fit2C # per-continent spatial fit to use for remaining countries (including for countries not in GISAID)

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

extrapolate = 30 # nr of days to extrapolate predictions for in plots of model predictions
date.from = min(GISAID_ftable$date)
date.to = max(GISAID_ftable$date)+extrapolate
date.from_num = as.numeric(date.from)
date.to_num = as.numeric(date.to)
dateseq = as.numeric(seq(date.from, date.to, by=1))

# predictions of per-country fit
predgrid1 = expand.grid(list(date_num=dateseq,
                            country=countries_levels))
predgrid1$longitude = countries_long_lat$longitude[match(predgrid1$country, countries_long_lat$country)]
predgrid1$latitude = countries_long_lat$latitude[match(predgrid1$country, countries_long_lat$country)]
predgrid1$continent = countries_long_lat$continent[match(interaction(predgrid1$longitude, predgrid1$latitude),
                                                         interaction(countries_long_lat$longitude, countries_long_lat$latitude))]

fit_preds1 = data.frame(predgrid1, as.data.frame(predict(fit_percountry, newdata=predgrid1, type="prob")),check.names=F)
fit_preds1 = gather(fit_preds1, VARIANT, prob, all_of(VARIANT_levels), factor_key=TRUE)
fit_preds1$date = as.Date(fit_preds1$date_num, origin="1970-01-01")
fit_preds1$VARIANT = factor(fit_preds1$VARIANT, levels=VARIANT_levels)
# fit_preds1$country = countries_long_lat$country[match(interaction(fit_preds1$latitude,fit_preds1$longitude),
#                                                      interaction(countries_long_lat$latitude,countries_long_lat$longitude))]
fit_preds1$country = factor(fit_preds1$country, levels=sort(unique(fit_preds1$country)))
# fit_preds1$continent = countries_long_lat$continent[match(fit_preds1$country,countries_long_lat$country)]
fit_preds1$continent = factor(fit_preds1$continent, levels=sort(unique(fit_preds1$continent)))

fit_preds1$iso3c = countries_long_lat$iso3c[match(fit_preds1$country,countries_long_lat$country)]


# prediction of per-continent fit
predgrid2 = expand.grid(list(date_num=dateseq,
                             country=countries_remaining))
predgrid2$iso3c = countrycode(predgrid2$country, "country.name", "iso3c")
predgrid2$longitude = countries_long_lat$longitude[match(predgrid2$iso3c, countries_long_lat$iso3c)]
predgrid2$latitude = countries_long_lat$latitude[match(predgrid2$iso3c, countries_long_lat$iso3c)]
predgrid2$continent = countries_long_lat$continent[match(predgrid2$iso3c,countries_long_lat$iso3c)] # countrycode(predgrid2$iso3c, "iso3c", "continent")

fit_preds2 = data.frame(predgrid2, as.data.frame(predict(fit_percontinent, newdata=predgrid2, type="prob")),check.names=F)
fit_preds2 = gather(fit_preds2, VARIANT, prob, all_of(VARIANT_levels), factor_key=TRUE)
fit_preds2$date = as.Date(fit_preds2$date_num, origin="1970-01-01")
fit_preds2$VARIANT = factor(fit_preds2$VARIANT, levels=VARIANT_levels)
# fit_preds2$country = countries_long_lat$country[match(interaction(fit_preds2$latitude,fit_preds2$longitude),
#                                                 interaction(countries_long_lat$latitude,countries_long_lat$longitude))]
fit_preds2$country = factor(fit_preds2$country, levels=sort(unique(fit_preds2$country)))
# fit_preds2$continent = countries_long_lat$continent[match(fit_preds$country,countries_long_lat$country)]
fit_preds2$continent = factor(fit_preds2$continent, levels=sort(unique(fit_preds2$continent)))

# combined predictions for all countries
fit_preds = rbind(fit_preds1[fit_preds1$iso3c %in% iso3c_sufficient_data,], fit_preds2)
continents_levels = sort(unique(fit_preds$continent))

# line plot of multinomial lineage fit by continent & country (on logit scale) ####
lineplots_logit = list(length(continents_levels))
countries_continent = list(length(continents_levels))
for (i in 1:length(continents_levels)) {  
  continent = continents_levels[[i]]
  countries_continent[[i]] = unique(fit_preds[fit_preds$continent==continent,"country"])
  lineplots_logit[[i]] = qplot(data=fit_preds[fit_preds$continent==continent,], x=date, y=prob, geom="blank") +
      facet_wrap(~ country, ncol=5) +
      # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
      #                 fill=VARIANT
      # ), alpha=I(0.3)) +
      geom_line(aes(y=prob,
                    colour=VARIANT
      ), alpha=I(1)) +
      ylab("Share (%)") +
      theme_hc() + 
      xlab("") +
      ggtitle(paste0("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\nIN ",toupper(continent)), 
              "(GISAID data, spatiotemporal multinomial fit)") +
      xaxis +
      scale_y_continuous( trans="logit", 
                          breaks=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5,0.9,0.99,0.999),
                          labels = c("0.001","0.01","0.1","1","10","50","90","99","99.9")
      ) +
      scale_fill_manual(values=VARIANT_colours) +
      scale_colour_manual(values=VARIANT_colours) +
      geom_point(data=GISAID_ftable[GISAID_ftable$continent==continent,],
                 aes(x=week_midpoint, y=Prop, size=Total,
                     colour=VARIANT
                 ),
                 alpha=I(1)) +
      scale_size_continuous("total number\nsequenced", trans="sqrt",
                            range=c(0.1, 2), limits=c(1,max(GISAID_ftable$Total)), breaks=c(10,100,1000,10000)) +
      # guides(fill=FALSE) +
      # guides(colour=FALSE) +
      theme(legend.position = "right") +
      xlab("Collection date") +
      coord_cartesian(ylim=c(0.001, 0.9901), expand=c(0,0))
  
  # ggsave(file=paste0(".\\plots\\",plotdir,"\\line plot_fit_",continent,"_logit scale.png"), width=10, height=6)
}

# save.image("~/GitHub/newcovid_belgium/environment_24_1_2022.RData")

lineplots_logit[[1]] # not OK: Cabo Verde= Cape Verde, Cote d'Ivoire, Democratic Republic of the Congo, Ethiopia, Niger, Nigeria, South Sudan, Chad
lineplots_logit[[2]] # not OK: Myanmar = Myanmar (Burma), Palestine = Palestinian Territories 
lineplots_logit[[3]] # not OK: Azerbaijan, Canary Islands, Crimea, Kosovo, Turkey, 
lineplots_logit[[4]] # not OK: Antigua and Barbuda = Antigua & Barbuda, Bahames = The Bahamas, Aruba, Barbados, Bermuda, Bonaire, Curacao, Domincican Republic, Haiti, Saint Barthelemy, Saint Kitts and Nevis = St. Kitts & Nevis, Saint Lucia = St. Lucia, Saint Martin, Saint Vincent and the Grenadines, Sint Eustatius, The Bahamas, Trinidad and Tobago, Turks and Caicos Islands, St. Vincent & Grenadines
lineplots_logit[[5]] # not OK: Fiji, Guam, Northern Mariana Islands, Vanuatu, Wallis and Futuna Islands = Wallis & Futuna, Marshall Islands, Nauru, Tuvalu, 
lineplots_logit[[6]] # not OK: Aruba, Bonaire, Curacao, Trinidad and Tobago, 






# AFRICA
countries_continent[[1]]
# Botswana                        Gambia                           Ghana      Kenya                            
# Malawi                           Mauritius                       
# [31] Mayotte                          Morocco                          Mozambique                       Republic of the Congo            Reunion  
# Senegal South Africa                     Uganda

# ASIA
countries_continent[[2]]
# Bangladesh           Cambodia             China               
# Georgia              Hong Kong            India                Indonesia            Iran                 Iraq                 Israel              
# Japan                Jordan               Kazakhstan           Malaysia            
# Nepal                Oman                 Pakistan                     
# Singapore            South Korea          Sri Lanka            Thailand

# EUROPE
countries_continent[[3]]
# Austria                Belgium                Bosnia and Herzegovina
# Bulgaria               Croatia                Czech Republic         Denmark               
# Estonia                Finland                France                 Germany                Gibraltar              Greece                
# Iceland                Ireland                Italy                      
# Lithuania              Luxembourg             Malta                  Netherlands           
# North Macedonia        Norway                 Poland                 Portugal               Romania                Russia                 Serbia                
# Slovakia               Slovenia               Spain                  Sweden                 Switzerland            Turkey                                
# United Kingdom

# NORTH AMERICA
countries_continent[[4]]
# Canada                           Costa Rica                       Guadeloupe                      
# Martinique                      
# Mexico                           Puerto Rico                      Sint Maarten                     USA 

# OCEANIA
countries_continent[[5]]
# Australia                 New Zealand               

# SOUTH AMERICA
countries_continent[[6]]
# Argentina           Aruba               Bonaire             Brazil              Chile               Colombia            Curacao            
# Ecuador             French Guiana       Paraguay            Peru                Suriname            Trinidad and Tobago 



# # Muller plot of multinomial lineage fit by continent & country
# for (continent in continents_levels) {  
#   print(
#     ggplot(data=fit_preds[fit_preds$continent==continent,], 
#                           aes(x=date, y=prob, group=VARIANT)) + 
#   facet_wrap(~ country, ncol=5) +
#   scale_y_continuous(expand=c(0,0)) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=VARIANT, group=VARIANT), position="stack") +
#   scale_fill_manual("", values=VARIANT_colours) +
#   # annotate("rect", xmin=max(GISAID_sel$date_num)+1, 
#   #          xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.2, fill="white") + # extrapolated part
#   xaxis +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("Share") +
#   ggtitle(paste0("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\nIN ",toupper(continent)), 
#           "(GISAID data, spatiotemporal multinomial fit)") 
#   )
#   
#   # ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_fit_",continent,".png"), width=10, height=6)
# }








# 4. multinomial fit fit3_nnet ifo latitude & longitude & date ####
date.from_num = min(GISAID_ftable$date_num)
date.to_num = max(GISAID_ftable$date_num)

predgrid = expand.grid(list(date_num=seq(date.from_num, date.to_num),
                            country=countries_levels))
predgrid$longitude = GISAID_ftable$longitude[match(predgrid$country, GISAID_ftable$country)]
predgrid$latitude = GISAID_ftable$latitude[match(predgrid$country, GISAID_ftable$country)]
predgrid$continent = GISAID_ftable$continent[match(predgrid$country, GISAID_ftable$country)]

fit3_nnet_bycountry = data.frame(predgrid, as.data.frame(predict(fit3_nnet, newdata=predgrid, type="prob")),check.names=F)
fit3_nnet_bycountry = gather(fit3_nnet_bycountry, VARIANT, prob, all_of(VARIANT_levels), factor_key=TRUE)
fit3_nnet_bycountry$date = as.Date(fit3_nnet_bycountry$date_num, origin="1970-01-01")
fit3_nnet_bycountry$VARIANT = factor(fit3_nnet_bycountry$VARIANT, levels=VARIANT_levels)
fit3_nnet_bycountry$country = GISAID_ftable$country[match(interaction(fit3_nnet_bycountry$latitude,fit3_nnet_bycountry$longitude),
                                                          interaction(GISAID_ftable$latitude,GISAID_ftable$longitude))]
fit3_nnet_bycountry$country = factor(fit3_nnet_bycountry$country, levels=countries_levels)
fit3_nnet_bycountry$continent = GISAID_ftable$continent[match(fit3_nnet_bycountry$country,GISAID_ftable$country)]
fit3_nnet_bycountry$continent = factor(fit3_nnet_bycountry$continent, levels=continents_levels)

# Muller plot of multinomial lineage fit by continent & country

for (continent in continents_levels) {  
  print(
    ggplot(data=fit1_nnet_bycountry[fit1_nnet_bycountry$continent==continent,], 
           aes(x=date, y=prob, group=VARIANT)) + 
      facet_wrap(~ country, ncol=5) +
      scale_y_continuous(expand=c(0,0)) +
      geom_area(aes(lwd=I(1.2), colour=NULL, fill=VARIANT, group=VARIANT), position="stack") +
      scale_fill_manual("", values=VARIANT_colours) +
      # annotate("rect", xmin=max(GISAID_sel$date_num)+1, 
      #          xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.2, fill="white") + # extrapolated part
      xaxis +
      theme_hc() + theme(legend.position="right", 
                         axis.title.x=element_blank()) + 
      ylab("Share") +
      ggtitle(paste0("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\nIN ",toupper(continent)), 
              "(GISAID data, spatiotemporal multinomial fit)") 
  )
  
  # ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_fit_",continent,".png"), width=10, height=6)
}

# line plot of multinomial lineage fit by continent & country (on logit scale)
for (continent in continents_levels) {  
  print(
    qplot(data=fit1_nnet_bycountry[fit1_nnet_bycountry$continent==continent,], x=date, y=prob, geom="blank") +
      facet_wrap(~ country, ncol=5) +
      # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
      #                 fill=VARIANT
      # ), alpha=I(0.3)) +
      geom_line(aes(y=prob,
                    colour=VARIANT
      ), alpha=I(1)) +
      ylab("Share (%)") +
      theme_hc() + xlab("") +
      ggtitle(paste0("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\nIN ",toupper(continent)), 
              "(GISAID data, spatiotemporal multinomial fit)") +
      xaxis +
      scale_y_continuous( trans="logit", 
                          breaks=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5,0.9,0.99,0.999),
                          labels = c("0.001","0.01","0.1","1","10","50","90","99","99.9")
      ) +
      scale_fill_manual(values=VARIANT_colours) +
      scale_colour_manual(values=VARIANT_colours) +
      geom_point(data=GISAID_ftable[GISAID_ftable$continent==continent,],
                 aes(x=week_midpoint, y=Prop, size=Total,
                     colour=VARIANT
                 ),
                 alpha=I(1)) +
      scale_size_continuous("total number\nsequenced", trans="sqrt",
                            range=c(0.1, 2), limits=c(1,max(GISAID_ftable$Total)), breaks=c(10,100,1000,10000)) +
      # guides(fill=FALSE) +
      # guides(colour=FALSE) +
      theme(legend.position = "right") +
      xlab("Collection date")+
      coord_cartesian(ylim=c(0.001, 0.9901), expand=c(0,0)) 
  )
  
  # ggsave(file=paste0(".\\plots\\",plotdir,"\\line plot_fit_",continent,"_logit scale.png"), width=10, height=6)
}





# growth rate advantage compared to Omicron (BA.1) (difference in growth rate per day) 
emtrbelgium = emtrends(fit3_belgium_multi, trt.vs.ctrl ~ LINEAGE,  
                       var="DATE_NUM",  mode="latent",
                       at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium = data.frame(confint(emtrbelgium, 
                                     adjust="none", df=NA)$contrasts, 
                             p.value=as.data.frame(emtrbelgium$contrasts)$p.value)
delta_r_belgium
# contrast    estimate          SE df  asymp.LCL   asymp.UCL      p.value
# 1          Alpha - Omicron (BA.1) -0.10413057 0.015987065 NA -0.1354646 -0.07279650 3.071108e-08
# 2           Beta - Omicron (BA.1) -0.04633638 0.031304273 NA -0.1076916  0.01501886 4.787073e-01
# 3          Gamma - Omicron (BA.1) -0.07100060 0.014913581 NA -0.1002307 -0.04177052 4.673001e-05
# 4          Delta - Omicron (BA.1) -0.14457236 0.007819939 NA -0.1598992 -0.12924556 1.318344e-10
# 5 Omicron (BA.2) - Omicron (BA.1)  0.20010030 0.049666145 NA  0.1027564  0.29744415 7.031223e-04
# 6          Other - Omicron (BA.1) -0.09242247 0.008910148 NA -0.1098860 -0.07495890 1.318609e-10

# pairwise growth rate difference (differences in growth rate per day) 
emtrbelgium_pairw = emtrends(fit3_belgium_multi, pairwise ~ LINEAGE,  
                             var="DATE_NUM",  mode="latent",
                             at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                           adjust="none", df=NA)$contrasts, 
                                   p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#                         contrast    estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1           Omicron (BA.1) - Alpha  0.10413057 0.015987065 NA  0.072796503  0.135464646 1.064215e-07
# 2            Omicron (BA.1) - Beta  0.04633638 0.031304273 NA -0.015018865  0.107691630 7.556511e-01
# 3           Omicron (BA.1) - Gamma  0.07100060 0.014913581 NA  0.041770516  0.100230679 1.574078e-04
# 4           Omicron (BA.1) - Delta  0.14457236 0.007819939 NA  0.129245558  0.159899155 1.582013e-10
# 5  Omicron (BA.1) - Omicron (BA.2) -0.20010030 0.049666145 NA -0.297444151 -0.102756441 2.273156e-03
# 6           Omicron (BA.1) - Other  0.09242247 0.008910148 NA  0.074958905  0.109886043 1.582575e-10
# 7                     Alpha - Beta -0.05779419 0.032781624 NA -0.122044993  0.006456610 5.763410e-01
# 8                    Alpha - Gamma -0.03312998 0.017894463 NA -0.068202480  0.001942527 5.178721e-01
# 9                    Alpha - Delta  0.04044178 0.013982038 NA  0.013037492  0.067846073 6.967950e-02
# 10          Alpha - Omicron (BA.2) -0.30423087 0.052111149 NA -0.406366846 -0.202094895 1.984969e-06
# 11                   Alpha - Other -0.01170810 0.013844020 NA -0.038841881  0.015425680 9.792232e-01
# 12                    Beta - Gamma  0.02466421 0.032651019 NA -0.039330607  0.088659036 9.884092e-01
# 13                    Beta - Delta  0.09823597 0.030359358 NA  0.038732725  0.157739223 2.773469e-02
# 14           Beta - Omicron (BA.2) -0.24643668 0.058646069 NA -0.361380862 -0.131492495 1.241968e-03
# 15                    Beta - Other  0.04608609 0.030459977 NA -0.013614368  0.105786550 7.363156e-01
# 16                   Gamma - Delta  0.07357176 0.012759046 NA  0.048564489  0.098579029 2.691296e-06
# 17          Gamma - Omicron (BA.2) -0.27110089 0.051791746 NA -0.372610850 -0.169590937 2.429803e-05
# 18                   Gamma - Other  0.02142188 0.013040634 NA -0.004137296  0.046981049 6.553355e-01
# 19          Delta - Omicron (BA.2) -0.34467265 0.050211796 NA -0.443085965 -0.246259340 2.247302e-08
# 20                   Delta - Other -0.05214988 0.004431786 NA -0.060836024 -0.043463742 1.582182e-10
# 21          Omicron (BA.2) - Other  0.29252277 0.050391942 NA  0.193756378  0.391289162 2.284758e-06


# estimated proportion of different LINEAGES among new lab diagnoses in Belgium today
today # "2021-10-25"
multinom_preds_today_avg = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
