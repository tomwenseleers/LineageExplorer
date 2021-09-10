# CALCULATION OF GLOBAL DEATH TOLL CAUSED BY DIFFERENT VARIANTS 
# CALCULATED FROM GISAID LINEAGE FREQUENCIES & ANALYSIS OF TOTAL DEATH TOLL BY THE ECONOMIST

# PLUS CALCULATION OF CUMULATIVE TOTAL OF PREVIOUSLY INFECTED PEOPLE

# T. Wenseleers, 10 Sept. 2021

fit_economist = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/output-for-interactive/second_map.csv")
OWID = read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv")
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

# estimated global % that was already infected by the virus on 8 Sept 2021
# 57% [25-67%]
weighted.mean(fit_economist$implied_infections_per_100_persons_clipped,
              w=fit_economist$population, na.rm=T)
weighted.mean(fit_economist$implied_infections_per_100_persons_bot_95_clipped,
              w=fit_economist$population, na.rm=T)
weighted.mean(fit_economist$implied_infections_per_100_persons_top_95_clipped,
              w=fit_economist$population, na.rm=T)

# estimated % already infected by the virus on 8 Sept 2021
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



