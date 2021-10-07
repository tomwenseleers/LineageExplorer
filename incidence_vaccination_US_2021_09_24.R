# ANALYSIS RELATIONSHIP CASES, HOSPITALISATIONS & VACCINATION COVERAGE IN THE US
# DATA NYT, https://www.nytimes.com/interactive/2021/us/covid-cases.html
# downloaded 24 SEPT 2021
# T. Wenseleers, 24 Sept 2021

# library(devtools)
# devtools::install_github("tomwenseleers/export")
library(export)

data = read.csv(".//data//us/vaccination_cases_hosps_NYT_24_sept_2021.csv", encoding="UTF-8")
data$STATE = as.factor(data$STATE)
data$CASES_DAILY = as.numeric(data$CASES_DAILY)
data$CASES_PER100K = as.numeric(data$CASES_PER100K)
data$HOSPITALISED_DAILY = as.numeric(data$HOSPITALISED_DAILY)
data$HOSPITALISED_PER100K = as.numeric(data$HOSPITALISED_PER100K)
data$POPULATION_2021 = as.numeric(data$POPULATION_2021)

data$log1pCASES_PER100K = log1p(data$CASES_PER100K)
data$log1pHOSPITALISED_PER100K = log1p(data$HOSPITALISED_PER100K)
data$log1pDEATHS_PER_100K = log1p(data$DEATHS_PER_100K)
data$logPOP = log(data$POPULATION_2021)
data$logPOPDENS = log(data$POPDENS)

data = data[complete.cases(data),] # normally not necessary

head(data)
str(data)


# regression model FULLY_VAXXED ifo vote

fit_vaxxed = lm(qlogis(FULLY_VAXXED) ~ MEDIAN_AGE+TRUMP_VOTE_2020+POPDENS, data=data)
summary(fit_vaxxed)

plot(effect(mod=fit_vaxxed, term="TRUMP_VOTE_2020", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion that voted Trump in 2020", ylab="Proportion fully vaccinated",
     axes=list(y=list(transform=list(trans=qlogis, inverse=plogis))), type="rescale",
     main="Proportion fully vaccinated i.f.o. 2020 vote in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX TRUMP VOTE.png", width=8, height=6)


# regression model log1p(DEATHS PER 100K)

fit_deaths = lm(log1pDEATHS_PER_100K ~ FULLY_VAXXED+MEDIAN_AGE, data=data)
summary(fit_deaths)

plot(effect(mod=fit_deaths, term="FULLY_VAXXED", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion fully vaccinated", ylab="Daily deaths per 100K",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Daily deaths per 100K i.f.o. prop. fully vaccinated in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX DEATHS.png", width=8, height=6)


# regression model log1p(DEATHS PER 100K) ifo vote

fit_deaths2 = lm(log1pDEATHS_PER_100K ~ MEDIAN_AGE+TRUMP_VOTE_2020, data=data)
summary(fit_deaths2)

plot(effect(mod=fit_deaths2, term="TRUMP_VOTE_2020", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion that voted Trump in 2020", ylab="Daily deaths per 100K",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Daily deaths per 100K i.f.o. 2020 vote in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX DEATHS TRUMP VOTE.png", width=8, height=6)


# regression model log1p(HOSPS PER 100K)

fit_hosps = lm(log1pHOSPITALISED_PER100K ~ FULLY_VAXXED, data=data) # adding +MEDIAN_AGE doesn't improve fit
summary(fit_hosps)

plot(effect(mod=fit_hosps, term="FULLY_VAXXED", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion fully vaccinated", ylab="Daily hospitalisations per 100K",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Daily hospitalisations per 100K i.f.o. prop. fully vaccinated in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX HOSPS.png", width=8, height=6)


# regression model log1p(HOSPS PER 100K) ifo vote

fit_hosps2 = lm(log1pHOSPITALISED_PER100K ~ MEDIAN_AGE+TRUMP_VOTE_2020, data=data)
summary(fit_hosps2)

plot(effect(mod=fit_hosps2, term="TRUMP_VOTE_2020", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion that voted Trump in 2020", ylab="Daily hospitalisations per 100K",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Daily hospitalisations per 100K i.f.o. 2020 vote in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX HOSPS TRUMP VOTE.png", width=8, height=6)


# regression model log1p(CASES PER 100K)

fit_cases = lm(log1pCASES_PER100K ~ FULLY_VAXXED, data=data)
summary(fit_cases)

plot(effect(mod=fit_cases, term="FULLY_VAXXED", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion fully vaccinated", ylab="Daily cases per 100K",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Daily cases per 100K i.f.o. prop. fully vaccinated in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX CASES.png", width=8, height=6)


# regression model log1p(CASES PER 100K) ifo vote

fit_cases2 = lm(log1pCASES_PER100K ~ MEDIAN_AGE+TRUMP_VOTE_2020, data=data)
summary(fit_cases2)

plot(effect(mod=fit_cases2, term="TRUMP_VOTE_2020", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion that voted Trump in 2020", ylab="Daily cases per 100K",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Daily cases per 100K i.f.o. 2020 vote in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX CASES TRUMP VOTE.png", width=8, height=6)


# regression model HOSPITALISED_14D_CHANGE

fit_hosps_change = lm(HOSPITALISED_14D_CHANGE ~ FULLY_VAXXED, data=data)
summary(fit_hosps_change)

plot(effect(mod=fit_hosps_change, term="FULLY_VAXXED", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion fully vaccinated", ylab="Daily hospitalisations (14D change)",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Hospitalised 14D change i.f.o. prop. fully vaccinated in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX HOSPS 14D CHANGE.png", width=8, height=6)

# regression model CASES_14D_CHANGE

fit_cases_change = lm(CASES_14D_CHANGE ~ FULLY_VAXXED, data=data)
summary(fit_cases_change)

plot(effect(mod=fit_cases_change, term="FULLY_VAXXED", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black",
     xlab="Proportion fully vaccinated", ylab="Daily cases (14D change)",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale",
     main="Cases 14D change i.f.o. prop. fully vaccinated in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on semilog scale
graph2png(file=".//plots/us/effect VAX CASES 14D CHANGE.png", width=8, height=6)


# regression model CFR

fit_CFR = lm(CFR ~ FULLY_VAXXED, data=data)
summary(fit_CFR)

plot(effect(mod=fit_CFR, term="FULLY_VAXXED", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.75), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="Proportion fully vaccinated", ylab="Deaths/cases (no lag)", 
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="response",
     main="Deaths per cases (no lag) in the US by state",
     id=list(n=nrow(data), labels=data$STATE)) # on backtransformed scale
graph2png(file=".//plots/us/effect VAX CFR.png", width=8, height=6)


