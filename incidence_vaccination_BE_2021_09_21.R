# ANALYSIS RELATIONSHIP INCIDENCE & VACCINATION COVERAGE IN BELGIUM
# T. Wenseleers, 23 Sept 2021

# data provided by Dries de Smet:
# NONVAC = % not vaccinated
# INC = cumulative incidence the last 14 days per 100.000 inhabitants
# INC2020 = cumulative incidence in 2020 per 100.000 inhabitants
# TAXINC = mean net income per inhabitant (in 2018)
# POPDENS = population density
# RAND = 19 municipalities in Vlaams rand
# PCTMIGR = percentage inhabitants that are not Belgian or with non-Belgian background
# PCTNOTBORNINB = percentage inhabitants not born in Belgium
# CUMULCASES = total cumulative nr of infections since March 1 2020

library(devtools)
# devtools::install_github("tomwenseleers/export")
library(export)

data = read.csv(".//data//incidence_vaccination_BE_2021_09_21.csv", encoding="UTF-8")
data$Gemeente = iconv(data$Gemeente, from="ISO_8859-2", to="UTF-8")
# data$PROV = iconv(data$PROV, from="UTF-8", to="ISO_8859-2")
data$Gemeente = as.factor(data$Gemeente)
data$Regio = as.factor(data$Regio)
data$PROV = as.factor(data$PROV)
data$POP = as.numeric(data$POP)
data$TAXINC = as.numeric(data$TAXINC)
data$CUMULCASES = as.numeric(data$CUMULCASES)
data$INC2020 = as.numeric(data$INC2020) # PS last value was NA - fix?

head(data)
str(data)

data = data[complete.cases(data),] # remove lines with NAs

data$log1pINC = log1p(data$INC)
data$log1pINC2020 = log1p(data$INC2020)
data$log1pTAXINC = log1p(data$TAXINC)
data$logPOP = log(data$POP)
data$logAGE = log(data$MEAN_AGE)
data$sqrtAGE = sqrt(data$MEAN_AGE)
data$expAGE = exp(data$MEAN_AGE)



# multiple regression model for incidence INC
library(bestglm)
bestglmfit = bestglm(Xy=data[,c("NONVAC","log1pTAXINC","logPOP","log1pINC2020","PCTMIGR","MEAN_AGE","log1pINC")], family=gaussian, IC="AIC")
# PS using glmulti is another option, that one also allows all possible 1st order interaction effects to be considered, see https://rstudio-pubs-static.s3.amazonaws.com/2897_9220b21cfc0c43a396ff9abf122bb351.html for an example
bestfit = bestglmfit$BestModel
summary(bestfit) # best model based on AIC value

bestfit = lm(log1p(INC) ~ NONVAC+log(POP)+log1p(TAXINC)+log1p(INC2020), data=data)
summary(bestfit)
# Coefficients:
#         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    -15.87078    2.42714  -6.539 1.39e-10 ***
#         NONVAC           6.37076    0.48102  13.244  < 2e-16 ***
#         log(POP)         0.13610    0.03052   4.460 9.90e-06 ***
#         log1p(TAXINC)    1.47726    0.22432   6.585 1.04e-10 ***
#         log1p(INC2020)   0.39220    0.08576   4.573 5.90e-06 ***
#         ---
#         Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5838 on 564 degrees of freedom
# Multiple R-squared:  0.4745,	Adjusted R-squared:  0.4708 
# F-statistic: 127.3 on 4 and 564 DF,  p-value: < 2.2e-16


library(car)
library(MASS)
vif(bestfit) # variance inflation factors not crazy high (none >10), so collinearity OK
# NONVAC       log(POP)  log1p(TAXINC) log1p(INC2020) 
# 2.640562       1.089854       1.643984       2.071303 

# effect plots (partial effects of each covariate, with partial residuals)
# for all terms
plot(allEffects(bestfit, residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.2), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black")
graph2ppt(file=".//plots/BE/all effects.pptx", width=8, height=6)
graph2pdf(file=".//plots/BE/all effects.pdf", width=8, height=6)
graph2png(file=".//plots/BE/all effects.png", width=8, height=6)


# for NONVAC
plot(effect(mod=bestfit, term="NONVAC", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="% non-vaccinated", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="response") # on backtransformed scale
plot(effect(mod=bestfit, term="NONVAC", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="% non-vaccinated", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale") # on semilog scale
graph2ppt(file=".//plots/BE/effect NONVAC.pptx", width=8, height=6)
graph2pdf(file=".//plots/BE/effect NONVAC.pdf", width=8, height=6)
graph2png(file=".//plots/BE/effect NONVAC.png", width=8, height=6)

# for log(POP)
plot(effect(mod=bestfit, term="log(POP)", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="Population size", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="response") # on backtransformed scale
plot(effect(mod=bestfit, term="log(POP)", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="Population size", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale") # on semilog scale
graph2ppt(file=".//plots/BE/effect POP.pptx", width=8, height=6)
graph2pdf(file=".//plots/BE/effect POP.pdf", width=8, height=6)
graph2png(file=".//plots/BE/effect POP.png", width=8, height=6)

# for log1p(TAXINC)
plot(effect(mod=bestfit, term="log1p(TAXINC)", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="Mean income per inhabitant", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="response") # on backtransformed scale
plot(effect(mod=bestfit, term="log1p(TAXINC)", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="Mean income per inhabitant", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale") # on semilog scale
graph2ppt(file=".//plots/BE/effect TAXINC.pptx", width=8, height=6)
graph2pdf(file=".//plots/BE/effect TAXINC.pdf", width=8, height=6)
graph2png(file=".//plots/BE/effect TAXINC.png", width=8, height=6)

# for log1p(INC2020)
plot(effect(mod=bestfit, term="log1p(INC2020)", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="Mean incidence in 2020", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="response") # on backtransformed scale
plot(effect(mod=bestfit, term="log1p(INC2020)", residuals=T), residuals.pch=16, residuals.color=alpha("steelblue", 0.3), smooth.residuals=FALSE, 
     ci.style="bands", band.colors="black", band.transparency=0.2, colors="black", 
     xlab="Mean incidence in 2020", ylab="Incidence over last 2 weeks", main="",
     axes=list(y=list(transform=list(trans=log1p, inverse=function(yt) exp(yt)-1))), type="rescale") # on semilog scale
graph2ppt(file=".//plots/BE/effect INC2020.pptx", width=8, height=6)
graph2pdf(file=".//plots/BE/effect INC2020.pdf", width=8, height=6)
graph2png(file=".//plots/BE/effect INC2020.png", width=8, height=6)




