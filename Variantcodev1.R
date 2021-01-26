# Data VOC Emmanuel André
# Analysis 
# Jan 2021
############################
library(ggplot2)
library(scales)
library(EpiEstim)
library(lme4)

# setwd("/Users/Niel/Dropbox (DSI)/Niel/CENSTAT/Projects/nCov2019/DataVariant")

# Read in testdata
testdata=read.csv(".//data//all valid PCRs since 1st of january.csv")
head(testdata)
testdata$date=as.Date(testdata$Analysis.created.at..UTC.)

# Read in sdropdata
sdropdata=read.csv(".//data//S dropouts raw data for Niel Hens.csv")
head(sdropdata)
sdropdata$date=as.Date(sdropdata$Analysis.created.at..UTC.)

# select subset based on last date
sel=min(max(sdropdata$date),max(testdata$date))
sdropdata=subset(sdropdata,date<=sel)
testdata=subset(testdata,date<=sel)

# Check sdrop being part of testdata
table(sdropdata$Sample.ID%in%testdata$Sample.ID)

# Travellers: ULB & Namur until 8 Jan hence look at last 14 days
# Exclude Antwerp (maybe Gent): pro-active for outbreaks

# sdrop only analysis
######################
# Growth rate
tab=table(sdropdata$date)
attributes(tab)

delta=14
tmp=as.vector(tab[(length(tab)-delta+1):length(tab)])
fit=glm(tmp~c(1:delta),quasipoisson(link="log"))
summary(fit)

plot(c(1:delta),tmp,xlab="last x days",ylab="number of S gene dropouts")
lines(c(1:delta),predict(fit,type="response"))

# growth rate = log(2)/Td
# doubling time (take care of the order which can change)
log(2)/c(fit$coeff[2],rev(confint(fit)[2,]))
Gr.s=fit$coeff[2]

# R-estimate
res=estimate_R(incid=tmp,method="parametric_si",config = make_config(list(mean_si = 3.6, std_si = 3.1)))
res=estimate_R(incid=tmp,method="parametric_si",config = make_config(list(mean_si = 5.5, std_si = 3.1)))
res$R
R.s=res$R

# Conclusion depends on settings:
# Based on the last 14 days; we see an increase in S-gene dropouts. 
# Merely analysing the number of S-gene dropouts we inferred a doubling time of 3.7 days (95% CI: 2.8,4.8).
# Estimating the reproduction number on this limited time span, gives Rt(S gene) 1.77 (95% CI: 1.62,1.93)
# *note that testing bias could influence these numbers; careful interpretation is needed

# Looking at hospitals as fixed effects; note GOF hasn't been checked
tab.factors=table(sdropdata$Laboratory,as.Date(sdropdata$Analysis.created.at..UTC.))
tab.factors

facdata=as.data.frame(tab.factors)
names(facdata)=c("hosp","date","freq")
fitglm=glm(freq~hosp+as.numeric(date):hosp,data=facdata[as.Date(facdata$date)>="2021-01-10",],family=quasipoisson(link = "log"))
summary(fitglm)

# non s-drop only analysis
###########################
# Growth rate
tab=table(testdata$date[!(testdata$Sample.ID%in%sdropdata$Sample.ID)&(testdata$Outcome=="Detected"|testdata$Outcome=="Positive")])
attributes(tab)

delta=14
tmp=as.vector(tab[(length(tab)-delta+1):length(tab)])
fit=glm(tmp~c(1:delta),quasipoisson(link="log"))
summary(fit)
Gr.wt=fit$coeff[2]

plot(c(1:delta),tmp,xlab="last x days",ylab="number of tests")
lines(c(1:delta),predict(fit,type="response"))

# growth rate = log(2)/Td
# doubling time (take care of the order which can change)
log(2)/c(fit$coeff[2],(confint(fit)[2,]))

# R-estimate
res=estimate_R(incid=tmp,method="parametric_si",config = make_config(list(mean_si = 3.6, std_si = 3.1)))
res=estimate_R(incid=tmp,method="parametric_si",config = make_config(list(mean_si = 5.5, std_si = 3.1)))
res$R
R.wt=res$R

# Selection rate between the s-drop and wild-type variant
##########################################################
s=Gr.s-Gr.wt
sT=s*c(3.6,5.5)
exp(sT)

# Ratio based on Rt's estimated using Cori et al. (2013) - check which gen time used
R.s$"Mean(R)"[dim(R.s)[1]] / R.wt$"Mean(R)"[dim(R.s)[1]]

