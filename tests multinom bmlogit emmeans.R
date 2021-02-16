library(nnet)
library(mclogit)
library(broom)

# multinomial spline fit on test outcome data (negative / positive wild type / positive B.1.1.7
# to be able to estimate growth rate and Rt of B.1.1.7 and wild type separately

data_ag_long = read.csv("https://github.com/tomwenseleers/newcovid_belgium/raw/main/temp/data_ag_long.csv")
data_ag_long$outcome = factor(data_ag_long$outcome, levels=c("n_neg","n_spos","n_b117"))
data_ag_long$LABORATORY = factor(data_ag_long$LABORATORY)
data_ag_long$collection_date = as.Date(data_ag_long$collection_date)
set.seed(1)
date.from = as.Date("2021-01-14")
data_ag_long_subs = data_ag_long[(data_ag_long$collection_date>=date.from),]

# example nnet::multinom & mclogit::mblogit multinomial fit
mfit0 = nnet::multinom(outcome ~ scale(collection_date_num) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mblogitfit0 = mblogit(outcome ~ scale(collection_date_num) + LABORATORY, weights=count, data=data_ag_long_subs)

# mutlinomial coefs & conf intervals in function of time (collection_date_nu) of both models are
mfit0_coefs = data.frame(tidy(mfit0, conf.int=TRUE))[,c("term","estimate","conf.low","conf.high")]
rownames(mfit0_coefs) = paste0(rep(mfit0$lab[-1],each=nrow(mfit0_coefs)/length(mfit0$lab[-1])),"~",mfit0_coefs$term) 
mfit0_coefs[,c("term")] = NULL
mfit0_coefs = mfit0_coefs[grepl("date",rownames(mfit0_coefs)),]
mfit0_coefs = mfit0_coefs / attr(scale(data_ag_long_subs$collection_date_num),"scaled:scale")
mfit0_coefs
#                                      estimate    conf.low   conf.high
# n_spos~scale(collection_date_num) -0.03388852 -0.03642875 -0.03134829
# n_b117~scale(collection_date_num)  0.01453788  0.01036944  0.01870632

mblogitfit_coefs = data.frame(tidy(mblogitfit0, conf.int=TRUE))[,c("term","estimate","conf.low","conf.high")]
rownames(mblogitfit_coefs) = mblogitfit_coefs$term
mblogitfit_coefs[,c("term")] = NULL
mblogitfit_coefs = mblogitfit_coefs[grepl("date",rownames(mblogitfit_coefs)),]
mblogitfit_coefs = mblogitfit_coefs / attr(scale(data_ag_long_subs$collection_date_num),"scaled:scale")
mblogitfit_coefs
#                                    estimate    conf.low   conf.high
# n_spos~scale(collection_date_num) -0.03388029 -0.03642333 -0.03133724
# n_b117~scale(collection_date_num)  0.01454377  0.01037025  0.01871729

# these match, so both models agree, as should be the case


# the coefficients & confidence intervals above should match with 
# confint(emtrends(mfit0, trt.vs.ctrl~outcome|1, var="collection_date_num",  mode="latent", adjust="none"))$contrasts
# test: this is indeed the case for the development version (there was a bug in version 1.5.4 though)
# the emmeans result gives nonsensical results for the mblogit model though:

unloadNamespace("emmeans")
require(devtools)
remotes::install_github("rvlenth/emmeans")
library(emmeans)

confint(emtrends(mfit0, trt.vs.ctrl~outcome|1, var="collection_date_num",  mode="latent", adjust="none"))$contrasts
# contrast       estimate      SE df lower.CL upper.CL
# n_spos - n_neg  -0.0339 0.00130 16  -0.0366  -0.0311
# n_b117 - n_neg   0.0145 0.00213 16   0.0100   0.0190
# correct

confint(emtrends(mblogitfit0, trt.vs.ctrl~outcome|1, var="collection_date_num",  mode="latent", adjust="none"))$contrasts
# contrast       estimate    SE  df asymp.LCL asymp.UCL
# n_spos - n_neg   -19.64 1.37 Inf    -22.33     -16.9
# n_b117 - n_neg     8.43 2.80 Inf      2.94      13.9
# these numbers do not make sense - they have a totally different range than the ones above
# Russ Length is currently diagnosing the problem & making a bug fix to the emmeans and/or the mclogit package

rg = ref_grid(mfit0, mode = "latent")
RG = ref_grid(mblogitfit0, mode = "latent")
RG@linfct = rg@linfct # problem is that there are NAs in RG@linfct due to the terms() result not containing the "predvars" attribute needed by scale()
summary(RG)