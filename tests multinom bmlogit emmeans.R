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
mfit0 = nnet::multinom(outcome ~ scale(collection_date_num, center=TRUE, scale=FALSE) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mblogitfit0 = mblogit(outcome ~ scale(collection_date_num, center=TRUE, scale=FALSE) + LABORATORY, weights=count, data=data_ag_long_subs)
attr(mblogitfit0$terms, "predvars") = attr(mfit0$terms, "predvars")

# mutlinomial coefs & conf intervals in function of time (collection_date_nu) of both models are
mfit0_coefs = data.frame(tidy(mfit0, conf.int=TRUE))[,c("term","estimate","conf.low","conf.high")]
rownames(mfit0_coefs) = paste0(rep(mfit0$lab[-1],each=nrow(mfit0_coefs)/length(mfit0$lab[-1])),"~",mfit0_coefs$term) 
mfit0_coefs[,c("term")] = NULL
mfit0_coefs = mfit0_coefs[grepl("date",rownames(mfit0_coefs)),]
mfit0_coefs = mfit0_coefs
mfit0_coefs
#                                                                    estimate    conf.low   conf.high
# n_spos~scale(collection_date_num, center = TRUE, scale = FALSE) -0.03389104 -0.03643120 -0.03135088
# n_b117~scale(collection_date_num, center = TRUE, scale = FALSE)  0.01454299  0.01037375  0.01871222

mblogitfit_coefs = data.frame(tidy(mblogitfit0, conf.int=TRUE))[,c("term","estimate","conf.low","conf.high")]
rownames(mblogitfit_coefs) = mblogitfit_coefs$term
mblogitfit_coefs[,c("term")] = NULL
mblogitfit_coefs = mblogitfit_coefs[grepl("date",rownames(mblogitfit_coefs)),]
mblogitfit_coefs = mblogitfit_coefs
mblogitfit_coefs
#                                    estimate    conf.low   conf.high
# n_spos~scale(collection_date_num) -0.03388029 -0.03642333 -0.03133724
# n_b117~scale(collection_date_num)  0.01454377  0.01037025  0.01871729

# these match, so both models agree, as should be the case


# the coefficients & confidence intervals above should match with 
# confint(emtrends(mfit0, trt.vs.ctrl~outcome|1, var="collection_date_num",  mode="latent"))$contrasts

unloadNamespace("emmeans")
require(devtools)
remotes::install_github("rvlenth/emmeans", force=TRUE)
library(emmeans)

confint(emtrends(mfit0, trt.vs.ctrl~outcome|1, var="collection_date_num",  mode="latent"))$contrasts
# contrast       estimate      SE df lower.CL upper.CL
# n_spos - n_neg  -0.0339 0.00130 16  -0.0366  -0.0311
# n_b117 - n_neg   0.0145 0.00213 16   0.0100   0.0190
# correct

# these confidence intervals are not correct though:
confint(emtrends(mblogitfit0, trt.vs.ctrl~outcome|1, var="collection_date_num",  mode="latent", adjust="none"))$contrasts
# contrast       estimate     SE  df asymp.LCL asymp.UCL
# n_spos - n_neg  -0.0339 0.0587 Inf    -0.149    0.0812
# n_b117 - n_neg   0.0145 0.0900 Inf    -0.162    0.1909
