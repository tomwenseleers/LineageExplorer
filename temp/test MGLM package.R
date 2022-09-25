# TESTS MGLM PACKAGE ####

softMax <- function(eta){
  exp_eta <- exp(eta)
  return(sweep(exp_eta, 1, STATS=rowSums(exp_eta), FUN="/"))
}
library(splines)

# dataset: SARS-CoV2 variant frequencies in function of time ####
dat = read.csv("https://www.dropbox.com/s/u27cn44p5srievq/dat.csv?dl=1")
dat$variant = factor(dat$variant)
variant = dat$variant
count = dat$count
# dat$variant = relevel(dat$variant, ref="Other")
# observed counts:
Y = dat %>%
  pivot_wider(id_cols ="collection_date_num", names_from = "variant", values_from = "count")
Y[is.na(Y)] = 0
collection_date_num = Y$collection_date_num
Y = as.matrix(Y[,-which(colnames(Y)==c("collection_date_num"))])
# observed proportions:
Y_prop = (dat %>% 
            pivot_wider(id_cols ="collection_date_num", names_from = "variant", values_from = "prop"))[,-1]
Y_prop[is.na(Y_prop)] = 0
# model matrix:
X = model.matrix(formula(~ ns(collection_date_num, df=2)))
p = ncol(X)

# nnet::multinom fits OK ####
library(nnet)
set.seed(1)
fit_multinom = nnet::multinom(variant ~ ns(collection_date_num, df=2), 
                              weights=count, data=dat)
preds = softMax(X %*% t(rbind(0,coef(fit_multinom))))  # = fitted(fit_multinom)

matplot(preds, type="l", col=1:(1+ncol(Y)))
matplot(Y_prop, type="p", add=T, col=1:(1+ncol(Y)), pch=16)

# mblogit fit: return error ####
library(mclogit)
fit_mblogit = mblogit(formula=variant ~ ns(collection_date_num, df=2),
                      weights=count,
                      data=dat,
                      from.table=FALSE, dispersion=FALSE,
                      control=mclogit.control(maxit=100))
# returns error Error: no valid set of coefficients has been found: please supply starting values


# MGLM fits: return errors ####
library(MGLM)
fit_MGLM_MN = MGLMreg(Y ~ ns(collection_date_num, df=2), 
                      dist="MN", LRT=FALSE, maxiters = 1000) # multinomial
# gives error
# Error in if (is.nan(ll.Newton) | ll.Newton >= 0) ll.Newton <- NA : 
# missing value where TRUE/FALSE needed

fit_MGLM_DM = MGLMreg(Y ~ ns(collection_date_num, df=2), # Dirichlet multinomial
                      dist="DM", LRT=FALSE, maxiters = 1000)
# gives error
# Error in if (is.nan(ll.Newton) | ll.Newton >= 0) ll.Newton <- NA : 
# missing value where TRUE/FALSE needed
fit_MGLM_NegMN = MGLMreg(Y ~ ns(collection_date_num, df=2), # negative multinomial
                         dist="NegMN", LRT=FALSE, maxiters = 1000)
# gives error
# Error: inner loop 1; cannot correct step size
fit_MGLM_GDM = MGLMreg(Y ~ ns(collection_date_num, df=2), # generalized Dirichlet multinomial
                       dist="GDM", LRT=FALSE, maxiters = 1000)
# gives error
# Error in if (mean(gradients^2) > 1e-04) { : 
# missing value where TRUE/FALSE needed


# MGLM sparse fits: only 2 without errors, but bad fit ####
fit_MGLM_MN_sparse = MGLMsparsereg(Y ~ ns(collection_date_num, df=2), # multinomial
                                   dist="MN", maxiters = 1000, penidx=c(FALSE,rep(TRUE,p-1)), 
                                   lambda=1E-5,
                                   penalty="sweep")
# Error in B[b > lambda] <- B[b > lambda] - lambda : 
#  NAs are not allowed in subscripted assignments

fit_MGLM_DM_sparse = MGLMsparsereg(Y ~ ns(collection_date_num, df=2), # Dirichlet multinomial
                                   dist="DM", maxiters = 1000, penidx=c(FALSE,rep(TRUE,p-1)),
                                   lambda=1E-5,
                                   penalty="sweep")
# no errors, but bad fit to data:
preds = softMax(X %*% coef(fit_MGLM_DM_sparse))
matplot(preds, type="l", col=1:(1+ncol(Y)))
matplot(Y_prop, type="p", add=T, col=1:(1+ncol(Y)), pch=16)

fit_MGLM_NegMN_sparse = MGLMsparsereg(Y ~ 0+X, # Dirichlet multinomial
                                      dist="NegMN", maxiters = 1000, penidx=c(FALSE,rep(TRUE,p-1)),
                                      lambda=1E-5,
                                      penalty="sweep")
# Error in B[b > lambda] <- B[b > lambda] - lambda : 
#   NAs are not allowed in subscripted assignments
fit_MGLM_GDM_sparse = MGLMsparsereg(Y ~ 0+X, # Dirichlet multinomial
                                    dist="GDM", maxiters = 1000, penidx=c(FALSE,rep(TRUE,p-1)),
                                    lambda=1E-5,
                                    penalty="sweep")
# no errors, but bad fit to data:
preds = softMax(X %*% coef(fit_MGLM_GDM_sparse))
matplot(preds, type="l", col=1:(1+p))
matplot(Y_prop, type="p", add=T, col=1:(1+p), pch=16) 

# predict method for MGLMsparsereg is also missing:
predict(fit_MGLM_DM_sparse, newdata=X)

# as are vcov methods for MGLMreg & MGLMsparsereg objects

# any chance that those could be provided?

# I am asking this in the context of hoping to define methods for
# emmeans and marginaleffects, so that MGLMreg & MGLMsparsereg objects could be supported there
# (now these packages only support nnet::multinom & mclogit::mblogit objects)
# see
# https://github.com/rvlenth/emmeans/blob/master/R/multinom-support.R
# and
# https://github.com/tomwenseleers/marginaleffects/blob/main/R/methods_nnet.R
