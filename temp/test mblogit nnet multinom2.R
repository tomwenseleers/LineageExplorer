library(MASS)
library(marginaleffects)
library(nnet)
library(mclogit)
library(memisc)
library(mclogit)
library(insight)
library(dplyr)
library(tidyr)
library(emmeans)
library(splines)
library(MGLM)

# 1. Dataset: covid variant frequencies ("variant") through time ("collection_date_num") ####
# "count" = actual count per week
datfull = read.csv("https://www.dropbox.com/s/u27cn44p5srievq/dat.csv?dl=1")
# simplify example to just 3 categories:
# inverse_softMax <- function(mu) {
#   log_mu <- log(mu)
#   # we let the log(odds) sum to zero
#   # these predictions are referred to as type="latent" in the emmeans package
#   return(sweep(log_mu, 1, STATS=rowMeans(log_mu), FUN="-")) 
# }
dat = datfull[datfull$variant %in% c("Alpha", "Beta", "Other"),] 
dat$variant = factor(dat$variant, levels=c("Other", "Beta", "Alpha"))
datfull$variant = factor(datfull$variant)
datfull$variant = relevel(datfull$variant, ref="Other")


# # 2. fixed get_predict methods for nnet::multinom and mclogit::mblogit ####
# 
# # softMax <- function(eta){
# #   exp_eta <- exp(eta)
# #   return(sweep(exp_eta, 1, STATS=rowSums(exp_eta), FUN="/"))
# # }
# # this is the link function for multinomial
# # = generalized logit = inverse softmax (you could also call it glogit)
# # here assuming mu is a matrix with predictions on the response = probability scale
# # for ALL outcome levels as columns (including for the reference level)
# 
# # this is the inverse link function for multinomial 
# # = inverse generalized logit = softmax (you could also call it inv_glogit)
# # here assuming eta is a matrix with predictions on the latent scale
# # for ALL outcome levels as columns (including for the reference level)
# 
# # these allow prediction with either type="response" or type="latent"
# # and for both prediction types return predictions for all
# # outcome levels, including the reference level
# # (as is done in the emmeans package) - I don't think anything else
# # makes sense
# # default I think should be type="latent"
# 
# # pr_multinom <- function(model, newdata, type=c("latent", "link", "probs", "response")) # probs=response, latent=link (latent is how it is referred to in the emmeans package)
# # {
# #   type <- match.arg(type)
# #   if ("mblogit" %in% class(model)) { resptype="response" } else { resptype="probs" }
# #   if (type == "probs"|type == "response") return(predict(model, newdata, type=resptype))
# #   
# #   mu <- predict(model, newdata, type=resptype)
# #   return(inverse_softMax(mu)) # when type=="latent" | type=="link"
# #   
# # }
# # 
# # get_predict.multinom <- function(model,
# #                                  newdata = insight::get_data(model),
# #                                  vcov = FALSE,
# #                                  conf_level = 0.95,
# #                                  type = "latent",
# #                                  ...) {
# #   
# #   # type <- sanitize_type(model, type) # type_dictionary still needs to be updated as well, just unmarked this to test
# #   
# #   pred <- pr_multinom(model,
# #                   newdata = newdata,
# #                   type = type,
# #                   ...)
# #   
# #   # atomic vector means there is only one row in `newdata`
# #   if (isTRUE(checkmate::check_atomic_vector(pred))) {
# #     pred <- matrix(pred, nrow = 1, dimnames = list(NULL, names(pred)))
# #   }
# #   
# #   # matrix with outcome levels as columns
# #   out <- data.frame(
# #     group = rep(colnames(pred), each = nrow(pred)),
# #     predicted = c(pred))
# #   
# #   # usually when `newdata` is supplied by `comparisons`
# #   if ("rowid" %in% colnames(newdata)) {
# #     out$rowid <- rep(newdata$rowid, times = ncol(pred))
# #   } else {
# #     out$rowid <- rep(seq_len(nrow(pred)), times = ncol(pred))
# #   }
# #   
# #   return(out)
# # }
# # 
# # get_predict.mblogit <- function(model,
# #                                 newdata = insight::get_data(model),
# #                                 vcov = FALSE,
# #                                 conf_level = 0.95,
# #                                 type = "latent",
# #                                 ...) {
# #   
# #   if (!isTRUE(checkmate::check_flag(vcov, null.ok = TRUE)) &&
# #       isTRUE(list(...)$calling_function == "predictions")) {
# #     stop("The `vcov` argument is not supported for this model class.", call. = FALSE)
# #   }
# #   
# #   out <- suppressMessages(
# #     get_predict.multinom(model = model,
# #                          newdata = newdata,
# #                          type = type,
# #                          ...))
# #   return(out)
# # }


# 2. mclogit::mblogit multinomial fit ####
fit_mblogit = mblogit(formula=variant ~ ns(collection_date_num, df=2),
                      weights=count,
                      data=dat,
                      from.table=FALSE, dispersion=FALSE,
                      control=mclogit.control(maxit=100))

# test of mclogit::mblogit model prediction for row 160 in dataset
# prediction on response scale
mclogit:::predict.mblogit(fit_mblogit, type="response")[160,]
#      Other       Beta      Alpha 
# 0.86616413 0.01788666 0.11594921 

# prediction on link scale (implicitly here reference level Other=0, but dropped, better maybe to add explicit column Other=0 here?)
mclogit:::predict.mblogit(fit_mblogit, type="link")[160,]
#      Beta     Alpha 
# -3.880019 -2.010922 

# standard errors on link scale
mclogit:::predict.mblogit(fit_mblogit, type="link", se.fit=TRUE)$se.fit[160,]
#       Beta      Alpha 
# 0.15458046 0.06828786

# emmeans result (response/prob scale, with SEs and 95% CLs on response scale)
preds_emmeans_response = as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
        at=list(collection_date_num=unique(dat$collection_date_num)),
        mode="prob", level=0.95))
# confidence intervals go slightly outside allowed domain [0,1], which is undesirable
range(preds_emmeans_response$asymp.LCL) # -0.0006024099  1.0000000000
range(preds_emmeans_response$asymp.UCL) # 1.918237e-27 1.000169e+00

# emmeans result (latent scale, with SEs and 95% CLs on latent scale)
preds_emmeans_latentbacktransformed = as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
                      at=list(collection_date_num=unique(dat$collection_date_num)),
                      mode="latent", level=0.95), type="response")
# with type="response" it just takes exp(y_on_latent_scale) as backtransform,
# whereas correct backtransform would have been exp(y_on_latent_scale)/sum(exp(y_on_latent_scale))
preds_emmeans_latentbacktransformed[1,]
#   collection_date_num variant     e^y      SE  df asymp.LCL asymp.UCL
# 1               18414   Other 1949560 1860574 Inf  300322.5  12655681

# the correct backtransformed result would have been :
library(dplyr)
softmax <- function(x) exp(x) / sum(exp(x)) # 
preds_emmeans_latent = as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
                                                            at=list(collection_date_num=unique(dat$collection_date_num)),
                                                            mode="latent", level=0.95), type="latent")
preds_emmeans_latent_backtransformed = preds_emmeans_latent |>
                group_by(collection_date_num ) |> # softmax backtransform should be done grouped over outcome levels
                mutate_at(c("emmean", "asymp.LCL", "asymp.UCL"), softmax)
range(preds_emmeans_latent_backtransformed$asymp.LCL) # 1.073096e-28 1.000000e+00
range(preds_emmeans_latent_backtransformed$asymp.UCL) # 6.552742e-28 1.000000e+00

# the softmax backtransform of predictions on the latent scale to the response
# scale with the different outcome levels in different columns would be
softMax <- function(eta){ # if eta has different outcome levels in different columns
  exp_eta <- exp(eta)
  return(sweep(exp_eta, 1, STATS=rowSums(exp_eta), FUN="/"))
}









# predictions modified get_predict function for row 160: correct
# on response scale
preds_mblogit_response = marginaleffects:::get_predict.multinom(fit_mblogit, type="response")
preds_mblogit_response[preds_mblogit_response$rowid==160,]
#     group  predicted rowid
# 160 Other 0.86616413   160
# 487  Beta 0.01788666   160
# 814 Alpha 0.11594921   160

# on latent scale
preds_mblogit_latent = marginaleffects:::get_predict.multinom(fit_mblogit, type="clr")
preds_mblogit_latent[preds_mblogit_latent$rowid==160,]
#     group   predicted rowid
# 160 Other  1.96364709   160
# 487  Beta -1.91637203   160
# 814 Alpha -0.04727506   160


# test with new marginaleffects development version:
# marginaleffects predictions with type="response" :
predictions(fit_mblogit, 
            newdata = datagrid(collection_date_num=dat[160,'collection_date_num']),
            type="response" 
             ) 
# transform_post = function (eta) softMax(eta), NOTE: insight::link_inverse(fit_multinom) NOT CORRECT
# rowid     type group  predicted   std.error  statistic      p.value   conf.low
# 1     1 response Other 0.86616413 0.007360510 117.677182 0.000000e+00 0.85173779
# 2     1 response  Beta 0.01788666 0.002714336   6.589701 4.407117e-11 0.01256666
# 3     1 response Alpha 0.11594921 0.006996554  16.572332 1.104583e-61 0.10223622
# conf.high    count collection_date_num variant
# 1 0.88059046 73.10398               18764   Alpha
# 2 0.02320666 73.10398               18764   Alpha
# 3 0.12966221 73.10398               18764   Alpha

# this closely matches output of emmeans call (slight differences maybe
# due to difference in assumed df?, taken as Inf in emmeans)
as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
                      at=list(collection_date_num=dat[160,'collection_date_num']),
                      mode="prob", level=0.95))
#   collection_date_num variant       prob          SE  df  asymp.LCL  asymp.UCL
# 1               18764   Other 0.86616413 0.007357713 Inf 0.85174327 0.88058498
# 2               18764    Beta 0.01788666 0.002713279 Inf 0.01256873 0.02320459
# 3               18764   Alpha 0.11594921 0.006993861 Inf 0.10224150 0.12965693
# NOTE: these confidence intervals, which I would think cannot be correct

# TIMINGS FOR LARGER GRID

# timings adjusted predictions emmeans vs marginaleffects
dates = seq(min(dat$collection_date_num), max(dat$collection_date_num)) 
system.time(as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
                      at=list(collection_date_num=dates),
                      mode="prob", level=0.95))) # emmeans response scale 22s
system.time(as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
                                  at=list(collection_date_num=dates),
                                  mode="latent", level=0.95))) # emmeans latent scale 0.19s
system.time(predictions(fit_mblogit, 
            newdata = datagrid(collection_date_num=dates),
            type="response")) # marginaleffects response scale 0.03s
system.time(predictions(fit_mblogit, 
                        newdata = datagrid(collection_date_num=dates),
                        type="clr")) # marginaleffects clr=latent scale 0.09s

# timings emtrends vs emtrends in marginaleffects ####
library(microbenchmark)
microbenchmark(as.data.frame(confint(pairs(emtrends(fit_mblogit, ~ variant|collection_date_num,
                                     var = "collection_date_num",  
                                     mode = "latent",
                                     at = list(collection_date_num = 
                                                 dat[160,'collection_date_num'])),
                            reverse=TRUE), level=0.9))) # emtrends 83 ms

microbenchmark(comparisons(
  fit_mblogit,
  newdata = datagrid(collection_date_num = 
                       c(dat[160,'collection_date_num'])),
  variables = "collection_date_num",
  type = "clr",
  hypothesis = "pairwise")) # marginaleffects 83 ms


# all comparisons against the reference level (Other) can be obtained using type="link"="alr"
microbenchmark(marginaleffects(fit_mblogit, 
                           type = "link", 
                           variables = c("collection_date_num"),
                           by = c("group", "collection_date_num"), 
                           vcov = vcov(fit_mblogit),
                           newdata = datagrid(collection_date_num=dat[160,'collection_date_num']))) # 108 ms
# correct for Beta and Alpha vs Other
# TO DO: adapt type="link" to allow additional ref argument with ref = nr of reference group?









# marginaleffects predictions with type="clr" :
predictions(fit_mblogit, 
            newdata = datagrid(collection_date_num=dat[160,'collection_date_num']),
            type="clr" 
)
# rowid   type group   predicted  std.error   statistic       p.value   conf.low
# 1     1 latent Other  1.96364709 0.05718955  34.3357663 2.297228e-258  1.8515576
# 2     1 latent  Beta -1.91637203 0.10460975 -18.3192493  5.810889e-75 -2.1214034
# 3     1 latent Alpha -0.04727506 0.06732437  -0.7021983  4.825555e-01 -0.1792284
# conf.high    count collection_date_num variant
# 1  2.07573655 73.10398               18764   Alpha
# 2 -1.71134070 73.10398               18764   Alpha
# 3  0.08467828 73.10398               18764   Alpha

# predictions & conf intervals on backtransformed response scale would be
library(dplyr)
softmax <- function(x) exp(x) / sum(exp(x)) 
as.data.frame(predictions(fit_mblogit, type = "clr") |>
  group_by(rowid) |>
  mutate_at(c("predicted", "conf.low", "conf.high"), softmax) |>
  dplyr::filter(rowid == 160))

# rowid   type group  predicted  std.error   statistic       p.value   conf.low
# 1   160 latent Other 0.86616413 0.05718955  34.3357663 2.297228e-258 0.86952744
# 2   160 latent  Beta 0.01788666 0.10460975 -18.3192493  5.810889e-75 0.01636245
# 3   160 latent Alpha 0.11594921 0.06732437  -0.7021983  4.825555e-01 0.11411011
# conf.high variant collection_date_num count (weights)
# 1 0.86265454    Beta               18764     0         0
# 2 0.01954925    Beta               18764     0         0
# 3 0.11779620    Beta               18764     0         0

# this matches result of
library(dplyr)
softmax <- function(x) exp(x) / sum(exp(x)) 
as.data.frame(as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
                                             at=list(collection_date_num=dat[160,'collection_date_num']),
                                             mode="latent", level=0.95)) |>
  group_by(collection_date_num) |>
  mutate_at(c("emmean", "asymp.LCL", "asymp.UCL"), softmax) |>
  dplyr::filter(collection_date_num == 18764))
#   collection_date_num variant     emmean         SE  df  asymp.LCL  asymp.UCL
# 1               18764   Other 0.86616413 0.05718955 Inf 0.86952744 0.86265454
# 2               18764    Beta 0.01788666 0.10460975 Inf 0.01636245 0.01954925
# 3               18764   Alpha 0.11594921 0.06732437 Inf 0.11411011 0.11779620


# emmeans output with mode="latent" 
as.data.frame(emmeans(fit_mblogit, ~ collection_date_num, by="variant", 
at=list(collection_date_num=dat[160,'collection_date_num']),
mode="latent", level=0.95))
#   collection_date_num variant      emmean         SE  df  asymp.LCL   asymp.UCL
# 1               18764   Other  1.96364709 0.05718955 Inf  1.8515576  2.07573655
# 2               18764    Beta -1.91637203 0.10460975 Inf -2.1214034 -1.71134070
# 3               18764   Alpha -0.04727506 0.06732437 Inf -0.1792284  0.08467828


# pairwise contrasts in marginal trends on the latent link scale
# (this would measure pairwise differences in growth rate among
# SARS-CoV2 variants) calculated using emtrends :
as.data.frame(confint(pairs(emtrends(fit_mblogit, ~ variant|collection_date_num,
                                     var = "collection_date_num",  
                                     mode = "latent",
                                     at = list(collection_date_num = 
                                                 dat[160,'collection_date_num'])),
                            reverse=TRUE), level=0.9))
#        contrast collection_date_num    estimate          SE  df   asymp.LCL
# 1  Beta - Other               18764 -0.01606104 0.005195378 Inf -0.02672347
# 2 Alpha - Other               18764 -0.03421657 0.002180285 Inf -0.03869116
# 3  Alpha - Beta               18764 -0.01815554 0.005577478 Inf -0.02960216
# asymp.UCL
# 1 -0.005398602
# 2 -0.029741991
# 3 -0.006708918

# the same emtrends on latent scale calculated using marginaleffects
# pairwise contrasts in marginal trends = 
# growth rate differences between variants :
comparisons(
  fit_mblogit,
  newdata = datagrid(collection_date_num = 
                       c(dat[160,'collection_date_num'])),
  variables = "collection_date_num",
  type = "clr",
  hypothesis = "pairwise")
#     type          term  comparison   std.error  statistic      p.value
# 1 latent  Beta - Other -0.01583897 0.005161144  -3.068888 2.148571e-03
# 2 latent Alpha - Other -0.03392950 0.002167369 -15.654696 3.085834e-55
# 3 latent  Alpha - Beta -0.01809053 0.005541114  -3.264782 1.095484e-03
# conf.low    conf.high
# 1 -0.02595463 -0.005723316
# 2 -0.03817747 -0.029681536
# 3 -0.02895091 -0.007230145

# PS wasn't sure how to do this for multiple timepoints / collection_date_num
# without putting this in a loop - is that possible?
# If I pass two timepoints I get also all contrasts between the
# different timepoints, whereas I would just like to get the
# marginal trend contrasts at the given timepoints between the different variants

# solution: put it in a loop
do.call("rbind", lapply(seq(min(dat$collection_date_num),
                            max(dat$collection_date_num),
                            by=7), function(collection_date_num)
  transform(
    comparisons(
      fit_mblogit,
      newdata = datagrid(collection_date_num = collection_date_num, 
                         model = fit_mblogit),
      hypothesis = "pairwise",
      type = "clr",
      variables = "collection_date_num"),
    collection_date_num = collection_date_num)))



# example of emtrends (on latent scale, here this describes growth rate advantage of each variant)
# here also at timepoint of observation 160
as.data.frame(confint(pairs(emtrends(fit_mblogit, ~ variant|collection_date_num,
                       var = "collection_date_num",  
                       mode = "latent",
                       at = list(collection_date_num = 
                                   dat[160,'collection_date_num'])),
              reverse=TRUE), level=0.9))
#        contrast collection_date_num    estimate          SE  df   asymp.LCL
# 1  Beta - Other               18764 -0.01606104 0.005195378 Inf -0.02672347
# 2 Alpha - Other               18764 -0.03421657 0.002180285 Inf -0.03869116
# 3  Alpha - Beta               18764 -0.01815554 0.005577478 Inf -0.02960216
#      asymp.UCL
# 1 -0.005398602
# 2 -0.029741991
# 3 -0.006708918



# 3. nnet::multinom multinomial fit ####
set.seed(1)
fit_multinom = nnet::multinom(variant ~ ns(collection_date_num, df=2), 
                          weights=count, data=dat)
fit_multinom_small = nnet::multinom(variant ~ ns(collection_date_num, df=2), 
                                   weights=count, data=datfull, maxit=10000)

dim(coef(fit_multinom_small)) # 11 3
dim(nnet:::multinomHess(fit_multinom_small)) # 33 33
dim(vcov(fit_multinom_small)) # 33 33


library(Rcpp)
library(RcppArmadillo)
# sourceCpp("..//multinomHessRcpp.cpp") # basic Rcpp version
sourceCpp("..//multinomHessRcpp2.cpp") # RcppArmadillo version
multinomHessRcpp(probs = fitted(fit_multinom_small), 
                 Z = model.matrix(fit_multinom_small),
                 coefs = t(rbind(0,coef(fit_multinom_small))),
                 row_totals = fit_multinom_small$weights)
sourceCpp("..//multinomHessRcpp4.cpp")
hessian_multinom(X=model.matrix(fit_multinom_small), 
                 Y=fitted(fit_multinom_small), 
                 beta=t(rbind(0,coef(fit_multinom_small))), 
                 W=fit_multinom_small$weights, 
                 # n=as.integer(nrow(model.matrix(fit_multinom_small))), 
                 # p=as.integer(ncol(model.matrix(fit_multinom_small))), 
                 k=as.integer(length(fit_multinom_small$lev)))
# test hessian_multinom(arma::mat X, arma::mat Y, arma::mat beta, arma::mat W, int n, int p, int k)
# X = model.matrix
# beta = coef
# Y = observed
# W = counts/weights?
# n = nrow(X)
# p = ncol(X)
# k = nr of classes


# CHECK: how to calculate Hessian / information matrix more efficiently via Kronecker products
# see
# https://stackoverflow.com/questions/73811835/faster-way-to-calculate-the-hessian-fisher-information-matrix-of-a-nnetmulti
# https://stats.stackexchange.com/questions/589848/calculate-hessian-of-multinomial-regression-variance-covariance-matrix-via-kro

# code in nnet:
object = fit_multinom_small
multinomHess <- function(object, Z = model.matrix(object))
{
  probs <- object$fitted # avoid napredict from fitted.default
  coefs <- coef(object) # 11 3
  if (is.vector(coefs)){ # ie there are only 2 response categories
    coefs <- t(as.matrix(coefs))
    probs <- cbind(1 - probs, probs)
  }
  coefdim <- dim(coefs)
  p <- coefdim[2L] # 3
  k <- coefdim[1L] # 11
  ncoefs <- k * p # 33
  kpees <- rep(p, k) # 11 3's
  n <- dim(Z)[1L] # 3780
  ##  Now compute the observed (= expected, in this case) information,
  ##  e.g. as in T Amemiya "Advanced Econometrics" (1985) pp295-6.
  ##  Here i and j are as in Amemiya, and x, xbar are vectors
  ##  specific to (i,j) and to i respectively.
  info <- matrix(0, ncoefs, ncoefs)
  Names <- dimnames(coefs)
  if (is.null(Names[[1L]])) Names <- Names[[2L]]
  else Names <- as.vector(outer(Names[[2L]], Names[[1L]],
                                function(name2, name1)
                                  paste(name1, name2, sep = ":")))
  dimnames(info) <- list(Names, Names)
  x0 <- matrix(0, p, k+1L) # 3 x 12 0 matrix 
  row.totals <- object$weights # total counts
  for (i in seq_len(n)) { # loop over rows i
    Zi <- Z[i, ] # ith row of model matrix
    xbar <- rep(Zi, times=k) * rep(probs[i, -1, drop = FALSE], times=kpees) # ith row of model matrix with p cols replicated k times * row i probs with k (11) columns, each replicated p times
    for (j in seq_len(k+1)){ # loop over all outcome levels, from 1 to k+1
      x <- x0
      x[, j] <- Zi
      x <- x[, -1, drop = FALSE] # 1st column / outcome level is dropped
      x <- x - xbar
      dim(x) <- c(1, ncoefs)
      info <- info + (row.totals[i] * probs[i, j] * crossprod(x))
    }
  }
  info
}



# save(fit_multinom_small, datfull, file="C://TEMP/smallmodel.RData")
download.file("https://www.dropbox.com/s/gt0yennn2gkg3rd/smallmodel.RData?dl=1",
              "smallmodel.RData", 
              method = "auto", mode="wb")
load("smallmodel.RData")
length(fit_multinom_small$lev) # 12
dim(coef(fit_multinom_small)) # 11 x 3 = 33
system.time(hess <- nnet:::multinomHess(fit_multinom_small)) # 0.11s
dim(hess) # 33 33

# fast way to calculate Hessian using Kroncker delta
# see https://stats.stackexchange.com/questions/525042/derivation-of-hessian-for-multinomial-logistic-regression-in-b%C3%B6hning-1992?_gl=1*103i8gk*_ga*MTI5MzMyMzMyMi4xNjU1MjgxOTAy*_ga_WCZ03SZFCQ*MTY2MzkyMjA4OC4zNS4xLjE2NjM5MjU4MzMuMC4wLjA.
download.file("https://www.dropbox.com/s/gt0yennn2gkg3rd/smallmodel.RData?dl=1",
              "smallmodel.RData", 
              method = "auto", mode="wb")
load("smallmodel.RData")
length(fit_multinom_small$lev) # 12
dim(coef(fit_multinom_small)) # 11 x 3 = 33
fit = fit_multinom_small


object = fit_multinom_small


# Utility function to calculate observed Fisher information matrix
# of multinomial fit, with 
# probs=fitted probabilities (with 1st category/column dropped)
# Z = model matrix
# row_totals = row totals
# We do this using Kronecker products, as in
# https://ieeexplore.ieee.org/abstract/document/1424458
# B. Krishnapuram; L. Carin; M.A.T. Figueiredo; A.J. Hartemink
# Sparse multinomial logistic regression: fast algorithms and generalization bounds
# IEEE Transactions on Pattern Analysis and Machine Intelligence ( Volume: 27, Issue: 6, June 2005)
calc_infmatrix = function(probs, Z, row_totals) {
  require(fastmatrix) # for kronecker.prod Kronecker product function
  
  n <- nrow(Z)
  p <- ncol(Z)
  k <- ncol(probs)
  ncoefs <- k * p
  
  info <- matrix(0, ncoefs, ncoefs)
  
  for (i in 1:n) {
    info <- info +
#    kronecker.prod((diag(probs[i,]) - tcrossprod(probs[i,]))*sqrt(row_totals[i]), 
#                     tcrossprod(Z[i,])*sqrt(row_totals[i]) )
    kronecker.prod((diag(probs[i,]) - tcrossprod(probs[i,])), 
                   tcrossprod(Z[i,])*row_totals[i] ) # this is faster & numerically near-identical to the above, but not 100% sure which is correct
  }
  return(info)  
}

# RcppArmadillo utility function to calculate observed Fisher information matrix
# of multinomial fit, with 
# probs=fitted probabilities (with 1st category/column dropped)
# Z = model matrix
# row_totals = row totals
# We do this using Kronecker products, as in
# https://ieeexplore.ieee.org/abstract/document/1424458
# B. Krishnapuram; L. Carin; M.A.T. Figueiredo; A.J. Hartemink
# Sparse multinomial logistic regression: fast algorithms and generalization bounds
# IEEE Transactions on Pattern Analysis and Machine Intelligence ( Volume: 27, Issue: 6, June 2005)

# calc_infmatrix_RcppArma(probs = object$fitted[,-1], Z = model.matrix(object), row_totals = object$weights)

# (this function was faster than similar RcppEigen implementation)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("..//calc_infmatrix_arma.cpp")

fastmultinomHess <- function(object, Z = model.matrix(object)) {
  
  probs <- object$fitted # predicted probabilities, avoid napredict from fitted.default

  coefs <- coef(object)
  if (is.vector(coefs)){ # ie there are only 2 response categories
    coefs <- t(as.matrix(coefs))
    probs <- cbind(1 - probs, probs)
  }
  coefdim <- dim(coefs)
  p <- coefdim[2L] # nr of parameters
  k <- coefdim[1L] # nr out outcome categories-1
  ncoefs <- k * p # nr of coefficients
  n <- dim(Z)[1L] # nr of observations
  
  #  Now compute the Hessian = the observed (= expected, in this case) 
  #  Fisher information matrix info
  
  # info <- calc_infmatrix(probs = probs[, -1, drop=F], # pure R function
  #                        Z = Z, 
  #                        row_totals = object$weights) 
  
  info <- calc_infmatrix_RcppArma(probs = probs[, -1, drop=F], # using faster RcppArmadillo function
                                  Z = Z, 
                                  row_totals = object$weights)

  Names <- dimnames(coefs)
  if (is.null(Names[[1L]])) Names <- Names[[2L]] else Names <- as.vector(outer(Names[[2L]], Names[[1L]],
                                function(name2, name1)
                                  paste(name1, name2, sep = ":")))
  dimnames(info) <- list(Names, Names)

  return(info)
}



download.file("https://www.dropbox.com/s/mpz08jj7fmubd68/bigmodel.RData?dl=1",
              "bigmodel.RData", 
              method = "auto", mode="wb")
load("bigmodel.RData")
object = fit_global_multi_last3m # large model
system.time(info <- fastmultinomHess(object, Z = model.matrix(object))) # 103s
system.time(info <- nnet:::multinomHess(object, Z = model.matrix(object))) # 8127s = 2.25h

system.time(V <- ginv(info))
image(V)
range(V) # -629.0583 1095.7594





system.time(info <- calc_infmatrix(probs = object$fitted[,-1],
                           Z = model.matrix(object), 
                           row_totals = object$weights)) # 600s
system.time(V <- ginv(info))
image(V)

system.time(info <- multinomHess_10rows(object,
                                        Z = model.matrix(object))) # 21.5s for 10 rows, 
dim(model.matrix(object))[[1]]/10 # 378 * 21.5s = 8127s = 2.25h


sourceCpp("..//calc_infmatrix_eigen.cpp")
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
sourceCpp("..//calc_infmatrix_arma.cpp")
system.time(info <- calc_infmatrix_RcppArma(probs = object$fitted[,-1],
                                   Z = model.matrix(object), 
                                   row_totals = object$weights)) # 103s - fastest
system.time(V <- ginv(info)) # 167s
image(V)

system.time(info2 <- calc_infmatrix_RcppEigen(probs = object$fitted[,-1],
                                            Z = model.matrix(object), 
                                            row_totals = object$weights)) # 264s
system.time(V <- ginv(info2)) # 167s
image(V)


# still parallelize this with a Paral




image(ginv(fastmultinomHess(fit_multinom_small)))

system.time(info <- fastmultinomHess(fit_global_multi_last3m)) # 617s
system.time(V <- ginv(info)) # 167s
image(V)

library(sanic)
system.time(V1 <- solve_lu(info)) # 14s
system.time(V2 <- solve_chol(info)) # 15s 
system.time(V3 <- solve_qr(info)) # 10s 
system.time(V4 <- solve_cg(info, type = "CG")) # 




solve_chol(, b)


X = model.matrix(fit) # model matrix
n = nrow(X) # nr of observations
d = ncol(X) # nr of parameters
P = fitted(fit)[,-1] # fitted probabilities for all except reference category
k = ncol(P) # = nr of categories-1 = 11
row.totals = as.vector(fit$weights) # total counts per observation
library(fastmatrix)
hess <- matrix(0, d*k, d*k)
for (j in 1:n) {
       hess <- hess +
         kronecker.prod(diag(P[j,]) - 
                          tcrossprod(P[j,]), 
                        tcrossprod(X[j,]*sqrt(row.totals[j])))
}
library(MASS)
V = ginv(hess) # variance-covariance matrix
image(V)
range(V) # -629.0583 1095.7594

# matches multinomHess result
hess_multinom <- as.matrix(nnet:::multinomHess(fit))
V_multinom = ginv(hess_multinom) # variance-covariance matrix
image(V_multinom)
range(V_multinom) # -629.0583 1095.7594


  
p <- ncol(P) # 11
k <- ncol(X) # 3
ncoefs = k*p
n <- dim(X)[1L] # 1308
lambda_phat <- matrix(0, p, p)
D <- diag(p) # delta(p) 
for (i in seq_len(p)) {
  for (j in seq_len(p)) {
    lambda_phat[,j] = D[,j]%*%P[j,]
  }
}


    

dim(delta(k-1) - crossprod(p))
dim(crossprod(p)) # 19 19



# save(fit_multinom_small, datfull, file="C://TEMP/smallmodel.RData")
download.file("https://www.dropbox.com/s/mpz08jj7fmubd68/bigmodel.RData?dl=1",
              "bigmodel.RData", 
              method = "auto", mode="wb")
load("bigmodel.RData")
length(fit_global_multi_last3m$lev) # 20
dim(coef(fit_global_multi_last3m)) # 19 x 229 = 4351

# fast way to calculate Hessian using Kroncker delta
# see https://stats.stackexchange.com/questions/525042/derivation-of-hessian-for-multinomial-logistic-regression-in-b%C3%B6hning-1992?_gl=1*103i8gk*_ga*MTI5MzMyMzMyMi4xNjU1MjgxOTAy*_ga_WCZ03SZFCQ*MTY2MzkyMjA4OC4zNS4xLjE2NjM5MjU4MzMuMC4wLjA.
Z = model.matrix(fit_global_multi_last3m)
p = fitted(fit_global_multi_last3m)[,-1]
k = length(fit_global_multi_last3m$lev)
coefs = t(as.matrix(coef(fit_global_multi_last3m)))
dim(crossprod(Z)) # 229 229
dim(crossprod(p)) # 19 x 19
nrow(Z) # 3780
library(calculus)
library(fastmatrix)
library(microbenchmark)
microbenchmark(kronecker.prod(delta(k-1)*diag(p) - crossprod(p), crossprod(Z))) # 190ms
dim(kronecker.prod(delta(k-1)*diag(p) - crossprod(p), crossprod(Z))) # 4351 x 4351
dim(delta(k-1) - crossprod(p))
dim(crossprod(p)) # 19 19

library(RcppEigen)
library(Rcpp)
sourceCpp("..//kroneckerprod.cpp")
microbenchmark(kroneckerprod(delta(k-1)*diag(p) - crossprod(p), crossprod(Z))) # 180ms

# save(fit_global_multi_last3m, data_agbyweekcountry1, file="C://TEMP/bigmodel.RData")
download.file("https://www.dropbox.com/s/mpz08jj7fmubd68/bigmodel.RData?dl=1",
              "bigmodel.RData", 
              method = "auto", mode="wb")
load("bigmodel.RData")
length(coef(fit_global_multi_last3m)) # 4351
system.time(hess <- nnet:::multinomHess(fit_global_multi_last3m)) # takes forever


# test of nnet::multinom model prediction for row 160 in dataset
# prediction on response / probs scale
nnet:::predict.multinom(fit_multinom, type="probs")[160,]
#      Other       Beta      Alpha 
# 0.86628084 0.01791583 0.11580333 

# mu = nnet:::predict.multinom(fit_multinom, type="probs")
# 
# inverse_softMax_tolink <- function(mu) {
#   log_mu <- log(mu)
#   # we normalize log(odds) so that first outcome level comes out as zero & then drop that level as in predict.mblogit with type="link"
#   # these predictions are referred to as type="latent" in the emmeans package
#   return(sweep(log_mu, 1, STATS=log_mu[,1], FUN="-")[,-1]) 
# }
# 
# inverse_softMax_tolink(mu)[1:3,]
# #        Beta     Alpha
# # 1 -22.25909 -21.13826
# # 2 -21.65267 -20.44572
# # 3 -21.04682 -19.75393

# this matches (up to rounding/fitting differences)
mclogit:::predict.mblogit(fit_mblogit, type="link")[1:3,]
#        Beta     Alpha
# 1 -22.31522 -21.13413
# 2 -21.70682 -20.44180
# 3 -21.09900 -19.75022



# predictions on link or latent scale or calculation of SEs & CLs
# not supported by nnet but they are by emmeans:

# predictions modified get_predict function for row 160: correct
# on response scale
preds_multinom_response = marginaleffects:::get_predict.multinom(fit_multinom, type="probs")
preds_multinom_response[preds_multinom_response$rowid==160,]
#     group  predicted rowid
# 160 Other 0.86628084   160
# 487  Beta 0.01791583   160
# 814 Alpha 0.11580333   160

# on latent scale
preds_multi_latent = marginaleffects:::get_predict.multinom(fit_multinom, type="clr")
# TO DO : FIX
# Error in ns(collection_date_num, knots = c(`50%` = 18792), Boundary.knots = c(18414L,  : 
# object 'collection_date_num' not found

preds_multi_latent[preds_multi_latent$rowid==160,]
#     group   predicted rowid
# 160 Other  1.96361336   160
# 487  Beta -1.91491087   160
# 814 Alpha -0.04870248   160


# emmeans result (response/prob scale, with SEs and 95% CLs on response scale)
as.data.frame(emmeans(fit_multinom, ~ collection_date_num, by="variant", 
                      at=list(collection_date_num=dat[160,'collection_date_num']),
                      mode="prob", level=0.95, df=Inf))
#    collection_date_num variant       prob          SE df   lower.CL  upper.CL
# 1               18764   Other 0.86628084 0.007354173 Inf 0.8518669 0.88069475
# 2               18764    Beta 0.01791583 0.002712467 Inf 0.0125995 0.02323217
# 3               18764   Alpha 0.11580333 0.006990414 Inf 0.1021024 0.12950429

df = as.data.frame(emmeans(fit_multinom, ~ collection_date_num, by="variant", 
                           at=list(collection_date_num=unique(dat$collection_date_num)),
                           mode="prob", level=0.95, df=Inf))
range(df$asymp.UCL) # 1.843164e-27 1.000173e+00
range(df$asymp.LCL) # -0.0006037624  1.0000000000

# emmeans result (latent scale, with SEs and 95% CLs on latent scale)
softmax <- function(x) exp(x) / sum(exp(x)) 
as.data.frame(emmeans(fit_multinom, ~ collection_date_num, by="variant", 
                      at=list(collection_date_num=dat[160,'collection_date_num']),
                      mode="latent", level=0.95, df=Inf)) 
# NOTE: nnet::multinom takes into account df, but mclogit::mblogit does not & uses df=Inf
# change line 111 in https://github.com/rvlenth/emmeans/blob/master/R/multinom-support.R to object$edf = object$model.df
#   collection_date_num variant      emmean         SE  df  asymp.LCL   asymp.UCL
# 1               18764   Other  1.96361336 0.05710552 Inf  1.8516886  2.07553813
# 2               18764    Beta -1.91491087 0.10442145 Inf -2.1195731 -1.71024860
# 3               18764   Alpha -0.04870248 0.06727033 Inf -0.1805499  0.08314493


# predictions marginaleffects using type="probs"
preds_response = predictions(fit_multinom, 
            newdata = datagrid(collection_date_num=unique(dat$collection_date_num)),
            type = "probs") |> 
  transform(conf.low = predicted - 1.96 * std.error,
            conf.high = predicted + 1.96 * std.error)
range(preds_response$conf.low)  # -0.0006041339  1.0000000000
range(preds_response$conf.high) # 1.843259e-27   1.000174e+00

# predictions marginaleffects using type="latent" & backtransformed to response scale
softmax <- function(x) exp(x) / sum(exp(x)) 
preds_backtransformed_latent = predictions(fit_multinom, 
            newdata = datagrid(collection_date_num=unique(dat$collection_date_num)),
            type = "clr") |> 
  transform(conf.low = predicted - 1.96 * std.error,
            conf.high = predicted + 1.96 * std.error) |>
  group_by(rowid) |>
  mutate_at(c("predicted", "conf.low", "conf.high"), softmax)

range(preds_backtransformed_latent$conf.low) # 1.02708e-28 1.00000e+00
range(preds_backtransformed_latent$conf.high) # 6.316429e-28 1.000000e+00








#, 
#            transform_post = function (eta) softMax(eta) ) # NOTE: insight::link_inverse(fit_multinom) NOT CORRECT
# with type="latent" or "link": Error: The `type` argument for models of class `mblogit` must be an element of: link, response
# presumably R/type_dictionary.R still needs to be modified

# leaving out transform_post & using type="probs":
# SE matches emmeans result: 
# as.data.frame(emmeans(fit_multinom, ~ collection_date_num, by="variant", 
# at=list(collection_date_num=dat[160,'collection_date_num']),
# mode="prob", level=0.95))
# but in contrast to emmeans results leaves out confidence intervals
predictions(fit_multinom, 
            newdata = datagrid(collection_date_num=dat[160,'collection_date_num']),
            type="probs")
#   rowid  type group  predicted   std.error    count collection_date_num
# 1     1 probs Other 0.86628084 0.007356969 73.10398               18764
# 2     1 probs  Beta 0.01791583 0.002713518 73.10398               18764
# 3     1 probs Alpha 0.11580333 0.006993107 73.10398               18764


# example of emtrends (on latent scale, here this describes growth rate advantage of each variant)
# here also at timepoint of observation 160
as.data.frame(confint(pairs(emtrends(fit_multinom, ~ variant|collection_date_num,
                                     var = "collection_date_num",  
                                     mode = "latent",
                                     at = list(collection_date_num = 
                                                 dat[160,'collection_date_num'])),
                            reverse=TRUE), level=0.9))
# contrast collection_date_num    estimate          SE df    lower.CL
# 1  Beta - Other               18764 -0.01596339 0.005179619  6 -0.02899609
# 2 Alpha - Other               18764 -0.03425342 0.002181425  6 -0.03974221
# 3  Alpha - Beta               18764 -0.01829002 0.005563444  6 -0.03228847
# upper.CL
# 1 -0.002930702
# 2 -0.028764628
# 3 -0.004291572

microbenchmark(as.data.frame(confint(pairs(emtrends(fit_multinom, ~ variant|collection_date_num,
                                                    var = "collection_date_num",  
                                                    mode = "latent",
                                                    at = list(collection_date_num = 
                                                                dat[160,'collection_date_num'])),
                                           reverse=TRUE), level=0.9))) # 169 ms

# the same emtrends on latent scale calculated using marginaleffects
# pairwise contrasts in marginal trends = 
# growth rate differences between variants :
comparisons(
  fit_multinom,
  newdata = datagrid(collection_date_num = 
                       c(dat[160,'collection_date_num'])),
  variables = "collection_date_num",
  type = "clr", # type="link"=alr gives the same result, type="logit" would give incorrect result
  hypothesis = "pairwise")
# type          term  comparison   std.error  statistic      p.value    conf.low    conf.high
# 1  clr  Beta - Other -0.01574215 0.005145516  -3.059392 2.217869e-03 -0.02582717 -0.005657123
# 2  clr Alpha - Other -0.03396620 0.002168503 -15.663429 2.689945e-55 -0.03821639 -0.029716011
# 3  clr  Alpha - Beta -0.01822405 0.005527198  -3.297159 9.766826e-04 -0.02905716 -0.007390941


microbenchmark(comparisons(
  fit_multinom,
  newdata = datagrid(collection_date_num = 
                       c(dat[160,'collection_date_num'])),
  variables = "collection_date_num",
  type = "clr",
  hypothesis = "pairwise")) # marginaleffects 88 ms


# all comparisons against the reference level (Other) can be obtained using type="link"="alr"
microbenchmark(marginaleffects(fit_multinom, 
                               type = "link", 
                               variables = c("collection_date_num"),
                               by = c("group", "collection_date_num"), 
                               vcov = vcov(fit_multinom),
                               newdata = data.frame(collection_date_num=dat[160,'collection_date_num']))) # 105 ms
# correct for Beta and Alpha vs Other
# TO DO: adapt type="link" to allow additional ref argument with ref = nr of reference group?

microbenchmark(marginaleffects(fit_multinom, 
                               type = "link", 
                               variables = c("collection_date_num"),
                               by = c("group", "collection_date_num"), 
                               # vcov = FALSE, # vcov(fit_multinom),
                               newdata = datagrid(collection_date_num=dat[160,'collection_date_num'],
                                                  # group = c("Other","Beta"),
                                                  grid_type = "counterfactual"),
                               hypothesis = "reference")) # 114 ms




# FINAL NOTE: ILLUSTRATION THAT vcov MATRIX COMES OUT IN DIFFERENT ORDER
# IN nnet::multinom AND mclogit::mblogit
colnames(vcov(fit_multinom)) # (in row major order, this is standard)
# [1] "Beta:(Intercept)"                       "Beta:ns(collection_date_num, df = 2)1" 
# [3] "Beta:ns(collection_date_num, df = 2)2"  "Alpha:(Intercept)"                     
# [5] "Alpha:ns(collection_date_num, df = 2)1" "Alpha:ns(collection_date_num, df = 2)2"

colnames(vcov(fit_mblogit)) # (in column major order, this is nonstandard)
# [1] "Beta~(Intercept)"                       "Alpha~(Intercept)"                     
# [3] "Beta~ns(collection_date_num, df = 2)1"  "Alpha~ns(collection_date_num, df = 2)1"
# [5] "Beta~ns(collection_date_num, df = 2)2"  "Alpha~ns(collection_date_num, df = 2)2"

# emmeans deals with this here:
# https://github.com/rvlenth/emmeans/blob/master/R/multinom-support.R
# lines 112-118:
# # we have to arrange the vcov elements in row-major order
#     if(missing(vcov.))
# vcov. = vcov(object)
# perm = matrix(seq_along(as.numeric(object$coefmat)), 
#               ncol = ncol(object$coefmat))
# perm = as.numeric(t(perm))
# vcov. = vcov.[perm, perm]


# 4. MGLM::MGLMreg multinomial fit (also support Dirichlet multinomial) ####
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
fit_mblogit = mblogit(formula=variant ~ ns(collection_date_num, df=2),
                      weights=count,
                      data=dat,
                      from.table=FALSE, dispersion=FALSE,
                      control=mclogit.control(maxit=100))
# returns error Error: no valid set of coefficients has been found: please supply starting values


# MGLM fits: return errors ####
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

