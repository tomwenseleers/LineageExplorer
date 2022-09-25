library(MASS) # For 'housing' data
library(nnet)
library(memisc)
library(mclogit)
library(insight)
library(dplyr)
library(tidyr)

fit_multinom = multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
fit_mblogit = mblogit(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

# MBLOGIT MULTINOMIAL MODEL PREDICTIONS
# ON RESPONSE SCALE
predict(fit_mblogit, type="response")[1,]
# 0.3955687 0.2601077 0.3443236

preds_mblogit = as.data.frame(insight::get_predicted(fit_mblogit, ci=T, predict="response"))
  preds_mblogit[preds_mblogit$Row==1,"Predicted"] 
# 0.3955687 0.2601077 0.3443236
# ci=T is ignored

# ON LINK SCALE
preds_link = c(0,  # PS reference category is dropped here so adding 0
  predict(fit_mblogit, type="link")[1,])
#        Low     Medium       High 
#  0.0000000 -0.4192287 -0.1387428 
names(preds_link)[[1]] = levels(insight::get_response(fit_mblogit))[[1]]
preds_link
#       Low     Medium       High 
# 0.0000000 -0.4192287 -0.1387428

predict(fit_multinom, type="probs")
preds_multinom =as.data.frame(insight::get_predicted(fit_multinom, ci=T))
sum(preds_multinom[preds_multinom$Row==1,"Predicted"])




summary(fit_multinom)
summary(fit_mblogit)

# coefficients come out in the same order in both
dim(coef(fit_multinom))
dim(fit_mblogit$coefmat) 
# note that coef(fit_mblogit)returns a vector rather than a matrix, so that's no good

# see https://github.com/easystats/insight/blob/main/R/link_function.R
# and https://github.com/easystats/insight/blob/main/R/link_inverse.R
insight::link_inverse(fit_mblogit) # function (eta) eta, not correct
insight::link_function(fit_mblogit) # function (mu) mu, not correct
insight::link_inverse(fit_multinom) # link_inverse.multinom <- function(x, ...) { stats::make.link("logit")$linkinv }, not correct - this should be the SoftMax function function(x){ expx <- exp(x); return(expx/sum(expx)) }
insight::link_function(fit_multinom) # link_function.multinom <- function(x, ...) { stats::make.link(link = "logit")$linkfun }

# see also
# https://stackoverflow.com/questions/69586966/what-is-the-scale-of-parameter-estimates-produced-by-nnetmultinom
# https://online.stat.psu.edu/stat504/book/export/html/788
# https://rpubs.com/beane/n4_2
# https://stackoverflow.com/questions/17283595/fitted-values-for-multinom-in-r-coefficients-for-reference-category?rq=1

predict(fit_mblogit, type="link") 
# this leaves out the reference category Low - that one would be zero
predict(fit_mblogit, type="response") # has all categories, matches multinom
library(emmeans)
emmeans(fit_mblogit, ~ Sat|Infl + Type + Cont, mode="prob")
emmeans(fit_mblogit, ~ Sat|Infl + Type + Cont, mode="latent")
vcov(fit_mblogit)
rownames(vcov(fit_mblogit))
rownames(emmeans:::.my.vcov(fit_mblogit))

predict(fit_multinom, type="probs") # has all categories, matches mblogit
stats::make.link("logit")$linkfun(predict(fit_multinom, type="probs")) # not correct
emmeans(fit_multinom, ~ Sat|Infl + Type + Cont, mode="prob", df=Inf)
emmeans(fit_multinom, ~ Sat|Infl + Type + Cont, mode="latent", df=Inf)
vcov(fit_multinom) 
# note: comes out in the wrong order, 
# see emm_basis.mblogit https://rdrr.io/cran/emmeans/src/R/multinom-support.R
rownames(vcov(fit_multinom)) # incorrect order
rownames(emmeans:::.my.vcov(fit_multinom))

# example where we would like predictions for our original data
newdata = housing 
sapply(1:nrow(newdata), function (row) { fit_multinom_refgrid = ref_grid(fit_multinom, at = newdata[row,])
                                         predict(fit_multinom_refgrid, type="probs") } )

# for predictions with mode="latent"
object = fit_multinom
bhat = t(coef(object))
V = emmeans:::.my.vcov(object) # variance-covariance matrix
# NOTE: entries in vcov(object) come out in same order as
# in as.numeric(bhat), even though latter has been transposed
k = ifelse(is.matrix(coef(object)), ncol(bhat), 1)
# m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
X = model.matrix(object) # trms, m, contrasts.arg = object$contrasts
dim(X) # 72 7
# recenter for latent predictions
pat = (rbind(0, diag(k + 1, k)) - 1) / (k + 1)
X = kronecker(pat, X)
# colSums(X) # now sums to zero
dim(X) # 216 14
linear.predictor = as.vector(rbind(0, coef(fit_multinom)) %*% X[1,])


# see https://stackoverflow.com/questions/17283595/fitted-values-for-multinom-in-r-coefficients-for-reference-category?rq=1
# this is the inverse link function for multinomial (here assuming eta is a matrix with different outcome levels as columns)
# = inverse generalized logit
softMax <- function(eta){
  exp_eta <- exp(eta)
  return(sweep(exp_eta, 1, STATS=rowSums(exp_eta), FUN="/"))
}

# see https://math.stackexchange.com/questions/2786600/invert-the-softmax-function
# this is the link function for multinomial (here assuming mu is a matrix with predicted probabilities, with different outcome levels as columns)
# = generalized logit
inverse_softMax <- function(mu) {
  log_mu <- log(mu)
  return(sweep(log_mu, 1, STATS=rowMeans(log_mu), FUN="-")) # we let the log(odds) sum to zero - these predictions are referred to as type="latent" in the emmeans package
}


# predict on response, link or latent scale for nnet::multinom models
object = fit_multinom
betahat = t(rbind(0, coef(object))) # 0 for reference category
colnames(betahat)[[1]] = levels(insight::get_response(object))[[1]]
X = model.matrix(delete.response(terms(object)), newdat) # or whatever newdat argument
dim(X) # 72 7
preds_link = X %*% betahat # eta = predictions on linear predictor scale
preds_link
dim(preds_link) # 72 3
preds_latent = sweep(preds_link, 1,  # predictions on latent centered linear predictor scale
                                  STATS=rowMeans(preds_link), 
                                  FUN="-")
dim(preds_latent) # 72 3
preds_response = softMax(preds_latent) # predictions on response / probability scale
preds_response

# the standard errors on the latent/link scale would be given by
V <- vcov(object)
Xaug <- cbind(matrix(0, ncol=ncol(X), nrow=nrow(X)), X)
var.fit <- rowSums((Xaug %*% V) * Xaug)  # point-wise variance for predicted mean, as we only need diag(Xp %*% V %*% t(Xp))
std.er <- sqrt(var.fit)

D <- fit_mblogit$D
XD <- X %x% D
rspmat <- function(x){
  y <- t(matrix(x,nrow=nrow(D)))
  colnames(y) <- rownames(D)
  y
}
se.eta <- rspmat(sqrt(rowSums(XD * (XD %*% V))))
# see https://github.com/melff/mclogit/blob/master/pkg/R/mblogit.R

predict(fit_mblogit, type="link", se.fit=TRUE)


stats::make.link("logit")$linkinv(preds_link) # this would not be correct
# correct would be softMax(preds_link)
stats::make.link("logit")$linkfun(preds_response) # this would not be correct
# this is correct
# for conversion to latent scale
inverse_softMax(preds_response)[1,]
softMax(inverse_softMax(preds_response))[1,]

# softMax(matrix(log(preds_response[1,])-log(preds_response[1,1]),nrow=1))

# 
