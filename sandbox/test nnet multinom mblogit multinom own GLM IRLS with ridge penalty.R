download.file("https://www.dropbox.com/s/o6iu51wu7x90omd/data.rds?dl=1",
              "../data.rds",
              method = "auto", mode="wb")
datsubs = readRDS(file = "../data.rds")

# 1. nnet::multinom fit ####
library(nnet)
library(splines)
fit_nnet <- nnet::multinom(variant ~ date_num + date_num:region + division,
                           weights = count,
                           data = datsubs,
                           maxit = 10000, MaxNWts = 100000)
hist(coef(fit_nnet))
range(coef(fit_nnet))

# 2. mblogit multinomial fit ####
# devtools::install_github("melff/mclogit",subdir="pkg")
library(mclogit)
fit_mblogit <- mblogit(formula=variant ~ date_num + date_num:region + division, # runs OK too
                                               weights=count,
                                               data=datsubs,
                                               from.table=FALSE,
                                               dispersion=FALSE, # fit overdispersion?
                                               control=mclogit.control(maxit=10000))
hist(coef(fit_mblogit))
range(coef(fit_mblogit))

plot(coef(fit_nnet),
     coef(fit_mblogit), pch=16)


# 3. regular minimal IRLS GLM algorithm, adapted from https://bwlewis.github.io/GLM/

glm_irls = function(X, y, weights=rep(1,nrow(X)), family=poisson(log), maxit=25, tol=1e-16) {
  if (!is(family, "family")) family = family()
  variance = family$variance
  linkinv = family$linkinv
  mu.eta = family$mu.eta
  etastart = NULL
  
  nobs = nrow(X)    # needed by the initialize expression below
  nvars = ncol(X)   # needed by the initialize expression below
  eval(family$initialize) # initializes n and fitted values on response scale mustart (set equal to y)
  
  eta = family$linkfun(mustart) # initialize eta = fitted values on link scale with this mustart
  dev.resids = family$dev.resids
  dev = sum(dev.resids(y, linkinv(eta), weights))
  devold = 0
  beta_old = rep(1, nvars)
  
  for(j in 1:maxit)
  {
    # E-step
    mu = linkinv(eta) # mu = fitted values on response scale
    varg = variance(mu) #  variance of the response variable given the mean mu
    gprime = mu.eta(eta) # derivative of the inverse link function with respect to the mean, also known as the "working response" or "score function"
    z = eta + (y - mu) / gprime # "working residual" or "adjusted z scale", which is calculated by adding the residuals, scaled by the inverse of the working response, to the linear predictor; potentially -offset if you would have an offset argument as well
    W = weights * as.vector(gprime^2 / varg) # weights in weighted least square regression below

    # M-step: update coefficients beta using a weighted least square regression (regressing z on X using weights W)
    beta = solve(crossprod(X, W*X), crossprod(X, W*z), tol=2*.Machine$double.eps)
  
    # update predictions for E-step
    eta = X %*% beta # predictions on link scale, potentially +offset if you would have an offset argument as well
    dev = sum(dev.resids(y, mu, weights)) # deviance
    if (abs(dev - devold) / (0.1 + abs(dev)) < tol) break
    devold = dev
    beta_old = beta
  }
  list(coefficients=t(beta), iterations=j)
}

# hessian = expected Fisher information matrix = XTWX = tcrossprod(X, W) %*% X = crossprod(t(X), W) %*% X
# vcov matrix = solve(hessian)
# SEs = sqrt(diag(vcov))


## Poisson GLM from Dobson (1990) Page 93: Randomized Controlled Trial
y <- counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
X <- model.matrix(counts ~ outcome + treatment)

coef(glm.fit(x=X, y=y, family=poisson(log))) 
# (Intercept)      outcome2      outcome3    treatment2    treatment3 
# 3.044522e+00 -4.542553e-01 -2.929871e-01 -7.635479e-16 -9.532452e-16

coef(glm_irls(X=X, y=y, family=poisson(log)))
#         (Intercept)   outcome2   outcome3    treatment2   treatment3
# [1,]    3.044522 -0.4542553 -0.2929871 -3.151689e-16 -8.24099e-16

sum(coef(glm.fit(x=X, y=y, family = poisson(log))) - coef(glm_irls(X=X, y=y, family=poisson(log)))) # 5E-15

## Binomial GLM
X = as.matrix(data.frame(Intercept=1, x=seq(0:100)))
beta = c(-1, 0.05)
y1 = sapply(plogis(X %*% beta), function (prob) rbinom(1, 10, prob))
y2 = 10-y1
weights = y1+y2
y = y1/weights
coef(glm.fit(x=X, y=y, weights=weights, family=binomial(logit))) 

coef(glm_irls(X=X, y=y, weights=weights, family=binomial(logit)))


# my own multinomial function with ridge penalty ####
# still have to finish this one


softMax_derivative <- function(eta) {
  mu <- softMax(eta)
  mu_matrix <- mu # -p(i)*p(j) if i <> j https://www.mldawn.com/the-derivative-of-softmaxz-function-w-r-t-z/, https://deepnotes.io/softmax-crossentropy#derivative-of-softmax, https://stats.stackexchange.com/questions/453539/softmax-derivative-implementation
  diag(mu_matrix) <- diag(mu * (1 - mu)) #  p(i)*(1-p(i)) if i=j
  return(mu_matrix)
}

d_softmax <- function(eta) {
  mu <- softMax(eta)
  delta_matrix <- diag(rep(1, nrow(eta)))
  mu %*% (delta_matrix - t(mu) %*% mu)
}


softMax_derivative <- function(eta) { # equivalent of mu.eta(eta), i.e. derivative of the inverse link function with respect to eta, e.g. exp(eta)/(1+exp(eta))^2 for the binomial case
  mu <- exp(eta) / rowSums(exp(eta)) # inverse link function
  
  expx <- exp(as.vector(eta))
  s <- sum(expx)^2
  l <- length(eta)
  ans <- (-1) * outer(expx,expx) / s
  diag(ans) <- sapply(1:l,function(i){return(expx[i]*sum(expx[-i]))}) / s
  return(ans)
  
  
  mu # https://deepnotes.io/softmax-crossentropy#derivative-of-softmax
  mu_matrix <- mu %*% t(mu)
  diag(mu_matrix) <- mu * (1 - mu)
  mu_matrix <- mu_matrix - t(mu_matrix)
  return(mu_matrix)
}

inverse_softMax_to_clr = function(mu) { log_mu = log(mu)
                                        eta = log_mu-rowMeans(log_mu) # convert predictions on response scale to centered log ratio link scale
                                        return(eta) }

inverse_softMax_to_alr = function(mu) { eta = log(mu)-log(mu[,1]) # link function for multinomial: converts predictions on response scale to additive log ratio scale
                                        return(eta) }

softMax = function(eta) { require(mclustAddons)   # inverse link function for multinomial, note mclustAddons also has softmax(eta)
                s = as.vector(logsumexp(eta))
                mu = exp(eta-s) # = exp(eta) / rowSums(exp(eta)), but numerically more stable
                return(mu)
                }


multinom_irls = function(X, # covariate matrix
                         Y, # outcome matrix with counts
                         lambda=0, # ridge penalty
                         w=rep(1, ncol(X)), # adaptive penalty weights
                         maxit=25, # max nr of IRLS iterations
                         tol=1e-08) { # convergence tolerance
  k = ncol(Y) # number of classes
  nobs = nrow(X) # number of observations
  p = ncol(X) # number of predictors
  
  # initialise mu & eta
  mustart = (Y + 0.5)/(rowSums(Y) + 1) # analogous to how binomial GLMs are initialised, cf. binomial(logit)$initialize
  mustart = mustart/rowSums(mustart) # make sure rows sums to 1
  eta = inverse_softMax_to_alr(mustart) # initialize eta = fitted values on additive logratio link scale with this mustart, implement as family$linkfun(mustart)?
  
  beta_old = matrix(1, nrow=p, ncol=k-1) # coefficients, with k-1 columns as we will normalize the first column of coefficients for the reference baseline category to zero
  
  for(j in 1:maxit) {
    # E-step
    mu = softMax(eta) # mu = fitted values on response scale, implement as linkinv(eta)?
    varg = mu*(1-mu)  #  variance of the response variable given the mean mu, implement as variance(mu)
    gprime = mu.eta(eta) # derivative of the inverse link function with respect to eta, also known as the "working response" or "score function", implement as mu.eta(eta)? for binomial derivative of inverse logit with respect to eta = FullSimplify[D[Exp[eta]/(1 + Exp[eta]), eta]] = Exp[eta]/(1+Exp[eta])^2
    z = eta + (y - mu) / gprime # "working residual" or "adjusted z scale", which is calculated by adding the residuals, scaled by the inverse of the working response, to the linear predictor; potentially -offset if you would have an offset argument as well
    W = weights * as.vector(gprime^2 / varg) # weights in weighted least square regression below
    
    # M-step: update coefficients beta using a weighted least square regression (regressing z on X using weights W)
    beta = solve(crossprod(X, W*X), crossprod(X, W*z), tol=2*.Machine$double.eps)
    
    # update predictions for E-step
    eta = X %*% beta # predictions on link scale, potentially +offset if you would have an offset argument as well
    
    
    
    eta = X %*% beta # predictions on generalized logit scale
    P = exp(eta) / rowSums(exp(eta)) # predicted probabilities (softmax transforms)
    z = eta + (Y - P) %*% t(P) / n
    W = diag(P %*% t(P))
    beta_old = beta
    penalty = lambda * diag(w)
    x = solve(crossprod(X, W*X) + penalty, crossprod(X, W*z), tol=2.Machine$double.eps)
    if(sqrt(crossprod(x-xold)) < tol) break
  }
  list(coefficients=x, iterations=j)
}
 
# generate some multinomial data
library(MGLM)
## Generate data
n <- 2000
p <- 5 # nr of variables in X (not counting intercept)
d <- 4 # nr of outcome categories
m <- rep(20, n)
set.seed(1234)
X <- 0.1* matrix(rnorm(n*p),n, p)
alpha <- matrix(1, p, d-1)
beta <- matrix(1, p, d-1)
Alpha <- exp(X %*% alpha)
Beta <- exp(X %*% beta)
Y <- rgdirmn(n, m, Alpha, Beta)

# solution returned by MGLM
gdm.reg <- MGLMreg(Y~X, dist="GDM", LRT=FALSE)
coef(gdm.reg)
dim(coef(gdm.reg)) # 6 x 6

# # solution returned by Rfast2::multinom.reg
# library(Rfast2)
# Rfast2::multinom.reg(ytrain, xtrain) # only accepts data in long format, so no good


# solution returned by MultinomialMutations - very fast - but check quality of solution
remotes::install_github("MultinomialMutations/MultinomialMutations")
library(MultinomialMutations)
colnames(Y) = paste0("OUT",1:ncol(Y))
colnames(X) = paste0("V",1:ncol(X))
df = data.frame(Y, X)
fit_multmut = fast_multinom(cbind(OUT1, OUT2, OUT3, OUT4) ~ V1+V2+V3+V4+V5, data=df, refLevel=1, loglik=T, predictions=T, VC=T)
coefs_multmut = do.call(cbind, coef(fit_multmut))
dim(coefs_multmut) # 6 x 3
head(predict.fast_multinom(fit_multmut))


# solution returned by glmnet with lambda close to zero
library(glmnet)
fit_glmnet = glmnet(x=X, y=Y, intercept=T, alpha=1, nlambda=1000, lambda.min= 1E-7, family="multinomial", 
                    type.multinomial="grouped") # ungrouped
plot(fit_glmnet, xvar="lambda")
coefs = coef(fit_glmnet, s=1E-10)
do.call(cbind, coefs)

# sparse multinomial regression worth checking
# https://ieeexplore.ieee.org/abstract/document/1424458
# https://github.com/inzxx/adastra/blob/master/src/Accord24/Accord.Statistics/Models/Regression/Nonlinear/Fitting/LowerBoundNewtonRaphson.cs
# https://github.com/slavomir-sidor/princeton-mvpa-toolbox/blob/master/core/learn/smlr.m
# 



