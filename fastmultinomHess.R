# FASTER ALTERNATIVE TO nnet::multinomHess function TO CALCULATE HESSIAN OF MULTINOMIAL FIT ####
# see https://stackoverflow.com/questions/73811835/faster-way-to-calculate-the-hessian-fisher-information-matrix-of-a-nnetmulti/73840453#73840453


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
sourceCpp(".//src//calc_infmatrix_arma.cpp")

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
  # note: this could still be parallelized either within Rcpp code with parallelReduce or
  # on R side
  
  Names <- dimnames(coefs)
  if (is.null(Names[[1L]])) Names <- Names[[2L]] else Names <- as.vector(outer(Names[[2L]], Names[[1L]],
                                                                               function(name2, name1)
                                                                                 paste(name1, name2, sep = ":")))
  dimnames(info) <- list(Names, Names)
  
  return(info)
}
