# see predict.mblogit
# https://github.com/melff/mclogit/blob/cf0a90fa566f7a5c7ef8bbe729dde6e8ad86a031/pkg/R/mblogit.R

predict2.mblogit <- function(object, 
                             newdata=NULL,
                             type=c("link","response"),
                             se.fit=FALSE,...){
  
  type <- match.arg(type)
  mt <- terms(object)
  rhs <- delete.response(mt)
  if(missing(newdata)){
    m <- object$model
    na.act <- object$na.action
  }
  else{
    m <- model.frame(rhs,data=newdata,na.action=na.exclude)
    na.act <- attr(m,"na.action")
  }
  X <- model.matrix(rhs,m,
                    contrasts.arg=object$contrasts,
                    xlev=object$xlevels
  )
  rn <- rownames(X)
  D <- object$D
  XD <- X%x%D
  rspmat <- function(x){
    y <- t(matrix(x,nrow=nrow(D)))
    colnames(y) <- rownames(D)
    y
  }
  
  eta <- c(XD %*% coef(object))
  eta <- rspmat(eta)
  rownames(eta) <- rn
  if(se.fit){
    V <- vcov(object)
    stopifnot(ncol(XD)==ncol(V))
  }
  
  if(type=="response") {
    exp.eta <- exp(eta)
    sum.exp.eta <- rowSums(exp.eta)
    p <- exp.eta/sum.exp.eta
    
    if(se.fit){
      p.long <- as.vector(t(p))
      s <- rep(1:nrow(X),each=nrow(D))
      
      wX <- p.long*(XD - rowsum(p.long*XD,s)[s,,drop=FALSE])
      se.p.long <- sqrt(rowSums(wX * (wX %*% V)))
      se.p <- rspmat(se.p.long)
      rownames(se.p) <- rownames(p)
      if(is.null(na.act))
        list(fit=p,se.fit=se.p)
      else
        list(fit=napredict(na.act,p),
             se.fit=napredict(na.act,se.p))
    }
    else {
      if(is.null(na.act)) p
      else napredict(na.act,p)
    }
  }
  else if(se.fit) {
    se.eta <- sqrt(rowSums(XD * (XD %*% V)))
    se.eta <- rspmat(se.eta)
    eta <- eta[,-1,drop=FALSE]
    se.eta <- se.eta[,-1,drop=FALSE]
    if(is.null(na.act))
      list(fit=eta,se.fit=se.eta) 
    else
      list(fit=napredict(na.act,eta),
           se.fit=napredict(na.act,se.eta))
  }
  else {
    eta <- eta[,-1,drop=FALSE]
    if(is.null(na.act)) eta
    else napredict(na.act,eta)
  }
}
