#-----------------------------------------------
# Add the totals of the groups mentioned
# x : the data
# values : the variable(s) with the values
# groups : the variable(s) with the groups
# along : the variable(s) along which the totals should be calculated
source("functions/helperfunctions.R")

add_totals <- function(x,values, along, groups, name = "All",
                       na.rm = FALSE, use.na = TRUE){
  # Check the groups
  dims <- length(groups)
  # Check the names
  if(length(name) == 1) {
    name <- rep(name,dims)
  }
  # make the tapply list
  applyby <- x[c(groups, along)]
  mforsum <- seq.int(dims)
  
  if(use.na){
    applyby <- lapply(applyby, add_unknown)
  }
  
  # How many values
  morevalues <- length(values) > 1
  
  # add margins
  msums <- function(v){
    tab <- tapply(x[[v]], applyby, sum, na.rm = na.rm)
    tab <- replaceby0(tab)
    addmargins(tab, margin = mforsum)
  }
  
  tmp <- msums(values[1])

  # reprocess the names
  dnames <- dimnames(tmp)
  
  for(i in seq.int(dims)){
    nnames <- length(dnames[[i]])
    dnames[[i]][nnames] <- name[i]
  }
  dimnames(tmp) <- dnames
  out <- as.data.frame(as.table(tmp),
                       responseName = values[1],
                       stringsAsFactors = FALSE)
  
  if(morevalues){
    
    for(i in values[-1]){
      extra <- msums(i)
      out[[i]] <- as.vector(extra)
    }
  }
  return(out)
}
