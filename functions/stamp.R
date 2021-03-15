#' Functions to create and check the timestamps.
#' 
stamp <- function(x, suffix = ".csv"){
  if(substring(suffix,1,1) != "."){
    suffix <- paste0(".",suffix)
  }
  # make sure the same format is used everywhere.
  thedate <- format(Sys.Date(), "%Y-%m-%d")
  paste0(x, thedate, suffix)
}

readstamp <- function(x, as_date = TRUE){
  tmp <- gsub(".*(\\d{4}-\\d{2}-\\d{2})","\\1",x)
  if(as_date){
    as.Date(tmp)
  } else {
    tmp
  }
}
