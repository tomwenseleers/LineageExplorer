#' Function to refresh the data 
#' It checks the timestamp, compares to today, and if the timestamp
#' is earlier, it removes the old file and reads the new one.
#' 
#' @param name the name of the dataset (NOT the file)
#' @param link a link from where the data is downloaded
#' @param suffix the suffix that should be added (normally .csv)
#' @param dir the directory where the file should be looked up
#' @param force if TRUE, download regardless of timestamp
#' @param process a function used to process the data before saving. It should take a data frame as first argument and return the data frame.
#' 
#' @return the (refreshed) data
#' 
source("functions/stamp.R")

refresh <- function(name, file, 
                    loc = "https://epistat.sciensano.be/Data/", 
                    suffix = ".csv", dir = "Data",
                    force = FALSE,
        process = function(x) x ) {
  # stop()
  
  oldfile <- grep(name,dir(dir), value = TRUE)
  link <- file.path(loc,file)
  # Check whether the old file is there.
  if(length(oldfile)){
    thestamp <- readstamp(oldfile)
    go <- thestamp < Sys.Date()
    
    if(go){
      unlink(file.path(dir,oldfile))
    }
  } else {
    go <- TRUE
  }
  
  if(go || force){
    thefile <- stamp(name)
    tmp <- tryCatch(
      read.csv(link, encoding = "UTF-8"),
      error = function(e){}, warning = function(w){
        stop(paste("The data from",link,"could not be downloaded. The server of epistat might be temporarily down. Check whether you can access\nhttps://epistat.wiv-isp.be/Covid/"),
           call. = FALSE)
      }
    )
    
    write.csv(tmp, file = file.path(dir,thefile), 
              row.names = FALSE)
  } else {
    tmp <- read.csv(file.path(dir,oldfile))
  }
  tmp <- process(tmp)
  return(tmp)
}
