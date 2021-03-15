#' Replace a RDS file with a new one
#' 
#' This function checks the time stamp on a file, and if the
#' time stamp is not today, it replaces the file with a new
#' file and a time stamp. If the file doesn't exist, then
#' it just creates the file
#' 
#' @param x the data object to be stored
#' @param prefix the prefix of the file
#' @param extension the extension of the file
#' @param path the path to the file
#' 
#' @return the file name (silently)
#' 
replacefile <- function(x, prefix, extension, path){
  
  fpattern <- paste0(prefix,"\\d{4}-\\d{2}-\\d{2}\\.",extension)
  oldfile <- dir(path,pattern = fpattern,
                 ignore.case = TRUE)
  
  hasfile <- length(oldfile) > 0
  
  if(hasfile){
    tstamp <- readstamp(oldfile)
    id <- which(as.Date(tstamp) < Sys.Date())
    unlink(file.path(path,oldfile[id]))
    
    if(length(oldfile[-id])){
      newfile <- oldfile[-id]
      needsave <- FALSE
    } else {
      needsave <- TRUE
    }
  } else {
    needsave <- TRUE
  }
  
  if(needsave){
    if(substring(extension,1,1) != "."){
      extension <- paste0(".", extension)
    }
    newfile <- stamp(prefix, extension)
    extension <- tolower(extension)
    
    if(extension == ".csv"){
      write.csv(x,
                file = file.path(path, newfile),
                row.names = FALSE)
    } else if(extension == ".rds"){
      saveRDS(x,
              file = file.path(path, newfile))
    }
  } # END if needsave
  return(invisible(newfile))
}
