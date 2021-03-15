#' This file contains a number of convenience functions for processing
#' 

#----------------------------------
# Link provinces to regions
regions <- c(
  "Brussels"= "Brussels",
  "Antwerpen"= "Flanders",
  "BrabantWallon" = "Wallonia",
  "Hainaut" = "Wallonia",
  "Liège" = "Wallonia",
  "Limburg" = "Flanders",
  "Luxembourg" = "Wallonia",
  "Namur" = "Wallonia",
  "VlaamsBrabant" = "Flanders",
  "OostVlaanderen" = "Flanders",
  "WestVlaanderen" = "Flanders",
  "unknown" = "unknown"
)

#----------------------------------
# Define bank holidays.
bankholidays = as.Date(c("2020-04-13", # public holidays (excluding the ones that fell in a weekend), regional holidays excluded
                         "2020-05-01",
                         "2020-05-21",
                         "2020-06-01",
                         "2020-07-21",
                         "2020-11-01",
                         "2020-11-11",
                         "2020-12-24", # not an official one though
                         "2020-12-25"))
bankholiday = function(date) { date=as.Date(date)
                               as.factor(c("no","yes")[(date %in% bankholidays)*1+1]) }

#----------------------------------
# Translates provinces to regions.
add_region <- function(x){
  regions[x]
}

#----------------------------------
# Replace NA by "unknown"
add_unknown <- function(x){
  x[is.na(x)] <- "unknown"
  return(x)
}

#----------------------------------
# Convert the date and cut the last day

cleanDate <- function(x, var = "DATE"){
  tmp <- as.Date(x[[var]])
  x[[var]] <- tmp
  id <- !is.na(tmp) & tmp < (Sys.Date() - 1)
  x[id,]
}

#---------------------------------
# Replace missing values by 0
replaceby0 <- function(x){
  x[is.na(x)] <- 0
  x
}

#---------------------------------
# Calculate the change compared to n values ago

calculate_change <- function(x,n=7,type=c("abs","rel")){
  type <- match.arg(type)
  nx <- length(x)
  if(n >= nx) stop("x isn't large enough.")
  tmp <- x[(n+1):nx] - x[1:(nx - n)]
  if(type == "rel"){
    tmp <- tmp / x[1:(nx - n)]
  }
  return(c(rep(NA,n),tmp))
}

#-------------------------------
# Clean up impossible combinations

keep_combs<- function(region,province){
  poss <- paste(names(regions), regions, sep = "-")
  tmp <- paste(province, region, sep = "-")
  id <- tmp %in% poss | grepl("All-",tmp)
}

#-------------------------------
# categorize cases
bin_cases <- function(x, breaks = c(0,5,10,15,20,30,50,100,+Inf),
                      labels = c("<5","5-9","10-14","15-19",
                                 "20-29","30-49","50-99","100+")){
  idless <- x == "<5"
  x[idless] <- NA
  x <- cut(as.numeric(x), 
           breaks = breaks, labels = labels,
           right = FALSE)
  x[idless] <- "<5"
  x
}

# Calculate the change
changeabsolute <- function(x){
  n <- length(x)
  c(rep(NA,7), (x[8:n] - x[1:(n-7)]))
}

changepercent <- function(x){
  n <- length(x)
  c(rep(NA,7), (x[8:n] - x[1:(n-7)])/ x[1:(n-7)] * 100 )
}
