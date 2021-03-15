#----------------------------------------
# Functions to convert the descriptive of 
# the administrative units. Use for processing
# the shape info and making names consistent.

library(stringr)

translate_regions <- function(x){
  recode(x,
         "Vlaams Gewest" = "Flanders",
         "Brussels Hoofdstedelijk Gewest" = "Brussels",
         "Waals Gewest" = "Wallonia")
}

translate_districts <- function(x_nl, x_fr,region){
  case_when(
    region == "Flanders" ~ str_remove(x_nl, "Arrondissement "),
    region == "Wallonia" ~ str_remove(x_fr, "Arrondissement d[eu']\\s?"),
    region == "Brussels" ~ str_remove(x_nl,"Arrondissement ")
  )
}

translate_provinces <- function(x_nl, x_fr, region){
  case_when(
    region == "Flanders" ~ str_remove(x_nl,"Provincie "),
    region == "Wallonia" ~ str_remove(x_fr,"Province d[eu] "),
    region == "Brussels" ~ "Brussel"
  )
}

translate_munty <- function(x_nl, x_fr,region){
  case_when(
    region == "Flanders" ~ x_nl,
    region == "Wallonia" ~ x_fr,
    region == "Brussels" ~ x_nl
  )
}

# Remove dashes and spaces, and capitalize every first letter
clean_prov <- function(x){
  str_replace(x, "-", " ") %>%
    gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", ., perl=TRUE) %>%
    str_remove(" ")
}
