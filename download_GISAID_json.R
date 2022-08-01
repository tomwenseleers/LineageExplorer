# DOWNLOAD GISAID json METADATA 
# T. Wenseleers, 29 JULY 2022

library(curl)
library(glue)

# note: set GISAID JSON access credentials first using ####
# Sys.setenv(GISAIDJSON_USERNAME = "XXX")
# Sys.setenv(GISAIDJSON_PASSWORD = "XXX")
# Sys.setenv(GISAIDJSON_STREAM = "XXX")
source(".//set_GISAID_credentials.R")

custom_curl <- function(user, pass, feed, dest) {
  custom_handle <- curl::new_handle()
  curl::handle_setopt(
    custom_handle,
    username = user,
    password = pass
  )
  
  url <- glue::glue(
    "https://www.epicov.org/epi3/3p/{feed}/export/provision.json.xz"
  )
  
  curl::curl_download(url, dest, handle = custom_handle)
}

message("Downloading GISAID JSON...")
custom_curl(Sys.getenv("GISAIDJSON_USERNAME"), 
            Sys.getenv("GISAIDJSON_PASSWORD"), 
            Sys.getenv("GISAIDJSON_STREAM"), 
            'data/GISAID_json/provision.json.xz')
message("GISAID JSON download complete")

