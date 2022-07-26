# DOWNLOAD GISAID json METADATA 
# T. Wenseleers, 12 JULY 2022

library(curl)
library(glue)

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
custom_curl('u0008355', 'f32oF9ASmWAsX', 'kuleuven', 
            'data/GISAID_json/provision.json.xz')
message("GISAID JSON download complete")

