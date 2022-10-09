# DOWNLOAD GISAID json METADATA 
# T. Wenseleers, 9 OCTOBER 2022

# TO DO: I should wrap this in a function like in download_GISAID.R
# make the names of the fields consistent with those returned by 
# download_GISAID.R
# and maybe apply some filtering to filter out targeted sequencing & sequencing
# of travel related cases?

library(curl)
library(glue)

target_dir = "C:/Users/bherr/OneDrive - KU Leuven/Documents/Github/LineageExplorer/LineageExplorer/data/GISAID" # target download directory GISAID data

# note: set GISAID JSON access credentials first using ####
# Sys.setenv(GISAIDJSON_USERNAME = "XXX")
# Sys.setenv(GISAIDJSON_PASSWORD = "XXX")
# Sys.setenv(GISAIDJSON_STREAM = "XXX")
source("..//set_GISAID_credentials.R")

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
            file.path(target_dir, 'provision.json.xz') )
message("GISAID JSON download complete")

message("Parsing JSON GISAID records...")
library(dplyr)

library(jsonlite)
# GISAID_json = jsonlite::stream_in(gzfile(file.path(target_dir, 'provision.json.xz'))) 
# this runs out of memory on my 64 Gb laptop :-(

# trying instead like this
con_in = gzfile(file.path(target_dir, 'provision.json.xz'))
outfile = file.path(target_dir, 'provision.json.bin')
con_out = file(out <- outfile, open = "wb")
stream_in(con_in, handler = function(df){
  # df <- dplyr::filter(df, XXX) # drop any unneeded columns
  stream_out(df, con_out, pagesize = 100000)
}, pagesize = 100000)
close(con_out)
# stream it back in
GISAID_JSON = stream_in(file(out))
nrow(GISAID_JSON)


