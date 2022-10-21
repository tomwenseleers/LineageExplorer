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

# here instead loading country by country
country = "Austria" # example country
con_in = gzfile(file.path(target_dir, 'provision.json.xz'))
outfile = file.path(target_dir, paste0("provision_",country,".json.gz"))
con_out = gzfile(outfile, open="wb") # file(outfile, open = "wb")
system.time(stream_in(con_in, handler = function(df){
  df <- dplyr::filter(df, grepl(country,covv_location)) # keep only records from specified country
  # TO DO: do all processing in chunks & return processed result - i.e. removing non-baseline surveillance records, lineage calling based on pangolin lineage & AA substitutions & aggregation per week
  stream_out(df, con_out, pagesize = 10000)
}, pagesize = 100000)) # 1931s=32 mins
close(con_out)
# stream it back in
system.time(GISAID_JSON <- stream_in(file(outfile))) # 192410 records
nrow(GISAID_JSON) # 192 410

# note: for Austrian data records from Lifebrain are biased by pre-selected sample sequencing
# covv_virus_name for those contains LB- in the name

head(GISAID_JSON)

names(GISAID_JSON)
# [1] "covsurver_existingmutlist" "covsurver_prot_mutations"  "covsurver_uniquemutlist"  
# [4] "covv_accession_id"         "covv_add_host_info"        "covv_add_location"        
# [7] "covv_clade"                "covv_collection_date"      "covv_gender"              
# [10] "covv_host"                 "covv_last_vaccinated"      "covv_lineage"             
# [13] "pangolin_lineages_version" "covv_location"             "covv_passage"             
# [16] "covv_patient_age"          "covv_patient_status"       "covv_sampling_strategy"   
# [19] "covv_specimen"             "covv_subm_date"            "covv_treatment"           
# [22] "covv_type"                 "covv_virus_name"           "gc_content"               
# [25] "is_complete"               "is_reference"              "is_high_coverage"         
# [28] "n_content"                 "sequence_length"           "covv_variant" 

unique(GISAID_JSON$covv_sampling_strategy)
# [1] ""                               "unknown"                       
# [3] "Baseline surveillance"          "Preselected sample"            
# [5] "Random sample"                  "random"                        
# [7] "Active surveillance"            "Other"                         
# [9] "Random samples"                 "Sentinel surveillance (ILI)"   
# [11] "Randomly selected"              "Sentinel surveillance ILI"     
# [13] "Same-patient sampling strategy" "SURVEILLANCE"                  
# [15] "Outbreak investigation" 

# GISAID_JSON$covv_virus_name

GISAID_JSON_baseline = GISAID_JSON[GISAID_JSON$covv_sampling_strategy %in% c("Baseline surveillance",
                                                                             "Random sample", 
                                                                             "random",
                                                                             "Random samples", 
                                                                             "Randomly selected",
                                                                             "SURVEILLANCE"),]
head(GISAID_JSON_baseline$covv_virus_name)
nrow(GISAID_JSON_baseline) # 7507
GISAID_JSON_baseline$covv_virus_name # IMBA in virus name
GISAID_JSON_baseline$covv_virus_name[!grepl("IMBA|LB",GISAID_JSON_baseline$covv_virus_name)]
GISAID_JSON_baseline$covv_virus_name[!grepl("IMBA",GISAID_JSON_baseline$covv_virus_name)]

sum(grepl("IMBA",GISAID_JSON_baseline$covv_virus_name)) # 4871
sum(grepl("IMBA",GISAID_JSON$covv_virus_name)) # 143367
nrow(GISAID_JSON) # 192410

GISAID_IMBA = GISAID_JSON[grepl("IMBA",GISAID_JSON$covv_virus_name),]
unique(GISAID_IMBA$covv_sampling_strategy)
sum(GISAID_IMBA$covv_sampling_strategy=="random") # 4871
tail(GISAID_IMBA$covv_virus_name[GISAID_IMBA$covv_sampling_strategy=="random"])
tail(GISAID_IMBA$covv_virus_name[GISAID_IMBA$covv_sampling_strategy!="random"])

GISAID_JSON_baseline[!grepl("IMBA|LB", GISAID_JSON_baseline$covv_virus_name),
                     "covv_virus_name"]


# IMBA can be filtered as "Elling laboratory/IMBA/OEAW and Institut für Infektionsepidemiologie & Surveillance/AGES" in the field 'submitting lab'. 
# The 'sampling strategy' should be also added to metadata, but it often isn't. 
# IMBA puts "random" in that field, it seems.
unique(GISAID_JSON$)
