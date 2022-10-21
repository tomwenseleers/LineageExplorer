

# PS best to run this on a workstation with a lot of memory
# takes ca half an hour to parse
library(jsonlite)
library(dplyr)

message("Parsing JSON GISAID records...")

start_time <- Sys.time()

# If shell commands are used:
# cmd <- "xz -d -c ./data/gisaid-austria.txt.xz | grep \"Austria / Vienna\" "
# then needed to pass the command as pipe:     stream_in( pipe( cmd ),

#filename <- "./data/provision.json.xz"
filename <- "./data/gisaid-austria.txt.xz"

handle_one_df_page <- function (df) {
  idx <- as.character(length(GISAID_json) + 1)
  # here a chance to create different
  # (smaller) data than original
  t <- df[
    # ==================rows
	df$covv_host=="Human"
    ,  # ============ columns
	c(
	"covv_accession_id",
	"covv_virus_name",
	"covv_sampling_strategy",
	"covv_add_host_info",
	"covv_add_location",
	"covv_lineage"
	)
  ]
  GISAID_json[[ idx ]] <- t
}

GISAID_json <- new.env()
stream_in( gzfile( filename ),
	handler = handle_one_df_page,
	pagesize = 10000 )
GISAID_json <- bind_rows(as.list(( GISAID_json )))

str( GISAID_json )
names( GISAID_json )

end_time <- Sys.time()
end_time - start_time

nrow(GISAID_json) # data frame with 12025647 rows

# BEWARE as already filtered before and no covv_host,
# nothing will remain if uncommented!
# only keep human samples
# GISAID_json = GISAID_json[GISAID_json$covv_host=="Human",]
print( "only humans" )
nrow(GISAID_json) # data frame with 11853650  rows
print( "Names" )
names(GISAID_json)
print( "Head covv_location" )
head(GISAID_json$covv_location)
#print( "H covv_collection_date" )
#head(GISAID_json$covv_collection_date)
#print( "H covv_lineage" )
head(GISAID_json$covv_lineage)
print( "Unique covv_lineage" )
unique(GISAID_json$covv_lineage)
#head(GISAID_json$covsurver_existingmutlist)
#head(GISAID_json$covsurver_prot_mutations)
#head(GISAID_json$covsurver_uniquemutlist)
unique(GISAID_json$covv_sampling_strategy)

print( "Object size" )
object.size(GISAID_json)

#library(dplyr)
GISAID_json = mutate_at(GISAID_json,
                        "covv_sampling_strategy", .funs=toupper)
GISAID_json = mutate_at(GISAID_json,
                        "covv_add_host_info", .funs=toupper)
GISAID_json = mutate_at(GISAID_json,
                        "covv_add_location", .funs=toupper)
