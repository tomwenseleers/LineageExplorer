# IMPORT GISAID json METADATA 
# T. Wenseleers, 28 January 2022

# to download latest GISAID JSON run
# source(".//download_GISAID_json.R")

# PS takes ca an hour to load on my laptop, better run this on a workstation with a lot of memory

library(jsonlite)
message("Parsing JSON GISAID records...")
GISAID_json = jsonlite::stream_in(gzfile(".//data//GISAID_json//provision.json.xz")) 
saveRDS(GISAID_json, ".//data/GISAID_json//GISAID_json.rds")
nrow(GISAID_json) # data frame with 4425000 rows
library(inspectdf)
memcons = data.frame(inspect_mem(GISAID_json))
memcons
# PS most memory taken up by these columns
#                     col_name     bytes      size      pcnt
# 1   covsurver_prot_mutations 665094432 634.28 Mb 23.158827
# 2  covsurver_existingmutlist 627688440 598.61 Mb 21.856337
# 3            covv_virus_name 427700024 407.89 Mb 14.892668
# 4          covv_accession_id 308555832 294.26 Mb 10.744025


# only keep human samples
GISAID_json = GISAID_json[GISAID_json$covv_host=="Human",] 
nrow(GISAID_json) # data frame with 4285497  rows
names(GISAID_json)
head(GISAID_json$covv_location)
head(GISAID_json$covv_collection_date)
head(GISAID_json$covv_lineage)
unique(GISAID_json$covv_lineage)
head(GISAID_json$covsurver_existingmutlist)
head(GISAID_json$covsurver_prot_mutations)
head(GISAID_json$covsurver_uniquemutlist)
unique(GISAID_json$covv_sampling_strategy)

library(dplyr)
GISAID_json = mutate_at(GISAID_json, "covv_sampling_strategy", .funs=toupper)
GISAID_json = mutate_at(GISAID_json, "covv_add_host_info", .funs=toupper)
GISAID_json = mutate_at(GISAID_json, "covv_add_location", .funs=toupper)

# remove travel-related cases, active surveillance, surge testing etc
GISAID_json = GISAID_json[!grepl("ACTIVE|S GENE|COVIGEN|S-GENE|SUSPECT|LONGITUDINAL|ASSAY|DROPOUT|CLUSTER|OUTBREAK|LONGITUDINAL|S GENE|SAME-PATIENT|TRAVEL|VOC|VARIANT|CONTACT|TRACING|PCR TEST|RETURNING|SGTF|TRIAL|FAMILY|DROPOUT|TAQPATH|WEEKLY|NON-RANDOM|SAME PATIENT|TIME-SERIES|TARGET|QUARANTINE|PASSAGE",
                                 GISAID_json$covv_sampling_strategy),]
nrow(GISAID_json) # 4211029 rows
GISAID_json = GISAID_json[!grepl("ACTIVE|S GENE|COVIGEN|S-GENE|SUSPECT|LONGITUDINAL|ASSAY|DROPOUT|CLUSTER|OUTBREAK|LONGITUDINAL|S GENE|SAME-PATIENT|TRAVEL|VOC|VARIANT|CONTACT|TRACING|PCR TEST|RETURNING|SGTF|TRIAL|FAMILY|DROPOUT|TAQPATH|WEEKLY|NON-RANDOM|SAME PATIENT|TIME-SERIES|TARGET|QUARANTINE|PASSAGE|PASSENGER|INFECTED IN|HOLIDAY|CONTACT|CAME FROM|RETURN|FROM|VISIT|FOREIGN|SKIING|IMPORT|RELATED WITH|ENTERED FROM|CASE FROM",
                                 GISAID_json$covv_add_host_info),]
nrow(GISAID_json) # 4202123 rows
GISAID_json = GISAID_json[!grepl("ACTIVE|S GENE|COVIGEN|S-GENE|SUSPECT|LONGITUDINAL|ASSAY|DROPOUT|CLUSTER|OUTBREAK|LONGITUDINAL|S GENE|SAME-PATIENT|TRAVEL|VOC|VARIANT|CONTACT|TRACING|PCR TEST|RETURNING|SGTF|TRIAL|FAMILY|DROPOUT|TAQPATH|WEEKLY|NON-RANDOM|SAME PATIENT|TIME-SERIES|TARGET|QUARANTINE|PASSAGE|PASSENGER|INFECTED IN|HOLIDAY|CONTACT|CAME FROM|RETURN|FROM|VISIT|FOREIGN|SKIING|IMPORT|RELATED WITH|ENTERED FROM|CASE FROM|HISTOR|SAMPLED AT|TOURIST|VISIT",
                                 GISAID_json$covv_add_location),]
nrow(GISAID_json) # 4195414 rows

# parse dates & remove records with invalid dates
library(stringr)
date_isvalid = sapply(GISAID_json$covv_collection_date, function (s) str_count(s, pattern = "-")==2)
GISAID_json = GISAID_json[date_isvalid,]
GISAID_json$date = as.Date(GISAID_json$covv_collection_date) 
GISAID_json = GISAID_json[!is.na(GISAID_json$date),]
nrow(GISAID_json) # 4068724

# add numeric version of date, week midpoint & week & year

# GISAID_json = GISAID_json[GISAID_json$date>=as.Date("2021-01-01"),]
range(GISAID_json$date) # "2019-12-24" "2021-10-11"

GISAID_json$Week = lubridate::week(GISAID_json$date)
GISAID_json$Year = lubridate::year(GISAID_json$date)
GISAID_json$Year_Week = interaction(GISAID_json$Year,GISAID_json$Week)
library(lubridate)
GISAID_json$week_midpoint = as.Date(as.character(cut(GISAID_json$date, "week")))+3.5 # week midpoint date
GISAID_json$DATE_NUM = as.numeric(GISAID_json$date)
names(GISAID_json)


# parse continent / country / location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID_json$covv_location), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
unique(loc[,1]) # continent
unique(loc[,2]) # country
unique(loc[,3]) # city or province

levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
GISAID_json$continent = factor(loc[,1], levels=levels_continents)
GISAID_json$country = factor(loc[,2])
levels_countries = levels(GISAID_json$country)
GISAID_json$location = factor(loc[,3])
levels_locations = levels(GISAID_json$location)



# recode pango lineages as desired
GISAID_json = GISAID_json[!(GISAID_json$covv_lineage=="None"|GISAID_json$covv_lineage==""|is.na(GISAID_json$covv_lineage)),]
nrow(GISAID_json) # 3981738
options(max.print=1000000)
unique(GISAID_json$covv_lineage)
length(unique(GISAID_json$covv_lineage)) # 1393 pango lineages
# TO DO: allow user to define custom grepl patterns to recode & name lineages (based on combinations of pango lineages & mutations?)

# GISAID$Lineage[grepl("B.1.177",GISAID$Lineage,fixed=T)] = "B.1.177+"
# GISAID$Lineage[grepl("B.1.36\\>",GISAID$Lineage)] = "B.1.36+"
# 
# sel_target_VOC = "B.1.617"
# GISAID$LINEAGE1 = GISAID$Lineage
# GISAID$LINEAGE2 = GISAID$Lineage
# GISAID[grepl(sel_target_VOC, GISAID$LINEAGE1, fixed=TRUE),"LINEAGE1"] = paste0(sel_target_VOC,"+") # in LINEAGE1 we recode B.1.617.1,2&3 all as B.1.617+
# 
# # subset to selected countries
# GISAID_json = GISAID[GISAID$Country %in% sel_countries,]
# nrow(GISAID_json) # 7691
# 
# rowSums(table(GISAID_json$LINEAGE1,GISAID_json$Country))
# 
# main_lineages = names(table(GISAID_json$LINEAGE1))[100*table(GISAID_json$LINEAGE1)/sum(table(GISAID_json$LINEAGE1)) > 1]
# main_lineages
# # "A.23.1"    "B.1"       "B.1.1.318" "B.1.1.7"   "B.1.351"   "B.1.525"   "B.1.617+"
# VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
#          "B.1.525","B.1.526","B.1.1.519")
# main_lineages = union(main_lineages, VOCs)
# GISAID_json$LINEAGE1[!(GISAID_json$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
# GISAID_json$LINEAGE2[!(GISAID_json$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
# remove1 = names(table(GISAID_json$LINEAGE1))[table(GISAID_json$LINEAGE1)/sum(table(GISAID_json$LINEAGE1)) < 0.01]
# remove1 = remove1[!(remove1 %in% c("B.1.351","B.1.1.7","P.1","B.1.617+","B.1.1.519"))]
# remove2 = names(table(GISAID_json$LINEAGE2))[table(GISAID_json$LINEAGE2)/sum(table(GISAID_json$LINEAGE2)) < 0.01]
# remove2 = remove2[!(remove2 %in% c("B.1.351","B.1.1.7","P.1","B.1.617.2","B.1.617.1","B.1.1.519"))]
# GISAID_json$LINEAGE1[(GISAID_json$LINEAGE1 %in% remove1)] = "other" # minority VOCs
# GISAID_json$LINEAGE2[(GISAID_json$LINEAGE2 %in% remove2)] = "other" # minority VOCs
# table(GISAID_json$LINEAGE1)
# GISAID_json$LINEAGE1 = factor(GISAID_json$LINEAGE1)
# GISAID_json$LINEAGE1 = relevel(GISAID_json$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
# levels(GISAID_json$LINEAGE1)
# levels_LINEAGE1 = c("B.1.1.7",levels(GISAID_json$LINEAGE1)[!levels(GISAID_json$LINEAGE1) %in% c("B.1.1.7","B.1.617+","B.1.617.1","B.1.617.2","other")],
#                     "B.1.617+","other")
# GISAID_json$LINEAGE1 = factor(GISAID_json$LINEAGE1, levels=levels_LINEAGE1)
# 
# GISAID_json$LINEAGE2 = factor(GISAID_json$LINEAGE2)
# GISAID_json$LINEAGE2 = relevel(GISAID_json$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
# levels(GISAID_json$LINEAGE2)
# levels_LINEAGE2 = c("B.1.1.7",levels(GISAID_json$LINEAGE2)[!levels(GISAID_json$LINEAGE2) %in% c("B.1.1.7","B.1.617+","B.1.617.1","B.1.617.2","other")],
#                     "B.1.617.1","B.1.617.2","other")
# GISAID_json$LINEAGE2 = factor(GISAID_json$LINEAGE2, levels=levels_LINEAGE2)
# 
# GISAID_json$country = factor(GISAID_json$Country, levels=levels_countries)
# 
# table(GISAID_json$Country)


# code main clades & VOCs only
unique(GISAID_json$covv_clade)
GISAID_json$VOC = GISAID_json$covv_variant
GISAID_json$VOC[GISAID_json$covv_lineage=="A"] = "A" 
GISAID_json$VOC[GISAID_json$VOC==""&(grepl("S",GISAID_json$covv_clade)|grepl("A.",GISAID_json$covv_lineage,fixed=T))] = "A*" # = 19B, cf https://www.news-medical.net/health/Viral-Clades-of-SARS-CoV-2.aspx
GISAID_json$VOC[GISAID_json$covv_lineage=="B"] = "B"
GISAID_json$VOC[GISAID_json$VOC==""&(grepl("L|O|V|G",GISAID_json$covv_clade)|grepl("B.",GISAID_json$covv_lineage,fixed=T))&(!grepl("D614G",GISAID_json$covsurver_prot_mutations))] = "B*" # B* without D614G
GISAID_json$VOC[GISAID_json$VOC==""&(grepl("G",GISAID_json$covv_clade)|grepl("B.",GISAID_json$covv_lineage,fixed=T)|grepl("AE.",GISAID_json$covv_lineage,fixed=T)|grepl("AH.",GISAID_json$covv_lineage,fixed=T)|grepl("AK.",GISAID_json$covv_lineage,fixed=T)|grepl("AN.",GISAID_json$covv_lineage,fixed=T)|grepl("AS.",GISAID_json$covv_lineage,fixed=T)|grepl("AU.",GISAID_json$covv_lineage,fixed=T)|grepl("AV.",GISAID_json$covv_lineage,fixed=T)|grepl("AZ.",GISAID_json$covv_lineage,fixed=T))&(grepl("D614G",GISAID_json$covsurver_prot_mutations))] = "B* (+S:D614G)" # B* with D614G
GISAID_json$VOC[grepl("Alpha",GISAID_json$VOC)] = "Alpha"
GISAID_json$VOC[grepl("Beta",GISAID_json$VOC)] = "Beta"
GISAID_json$VOC[grepl("Gamma",GISAID_json$VOC)|grepl("P.1.",GISAID_json$covv_lineage,fixed=T)] = "Gamma"
GISAID_json$VOC[grepl("Delta",GISAID_json$VOC)|grepl("AY.",GISAID_json$covv_lineage,fixed=T)] = "Delta"
GISAID_json$VOC[grepl("Lambda",GISAID_json$VOC)] = "Lambda"
GISAID_json$VOC[grepl("Mu",GISAID_json$VOC)] = "Mu"
table(GISAID_json$VOC)
table(GISAID_json$covv_lineage[GISAID_json$VOC==""])  
# sum(GISAID_json$VOC=="") # 0 remaining unassigned
# GISAID_json = GISAID_json[GISAID_json$VOC!="",] # we remove these samples; instead we could also recode them as category "other"
nrow(GISAID_json) # 3981577

table(GISAID_json$covv_lineage[GISAID_json$covv_location=="Asia / China / Hubei / Wuhan"])

# use main clades/VOCs as my LINEAGE variable
levels_LINEAGE = c("A", "A*", "B", "B*", "B* (+S:D614G)", "Alpha", "Beta", "Gamma", "Delta", "Lambda", "Mu")
GISAID_json$LINEAGE = factor(GISAID_json$VOC, 
                             levels=levels_LINEAGE)

levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
levels_countries = levels(GISAID_json$country)


# saveRDS(object=GISAID_json, file=".//data//GISAID_json//GISAID_json.rds")
# GISAID_json = readRDS(file=".//data//GISAID_json//GISAID_json.rds")

# ANALYSIS OF MAIN CLADES & VOCs ####

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek = as.data.frame(table(GISAID_json$week_midpoint, GISAID_json$LINEAGE))
colnames(data_agbyweek) = c("week_midpoint", "LINEAGE", "count")
data_agbyweek_sum = aggregate(count ~ week_midpoint, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$week_midpoint, data_agbyweek_sum$week_midpoint)]
sum(data_agbyweek[data_agbyweek$LINEAGE=="Alpha","total"]) == nrow(GISAID_json) # correct
data_agbyweek$collection_date = as.Date(as.character(data_agbyweek$week_midpoint))
data_agbyweek$LINEAGE = factor(data_agbyweek$LINEAGE, levels=levels_LINEAGE)
data_agbyweek$collection_date_num = as.numeric(data_agbyweek$collection_date)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
data_agbyweek$week_midpoint = NULL

# aggregated by week & continent for selected variant lineages
data_agbyweek_bycontinent = as.data.frame(table(GISAID_json$week_midpoint, GISAID_json$continent, GISAID_json$LINEAGE))
colnames(data_agbyweek_bycontinent) = c("week_midpoint", "continent", "LINEAGE", "count")
data_agbyweek_bycontinent_sum = aggregate(count ~ week_midpoint+continent, data=data_agbyweek_bycontinent, sum)
data_agbyweek_bycontinent$total = data_agbyweek_bycontinent_sum$count[match(interaction(data_agbyweek_bycontinent$week_midpoint,data_agbyweek_bycontinent$continent), 
                                                                            interaction(data_agbyweek_bycontinent_sum$week_midpoint,data_agbyweek_bycontinent_sum$continent))]
sum(data_agbyweek_bycontinent[data_agbyweek_bycontinent$LINEAGE=="Alpha","total"]) == nrow(GISAID_json) # correct
data_agbyweek_bycontinent$collection_date = as.Date(as.character(data_agbyweek_bycontinent$week_midpoint))
data_agbyweek_bycontinent$LINEAGE = factor(data_agbyweek_bycontinent$LINEAGE, levels=levels_LINEAGE)
data_agbyweek_bycontinent$collection_date_num = as.numeric(data_agbyweek_bycontinent$collection_date)
data_agbyweek_bycontinent$prop = data_agbyweek_bycontinent$count/data_agbyweek_bycontinent$total
data_agbyweek_bycontinent$week_midpoint = NULL
data_agbyweek_bycontinent$continent = factor(data_agbyweek_bycontinent$continent, levels=levels_continents)

# aggregated by week, continent & country for selected variant lineages
data_agbyweek_bycontinent_country = as.data.frame(table(GISAID_json$week_midpoint, 
                                                        GISAID_json$continent, 
                                                        GISAID_json$country, 
                                                        GISAID_json$LINEAGE))
colnames(data_agbyweek_bycontinent_country) = c("week_midpoint", "continent", "country", "LINEAGE", "count")
data_agbyweek_bycontinent_country_sum = aggregate(count ~ week_midpoint+continent+country, 
                                                  data=data_agbyweek_bycontinent_country, sum)
data_agbyweek_bycontinent_country$total = data_agbyweek_bycontinent_country_sum$count[
  match(interaction(data_agbyweek_bycontinent_country$week_midpoint,data_agbyweek_bycontinent_country$continent,data_agbyweek_bycontinent_country$country), 
        interaction(data_agbyweek_bycontinent_country_sum$week_midpoint,data_agbyweek_bycontinent_country_sum$continent,data_agbyweek_bycontinent_country_sum$country))]
sum(data_agbyweek_bycontinent_country[data_agbyweek_bycontinent_country$LINEAGE=="Alpha","total"]) == nrow(GISAID_json) # correct
data_agbyweek_bycontinent_country$collection_date = as.Date(as.character(data_agbyweek_bycontinent_country$week_midpoint))
data_agbyweek_bycontinent_country$LINEAGE = factor(data_agbyweek_bycontinent_country$LINEAGE, levels=levels_LINEAGE)
data_agbyweek_bycontinent_country$collection_date_num = as.numeric(data_agbyweek_bycontinent_country$collection_date)
data_agbyweek_bycontinent_country$prop = data_agbyweek_bycontinent_country$count/data_agbyweek_bycontinent_country$total
data_agbyweek_bycontinent_country$week_midpoint = NULL
data_agbyweek_bycontinent_country$continent = factor(data_agbyweek_bycontinent_country$continent, levels=levels_continents)
data_agbyweek_bycontinent_country$country = factor(data_agbyweek_bycontinent_country$country, levels=levels_countries)
