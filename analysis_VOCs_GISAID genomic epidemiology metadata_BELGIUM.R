# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN BELGIUM BASED ON ANALYSIS OF GISAID PATIENT METADATA
# T. Wenseleers
# last update 24 MAY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-24")
today_num = as.numeric(today)
today # "2021-05-24"
plotdir = "VOCs_belgium"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import & parse manually downloaded GISAID patient metadata for Belgium (downloaded 24/5/2021) ####
d1 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_selection_2020.tsv", col_types = cols(.default = "c"))
d2 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_selection_jan_feb_2021.tsv", col_types = cols(.default = "c"))
d3 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_selection_mar_apr_2021.tsv", col_types = cols(.default = "c"))
d4 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_selection_may_2021.tsv", col_types = cols(.default = "c"))
d1 = as.data.frame(d1)
d2 = as.data.frame(d2)
d3 = as.data.frame(d3)
d4 = as.data.frame(d4)
GISAID_belgium1 = rbind(d1, d2, d3, d4)
GISAID_belgium1 = GISAID_belgium1[GISAID_belgium1[,"Host"]=="Human",]
GISAID_belgium1[,"Collection date"] = as.Date(GISAID_belgium1[,"Collection date"])
GISAID_belgium1 = GISAID_belgium1[!is.na(GISAID_belgium1[,"Collection date"]),]
nrow(GISAID_belgium1) # 21 449
unique(GISAID_belgium1[,"Virus name"])
unique(GISAID_belgium1[,"Location"])
unique(GISAID_belgium1[,"Host"])
unique(GISAID_belgium1[,"Additional location information"]) # ZIP codes # sometimes has additional info on whether it was a traveller

# anonymous records without location info are assumed to be active surveillance of cluster outbreaks & included in active surveillance (cf weekly Sciensano report)
GISAID_belgium1[,"Additional location information"][grepl("nknown",GISAID_belgium1[,"Additional location information"])] = NA
sum(is.na(GISAID_belgium1[,"Additional location information"])) # 6064
unique(GISAID_belgium1[,"Additional location information"])

library(stringr)
ZIP = str_extract(GISAID_belgium1[,"Additional location information"],"[0-9]{4}") # extract ZIP codes
unique(ZIP)
sum(!is.na(ZIP)) # 15312
sum(is.na(ZIP)) # 6137 - anonymous records, these are assumed to be active surveillance of cluster outbreaks

muni = read.csv(".//data//GISAID//Belgium//mapping_municips.csv") # post codes & cities & provinces
province = muni$province_label_nl[match(ZIP, muni$postcode)] # province
city = muni$municipality_label_nl[match(ZIP, muni$postcode)] # city
sum(is.na(city)) # 6164  
sum(is.na(ZIP))  # 6137
wrongZIP = which(is.na(city)&!is.na(ZIP))
ZIP[wrongZIP] = NA

anonym = grepl("Belgium$",GISAID_belgium1[,"Additional location information"])|
  is.na(GISAID_belgium1[,"Additional location information"])|
  grepl("CHU|infected|travel|Travel",GISAID_belgium1[,"Additional location information"]) |
  is.na(ZIP)
sum(anonym) # 6168
unique(GISAID_belgium1[!is.na(GISAID_belgium1[,"Sampling strategy"]),"Sampling strategy"])
traveller = grepl("travel|Travel|infected",GISAID_belgium1[,"Additional location information"])|grepl("travel|Travel",GISAID_belgium1[,"Sampling strategy"])|
  grepl("travel|Travel|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(traveller) # 143
Sdropout = grepl("S-gene dropout",GISAID_belgium1[,"Sampling strategy"])|
  grepl("S-gene dropout",GISAID_belgium1[,"Additional host information"])  
sum(Sdropout) # 61 - should be much higher though...
actsurv = anonym|traveller|grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster",GISAID_belgium1[,"Sampling strategy"])|
  grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster|HR contact|Outbreak|Nurse|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(actsurv) # 6797
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")]) # 4135 # In Sciensano weekly report of 21st of May 2021 this is 4710
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")&GISAID_belgium1[,"Collection date"]<=as.Date("2021-01-31")]) # 1172 # In Sciensano weekly report of 21st of May 2021 this is 1975

reinfection = grepl("Breakthrough",GISAID_belgium1[,"Sampling strategy"])
sum(reinfection) # 7
infectionpostvaccination = grepl("Vaccination|vaccination",GISAID_belgium1[,"Sampling strategy"])
sum(infectionpostvaccination) # 67

asymptomatic = grepl("Asymptomatic|asympto",GISAID_belgium1[,"Patient status"])
symptomatic = grepl("Ho|Out|Am|Dec|dec| sympto|Not hosp|Relea",GISAID_belgium1[,"Patient status"])|grepl(" symptoms|ILI",GISAID_belgium1[,"Additional host information"])
outpatient = grepl("Out|Am|Not hosp",GISAID_belgium1[,"Patient status"])
hospitalised = grepl("Ho|Dec|dec|Relea",GISAID_belgium1[,"Patient status"])
died = grepl("Dec|dec",GISAID_belgium1[,"Patient status"])
unknown = grepl("unknown",GISAID_belgium1[,"Patient status"])
asymptomatic[unknown] = NA
symptomatic[unknown] = NA
outpatient[unknown] = NA
hospitalised[unknown] = NA
died[unknown] = NA

# bassurv = grepl("Baseline|baseline|Basic|^surveillance$|^Surveillance$|entinel",GISAID_belgium1[,"Sampling strategy"])|
#   grepl("Baseline",GISAID_belgium1[,"Specimen"])|
#   grepl("unbiased|Baseline|entinel",GISAID_belgium1[,"Additional host information"])
# # Sentinel surveillance, Sentinel surveillance (ILI), Sentinel surveillance (SARI), Non-sentinel-surveillance (GP network),
# # Non-sentinel-surveillance, Non-sentinel-surveillance (hospital), Non-sentinel-surveillance (Hospital) are all baseline surveillance
# sum(bassurv) # 8856

bassurv = !actsurv
sum(bassurv) # 14652
purpose_of_sequencing = rep("active_surveillance", nrow(GISAID_belgium1))
purpose_of_sequencing[bassurv] = "baseline_surveillance"

accessionIDs_bassurv = GISAID_belgium1[,"Accession ID"][bassurv]

GISAID_belgium1[,"Gender"][grepl("F|V",GISAID_belgium1[,"Gender"])] = "F"
GISAID_belgium1[,"Gender"][grepl("M|m",GISAID_belgium1[,"Gender"])] = "M"
GISAID_belgium1[,"Gender"][grepl("unknown",GISAID_belgium1[,"Gender"])] = NA
unique(GISAID_belgium1[,"Gender"]) # "F" "M" NA 

sum(!is.na(GISAID_belgium1[,"Gender"] )) # 18984 with gender info

GISAID_belgium1[,"Patient age"][grepl("nknown",GISAID_belgium1[,"Patient age"])] = NA
GISAID_belgium1[,"Patient age"][grepl("month",GISAID_belgium1[,"Patient age"])] = 0
GISAID_belgium1[,"Patient age"][grepl("day",GISAID_belgium1[,"Patient age"])] = 0
GISAID_belgium1[,"Patient age"][grepl("2020-1954",GISAID_belgium1[,"Patient age"])] = as.character(2021-1954)
GISAID_belgium1[,"Patient age"][grepl("2020-1991",GISAID_belgium1[,"Patient age"])] = as.character(2021-1991)
GISAID_belgium1[,"Patient age"][grepl("2020-1984",GISAID_belgium1[,"Patient age"])] = as.character(2021-1984)
GISAID_belgium1[,"Patient age"][grepl("2020-1979",GISAID_belgium1[,"Patient age"])] = as.character(2021-1979)
GISAID_belgium1[,"Patient age"][grepl("2020-1994",GISAID_belgium1[,"Patient age"])] = as.character(2021-1994)
GISAID_belgium1[,"Patient age"][grepl("2020-1985",GISAID_belgium1[,"Patient age"])] = as.character(2021-1985)
GISAID_belgium1[,"Patient age"][grepl("2020-1972",GISAID_belgium1[,"Patient age"])] = as.character(2021-1972)
GISAID_belgium1[,"Patient age"] = as.numeric(GISAID_belgium1[,"Patient age"])
sum(!is.na(GISAID_belgium1[,"Patient age"] )) # 18995 with age info

GISAID_belgium1[,"Virus name"] = gsub("hCoV-19/","",GISAID_belgium1[,"Virus name"])

colnames(GISAID_belgium1) = c("strain",
                              "gisaid_epi_isl",
                              "date",
                              "location",
                              "host",
                              "additional_location_information",
                              "sampling_strategy",
                              "sex",                         
                              "age",
                              "patient_status",
                              "last_vaccinated",
                              "passage",
                              "specimen",
                              "additional_host_information",
                              "pango_lineage",
                              "GISAID_clade",
                              "AA_substitutions") 
# add some extra parsed fields
GISAID_belgium1$travel_history = c("no","yes")[traveller*1+1] # travel history?
GISAID_belgium1$reinfection = c("no","yes")[reinfection*1+1]
GISAID_belgium1$infectionpostvaccination = c("no","yes")[infectionpostvaccination*1+1]
GISAID_belgium1$asymptomatic = c("no","yes")[asymptomatic*1+1]
GISAID_belgium1$symptomatic = c("no","yes")[symptomatic*1+1]
GISAID_belgium1$outpatient = c("no","yes")[outpatient*1+1]
GISAID_belgium1$hospitalised = c("no","yes")[hospitalised*1+1]
GISAID_belgium1$died = c("no","yes")[died*1+1]
GISAID_belgium1$Sdropout_sequencing = c("no","yes")[Sdropout*1+1]
# add postcode, city & province
GISAID_belgium1$region = "Europe"
GISAID_belgium1$country = "Belgium"
GISAID_belgium1$postcode = as.integer(ZIP)
GISAID_belgium1$city = city
GISAID_belgium1$province = province
GISAID_belgium1$province = factor(GISAID_belgium1$province)
levels(GISAID_belgium1$province)
# "Antwerpen"                      "Brussels Hoofdstedelijk Gewest" "Henegouwen"                     "Limburg"                       
# "Luik"                           "Luxemburg"                      "Namen"                          "Oost-Vlaanderen"               
# "Vlaams-Brabant"                 "Waals-Brabant"                  "West-Vlaanderen"  
levels_PROVINCES = levels(GISAID_belgium1$province)

# add parsed purpose of sequencing
GISAID_belgium1$purpose_of_sequencing = purpose_of_sequencing
library(lubridate)
GISAID_belgium1$Week = lubridate::week(GISAID_belgium1$date)
GISAID_belgium1$Year = lubridate::year(GISAID_belgium1$date)
GISAID_belgium1$Year_Week = interaction(GISAID_belgium1$Year,GISAID_belgium1$Week)
GISAID_belgium1$floor_date = as.Date(as.character(cut(GISAID_belgium1$date, "week")))+3.5 # week midpoint date
GISAID_belgium1$DATE_NUM = as.numeric(GISAID_belgium1$date)

# recode some lineages
GISAID_belgium1$pango_lineage_simplified = GISAID_belgium1$pango_lineage
GISAID_belgium1$pango_lineage_simplified[grepl("B.1.177", GISAID_belgium1$pango_lineage_simplified,fixed=T)] = "B.1.177+"
GISAID_belgium1$pango_lineage_simplified[grepl("B.1.36\\>",GISAID_belgium1$pango_lineage_simplified)] = "B.1.36+"
GISAID_belgium1$LINEAGE1 = GISAID_belgium1$pango_lineage_simplified
GISAID_belgium1$LINEAGE2 = GISAID_belgium1$pango_lineage_simplified
GISAID_belgium1[grepl("B.1.617", GISAID_belgium1$LINEAGE1, fixed=TRUE),"LINEAGE1"] = paste0("B.1.617","+") # in LINEAGE1 we recode B.1.617.1,2&3 all as B.1.617+

# write parsed & cleaned up file to csv
write.csv(GISAID_belgium1, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_ALL PARSED.csv",row.names=F)


# PS fields "genbank_accession", "Nextstrain_clade", "originating_lab", "submitting_lab", "authors", "url", "title", "paper_url"
# also still given in GISAID genomic epidemiology file, but we don't need those right now
# (unless some labs were only doing active surveillance or only baseline surveillance, check with Emmanuel)
# genomic epidemiology GISAID data file version metadata_2021-05-21_13-00.tsv.gz with data for all countries :
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-05-21_13-00.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID_genepi_belgium = GISAID[GISAID$country=="Belgium",]
nrow(GISAID_genepi_belgium) # 20590
uniqueoriginatinglabs = unique(GISAID_genepi_belgium$originating_lab)
uniquesubmittinglabs = unique(GISAID_genepi_belgium$submitting_lab)
write.csv(uniqueoriginatinglabs,".//data//GISAID//Belgium//unique_originating_labs.csv")
write.csv(uniquesubmittinglabs,".//data//GISAID//Belgium//unique_submitting_labs.csv")



# ANALYSIS OF VOC LINEAGE FREQUENCIES IN BELGIUM USING MULTINOMIAL FITS ####

# GISAID_belgium = GISAID_sel[GISAID_sel$country=="Belgium",]
GISAID_belgium = GISAID_belgium1[GISAID_belgium1$purpose_of_sequencing=="baseline_surveillance",]
nrow(GISAID_belgium) # 14652
sum(!is.na(GISAID_belgium$age)) # 13612 with age info

# unique(GISAID_belgium$originating_lab) # 133 labs # maybe some labs did more active surveillance & could be excluded?
# unique(GISAID_belgium$submitting_lab) # 41 labs
# sum(GISAID_belgium$country_exposure!="Belgium") # only 23 are indicated as travellers

# # we remove travellers
# GISAID_belgium = GISAID_belgium[GISAID_belgium$country_exposure=="Belgium",]
# nrow(GISAID_belgium) # 19283
# sum(GISAID_belgium$division_exposure==GISAID_belgium$division) # 19283

unique(GISAID_belgium$province) # 
unique(GISAID_belgium$province[GISAID_belgium$LINEAGE1=="B.1.617+"]) # "Belgium"         "Hasselt"         "Brussels"        "Brugge"          "Gent"            "Liège"           "Halle-Vilvoorde" "Antwerpen"       "Mechelen"
sum(GISAID_belgium$LINEAGE1=="B.1.617+") # 19 among baseline surveillance
unique(GISAID_belgium$province[GISAID_belgium$LINEAGE1=="B.1.1.7"])
sum(GISAID_belgium$LINEAGE1=="B.1.1.7") # 9821 among baseline surveillance

table(GISAID_belgium$LINEAGE1)
table(GISAID_belgium$LINEAGE2)

main_lineages = names(table(GISAID_belgium$LINEAGE1))[100*table(GISAID_belgium$LINEAGE1)/sum(table(GISAID_belgium$LINEAGE1)) > 3]
main_lineages
# "B.1.1.7"  "B.1.160"  "B.1.177+" "B.1.221"  "B.1.351"  "P.1"
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.207","B.1.429", # "B.1.1.318",
         "B.1.214.2") # I added B.1.214.2 here # cut "B.1.525","B.1.526",
main_lineages = union(main_lineages, VOCs)
main_lineages = c(main_lineages,"B.1","B.1.1","B.1.258")
GISAID_belgium$LINEAGE1[!(GISAID_belgium$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_belgium$LINEAGE2[!(GISAID_belgium$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
# remove = names(table(GISAID_belgium$LINEAGE1))[table(GISAID_belgium$LINEAGE1) < 10]
# GISAID_belgium$LINEAGE1[(GISAID_belgium$LINEAGE1 %in% remove)] = "other" # minority VOCs
# GISAID_belgium$LINEAGE2[(GISAID_belgium$LINEAGE2 %in% remove)] = "other" # minority VOCs
table(GISAID_belgium$LINEAGE1)
GISAID_belgium$LINEAGE1 = factor(GISAID_belgium$LINEAGE1)
GISAID_belgium$LINEAGE1 = relevel(GISAID_belgium$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_belgium$LINEAGE1)
# [1] "B.1.1.7"   "B.1.1.318" "B.1.160"   "B.1.177+"  "B.1.214.2" "B.1.221"   "B.1.351"   "B.1.525"   "B.1.617+"  "other"     "P.1"   
levels_LINEAGE1 = c("B.1.1.7","B.1","B.1.1","B.1.221","B.1.258","B.1.160","B.1.177+","B.1.214.2",
                    "B.1.351","P.1","B.1.617+","other")
GISAID_belgium$LINEAGE1 = factor(GISAID_belgium$LINEAGE1, levels=levels_LINEAGE1)

GISAID_belgium$LINEAGE2 = factor(GISAID_belgium$LINEAGE2)
GISAID_belgium$LINEAGE2 = relevel(GISAID_belgium$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_belgium$LINEAGE2)
# "B.1.1.7"   "B.1.1.318" "B.1.160"   "B.1.177+"  "B.1.214.2" "B.1.221"   "B.1.351"   "B.1.525"   "B.1.617.1" "B.1.617.2" "other"     "P.1" 
levels_LINEAGE2 = c("B.1.1.7","B.1","B.1.1","B.1.221","B.1.258","B.1.160","B.1.177+","B.1.214.2",
                    "B.1.351","P.1","B.1.617.1","B.1.617.2","other")
GISAID_belgium$LINEAGE2 = factor(GISAID_belgium$LINEAGE2, levels=levels_LINEAGE2)

# # pango lineages since May 1st
# tbl = 100*table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01"),]$pango_lineage)/sum(table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01"),]$pango_lineage))
# tbl[order(tbl,decreasing=T)]
# 
# # pango lineages in category other since May 1st
# tbl = 100*table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01")&GISAID_belgium$LINEAGE1=="other",]$pango_lineage)/sum(table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01")&GISAID_belgium$LINEAGE1=="other",]$pango_lineage))
# tbl[order(tbl,decreasing=T)]

# clip data from last few weeks to avoid deposition biases
GISAID_belgium = GISAID_belgium[GISAID_belgium$date<=(max(GISAID_belgium$date)-14),]

# use data from Nov 1 onwards
GISAID_belgium = GISAID_belgium[GISAID_belgium$date>=as.Date("2020-11-01"),]
nrow(GISAID_belgium) # 13658
range(GISAID_belgium$date) # "2020-11-03" "2021-05-04"

write.csv(GISAID_belgium, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_ALL PARSED_BASELINE SELECTION.csv",row.names=F)


# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of Belgium
data_agbyweek1 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$LINEAGE1))
colnames(data_agbyweek1) = c("floor_date", "LINEAGE1", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_belgium) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$LINEAGE1 = factor(data_agbyweek1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

data_agbyweek2 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_belgium) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL


# aggregated by week and province for selected variant lineages
data_agbyweekregion1 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$province, GISAID_belgium$LINEAGE1))
colnames(data_agbyweekregion1) = c("floor_date", "province", "LINEAGE1", "count")
data_agbyweekregion1_sum = aggregate(count ~ floor_date + province, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$floor_date,data_agbyweekregion1$province),
                                                                  interaction(data_agbyweekregion1_sum$floor_date,data_agbyweekregion1_sum$province))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_belgium) # correct
data_agbyweekregion1$collection_date = as.Date(as.character(data_agbyweekregion1$floor_date))
data_agbyweekregion1$LINEAGE1 = factor(data_agbyweekregion1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweekregion1$province = factor(data_agbyweekregion1$province, levels=levels_PROVINCES)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]
data_agbyweekregion1$floor_date = NULL

data_agbyweekregion2 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$province, GISAID_belgium$LINEAGE2))
colnames(data_agbyweekregion2) = c("floor_date", "province", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + province, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$province),
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$province))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_belgium) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$province = factor(data_agbyweekregion2$province, levels=levels_PROVINCES)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]
data_agbyweekregion2$floor_date = NULL


# MULLER PLOT (RAW DATA)
library(scales)
n1 = length(levels(GISAID_belgium$LINEAGE1))
lineage_cols1 = hcl(h = seq(15, 320, length = n1), l = 65, c = 200)
lineage_cols1[which(levels(GISAID_belgium$LINEAGE1)=="B.1.617+")] = "magenta"
lineage_cols1[which(levels(GISAID_belgium$LINEAGE1)=="other")] = "grey75"

n2 = length(levels(GISAID_belgium$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_belgium$LINEAGE2)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(GISAID_belgium$LINEAGE2)=="B.1.617.2")] = "magenta"
lineage_cols2[which(levels(GISAID_belgium$LINEAGE2)=="other")] = "grey75"

muller_belgium_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE1)) + 
  # facet_wrap(~ province, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1, group=LINEAGE1), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_belgium$Week)))
muller_belgium_raw1

muller_belgium_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ province, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nRaw GISAID data") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_belgium$Week)))
muller_belgium_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_raw data.png"), width=8, height=5)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_raw data.pdf"), width=8, height=5)


muller_belgiumbyprovince_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) +
  facet_wrap(~ province, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nRaw GISAID data")
# +
# coord_cartesian(xlim=c(1,max(GISAID_belgium$Week)))
muller_belgiumbyprovince_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by state_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by state_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_belgium_multi = nnet::multinom(LINEAGE2 ~ province + DATE_NUM, data=GISAID_belgium, maxit=1000)
fit2_belgium_multi = nnet::multinom(LINEAGE2 ~ province + ns(DATE_NUM, df=2), data=GISAID_belgium, maxit=1000)
fit3_belgium_multi = nnet::multinom(LINEAGE2 ~ province + ns(DATE_NUM, df=3), data=GISAID_belgium, maxit=1000)
fit4_belgium_multi = nnet::multinom(LINEAGE2 ~ province * DATE_NUM, data=GISAID_belgium, maxit=1000)
fit5_belgium_multi = nnet::multinom(LINEAGE2 ~ province * ns(DATE_NUM, df=2), data=GISAID_belgium, maxit=1000)
fit6_belgium_multi = nnet::multinom(LINEAGE2 ~ province * ns(DATE_NUM, df=3), data=GISAID_belgium, maxit=1000)
BIC(fit1_belgium_multi, fit2_belgium_multi, fit3_belgium_multi, fit4_belgium_multi, fit5_belgium_multi, fit6_belgium_multi) 
# fit2_belgium_multi fits best (lowest BIC)

# equivalent fit with B.1.617.1,2&3 all recoded to B.1.617+
fit2_belgium_multi1 = nnet::multinom(LINEAGE1 ~ province + ns(DATE_NUM, df=2), data=GISAID_belgium, maxit=1000)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrbelgium = emtrends(fit2_belgium_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium = data.frame(confint(emtrbelgium, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrbelgium$contrasts)$p.value)
delta_r_belgium
#               contrast     estimate           SE df    asymp.LCL    asymp.UCL p.value
# 1         B.1 - B.1.1.7 -0.05379284 0.013930636 NA -0.081096390 -0.02648930 1.831004e-03
# 2       B.1.1 - B.1.1.7 -0.00228145 0.009427356 NA -0.020758728  0.01619583 9.996855e-01
# 3     B.1.221 - B.1.1.7 -0.10612712 0.008439138 NA -0.122667531 -0.08958672 2.109424e-15
# 4     B.1.258 - B.1.1.7 -0.06886362 0.019734355 NA -0.107542245 -0.03018500 6.705334e-03
# 5     B.1.160 - B.1.1.7 -0.09073316 0.008173490 NA -0.106752904 -0.07471341 2.609024e-14
# 6  (B.1.177+) - B.1.1.7 -0.12598536 0.009664864 NA -0.144928142 -0.10704257 1.443290e-15
# 7   B.1.214.2 - B.1.1.7 -0.02899817 0.004862009 NA -0.038527533 -0.01946881 1.895263e-07
# 8     B.1.351 - B.1.1.7 -0.05573465 0.004665561 NA -0.064878984 -0.04659032 6.661338e-15
# 9         P.1 - B.1.1.7  0.00383270 0.003808290 NA -0.003631411  0.01129681 8.919160e-01
# 10  B.1.617.1 - B.1.1.7 -0.59778114 0.244049231 NA -1.076108847 -0.11945344 1.277031e-01
# 11  B.1.617.2 - B.1.1.7  0.10972308 0.039282042 NA  0.032731696  0.18671447 5.448279e-02
# 12      other - B.1.1.7 -0.01733686 0.003577685 NA -0.024348995 -0.01032473 3.540335e-05

# avg growth advantage of B.1.617+ over B.1.1.7 :
emtrbelgium1 = emtrends(fit2_belgium_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium1 = data.frame(confint(emtrbelgium1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrbelgium1$contrasts)$p.value)
delta_r_belgium1
#                contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1         B.1 - B.1.1.7 -0.053769904 0.013925388 NA -0.081063162 -0.02647665 1.740201e-03
# 2       B.1.1 - B.1.1.7 -0.002279692 0.009427428 NA -0.020757111  0.01619773 9.995391e-01
# 3     B.1.221 - B.1.1.7 -0.106115216 0.008438748 NA -0.122654858 -0.08957557 0.000000e+00
# 4     B.1.258 - B.1.1.7 -0.068953818 0.019750919 NA -0.107664907 -0.03024273 6.287506e-03
# 5     B.1.160 - B.1.1.7 -0.090729280 0.008173688 NA -0.106749415 -0.07470915 1.598721e-14
# 6  (B.1.177+) - B.1.1.7 -0.125989836 0.009665323 NA -0.144933522 -0.10704615 0.000000e+00
# 7   B.1.214.2 - B.1.1.7 -0.028971670 0.004861419 NA -0.038499876 -0.01944346 2.058238e-07
# 8     B.1.351 - B.1.1.7 -0.055737277 0.004665266 NA -0.064881029 -0.04659352 0.000000e+00
# 9         P.1 - B.1.1.7  0.003832220 0.003808423 NA -0.003632153  0.01129659 8.793448e-01
# 10 (B.1.617+) - B.1.1.7  0.115914933 0.032795271 NA  0.051637384  0.18019248 5.441088e-03
# 11      other - B.1.1.7 -0.017340539 0.003577724 NA -0.024352750 -0.01032833 3.473688e-05


# pairwise growth rate difference (differences in growth rate per day) 
emtrbelgium_pairw = emtrends(fit2_belgium_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1           B.1.1.7 - B.1  0.053792844 0.013930636 NA  0.026489298  8.109639e-02 1.017531e-02
# 2         B.1.1.7 - B.1.1  0.002281450 0.009427356 NA -0.016195828  2.075873e-02 1.000000e+00
# 3       B.1.1.7 - B.1.221  0.106127124 0.008439138 NA  0.089586717  1.226675e-01 3.663736e-15
# 4       B.1.1.7 - B.1.258  0.068863620 0.019734355 NA  0.030184995  1.075422e-01 3.439446e-02
# 5       B.1.1.7 - B.1.160  0.090733159 0.008173490 NA  0.074713413  1.067529e-01 6.616929e-14
# 6    B.1.1.7 - (B.1.177+)  0.125985356 0.009664864 NA  0.107042570  1.449281e-01 1.887379e-15
# 7     B.1.1.7 - B.1.214.2  0.028998169 0.004862009 NA  0.019468806  3.852753e-02 1.217066e-06
# 8       B.1.1.7 - B.1.351  0.055734651 0.004665561 NA  0.046590319  6.487898e-02 1.521006e-14
# 9           B.1.1.7 - P.1 -0.003832700 0.003808290 NA -0.011296810  3.631411e-03 9.985381e-01
# 10    B.1.1.7 - B.1.617.1  0.597781143 0.244049231 NA  0.119453439  1.076109e+00 4.156712e-01
# 11    B.1.1.7 - B.1.617.2 -0.109723083 0.039282042 NA -0.186714470 -3.273170e-02 2.161919e-01
# 12        B.1.1.7 - other  0.017336861 0.003577685 NA  0.010324727  2.434899e-02 2.190358e-04
# 13            B.1 - B.1.1 -0.051511394 0.016523757 NA -0.083897362 -1.912542e-02 9.939048e-02
# 14          B.1 - B.1.221  0.052334280 0.015893691 NA  0.021183219  8.348534e-02 6.156404e-02
# 15          B.1 - B.1.258  0.015070776 0.023785768 NA -0.031548473  6.169002e-02 9.999884e-01
# 16          B.1 - B.1.160  0.036940315 0.015786917 NA  0.005998526  6.788210e-02 4.919485e-01
# 17       B.1 - (B.1.177+)  0.072192512 0.016573674 NA  0.039708707  1.046763e-01 1.636898e-03
# 18        B.1 - B.1.214.2 -0.024794675 0.014720286 NA -0.053645904  4.056555e-03 8.980146e-01
# 19          B.1 - B.1.351  0.001941807 0.014635966 NA -0.026744158  3.062777e-02 1.000000e+00
# 20              B.1 - P.1 -0.057625544 0.014436191 NA -0.085919959 -2.933113e-02 6.425052e-03
# 21        B.1 - B.1.617.1  0.543988299 0.244445676 NA  0.064883579  1.023093e+00 5.744133e-01
# 22        B.1 - B.1.617.2 -0.163515927 0.041680155 NA -0.245207530 -8.182432e-02 8.202692e-03
# 23            B.1 - other -0.036455983 0.014186062 NA -0.064260155 -8.651812e-03 3.376643e-01
# 24        B.1.1 - B.1.221  0.103845674 0.012179929 NA  0.079973452  1.277179e-01 1.089018e-12
# 25        B.1.1 - B.1.258  0.066582170 0.021184131 NA  0.025062037  1.081023e-01 9.290375e-02
# 26        B.1.1 - B.1.160  0.088451708 0.011976109 NA  0.064978966  1.119245e-01 6.663303e-10
# 27     B.1.1 - (B.1.177+)  0.123703905 0.013042493 NA  0.098141088  1.492667e-01 1.234568e-13
# 28      B.1.1 - B.1.214.2  0.026716719 0.010561951 NA  0.006015675  4.741776e-02 3.629423e-01
# 29        B.1.1 - B.1.351  0.053453201 0.010457351 NA  0.032957169  7.394923e-02 6.841996e-05
# 30            B.1.1 - P.1 -0.006114150 0.010146784 NA -0.026001482  1.377318e-02 9.999934e-01
# 31      B.1.1 - B.1.617.1  0.595499693 0.244230437 NA  0.116816833  1.074183e+00 4.232512e-01
# 32      B.1.1 - B.1.617.2 -0.112004534 0.040391955 NA -0.191171310 -3.283776e-02 2.258260e-01
# 33          B.1.1 - other  0.015055410 0.009785346 NA -0.004123515  3.423434e-02 9.448768e-01
# 34      B.1.221 - B.1.258 -0.037263504 0.020850451 NA -0.078129638  3.602630e-03 8.531296e-01
# 35      B.1.221 - B.1.160 -0.015393966 0.010992821 NA -0.036939500  6.151569e-03 9.728265e-01
# 36   B.1.221 - (B.1.177+)  0.019858231 0.012074542 NA -0.003807436  4.352390e-02 9.127179e-01
# 37    B.1.221 - B.1.214.2 -0.077128955 0.009698877 NA -0.096138405 -5.811950e-02 2.713474e-11
# 38      B.1.221 - B.1.351 -0.050392473 0.009569915 NA -0.069149162 -3.163578e-02 3.409422e-05
# 39          B.1.221 - P.1 -0.109959824 0.009273057 NA -0.128134682 -9.178497e-02 1.942890e-14
# 40    B.1.221 - B.1.617.1  0.491654019 0.244195101 NA  0.013040416  9.702676e-01 7.225767e-01
# 41    B.1.221 - B.1.617.2 -0.215850207 0.040179050 NA -0.294599699 -1.371007e-01 2.089087e-05
# 42        B.1.221 - other -0.088790264 0.008707872 NA -0.105857379 -7.172315e-02 1.199041e-13
# 43      B.1.258 - B.1.160  0.021869539 0.020782560 NA -0.018863531  6.260261e-02 9.977580e-01
# 44   B.1.258 - (B.1.177+)  0.057121736 0.021481588 NA  0.015018597  9.922487e-02 2.851740e-01
# 45    B.1.258 - B.1.214.2 -0.039865451 0.020313316 NA -0.079678819 -5.208295e-05 7.551983e-01
# 46      B.1.258 - B.1.351 -0.013128969 0.020257876 NA -0.052833676  2.657574e-02 9.999851e-01
# 47          B.1.258 - P.1 -0.072696320 0.020099136 NA -0.112089903 -3.330274e-02 2.304473e-02
# 48    B.1.258 - B.1.617.1  0.528917523 0.244845989 NA  0.049028203  1.008807e+00 6.212922e-01
# 49    B.1.258 - B.1.617.2 -0.178586703 0.043960288 NA -0.264747284 -9.242612e-02 4.972399e-03
# 50        B.1.258 - other -0.051526759 0.019774550 NA -0.090284164 -1.276935e-02 3.159758e-01
# 51   B.1.160 - (B.1.177+)  0.035252197 0.011955114 NA  0.011820605  5.868379e-02 1.516808e-01
# 52    B.1.160 - B.1.214.2 -0.061734989 0.009466137 NA -0.080288278 -4.318170e-02 7.096377e-08
# 53      B.1.160 - B.1.351 -0.034998508 0.009331012 NA -0.053286956 -1.671006e-02 1.484392e-02
# 54          B.1.160 - P.1 -0.094565858 0.009027165 NA -0.112258776 -7.687294e-02 1.075806e-13
# 55    B.1.160 - B.1.617.1  0.507047985 0.244185730 NA  0.028452749  9.856432e-01 6.800920e-01
# 56    B.1.160 - B.1.617.2 -0.200456242 0.040124119 NA -0.279098071 -1.218144e-01 1.141996e-04
# 57        B.1.160 - other -0.073396298 0.008512633 NA -0.090080753 -5.671184e-02 6.654677e-13
# 58 (B.1.177+) - B.1.214.2 -0.096987186 0.010785756 NA -0.118126879 -7.584749e-02 1.860734e-13
# 59   (B.1.177+) - B.1.351 -0.070250705 0.010667007 NA -0.091157654 -4.934376e-02 5.070585e-08
# 60       (B.1.177+) - P.1 -0.129818056 0.010406654 NA -0.150214722 -1.094214e-01 5.440093e-15
# 61 (B.1.177+) - B.1.617.1  0.471795788 0.244241210 NA -0.006908188  9.504998e-01 7.741704e-01
# 62 (B.1.177+) - B.1.617.2 -0.235708439 0.040454792 NA -0.314998373 -1.564185e-01 2.399342e-06
# 63     (B.1.177+) - other -0.108648495 0.009929131 NA -0.128109235 -8.918776e-02 7.516210e-14
# 64    B.1.214.2 - B.1.351  0.026736482 0.006664167 NA  0.013674954  3.979801e-02 5.973530e-03
# 65        B.1.214.2 - P.1 -0.032830869 0.006050997 NA -0.044690605 -2.097113e-02 1.629213e-05
# 66  B.1.214.2 - B.1.617.1  0.568782974 0.244089940 NA  0.090375482  1.047190e+00 4.988790e-01
# 67  B.1.214.2 - B.1.617.2 -0.138721253 0.039571006 NA -0.216279000 -6.116351e-02 3.272934e-02
# 68      B.1.214.2 - other -0.011661309 0.005937101 NA -0.023297812 -2.480514e-05 7.541877e-01
# 69          B.1.351 - P.1 -0.059567351 0.005965060 NA -0.071258654 -4.787605e-02 1.294520e-13
# 70    B.1.351 - B.1.617.1  0.542046492 0.244092917 NA  0.063633167  1.020460e+00 5.778363e-01
# 71    B.1.351 - B.1.617.2 -0.165457734 0.039557948 NA -0.242989888 -8.792558e-02 3.182754e-03
# 72        B.1.351 - other -0.038397790 0.005758749 NA -0.049684731 -2.711085e-02 3.291650e-08
# 73        P.1 - B.1.617.1  0.601613843 0.244076889 NA  0.123231931  1.079996e+00 4.052804e-01
# 74        P.1 - B.1.617.2 -0.105890383 0.039433205 NA -0.183178045 -2.860272e-02 2.707331e-01
# 75            P.1 - other  0.021169560 0.005161051 NA  0.011054087  3.128503e-02 4.303095e-03
# 76  B.1.617.1 - B.1.617.2 -0.707504226 0.247195848 NA -1.191999186 -2.230093e-01 1.855697e-01
# 77      B.1.617.1 - other -0.580444283 0.244074086 NA -1.058820700 -1.020679e-01 4.648976e-01
# 78      B.1.617.2 - other  0.127059944 0.039440163 NA  0.049758645  2.043612e-01 7.512879e-02

# estimated proportion of different LINEAGES in Belgium today
today # "2021-05-24"
# 62% [58%-65%] 95% CLs now estimated to be B.1.617.2 across all provinces
multinom_preds_today_avg = data.frame(emmeans(fit2_belgium_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1    B.1.1.7 8.228281e-01 5.760338e-02 NA  7.099276e-01 9.357287e-01
# 2        B.1 8.178679e-05 8.019710e-05 NA -7.539663e-05 2.389702e-04
# 3      B.1.1 7.096843e-04 4.793703e-04 NA -2.298642e-04 1.649233e-03
# 4    B.1.221 2.119297e-05 1.286224e-05 NA -4.016561e-06 4.640249e-05
# 5    B.1.258 1.216902e-05 1.903884e-05 NA -2.514641e-05 4.948445e-05
# 6    B.1.160 4.913655e-05 2.863917e-05 NA -6.995193e-06 1.052683e-04
# 7   B.1.177+ 4.331101e-06 3.012075e-06 NA -1.572458e-06 1.023466e-05
# 8  B.1.214.2 6.974436e-03 1.995107e-03 NA  3.064098e-03 1.088477e-02
# 9    B.1.351 2.845962e-03 7.714444e-04 NA  1.333959e-03 4.357965e-03
# 10       P.1 7.933595e-02 1.376814e-02 NA  5.235089e-02 1.063210e-01
# 11 B.1.617.1 1.493088e-16 2.077326e-15 NA -3.922176e-15 4.220793e-15
# 12 B.1.617.2 7.580888e-02 6.503611e-02 NA -5.165955e-02 2.032773e-01
# 13     other 1.132834e-02 2.546631e-03 NA  6.337037e-03 1.631965e-02


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-11-01"))
date.to = as.numeric(as.Date("2021-06-14")) # max(GISAID_belgium$DATE_NUM)+extrapolate

fit_belgium_multi_predsbyprovince = data.frame(emmeans(fit2_belgium_multi,
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "province"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=3)), # by=3 just to speed up things a bit
                                                  mode="prob", df=NA))
fit_belgium_multi_predsbyprovince$collection_date = as.Date(fit_belgium_multi_predsbyprovince$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_predsbyprovince$LINEAGE2 = factor(fit_belgium_multi_predsbyprovince$LINEAGE2, levels=levels_LINEAGE2)
fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=levels_PROVINCES)

fit_belgium_multi_preds = data.frame(emmeans(fit2_belgium_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=3)), # by=3 just to speed up things a bit
                                           mode="prob", df=NA))
fit_belgium_multi_preds$collection_date = as.Date(fit_belgium_multi_preds$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_preds$LINEAGE2 = factor(fit_belgium_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 

muller_belgium_mfit = ggplot(data=fit_belgium_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ province, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_belgium$DATE_NUM)+1, 
           xmax=as.Date("2021-06-14"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(multinomial fit)")
muller_belgium_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_multinom fit.png"), width=8, height=5)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_multinom fit.pdf"), width=8, height=5)


library(ggpubr)
ggarrange(muller_belgium_raw2+coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date("2021-06-14")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_belgium_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots multipanel_multinom fit.png"), width=8, height=8)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots multipanel_multinom fit.pdf"), width=8, height=8)


muller_belgiumbyprovince_mfit = ggplot(data=fit_belgium_multi_predsbyprovince,
                                  aes(x=collection_date, y=prob, group=LINEAGE2)) +
  facet_wrap(~ province, ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_belgium$DATE_NUM)+1,
           xmax=as.Date("2021-06-14"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("Share") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(multinomial fit)")
muller_belgiumbyprovince_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_multinom fit.pdf"), width=8, height=6)

ggarrange(muller_belgiumbyprovince_raw2 +
            coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date("2021-06-14")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0)))+
            ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nRaw GISAID data"),
          muller_belgiumbyprovince_mfit+ggtitle("\nMultinomial fit")+
            coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date("2021-06-14"))), nrow=2)

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province multipanel_multinom fit.png"), width=8, height=12)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province multipanel_multinom fit.pdf"), width=8, height=12)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS ####

fit_belgium_multi_preds2 = fit_belgium_multi_preds
fit_belgium_multi_preds2$LINEAGE2 = factor(fit_belgium_multi_preds2$LINEAGE2, levels=levels_LINEAGE1)
fit_belgium_multi_preds2$LINEAGE1 = fit_belgium_multi_preds2$LINEAGE2
levels(fit_belgium_multi_preds2$LINEAGE1)

# on logit scale:

fit_belgium_multi_preds2 = fit_belgium_multi_preds
ymin = 0.001
ymax = 0.998
fit_belgium_multi_preds2$asymp.LCL[fit_belgium_multi_preds2$asymp.LCL<ymin] = ymin
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$asymp.UCL<ymin] = ymin
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$asymp.UCL>ymax] = ymax
fit_belgium_multi_preds2$prob[fit_belgium_multi_preds2$prob<ymin] = ymin

plot_belgium_mfit_logit = qplot(data=fit_belgium_multi_preds2, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-14",NA)), expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date("2021-06-14")), ylim=c(0.005, 0.95), expand=c(0,0))
plot_belgium_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_belgium_mfit = qplot(data=fit_belgium_multi_preds, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share among new infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-14",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-11-01","2021-06-14")),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_belgium_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_response scale.pdf"), width=8, height=6)


# # project multinomial fit onto case data ####
# cases_belgium_byprovince = XXX # cumulative cases
# cases_belgium_byprovince$Date = as.Date(cases_belgium_byprovince$Date)
# cases_belgium_byprovince = cases_belgium_byprovince[cases_belgium_byprovince$Date >= as.Date("2020-06-01"),]
# head(cases_belgium_byprovince)
# levels_PROVINCES
# 
# # plot new cases per day by province
# ggplot(data=cases_belgium_byprovince,
#        aes(x=Date, y=newcases, 
#            group=State)) +
#   facet_wrap(~ State, scale="free", ncol=5) +
#   geom_smooth(aes(lwd=I(1), colour=State), method="loess", span=0.3, se=FALSE) +
#   # geom_line(aes(lwd=I(1), colour=State)) +
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                      limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right",
#                      axis.title.x=element_blank()) +
#   ylab("New confirmed cases per day") +
#   ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY province IN INDIA") +
#   scale_y_log10() +
#   theme(legend.position = "none") # +
# #  coord_cartesian(ylim=c(1,NA)) # +
# # coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_belgium_multi_predsbyprovince2$collection_date)-20))
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_cases per day by province.png"), width=12, height=12)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_cases per day by province.pdf"), width=12, height=12)
# 
# 
# cases_belgium_byprovince2 = cases_belgium_byprovince[cases_belgium_byprovince$State %in% levels_PROVINCES,]
# colnames(cases_belgium_byprovince2)[2]="province"
# 
# newdat = expand.grid(DATE_NUM=seq(as.numeric(min(cases_belgium_byprovince2$Date)),as.numeric(max(cases_belgium_byprovince2$Date))),
#                      province=unique(as.character(cases_belgium_byprovince2$province)))
# fit_belgium_multi_predsbyprovince = data.frame(newdat,
#                                           predict(fit5_belgium_multi, 
#                                                   newdata = newdat,
#                                                   type = "prob"), check.names=F)  
# fit_belgium_multi_predsbyprovince = gather(fit_belgium_multi_predsbyprovince, LINEAGE2, prob, all_of(levels_LINEAGE2))
# fit_belgium_multi_predsbyprovince$collection_date = as.Date(fit_belgium_multi_predsbyprovince$DATE_NUM, origin="1970-01-01")
# fit_belgium_multi_predsbyprovince$LINEAGE2 = factor(fit_belgium_multi_predsbyprovince$LINEAGE2, levels=levels_LINEAGE2)
# colnames(fit_belgium_multi_predsbyprovince)[2] = "province"
# fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
#                                                                                          "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))
# fit_belgium_multi_predsbyprovince$totnewcases = cases_belgium_byprovince2$newcases[match(interaction(fit_belgium_multi_predsbyprovince$province,fit_belgium_multi_predsbyprovince$collection_date),
#                                                                                interaction(cases_belgium_byprovince2$province,cases_belgium_byprovince2$Date))]
# fit_belgium_multi_predsbyprovince$cases = fit_belgium_multi_predsbyprovince$totnewcases*fit_belgium_multi_predsbyprovince$prob
# fit_belgium_multi_predsbyprovince$cases[fit_belgium_multi_predsbyprovince$cases==0] = NA
# fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province,
#                                             levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
#                                                      "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))
# 
# fit_belgium_multi_predsbyprovince2 = fit_belgium_multi_predsbyprovince
# fit_belgium_multi_predsbyprovince2$cases[fit_belgium_multi_predsbyprovince2$cases==0] = NA
# fit_belgium_multi_predsbyprovince2$cases[fit_belgium_multi_predsbyprovince2$cases<=1] = NA
# fit_belgium_multi_predsbyprovince2$province = factor(fit_belgium_multi_predsbyprovince2$province,
#                                              levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
#                                                       "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))
# cases_belgium_byprovince2$province = factor(cases_belgium_byprovince2$province,
#                                     levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
#                                              "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))
#
# ggplot(data=fit_belgium_multi_predsbyprovince2, 
#        aes(x=collection_date, y=cases)) + 
#   facet_wrap(~ province, scale="free", ncol=3) +
#   geom_smooth(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2), method="loess", span=0.3, se=FALSE) +
#   geom_smooth(data=cases_belgium_byprovince2, aes(x=Date, y=newcases, lwd=I(1.5)), method="loess", span=0.3, se=FALSE, colour="black") +
#   # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                      limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("New confirmed cases per day") +
#   ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN BELGIUM\n(multinomial fit)") +
#   scale_colour_manual("lineage", values=lineage_cols2) +
#   scale_y_log10() +
#   coord_cartesian(ylim=c(1,NA)) # +
# # coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_belgium_multi_predsbyprovince2$collection_date)-20))
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_confirmed cases multinomial fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_confirmed cases multinomial fit.pdf"), width=8, height=6)
# 
# ggplot(data=fit_belgium_multi_predsbyprovince2, 
#        aes(x=collection_date, y=cases, group=LINEAGE2)) + 
#   facet_wrap(~ province, scale="free", ncol=3) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
#   scale_fill_manual("", values=lineage_cols2) +
#   annotate("rect", xmin=max(GISAID_belgium$DATE_NUM)+1, 
#            xmax=as.Date("2021-06-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                      limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("New confirmed cases per day") +
#   ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN BELGIUM\n(multinomial fit)")
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_confirmed cases stacked area multinomial fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_confirmed cases stacked area multinomial fit.pdf"), width=8, height=6)
