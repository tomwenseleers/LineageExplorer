# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN BELGIUM BASED ON ANALYSIS OF GISAID PATIENT METADATA
# T. Wenseleers
# last update 31 MAY 2021

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
today = as.Date("2021-05-31")
today_num = as.numeric(today)
today # "2021-05-31"
plotdir = "VOCs_belgium"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import & parse manually downloaded GISAID patient metadata for Belgium (downloaded 24/5/2021) ####
d1 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_08_belgium_subm_jan_dec_2020.tsv", col_types = cols(.default = "c"))
d2 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_08_belgium_subm_jan_feb_2021.tsv", col_types = cols(.default = "c"))
d3 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_09_belgium_subm_mar_apr_2021.tsv", col_types = cols(.default = "c"))
d4 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_06_01_21_belgium_subm_may_2021.tsv", col_types = cols(.default = "c")) # downloaded 1/6/2021
d1 = as.data.frame(d1)
d2 = as.data.frame(d2)
d3 = as.data.frame(d3)
d4 = as.data.frame(d4)
GISAID_belgium1 = rbind(d1, d2, d3, d4)
GISAID_belgium1 = GISAID_belgium1[GISAID_belgium1[,"Host"]=="Human",]
GISAID_belgium1[,"Collection date"] = as.Date(GISAID_belgium1[,"Collection date"])
GISAID_belgium1 = GISAID_belgium1[!is.na(GISAID_belgium1[,"Collection date"]),]
nrow(GISAID_belgium1) # 23174
unique(GISAID_belgium1[,"Virus name"])
unique(GISAID_belgium1[,"Location"])
unique(GISAID_belgium1[,"Host"])
unique(GISAID_belgium1[,"Additional location information"]) # ZIP codes # sometimes has additional info on whether it was a traveller

# anonymous records without location info are assumed to be active surveillance of cluster outbreaks & included in active surveillance (cf weekly Sciensano report)
GISAID_belgium1[,"Additional location information"][grepl("nknown",GISAID_belgium1[,"Additional location information"])] = NA
sum(is.na(GISAID_belgium1[,"Additional location information"])) # 6156
unique(GISAID_belgium1[,"Additional location information"])

library(stringr)
ZIP = str_extract(GISAID_belgium1[,"Additional location information"],"[0-9]{4}") # extract ZIP codes
unique(ZIP)
sum(!is.na(ZIP)) # 16943
sum(is.na(ZIP)) # 6231 - anonymous records, these are assumed to be active surveillance of cluster outbreaks

muni = read.csv(".//data//GISAID//Belgium//mapping_municips.csv") # post codes & cities & provinces
province = muni$province_label_nl[match(ZIP, muni$postcode)] # province
city = muni$municipality_label_nl[match(ZIP, muni$postcode)] # city
sum(is.na(city)) # 6261
sum(is.na(ZIP))  # 6231
wrongZIP = which(is.na(city)&!is.na(ZIP))
ZIP[wrongZIP] = NA

anonym = grepl("Belgium$",GISAID_belgium1[,"Additional location information"])|
  is.na(GISAID_belgium1[,"Additional location information"])|
  grepl("CHU|infected|travel|Travel",GISAID_belgium1[,"Additional location information"]) |
  is.na(ZIP)
sum(anonym) # 6265
unique(GISAID_belgium1[!is.na(GISAID_belgium1[,"Sampling strategy"]),"Sampling strategy"])
traveller = grepl("travel|Travel|infected",GISAID_belgium1[,"Additional location information"])|grepl("travel|Travel",GISAID_belgium1[,"Sampling strategy"])|
  grepl("travel|Travel|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(traveller) # 143
Sdropout = grepl("S-gene dropout",GISAID_belgium1[,"Sampling strategy"])|
  grepl("S-gene dropout",GISAID_belgium1[,"Additional host information"])  
sum(Sdropout) # 61 - should be much higher though...
actsurv = anonym|traveller|grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster",GISAID_belgium1[,"Sampling strategy"])|
  grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster|HR contact|Outbreak|Nurse|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(actsurv) # 6954
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")]) # 4292 # In Sciensano weekly report of 21st of May 2021 this is 4710
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")&GISAID_belgium1[,"Collection date"]<=as.Date("2021-01-31")]) # 1172 # In Sciensano weekly report of 21st of May 2021 this is 1975

reinfection = grepl("Breakthrough",GISAID_belgium1[,"Sampling strategy"])
sum(reinfection) # 7
infectionpostvaccination = grepl("Vaccination|vaccination",GISAID_belgium1[,"Sampling strategy"])
sum(infectionpostvaccination) # 74

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
sum(bassurv) # 16220
purpose_of_sequencing = rep("active_surveillance", nrow(GISAID_belgium1))
purpose_of_sequencing[bassurv] = "baseline_surveillance"

accessionIDs_bassurv = GISAID_belgium1[,"Accession ID"][bassurv]

GISAID_belgium1[,"Gender"][grepl("F|V",GISAID_belgium1[,"Gender"])] = "F"
GISAID_belgium1[,"Gender"][grepl("M|m",GISAID_belgium1[,"Gender"])] = "M"
GISAID_belgium1[,"Gender"][grepl("unknown",GISAID_belgium1[,"Gender"])] = NA
unique(GISAID_belgium1[,"Gender"]) # "F" "M" NA 

sum(!is.na(GISAID_belgium1[,"Gender"] )) # 20621 with gender info

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
sum(!is.na(GISAID_belgium1[,"Patient age"] )) # 20633 with age info

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



# get extra fields "genbank_accession", "Nextstrain_clade", "originating_lab", "submitting_lab", "authors", "url", "title", "paper_url"
# from genomic epidemiology GISAID metadata
# (check with Emmanuel André & Lize Cuypers which labs were doing active surveillance vs baseline surveillance)
# we use genomic epidemiology GISAID data file version metadata_2021-05-31_10-15.tsv.gz with data for all countries :
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-05-31_10-15.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID_genepi_belgium = GISAID[GISAID$country=="Belgium",]
nrow(GISAID_genepi_belgium) # 22695

# add (genbank_accession,) Nextstrain_clade, originating_lab, submitting_lab, authors, url, title & paper_url to GISAID_belgium1
# GISAID_belgium1$genbank_accession = GISAID_genepi_belgium$genbank_accession[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
# GISAID_belgium1$genbank_accession[GISAID_belgium1$genbank_accession=="?"] = NA
# sum(!is.na(GISAID_belgium1$genbank_accession)) # 0 - all NA unfortunately
GISAID_belgium1$Nextstrain_clade = GISAID_genepi_belgium$Nextstrain_clade[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
GISAID_belgium1$originating_lab = GISAID_genepi_belgium$originating_lab[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
GISAID_belgium1$submitting_lab = GISAID_genepi_belgium$submitting_lab[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
GISAID_belgium1$authors = GISAID_genepi_belgium$authors[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
GISAID_belgium1$url = GISAID_genepi_belgium$url[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
GISAID_belgium1$title = GISAID_genepi_belgium$title[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
GISAID_belgium1$paper_url = GISAID_genepi_belgium$paper_url[match(GISAID_belgium1$gisaid_epi_isl,GISAID_genepi_belgium$gisaid_epi_isl)]
GISAID_belgium1$paper_url[GISAID_belgium1$paper_url=="?"] = NA


# write unique originating & submitting labs to file
uniqueoriginatinglabs = unique(interaction(GISAID_belgium1$province,GISAID_belgium1$originating_lab))
uniquesubmittinglabs = unique(interaction(GISAID_belgium1$province,GISAID_belgium1$submitting_lab))
uniquelabs = unique(c(GISAID_belgium1$originating_lab,GISAID_belgium1$submitting_lab))
write.csv(uniqueoriginatinglabs,".//data//GISAID//Belgium//unique_originating_labs.csv")
write.csv(uniquesubmittinglabs,".//data//GISAID//Belgium//unique_submitting_labs.csv")
write.csv(uniquelabs,".//data//GISAID//Belgium//unique_labs.csv")
write.csv(data.frame(table(GISAID_belgium1$originating_lab)), ".//data//GISAID//Belgium//counts_originating_labs.csv")
write.csv(data.frame(table(GISAID_belgium1$submitting_lab)), ".//data//GISAID//Belgium//counts_submitting_labs.csv")

# fix lab names & also include under baseline surveillance only labs that actually did baseline surveillance
labnames = read.csv(".//data//GISAID//Belgium//lab_names.csv")
GISAID_belgium1$originating_lab_cleaned = labnames$lab_name[match(GISAID_belgium1$originating_lab,labnames$lab)]
GISAID_belgium1$submitting_lab_cleaned = labnames$lab_name[match(GISAID_belgium1$submitting_lab,labnames$lab)]
GISAID_belgium1$originating_lab_does_baseline_surveillance = labnames$do_baseline_surveillance[match(GISAID_belgium1$originating_lab,labnames$lab)]
GISAID_belgium1$purpose_of_sequencing[GISAID_belgium1$originating_lab_does_baseline_surveillance=="no"] = "active_surveillance"

# write parsed & cleaned up file to csv
write.csv(GISAID_belgium1, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_06_01_ALL PARSED.csv",row.names=F)




# ANALYSIS OF VOC LINEAGE FREQUENCIES IN BELGIUM USING MULTINOMIAL FITS ####

# GISAID_belgium = GISAID_sel[GISAID_sel$country=="Belgium",]
GISAID_belgium = GISAID_belgium1[GISAID_belgium1$purpose_of_sequencing=="baseline_surveillance",]
nrow(GISAID_belgium) # 15092
sum(!is.na(GISAID_belgium$age)) # 13971 with age info

# plot age distribution of B.1.617.2 & B.1.1.7 cases in Belgium
df = GISAID_belgium[GISAID_belgium$date>=as.Date("2021-04-01")&GISAID_belgium$pango_lineage %in% c("B.1.1.7","B.1.617.2"),]
qplot(data=df, binwidth=10.0, x=age, y = ..density.., geom="histogram", fill=pango_lineage) + facet_wrap(~pango_lineage, ncol=1)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_B1672_B117_age distribution.png"), width=8, height=5)


# unique(GISAID_belgium$originating_lab) # 133 labs # maybe some labs did more active surveillance & could be excluded?
# unique(GISAID_belgium$submitting_lab) # 41 labs
# sum(GISAID_belgium$country_exposure!="Belgium") # only 23 are indicated as travellers

# # we remove travellers
# GISAID_belgium = GISAID_belgium[GISAID_belgium$country_exposure=="Belgium",]
# nrow(GISAID_belgium) # 19283
# sum(GISAID_belgium$division_exposure==GISAID_belgium$division) # 19283

unique(GISAID_belgium$province) # 
unique(GISAID_belgium$province[GISAID_belgium$LINEAGE1=="B.1.617+"]) # "Belgium"         "Hasselt"         "Brussels"        "Brugge"          "Gent"            "Liège"           "Halle-Vilvoorde" "Antwerpen"       "Mechelen"
sum(GISAID_belgium$LINEAGE1=="B.1.617+") # 27 among baseline surveillance
unique(GISAID_belgium$province[GISAID_belgium$LINEAGE1=="B.1.1.7"])
sum(GISAID_belgium$LINEAGE1=="B.1.1.7") # 10630 among baseline surveillance

table(GISAID_belgium$LINEAGE1)
table(GISAID_belgium$LINEAGE2)

main_lineages = names(table(GISAID_belgium$LINEAGE1))[100*table(GISAID_belgium$LINEAGE1)/sum(table(GISAID_belgium$LINEAGE1)) > 3]
main_lineages
# "B.1.1.7"  "B.1.160"  "B.1.177+" "B.1.221"  "B.1.351"  "P.1"
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.620","B.1.1.7","B.1.351","P.1","B.1.1.207","B.1.429", # "B.1.1.318",
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
# [1] "B.1.1.7"   "B.1"       "B.1.1"     "B.1.160"   "B.1.177+"  "B.1.214.2" "B.1.221"   "B.1.258"   "B.1.351"   "B.1.617+"  "B.1.620"   "other"     "P.1" 
levels_LINEAGE1 = c("B.1.1.7","B.1","B.1.1","B.1.221","B.1.258","B.1.160","B.1.177+","B.1.214.2",
                    "B.1.351","P.1","B.1.620","B.1.617+","other")
GISAID_belgium$LINEAGE1 = factor(GISAID_belgium$LINEAGE1, levels=levels_LINEAGE1)

GISAID_belgium$LINEAGE2 = factor(GISAID_belgium$LINEAGE2)
GISAID_belgium$LINEAGE2 = relevel(GISAID_belgium$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_belgium$LINEAGE2)
# "B.1.1.7"   "B.1"       "B.1.1"     "B.1.160"   "B.1.177+"  "B.1.214.2" "B.1.221"   "B.1.258"   "B.1.351"   "B.1.617.1" "B.1.617.2" "B.1.620"   "other"    
# "P.1"
levels_LINEAGE2 = c("B.1.1.7","B.1","B.1.1","B.1.221","B.1.258","B.1.160","B.1.177+","B.1.214.2",
                    "B.1.351","B.1.620","P.1","B.1.617.1","B.1.617.2","other")
GISAID_belgium$LINEAGE2 = factor(GISAID_belgium$LINEAGE2, levels=levels_LINEAGE2)

# # pango lineages since May 1st
# tbl = 100*table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01"),]$pango_lineage)/sum(table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01"),]$pango_lineage))
# tbl[order(tbl,decreasing=T)]
# 
# # pango lineages in category other since May 1st
# tbl = 100*table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01")&GISAID_belgium$LINEAGE1=="other",]$pango_lineage)/sum(table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01")&GISAID_belgium$LINEAGE1=="other",]$pango_lineage))
# tbl[order(tbl,decreasing=T)]

# clip data from last few weeks to avoid deposition biases
# GISAID_belgium = GISAID_belgium[GISAID_belgium$date<=(max(GISAID_belgium$date)-14),]

# use data from Nov 1 onwards
GISAID_belgium = GISAID_belgium[GISAID_belgium$date>=as.Date("2020-11-01"),]
nrow(GISAID_belgium) # 14807
range(GISAID_belgium$date) # "2020-11-03" "2021-05-26"

# write.csv(GISAID_belgium, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_ALL PARSED_BASELINE SELECTION.csv",row.names=F)


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

# write.csv(data_agbyweek2, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_BASELINE SELECTION_aggregated counts by week.csv",row.names=F)


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
# 1         B.1 - B.1.1.7 -4.377182e-02 0.013022512 NA -6.929547e-02 -1.824816e-02 1.080461e-02
# 2       B.1.1 - B.1.1.7  3.562279e-02 0.007277365 NA  2.135941e-02  4.988616e-02 2.902215e-05
# 3     B.1.221 - B.1.1.7 -1.106480e-01 0.010234930 NA -1.307081e-01 -9.058794e-02 0.000000e+00
# 4     B.1.258 - B.1.1.7 -4.607312e-02 0.023332610 NA -9.180420e-02 -3.420451e-04 3.440987e-01
# 5     B.1.160 - B.1.1.7 -1.082899e-01 0.011354366 NA -1.305441e-01 -8.603576e-02 1.132427e-14
# 6  (B.1.177+) - B.1.1.7 -1.455738e-01 0.012904753 NA -1.708667e-01 -1.202810e-01 0.000000e+00
# 7   B.1.214.2 - B.1.1.7 -3.718344e-02 0.004867989 NA -4.672452e-02 -2.764236e-02 2.023204e-11
# 8     B.1.351 - B.1.1.7 -6.012045e-02 0.005214107 NA -7.033991e-02 -4.990099e-02 0.000000e+00
# 9     B.1.620 - B.1.1.7 -8.185889e-02 0.044970812 NA -1.700001e-01  6.282286e-03 4.373267e-01
# 10        P.1 - B.1.1.7 -4.289362e-04 0.003258898 NA -6.816259e-03  5.958387e-03 9.999867e-01
# 11  B.1.617.1 - B.1.1.7 -8.700874e+02 0.001255733 NA -8.700899e+02 -8.700849e+02 0.000000e+00
# 12  B.1.617.2 - B.1.1.7  1.194345e-02 0.029225663 NA -4.533780e-02  6.922470e-02 9.976922e-01
# 13      other - B.1.1.7  2.586467e-03 0.002775713 NA -2.853831e-03  8.026764e-03 9.270410e-01


# avg growth advantage of B.1.617+ over B.1.1.7 :
emtrbelgium1 = emtrends(fit2_belgium_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium1 = data.frame(confint(emtrbelgium1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrbelgium1$contrasts)$p.value)
delta_r_belgium1
#                contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1         B.1 - B.1.1.7 -0.0437806377 0.013021617 NA -0.069302537 -0.0182587380 1.017726e-02
# 2       B.1.1 - B.1.1.7  0.0356316113 0.007277606 NA  0.021367766  0.0498954568 2.843503e-05
# 3     B.1.221 - B.1.1.7 -0.1106484536 0.010235159 NA -0.130708997 -0.0905879100 3.474998e-14
# 4     B.1.258 - B.1.1.7 -0.0460365261 0.023321329 NA -0.091745491 -0.0003275608 3.292957e-01
# 5     B.1.160 - B.1.1.7 -0.1082878147 0.011354239 NA -0.130541713 -0.0860339161 4.962697e-14
# 6  (B.1.177+) - B.1.1.7 -0.1455684170 0.012904912 NA -0.170861579 -0.1202752549 2.020606e-14
# 7   B.1.214.2 - B.1.1.7 -0.0371830618 0.004868581 NA -0.046725305 -0.0276408183 2.507772e-11
# 8     B.1.351 - B.1.1.7 -0.0601310708 0.005214311 NA -0.070350932 -0.0499112091 1.465494e-14
# 9         P.1 - B.1.1.7 -0.0004229302 0.003259232 NA -0.006810907  0.0059650471 9.999796e-01
# 10    B.1.620 - B.1.1.7 -0.0818412212 0.044958082 NA -0.169957444  0.0062750013 4.199212e-01
# 11 (B.1.617+) - B.1.1.7  0.0324161878 0.023018199 NA -0.012698653  0.0775310290 6.859237e-01
# 12      other - B.1.1.7  0.0025843956 0.002775843 NA -0.002856158  0.0080249489 9.180671e-01


# pairwise growth rate difference (differences in growth rate per day) 
emtrbelgium_pairw = emtrends(fit2_belgium_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1           B.1.1.7 - B.1  4.377182e-02 0.013022512 NA  1.824816e-02  6.929547e-02 5.656862e-02
# 2         B.1.1.7 - B.1.1 -3.562279e-02 0.007277365 NA -4.988616e-02 -2.135941e-02 1.935513e-04
# 3       B.1.1.7 - B.1.221  1.106480e-01 0.010234930 NA  9.058794e-02  1.307081e-01 5.229150e-14
# 4       B.1.1.7 - B.1.258  4.607312e-02 0.023332610 NA  3.420451e-04  9.180420e-02 7.810524e-01
# 5       B.1.1.7 - B.1.160  1.082899e-01 0.011354366 NA  8.603576e-02  1.305441e-01 9.525714e-14
# 6    B.1.1.7 - (B.1.177+)  1.455738e-01 0.012904753 NA  1.202810e-01  1.708667e-01 9.992007e-15
# 7     B.1.1.7 - B.1.214.2  3.718344e-02 0.004867989 NA  2.764236e-02  4.672452e-02 1.414805e-10
# 8       B.1.1.7 - B.1.351  6.012045e-02 0.005214107 NA  4.990099e-02  7.033991e-02 0.000000e+00
# 9       B.1.1.7 - B.1.620  8.185889e-02 0.044970812 NA -6.282286e-03  1.700001e-01 8.639409e-01
# 10          B.1.1.7 - P.1  4.289362e-04 0.003258898 NA -5.958387e-03  6.816259e-03 1.000000e+00
# 11    B.1.1.7 - B.1.617.1  8.700874e+02 0.001255733 NA  8.700849e+02  8.700899e+02 0.000000e+00
# 12    B.1.1.7 - B.1.617.2 -1.194345e-02 0.029225663 NA -6.922470e-02  4.533780e-02 1.000000e+00
# 13        B.1.1.7 - other -2.586467e-03 0.002775713 NA -8.026764e-03  2.853831e-03 9.996471e-01
# 14            B.1 - B.1.1 -7.939460e-02 0.014674663 NA -1.081564e-01 -5.063279e-02 1.863983e-05
# 15          B.1 - B.1.221  6.687622e-02 0.016098261 NA  3.532421e-02  9.842823e-02 3.953387e-03
# 16          B.1 - B.1.258  2.301305e-03 0.026411024 NA -4.946335e-02  5.406596e-02 1.000000e+00
# 17          B.1 - B.1.160  6.451809e-02 0.016854811 NA  3.148327e-02  9.755291e-02 1.280873e-02
# 18       B.1 - (B.1.177+)  1.018020e-01 0.017911767 NA  6.669560e-02  1.369084e-01 5.029474e-06
# 19        B.1 - B.1.214.2 -6.588376e-03 0.013862546 NA -3.375847e-02  2.058172e-02 9.999999e-01
# 20          B.1 - B.1.351  1.634863e-02 0.013951380 NA -1.099557e-02  4.369284e-02 9.962402e-01
# 21          B.1 - B.1.620  3.808707e-02 0.046819114 NA -5.367671e-02  1.298508e-01 9.999218e-01
# 22              B.1 - P.1 -4.334288e-02 0.013414510 NA -6.963484e-02 -1.705092e-02 8.175594e-02
# 23        B.1 - B.1.617.1  8.700436e+02 0.013082996 NA  8.700180e+02  8.700693e+02 0.000000e+00
# 24        B.1 - B.1.617.2 -5.571527e-02 0.031995042 NA -1.184244e-01  6.993862e-03 8.978699e-01
# 25            B.1 - other -4.635828e-02 0.013153820 NA -7.213930e-02 -2.057727e-02 3.462706e-02
# 26        B.1.1 - B.1.221  1.462708e-01 0.012045859 NA  1.226614e-01  1.698803e-01 0.000000e+00
# 27        B.1.1 - B.1.258  8.169591e-02 0.023822809 NA  3.500406e-02  1.283878e-01 4.626866e-02
# 28        B.1.1 - B.1.160  1.439127e-01 0.012968216 NA  1.184955e-01  1.693299e-01 2.731149e-14
# 29     B.1.1 - (B.1.177+)  1.811966e-01 0.014349008 NA  1.530731e-01  2.093202e-01 0.000000e+00
# 30      B.1.1 - B.1.214.2  7.280623e-02 0.008701746 NA  5.575112e-02  8.986133e-02 2.034484e-12
# 31        B.1.1 - B.1.351  9.574324e-02 0.008871193 NA  7.835602e-02  1.131305e-01 5.262457e-14
# 32        B.1.1 - B.1.620  1.174817e-01 0.045554205 NA  2.819707e-02  2.067663e-01 3.619539e-01
# 33            B.1.1 - P.1  3.605172e-02 0.007940478 NA  2.048867e-02  5.161477e-02 8.679128e-04
# 34      B.1.1 - B.1.617.1  8.701230e+02 0.007384666 NA  8.701086e+02  8.701375e+02 0.000000e+00
# 35      B.1.1 - B.1.617.2  2.367933e-02 0.030098340 NA -3.531233e-02  8.267099e-02 9.999466e-01
# 36          B.1.1 - other  3.303632e-02 0.007529971 NA  1.827785e-02  4.779479e-02 1.607042e-03
# 37      B.1.221 - B.1.258 -6.457492e-02 0.024753003 NA -1.130899e-01 -1.605992e-02 3.429806e-01
# 38      B.1.221 - B.1.160 -2.358129e-03 0.014231310 NA -3.025098e-02  2.553473e-02 1.000000e+00
# 39   B.1.221 - (B.1.177+)  3.492580e-02 0.015455331 NA  4.633908e-03  6.521769e-02 5.867057e-01
# 40    B.1.221 - B.1.214.2 -7.346460e-02 0.011296257 NA -9.560485e-02 -5.132434e-02 7.686933e-08
# 41      B.1.221 - B.1.351 -5.052759e-02 0.011377714 NA -7.282750e-02 -2.822768e-02 1.297753e-03
# 42      B.1.221 - B.1.620 -2.878915e-02 0.046128076 NA -1.191985e-01  6.162022e-02 9.999965e-01
# 43          B.1.221 - P.1 -1.102191e-01 0.010751767 NA -1.312922e-01 -8.914603e-02 8.504308e-14
# 44    B.1.221 - B.1.617.1  8.699768e+02 0.010312187 NA  8.699565e+02  8.699970e+02 0.000000e+00
# 45    B.1.221 - B.1.617.2 -1.225915e-01 0.030961361 NA -1.832746e-01 -6.190834e-02 8.073488e-03
# 46        B.1.221 - other -1.132345e-01 0.010187773 NA -1.332022e-01 -9.326684e-02 2.819966e-14
# 47      B.1.258 - B.1.160  6.221679e-02 0.025285252 NA  1.265861e-02  1.117750e-01 4.417097e-01
# 48   B.1.258 - (B.1.177+)  9.950072e-02 0.026106016 NA  4.833387e-02  1.506676e-01 1.355311e-02
# 49    B.1.258 - B.1.214.2 -8.889681e-03 0.023822922 NA -5.558175e-02  3.780239e-02 1.000000e+00
# 50      B.1.258 - B.1.351  1.404733e-02 0.023875085 NA -3.274698e-02  6.084164e-02 9.999983e-01
# 51      B.1.258 - B.1.620  3.578577e-02 0.050666259 NA -6.351828e-02  1.350898e-01 9.999847e-01
# 52          B.1.258 - P.1 -4.564418e-02 0.023556572 NA -9.181422e-02  5.258489e-04 8.027680e-01
# 53    B.1.258 - B.1.617.1  8.700413e+02 0.023366617 NA  8.699955e+02  8.700871e+02 0.000000e+00
# 54    B.1.258 - B.1.617.2 -5.801657e-02 0.037392032 NA -1.313036e-01  1.527046e-02 9.556076e-01
# 55        B.1.258 - other -4.865959e-02 0.023297020 NA -9.432091e-02 -2.998266e-03 7.079622e-01
# 56   B.1.160 - (B.1.177+)  3.728393e-02 0.016250662 NA  5.433218e-03  6.913464e-02 5.615043e-01
# 57    B.1.160 - B.1.214.2 -7.110647e-02 0.012316698 NA -9.524675e-02 -4.696618e-02 3.240818e-06
# 58      B.1.160 - B.1.351 -4.816946e-02 0.012388026 NA -7.244954e-02 -2.388937e-02 1.038237e-02
# 59      B.1.160 - B.1.620 -2.643102e-02 0.046386106 NA -1.173461e-01  6.448407e-02 9.999988e-01
# 60          B.1.160 - P.1 -1.078610e-01 0.011823234 NA -1.310341e-01 -8.468786e-02 1.110223e-13
# 61    B.1.160 - B.1.617.1  8.699791e+02 0.011423726 NA  8.699567e+02  8.700015e+02 0.000000e+00
# 62    B.1.160 - B.1.617.2 -1.202334e-01 0.031349864 NA -1.816780e-01 -5.878876e-02 1.248918e-02
# 63        B.1.160 - other -1.108764e-01 0.011360521 NA -1.331426e-01 -8.861016e-02 9.303669e-14
# 64 (B.1.177+) - B.1.214.2 -1.083904e-01 0.013760867 NA -1.353612e-01 -8.141960e-02 3.549172e-11
# 65   (B.1.177+) - B.1.351 -8.545339e-02 0.013824164 NA -1.125483e-01 -5.835852e-02 4.139551e-07
# 66   (B.1.177+) - B.1.620 -6.371495e-02 0.046804679 NA -1.554504e-01  2.802053e-02 9.849325e-01
# 67       (B.1.177+) - P.1 -1.451449e-01 0.013322548 NA -1.712566e-01 -1.190332e-01 4.218847e-14
# 68 (B.1.177+) - B.1.617.1  8.699418e+02 0.012966034 NA  8.699164e+02  8.699672e+02 0.000000e+00
# 69 (B.1.177+) - B.1.617.2 -1.575173e-01 0.031943086 NA -2.201246e-01 -9.490999e-02 1.652123e-04
# 70     (B.1.177+) - other -1.481603e-01 0.012902163 NA -1.734481e-01 -1.228725e-01 0.000000e+00
# 71    B.1.214.2 - B.1.351  2.293701e-02 0.007064298 NA  9.091240e-03  3.678278e-02 7.826286e-02
# 72    B.1.214.2 - B.1.620  4.467545e-02 0.045227105 NA -4.396805e-02  1.333189e-01 9.993402e-01
# 73        B.1.214.2 - P.1 -3.675450e-02 0.005768675 NA -4.806090e-02 -2.544811e-02 1.542729e-07
# 74  B.1.214.2 - B.1.617.1  8.700502e+02 0.005023681 NA  8.700404e+02  8.700601e+02 0.000000e+00
# 75  B.1.214.2 - B.1.617.2 -4.912689e-02 0.029627233 NA -1.071952e-01  8.941417e-03 9.272149e-01
# 76      B.1.214.2 - other -3.976991e-02 0.005507321 NA -5.056406e-02 -2.897576e-02 1.519776e-09
# 77      B.1.351 - B.1.620  2.173844e-02 0.045266075 NA -6.698144e-02  1.104583e-01 9.999999e-01
# 78          B.1.351 - P.1 -5.969151e-02 0.006110951 NA -7.166876e-02 -4.771427e-02 9.459100e-14
# 79    B.1.351 - B.1.617.1  8.700273e+02 0.005362831 NA  8.700168e+02  8.700378e+02 0.000000e+00
# 80    B.1.351 - B.1.617.2 -7.206390e-02 0.029688859 NA -1.302530e-01 -1.387481e-02 4.651908e-01
# 81        B.1.351 - other -6.270692e-02 0.005780114 NA -7.403573e-02 -5.137810e-02 4.951595e-14
# 82          B.1.620 - P.1 -8.142995e-02 0.045076019 NA -1.697773e-01  6.917423e-03 8.702881e-01
# 83    B.1.620 - B.1.617.1  8.700055e+02 0.044987579 NA  8.699174e+02  8.700937e+02 0.000000e+00
# 84    B.1.620 - B.1.617.2 -9.380234e-02 0.053643227 NA -1.989411e-01  1.133645e-02 8.949963e-01
# 85        B.1.620 - other -8.444535e-02 0.045052820 NA -1.727473e-01  3.856552e-03 8.372805e-01
# 86        P.1 - B.1.617.1  8.700870e+02 0.003487909 NA  8.700801e+02  8.700938e+02 0.000000e+00
# 87        P.1 - B.1.617.2 -1.237239e-02 0.029374543 NA -6.994543e-02  4.520066e-02 1.000000e+00
# 88            P.1 - other -3.015403e-03 0.004200716 NA -1.124865e-02  5.217849e-03 9.999815e-01
# 89  B.1.617.1 - B.1.617.2 -8.700993e+02 0.029252681 NA -8.701567e+02 -8.700420e+02 0.000000e+00
# 90      B.1.617.1 - other -8.700900e+02 0.003045595 NA -8.700960e+02 -8.700840e+02 0.000000e+00
# 91      B.1.617.2 - other  9.356986e-03 0.029347232 NA -4.816253e-02  6.687650e-02 1.000000e+00

# estimated proportion of different LINEAGES in Belgium today
today # "2021-05-31"
# 62% [58%-65%] 95% CLs now estimated to be B.1.617.2 across all provinces
multinom_preds_today_avg = data.frame(emmeans(fit2_belgium_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1   B.1.1.7 2.157249e-01 1.037562e-02 NA  1.953891e-01 2.360608e-01
# 2     other 2.031711e-02 1.928644e-03 NA  1.653704e-02 2.409718e-02
# 3  B.1.177+ 6.429385e-07 3.788120e-07 NA -9.951926e-08 1.385396e-06
# 4   B.1.351 1.992858e-03 3.846646e-04 NA  1.238929e-03 2.746787e-03
# 5   B.1.525 8.276596e-04 2.361210e-04 NA  3.648711e-04 1.290448e-03
# 6 B.1.617.1 3.529633e-04 1.420207e-04 NA  7.460772e-05 6.313188e-04
# 7 B.1.617.2 7.596203e-01 1.137947e-02 NA  7.373170e-01 7.819237e-01
# 8       P.1 1.163470e-03 4.123223e-04 NA  3.553333e-04 1.971607e-03


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-11-01"))
date.to = as.numeric(as.Date("2021-06-01")) # max(GISAID_belgium$DATE_NUM)+extrapolate

fit_belgium_multi_predsbyprovince = data.frame(emmeans(fit2_belgium_multi,
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "province"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                                  mode="prob", df=NA))
fit_belgium_multi_predsbyprovince$collection_date = as.Date(fit_belgium_multi_predsbyprovince$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_predsbyprovince$LINEAGE2 = factor(fit_belgium_multi_predsbyprovince$LINEAGE2, levels=levels_LINEAGE2)
fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=levels_PROVINCES)

fit_belgium_multi_preds = data.frame(emmeans(fit2_belgium_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
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
