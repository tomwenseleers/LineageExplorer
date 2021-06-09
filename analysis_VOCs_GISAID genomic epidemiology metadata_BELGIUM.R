# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN BELGIUM BASED ON ANALYSIS OF GISAID PATIENT METADATA
# T. Wenseleers
# last update 3 JUNE 2021

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
today = as.Date("2021-06-03")
today_num = as.numeric(today)
today # "2021-06-03"
plotdir = "VOCs_belgium"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import & parse manually downloaded GISAID patient metadata for Belgium (downloaded 24/5/2021) ####
d1 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_08_belgium_subm_jan_dec_2020.tsv", col_types = cols(.default = "c"))
d2 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_08_belgium_subm_jan_feb_2021.tsv", col_types = cols(.default = "c"))
d3 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_09_belgium_subm_mar_apr_2021.tsv", col_types = cols(.default = "c"))
d4 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_06_03_08_belgium_subm_may_2021.tsv", col_types = cols(.default = "c")) # downloaded 3/6/2021
d1 = as.data.frame(d1)
d2 = as.data.frame(d2)
d3 = as.data.frame(d3)
d4 = as.data.frame(d4)
GISAID_belgium1 = rbind(d1, d2, d3, d4)
GISAID_belgium1 = GISAID_belgium1[GISAID_belgium1[,"Host"]=="Human",]
GISAID_belgium1[,"Collection date"] = as.Date(GISAID_belgium1[,"Collection date"])
GISAID_belgium1 = GISAID_belgium1[!is.na(GISAID_belgium1[,"Collection date"]),]
nrow(GISAID_belgium1) # 23497
unique(GISAID_belgium1[,"Virus name"])
unique(GISAID_belgium1[,"Location"])
unique(GISAID_belgium1[,"Host"])
unique(GISAID_belgium1[,"Additional location information"]) # ZIP codes # sometimes has additional info on whether it was a traveller

# anonymous records without location info are assumed to be active surveillance of cluster outbreaks & included in active surveillance (cf weekly Sciensano report)
GISAID_belgium1[,"Additional location information"][grepl("nknown",GISAID_belgium1[,"Additional location information"])] = NA
sum(is.na(GISAID_belgium1[,"Additional location information"])) # 6193
unique(GISAID_belgium1[,"Additional location information"])

library(stringr)
ZIP = str_extract(GISAID_belgium1[,"Additional location information"],"[0-9]{4}") # extract ZIP codes
unique(ZIP)
sum(!is.na(ZIP)) # 17229
sum(is.na(ZIP)) # 6268 - anonymous records, these are assumed to be active surveillance of cluster outbreaks

muni = read.csv(".//data//GISAID//Belgium//mapping_municips.csv") # post codes & cities & provinces
province = muni$province_label_nl[match(ZIP, muni$postcode)] # province
city = muni$municipality_label_nl[match(ZIP, muni$postcode)] # city
sum(is.na(city)) # 6298
sum(is.na(ZIP))  # 6268
wrongZIP = which(is.na(city)&!is.na(ZIP))
ZIP[wrongZIP] = NA

anonym = grepl("Belgium$",GISAID_belgium1[,"Additional location information"])|
  is.na(GISAID_belgium1[,"Additional location information"])|
  grepl("CHU|infected|travel|Travel",GISAID_belgium1[,"Additional location information"]) |
  is.na(ZIP)
sum(anonym) # 6302
unique(GISAID_belgium1[!is.na(GISAID_belgium1[,"Sampling strategy"]),"Sampling strategy"])
traveller = grepl("travel|Travel|infected",GISAID_belgium1[,"Additional location information"])|grepl("travel|Travel",GISAID_belgium1[,"Sampling strategy"])|
  grepl("travel|Travel|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(traveller) # 143
Sdropout = grepl("S-gene dropout",GISAID_belgium1[,"Sampling strategy"])|
  grepl("S-gene dropout",GISAID_belgium1[,"Additional host information"])  
sum(Sdropout) # 61 - should be much higher though...
actsurv = anonym|traveller|grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster",GISAID_belgium1[,"Sampling strategy"])|
  grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster|HR contact|Outbreak|Nurse|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(actsurv) # 6995
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")]) # 4333 # In Sciensano weekly report of 21st of May 2021 this is 4710
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")&GISAID_belgium1[,"Collection date"]<=as.Date("2021-01-31")]) # 1172 # In Sciensano weekly report of 21st of May 2021 this is 1975

reinfection = grepl("Breakthrough",GISAID_belgium1[,"Sampling strategy"])
sum(reinfection) # 14
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
sum(bassurv) # 16502
purpose_of_sequencing = rep("active_surveillance", nrow(GISAID_belgium1))
purpose_of_sequencing[bassurv] = "baseline_surveillance"

accessionIDs_bassurv = GISAID_belgium1[,"Accession ID"][bassurv]

GISAID_belgium1[,"Gender"][grepl("F|V",GISAID_belgium1[,"Gender"])] = "F"
GISAID_belgium1[,"Gender"][grepl("M|m",GISAID_belgium1[,"Gender"])] = "M"
GISAID_belgium1[,"Gender"][grepl("unknown",GISAID_belgium1[,"Gender"])] = NA
unique(GISAID_belgium1[,"Gender"]) # "F" "M" NA 

sum(!is.na(GISAID_belgium1[,"Gender"] )) # 20830 with gender info

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
sum(!is.na(GISAID_belgium1[,"Patient age"] )) # 20832 with age info

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
# we use genomic epidemiology GISAID data file version metadata_2021-06-01_08-07.tsv.gz with data for all countries :
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-06-01_08-07.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID_genepi_belgium = GISAID[GISAID$country=="Belgium",]
nrow(GISAID_genepi_belgium) # 22822

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
write.csv(GISAID_belgium1, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_06_03_ALL PARSED.csv",row.names=F)




# ANALYSIS OF VOC LINEAGE FREQUENCIES IN BELGIUM USING MULTINOMIAL FITS ####

# GISAID_belgium = GISAID_sel[GISAID_sel$country=="Belgium",]
GISAID_belgium = GISAID_belgium1[GISAID_belgium1$purpose_of_sequencing=="baseline_surveillance",]
nrow(GISAID_belgium) # 15374
sum(!is.na(GISAID_belgium$age)) # 14129 with age info

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
sum(GISAID_belgium$LINEAGE1=="B.1.617+") # 44 among baseline surveillance
unique(GISAID_belgium$province[GISAID_belgium$LINEAGE1=="B.1.1.7"])
sum(GISAID_belgium$LINEAGE1=="B.1.1.7") # 10632 among baseline surveillance

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
GISAID_belgium = GISAID_belgium[GISAID_belgium$date>=as.Date("2021-02-01"),]
nrow(GISAID_belgium) # 13233
range(GISAID_belgium$date) # "2021-02-01" "2021-05-30"

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

write.csv(data_agbyweek2, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_06_03_BASELINE SELECTION_aggregated counts by week.csv",row.names=F)


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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +
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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +
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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_raw data.pdf"), width=8, height=6)



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
# 1         B.1 - B.1.1.7 -0.025571153 0.022693021 NA  -0.070048658  0.01890635 8.548631e-01
# 2       B.1.1 - B.1.1.7  0.046574141 0.015848318 NA   0.015512009  0.07763627 3.877118e-02
# 3     B.1.221 - B.1.1.7 -0.052429661 0.019347071 NA  -0.090349223 -0.01451010 7.179123e-02
# 4     B.1.258 - B.1.1.7 -1.569758630 1.441235173 NA  -4.394527662  1.25501040 8.710519e-01
# 5     B.1.160 - B.1.1.7 -0.065320057 0.024071753 NA  -0.112499827 -0.01814029 7.112924e-02
# 6  (B.1.177+) - B.1.1.7 -0.092408925 0.030863963 NA  -0.152901180 -0.03191667 3.312954e-02
# 7   B.1.214.2 - B.1.1.7 -0.027649085 0.006442588 NA  -0.040276326 -0.01502184 3.696292e-04
# 8     B.1.351 - B.1.1.7 -0.058741562 0.007234555 NA  -0.072921029 -0.04456210 1.224243e-12
# 9     B.1.620 - B.1.1.7 -0.091846291 0.052142084 NA  -0.194042897  0.01035031 4.749712e-01
# 10        P.1 - B.1.1.7  0.005288836 0.003857032 NA  -0.002270807  0.01284848 7.264389e-01
# 11  B.1.617.1 - B.1.1.7 -2.496613492 3.960172621 NA -10.258409201  5.26518222 9.848196e-01
# 12  B.1.617.2 - B.1.1.7  0.020494209 0.028238384 NA  -0.034852008  0.07584043 9.727602e-01
# 13      other - B.1.1.7  0.026536966 0.004282854 NA   0.018142727  0.03493120 5.533290e-08


# avg growth advantage of B.1.617+ over B.1.1.7 :
emtrbelgium1 = emtrends(fit2_belgium_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium1 = data.frame(confint(emtrbelgium1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrbelgium1$contrasts)$p.value)
delta_r_belgium1
#                contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1         B.1 - B.1.1.7 -0.025490836 0.022672610 NA -0.069928334  0.01894666 8.424097e-01
# 2       B.1.1 - B.1.1.7  0.046589653 0.015850809 NA  0.015522638  0.07765667 3.651131e-02
# 3     B.1.221 - B.1.1.7 -0.052417537 0.019345943 NA -0.090334888 -0.01450019 6.784953e-02
# 4     B.1.258 - B.1.1.7 -1.123758508 1.181130072 NA -3.438730911  1.19121389 9.114190e-01
# 5     B.1.160 - B.1.1.7 -0.065317572 0.024074736 NA -0.112503187 -0.01813196 6.721726e-02
# 6  (B.1.177+) - B.1.1.7 -0.092308151 0.030845956 NA -0.152765114 -0.03185119 3.137967e-02
# 7   B.1.214.2 - B.1.1.7 -0.027628785 0.006442662 NA -0.040256170 -0.01500140 3.602238e-04
# 8     B.1.351 - B.1.1.7 -0.058749336 0.007234548 NA -0.072928790 -0.04456988 1.626699e-12
# 9         P.1 - B.1.1.7  0.005299635 0.003857641 NA -0.002261201  0.01286047 7.071449e-01
# 10    B.1.620 - B.1.1.7 -0.091866930 0.052134417 NA -0.194048510  0.01031465 4.565968e-01
# 11 (B.1.617+) - B.1.1.7  0.038481841 0.022408932 NA -0.005438858  0.08240254 4.854358e-01
# 12      other - B.1.1.7  0.026543835 0.004283296 NA  0.018148728  0.03493894 5.869722e-08


# pairwise growth rate difference (differences in growth rate per day) 
emtrbelgium_pairw = emtrends(fit2_belgium_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1           B.1.1.7 - B.1  0.0255711528 0.022693021 NA  -0.018906352  0.070048658 9.974428e-01
# 2         B.1.1.7 - B.1.1 -0.0465741413 0.015848318 NA  -0.077636273 -0.015512009 1.721120e-01
# 3       B.1.1.7 - B.1.221  0.0524296610 0.019347071 NA   0.014510099  0.090349223 2.826766e-01
# 4       B.1.1.7 - B.1.258  1.5697586302 1.441235173 NA  -1.255010401  4.394527662 9.981823e-01
# 5       B.1.1.7 - B.1.160  0.0653200575 0.024071753 NA   0.018140288  0.112499827 2.806522e-01
# 6    B.1.1.7 - (B.1.177+)  0.0924089245 0.030863963 NA   0.031916669  0.152901180 1.508698e-01
# 7     B.1.1.7 - B.1.214.2  0.0276490847 0.006442588 NA   0.015021844  0.040276326 2.339247e-03
# 8       B.1.1.7 - B.1.351  0.0587415621 0.007234555 NA   0.044562095  0.072921029 8.572254e-12
# 9       B.1.1.7 - B.1.620  0.0918462914 0.052142084 NA  -0.010350315  0.194042897 8.897959e-01
# 10          B.1.1.7 - P.1 -0.0052888364 0.003857032 NA  -0.012848480  0.002270807 9.839458e-01
# 11    B.1.1.7 - B.1.617.1  2.4966134923 3.960172621 NA  -5.265182216 10.258409201 9.999960e-01
# 12    B.1.1.7 - B.1.617.2 -0.0204942090 0.028238384 NA  -0.075840426  0.034852008 9.999789e-01
# 13        B.1.1.7 - other -0.0265369658 0.004282854 NA  -0.034931205 -0.018142727 3.838710e-07
# 14            B.1 - B.1.1 -0.0721452941 0.027620875 NA  -0.126281215 -0.018009373 3.409674e-01
# 15          B.1 - B.1.221  0.0268585082 0.029692338 NA  -0.031337406  0.085054422 9.997447e-01
# 16          B.1 - B.1.258  1.5441874774 1.441376152 NA  -1.280857868  4.369232823 9.984635e-01
# 17          B.1 - B.1.160  0.0397489046 0.032943744 NA  -0.024819647  0.104317456 9.950164e-01
# 18       B.1 - (B.1.177+)  0.0668377717 0.038131037 NA  -0.007897687  0.141573230 8.933067e-01
# 19        B.1 - B.1.214.2  0.0020779318 0.023528252 NA  -0.044036594  0.048192457 1.000000e+00
# 20          B.1 - B.1.351  0.0331704093 0.023719500 NA  -0.013318956  0.079659774 9.809759e-01
# 21          B.1 - B.1.620  0.0662751386 0.056861163 NA  -0.045170694  0.177720971 9.964317e-01
# 22              B.1 - P.1 -0.0308599892 0.022988107 NA  -0.075915852  0.014195873 9.866750e-01
# 23        B.1 - B.1.617.1  2.4710423395 3.960237502 NA  -5.290880535 10.232965214 9.999965e-01
# 24        B.1 - B.1.617.2 -0.0460653618 0.036221551 NA  -0.117058297  0.024927573 9.918293e-01
# 25            B.1 - other -0.0521081186 0.023032485 NA  -0.097250959 -0.006965279 5.848205e-01
# 26        B.1.1 - B.1.221  0.0990038023 0.024921692 NA   0.050158184  0.147849421 7.703817e-03
# 27        B.1.1 - B.1.258  1.6163327715 1.441302659 NA  -1.208568532  4.441234075 9.975621e-01
# 28        B.1.1 - B.1.160  0.1118941987 0.028691573 NA   0.055659748  0.168128649 9.971634e-03
# 29     B.1.1 - (B.1.177+)  0.1389830658 0.034587875 NA   0.071192076  0.206774055 6.533728e-03
# 30      B.1.1 - B.1.214.2  0.0742232259 0.017029001 NA   0.040846997  0.107599455 1.799868e-03
# 31        B.1.1 - B.1.351  0.1053157034 0.017362086 NA   0.071286641  0.139344766 7.480827e-07
# 32        B.1.1 - B.1.620  0.1384204327 0.054484922 NA   0.031631948  0.245208917 3.871089e-01
# 33            B.1.1 - P.1  0.0412853049 0.016248097 NA   0.009439620  0.073130990 3.868401e-01
# 34      B.1.1 - B.1.617.1  2.5431876336 3.960204257 NA  -5.218670081 10.305045348 9.999950e-01
# 35      B.1.1 - B.1.617.2  0.0260799323 0.032333872 NA  -0.037293293  0.089453157 9.999290e-01
# 36          B.1.1 - other  0.0200371755 0.016343716 NA  -0.011995920  0.052070271 9.941993e-01
# 37      B.1.221 - B.1.258  1.5173289692 1.441282583 NA  -1.307532986  4.342190924 9.987158e-01
# 38      B.1.221 - B.1.160  0.0128903965 0.030549062 NA  -0.046984664  0.072765457 1.000000e+00
# 39   B.1.221 - (B.1.177+)  0.0399792635 0.035996835 NA  -0.030573237  0.110531764 9.977872e-01
# 40    B.1.221 - B.1.214.2 -0.0247805763 0.020297814 NA  -0.064563561  0.015002408 9.944255e-01
# 41      B.1.221 - B.1.351  0.0063119011 0.020526902 NA  -0.033920088  0.046543891 1.000000e+00
# 42      B.1.221 - B.1.620  0.0394166304 0.055612756 NA  -0.069582369  0.148415630 9.999840e-01
# 43          B.1.221 - P.1 -0.0577184974 0.019698282 NA  -0.096326420 -0.019110575 1.756083e-01
# 44    B.1.221 - B.1.617.1  2.4441838313 3.960219944 NA  -5.317704630 10.206072292 9.999969e-01
# 45    B.1.221 - B.1.617.2 -0.0729238700 0.034229950 NA  -0.140013339 -0.005834401 6.793243e-01
# 46        B.1.221 - other -0.0789666268 0.019701927 NA  -0.117581694 -0.040351560 6.779665e-03
# 47      B.1.258 - B.1.160 -1.5044385728 1.441321281 NA  -4.329376375  1.320499229 9.988245e-01
# 48   B.1.258 - (B.1.177+) -1.4773497057 1.441352877 NA  -4.302349433  1.347650022 9.990269e-01
# 49    B.1.258 - B.1.214.2 -1.5421095456 1.441247650 NA  -4.366903032  1.282683941 9.984832e-01
# 50      B.1.258 - B.1.351 -1.5110170681 1.441246370 NA  -4.335808046  1.313773910 9.987696e-01
# 51      B.1.258 - B.1.620 -1.4779123388 1.442177604 NA  -4.304528501  1.348703824 9.990288e-01
# 52          B.1.258 - P.1 -1.5750474666 1.441240429 NA  -4.399826801  1.249731868 9.981194e-01
# 53    B.1.258 - B.1.617.1  0.9268548621 4.214274012 NA  -7.332970423  9.186680147 1.000000e+00
# 54    B.1.258 - B.1.617.2 -1.5902528392 1.441511292 NA  -4.415563055  1.235057377 9.979318e-01
# 55        B.1.258 - other -1.5962955960 1.441229633 NA  -4.421053770  1.228462578 9.978472e-01
# 56   B.1.160 - (B.1.177+)  0.0270888671 0.038721637 NA  -0.048804147  0.102981881 9.999863e-01
# 57    B.1.160 - B.1.214.2 -0.0376709728 0.024838624 NA  -0.086353782  0.011011837 9.628575e-01
# 58      B.1.160 - B.1.351 -0.0065784953 0.024997881 NA  -0.055573442  0.042416451 1.000000e+00
# 59      B.1.160 - B.1.620  0.0265262340 0.057427554 NA  -0.086029703  0.139082171 9.999999e-01
# 60          B.1.160 - P.1 -0.0706088938 0.024357414 NA  -0.118348549 -0.022869239 1.887295e-01
# 61    B.1.160 - B.1.617.1  2.4312934349 3.960245885 NA  -5.330645869 10.193232739 9.999971e-01
# 62    B.1.160 - B.1.617.2 -0.0858142665 0.037105587 NA  -0.158539881 -0.013088652 5.480569e-01
# 63        B.1.160 - other -0.0918570233 0.024355932 NA  -0.139593772 -0.044120274 1.552666e-02
# 64 (B.1.177+) - B.1.214.2 -0.0647598398 0.031454770 NA  -0.126410056 -0.003109624 7.278787e-01
# 65   (B.1.177+) - B.1.351 -0.0336673624 0.031557415 NA  -0.095518760  0.028184035 9.985278e-01
# 66   (B.1.177+) - B.1.620 -0.0005626331 0.060615153 NA  -0.119366150  0.118240884 1.000000e+00
# 67       (B.1.177+) - P.1 -0.0976977609 0.031091703 NA  -0.158636380 -0.036759142 1.037820e-01
# 68 (B.1.177+) - B.1.617.1  2.4042045678 3.960293180 NA  -5.357827433 10.166236569 9.999975e-01
# 69 (B.1.177+) - B.1.617.2 -0.1129031335 0.041831553 NA  -0.194891471 -0.030914796 2.888878e-01
# 70     (B.1.177+) - other -0.1189458903 0.031067821 NA  -0.179837700 -0.058054081 1.277731e-02
# 71    B.1.214.2 - B.1.351  0.0310924775 0.009581814 NA   0.012312468  0.049872487 7.868469e-02
# 72    B.1.214.2 - B.1.620  0.0641972068 0.052523704 NA  -0.038747361  0.167141774 9.943645e-01
# 73        B.1.214.2 - P.1 -0.0329379210 0.007368181 NA  -0.047379291 -0.018496551 1.153193e-03
# 74  B.1.214.2 - B.1.617.1  2.4689644077 3.960176786 NA  -5.292839465 10.230768280 9.999965e-01
# 75  B.1.214.2 - B.1.617.2 -0.0481432937 0.028949302 NA  -0.104882882  0.008596295 9.256815e-01
# 76      B.1.214.2 - other -0.0541860505 0.007590943 NA  -0.069064026 -0.039308075 2.419172e-09
# 77      B.1.351 - B.1.620  0.0331047293 0.052623924 NA  -0.070036266  0.136245724 9.999961e-01
# 78          B.1.351 - P.1 -0.0640303985 0.008132429 NA  -0.079969667 -0.048091130 3.617107e-11
# 79    B.1.351 - B.1.617.1  2.4378719302 3.960179067 NA  -5.323936413 10.199680274 9.999970e-01
# 80    B.1.351 - B.1.617.2 -0.0792357711 0.029151155 NA  -0.136370985 -0.022100558 2.781114e-01
# 81        B.1.351 - other -0.0852785279 0.008268975 NA  -0.101485422 -0.069071634 7.571721e-14
# 82          B.1.620 - P.1 -0.0971351278 0.052268805 NA  -0.199580103  0.005309848 8.454401e-01
# 83    B.1.620 - B.1.617.1  2.4047672009 3.960513965 NA  -5.357697530 10.167231932 9.999975e-01
# 84    B.1.620 - B.1.617.2 -0.1123405004 0.059308578 NA  -0.228583177  0.003902176 8.268524e-01
# 85        B.1.620 - other -0.1183832572 0.052307778 NA  -0.220904618 -0.015861896 5.842147e-01
# 86        P.1 - B.1.617.1  2.5019023287 3.960174251 NA  -5.259896575 10.263701232 9.999959e-01
# 87        P.1 - B.1.617.2 -0.0152053726 0.028445164 NA  -0.070956870  0.040546125 9.999995e-01
# 88            P.1 - other -0.0212481294 0.005605752 NA  -0.032235202 -0.010261057 1.455939e-02
# 89  B.1.617.1 - B.1.617.2 -2.5171077013 3.960272718 NA -10.279099597  5.244884195 9.999956e-01
# 90      B.1.617.1 - other -2.5231504581 3.960174822 NA -10.284950481  5.238649565 9.999955e-01
# 91      B.1.617.2 - other -0.0060427568 0.028525655 NA  -0.061952013  0.049866500 1.000000e+00

# estimated proportion of different LINEAGES in Belgium today
today # "2021-06-03"
# 62% [58%-65%] 95% CLs now estimated to be B.1.617.2 across all provinces
multinom_preds_today_avg = data.frame(emmeans(fit2_belgium_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1    B.1.1.7 8.135461e-01 1.694109e-02 NA  7.803422e-01 8.467500e-01
# 2        B.1 2.271596e-04 2.522467e-04 NA -2.672349e-04 7.215540e-04
# 3      B.1.1 7.730407e-03 4.567008e-03 NA -1.220763e-03 1.668158e-02
# 4    B.1.221 1.832953e-04 1.937086e-04 NA -1.963666e-04 5.629572e-04
# 5    B.1.258 9.843938e-48 8.823933e-46 NA -1.719615e-45 1.739303e-45
# 6    B.1.160 6.493338e-05 8.708438e-05 NA -1.057489e-04 2.356156e-04
# 7   B.1.177+ 9.647403e-06 1.696255e-05 NA -2.359858e-05 4.289339e-05
# 8  B.1.214.2 5.622821e-03 1.579650e-03 NA  2.526763e-03 8.718879e-03
# 9    B.1.351 1.350115e-03 4.631194e-04 NA  4.424177e-04 2.257812e-03
# 10   B.1.620 2.291534e-05 4.999045e-05 NA -7.506414e-05 1.208948e-04
# 11       P.1 8.338785e-02 1.017312e-02 NA  6.344891e-02 1.033268e-01
# 12 B.1.617.1 1.410855e-55 2.658013e-53 NA -5.195501e-53 5.223718e-53
# 13 B.1.617.2 1.900155e-02 1.103772e-02 NA -2.631991e-03 4.063509e-02
# 14     other 6.885321e-02 1.056597e-02 NA  4.814430e-02 8.956213e-02


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2021-02-01"))
date.to = as.numeric(as.Date("2021-06-01")) # max(GISAID_belgium$DATE_NUM)+extrapolate

fit_belgium_multi_predsbyprovince = data.frame(emmeans(fit2_belgium_multi,
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "province"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=2)), # by=2 just to speed up things a bit
                                                  mode="prob", df=NA))
fit_belgium_multi_predsbyprovince$collection_date = as.Date(fit_belgium_multi_predsbyprovince$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_predsbyprovince$LINEAGE2 = factor(fit_belgium_multi_predsbyprovince$LINEAGE2, levels=levels_LINEAGE2)
fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=levels_PROVINCES)

fit_belgium_multi_preds = data.frame(emmeans(fit2_belgium_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=2)), # by=2 just to speed up things a bit
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
