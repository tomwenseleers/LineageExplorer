# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN BELGIUM BASED ON ANALYSIS OF GISAID PATIENT METADATA
# T. Wenseleers
# last update 9 Juni 2022

library(remotes)
# remotes::install_github("rvlenth/emmeans")
# remotes::install_github("vincentarelbundock/marginaleffects")
# devtools::install_github("melff/mclogit",subdir="pkg")

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(marginaleffects)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
# today = as.Date("2021-08-31")
today # "2022-01-17"
today_num = as.numeric(today)
plotdir = "BE_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))



# import & parse manually downloaded GISAID patient metadata for Belgium (downloaded 30/6/2021) ####
d1 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_jan_dec_2020.tsv", col_types = cols(.default = "c"))
d2 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_jan_feb_2021.tsv", col_types = cols(.default = "c"))
d3 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_mar_apr_2021.tsv", col_types = cols(.default = "c"))
d4 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_may_2021.tsv", col_types = cols(.default = "c")) 
d5 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_june_2021.tsv", col_types = cols(.default = "c"))
d6 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_july_2021.tsv", col_types = cols(.default = "c"))
d7 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_aug_2021.tsv", col_types = cols(.default = "c"))
d8 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_sept_2021.tsv", col_types = cols(.default = "c"))
d9 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_oct_1_2021.tsv", col_types = cols(.default = "c"))
d10 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_12_01_17_belgium_subm_oct_2_2021.tsv", col_types = cols(.default = "c"))
d11 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_01_17_08_belgium_subm_nov_2021.tsv", col_types = cols(.default = "c"))
d12 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_01_17_08_belgium_subm_dec_2021.tsv", col_types = cols(.default = "c"))
d13 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_04_13_10_belgium_subm_jan_2022.tsv", col_types = cols(.default = "c"))
d14 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_04_13_10_belgium_subm_feb_2022.tsv", col_types = cols(.default = "c"))
d15 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_04_13_10_belgium_subm_mar_2022.tsv", col_types = cols(.default = "c"))
d16 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_05_16_00_subm_apr_2022.tsv", col_types = cols(.default = "c"))
d17 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_06_09_07_subm_may_june_2022.tsv", col_types = cols(.default = "c"))

d1 = as.data.frame(d1)
d2 = as.data.frame(d2)
d3 = as.data.frame(d3)
d4 = as.data.frame(d4)
d5 = as.data.frame(d5)
d6 = as.data.frame(d6)
d7 = as.data.frame(d7)
d8 = as.data.frame(d8)
d9 = as.data.frame(d9)
d10 = as.data.frame(d10)
d11 = as.data.frame(d11)
d12 = as.data.frame(d12)
d13 = as.data.frame(d13)
d14 = as.data.frame(d14)
d15 = as.data.frame(d15)
d16 = as.data.frame(d16)
d17 = as.data.frame(d17)
GISAID_belgium1 = rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17)
GISAID_belgium1 = GISAID_belgium1[GISAID_belgium1[,"Host"]=="Human",]
GISAID_belgium1[,"Collection date"] = as.Date(GISAID_belgium1[,"Collection date"])
GISAID_belgium1 = GISAID_belgium1[!is.na(GISAID_belgium1[,"Collection date"]),]
nrow(GISAID_belgium1) # 110967
unique(GISAID_belgium1[,"Virus name"])
unique(GISAID_belgium1[,"Location"])
unique(GISAID_belgium1[,"Host"])
unique(GISAID_belgium1[,"Additional location information"]) # ZIP codes # sometimes has additional info on whether it was a traveller
range(GISAID_belgium1$`Collection date`)

# anonymous records without location info are assumed to be active surveillance of cluster outbreaks & included in active surveillance (cf weekly Sciensano report)
GISAID_belgium1[,"Additional location information"][grepl("nknown",GISAID_belgium1[,"Additional location information"])] = NA
sum(is.na(GISAID_belgium1[,"Additional location information"])) # 10891
unique(GISAID_belgium1[,"Additional location information"])

library(stringr)
ZIP = str_extract(GISAID_belgium1[,"Additional location information"],"[0-9]{4}") # extract ZIP codes
unique(ZIP)
sum(!is.na(ZIP)) # 95906
sum(is.na(ZIP)) # 10341 - anonymous records, these are assumed to be active surveillance of cluster outbreaks

muni = read.csv(".//data//GISAID//Belgium//mapping_municips.csv") # post codes & cities & provinces
province = muni$province_label_nl[match(ZIP, muni$postcode)] # province
city = muni$municipality_label_nl[match(ZIP, muni$postcode)] # city
sum(is.na(city)) # 10490
sum(is.na(ZIP))  # 10341
wrongZIP = which(is.na(city)&!is.na(ZIP))
ZIP[wrongZIP] = NA

anonym = grepl("Belgium$",GISAID_belgium1[,"Additional location information"])|
  is.na(GISAID_belgium1[,"Additional location information"])|
  grepl("CHU|infected|travel|Travel",GISAID_belgium1[,"Additional location information"]) |
  is.na(ZIP)
sum(anonym) # 10630
unique(GISAID_belgium1[!is.na(GISAID_belgium1[,"Sampling strategy"]),"Sampling strategy"])
traveller = grepl("travel|Travel|infected",GISAID_belgium1[,"Additional location information"])|grepl("travel|Travel",GISAID_belgium1[,"Sampling strategy"])|
  grepl("travel|Travel|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(traveller) # 3391
Sdropout = grepl("S-gene dropout",GISAID_belgium1[,"Sampling strategy"])|
  grepl("S-gene dropout",GISAID_belgium1[,"Additional host information"])  
sum(Sdropout) # 21 - should be much higher though...
actsurv = anonym|traveller|grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster",GISAID_belgium1[,"Sampling strategy"])|
  grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster|HR contact|Outbreak|Nurse|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(actsurv) # 23635
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")]) # 19449 # In Sciensano weekly report of 21st of May 2021 this is 4710
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")&GISAID_belgium1[,"Collection date"]<=as.Date("2021-01-31")]) # 1204 # In Sciensano weekly report of 21st of May 2021 this is 1975

reinfection = grepl("Breakthrough",GISAID_belgium1[,"Sampling strategy"])
sum(reinfection) # 354
infectionpostvaccination = grepl("Vaccination|vaccination",GISAID_belgium1[,"Sampling strategy"])
sum(infectionpostvaccination) # 5854

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
sum(bassurv) # 82612
purpose_of_sequencing = rep("active_surveillance", nrow(GISAID_belgium1))
purpose_of_sequencing[bassurv] = "baseline_surveillance"

accessionIDs_bassurv = GISAID_belgium1[,"Accession ID"][bassurv]

GISAID_belgium1[,"Gender"][grepl("F|V",GISAID_belgium1[,"Gender"])] = "F"
GISAID_belgium1[,"Gender"][grepl("M|m",GISAID_belgium1[,"Gender"])] = "M"
GISAID_belgium1[,"Gender"][grepl("unknown",GISAID_belgium1[,"Gender"])] = NA
GISAID_belgium1[,"Gender"][!grepl("F|M",GISAID_belgium1[,"Gender"])] = NA
unique(GISAID_belgium1[,"Gender"]) # "F" "M" NA 

sum(!is.na(GISAID_belgium1[,"Gender"] )) # 93474 with gender info

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
GISAID_belgium1[,"Patient age"][GISAID_belgium1[,"Patient age"]>200] = NA
sum(!is.na(GISAID_belgium1[,"Patient age"] )) # 93225 with age info

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
# GISAID_belgium1$pango_lineage_simplified[grepl("B.1.617.2|AY",GISAID_belgium1$pango_lineage_simplified)] = "Delta"
GISAID_belgium1$LINEAGE2 = GISAID_belgium1$pango_lineage_simplified


# get extra fields "genbank_accession", "Nextstrain_clade", "originating_lab", "submitting_lab", "authors", "url", "title", "paper_url"
# from genomic epidemiology GISAID metadata
# (check with Emmanuel André & Lize Cuypers which labs were doing active surveillance vs baseline surveillance)
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2022-05-12_23-45.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID_genepi_belgium = GISAID[GISAID$country=="Belgium",]
nrow(GISAID_genepi_belgium) # 103271

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
write.csv(GISAID_belgium1, ".//data//GISAID//Belgium//gisaid_hcov-19_2022_05_16_00_ALL PARSED.csv",row.names=F)
range(GISAID_belgium1$date) # "2020-02-03" "2022-05-10"




# ANALYSIS OF VOC LINEAGE FREQUENCIES IN BELGIUM USING MULTINOMIAL FITS ####

GISAID_belgium1 = read.csv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_05_16_00_ALL PARSED.csv")
GISAID_belgium1$date = as.Date(GISAID_belgium1$date)
GISAID_belgium1$province = as.factor(GISAID_belgium1$province)
str(GISAID_belgium1)

library(dplyr)

GISAID_belgium1$LINEAGE = case_when(
  # grepl("N679K", AA_substitutions) & grepl("H655Y", AA_substitutions) & grepl("P681H", AA_substitutions) ~ "B.1.1.529",
  (grepl("^BA\\.1\\.1$|BA\\.1\\.1\\.", GISAID_belgium1$pango_lineage)&grepl("R346K", GISAID_belgium1$AA_substitutions)) ~ "Omicron (BA.1.1, with S:R346K)",
  (grepl("B\\.1\\.1\\.529|^BA\\.1$|BA\\.1\\.", GISAID_belgium1$pango_lineage)) ~ "Omicron (BA.1)",
  (grepl("^BA\\.3$|BA\\.3\\.", GISAID_belgium1$pango_lineage)) ~ "Omicron (BA.3)",
  (((grepl("BA\\.4",GISAID_belgium1$pango_lineage)))) ~ "Omicron (BA.4)",
  (((grepl("BA\\.5",GISAID_belgium1$pango_lineage)))) ~ "Omicron (BA.5)",
  (grepl("^BA\\.2",GISAID_belgium1$pango_lineage)) ~ "Omicron (BA.2)",
  grepl("AY", GISAID_belgium1$pango_lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", GISAID_belgium1$pango_lineage) ~ "Alpha",
  grepl("B.1.351", GISAID_belgium1$pango_lineage, fixed=T) ~ "Beta",
  grepl("P.1", GISAID_belgium1$pango_lineage, fixed=T) ~ "Gamma",
  # grepl("C.1.2", GISAID_belgium1$pango_lineage, fixed=T) ~ "C.1.2",
  T ~ "Other"
)

# GISAID_belgium1$LINEAGE[GISAID_belgium1$LINEAGE=="Alpha"&GISAID_belgium1$date<=as.Date("2020-10-01")] = "Other"
table(GISAID_belgium1$LINEAGE)

# GISAID_belgium1$LINEAGE = case_when(
#   grepl("B.1.1.529|BA.1", GISAID_belgium1$pango_lineage) ~ "Omicron (BA.1)", # B.1.1.529|
#   # (GISAID_belgium1$pango_lineage=="BA.3") ~ "Omicron (BA.3)",
#   # ((grepl("BA.2",GISAID_belgium1$pango_lineage))&
#   #    ((grepl("L452R", GISAID_belgium1$AA_substitutions)&
#   #        grepl("486V", GISAID_belgium1$AA_substitutions)&
#   #        grepl("11F", GISAID_belgium1$AA_substitutions)&
#   #        (!grepl("D3N",GISAID_belgium1$AA_substitutions)) ))) ~ "Omicron (BA.4)",
#   # ((grepl("BA.2",GISAID_belgium1$pango_lineage))&
#   #    ((grepl("L452R",GISAID_belgium1$AA_substitutions)&
#   #        grepl("486V",GISAID_belgium1$AA_substitutions)&
#   #        (!grepl("11F", GISAID_belgium1$AA_substitutions))&
#   #        grepl("D3N",GISAID_belgium1$AA_substitutions)))) ~ "Omicron (BA.5)",
#   (grepl("BA.2",GISAID_belgium1$pango_lineage)) ~ "Omicron (BA.2)",
#   grepl("B.1.617.2", GISAID_belgium1$pango_lineage, fixed=T) | grepl("AY", GISAID_belgium1$pango_lineage)  ~ "Delta",
#   grepl("B.1.1.7", GISAID_belgium1$pango_lineage, fixed=T) ~ "Alpha",
#   grepl("B.1.351", GISAID_belgium1$pango_lineage, fixed=T) ~ "Beta",
#   grepl("P.1", GISAID_belgium1$pango_lineage, fixed=T) ~ "Gamma",
#   T ~ "Other"
# )



# GISAID_belgium1$VOC = GISAID_json$covv_variant
# GISAID_belgium1$VOC[GISAID_json$covv_lineage=="A"] = "A" 
# GISAID_belgium1$VOC[GISAID_json$VOC==""&(grepl("S",GISAID_json$covv_clade)|grepl("A.",GISAID_json$covv_lineage,fixed=T))] = "A*" # = 19B, cf https://www.news-medical.net/health/Viral-Clades-of-SARS-CoV-2.aspx
# GISAID_json$VOC[GISAID_json$covv_lineage=="B"] = "B"
# GISAID_json$VOC[GISAID_json$VOC==""&(grepl("L|O|V|G",GISAID_json$covv_clade)|grepl("B.",GISAID_json$covv_lineage,fixed=T))&(!grepl("D614G",GISAID_json$covsurver_prot_mutations))] = "B*" # B* without D614G
# GISAID_json$VOC[GISAID_json$VOC==""&(grepl("G",GISAID_json$covv_clade)|grepl("B.",GISAID_json$covv_lineage,fixed=T)|grepl("AE.",GISAID_json$covv_lineage,fixed=T)|grepl("AH.",GISAID_json$covv_lineage,fixed=T)|grepl("AK.",GISAID_json$covv_lineage,fixed=T)|grepl("AN.",GISAID_json$covv_lineage,fixed=T)|grepl("AS.",GISAID_json$covv_lineage,fixed=T)|grepl("AU.",GISAID_json$covv_lineage,fixed=T)|grepl("AV.",GISAID_json$covv_lineage,fixed=T)|grepl("AZ.",GISAID_json$covv_lineage,fixed=T))&(grepl("D614G",GISAID_json$covsurver_prot_mutations))] = "B* (+S:D614G)" # B* with D614G
# GISAID_json$VOC[grepl("Alpha",GISAID_json$VOC)] = "Alpha"
# GISAID_json$VOC[grepl("Beta",GISAID_json$VOC)] = "Beta"
# GISAID_json$VOC[grepl("Gamma",GISAID_json$VOC)|grepl("P.1.",GISAID_json$covv_lineage,fixed=T)] = "Gamma"
# GISAID_json$VOC[grepl("Delta",GISAID_json$VOC)|grepl("AY.",GISAID_json$covv_lineage,fixed=T)] = "Delta"
# GISAID_json$VOC[grepl("Lambda",GISAID_json$VOC)] = "Lambda"
# GISAID_json$VOC[grepl("Mu",GISAID_json$VOC)] = "Mu"




# selected lineages
main_lineages = c("Alpha","Beta","Gamma","Delta",
                  "Omicron (BA.1)", "Omicron (BA.2)",
                  "Omicron (BA.3)", "Omicron (BA.1.1, with S:R346K)", "Omicron (BA.4)", "Omicron (BA.5)") 

levels_LINEAGE = c(main_lineages, "Other")
levels_LINEAGE_plot = c("Other", "Alpha", "Beta", "Gamma", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.1.1, with S:R346K)", "Omicron (BA.4)", "Omicron (BA.5)")


# lineage_cols = c("#0085FF","#9A9D00","cyan3","magenta",
#                   muted("red",c=150,l=65),muted("red",c=150,l=25),
#                   "grey70")
# lineage_cols = c("#0085FF","#9A9D00","cyan3","magenta",
#                   colorRampPalette(c("red", "orange", "blue"))(2),
#                   "grey70")
n = length(levels_LINEAGE)
lineage_cols = hcl(h = seq(0, 180, length = n), l = 60, c = 180)
lineage_cols[which(levels_LINEAGE=="Alpha")] = "#0085FF"
lineage_cols[which(levels_LINEAGE=="Beta")] = "green4"
lineage_cols[which(levels_LINEAGE=="Delta")] = "mediumorchid"
# lineage_cols[which(levels_LINEAGE=="C.1.2")] = "darkorange"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.1)")] = "red" # "magenta"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.1.1, with S:R346K)")] = "black"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.2)")] = "red3"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.3)")] = "red4" 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.4)")] = "darkorange" 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.5)")] = "darkorange3" 
lineage_cols[which(levels_LINEAGE=="Other")] = "grey65"

lineage_cols_plot = lineage_cols[match(levels_LINEAGE_plot,levels_LINEAGE)]


unique(GISAID_belgium1$LINEAGE)
GISAID_belgium1$LINEAGE[!(GISAID_belgium1$LINEAGE %in% main_lineages)] = "Other" # minority lineages & non-VOCs
table(GISAID_belgium1$LINEAGE)
GISAID_belgium1$LINEAGE = factor(GISAID_belgium1$LINEAGE, 
                                  levels=levels_LINEAGE, 
                                  labels=levels_LINEAGE)
table(GISAID_belgium1$LINEAGE)
sum(table(GISAID_belgium1$LINEAGE)) # 106247

# GISAID_belgium = GISAID_sel[GISAID_sel$country=="Belgium",]
GISAID_belgium = GISAID_belgium1 # [GISAID_belgium1$purpose_of_sequencing=="baseline_surveillance",] # if desired subset to baseline surveillance
# GISAID_belgium = GISAID_belgium[-which((GISAID_belgium$LINEAGE=="B.1.617.2")&(GISAID_belgium$date<="2021-05-01")),]
nrow(GISAID_belgium) # 106247
sum(!is.na(GISAID_belgium$age)) # 93225 with age info

# use data from Feb 1 onwards
GISAID_belgium_all = GISAID_belgium
# GISAID_belgium = GISAID_belgium_all[GISAID_belgium_all$date>=as.Date("2021-02-01"),]
nrow(GISAID_belgium) # 106247
range(GISAID_belgium$date) # "2020-02-03" "2022-05-10"


# # plot age distribution of B.1.617.2 & B.1.1.7 cases in Belgium though time
# df = GISAID_belgium[(GISAID_belgium$date>=as.Date("2021-04-01")&GISAID_belgium$pango_lineage %in% c("B.1.1.7","B.1.617.2")),]
# df$variant = factor(df$pango_lineage, levels=c("B.1.1.7","B.1.617.2"), labels=c("B.1.1.7 (alfa)","B.1.617.2 (delta)"))
# qplot(data=df[df$Week>18&df$Week<26,], binwidth=10.0, x=age, y = ..density..*100, geom="histogram", fill=variant) + 
#   facet_wrap(~floor_date+variant, ncol=2) + # , ncol=1
#   ylab("%") + scale_fill_manual(values=lineage_cols[c(1,4)]) + theme(legend.position="none")
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_B1672_B117_age distribution.png"), width=8, height=10)


# unique(GISAID_belgium$originating_lab) # 133 labs # maybe some labs did more active surveillance & could be excluded?
# unique(GISAID_belgium$submitting_lab) # 41 labs
# sum(GISAID_belgium$country_exposure!="Belgium") # only 23 are indicated as travellers

# # we remove travellers
# GISAID_belgium = GISAID_belgium[GISAID_belgium$country_exposure=="Belgium",]
# nrow(GISAID_belgium) # 19283
# sum(GISAID_belgium$division_exposure==GISAID_belgium$division) # 19283

# unique(GISAID_belgium$province) # 
# unique(GISAID_belgium$province[GISAID_belgium$LINEAGE=="B.1.617.2"]) # "Belgium"         "Hasselt"         "Brussels"        "Brugge"          "Gent"            "Liège"           "Halle-Vilvoorde" "Antwerpen"       "Mechelen"
# sum(GISAID_belgium$LINEAGE=="B.1.617.2") # 171 among baseline surveillance
# unique(GISAID_belgium$province[GISAID_belgium$LINEAGE1=="B.1.1.7"])
# sum(GISAID_belgium$LINEAGE1=="B.1.1.7") # 13339 among baseline surveillance
# 
# table(GISAID_belgium$LINEAGE1)
# table(GISAID_belgium$LINEAGE)

# main_lineages = names(table(GISAID_belgium$LINEAGE1))[100*table(GISAID_belgium$LINEAGE1)/sum(table(GISAID_belgium$LINEAGE1)) > 3]
# main_lineages
# # "B.1.1.7" "B.1.351" "P.1" 
# VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.620","B.1.1.7","B.1.351","P.1","B.1.1.207","B.1.429", # "B.1.1.318",
#          "B.1.214.2") # I added B.1.214.2 here # cut "B.1.525","B.1.526",
# main_lineages = union(main_lineages, VOCs)
# # main_lineages = c(main_lineages,"B.1","B.1.1","B.1.258")


# remove = names(table(GISAID_belgium$LINEAGE1))[table(GISAID_belgium$LINEAGE1) < 10]
# GISAID_belgium$LINEAGE1[(GISAID_belgium$LINEAGE1 %in% remove)] = "other" # minority VOCs
# GISAID_belgium$LINEAGE[(GISAID_belgium$LINEAGE %in% remove)] = "other" # minority VOCs
# table(GISAID_belgium$LINEAGE1)
# GISAID_belgium$LINEAGE1 = factor(GISAID_belgium$LINEAGE1)
# GISAID_belgium$LINEAGE1 = relevel(GISAID_belgium$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
# levels(GISAID_belgium$LINEAGE1)
# "B.1.1.7"   "B.1.214.2" "B.1.351"   "B.1.617+"  "B.1.620"   "other"     "P.1"  
# levels_LINEAGE1 = c("B.1.1.7","B.1.214.2",
#                    "B.1.351","P.1","B.1.620","B.1.617+","other")
# GISAID_belgium$LINEAGE1 = factor(GISAID_belgium$LINEAGE1, levels=levels_LINEAGE1)
# 
# GISAID_belgium$LINEAGE = factor(GISAID_belgium$LINEAGE)
# GISAID_belgium$LINEAGE = relevel(GISAID_belgium$LINEAGE, ref="B.1.1.7") # we code UK strain as the reference level
# levels(GISAID_belgium$LINEAGE)
# # "B.1.1.7"   "B.1.214.2" "B.1.351"   "B.1.617.1" "B.1.617.2" "B.1.620"   "other"     "P.1"  
# levels_LINEAGE = c("B.1.1.7","B.1.214.2",
#                     "B.1.351","P.1","B.1.620","B.1.617.1","B.1.617.2","other")
# GISAID_belgium$LINEAGE = factor(GISAID_belgium$LINEAGE, levels=levels_LINEAGE)

# # pango lineages since May 1st
# tbl = 100*table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01"),]$pango_lineage)/sum(table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01"),]$pango_lineage))
# tbl[order(tbl,decreasing=T)]
# 
# # pango lineages in category other since May 1st
# tbl = 100*table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01")&GISAID_belgium$LINEAGE1=="other",]$pango_lineage)/sum(table(GISAID_belgium[GISAID_belgium$date>=as.Date("2020-05-01")&GISAID_belgium$LINEAGE1=="other",]$pango_lineage))
# tbl[order(tbl,decreasing=T)]

# clip data from last few weeks to avoid deposition biases
# GISAID_belgium = GISAID_belgium[GISAID_belgium$date<=(max(GISAID_belgium$date)-14),]



# write.csv(GISAID_belgium, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_24_10_ALL PARSED_BASELINE SELECTION.csv",row.names=F)


write.csv(GISAID_belgium[,"gisaid_epi_isl"],
          paste0(".\\plots\\",plotdir,
                 "\\accession IDs Belgium.csv"), 
          row.names=F)


# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of Belgium
data_agbyweek2 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$LINEAGE))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE=="Delta","total"]) == nrow(GISAID_belgium) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE = factor(data_agbyweek2$LINEAGE, levels=levels_LINEAGE)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL

write.csv(data_agbyweek2, ".//data//GISAID//Belgium//gisaid_hcov-19_2022_05_16_00_BASELINE SELECTION_aggregated counts by week.csv",row.names=F)


# aggregated by week and province for selected variant lineages
data_agbyweekregion2 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$province, GISAID_belgium$LINEAGE))
colnames(data_agbyweekregion2) = c("floor_date", "province", "LINEAGE", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + province, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$province),
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$province))]
# sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE=="B.1.617.2","total"]) == nrow(GISAID_belgium) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE = factor(data_agbyweekregion2$LINEAGE, levels=levels_LINEAGE)
data_agbyweekregion2$province = factor(data_agbyweekregion2$province, levels=levels_PROVINCES)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]
data_agbyweekregion2$floor_date = NULL
data_agbyweekregion2$DATE_NUM = as.numeric(data_agbyweekregion2$collection_date)

# data_agbyweekregion2[data_agbyweekregion2$LINEAGE=="AY.5.2"&data_agbyweekregion2$count>0,]

# MULLER PLOT (RAW DATA)
# library(scales)
# n1 = length(levels(data_agbyweek1$LINEAGE1))
# lineage_cols1 = hcl(h = seq(15, 320, length = n1), l = 65, c = 200)
# lineage_cols1[which(levels(data_agbyweek1$LINEAGE1)=="B.1.617+")] = "magenta"
# lineage_cols1[which(levels(data_agbyweek1$LINEAGE1)=="other")] = "grey75"
# 
# n2 = length(levels(data_agbyweek2$LINEAGE))
# lineage_cols = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
# lineage_cols[which(levels(data_agbyweek2$LINEAGE)=="B.1.617.1")] = muted("magenta")
# lineage_cols[which(levels(data_agbyweek2$LINEAGE)=="B.1.617.2")] = "magenta"
# lineage_cols[which(levels(data_agbyweek2$LINEAGE)=="other")] = "grey75"

data_agbyweek2$LINEAGE = factor(data_agbyweek2$LINEAGE, levels=levels_LINEAGE_plot)
muller_belgium_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE)) + 
  # facet_wrap(~ province, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + xlab("Collection date") +
  theme(legend.position="right") + # ,   axis.title.x=element_blank()
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nGISAID data") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_belgium$Week)))
muller_belgium_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_raw data.png"), width=8, height=5)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_raw data.pdf"), width=8, height=5)

data_agbyweekregion2$LINEAGE = factor(data_agbyweekregion2$LINEAGE, levels=levels_LINEAGE_plot)
muller_belgiumbyprovince_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE)) +
  facet_wrap(~ province, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + xlab("Collection date") +
  theme(legend.position="right") + # axis.title.x=element_blank()
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nGISAID data")
# +
# coord_cartesian(xlim=c(1,max(GISAID_belgium$Week)))
muller_belgiumbyprovince_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
data_agbyweekregion2$LINEAGE = relevel(data_agbyweekregion2$LINEAGE, ref="Omicron (BA.1)") 
fit0_belgium_multi = nnet::multinom(LINEAGE ~ ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
fit1_belgium_multi = nnet::multinom(LINEAGE ~ province + scale(DATE_NUM), weights=count, data=data_agbyweekregion2, maxit=1000)
fit2_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
fit3_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
# fit4_belgium_multi = nnet::multinom(LINEAGE ~ province * DATE_NUM, weights=count, data=data_agbyweekregion2, maxit=1000)
# fit5_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
# fit6_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
BIC(fit1_belgium_multi, fit2_belgium_multi, fit3_belgium_multi) # fit4_belgium_multi, fit5_belgium_multi, fit6_belgium_multi 
# fit3_belgium_multi fits best (lowest BIC)

# library(mclogit)
# mblogit_fit1 = mblogit(LINEAGE ~ scale(DATE_NUM), random=~1|province, weights=count, data=data_agbyweekregion2)
# mblogit_fit2 = mblogit(LINEAGE ~ ns(DATE_NUM, df=2), random=~1|province, weights=count, data=data_agbyweekregion2)
# mblogit_fit3 = mblogit(LINEAGE ~ ns(DATE_NUM, df=3), random=~1|province, weights=count, data=data_agbyweekregion2)
  

# growth rate advantage compared to Omicron (BA.1) (difference in growth rate per day) 
emtrbelgium = emtrends(fit3_belgium_multi, trt.vs.ctrl ~ LINEAGE,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=seq(today_num, today_num, by=1)))
delta_r_belgium = data.frame(confint(emtrbelgium, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrbelgium$contrasts)$p.value)
delta_r_belgium
# contrast    estimate          SE df   asymp.LCL    asymp.UCL      p.value
# 1                           Other - Omicron (BA.1)  0.10885499 0.004193745 NA  0.10063540  0.117074578 1.143530e-14
# 2                           Alpha - Omicron (BA.1) -0.68265728 0.050369879 NA -0.78138043 -0.583934135 1.143530e-14
# 3                            Beta - Omicron (BA.1)  0.18468291 0.018608941 NA  0.14821005  0.221155761 5.151435e-14
# 4                           Gamma - Omicron (BA.1)  0.03794555 0.082949639 NA -0.12463276  0.200523852 9.917987e-01
# 5                           Delta - Omicron (BA.1)  0.06965813 0.005109798 NA  0.05964311  0.079673145 1.143530e-14
# 6                  Omicron (BA.2) - Omicron (BA.1)  0.07187356 0.003662782 NA  0.06469464  0.079052484 1.143530e-14
# 7                  Omicron (BA.3) - Omicron (BA.1) -0.05661763 0.019125829 NA -0.09410356 -0.019131690 2.986852e-02
# 8  Omicron (BA.1.1, with S:R346K) - Omicron (BA.1) -0.01163099 0.002892810 NA -0.01730079 -0.005961185 8.921202e-04
# 9                  Omicron (BA.4) - Omicron (BA.1)  0.15905371 0.007319791 NA  0.14470718  0.173400239 1.143530e-14
# 10                 Omicron (BA.5) - Omicron (BA.1)  0.21311718 0.008893112 NA  0.19568701  0.230547363 1.143530e-14

# pairwise growth rate difference (differences in growth rate per day) 
# based on multinomial 3 df spline model with province included as a main effect
emtrbelgium_pairw = emtrends(fit3_belgium_multi, pairwise ~ LINEAGE,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=seq(today_num, today_num, by=1)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#                                           contrast     estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1                           Omicron (BA.1) - Other -0.055341857 0.007326450 NA -0.069701436 -0.040982279 2.734357e-10
# 2                           Omicron (BA.1) - Alpha  0.663614519 0.040968706 NA  0.583317331  0.743911707 1.265654e-14
# 3                            Omicron (BA.1) - Beta  0.355824450 0.315372494 NA -0.262294280  0.973943180 9.883512e-01
# 4                           Omicron (BA.1) - Gamma  0.188282852 0.076755373 NA  0.037845085  0.338720619 3.414799e-01
# 5                           Omicron (BA.1) - Delta -0.022792344 0.007386238 NA -0.037269104 -0.008315584 8.366611e-02
# 6                  Omicron (BA.1) - Omicron (BA.2) -0.080677388 0.003381889 NA -0.087305768 -0.074049007 1.265654e-14
# 7                  Omicron (BA.1) - Omicron (BA.3) -0.050839813 0.025102785 NA -0.100040368 -0.001639257 6.317714e-01
# 8  Omicron (BA.1) - Omicron (BA.1.1, with S:R346K)  0.003353692 0.002598869 NA -0.001739998  0.008447381 9.690689e-01
# 9                  Omicron (BA.1) - Omicron (BA.4) -0.192680513 0.015582732 NA -0.223222107 -0.162138919 1.776357e-14
# 10                 Omicron (BA.1) - Omicron (BA.5) -0.206344376 0.039306677 NA -0.283384046 -0.129304705 2.926338e-05
# 11                                   Other - Alpha  0.718956376 0.040225877 NA  0.640115107  0.797797645 1.265654e-14
# 12                                    Other - Beta  0.411166308 0.315282551 NA -0.206776136  1.029108752 9.667311e-01
# 13                                   Other - Gamma  0.243624710 0.076339658 NA  0.094001730  0.393247689 6.287694e-02
# 14                                   Other - Delta  0.032549513 0.002660800 NA  0.027334441  0.037764586 1.887379e-14
# 15                          Other - Omicron (BA.2) -0.025335530 0.007918085 NA -0.040854692 -0.009816369 6.143212e-02
# 16                          Other - Omicron (BA.3)  0.004502045 0.026102567 NA -0.046658047  0.055662137 1.000000e+00
# 17          Other - Omicron (BA.1.1, with S:R346K)  0.058695549 0.007636251 NA  0.043728772  0.073662326 1.323133e-10
# 18                          Other - Omicron (BA.4) -0.137338655 0.017167725 NA -0.170986778 -0.103690532 2.347611e-11
# 19                          Other - Omicron (BA.5) -0.151002518 0.039954097 NA -0.229311109 -0.072693927 1.020614e-02
# 20                                    Alpha - Beta -0.307790068 0.314768087 NA -0.924724182  0.309144045 9.962229e-01
# 21                                   Alpha - Gamma -0.475331667 0.078054634 NA -0.628315939 -0.322347394 5.588734e-07
# 22                                   Alpha - Delta -0.686406863 0.040505306 NA -0.765795804 -0.607017921 1.265654e-14
# 23                          Alpha - Omicron (BA.2) -0.744291907 0.041076226 NA -0.824799830 -0.663783983 1.265654e-14
# 24                          Alpha - Omicron (BA.3) -0.714454331 0.048020223 NA -0.808572239 -0.620336424 1.265654e-14
# 25          Alpha - Omicron (BA.1.1, with S:R346K) -0.660260827 0.041023320 NA -0.740665056 -0.579856598 1.265654e-14
# 26                          Alpha - Omicron (BA.4) -0.856295031 0.043810625 NA -0.942162278 -0.770427784 1.265654e-14
# 27                          Alpha - Omicron (BA.5) -0.869958894 0.056752971 NA -0.981192674 -0.758725115 1.265654e-14
# 28                                    Beta - Gamma -0.167541598 0.322992009 NA -0.800594304  0.465511107 9.999867e-01
# 29                                    Beta - Delta -0.378616794 0.315306455 NA -0.996606091  0.239372503 9.814788e-01
# 30                           Beta - Omicron (BA.2) -0.436501838 0.315386674 NA -1.054648361  0.181644685 9.503932e-01
# 31                           Beta - Omicron (BA.3) -0.406664263 0.316365966 NA -1.026730162  0.213401637 9.698952e-01
# 32           Beta - Omicron (BA.1.1, with S:R346K) -0.352470759 0.315379727 NA -0.970603665  0.265662148 9.891623e-01
# 33                           Beta - Omicron (BA.4) -0.548504963 0.315758582 NA -1.167380412  0.070370486 8.136764e-01
# 34                           Beta - Omicron (BA.5) -0.562168826 0.317808746 NA -1.185062521  0.060724870 7.961173e-01
# 35                                   Gamma - Delta -0.211075196 0.076527444 NA -0.361066230 -0.061084163 1.856757e-01
# 36                          Gamma - Omicron (BA.2) -0.268960240 0.076812382 NA -0.419509743 -0.118410737 2.521690e-02
# 37                          Gamma - Omicron (BA.3) -0.239122665 0.080739310 NA -0.397368804 -0.080876525 1.150363e-01
# 38          Gamma - Omicron (BA.1.1, with S:R346K) -0.184929160 0.076784194 NA -0.335423415 -0.034434906 3.689936e-01
# 39                          Gamma - Omicron (BA.4) -0.380963365 0.078309469 NA -0.534447103 -0.227479626 1.571515e-04
# 40                          Gamma - Omicron (BA.5) -0.394627228 0.086219423 NA -0.563614192 -0.225640264 5.188562e-04
# 41                          Delta - Omicron (BA.2) -0.057885044 0.007984604 NA -0.073534581 -0.042235507 1.416174e-09
# 42                          Delta - Omicron (BA.3) -0.028047469 0.026122212 NA -0.079246063  0.023151126 9.920423e-01
# 43          Delta - Omicron (BA.1.1, with S:R346K)  0.026146036 0.007702808 NA  0.011048808  0.041243263 3.499471e-02
# 44                          Delta - Omicron (BA.4) -0.169888169 0.017198656 NA -0.203596914 -0.136179423 1.174616e-13
# 45                          Delta - Omicron (BA.5) -0.183552032 0.039967228 NA -0.261886358 -0.105217705 4.871459e-04
# 46                 Omicron (BA.2) - Omicron (BA.3)  0.029837575 0.025029274 NA -0.019218901  0.078894052 9.824367e-01
# 47 Omicron (BA.2) - Omicron (BA.1.1, with S:R346K)  0.084031080 0.002903835 NA  0.078339668  0.089722491 1.265654e-14
# 48                 Omicron (BA.2) - Omicron (BA.4) -0.112003125 0.015317251 NA -0.142024386 -0.081981864 1.011542e-09
# 49                 Omicron (BA.2) - Omicron (BA.5) -0.125666988 0.039178583 NA -0.202455599 -0.048878377 6.010595e-02
# 50 Omicron (BA.3) - Omicron (BA.1.1, with S:R346K)  0.054193504 0.025059007 NA  0.005078752  0.103308256 5.350792e-01
# 51                 Omicron (BA.3) - Omicron (BA.4) -0.141840700 0.029323597 NA -0.199313893 -0.084367507 1.767716e-04
# 52                 Omicron (BA.3) - Omicron (BA.5) -0.155504563 0.046483893 NA -0.246611320 -0.064397806 4.048364e-02
# 53 Omicron (BA.1.1, with S:R346K) - Omicron (BA.4) -0.196034204 0.015518776 NA -0.226450447 -0.165617962 1.487699e-14
# 54 Omicron (BA.1.1, with S:R346K) - Omicron (BA.5) -0.209698067 0.039273949 NA -0.286673592 -0.132722542 1.951839e-05
# 55                 Omicron (BA.4) - Omicron (BA.5) -0.013663863 0.042039907 NA -0.096060566  0.068732840 9.999999e-01


# pairwise growth rate difference (differences in growth rate per day) 
# based on plain multinomial model with province included as a main effect
# (here growth rate differences are not time dependent)
emtrbelgium_pairw1 = emtrends(fit1_belgium_multi, pairwise ~ LINEAGE,  
                             var="DATE_NUM",  mode="latent",
                             at=list(DATE_NUM=seq(today_num, today_num, by=1)))
delta_r_belgium_pairw1 = data.frame(confint(emtrbelgium_pairw1, 
                                           adjust="none", df=NA)$contrasts, 
                                   p.value=as.data.frame(emtrbelgium_pairw1$contrasts)$p.value)
delta_r_belgium_pairw1
#                                           contrast     estimate           SE df    asymp.LCL    asymp.UCL      p.value
# 1                           Omicron (BA.1) - Other  0.160756314 1.555898e-06 NA  0.160753265  0.160759364 2.409184e-14
# 2                           Omicron (BA.1) - Alpha  0.165408237 1.508837e-06 NA  0.165405279  0.165411194 2.409184e-14
# 3                            Omicron (BA.1) - Beta  0.173394372 2.365307e-06 NA  0.173389736  0.173399008 2.409184e-14
# 4                           Omicron (BA.1) - Gamma  0.161864160 1.987507e-06 NA  0.161860265  0.161868056 2.409184e-14
# 5                           Omicron (BA.1) - Delta  0.125807841 1.230642e-06 NA  0.125805429  0.125810253 2.409184e-14
# 6                  Omicron (BA.1) - Omicron (BA.2) -0.103408811 1.172682e-06 NA -0.103411109 -0.103406512 2.409184e-14
# 7                  Omicron (BA.1) - Omicron (BA.3) -0.051589063 1.585422e-05 NA -0.051620137 -0.051557990 2.409184e-14
# 8  Omicron (BA.1) - Omicron (BA.1.1, with S:R346K) -0.026540071 8.321634e-07 NA -0.026541702 -0.026538440 2.409184e-14
# 9                  Omicron (BA.1) - Omicron (BA.4) -0.246294081 9.187242e-06 NA -0.246312088 -0.246276074 2.409184e-14
# 10                 Omicron (BA.1) - Omicron (BA.5) -0.138039266 1.856269e-05 NA -0.138075648 -0.138002884 2.409184e-14
# 11                                   Other - Alpha  0.004651922 7.785134e-07 NA  0.004650397  0.004653448 2.409184e-14
# 12                                    Other - Beta  0.012638058 1.962838e-06 NA  0.012634211  0.012641905 2.409184e-14
# 13                                   Other - Gamma  0.001107846 1.527865e-06 NA  0.001104852  0.001110841 2.409184e-14
# 14                                   Other - Delta -0.034948473 9.534058e-07 NA -0.034950342 -0.034946605 2.409184e-14
# 15                          Other - Omicron (BA.2) -0.264165125 1.869910e-06 NA -0.264168790 -0.264161460 2.409184e-14
# 16                          Other - Omicron (BA.3) -0.212345377 1.592044e-05 NA -0.212376581 -0.212314174 2.409184e-14
# 17          Other - Omicron (BA.1.1, with S:R346K) -0.187296385 1.641456e-06 NA -0.187299602 -0.187293168 2.409184e-14
# 18                          Other - Omicron (BA.4) -0.407050395 9.301946e-06 NA -0.407068627 -0.407032164 2.409184e-14
# 19                          Other - Omicron (BA.5) -0.298795580 1.861971e-05 NA -0.298832074 -0.298759086 2.409184e-14
# 20                                    Alpha - Beta  0.007986136 1.885016e-06 NA  0.007982441  0.007989830 2.409184e-14
# 21                                   Alpha - Gamma -0.003544076 1.448181e-06 NA -0.003546915 -0.003541238 2.409184e-14
# 22                                   Alpha - Delta -0.039600396 8.743148e-07 NA -0.039602109 -0.039598682 2.409184e-14
# 23                          Alpha - Omicron (BA.2) -0.268817047 1.830936e-06 NA -0.268820636 -0.268813459 2.409184e-14
# 24                          Alpha - Omicron (BA.3) -0.216997300 1.591591e-05 NA -0.217028494 -0.216966105 2.409184e-14
# 25          Alpha - Omicron (BA.1.1, with S:R346K) -0.191948307 1.596916e-06 NA -0.191951437 -0.191945178 2.409184e-14
# 26                          Alpha - Omicron (BA.4) -0.411702318 9.294190e-06 NA -0.411720534 -0.411684101 2.409184e-14
# 27                          Alpha - Omicron (BA.5) -0.303447502 1.861584e-05 NA -0.303483989 -0.303411016 2.409184e-14
# 28                                    Beta - Gamma -0.011530212 2.322410e-06 NA -0.011534764 -0.011525660 2.409184e-14
# 29                                    Beta - Delta -0.047586532 2.020507e-06 NA -0.047590492 -0.047582571 2.409184e-14
# 30                           Beta - Omicron (BA.2) -0.276803183 2.582715e-06 NA -0.276808245 -0.276798121 2.409184e-14
# 31                           Beta - Omicron (BA.3) -0.224983436 1.601981e-05 NA -0.225014834 -0.224952037 2.409184e-14
# 32           Beta - Omicron (BA.1.1, with S:R346K) -0.199934443 2.422443e-06 NA -0.199939191 -0.199929695 2.409184e-14
# 33                           Beta - Omicron (BA.4) -0.419688453 9.471011e-06 NA -0.419707016 -0.419669891 2.409184e-14
# 34                           Beta - Omicron (BA.5) -0.311433638 1.870475e-05 NA -0.311470299 -0.311396978 2.409184e-14
# 35                                   Gamma - Delta -0.036056320 1.561488e-06 NA -0.036059380 -0.036053259 2.409184e-14
# 36                          Gamma - Omicron (BA.2) -0.265272971 2.241858e-06 NA -0.265277365 -0.265268577 2.409184e-14
# 37                          Gamma - Omicron (BA.3) -0.213453224 1.596840e-05 NA -0.213484521 -0.213421926 2.409184e-14
# 38          Gamma - Omicron (BA.1.1, with S:R346K) -0.188404231 2.055174e-06 NA -0.188408259 -0.188400203 2.409184e-14
# 39                          Gamma - Omicron (BA.4) -0.408158241 9.383793e-06 NA -0.408176633 -0.408139850 2.409184e-14
# 40                          Gamma - Omicron (BA.5) -0.299903426 1.866073e-05 NA -0.299940000 -0.299866852 2.409184e-14
# 41                          Delta - Omicron (BA.2) -0.229216652 1.609405e-06 NA -0.229219806 -0.229213497 2.409184e-14
# 42                          Delta - Omicron (BA.3) -0.177396904 1.589195e-05 NA -0.177428052 -0.177365756 2.409184e-14
# 43          Delta - Omicron (BA.1.1, with S:R346K) -0.152347912 1.337157e-06 NA -0.152350532 -0.152345291 2.409184e-14
# 44                          Delta - Omicron (BA.4) -0.372101922 9.253098e-06 NA -0.372120058 -0.372083786 2.409184e-14
# 45                          Delta - Omicron (BA.5) -0.263847107 1.859536e-05 NA -0.263883553 -0.263810660 2.409184e-14
# 46                 Omicron (BA.2) - Omicron (BA.3)  0.051819748 1.585514e-05 NA  0.051788672  0.051850823 2.409184e-14
# 47 Omicron (BA.2) - Omicron (BA.1.1, with S:R346K)  0.076868740 1.064368e-06 NA  0.076866654  0.076870826 2.409184e-14
# 48                 Omicron (BA.2) - Omicron (BA.4) -0.142885270 9.112874e-06 NA -0.142903131 -0.142867409 2.409184e-14
# 49                 Omicron (BA.2) - Omicron (BA.5) -0.034630455 1.852870e-05 NA -0.034666771 -0.034594139 2.409184e-14
# 50 Omicron (BA.3) - Omicron (BA.1.1, with S:R346K)  0.025048992 1.585011e-05 NA  0.025017927  0.025080058 2.409184e-14
# 51                 Omicron (BA.3) - Omicron (BA.4) -0.194705018 1.828596e-05 NA -0.194740858 -0.194669178 2.409184e-14
# 52                 Omicron (BA.3) - Omicron (BA.5) -0.086450202 2.438355e-05 NA -0.086497993 -0.086402412 2.409184e-14
# 53 Omicron (BA.1.1, with S:R346K) - Omicron (BA.4) -0.219754010 9.173918e-06 NA -0.219771991 -0.219736030 2.409184e-14
# 54 Omicron (BA.1.1, with S:R346K) - Omicron (BA.5) -0.111499195 1.855600e-05 NA -0.111535564 -0.111462826 2.409184e-14
# 55                 Omicron (BA.4) - Omicron (BA.5)  0.108254815 2.063603e-05 NA  0.108214369  0.108295261 2.409184e-14


# estimated proportion of different LINEAGES among new lab diagnoses in Belgium today
# based on multinomial spline model
today # "2022-05-16"
multinom_preds_today_avg = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg


# estimated proportion of different LINEAGES among new infections in Belgium today
today+7 # "2021-07-29"
multinom_preds_infections_today_avg = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num+7), 
                                              mode="prob", df=NA))
multinom_preds_infections_today_avg


# estimates proportion of lab diagnoses that is BA.4/5 by province
multinom_preds_today_byprovince = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1, by="province",
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_delta_today_byprovince_BA4 = multinom_preds_today_byprovince[multinom_preds_today_byprovince$LINEAGE=="Omicron (BA.4)",]
multinom_delta_today_byprovince_BA4 = multinom_delta_today_byprovince_BA4[order(multinom_delta_today_byprovince_BA4$prob, decreasing=T),]
multinom_delta_today_byprovince_BA4


multinom_delta_today_byprovince_BA5 = multinom_preds_today_byprovince[multinom_preds_today_byprovince$LINEAGE=="Omicron (BA.5)",]
multinom_delta_today_byprovince_BA5 = multinom_delta_today_byprovince_BA5[order(multinom_delta_today_byprovince_BA5$prob, decreasing=T),]
multinom_delta_today_byprovince_BA5



# reorder provinces by incidence of BA.4
# levels_PROVINCES = as.character(multinom_delta_today_byprovince$province[order(multinom_delta_today_byprovince$prob, decreasing=T)])
# data_agbyweekregion2$province = factor(data_agbyweekregion2$province, levels=levels_PROVINCES)
data_agbyweekregion2$LINEAGE = factor(data_agbyweekregion2$LINEAGE, levels=levels_LINEAGE_plot)

# redo plot of raw data by province using this order
muller_belgiumbyprovince_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE)) +
  facet_wrap(~ province, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + xlab("Collection date") +
  theme(legend.position="right") + # axis.title.x=element_blank()
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nGISAID data")
# +
# coord_cartesian(xlim=c(1,max(GISAID_belgium$Week)))
muller_belgiumbyprovince_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_raw data.pdf"), width=8, height=6)


# PLOT MULTINOMIAL FIT

extrapolate = 30
date.from = as.numeric(as.Date("2020-03-01"))
date.to = today_num+extrapolate # max(GISAID_belgium$DATE_NUM)+extrapolate

fit_belgium_multi_predsbyprovince = data.frame(emmeans(fit3_belgium_multi,
                                                  ~ LINEAGE,
                                                  by=c("DATE_NUM", "province"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                                  mode="prob", df=NA))
fit_belgium_multi_predsbyprovince$collection_date = as.Date(fit_belgium_multi_predsbyprovince$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_predsbyprovince$prob[fit_belgium_multi_predsbyprovince$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_predsbyprovince$LINEAGE!="Other"] = 0
fit_belgium_multi_predsbyprovince$prob[fit_belgium_multi_predsbyprovince$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_predsbyprovince$LINEAGE=="Other"] = 1
fit_belgium_multi_predsbyprovince$asymp.LCL[fit_belgium_multi_predsbyprovince$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_predsbyprovince$LINEAGE!="Other"] = 0
fit_belgium_multi_predsbyprovince$asymp.LCL[fit_belgium_multi_predsbyprovince$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_predsbyprovince$LINEAGE=="Other"] = 1
fit_belgium_multi_predsbyprovince$asymp.UCL[fit_belgium_multi_predsbyprovince$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_predsbyprovince$LINEAGE!="Other"] = 0
fit_belgium_multi_predsbyprovince$asymp.UCL[fit_belgium_multi_predsbyprovince$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_predsbyprovince$LINEAGE=="Other"] = 1
fit_belgium_multi_predsbyprovince$LINEAGE = factor(fit_belgium_multi_predsbyprovince$LINEAGE, levels=levels_LINEAGE_plot)
fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=levels_PROVINCES)

# multinomial model predictions overall for Belgium (with model without province) (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to)))
fit_belgium_multi_preds0 = data.frame(predgrid, as.data.frame(predict(fit0_belgium_multi, newdata=predgrid, type="prob")),check.names=F)
fit_belgium_multi_preds0 = gather(fit_belgium_multi_preds0, LINEAGE, prob, all_of(levels_LINEAGE), factor_key=TRUE)
fit_belgium_multi_preds0$collection_date = as.Date(fit_belgium_multi_preds0$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_preds0$prob[fit_belgium_multi_preds0$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds0$LINEAGE!="Other"] = 0
fit_belgium_multi_preds0$prob[fit_belgium_multi_preds0$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds0$LINEAGE=="Other"] = 1
fit_belgium_multi_preds0$LINEAGE = factor(fit_belgium_multi_preds0$LINEAGE, levels=levels_LINEAGE_plot)

# multinomial model predictions overall for Belgium, here with confidence intervals, 
# using multinomial 3 df spline model with province included
# fit looks OK
fit_belgium_multi_preds = data.frame(emmeans(fit3_belgium_multi, 
                                           ~ LINEAGE,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 to speed up things a bit
                                           mode="prob", df=NA))
fit_belgium_multi_preds$collection_date = as.Date(fit_belgium_multi_preds$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_preds$prob[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE!="Other"] = 0
fit_belgium_multi_preds$prob[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE=="Other"] = 1
fit_belgium_multi_preds$asymp.LCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE!="Other"] = 0
fit_belgium_multi_preds$asymp.LCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE=="Other"] = 1
fit_belgium_multi_preds$asymp.UCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE!="Other"] = 0
fit_belgium_multi_preds$asymp.UCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE=="Other"] = 1
fit_belgium_multi_preds$LINEAGE = factor(fit_belgium_multi_preds$LINEAGE, levels=levels_LINEAGE_plot) 

# multinomial model predictions overall for Belgium, here with confidence intervals, ####
# using multinomial 2 df spline model with province included
# fit looks OK - I WILL USE THIS MODEL
fit_belgium_multi_preds2 = data.frame(emmeans(fit2_belgium_multi, 
                                             ~ LINEAGE,
                                             by=c("DATE_NUM"),
                                             at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 to speed up things a bit
                                             mode="prob", df=NA))
fit_belgium_multi_preds2$collection_date = as.Date(fit_belgium_multi_preds2$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_preds2$prob[fit_belgium_multi_preds2$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds2$LINEAGE!="Other"] = 0
fit_belgium_multi_preds2$prob[fit_belgium_multi_preds2$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds2$LINEAGE=="Other"] = 1
fit_belgium_multi_preds2$asymp.LCL[fit_belgium_multi_preds2$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds2$LINEAGE!="Other"] = 0
fit_belgium_multi_preds2$asymp.LCL[fit_belgium_multi_preds2$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds2$LINEAGE=="Other"] = 1
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds2$LINEAGE!="Other"] = 0
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds2$LINEAGE=="Other"] = 1
fit_belgium_multi_preds2$LINEAGE = factor(fit_belgium_multi_preds2$LINEAGE, levels=levels_LINEAGE_plot)

# CALCULATE PAIRWISE GROWTH RATE ADVANTAGES OF EACH VARIANT AT TIMEPOINT WHERE IT FIRST EXCEEDED 1% OF THE CASES ####
# BASED ON MULTINOMIAL MULTINOMIAL 2 DF SPLINE FIT WITH PROVINCE INCLUDED AS MAIN EFFECT
# (PS should still be subsetted to LINEAGES THAT THEN CO-OCCURRED AT REASONABLE FREQUENCY)
library(dplyr)
# fit_belgium_multi_preds2_maxima = as.data.frame(fit_belgium_multi_preds2[fit_belgium_multi_preds2$collection_date<=today,] 
#                                                 %>% group_by(LINEAGE) %>% top_n(1, prob))
# fit_belgium_multi_preds2_maxima = fit_belgium_multi_preds2_maxima[fit_belgium_multi_preds2_maxima$LINEAGE!="Other",]
date_at1perc = as.data.frame(fit_belgium_multi_preds2[fit_belgium_multi_preds2$collection_date<=today,]
                                                %>% group_by(LINEAGE) %>% summarise(DATE_1PERCENT=first(collection_date[prob>=0.01])))
date_at1perc = date_at1perc[date_at1perc$LINEAGE!="Other"&!is.na(date_at1perc$DATE_1PERCENT),]

growth_advantages_at1perc = do.call(rbind, lapply(date_at1perc$LINEAGE,
       function (LINEAGE) { 
         t = as.numeric(date_at1perc[date_at1perc$LINEAGE==LINEAGE,"DATE_1PERCENT"]) 
         emtrbelgium_pw = emtrends(fit2_belgium_multi, pairwise ~ LINEAGE,  
                                       var="DATE_NUM",  mode="latent",
                                       at=list(DATE_NUM=t))
         delta_r_belgium_pw = data.frame(LINEAGE_AT1PERCENT = LINEAGE,
                                         DATE_NUM_AT1PERCENT=t, 
                                         DATE_AT1PERCENT=as.Date(t, origin="1970-01-01"),
                                         confint(emtrbelgium_pw, 
                                                     adjust="none", df=NA)$contrasts, 
                                         p.value=as.data.frame(emtrbelgium_pw$contrasts)$p.value) 
         return(delta_r_belgium_pw)
         }) )

write.csv(growth_advantages_at1perc, ".//data//GISAID//Belgium//pairwise daily growth rate advantage when lineages at 1 percent_belgium_multinomial 2 df spline fit with province.csv",
          row.names=F)


# PLOT OF MULTINOMIAL FIT ON LOGIT SCALE ####
fit_belgium_multi_pr = fit_belgium_multi_preds2
ymin = 0.001
ymax = 0.99
fit_belgium_multi_pr$asymp.LCL[fit_belgium_multi_pr$asymp.LCL<ymin] = ymin
fit_belgium_multi_pr$asymp.UCL[fit_belgium_multi_pr$asymp.UCL<ymin] = ymin
fit_belgium_multi_pr$asymp.UCL[fit_belgium_multi_pr$asymp.UCL>ymax] = ymax
fit_belgium_multi_pr$prob[fit_belgium_multi_pr$prob<ymin] = ymin
data_agbyweek2$prop[data_agbyweek2$collection_date<=as.Date("2021-09-01")&grepl("Omicron", data_agbyweek2$LINEAGE)] = 0 # fix above

plot_belgium_mfit_logit = qplot(data=fit_belgium_multi_pr, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial spline fit VARIANT ~ ns(DATE, df=2)+PROVINCE)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-03-01"),NA), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_belgium_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom 2 df spline fit with province_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_logit scale.pdf"), width=8, height=6)









# multinomial model predictions overall for Belgium, here with confidence intervals, 
# using plain multinomial model with province included
# fit not OK
fit_belgium_multi_preds1 = data.frame(emmeans(fit1_belgium_multi, 
                                             ~ LINEAGE,
                                             by=c("DATE_NUM"),
                                             at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 to speed up things a bit
                                             mode="prob", df=NA))
fit_belgium_multi_preds1$collection_date = as.Date(fit_belgium_multi_preds1$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_preds1$prob[fit_belgium_multi_preds1$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds1$LINEAGE!="Other"] = 0
fit_belgium_multi_preds1$prob[fit_belgium_multi_preds1$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds1$LINEAGE=="Other"] = 1
fit_belgium_multi_preds1$asymp.LCL[fit_belgium_multi_preds1$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds1$LINEAGE!="Other"] = 0
fit_belgium_multi_preds1$asymp.LCL[fit_belgium_multi_preds1$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds1$LINEAGE=="Other"] = 1
fit_belgium_multi_preds1$asymp.UCL[fit_belgium_multi_preds1$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds1$LINEAGE!="Other"] = 0
fit_belgium_multi_preds1$asymp.UCL[fit_belgium_multi_preds1$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds1$LINEAGE=="Other"] = 1
fit_belgium_multi_preds1$LINEAGE = factor(fit_belgium_multi_preds1$LINEAGE, levels=levels_LINEAGE_plot) 


muller_belgium_mfit = ggplot(data=fit_belgium_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE)) + 
  # facet_wrap(~ province, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  annotate("rect", xmin=max(GISAID_belgium$DATE_NUM)+1, 
           xmax=date.to, ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +  # , axis.title.x=element_blank()
  ylab("Share") + xlab("Collection date") + 
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)")
muller_belgium_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_multinom fit.png"), width=8, height=5)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_multinom fit.pdf"), width=8, height=5)


library(ggpubr)
ggarrange(muller_belgium_raw2+coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to,origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))) + theme(axis.title.x=element_blank()), 
          muller_belgium_mfit+ggtitle("Multinomial fit")+coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to,origin="1970-01-01"))), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots multipanel_multinom fit.png"), width=8, height=8)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots multipanel_multinom fit.pdf"), width=8, height=8)


muller_belgiumbyprovince_mfit = ggplot(data=fit_belgium_multi_predsbyprovince,
                                  aes(x=collection_date, y=prob, group=LINEAGE)) +
  facet_wrap(~ province, ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  annotate("rect", xmin=max(GISAID_belgium$DATE_NUM)+1, 
           xmax=date.to, ymin=0, ymax=1, alpha=0.2, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("Share") + xlab("Collection date") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)")
muller_belgiumbyprovince_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_multinom fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_multinom fit.pdf"), width=8, height=6)

ggarrange(muller_belgiumbyprovince_raw2 +
            coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to,origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0)))+
            ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nRaw GISAID data") + theme(axis.title.x=element_blank()),
          muller_belgiumbyprovince_mfit+ggtitle("\nMultinomial fit")+
            coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date(date.to,origin="1970-01-01"))), nrow=2)

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province multipanel_multinom fit.png"), width=8, height=12)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province multipanel_multinom fit.pdf"), width=8, height=12)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS ####


# MULTINOMIAL SPLINE FIT
# OVERALL FOR THE WHOLE OF BELGIUM 

# on logit scale:

fit_belgium_multi_preds2 = fit_belgium_multi_preds
ymin = 0.001
ymax = 0.999
fit_belgium_multi_preds2$asymp.LCL[fit_belgium_multi_preds2$asymp.LCL<ymin] = ymin
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$asymp.UCL<ymin] = ymin
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$asymp.UCL>ymax] = ymax
fit_belgium_multi_preds2$prob[fit_belgium_multi_preds2$prob<ymin] = ymin
data_agbyweek2$prop[data_agbyweek2$collection_date<=as.Date("2021-09-01")&grepl("Omicron", data_agbyweek2$LINEAGE)] = 0 # fix above

plot_belgium_mfit_logit = qplot(data=fit_belgium_multi_preds2, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-03-01"),NA), ylim=c(0.001, 0.999), expand=c(0,0))
plot_belgium_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom spline fit with province_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_belgium_mfit = qplot(data=fit_belgium_multi_preds, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share among new infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-03-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_belgium_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom spline fit with province_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_response scale.pdf"), width=8, height=6)


# MULTINOMIAL FIT
# OVERALL FOR THE WHOLE OF BELGIUM 

# on logit scale:

fit_belgium_multi_pr = fit_belgium_multi_preds2
ymin = 0.001
ymax = 0.999
fit_belgium_multi_pr$asymp.LCL[fit_belgium_multi_pr$asymp.LCL<ymin] = ymin
fit_belgium_multi_pr$asymp.UCL[fit_belgium_multi_pr$asymp.UCL<ymin] = ymin
fit_belgium_multi_pr$asymp.UCL[fit_belgium_multi_pr$asymp.UCL>ymax] = ymax
fit_belgium_multi_pr$prob[fit_belgium_multi_pr$prob<ymin] = ymin
data_agbyweek2$prop[data_agbyweek2$collection_date<=as.Date("2021-09-01")&grepl("Omicron", data_agbyweek2$LINEAGE)] = 0 # fix above

plot_belgium_mfit_logit = qplot(data=fit_belgium_multi_pr, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-03-01"),NA), ylim=c(0.001, 0.999), expand=c(0,0))
plot_belgium_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_belgium_mfit = qplot(data=fit_belgium_multi_preds, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share among new infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-03-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE
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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_response scale.pdf"), width=8, height=6)




# BY PROVINCE

# on logit scale:

fit_belgium_multi_preds3 = fit_belgium_multi_predsbyprovince
ymin = 0.001
ymax = 0.999
fit_belgium_multi_preds3$asymp.LCL[fit_belgium_multi_preds3$asymp.LCL<ymin] = ymin
fit_belgium_multi_preds3$asymp.UCL[fit_belgium_multi_preds3$asymp.UCL<ymin] = ymin
fit_belgium_multi_preds3$asymp.UCL[fit_belgium_multi_preds3$asymp.UCL>ymax] = ymax
fit_belgium_multi_preds3$prob[fit_belgium_multi_preds3$prob<ymin] = ymin

plot_belgium_mfit_logit = qplot(data=fit_belgium_multi_preds3, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~province) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweekregion2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-03-01"),NA), ylim=c(0.001, 0.999), expand=c(0,0))
plot_belgium_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit by province_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_belgium_mfit = qplot(data=fit_belgium_multi_preds, x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~province) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-03-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweekregion2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_belgium_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit by province_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_response scale.pdf"), width=8, height=6)




# project multinomial fit onto case data & infections as estimated by IHME ####

source("scripts/downloadData.R") # download latest data with new confirmed cases per day from Sciensano website, code adapted from https://github.com/JoFAM/covidBE_analysis by Joris Meys
tail(cases_tot)
tail(cases_prov)
# cases_tot = cases_tot[cases_tot$DATE<(as.Date(Sys.time())-3),] # we leave out last 3 days since they are incomplete
range(cases_tot$DATE) # "2020-03-01" "2022-04-27"

# also add smoothed case, growth in cases, smoothed tests, growth in tests, smoothed positivity ratio & growth in positivity ratio
library(mgcv)
k=50
fit_cases = gam(CASES ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F) +
                  WEEKDAY, # + # + offset(log(testspercase)),
                # BANKHOLIDAY,
                # + s(tests_new, bs="cs", k=8, fx=F),
                family=nb(), data=cases_tot,
                method = "REML",
                knots = list(DATE_NUM = c(min(cases_tot$DATE_NUM)-14,
                                          seq(min(cases_tot$DATE_NUM)+1*diff(range(cases_tot$DATE_NUM))/(k-2), 
                                              max(cases_tot$DATE_NUM)-1*diff(range(cases_tot$DATE_NUM))/(k-2), length.out=k-2),
                                          max(cases_tot$DATE_NUM)+14))
) 
BIC(fit_cases) # 9953.436

date.from = as.numeric(min(cases_tot$DATE)-14) # as.numeric(min(GISAID_sel$date_num))
date.to = today_num+extrapolate # max(GISAID_sel$date_num)+extrapolate

growth_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", 
                                      at=list(DATE_NUM = seq(date.from,
                                                             date.to, by=1)),
                                      type="link"))
colnames(growth_cases)[2] = "growth_cases"
colnames(growth_cases)[5] = "growth_cases_LOWER"
colnames(growth_cases)[6] = "growth_cases_UPPER"
growth_cases$DATE = as.Date(growth_cases$DATE_NUM, origin="1970-01-01")

cases_smoothed = as.data.frame(emmeans(fit_cases, ~ DATE_NUM, 
                                       at=list(DATE_NUM = seq(date.from,
                                                              date.to, by=1)),
                                       type="response"))
colnames(cases_smoothed)[2] = "cases_smoothed"
colnames(cases_smoothed)[5] = "cases_smootheds_LOWER"
colnames(cases_smoothed)[6] = "cases_smoothed_UPPER"
cases_smoothed$DATE = as.Date(cases_smoothed$DATE_NUM, origin="1970-01-01")


fit_tests = gam(TESTS_ALL ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F) +
                  WEEKDAY, # + # + offset(log(testspercase)),
                # BANKHOLIDAY,
                # + s(tests_new, bs="cs", k=8, fx=F),
                family=nb(), data=cases_tot,
                method = "REML",
                knots = list(DATE_NUM = c(min(cases_tot$DATE_NUM)-14,
                                          seq(min(cases_tot$DATE_NUM)+1*diff(range(cases_tot$DATE_NUM))/(k-2), 
                                              max(cases_tot$DATE_NUM)-1*diff(range(cases_tot$DATE_NUM))/(k-2), length.out=k-2),
                                          max(cases_tot$DATE_NUM)+14))
) 
BIC(fit_tests) # 9953.436

growth_tests = as.data.frame(emtrends(fit_tests, ~ DATE_NUM, var="DATE_NUM", 
                                      at=list(DATE_NUM = seq(date.from,
                                                             date.to, by=1)),
                                      type="link"))
colnames(growth_tests)[2] = "growth_tests"
colnames(growth_tests)[5] = "growth_tests_LOWER"
colnames(growth_tests)[6] = "growth_testss_UPPER"
growth_tests$DATE = as.Date(growth_tests$DATE_NUM, origin="1970-01-01")

tests_smoothed = as.data.frame(emmeans(fit_tests, ~ DATE_NUM, 
                                       at=list(DATE_NUM = seq(date.from,
                                                              date.to, by=1)),
                                       type="response"))
colnames(tests_smoothed)[2] = "tests_smoothed"
colnames(tests_smoothed)[5] = "tests_smoothed_LOWER"
colnames(tests_smoothed)[6] = "tests_smoothed_UPPER"
tests_smoothed$DATE = as.Date(tests_smoothed$DATE_NUM, origin="1970-01-01")

fit_posratio = gam(cbind(CASES, TESTS_ALL-CASES) ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F) +
                     WEEKDAY, # + # + offset(log(testspercase)),
                   # BANKHOLIDAY,
                   # + s(tests_new, bs="cs", k=8, fx=F),
                   family=binomial, data=cases_tot,
                   method = "REML",
                   knots = list(DATE_NUM = c(min(cases_tot$DATE_NUM)-14,
                                             seq(min(cases_tot$DATE_NUM)+1*diff(range(cases_tot$DATE_NUM))/(k-2), 
                                                 max(cases_tot$DATE_NUM)-1*diff(range(cases_tot$DATE_NUM))/(k-2), length.out=k-2),
                                             max(cases_tot$DATE_NUM)+14))
) 
BIC(fit_posratio) 

growth_posratio = as.data.frame(emtrends(fit_posratio, ~ DATE_NUM, var="DATE_NUM", 
                                         at=list(DATE_NUM = seq(date.from,
                                                                date.to, by=1)),
                                         type="link"))
colnames(growth_posratio)[2] = "growth_posratio"
colnames(growth_posratio)[5] = "growth_posratio_LOWER"
colnames(growth_posratio)[6] = "growth_posratio_UPPER"
growth_posratio$DATE = as.Date(growth_posratio$DATE_NUM, origin="1970-01-01")

posratio_smoothed = as.data.frame(emmeans(fit_posratio, ~ DATE_NUM, 
                                          at=list(DATE_NUM = seq(date.from,
                                                                 date.to, by=1)),
                                          type="response"))
colnames(posratio_smoothed)[2] = "posratio_smoothed"
colnames(posratio_smoothed)[5] = "posratio_smoothed_LOWER"
colnames(posratio_smoothed)[6] = "posratio_smoothed_UPPER"
posratio_smoothed$DATE = as.Date(posratio_smoothed$DATE_NUM, origin="1970-01-01")

cases_tot$growth_cases = growth_cases$growth_cases[match(cases_tot$DATE, growth_cases$DATE)]
cases_tot$cases_smoothed = cases_smoothed$cases_smoothed[match(cases_tot$DATE, cases_smoothed$DATE)]
cases_tot$growth_tests = growth_tests$growth_tests[match(cases_tot$DATE, growth_tests$DATE)]
cases_tot$tests_smoothed = tests_smoothed$tests_smoothed[match(cases_tot$DATE, tests_smoothed$DATE)]
cases_tot$growth_posratio = growth_posratio$growth_posratio[match(cases_tot$DATE, growth_posratio$DATE)]
cases_tot$posratio_smoothed = posratio_smoothed$posratio_smoothed[match(cases_tot$DATE, posratio_smoothed$DATE)]

# estimated true infections (IHME)
# (https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00484-6/fulltext?utm_campaign=lancetcovid22&utm_source=twitter&utm_medium=social#sec1, https://github.com/ihmeuw/covid-historical-model, https://github.com/ihmeuw/covid-model-infections, downloaded from https://ourworldindata.org/grapher/daily-new-estimated-covid-19-infections-ihme-model?country=~BEL)
infections_IHME = read.csv(".//data//daily-new-estimated-covid-19-infections-ihme-model.csv")
colnames(infections_IHME) = c("COUNTRY","ISO3","DATE","infections_IHME","infections_IHME_lower","infections_IHME_upper","infections_IHME_rolling_7d_mean")
infections_IHME$DATE = as.Date(infections_IHME$DATE)
infections_IHME$DATE_NUM = as.numeric(infections_IHME$DATE)
infections_IHME_BE = infections_IHME[infections_IHME$ISO3=="BEL",]
cases_tot$infections_IHME = infections_IHME_BE$infections_IHME[match(cases_tot$DATE, infections_IHME$DATE)]

# estimate new infections from cases & pos ratio
library(Hmisc)
library(splines)
fit_infections = glm(Lag(infections_IHME, shift=8) ~ ns(DATE_NUM, df=3) * log(cases_smoothed), # + 
                     # log(tests_smoothed),
                     # ns(log(tests_smoothed), df=2), 
                     family=poisson, data=cases_tot[!is.na(cases_tot$infections_IHME),]) # good model
# fit_infections = glm(Lag(infections_IHME, shift=8) ~ ns(date_num, df=3) + 
#                        log(cases_smoothed) + 
#                        ns(growth_posratio, df=2),  
#                      family=poisson, data=cases_tot[!is.na(cases_tot$infections_IHME),])


plot(model.frame(fit_infections)[,1], col="black", type="l")
lines(fitted(fit_infections), col="red")
sum(fitted(fit_infections)) # 10 549 856, BE population is 11,681,636 as of Sunday, May 1, 2022, based on Worldometer elaboration of the latest United Nations data
sum(cases_tot$infections_IHME, na.rm=T) # 10 915 123

plot(predict(fit_infections, newdata=cases_tot, type="response"), type="l", col="red")
sum(predict(fit_infections, newdata=cases_tot, type="response"), na.rm=T) # 12 785 249

cases_tot$est_infections = predict(fit_infections, newdata=cases_tot, type="response")
cases_tot$est_infections[is.na(cases_tot$est_infections)] = 0
cases_tot$est_infections_cumsum = cumsum(cases_tot$est_infections)
populationBE = 10549856
cases_tot$est_infections_cumsum_perpop = 100*cases_tot$est_infections_cumsum/populationBE
plot(cases_tot$est_infections_cumsum_perpop, type="l", ylab="Cumulative infections (IHME & Sciensano case data) (% of population)")

# also add variable with estimated excess deaths
# TO DO : plot infection fatality rate over time based on IHME true infection estimates & excess mortality ests of the Economist
excess_mort = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/export_country_cumulative.csv")
excess_mort$date = as.Date(excess_mort$date)
excess_mort_BE = excess_mort[excess_mort$iso3c=="BEL",]
# TO DO: interpolate, as data now only given per week
cases_tot$cumulative_estimated_daily_excess_deaths_7d_later = excess_mort_BE$cumulative_estimated_daily_excess_deaths[match(cases_tot$DATE+7,excess_mort_BE$date)]      
cases_tot$infection_fatality_rate = 100*cases_tot$cumulative_estimated_daily_excess_deaths_7d_later/cases_tot$est_infections_cumsum
plot(x=cases_tot$DATE, y=cases_tot$infection_fatality_rate, type="p", ylim=c(0,1))
sort(unique(excess_mort$iso3c))
plot(x=excess_mort[excess_mort$iso3c=="BEL",]$date,
     y=excess_mort[excess_mort$iso3c=="BEL",]$cumulative_estimated_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="BEL",]$date,
     y=excess_mort[excess_mort$iso3c=="BEL",]$cumulative_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="BEL",]$date,
     y=excess_mort[excess_mort$iso3c=="BEL",]$cumulative_daily_covid_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="ZAF",]$date,
     y=excess_mort[excess_mort$iso3c=="ZAF",]$cumulative_estimated_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="ZAF",]$date,
     y=excess_mort[excess_mort$iso3c=="ZAF",]$cumulative_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="ZAF",]$date,
     y=excess_mort[excess_mort$iso3c=="ZAF",]$cumulative_daily_covid_deaths,
     type="l")

# plot of cum in excess mortality in Belgium
ggplot(data=excess_mort[excess_mort$iso3c=="BEL",], 
       aes(x=date, y=100*cumulative_estimated_daily_excess_deaths/populationBE)) + 
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=I("red3"), width=I(1.1))) +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("Cumulative exces mortality (% of population)") + xlab("Date") +
  ggtitle("CUMULATIVE EXCESS MORTALITY IN BELGIUM",
          "(estimates by The Economist)") +
  # scale_fill_manual("variant", values=lineage_cols_plot) +
  # scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(0,0.5))
# coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cumulative excess mortality Belgium.png"), width=8, height=6)


# plot of cum in excess mortality in South Africa
ggplot(data=excess_mort[excess_mort$iso3c=="ZAF",], 
       aes(x=date, y=100*cumulative_estimated_daily_excess_deaths/populationSA)) + 
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=I("red3"), width=I(1.1))) +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("Cumulative exces mortality (% of population)") + xlab("Date") +
  ggtitle("CUMULATIVE EXCESS MORTALITY IN SOUTH AFRICA",
          "(estimates by The Economist)") +
  # scale_fill_manual("variant", values=lineage_cols_plot) +
  # scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(0,0.5))
# coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cumulative excess mortality South Africa.png"), width=8, height=6)




# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####


# model with correction for variable testing intensity
fit_cases_testing = gam(CASES ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F) +
                          WEEKDAY + # + # + offset(log(testspercase)),
                          # BANKHOLIDAY,
                          s(TESTS_ALL, bs="cs", k=8, fx=F),
                        family=nb(), data=cases_tot,
                        method = "REML",
                        knots = list(date_num = c(min(cases_tot$DATE_NUM)-14,
                                                  seq(min(cases_tot$DATE_NUM)+0.75*diff(range(cases_tot$DATE_NUM))/(k-2), 
                                                      max(cases_tot$DATE_NUM)-0.75*diff(range(cases_tot$DATE_NUM))/(k-2), length.out=k-2),
                                                  max(cases_tot$DATE_NUM)+14))
)

                        
fit_belgium_multi_preds0$totcases = cases_tot$CASES[match(round(fit_belgium_multi_preds0$DATE_NUM),cases_tot$DATE_NUM)]
fit_belgium_multi_preds0$cases = fit_belgium_multi_preds0$totcases * fit_belgium_multi_preds0$prob
fit_belgium_multi_preds0$cases[fit_belgium_multi_preds0$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases_testing, ~ DATE_NUM, at=list(DATE_NUM=seq(date.from, date.to, by=1),
                                                                             testspercase=20,
                                                                             TESTS_ALL=max(cases_tot$TESTS_ALL, na.rm=T)), 
                                      type="response"))
fit_belgium_multi_preds0$smoothed_totcases = cases_emmeans$response[match(fit_belgium_multi_preds0$DATE_NUM,cases_emmeans$DATE_NUM)]
fit_belgium_multi_preds0$smoothed_cases = fit_belgium_multi_preds0$smoothed_totcases * fit_belgium_multi_preds0$prob
fit_belgium_multi_preds0$smoothed_cases[fit_belgium_multi_preds0$smoothed_cases<=0.001] = NA
fit_belgium_multi_preds0$LINEAGE = factor(fit_belgium_multi_preds0$LINEAGE, levels=levels_LINEAGE_plot)

qplot(data=fit_belgium_multi_preds0, x=collection_date, y=cases, geom="line", colour=LINEAGE) +
  xlim(c(as.Date("2021-10-01"), today)) + scale_y_log10() + 
  geom_line(aes(y=smoothed_cases), lwd=I(2)) +
  ylab("new Omicron cases per day") +
  scale_colour_manual(values=lineage_cols_plot) +
  ggtitle("Inferred nr. of new Omicron cases per day in Belgium", "(black=observed,\nred=GAM negative binomial spline smooth with correction for\nvariable testing intensity & weekday effects & marginal means\ncalculated at uniform, maximal testing effort)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\variant cases observed and fitted.png"), width=8, height=6)


ggplot(data=fit_belgium_multi_preds0, 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE, width=I(1.1)), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN BELGIUM","(case data Sciensano plus multinomial spline fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
# coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)
write.csv(fit_belgium_multi_preds0, file=paste0(".\\plots\\",plotdir,"\\cases per day by variant Belgium.csv"), row.names=F)

ggplot(data=fit_belgium_multi_preds0[fit_belgium_multi_preds0$collection_date<=today,],  # or (today-1) or (today-2) to be on the safe side
       aes(x=collection_date, y=smoothed_cases, group=LINEAGE)) + 
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN BELGIUM","(negative binomial fit to Sciensano case data, with correction for weekday effects &\nvariable testing intensity; marginal means calculated at constant, maximal testing effort\nvariant frequencies based on multinomial spline fit to GISAID data)") + #  
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit case data.png"), width=8, height=6)



# STACKED AREA CHART OF NEW INFECTIONS BY VARIANT (MULTINOMIAL FIT MAPPED ONTO ESTIMATED INFECTIONS) ####

fit_belgium_multi_preds0$totinfections = cases_tot$est_infections[match(round(fit_belgium_multi_preds0$DATE_NUM),cases_tot$DATE_NUM)]
fit_belgium_multi_preds0$totinfections[is.na(fit_belgium_multi_preds$totinfections)] = 0
fit_belgium_multi_preds0$infections = fit_belgium_multi_preds0$totinfections * fit_belgium_multi_preds0$prob
fit_belgium_multi_preds0$infections[fit_belgium_multi_preds0$infections<=0.001] = 0

fit_belgium_multi_preds0$LINEAGE = factor(fit_belgium_multi_preds0$LINEAGE, levels=levels_LINEAGE_plot)

ggplot(data=fit_belgium_multi_preds0, 
       aes(x=collection_date, y=infections*100/populationBE, group=LINBEAGE)) + 
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE, width=I(1.1)), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("Estimated infections per day (% of population)") + xlab("Date of diagnosis") +
  ggtitle("ESTIMATED SARS-CoV2 INFECTIONS PER DAY BY VARIANT IN BELGIUM","(case data Sciensano mapped to true number of daily infections as estimated by IHME using model\nglm(Lag(infections_IHME, shift=8) ~ ns(date, df=3)*log(cases_smoothed), family=poisson(log)),\ncombined with multinomial 3 df spline fit to GISAID lineage frequency data)") +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
# coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\estimated infections IHME per day by variant_stacked area multinomial fit.png"), width=8, height=6)

# plot of cumulative nr of infections by variant
library(dplyr)
fit_belgium_multi_preds0 = fit_belgium_multi_preds0 %>% 
  group_by(LINEAGE) %>% 
  arrange(collection_date) %>% 
  mutate(cuminfections = cumsum(ifelse(is.na(infections), 0, infections)) + infections*0,
         cuminfectionsperpop = 100*(cumsum(ifelse(is.na(infections), 0, infections)) + infections*0)/populationBE)

ggplot(data=fit_belgium_multi_preds0[fit_belgium_multi_preds0$collection_date<=(today-9),], 
       aes(x=collection_date, y=cuminfectionsperpop, group=LINEAGE)) + 
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE, width=I(1.1)), 
           position = position_stack(reverse = T)) +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("Cumulative infections (% of population)") + xlab("Date of diagnosis") +
  ggtitle("CUMULATIVE SARS-CoV2 INFECTIONS BY VARIANT IN BELGIUM","(case data Sciensano mapped to true number of daily infections as estimated by IHME using model\nglm(Lag(infections_IHME, shift=8) ~ ns(date, df=3)*log(cases_smoothed), family=poisson(log)),\ncombined with multinomial 3 df spline fit to GISAID lineage frequency data)") +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(0,150))
# coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cumulative estimated infections IHME by variant_stacked area multinomial fit raw case data.png"), width=8, height=6)
write.csv(fit_belgium_multi_preds, file=paste0(".\\plots\\",plotdir,"\\cases infections per day by variant South Africa 5 Jan 2022.csv"), row.names=F)


