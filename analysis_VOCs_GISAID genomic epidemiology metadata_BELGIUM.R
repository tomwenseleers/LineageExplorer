# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN BELGIUM BASED ON ANALYSIS OF GISAID PATIENT METADATA
# T. Wenseleers
# last update 14 April 2022

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
d16 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_04_13_10_belgium_subm_apr_2022.tsv", col_types = cols(.default = "c"))

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
GISAID_belgium1 = rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16)
GISAID_belgium1 = GISAID_belgium1[GISAID_belgium1[,"Host"]=="Human",]
GISAID_belgium1[,"Collection date"] = as.Date(GISAID_belgium1[,"Collection date"])
GISAID_belgium1 = GISAID_belgium1[!is.na(GISAID_belgium1[,"Collection date"]),]
nrow(GISAID_belgium1) # 95418
unique(GISAID_belgium1[,"Virus name"])
unique(GISAID_belgium1[,"Location"])
unique(GISAID_belgium1[,"Host"])
unique(GISAID_belgium1[,"Additional location information"]) # ZIP codes # sometimes has additional info on whether it was a traveller
range(GISAID_belgium1$`Collection date`)

# anonymous records without location info are assumed to be active surveillance of cluster outbreaks & included in active surveillance (cf weekly Sciensano report)
GISAID_belgium1[,"Additional location information"][grepl("nknown",GISAID_belgium1[,"Additional location information"])] = NA
sum(is.na(GISAID_belgium1[,"Additional location information"])) # 8718
unique(GISAID_belgium1[,"Additional location information"])

library(stringr)
ZIP = str_extract(GISAID_belgium1[,"Additional location information"],"[0-9]{4}") # extract ZIP codes
unique(ZIP)
sum(!is.na(ZIP)) # 85672
sum(is.na(ZIP)) # 9746 - anonymous records, these are assumed to be active surveillance of cluster outbreaks

muni = read.csv(".//data//GISAID//Belgium//mapping_municips.csv") # post codes & cities & provinces
province = muni$province_label_nl[match(ZIP, muni$postcode)] # province
city = muni$municipality_label_nl[match(ZIP, muni$postcode)] # city
sum(is.na(city)) # 9869
sum(is.na(ZIP))  # 9746
wrongZIP = which(is.na(city)&!is.na(ZIP))
ZIP[wrongZIP] = NA

anonym = grepl("Belgium$",GISAID_belgium1[,"Additional location information"])|
  is.na(GISAID_belgium1[,"Additional location information"])|
  grepl("CHU|infected|travel|Travel",GISAID_belgium1[,"Additional location information"]) |
  is.na(ZIP)
sum(anonym) # 10087
unique(GISAID_belgium1[!is.na(GISAID_belgium1[,"Sampling strategy"]),"Sampling strategy"])
traveller = grepl("travel|Travel|infected",GISAID_belgium1[,"Additional location information"])|grepl("travel|Travel",GISAID_belgium1[,"Sampling strategy"])|
  grepl("travel|Travel|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(traveller) # 2816
Sdropout = grepl("S-gene dropout",GISAID_belgium1[,"Sampling strategy"])|
  grepl("S-gene dropout",GISAID_belgium1[,"Additional host information"])  
sum(Sdropout) # 21 - should be much higher though...
actsurv = anonym|traveller|grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster",GISAID_belgium1[,"Sampling strategy"])|
  grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster|HR contact|Outbreak|Nurse|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(actsurv) # 22150
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")]) # 19449 # In Sciensano weekly report of 21st of May 2021 this is 4710
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")&GISAID_belgium1[,"Collection date"]<=as.Date("2021-01-31")]) # 1204 # In Sciensano weekly report of 21st of May 2021 this is 1975

reinfection = grepl("Breakthrough",GISAID_belgium1[,"Sampling strategy"])
sum(reinfection) # 464
infectionpostvaccination = grepl("Vaccination|vaccination",GISAID_belgium1[,"Sampling strategy"])
sum(infectionpostvaccination) # 4971

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
sum(bassurv) # 73268
purpose_of_sequencing = rep("active_surveillance", nrow(GISAID_belgium1))
purpose_of_sequencing[bassurv] = "baseline_surveillance"

accessionIDs_bassurv = GISAID_belgium1[,"Accession ID"][bassurv]

GISAID_belgium1[,"Gender"][grepl("F|V",GISAID_belgium1[,"Gender"])] = "F"
GISAID_belgium1[,"Gender"][grepl("M|m",GISAID_belgium1[,"Gender"])] = "M"
GISAID_belgium1[,"Gender"][grepl("unknown",GISAID_belgium1[,"Gender"])] = NA
GISAID_belgium1[,"Gender"][!grepl("F|M",GISAID_belgium1[,"Gender"])] = NA
unique(GISAID_belgium1[,"Gender"]) # "F" "M" NA 

sum(!is.na(GISAID_belgium1[,"Gender"] )) # 84102 with gender info

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
sum(!is.na(GISAID_belgium1[,"Patient age"] )) # 83690 with age info

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
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2022-04-11_23-09.tsv.gz"), col_types = cols(.default = "c")) 
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
write.csv(GISAID_belgium1, ".//data//GISAID//Belgium//gisaid_hcov-19_2022_04_13_10_ALL PARSED.csv",row.names=F)
range(GISAID_belgium1$date)



# ANALYSIS OF VOC LINEAGE FREQUENCIES IN BELGIUM USING MULTINOMIAL FITS ####


library(dplyr)
GISAID_belgium1$LINEAGE = case_when(
  grepl("B.1.1.529|BA.1", GISAID_belgium1$pango_lineage) ~ "Omicron (BA.1*)", # B.1.1.529|
  # (GISAID_belgium1$pango_lineage=="BA.3") ~ "Omicron (BA.3)",
  # ((grepl("BA.2",GISAID_belgium1$pango_lineage))&
  #    ((grepl("L452R", GISAID_belgium1$AA_substitutions)&
  #        grepl("486V", GISAID_belgium1$AA_substitutions)&
  #        grepl("11F", GISAID_belgium1$AA_substitutions)&
  #        (!grepl("D3N",GISAID_belgium1$AA_substitutions)) ))) ~ "Omicron (BA.4)",
  # ((grepl("BA.2",GISAID_belgium1$pango_lineage))&
  #    ((grepl("L452R",GISAID_belgium1$AA_substitutions)&
  #        grepl("486V",GISAID_belgium1$AA_substitutions)&
  #        (!grepl("11F", GISAID_belgium1$AA_substitutions))&
  #        grepl("D3N",GISAID_belgium1$AA_substitutions)))) ~ "Omicron (BA.5)",
  (grepl("BA.2",GISAID_belgium1$pango_lineage)) ~ "Omicron (BA.2)",
  grepl("B.1.617.2", GISAID_belgium1$pango_lineage, fixed=T) | grepl("AY", GISAID_belgium1$pango_lineage)  ~ "Delta",
  grepl("B.1.1.7", GISAID_belgium1$pango_lineage, fixed=T) ~ "Alpha",
  grepl("B.1.351", GISAID_belgium1$pango_lineage, fixed=T) ~ "Beta",
  grepl("P.1", GISAID_belgium1$pango_lineage, fixed=T) ~ "Gamma",
  T ~ "Other"
)

GISAID_belgium1$LINEAGE[GISAID_belgium1$LINEAGE=="Alpha"&GISAID_belgium1$date<=as.Date("2020-10-01")] = "Other"
table(GISAID_belgium1$LINEAGE)

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
                  "Omicron (BA.1*)", "Omicron (BA.2)") 
levels_LINEAGE = c(main_lineages, "Other")
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
# lineage_cols[which(levels_VARIANTS=="C.1.2")] = "darkorange"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.1*)")] = "red" # "magenta"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.2)")] = "red3"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.3)")] = "red4" 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.4)")] = "darkorange" 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.5)")] = "darkorange3" 
lineage_cols[which(levels_LINEAGE=="Other")] = "grey65"


unique(GISAID_belgium1$LINEAGE)
GISAID_belgium1$LINEAGE[!(GISAID_belgium1$LINEAGE %in% main_lineages)] = "Other" # minority lineages & non-VOCs
table(GISAID_belgium1$LINEAGE)
GISAID_belgium1$LINEAGE = factor(GISAID_belgium1$LINEAGE, 
                                  levels=levels_LINEAGE, 
                                  labels=levels_LINEAGE)
table(GISAID_belgium1$LINEAGE)
# Alpha            Beta           Gamma           Delta Omicron (BA.1*)  Omicron (BA.2)           Other 
# 18463            1047            1725           43418           16348            7884            8529 
sum(table(GISAID_belgium1$LINEAGE)) # 95418

# GISAID_belgium = GISAID_sel[GISAID_sel$country=="Belgium",]
GISAID_belgium = GISAID_belgium1 # [GISAID_belgium1$purpose_of_sequencing=="baseline_surveillance",] # if desired subset to baseline surveillance
# GISAID_belgium = GISAID_belgium[-which((GISAID_belgium$LINEAGE=="B.1.617.2")&(GISAID_belgium$date<="2021-05-01")),]
nrow(GISAID_belgium) # 95418
sum(!is.na(GISAID_belgium$age)) # 83690 with age info

# use data from Feb 1 onwards
GISAID_belgium_all = GISAID_belgium
# GISAID_belgium = GISAID_belgium_all[GISAID_belgium_all$date>=as.Date("2021-02-01"),]
nrow(GISAID_belgium) # 89143
range(GISAID_belgium$date) # "2021-02-01" "2022-01-25"


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

write.csv(data_agbyweek2, ".//data//GISAID//Belgium//gisaid_hcov-19_2022_04_13_10_BASELINE SELECTION_aggregated counts by week.csv",row.names=F)


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

muller_belgium_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE)) + 
  # facet_wrap(~ province, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
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


muller_belgiumbyprovince_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE)) +
  facet_wrap(~ province, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
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
data_agbyweekregion2$LINEAGE = data_agbyweekregion2$LINEAGE
data_agbyweekregion2$LINEAGE = relevel(data_agbyweekregion2$LINEAGE, ref="Omicron (BA.1*)") 
fit1_belgium_multi = nnet::multinom(LINEAGE ~ province + DATE_NUM, weights=count, data=data_agbyweekregion2, maxit=1000)
fit2_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
fit3_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
fit4_belgium_multi = nnet::multinom(LINEAGE ~ province * DATE_NUM, weights=count, data=data_agbyweekregion2, maxit=1000)
fit5_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
fit6_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
BIC(fit1_belgium_multi, fit2_belgium_multi, fit3_belgium_multi, fit4_belgium_multi, fit5_belgium_multi, fit6_belgium_multi) 
# fit3_belgium_multi fits best (lowest BIC)

# growth rate advantage compared to Omicron (BA.1*) (difference in growth rate per day) 
emtrbelgium = emtrends(fit3_belgium_multi, trt.vs.ctrl ~ LINEAGE,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium = data.frame(confint(emtrbelgium, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrbelgium$contrasts)$p.value)
delta_r_belgium
#                             contrast    estimate          SE df   asymp.LCL   asymp.UCL      p.value
# 1          Alpha - (Omicron (BA.1*))  0.12988502 0.010838887 NA  0.10864119  0.15112885 1.318404e-10
# 2           Beta - (Omicron (BA.1*)) -0.78357420 0.030641797 NA -0.84363102 -0.72351739 1.318344e-10
# 3          Gamma - (Omicron (BA.1*))  0.14340872 0.016450397 NA  0.11116654  0.17565091 1.331841e-10
# 4          Delta - (Omicron (BA.1*)) -0.05117165 0.008293990 NA -0.06742757 -0.03491573 1.380875e-07
# 5 Omicron (BA.2) - (Omicron (BA.1*))  0.09099401 0.002612576 NA  0.08587346  0.09611457 1.318344e-10
# 6          Other - (Omicron (BA.1*))  0.08397853 0.004807300 NA  0.07455639  0.09340066 1.318344e-10

# pairwise growth rate difference (differences in growth rate per day) 
emtrbelgium_pairw = emtrends(fit3_belgium_multi, pairwise ~ LINEAGE,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#                              contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1           (Omicron (BA.1*)) - Alpha -0.129885018 0.010838887 NA -0.151128847 -0.10864119 1.582136e-10
# 2            (Omicron (BA.1*)) - Beta  0.783574205 0.030641797 NA  0.723517386  0.84363102 1.582013e-10
# 3           (Omicron (BA.1*)) - Gamma -0.143408725 0.016450397 NA -0.175650911 -0.11116654 1.628744e-10
# 4           (Omicron (BA.1*)) - Delta  0.051171653 0.008293990 NA  0.034915731  0.06742757 4.782837e-07
# 5  (Omicron (BA.1*)) - Omicron (BA.2) -0.090994011 0.002612576 NA -0.096114566 -0.08587346 1.582013e-10
# 6           (Omicron (BA.1*)) - Other -0.083978528 0.004807300 NA -0.093400662 -0.07455639 1.582013e-10
# 7                        Alpha - Beta  0.913459223 0.030883830 NA  0.852928028  0.97399042 1.582013e-10
# 8                       Alpha - Gamma -0.013523707 0.018083006 NA -0.048965747  0.02191833 9.890034e-01
# 9                       Alpha - Delta  0.181056671 0.011126435 NA  0.159249258  0.20286408 1.582013e-10
# 10             Alpha - Omicron (BA.2)  0.038891007 0.010999129 NA  0.017333110  0.06044890 1.138742e-02
# 11                      Alpha - Other  0.045906490 0.010102509 NA  0.026105936  0.06570704 3.576883e-04
# 12                       Beta - Gamma -0.926982929 0.033163491 NA -0.991982178 -0.86198368 1.582013e-10
# 13                       Beta - Delta -0.732402552 0.029727109 NA -0.790666615 -0.67413849 1.582013e-10
# 14              Beta - Omicron (BA.2) -0.874568215 0.030597230 NA -0.934537684 -0.81459875 1.582013e-10
# 15                       Beta - Other -0.867552733 0.029908246 NA -0.926171817 -0.80893365 1.582013e-10
# 16                      Gamma - Delta  0.194580377 0.016667150 NA  0.161913364  0.22724739 1.582209e-10
# 17             Gamma - Omicron (BA.2)  0.052414714 0.016553483 NA  0.019970484  0.08485894 3.370469e-02
# 18                      Gamma - Other  0.059430197 0.016009743 NA  0.028051676  0.09080872 6.529003e-03
# 19             Delta - Omicron (BA.2) -0.142165664 0.008554082 NA -0.158931357 -0.12539997 1.582013e-10
# 20                      Delta - Other -0.135150181 0.005662904 NA -0.146249269 -0.12405109 1.582013e-10
# 21             Omicron (BA.2) - Other  0.007015483 0.005196824 NA -0.003170104  0.01720107 8.261329e-01


# estimated proportion of different LINEAGES among new lab diagnoses in Belgium today
today # "2022-04-13"
multinom_preds_today_avg = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1*) 1.669509e-02 1.669078e-03 NA  1.342376e-02 1.996643e-02
# 2           Alpha 3.075379e-04 4.101230e-04 NA -4.962884e-04 1.111364e-03
# 3            Beta 4.418289e-66 1.607583e-65 NA -2.708976e-65 3.592634e-65
# 4           Gamma 1.902981e-04 3.994613e-04 NA -5.926316e-04 9.732277e-04
# 5           Delta 9.428180e-07 5.760644e-07 NA -1.862474e-07 2.071883e-06
# 6  Omicron (BA.2) 9.785669e-01 2.323917e-03 NA  9.740121e-01 9.831217e-01
# 7           Other 4.239193e-03 1.176099e-03 NA  1.934081e-03 6.544305e-03


# estimated proportion of different LINEAGES among new infections in Belgium today
today+7 # "2021-07-29"
multinom_preds_infections_today_avg = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num+7), 
                                              mode="prob", df=NA))
multinom_preds_infections_today_avg
#           LINEAGE         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1*) 8.907198e-03 1.047559e-03 NA  6.854020e-03 1.096038e-02
# 2           Alpha 4.070824e-04 5.705654e-04 NA -7.112052e-04 1.525370e-03
# 3            Beta 9.760625e-69 3.732763e-68 NA -6.340019e-68 8.292144e-68
# 4           Gamma 2.772463e-04 6.109594e-04 NA -9.202121e-04 1.474705e-03
# 5           Delta 3.514871e-07 2.352967e-07 NA -1.096860e-07 8.126602e-07
# 6  Omicron (BA.2) 9.863387e-01 2.049295e-03 NA  9.823222e-01 9.903553e-01
# 7           Other 4.069406e-03 1.255910e-03 NA  1.607868e-03 6.530943e-03

# estimates proportion of lab diagnoses that is B.1.617.2 by province
multinom_preds_today_byprovince = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1, by="province",
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_delta_today_byprovince_BA1 = multinom_preds_today_byprovince[multinom_preds_today_byprovince$LINEAGE=="Omicron (BA.1*)",]
multinom_delta_today_byprovince_BA1 = multinom_delta_today_byprovince_BA1[order(multinom_delta_today_byprovince_BA1$prob, decreasing=T),]
multinom_delta_today_byprovince_BA1
# LINEAGE                       province        prob           SE df   asymp.LCL   asymp.UCL
# 29 Omicron (BA.1*)                           Luik 0.025668180 0.0032988626 NA 0.019202528 0.032133832
# 15 Omicron (BA.1*)                     Henegouwen 0.025452057 0.0044249234 NA 0.016779366 0.034124747
# 64 Omicron (BA.1*)                  Waals-Brabant 0.018259187 0.0030849341 NA 0.012212827 0.024305547
# 36 Omicron (BA.1*)                      Luxemburg 0.018093994 0.0034685050 NA 0.011295849 0.024892139
# 43 Omicron (BA.1*)                          Namen 0.017939154 0.0020305780 NA 0.013959294 0.021919013
# 57 Omicron (BA.1*)                 Vlaams-Brabant 0.017117370 0.0025148356 NA 0.012188383 0.022046358
# 1  Omicron (BA.1*)                      Antwerpen 0.015015119 0.0015464575 NA 0.011984118 0.018046120
# 50 Omicron (BA.1*)                Oost-Vlaanderen 0.014431209 0.0018891189 NA 0.010728604 0.018133814
# 22 Omicron (BA.1*)                        Limburg 0.012625804 0.0015095171 NA 0.009667204 0.015584403
# 8  Omicron (BA.1*) Brussels Hoofdstedelijk Gewest 0.012204270 0.0013423685 NA 0.009573276 0.014835264
# 71 Omicron (BA.1*)                West-Vlaanderen 0.006839692 0.0007432037 NA 0.005383040 0.008296345

multinom_delta_today_byprovince_BA2 = multinom_preds_today_byprovince[multinom_preds_today_byprovince$LINEAGE=="Omicron (BA.2)",]
multinom_delta_today_byprovince_BA2 = multinom_delta_today_byprovince_BA2[order(multinom_delta_today_byprovince_BA2$prob, decreasing=T),]
multinom_delta_today_byprovince_BA2
# LINEAGE                       province      prob          SE df asymp.LCL asymp.UCL
# 76 Omicron (BA.2)                West-Vlaanderen 0.9881878 0.001732267 NA 0.9847926 0.9915830
# 13 Omicron (BA.2) Brussels Hoofdstedelijk Gewest 0.9861018 0.001544069 NA 0.9830755 0.9891281
# 55 Omicron (BA.2)                Oost-Vlaanderen 0.9838147 0.002132153 NA 0.9796358 0.9879937
# 6  Omicron (BA.2)                      Antwerpen 0.9814586 0.002010690 NA 0.9775177 0.9853994
# 62 Omicron (BA.2)                 Vlaams-Brabant 0.9813439 0.002736932 NA 0.9759796 0.9867082
# 69 Omicron (BA.2)                  Waals-Brabant 0.9803224 0.003327224 NA 0.9738012 0.9868437
# 48 Omicron (BA.2)                          Namen 0.9779121 0.002744749 NA 0.9725325 0.9832917
# 27 Omicron (BA.2)                        Limburg 0.9776690 0.003509697 NA 0.9707901 0.9845479
# 41 Omicron (BA.2)                      Luxemburg 0.9721262 0.006463440 NA 0.9594581 0.9847943
# 20 Omicron (BA.2)                     Henegouwen 0.9682855 0.005832492 NA 0.9568540 0.9797169
# 34 Omicron (BA.2)                           Luik 0.9670143 0.004583982 NA 0.9580299 0.9759987


# reorder provinces by incidence of delta
levels_PROVINCES = as.character(multinom_delta_today_byprovince$province[order(multinom_delta_today_byprovince$prob, decreasing=T)])
data_agbyweekregion2$province = factor(data_agbyweekregion2$province, levels=levels_PROVINCES)
data_agbyweekregion2$LINEAGE = factor(data_agbyweekregion2$LINEAGE, levels=levels_LINEAGE)

# redo plot of raw data by province using this order
muller_belgiumbyprovince_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE)) +
  facet_wrap(~ province, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
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
fit_belgium_multi_predsbyprovince$LINEAGE = factor(fit_belgium_multi_predsbyprovince$LINEAGE, levels=levels_LINEAGE)
fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=levels_PROVINCES)

fit_belgium_multi_preds = data.frame(emmeans(fit3_belgium_multi, 
                                           ~ LINEAGE,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                           mode="prob", df=NA))
fit_belgium_multi_preds$collection_date = as.Date(fit_belgium_multi_preds$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_preds$prob[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE!="Other"] = 0
fit_belgium_multi_preds$prob[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE=="Other"] = 1
fit_belgium_multi_preds$asymp.LCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE!="Other"] = 0
fit_belgium_multi_preds$asymp.LCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE=="Other"] = 1
fit_belgium_multi_preds$asymp.UCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE!="Other"] = 0
fit_belgium_multi_preds$asymp.UCL[fit_belgium_multi_preds$collection_date<=as.Date("2020-09-01")&fit_belgium_multi_preds$LINEAGE=="Other"] = 1
fit_belgium_multi_preds$LINEAGE = factor(fit_belgium_multi_preds$LINEAGE, levels=levels_LINEAGE) 

muller_belgium_mfit = ggplot(data=fit_belgium_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE)) + 
  # facet_wrap(~ province, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols) +
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
  scale_fill_manual("", values=lineage_cols) +
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

# OVERALL FOR THE WHOLE OF BELGIUM 

fit_belgium_multi_preds2 = fit_belgium_multi_preds
fit_belgium_multi_preds2$LINEAGE = factor(fit_belgium_multi_preds2$LINEAGE, levels=levels_LINEAGE)
levels(fit_belgium_multi_preds2$LINEAGE)

# on logit scale:

fit_belgium_multi_preds2 = fit_belgium_multi_preds
ymin = 0.001
ymax = 0.999
fit_belgium_multi_preds2$asymp.LCL[fit_belgium_multi_preds2$asymp.LCL<ymin] = ymin
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$asymp.UCL<ymin] = ymin
fit_belgium_multi_preds2$asymp.UCL[fit_belgium_multi_preds2$asymp.UCL>ymax] = ymax
fit_belgium_multi_preds2$prob[fit_belgium_multi_preds2$prob<ymin] = ymin

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
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
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
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
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
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
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
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
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




# project multinomial fit onto case data ####

source("scripts/downloadData.R") # download latest data with new confirmed cases per day from Sciensano website, code adapted from https://github.com/JoFAM/covidBE_analysis by Joris Meys
tail(cases_tot)
tail(cases_prov)



# cases_tot = cases_tot[cases_tot$DATE<(as.Date(Sys.time())-3),] # we leave out last 3 days since they are incomplete
range(cases_tot$DATE) # "2020-03-01" "2021-07-18"



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
# fit_belgium_multi_predsbyprovince = gather(fit_belgium_multi_predsbyprovince, LINEAGE, prob, all_of(levels_LINEAGE))
# fit_belgium_multi_predsbyprovince$collection_date = as.Date(fit_belgium_multi_predsbyprovince$DATE_NUM, origin="1970-01-01")
# fit_belgium_multi_predsbyprovince$LINEAGE = factor(fit_belgium_multi_predsbyprovince$LINEAGE, levels=levels_LINEAGE)
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
#   geom_smooth(aes(lwd=I(1), colour=LINEAGE, group=LINEAGE), method="loess", span=0.3, se=FALSE) +
#   geom_smooth(data=cases_belgium_byprovince2, aes(x=Date, y=newcases, lwd=I(1.5)), method="loess", span=0.3, se=FALSE, colour="black") +
#   # geom_line(aes(lwd=I(1), colour=LINEAGE, group=LINEAGE)) +
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                      limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("New confirmed cases per day") +
#   ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN BELGIUM\n(multinomial fit)") +
#   scale_colour_manual("lineage", values=lineage_cols) +
#   scale_y_log10() +
#   coord_cartesian(ylim=c(1,NA)) # +
# # coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_belgium_multi_predsbyprovince2$collection_date)-20))
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_confirmed cases multinomial fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_confirmed cases multinomial fit.pdf"), width=8, height=6)
# 
# ggplot(data=fit_belgium_multi_predsbyprovince2, 
#        aes(x=collection_date, y=cases, group=LINEAGE)) + 
#   facet_wrap(~ province, scale="free", ncol=3) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
#   scale_fill_manual("", values=lineage_cols) +
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
