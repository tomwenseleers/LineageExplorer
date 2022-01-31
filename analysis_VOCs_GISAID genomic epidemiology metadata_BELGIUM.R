# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN BELGIUM BASED ON ANALYSIS OF GISAID PATIENT METADATA
# T. Wenseleers
# last update 17 January 2022

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
d13 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2022_01_17_08_belgium_subm_jan_2022.tsv", col_types = cols(.default = "c"))
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
GISAID_belgium1 = rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13)
GISAID_belgium1 = GISAID_belgium1[GISAID_belgium1[,"Host"]=="Human",]
GISAID_belgium1[,"Collection date"] = as.Date(GISAID_belgium1[,"Collection date"])
GISAID_belgium1 = GISAID_belgium1[!is.na(GISAID_belgium1[,"Collection date"]),]
nrow(GISAID_belgium1) # 73235
unique(GISAID_belgium1[,"Virus name"])
unique(GISAID_belgium1[,"Location"])
unique(GISAID_belgium1[,"Host"])
unique(GISAID_belgium1[,"Additional location information"]) # ZIP codes # sometimes has additional info on whether it was a traveller

# anonymous records without location info are assumed to be active surveillance of cluster outbreaks & included in active surveillance (cf weekly Sciensano report)
GISAID_belgium1[,"Additional location information"][grepl("nknown",GISAID_belgium1[,"Additional location information"])] = NA
sum(is.na(GISAID_belgium1[,"Additional location information"])) # 8718
unique(GISAID_belgium1[,"Additional location information"])

library(stringr)
ZIP = str_extract(GISAID_belgium1[,"Additional location information"],"[0-9]{4}") # extract ZIP codes
unique(ZIP)
sum(!is.na(ZIP)) # 64309
sum(is.na(ZIP)) # 8926 - anonymous records, these are assumed to be active surveillance of cluster outbreaks

muni = read.csv(".//data//GISAID//Belgium//mapping_municips.csv") # post codes & cities & provinces
province = muni$province_label_nl[match(ZIP, muni$postcode)] # province
city = muni$municipality_label_nl[match(ZIP, muni$postcode)] # city
sum(is.na(city)) # 9028
sum(is.na(ZIP))  # 8926
wrongZIP = which(is.na(city)&!is.na(ZIP))
ZIP[wrongZIP] = NA

anonym = grepl("Belgium$",GISAID_belgium1[,"Additional location information"])|
  is.na(GISAID_belgium1[,"Additional location information"])|
  grepl("CHU|infected|travel|Travel",GISAID_belgium1[,"Additional location information"]) |
  is.na(ZIP)
sum(anonym) # 9142
unique(GISAID_belgium1[!is.na(GISAID_belgium1[,"Sampling strategy"]),"Sampling strategy"])
traveller = grepl("travel|Travel|infected",GISAID_belgium1[,"Additional location information"])|grepl("travel|Travel",GISAID_belgium1[,"Sampling strategy"])|
  grepl("travel|Travel|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(traveller) # 2381
Sdropout = grepl("S-gene dropout",GISAID_belgium1[,"Sampling strategy"])|
  grepl("S-gene dropout",GISAID_belgium1[,"Additional host information"])  
sum(Sdropout) # 21 - should be much higher though...
actsurv = anonym|traveller|grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster",GISAID_belgium1[,"Sampling strategy"])|
  grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster|HR contact|Outbreak|Nurse|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(actsurv) # 19127
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")]) # 16426 # In Sciensano weekly report of 21st of May 2021 this is 4710
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")&GISAID_belgium1[,"Collection date"]<=as.Date("2021-01-31")]) # 1204 # In Sciensano weekly report of 21st of May 2021 this is 1975

reinfection = grepl("Breakthrough",GISAID_belgium1[,"Sampling strategy"])
sum(reinfection) # 307
infectionpostvaccination = grepl("Vaccination|vaccination",GISAID_belgium1[,"Sampling strategy"])
sum(infectionpostvaccination) # 4638

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
sum(bassurv) # 54108
purpose_of_sequencing = rep("active_surveillance", nrow(GISAID_belgium1))
purpose_of_sequencing[bassurv] = "baseline_surveillance"

accessionIDs_bassurv = GISAID_belgium1[,"Accession ID"][bassurv]

GISAID_belgium1[,"Gender"][grepl("F|V",GISAID_belgium1[,"Gender"])] = "F"
GISAID_belgium1[,"Gender"][grepl("M|m",GISAID_belgium1[,"Gender"])] = "M"
GISAID_belgium1[,"Gender"][grepl("unknown",GISAID_belgium1[,"Gender"])] = NA
GISAID_belgium1[,"Gender"][!grepl("F|M",GISAID_belgium1[,"Gender"])] = NA
unique(GISAID_belgium1[,"Gender"]) # "F" "M" NA 

sum(!is.na(GISAID_belgium1[,"Gender"] )) # 66105 with gender info

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
sum(!is.na(GISAID_belgium1[,"Patient age"] )) # 65963 with age info

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
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2022-01-14_18-50.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID_genepi_belgium = GISAID[GISAID$country=="Belgium",]
nrow(GISAID_genepi_belgium) # 79086

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
write.csv(GISAID_belgium1, ".//data//GISAID//Belgium//gisaid_hcov-19_2022_01_17_08_ALL PARSED.csv",row.names=F)




# ANALYSIS OF VOC LINEAGE FREQUENCIES IN BELGIUM USING MULTINOMIAL FITS ####


library(dplyr)
GISAID_belgium1$LINEAGE = case_when(
  grepl("B.1.1.529|BA.1", GISAID_belgium1$pango_lineage) ~ "Omicron (BA.1)", # B.1.1.529|
  grepl("BA.2", GISAID_belgium1$pango_lineage) ~ "Omicron (BA.2)", 
  grepl("B.1.617.2", GISAID_belgium1$pango_lineage, fixed=T) | grepl("AY", GISAID_belgium1$pango_lineage)  ~ "Delta",
  grepl("B.1.1.7", GISAID_belgium1$pango_lineage, fixed=T) ~ "Alpha",
  grepl("B.1.351", GISAID_belgium1$pango_lineage, fixed=T) ~ "Beta",
  grepl("P.1", GISAID_belgium1$pango_lineage, fixed=T) ~ "Gamma",
  T ~ "Other"
)
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
                  "Omicron (BA.1)", "Omicron (BA.2)") 
levels_LINEAGE = c(main_lineages, "Other")
# lineage_cols = c("#0085FF","#9A9D00","cyan3","magenta",
#                   muted("red",c=150,l=65),muted("red",c=150,l=25),
#                   "grey70")
lineage_cols = c("#0085FF","#9A9D00","cyan3","magenta",
                  colorRampPalette(c("red", "orange", "blue"))(2),
                  "grey70")
unique(GISAID_belgium1$LINEAGE)
GISAID_belgium1$LINEAGE[!(GISAID_belgium1$LINEAGE %in% main_lineages)] = "Other" # minority lineages & non-VOCs
table(GISAID_belgium1$LINEAGE)
GISAID_belgium1$LINEAGE = factor(GISAID_belgium1$LINEAGE, 
                                  levels=levels_LINEAGE, 
                                  labels=levels_LINEAGE)
table(GISAID_belgium1$LINEAGE)
# Alpha           Beta          Gamma          Delta Omicron (BA.1) Omicron (BA.2)          Other 
# 18461           1047           1725          41043           2516             18           8425 
sum(table(GISAID_belgium1$LINEAGE)) # 73235

# GISAID_belgium = GISAID_sel[GISAID_sel$country=="Belgium",]
GISAID_belgium = GISAID_belgium1 # [GISAID_belgium1$purpose_of_sequencing=="baseline_surveillance",] # if desired subset to baseline surveillance
# GISAID_belgium = GISAID_belgium[-which((GISAID_belgium$LINEAGE=="B.1.617.2")&(GISAID_belgium$date<="2021-05-01")),]
nrow(GISAID_belgium) # 73235
sum(!is.na(GISAID_belgium$age)) # 65963 with age info

# use data from Feb 1 onwards
GISAID_belgium = GISAID_belgium[GISAID_belgium$date>=as.Date("2021-02-01"),]
nrow(GISAID_belgium) #
range(GISAID_belgium$date) # "2021-02-01" "2022-01-10"


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

write.csv(data_agbyweek2, ".//data//GISAID//Belgium//gisaid_hcov-19_2022_01_17_08_BASELINE SELECTION_aggregated counts by week.csv",row.names=F)


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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
data_agbyweekregion2$LINEAGE = relevel(data_agbyweekregion2$LINEAGE, ref="Omicron (BA.1)") 
fit1_belgium_multi = nnet::multinom(LINEAGE ~ province + DATE_NUM, weights=count, data=data_agbyweekregion2, maxit=1000)
fit2_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
fit3_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
fit4_belgium_multi = nnet::multinom(LINEAGE ~ province * DATE_NUM, weights=count, data=data_agbyweekregion2, maxit=1000)
fit5_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
fit6_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
BIC(fit1_belgium_multi, fit2_belgium_multi, fit3_belgium_multi, fit4_belgium_multi, fit5_belgium_multi, fit6_belgium_multi) 
# fit3_belgium_multi fits best (lowest BIC)

# growth rate advantage compared to Omicron (BA.1) (difference in growth rate per day) 
emtrbelgium = emtrends(fit3_belgium_multi, trt.vs.ctrl ~ LINEAGE,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium = data.frame(confint(emtrbelgium, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrbelgium$contrasts)$p.value)
delta_r_belgium
# contrast    estimate          SE df  asymp.LCL   asymp.UCL      p.value
# 1          Alpha - Omicron (BA.1) -0.10413057 0.015987065 NA -0.1354646 -0.07279650 3.071108e-08
# 2           Beta - Omicron (BA.1) -0.04633638 0.031304273 NA -0.1076916  0.01501886 4.787073e-01
# 3          Gamma - Omicron (BA.1) -0.07100060 0.014913581 NA -0.1002307 -0.04177052 4.673001e-05
# 4          Delta - Omicron (BA.1) -0.14457236 0.007819939 NA -0.1598992 -0.12924556 1.318344e-10
# 5 Omicron (BA.2) - Omicron (BA.1)  0.20010030 0.049666145 NA  0.1027564  0.29744415 7.031223e-04
# 6          Other - Omicron (BA.1) -0.09242247 0.008910148 NA -0.1098860 -0.07495890 1.318609e-10

# pairwise growth rate difference (differences in growth rate per day) 
emtrbelgium_pairw = emtrends(fit3_belgium_multi, pairwise ~ LINEAGE,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#                         contrast    estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1           Omicron (BA.1) - Alpha  0.10413057 0.015987065 NA  0.072796503  0.135464646 1.064215e-07
# 2            Omicron (BA.1) - Beta  0.04633638 0.031304273 NA -0.015018865  0.107691630 7.556511e-01
# 3           Omicron (BA.1) - Gamma  0.07100060 0.014913581 NA  0.041770516  0.100230679 1.574078e-04
# 4           Omicron (BA.1) - Delta  0.14457236 0.007819939 NA  0.129245558  0.159899155 1.582013e-10
# 5  Omicron (BA.1) - Omicron (BA.2) -0.20010030 0.049666145 NA -0.297444151 -0.102756441 2.273156e-03
# 6           Omicron (BA.1) - Other  0.09242247 0.008910148 NA  0.074958905  0.109886043 1.582575e-10
# 7                     Alpha - Beta -0.05779419 0.032781624 NA -0.122044993  0.006456610 5.763410e-01
# 8                    Alpha - Gamma -0.03312998 0.017894463 NA -0.068202480  0.001942527 5.178721e-01
# 9                    Alpha - Delta  0.04044178 0.013982038 NA  0.013037492  0.067846073 6.967950e-02
# 10          Alpha - Omicron (BA.2) -0.30423087 0.052111149 NA -0.406366846 -0.202094895 1.984969e-06
# 11                   Alpha - Other -0.01170810 0.013844020 NA -0.038841881  0.015425680 9.792232e-01
# 12                    Beta - Gamma  0.02466421 0.032651019 NA -0.039330607  0.088659036 9.884092e-01
# 13                    Beta - Delta  0.09823597 0.030359358 NA  0.038732725  0.157739223 2.773469e-02
# 14           Beta - Omicron (BA.2) -0.24643668 0.058646069 NA -0.361380862 -0.131492495 1.241968e-03
# 15                    Beta - Other  0.04608609 0.030459977 NA -0.013614368  0.105786550 7.363156e-01
# 16                   Gamma - Delta  0.07357176 0.012759046 NA  0.048564489  0.098579029 2.691296e-06
# 17          Gamma - Omicron (BA.2) -0.27110089 0.051791746 NA -0.372610850 -0.169590937 2.429803e-05
# 18                   Gamma - Other  0.02142188 0.013040634 NA -0.004137296  0.046981049 6.553355e-01
# 19          Delta - Omicron (BA.2) -0.34467265 0.050211796 NA -0.443085965 -0.246259340 2.247302e-08
# 20                   Delta - Other -0.05214988 0.004431786 NA -0.060836024 -0.043463742 1.582182e-10
# 21          Omicron (BA.2) - Other  0.29252277 0.050391942 NA  0.193756378  0.391289162 2.284758e-06


# estimated proportion of different LINEAGES among new lab diagnoses in Belgium today
today # "2021-10-25"
multinom_preds_today_avg = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#       LINEAGE         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1) 9.122829e-01 3.616944e-02 NA  8.413921e-01 9.831737e-01
# 2          Alpha 1.315971e-05 1.990305e-05 NA -2.584956e-05 5.216898e-05
# 3           Beta 1.466360e-05 4.816685e-05 NA -7.974169e-05 1.090689e-04
# 4          Gamma 4.928127e-05 6.540275e-05 NA -7.890576e-05 1.774683e-04
# 5          Delta 2.489191e-02 4.188292e-03 NA  1.668300e-02 3.310081e-02
# 6 Omicron (BA.2) 6.211802e-02 3.692763e-02 NA -1.025880e-02 1.344948e-01
# 7          Other 6.301076e-04 2.213710e-04 NA  1.962285e-04 1.063987e-03


# estimated proportion of different LINEAGES among new infections in Belgium today
today+7 # "2021-07-29"
multinom_preds_infections_today_avg = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num+7), 
                                              mode="prob", df=NA))
multinom_preds_infections_today_avg
#       LINEAGE         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1) 8.388260e-01 7.377036e-02 NA  6.942387e-01 9.834132e-01
# 2          Alpha 5.498468e-06 8.911624e-06 NA -1.196799e-05 2.296493e-05
# 3           Beta 8.449307e-06 2.953691e-05 NA -4.944198e-05 6.634059e-05
# 4          Gamma 2.695007e-05 3.836011e-05 NA -4.823435e-05 1.021345e-04
# 5          Delta 8.043606e-03 1.995497e-03 NA  4.132503e-03 1.195471e-02
# 6 Omicron (BA.2) 1.528054e-01 7.468090e-02 NA  6.433506e-03 2.991772e-01
# 7          Other 2.841301e-04 1.210616e-04 NA  4.685376e-05 5.214064e-04

# estimates proportion of lab diagnoses that is B.1.617.2 by province
multinom_preds_today_byprovince = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1, by="province",
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_delta_today_byprovince_BA1 = multinom_preds_today_byprovince[multinom_preds_today_byprovince$LINEAGE=="Omicron (BA.1)",]
multinom_delta_today_byprovince_BA1 = multinom_delta_today_byprovince_BA1[order(multinom_delta_today_byprovince_BA1$prob, decreasing=T),]
multinom_delta_today_byprovince_BA1
# LINEAGE                       province      prob          SE df  asymp.LCL asymp.UCL
# 8  Omicron (BA.1) Brussels Hoofdstedelijk Gewest 0.9911264 0.001674592 NA 0.98784424 0.9944085
# 50 Omicron (BA.1)                Oost-Vlaanderen 0.9855645 0.003510507 NA 0.97868405 0.9924450
# 15 Omicron (BA.1)                     Henegouwen 0.9831935 0.004306700 NA 0.97475254 0.9916345
# 64 Omicron (BA.1)                  Waals-Brabant 0.9783258 0.005852965 NA 0.96685424 0.9897974
# 57 Omicron (BA.1)                 Vlaams-Brabant 0.9773491 0.005005659 NA 0.96753822 0.9871600
# 29 Omicron (BA.1)                           Luik 0.9719792 0.005284003 NA 0.96162278 0.9823357
# 36 Omicron (BA.1)                      Luxemburg 0.9552698 0.012129075 NA 0.93149723 0.9790423
# 43 Omicron (BA.1)                          Namen 0.9100533 0.065945089 NA 0.78080329 1.0393033
# 22 Omicron (BA.1)                        Limburg 0.8840068 0.053475317 NA 0.77919713 0.9888165
# 71 Omicron (BA.1)                West-Vlaanderen 0.7325691 0.078073205 NA 0.57954841 0.8855897
# 1  Omicron (BA.1)                      Antwerpen 0.6656739 0.294897396 NA 0.08768564 1.2436622

multinom_delta_today_byprovince_BA2 = multinom_preds_today_byprovince[multinom_preds_today_byprovince$LINEAGE=="Omicron (BA.2)",]
multinom_delta_today_byprovince_BA2 = multinom_delta_today_byprovince_BA2[order(multinom_delta_today_byprovince_BA2$prob, decreasing=T),]
multinom_delta_today_byprovince_BA2
#           LINEAGE                       province         prob           SE df     asymp.LCL    asymp.UCL
# 6  Omicron (BA.2)                      Antwerpen 3.200900e-01 3.011557e-01 NA -2.701644e-01 9.103443e-01
# 76 Omicron (BA.2)                West-Vlaanderen 2.256347e-01 8.218885e-02 NA  6.454753e-02 3.867219e-01
# 27 Omicron (BA.2)                        Limburg 7.354377e-02 5.540394e-02 NA -3.504595e-02 1.821335e-01
# 48 Omicron (BA.2)                          Namen 6.400461e-02 6.758194e-02 NA -6.845356e-02 1.964628e-01
# 62 Omicron (BA.2)                 Vlaams-Brabant 2.319392e-05 2.065104e-03 NA -4.024336e-03 4.070724e-03
# 20 Omicron (BA.2)                     Henegouwen 1.906813e-06 2.638731e-06 NA -3.265005e-06 7.078631e-06
# 34 Omicron (BA.2)                           Luik 3.081825e-13 4.264642e-13 NA -5.276720e-13 1.144037e-12
# 13 Omicron (BA.2) Brussels Hoofdstedelijk Gewest 1.300263e-13 1.799336e-13 NA -2.226371e-13 4.826897e-13
# 69 Omicron (BA.2)                  Waals-Brabant 2.984086e-14 4.129416e-14 NA -5.109420e-14 1.107759e-13
# 55 Omicron (BA.2)                Oost-Vlaanderen 5.821398e-19 8.055759e-19 NA -9.967599e-19 2.161039e-18
# 41 Omicron (BA.2)                      Luxemburg 4.227743e-25 5.850451e-25 NA -7.238929e-25 1.569442e-24


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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
date.from = as.numeric(as.Date("2021-02-01"))
date.to = today_num+extrapolate # max(GISAID_belgium$DATE_NUM)+extrapolate

fit_belgium_multi_predsbyprovince = data.frame(emmeans(fit3_belgium_multi,
                                                  ~ LINEAGE,
                                                  by=c("DATE_NUM", "province"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                                  mode="prob", df=NA))
fit_belgium_multi_predsbyprovince$collection_date = as.Date(fit_belgium_multi_predsbyprovince$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_predsbyprovince$LINEAGE = factor(fit_belgium_multi_predsbyprovince$LINEAGE, levels=levels_LINEAGE)
fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=levels_PROVINCES)

fit_belgium_multi_preds = data.frame(emmeans(fit3_belgium_multi, 
                                           ~ LINEAGE,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                           mode="prob", df=NA))
fit_belgium_multi_preds$collection_date = as.Date(fit_belgium_multi_preds$DATE_NUM, origin="1970-01-01")
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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +  # , axis.title.x=element_blank()
  ylab("Share") + xlab("Collection date") + 
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)")
muller_belgium_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_multinom fit.png"), width=8, height=5)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots_multinom fit.pdf"), width=8, height=5)


library(ggpubr)
ggarrange(muller_belgium_raw2+coord_cartesian(xlim=c(as.Date("2021-02-01"),as.Date(date.to,origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))) + theme(axis.title.x=element_blank()), 
          muller_belgium_mfit+ggtitle("Multinomial fit")+coord_cartesian(xlim=c(as.Date("2021-02-01"),as.Date(date.to,origin="1970-01-01"))), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots multipanel_multinom fit.png"), width=8, height=8)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots multipanel_multinom fit.pdf"), width=8, height=8)


muller_belgiumbyprovince_mfit = ggplot(data=fit_belgium_multi_predsbyprovince,
                                  aes(x=collection_date, y=prob, group=LINEAGE)) +
  facet_wrap(~ province, ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  annotate("rect", xmin=max(GISAID_belgium$DATE_NUM)+1, 
           xmax=date.to, ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("Share") + xlab("Collection date") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\n(GISAID data, multinomial fit)")
muller_belgiumbyprovince_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_multinom fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_muller plots by province_multinom fit.pdf"), width=8, height=6)

ggarrange(muller_belgiumbyprovince_raw2 +
            coord_cartesian(xlim=c(as.Date("2021-02-01"),as.Date(date.to,origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0)))+
            ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM\nRaw GISAID data") + theme(axis.title.x=element_blank()),
          muller_belgiumbyprovince_mfit+ggtitle("\nMultinomial fit")+
            coord_cartesian(xlim=c(as.Date("2021-02-01"),as.Date(date.to,origin="1970-01-01"))), nrow=2)

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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
  coord_cartesian(xlim=c(as.Date("2021-02-01"),NA), ylim=c(0.001, 0.999), expand=c(0,0))
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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-02-01",NA)),
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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
  coord_cartesian(xlim=c(as.Date("2021-02-01"),NA), ylim=c(0.001, 0.999), expand=c(0,0))
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
                     limits=as.Date(c("2021-02-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-02-01",NA)),
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
