# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN BELGIUM BASED ON ANALYSIS OF GISAID PATIENT METADATA
# T. Wenseleers
# last update 25 October 2021

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
today # 2021-10-25
today_num = as.numeric(today)
plotdir = "BE_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import & parse manually downloaded GISAID patient metadata for Belgium (downloaded 30/6/2021) ####
d1 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_08_belgium_subm_jan_dec_2020.tsv", col_types = cols(.default = "c"))
d2 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_08_belgium_subm_jan_feb_2021.tsv", col_types = cols(.default = "c"))
d3 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_05_29_09_belgium_subm_mar_apr_2021.tsv", col_types = cols(.default = "c"))
d4 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_06_18_22_belgium_subm_may_2021.tsv", col_types = cols(.default = "c")) 
d5 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_07_22_19_belgium_subm_june_2021.tsv", col_types = cols(.default = "c"))
d6 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_08_31_17_belgium_subm_july_2021.tsv", col_types = cols(.default = "c"))
d7 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_10_25_18_belgium_subm_aug_2021.tsv", col_types = cols(.default = "c"))
d8 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_10_25_18_belgium_subm_sept_2021.tsv", col_types = cols(.default = "c"))
d9 = read_tsv(".//data//GISAID//Belgium//gisaid_hcov-19_2021_10_25_18_belgium_subm_oct_2021.tsv", col_types = cols(.default = "c"))
d1 = as.data.frame(d1)
d2 = as.data.frame(d2)
d3 = as.data.frame(d3)
d4 = as.data.frame(d4)
d5 = as.data.frame(d5)
d6 = as.data.frame(d6)
d7 = as.data.frame(d7)
d8 = as.data.frame(d8)
d9 = as.data.frame(d9)
GISAID_belgium1 = rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9)
GISAID_belgium1 = GISAID_belgium1[GISAID_belgium1[,"Host"]=="Human",]
GISAID_belgium1[,"Collection date"] = as.Date(GISAID_belgium1[,"Collection date"])
GISAID_belgium1 = GISAID_belgium1[!is.na(GISAID_belgium1[,"Collection date"]),]
nrow(GISAID_belgium1) # 58321
unique(GISAID_belgium1[,"Virus name"])
unique(GISAID_belgium1[,"Location"])
unique(GISAID_belgium1[,"Host"])
unique(GISAID_belgium1[,"Additional location information"]) # ZIP codes # sometimes has additional info on whether it was a traveller

# anonymous records without location info are assumed to be active surveillance of cluster outbreaks & included in active surveillance (cf weekly Sciensano report)
GISAID_belgium1[,"Additional location information"][grepl("nknown",GISAID_belgium1[,"Additional location information"])] = NA
sum(is.na(GISAID_belgium1[,"Additional location information"])) # 7086
unique(GISAID_belgium1[,"Additional location information"])

library(stringr)
ZIP = str_extract(GISAID_belgium1[,"Additional location information"],"[0-9]{4}") # extract ZIP codes
unique(ZIP)
sum(!is.na(ZIP)) # 50118
sum(is.na(ZIP)) # 8203 - anonymous records, these are assumed to be active surveillance of cluster outbreaks

muni = read.csv(".//data//GISAID//Belgium//mapping_municips.csv") # post codes & cities & provinces
province = muni$province_label_nl[match(ZIP, muni$postcode)] # province
city = muni$municipality_label_nl[match(ZIP, muni$postcode)] # city
sum(is.na(city)) # 8281
sum(is.na(ZIP))  # 8203
wrongZIP = which(is.na(city)&!is.na(ZIP))
ZIP[wrongZIP] = NA

anonym = grepl("Belgium$",GISAID_belgium1[,"Additional location information"])|
  is.na(GISAID_belgium1[,"Additional location information"])|
  grepl("CHU|infected|travel|Travel",GISAID_belgium1[,"Additional location information"]) |
  is.na(ZIP)
sum(anonym) # 8331
unique(GISAID_belgium1[!is.na(GISAID_belgium1[,"Sampling strategy"]),"Sampling strategy"])
traveller = grepl("travel|Travel|infected",GISAID_belgium1[,"Additional location information"])|grepl("travel|Travel",GISAID_belgium1[,"Sampling strategy"])|
  grepl("travel|Travel|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(traveller) # 1403
Sdropout = grepl("S-gene dropout",GISAID_belgium1[,"Sampling strategy"])|
  grepl("S-gene dropout",GISAID_belgium1[,"Additional host information"])  
sum(Sdropout) # 61 - should be much higher though...
actsurv = anonym|traveller|grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster",GISAID_belgium1[,"Sampling strategy"])|
  grepl("Longitudinal|Outbreak|Active|active|dropout|Atypical|travel|Travel|Atypische|CLUSTER|cluster|Cluster|HR contact|Outbreak|Nurse|Holiday",GISAID_belgium1[,"Additional host information"])  
sum(actsurv) # 15677
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")]) # 12935 # In Sciensano weekly report of 21st of May 2021 this is 4710
sum(actsurv[GISAID_belgium1[,"Collection date"]>=as.Date("2020-11-30")&GISAID_belgium1[,"Collection date"]<=as.Date("2021-01-31")]) # 1207 # In Sciensano weekly report of 21st of May 2021 this is 1975

reinfection = grepl("Breakthrough",GISAID_belgium1[,"Sampling strategy"])
sum(reinfection) # 135
infectionpostvaccination = grepl("Vaccination|vaccination",GISAID_belgium1[,"Sampling strategy"])
sum(infectionpostvaccination) # 3357

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
sum(bassurv) # 42644
purpose_of_sequencing = rep("active_surveillance", nrow(GISAID_belgium1))
purpose_of_sequencing[bassurv] = "baseline_surveillance"

accessionIDs_bassurv = GISAID_belgium1[,"Accession ID"][bassurv]

GISAID_belgium1[,"Gender"][grepl("F|V",GISAID_belgium1[,"Gender"])] = "F"
GISAID_belgium1[,"Gender"][grepl("M|m",GISAID_belgium1[,"Gender"])] = "M"
GISAID_belgium1[,"Gender"][grepl("unknown",GISAID_belgium1[,"Gender"])] = NA
GISAID_belgium1[,"Gender"][!grepl("F|M",GISAID_belgium1[,"Gender"])] = NA
unique(GISAID_belgium1[,"Gender"]) # "F" "M" NA 

sum(!is.na(GISAID_belgium1[,"Gender"] )) # 52791 with gender info

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
sum(!is.na(GISAID_belgium1[,"Patient age"] )) # 52760 with age info

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
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-10-22_10-28.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID_genepi_belgium = GISAID[GISAID$country=="Belgium",]
nrow(GISAID_genepi_belgium) # 38994

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
write.csv(GISAID_belgium1, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_10_25_ALL PARSED.csv",row.names=F)




# ANALYSIS OF VOC LINEAGE FREQUENCIES IN BELGIUM USING MULTINOMIAL FITS ####

# selected lineages
main_lineages = c("B.1.1.7","B.1.351","P.1","B.1.617.2",
                  "AY.4","AY.5","AY.5.2","AY.9","AY.10","AY.33","AY.34","AY.39","AY.4.2") 
levels_LINEAGE2 = c(main_lineages, "other")
lineage_cols2 = c("#0085FF","#9A9D00","cyan3","magenta",
                  muted("red",c=150,l=65),muted("red",c=150,l=60),muted("red",c=150,l=55),muted("red",c=150,l=50),muted("red",c=150,l=45),
                  muted("red",c=150,l=40),muted("red",c=150,l=35),muted("red",c=150,l=30),muted("red",c=150,l=25),
                  "grey70")
lineage_cols2 = c("#0085FF","#9A9D00","cyan3","magenta",
                  colorRampPalette(c("red", "orange", "blue"))(9),
                  "grey70")
unique(GISAID_belgium1$LINEAGE2)
GISAID_belgium1$LINEAGE2[!(GISAID_belgium1$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
table(GISAID_belgium1$LINEAGE2)
GISAID_belgium1$LINEAGE2 = factor(GISAID_belgium1$LINEAGE2, 
                                  levels=levels_LINEAGE2, 
                                  labels=levels_LINEAGE2)
table(GISAID_belgium1$LINEAGE2)
# B.1.1.7   B.1.351       P.1 B.1.617.2      AY.4      AY.5    AY.5.2      AY.9     AY.10     AY.33     AY.34     AY.39    AY.4.2     other 
# 21327      1076      1851     14121      3552      1207        89      1254      1052      1293       395       946        74     10084 
sum(table(GISAID_belgium1$LINEAGE2)) # 58321

# GISAID_belgium = GISAID_sel[GISAID_sel$country=="Belgium",]
GISAID_belgium = GISAID_belgium1 # [GISAID_belgium1$purpose_of_sequencing=="baseline_surveillance",] # if desired subset to baseline surveillance
# GISAID_belgium = GISAID_belgium[-which((GISAID_belgium$LINEAGE2=="B.1.617.2")&(GISAID_belgium$date<="2021-05-01")),]
nrow(GISAID_belgium) # 58321
sum(!is.na(GISAID_belgium$age)) # 52760 with age info

# use data from Feb 1 onwards
GISAID_belgium = GISAID_belgium[GISAID_belgium$date>=as.Date("2021-02-01"),]
nrow(GISAID_belgium) # 52031
range(GISAID_belgium$date) # "2021-02-01" "2021-10-22"


# # plot age distribution of B.1.617.2 & B.1.1.7 cases in Belgium though time
# df = GISAID_belgium[(GISAID_belgium$date>=as.Date("2021-04-01")&GISAID_belgium$pango_lineage %in% c("B.1.1.7","B.1.617.2")),]
# df$variant = factor(df$pango_lineage, levels=c("B.1.1.7","B.1.617.2"), labels=c("B.1.1.7 (alfa)","B.1.617.2 (delta)"))
# qplot(data=df[df$Week>18&df$Week<26,], binwidth=10.0, x=age, y = ..density..*100, geom="histogram", fill=variant) + 
#   facet_wrap(~floor_date+variant, ncol=2) + # , ncol=1
#   ylab("%") + scale_fill_manual(values=lineage_cols2[c(1,4)]) + theme(legend.position="none")
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_B1672_B117_age distribution.png"), width=8, height=10)


# unique(GISAID_belgium$originating_lab) # 133 labs # maybe some labs did more active surveillance & could be excluded?
# unique(GISAID_belgium$submitting_lab) # 41 labs
# sum(GISAID_belgium$country_exposure!="Belgium") # only 23 are indicated as travellers

# # we remove travellers
# GISAID_belgium = GISAID_belgium[GISAID_belgium$country_exposure=="Belgium",]
# nrow(GISAID_belgium) # 19283
# sum(GISAID_belgium$division_exposure==GISAID_belgium$division) # 19283

# unique(GISAID_belgium$province) # 
# unique(GISAID_belgium$province[GISAID_belgium$LINEAGE2=="B.1.617.2"]) # "Belgium"         "Hasselt"         "Brussels"        "Brugge"          "Gent"            "Liège"           "Halle-Vilvoorde" "Antwerpen"       "Mechelen"
# sum(GISAID_belgium$LINEAGE2=="B.1.617.2") # 171 among baseline surveillance
# unique(GISAID_belgium$province[GISAID_belgium$LINEAGE1=="B.1.1.7"])
# sum(GISAID_belgium$LINEAGE1=="B.1.1.7") # 13339 among baseline surveillance
# 
# table(GISAID_belgium$LINEAGE1)
# table(GISAID_belgium$LINEAGE2)

# main_lineages = names(table(GISAID_belgium$LINEAGE1))[100*table(GISAID_belgium$LINEAGE1)/sum(table(GISAID_belgium$LINEAGE1)) > 3]
# main_lineages
# # "B.1.1.7" "B.1.351" "P.1" 
# VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.620","B.1.1.7","B.1.351","P.1","B.1.1.207","B.1.429", # "B.1.1.318",
#          "B.1.214.2") # I added B.1.214.2 here # cut "B.1.525","B.1.526",
# main_lineages = union(main_lineages, VOCs)
# # main_lineages = c(main_lineages,"B.1","B.1.1","B.1.258")


# remove = names(table(GISAID_belgium$LINEAGE1))[table(GISAID_belgium$LINEAGE1) < 10]
# GISAID_belgium$LINEAGE1[(GISAID_belgium$LINEAGE1 %in% remove)] = "other" # minority VOCs
# GISAID_belgium$LINEAGE2[(GISAID_belgium$LINEAGE2 %in% remove)] = "other" # minority VOCs
# table(GISAID_belgium$LINEAGE1)
# GISAID_belgium$LINEAGE1 = factor(GISAID_belgium$LINEAGE1)
# GISAID_belgium$LINEAGE1 = relevel(GISAID_belgium$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
# levels(GISAID_belgium$LINEAGE1)
# "B.1.1.7"   "B.1.214.2" "B.1.351"   "B.1.617+"  "B.1.620"   "other"     "P.1"  
# levels_LINEAGE1 = c("B.1.1.7","B.1.214.2",
#                    "B.1.351","P.1","B.1.620","B.1.617+","other")
# GISAID_belgium$LINEAGE1 = factor(GISAID_belgium$LINEAGE1, levels=levels_LINEAGE1)
# 
# GISAID_belgium$LINEAGE2 = factor(GISAID_belgium$LINEAGE2)
# GISAID_belgium$LINEAGE2 = relevel(GISAID_belgium$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
# levels(GISAID_belgium$LINEAGE2)
# # "B.1.1.7"   "B.1.214.2" "B.1.351"   "B.1.617.1" "B.1.617.2" "B.1.620"   "other"     "P.1"  
# levels_LINEAGE2 = c("B.1.1.7","B.1.214.2",
#                     "B.1.351","P.1","B.1.620","B.1.617.1","B.1.617.2","other")
# GISAID_belgium$LINEAGE2 = factor(GISAID_belgium$LINEAGE2, levels=levels_LINEAGE2)

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
data_agbyweek2 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617.2","total"]) == nrow(GISAID_belgium) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL

write.csv(data_agbyweek2, ".//data//GISAID//Belgium//gisaid_hcov-19_2021_10_25_BASELINE SELECTION_aggregated counts by week.csv",row.names=F)


# aggregated by week and province for selected variant lineages
data_agbyweekregion2 = as.data.frame(table(GISAID_belgium$floor_date, GISAID_belgium$province, GISAID_belgium$LINEAGE2))
colnames(data_agbyweekregion2) = c("floor_date", "province", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + province, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$province),
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$province))]
# sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617.2","total"]) == nrow(GISAID_belgium) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
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
# n2 = length(levels(data_agbyweek2$LINEAGE2))
# lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
# lineage_cols2[which(levels(data_agbyweek2$LINEAGE2)=="B.1.617.1")] = muted("magenta")
# lineage_cols2[which(levels(data_agbyweek2$LINEAGE2)=="B.1.617.2")] = "magenta"
# lineage_cols2[which(levels(data_agbyweek2$LINEAGE2)=="other")] = "grey75"

muller_belgium_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ province, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
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


muller_belgiumbyprovince_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) +
  facet_wrap(~ province, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
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
data_agbyweekregion2$LINEAGE = data_agbyweekregion2$LINEAGE2
data_agbyweekregion2$LINEAGE = relevel(data_agbyweekregion2$LINEAGE, ref="B.1.617.2") 
fit1_belgium_multi = nnet::multinom(LINEAGE ~ province + DATE_NUM, weights=count, data=data_agbyweekregion2, maxit=1000)
fit2_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
fit3_belgium_multi = nnet::multinom(LINEAGE ~ province + ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
fit4_belgium_multi = nnet::multinom(LINEAGE ~ province * DATE_NUM, weights=count, data=data_agbyweekregion2, maxit=1000)
fit5_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion2, maxit=1000)
fit6_belgium_multi = nnet::multinom(LINEAGE ~ province * ns(DATE_NUM, df=3), weights=count, data=data_agbyweekregion2, maxit=1000)
BIC(fit1_belgium_multi, fit2_belgium_multi, fit3_belgium_multi, fit4_belgium_multi, fit5_belgium_multi, fit6_belgium_multi) 
# fit3_belgium_multi fits best (lowest BIC)

# growth rate advantage compared to B.1.617.2 (difference in growth rate per day) 
emtrbelgium = emtrends(fit3_belgium_multi, trt.vs.ctrl ~ LINEAGE,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium = data.frame(confint(emtrbelgium, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrbelgium$contrasts)$p.value)
delta_r_belgium
#               contrast      estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1  B.1.1.7 - B.1.617.2 -0.2108503124 0.012860633 NA -0.236056690 -0.185643935 1.598721e-14
# 2  B.1.351 - B.1.617.2 -0.0174127064 0.024911848 NA -0.066239031  0.031413618 9.766758e-01
# 3      P.1 - B.1.617.2 -0.1583444153 0.022075942 NA -0.201612467 -0.115076364 2.326087e-10
# 4     AY.4 - B.1.617.2  0.0202646479 0.001624817 NA  0.017080066  0.023449230 1.720846e-14
# 5     AY.5 - B.1.617.2 -0.0237397508 0.002976283 NA -0.029573157 -0.017906344 2.177591e-12
# 6   AY.5.2 - B.1.617.2  0.1276004741 0.010654295 NA  0.106718440  0.148482508 2.009504e-14
# 7     AY.9 - B.1.617.2 -0.0035846717 0.003251084 NA -0.009956678  0.002787335 8.654360e-01
# 8    AY.10 - B.1.617.2 -0.0002414089 0.003343340 NA -0.006794235  0.006311417 9.999992e-01
# 9    AY.33 - B.1.617.2 -0.0433572132 0.003618741 NA -0.050449816 -0.036264611 1.998401e-14
# 10   AY.34 - B.1.617.2 -0.0182938768 0.005227735 NA -0.028540049 -0.008047704 6.754518e-03
# 11   AY.39 - B.1.617.2 -0.0014974694 0.003234962 NA -0.007837878  0.004842939 9.959993e-01
# 12  AY.4.2 - B.1.617.2 -0.0075062258 0.023000520 NA -0.052586417  0.037573965 9.991586e-01
# 13   other - B.1.617.2 -0.0129650721 0.002188249 NA -0.017253961 -0.008676183 1.981367e-07

# pairwise growth rate difference (differences in growth rate per day) 
emtrbelgium_pairw = emtrends(fit3_belgium_multi, pairwise ~ LINEAGE,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_belgium$DATE_NUM)))
delta_r_belgium_pairw = data.frame(confint(emtrbelgium_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrbelgium_pairw$contrasts)$p.value)
delta_r_belgium_pairw
#               contrast      estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1  B.1.617.2 - B.1.1.7  0.2108503124 0.012860633 NA  0.185643935  0.236056690 1.731948e-14
# 2  B.1.617.2 - B.1.351  0.0174127064 0.024911848 NA -0.031413618  0.066239031 9.999866e-01
# 3      B.1.617.2 - P.1  0.1583444153 0.022075942 NA  0.115076364  0.201612467 1.624548e-09
# 4     B.1.617.2 - AY.4 -0.0202646479 0.001624817 NA -0.023449230 -0.017080066 2.053913e-14
# 5     B.1.617.2 - AY.5  0.0237397508 0.002976283 NA  0.017906344  0.029573157 1.491340e-11
# 6   B.1.617.2 - AY.5.2 -0.1276004741 0.010654295 NA -0.148482508 -0.106718440 2.831069e-14
# 7     B.1.617.2 - AY.9  0.0035846717 0.003251084 NA -0.002787335  0.009956678 9.979587e-01
# 8    B.1.617.2 - AY.10  0.0002414089 0.003343340 NA -0.006311417  0.006794235 1.000000e+00
# 9    B.1.617.2 - AY.33  0.0433572132 0.003618741 NA  0.036264611  0.050449816 2.819966e-14
# 10   B.1.617.2 - AY.34  0.0182938768 0.005227735 NA  0.008047704  0.028540049 3.694665e-02
# 11   B.1.617.2 - AY.39  0.0014974694 0.003234962 NA -0.004842939  0.007837878 9.999999e-01
# 12  B.1.617.2 - AY.4.2  0.0075062258 0.023000520 NA -0.037573965  0.052586417 1.000000e+00
# 13   B.1.617.2 - other  0.0129650721 0.002188249 NA  0.008676183  0.017253961 1.370432e-06
# 14   B.1.1.7 - B.1.351 -0.1934376060 0.026457720 NA -0.245293783 -0.141581429 7.357426e-10
# 15       B.1.1.7 - P.1 -0.0525058971 0.021890104 NA -0.095409712 -0.009602082 4.853806e-01
# 16      B.1.1.7 - AY.4 -0.2311149603 0.012986591 NA -0.256568210 -0.205661710 1.731948e-14
# 17      B.1.1.7 - AY.5 -0.1871105616 0.013178851 NA -0.212940636 -0.161280488 1.731948e-14
# 18    B.1.1.7 - AY.5.2 -0.3384507865 0.016698632 NA -0.371179503 -0.305722070 1.731948e-14
# 19      B.1.1.7 - AY.9 -0.2072656407 0.013241308 NA -0.233218128 -0.181313153 1.731948e-14
# 20     B.1.1.7 - AY.10 -0.2106089035 0.013288235 NA -0.236653365 -0.184564442 1.731948e-14
# 21     B.1.1.7 - AY.33 -0.1674930992 0.013332422 NA -0.193624167 -0.141362032 1.909584e-14
# 22     B.1.1.7 - AY.34 -0.1925564356 0.013872153 NA -0.219745355 -0.165367516 1.743050e-14
# 23     B.1.1.7 - AY.39 -0.2093528430 0.013271487 NA -0.235364480 -0.183341206 1.731948e-14
# 24    B.1.1.7 - AY.4.2 -0.2033440866 0.026345902 NA -0.254981105 -0.151707068 6.871836e-11
# 25     B.1.1.7 - other -0.1978852403 0.012298604 NA -0.221990061 -0.173780419 1.731948e-14
# 26       B.1.351 - P.1  0.1409317089 0.032269085 NA  0.077685465  0.204177953 1.681409e-03
# 27      B.1.351 - AY.4 -0.0376773542 0.024963751 NA -0.086605408  0.011250699 9.644087e-01
# 28      B.1.351 - AY.5  0.0063270444 0.025076443 NA -0.042821881  0.055475970 1.000000e+00
# 29    B.1.351 - AY.5.2 -0.1450131805 0.027091459 NA -0.198111464 -0.091914897 2.268302e-05
# 30      B.1.351 - AY.9 -0.0138280347 0.025118163 NA -0.063058729  0.035402660 9.999992e-01
# 31     B.1.351 - AY.10 -0.0171712975 0.025133576 NA -0.066432200  0.032089605 9.999898e-01
# 32     B.1.351 - AY.33  0.0259445068 0.025166817 NA -0.023381547  0.075270561 9.989754e-01
# 33     B.1.351 - AY.34  0.0008811704 0.025450285 NA -0.049000471  0.050762812 1.000000e+00
# 34     B.1.351 - AY.39 -0.0159152370 0.025119444 NA -0.065148443  0.033317969 9.999958e-01
# 35    B.1.351 - AY.4.2 -0.0099064805 0.033902378 NA -0.076353921  0.056540960 1.000000e+00
# 36     B.1.351 - other -0.0044476342 0.024820950 NA -0.053095801  0.044200533 1.000000e+00
# 37          P.1 - AY.4 -0.1786090632 0.022150293 NA -0.222022839 -0.135195287 8.885781e-12
# 38          P.1 - AY.5 -0.1346046645 0.022263989 NA -0.178241282 -0.090968047 7.384530e-07
# 39        P.1 - AY.5.2 -0.2859448894 0.024514880 NA -0.333993171 -0.237896607 4.263256e-14
# 40          P.1 - AY.9 -0.1547597436 0.022308635 NA -0.198483865 -0.111035622 6.132607e-09
# 41         P.1 - AY.10 -0.1581030064 0.022335092 NA -0.201878983 -0.114327030 2.769076e-09
# 42         P.1 - AY.33 -0.1149872022 0.022360671 NA -0.158813312 -0.071161093 6.038650e-05
# 43         P.1 - AY.34 -0.1400505385 0.022686161 NA -0.184514597 -0.095586480 3.814322e-07
# 44         P.1 - AY.39 -0.1568469459 0.022321790 NA -0.200596849 -0.113097042 3.713901e-09
# 45        P.1 - AY.4.2 -0.1508381895 0.031877240 NA -0.213316432 -0.088359946 3.729087e-04
# 46         P.1 - other -0.1453793432 0.021805732 NA -0.188117792 -0.102640894 2.730498e-08
# 47         AY.4 - AY.5  0.0440043986 0.003232353 NA  0.037669102  0.050339695 1.743050e-14
# 48       AY.4 - AY.5.2 -0.1073358262 0.010709447 NA -0.128325957 -0.086345695 1.488809e-13
# 49         AY.4 - AY.9  0.0238493196 0.003459835 NA  0.017068168  0.030630471 7.839681e-09
# 50        AY.4 - AY.10  0.0205060568 0.003552461 NA  0.013543361  0.027468753 2.952416e-06
# 51        AY.4 - AY.33  0.0636218610 0.003814422 NA  0.056145732  0.071097990 1.731948e-14
# 52        AY.4 - AY.34  0.0385585246 0.005381499 NA  0.028010980  0.049106069 1.697284e-09
# 53        AY.4 - AY.39  0.0217621172 0.003446867 NA  0.015006382  0.028517852 1.826672e-07
# 54       AY.4 - AY.4.2  0.0277708737 0.023019269 NA -0.017346064  0.072887811 9.950562e-01 # AY.4.2 does not have a sign growth advantage over AY.4
# 55        AY.4 - other  0.0332297200 0.002625071 NA  0.028084675  0.038374765 1.865175e-14
# 56       AY.5 - AY.5.2 -0.1513402249 0.011018097 NA -0.172935298 -0.129745151 1.743050e-14
# 57         AY.5 - AY.9 -0.0201550791 0.004281422 NA -0.028546512 -0.011763647 4.136505e-04
# 58        AY.5 - AY.10 -0.0234983419 0.004343470 NA -0.032011386 -0.014985298 1.728261e-05
# 59        AY.5 - AY.33  0.0196174624 0.004573723 NA  0.010653130  0.028581795 2.289154e-03
# 60        AY.5 - AY.34 -0.0054458740 0.005935604 NA -0.017079444  0.006187696 9.997047e-01
# 61        AY.5 - AY.39 -0.0222422814 0.004261424 NA -0.030594518 -0.013890045 4.233111e-05
# 62       AY.5 - AY.4.2 -0.0162335249 0.023178564 NA -0.061662675  0.029195626 9.999863e-01
# 63        AY.5 - other -0.0107746786 0.003579865 NA -0.017791084 -0.003758273 1.444566e-01
# 64       AY.5.2 - AY.9  0.1311851458 0.011068950 NA  0.109490402  0.152879890 3.463896e-14
# 65      AY.5.2 - AY.10  0.1278418830 0.011107605 NA  0.106071376  0.149612390 5.362377e-14
# 66      AY.5.2 - AY.33  0.1709576873 0.011199003 NA  0.149008045  0.192907330 1.731948e-14
# 67      AY.5.2 - AY.34  0.1458943509 0.011826846 NA  0.122714159  0.169074543 2.131628e-14
# 68      AY.5.2 - AY.39  0.1290979435 0.011067463 NA  0.107406115  0.150789772 4.263256e-14
# 69     AY.5.2 - AY.4.2  0.1351067000 0.025216518 NA  0.085683233  0.184530167 2.213869e-05
# 70      AY.5.2 - other  0.1405655462 0.010856192 NA  0.119287801  0.161843292 1.776357e-14
# 71        AY.9 - AY.10 -0.0033432628 0.004503602 NA -0.012170160  0.005483635 9.999729e-01
# 72        AY.9 - AY.33  0.0397725415 0.004711624 NA  0.030537929  0.049007154 1.033063e-12
# 73        AY.9 - AY.34  0.0147092051 0.006062922 NA  0.002826097  0.026592313 4.656974e-01
# 74        AY.9 - AY.39 -0.0020872023 0.004443776 NA -0.010796844  0.006622439 9.999999e-01
# 75       AY.9 - AY.4.2  0.0039215541 0.023199438 NA -0.041548509  0.049391617 1.000000e+00
# 76        AY.9 - other  0.0093804004 0.003860596 NA  0.001813772  0.016947029 4.630712e-01
# 77       AY.10 - AY.33  0.0431158043 0.004790698 NA  0.033726209  0.052505400 1.771916e-13
# 78       AY.10 - AY.34  0.0180524679 0.006104347 NA  0.006088167  0.030016768 1.639906e-01
# 79       AY.10 - AY.39  0.0012560605 0.004506217 NA -0.007575963  0.010088084 1.000000e+00
# 80      AY.10 - AY.4.2  0.0072648170 0.023218162 NA -0.038241945  0.052771579 1.000000e+00
# 81       AY.10 - other  0.0127236633 0.003938751 NA  0.005003853  0.020443474 8.125549e-02
# 82       AY.33 - AY.34 -0.0250633364 0.006258416 NA -0.037329606 -0.012797067 6.704462e-03
# 83       AY.33 - AY.39 -0.0418597438 0.004730570 NA -0.051131491 -0.032587996 2.178258e-13
# 84      AY.33 - AY.4.2 -0.0358509873 0.023266170 NA -0.081451843  0.009749868 9.580832e-01
# 85       AY.33 - other -0.0303921410 0.004173875 NA -0.038572786 -0.022211496 8.724554e-10
# 86       AY.34 - AY.39 -0.0167964074 0.006058536 NA -0.028670920 -0.004921895 2.481522e-01
# 87      AY.34 - AY.4.2 -0.0107876509 0.023566190 NA -0.056976534  0.035401232 9.999999e-01
# 88       AY.34 - other -0.0053288046 0.005626292 NA -0.016356134  0.005698524 9.995830e-01
# 89      AY.39 - AY.4.2  0.0060087565 0.023199307 NA -0.039461051  0.051478563 1.000000e+00
# 90       AY.39 - other  0.0114676028 0.003838708 NA  0.003943874  0.018991332 1.525806e-01
# 91      AY.4.2 - other  0.0054588463 0.023096017 NA -0.039808514  0.050726207 1.000000e+00


# estimated proportion of different LINEAGES among new lab diagnoses in Belgium today
today # "2021-10-25"
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
#      LINEAGE         prob           SE df     asymp.LCL    asymp.UCL
# 1  B.1.617.2 4.821869e-01 1.438158e-02 NA  4.539995e-01 5.103743e-01
# 2    B.1.1.7 9.281472e-10 9.112552e-10 NA -8.578802e-10 2.714175e-09
# 3    B.1.351 2.618423e-05 5.146175e-05 NA -7.467896e-05 1.270474e-04
# 4        P.1 8.896540e-09 1.547357e-08 NA -2.143109e-08 3.922417e-08
# 5       AY.4 2.702892e-01 1.302664e-02 NA  2.447575e-01 2.958209e-01
# 6       AY.5 1.450725e-02 1.989080e-03 NA  1.060873e-02 1.840578e-02
# 7     AY.5.2 9.045981e-02 1.636622e-02 NA  5.838260e-02 1.225370e-01
# 8       AY.9 2.802057e-02 3.716797e-03 NA  2.073578e-02 3.530535e-02
# 9      AY.10 2.341439e-02 3.191549e-03 NA  1.715907e-02 2.966971e-02
# 10     AY.33 1.150431e-02 1.670597e-03 NA  8.230002e-03 1.477862e-02
# 11     AY.34 9.215505e-03 2.000468e-03 NA  5.294660e-03 1.313635e-02
# 12     AY.39 3.920751e-02 5.052464e-03 NA  2.930486e-02 4.911015e-02
# 13    AY.4.2 1.226491e-02 5.901764e-03 NA  6.976615e-04 2.383215e-02
# 14     other 1.890342e-02 2.229312e-03 NA  1.453405e-02 2.327279e-02

# estimates proportion of lab diagnoses that is B.1.617.2 by province
multinom_preds_today_byprovince = data.frame(emmeans(fit3_belgium_multi, ~ LINEAGE|1, by="province",
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_delta_today_byprovince = multinom_preds_today_byprovince[multinom_preds_today_byprovince$LINEAGE=="B.1.617.2",]
multinom_delta_today_byprovince = multinom_delta_today_byprovince[order(multinom_delta_today_byprovince$prob, decreasing=T),]
multinom_delta_today_byprovince

# reorder provinces by incidence of delta
levels_PROVINCES = as.character(multinom_delta_today_byprovince$province[order(multinom_delta_today_byprovince$prob, decreasing=T)])
data_agbyweekregion2$province = factor(data_agbyweekregion2$province, levels=levels_PROVINCES)

# redo plot of raw data by province using this order
muller_belgiumbyprovince_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) +
  facet_wrap(~ province, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
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
fit_belgium_multi_predsbyprovince$LINEAGE = factor(fit_belgium_multi_predsbyprovince$LINEAGE, levels=levels_LINEAGE2)
fit_belgium_multi_predsbyprovince$province = factor(fit_belgium_multi_predsbyprovince$province, levels=levels_PROVINCES)

fit_belgium_multi_preds = data.frame(emmeans(fit3_belgium_multi, 
                                           ~ LINEAGE,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                           mode="prob", df=NA))
fit_belgium_multi_preds$collection_date = as.Date(fit_belgium_multi_preds$DATE_NUM, origin="1970-01-01")
fit_belgium_multi_preds$LINEAGE = factor(fit_belgium_multi_preds$LINEAGE, levels=levels_LINEAGE2) 

muller_belgium_mfit = ggplot(data=fit_belgium_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE)) + 
  # facet_wrap(~ province, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
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
  scale_fill_manual("", values=lineage_cols2) +
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
fit_belgium_multi_preds2$LINEAGE = factor(fit_belgium_multi_preds2$LINEAGE, levels=levels_LINEAGE2)
levels(fit_belgium_multi_preds2$LINEAGE)

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
  coord_cartesian(xlim=c(as.Date("2021-02-01"),NA), ylim=c(0.005, 0.95), expand=c(0,0))
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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\belgium_multinom fit_response scale.pdf"), width=8, height=6)


# BY PROVINCE

# on logit scale:

fit_belgium_multi_preds3 = fit_belgium_multi_predsbyprovince
ymin = 0.001
ymax = 0.998
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
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweekregion2,
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
  coord_cartesian(xlim=c(as.Date("2021-02-01"),NA), ylim=c(0.005, 0.95), expand=c(0,0))
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
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweekregion2,
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
