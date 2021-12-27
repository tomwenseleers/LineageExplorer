library(readr)

# import GISAID metadata (file version metadata_tsv_2021_12_27.tar.xz)
message("Parsing GISAID metadata...")
GISAID = read_tsv(".//data//GISAID//metadata.tsv", col_types = cols(.default = "c"))
GISAID = as.data.frame(GISAID)
# only keep human samples
GISAID = GISAID[GISAID$Host=="Human",] 
GISAID = GISAID[!is.na(GISAID$Type),] 
nrow(GISAID) # data frame with 5967874  rows
names(GISAID)
head(GISAID$Location)
head(GISAID$"Collection date")
head(GISAID$"Pango lineage")
unique(GISAID$"Virus name")
unique(GISAID$"Type")
unique(GISAID$"Pango lineage")
unique(GISAID$"Pangolin version")
unique(GISAID$"Variant")
unique(GISAID$"AA Substitutions")
head(GISAID$covsurver_prot_mutations)
head(GISAID$covsurver_uniquemutlist)
unique(GISAID$covv_sampling_strategy)

# parse dates & remove records with invalid dates
library(stringr)
date_isvalid = sapply(GISAID$"Collection date", function (s) str_count(s, pattern = "-")==2)
GISAID = GISAID[date_isvalid,]
GISAID$date = as.Date(GISAID$"Collection date") 
GISAID = GISAID[!is.na(GISAID$date),]
nrow(GISAID) # 6230709

# only keep complete genomes
GISAID = GISAID[GISAID$"Is complete?"=="True",]
nrow(GISAID) # 6230709

# add numeric version of date, week midpoint & week & year
range(GISAID$date) # "2019-12-24" "2021-12-23"
# GISAID[GISAID$date==min(GISAID$date),] # first available Wuhan genome
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
library(lubridate)
GISAID$week_midpoint = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
names(GISAID)

# parse continent / country / location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID$Location), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
unique(loc[,1]) # continent
loc[,2][loc[,3]=="Scotland"] = "Scotland"
loc[,2][loc[,3]=="England"] = "England"
loc[,2][loc[,3]=="Wales"] = "Wales"
loc[,2][loc[,3]=="Northern Ireland"] = "Northern Ireland"
loc[,2][loc[,2]=="United Kingdom"] = "England"
sort(unique(loc[,2])) # country

sort(unique(loc[,3])) # city or province

levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
GISAID$continent = factor(loc[,1], levels=levels_continents)
GISAID$country = factor(loc[,2])
levels_countries = levels(GISAID$country)
GISAID$location = factor(loc[,3])
levels_locations = levels(GISAID$location)


# recode pango lineages as desired
GISAID = GISAID[!(GISAID$"Pango lineage"=="None"|GISAID$"Pango lineage"==""|is.na(GISAID$"Pango lineage")),]
nrow(GISAID) # 6178953
unique(GISAID$"Pango lineage")
length(unique(GISAID$"Pango lineage")) # 1542 pango lineages

# code Omicron lineage
GISAID$Omicron = "Other"
GISAID$Omicron[grepl("Omicron",GISAID$Variant)] = "Omicron"

# code Delta lineage
GISAID$Delta = "Other"
GISAID$Delta[grepl("Delta",GISAID$Variant)] = "Delta"

# code sgtf (S gene target failure) based on presence of Spike_H69del / Spike_V70del
GISAID$sgtf = grepl("Spike_H69del",GISAID$"AA Substitutions", fixed=T)|grepl("Spike_V70del",GISAID$"AA Substitutions", fixed=T)
GISAID$sgtf[GISAID$Omicron=="Omicron"] = TRUE
GISAID$sgtf[GISAID$Omicron=="Delta"] = TRUE
GISAID$sgtf[GISAID$sgtf==TRUE] = "sgtf"
GISAID$sgtf[GISAID$sgtf==FALSE] = "non_sgtf"

GISAID_sgtf_isomicron = as.data.frame(ftable(formula=Omicron~country+date, data=GISAID[GISAID$sgtf=="sgtf",]))
library(tidyr)
# proportion of sgtf samples that are omicron in function of time by country & date to run logistic regression on
GISAID_sgtf_isomicron = spread(GISAID_sgtf_isomicron, Omicron, Freq) 
GISAID_sgtf_isomicron$date = as.Date(GISAID_sgtf_isomicron$date)
GISAID_sgtf_isomicron$DATE_NUM = as.numeric(GISAID_sgtf_isomicron$date)
saveRDS(GISAID_sgtf_isomicron, file=".//data//GISAID//GISAID_sgtf_isomicron.rds")

# # code main clades & VOCs only # TO DO: FINISH THIS PART
# unique(GISAID$Clade)
# unique(GISAID$Variant)
# GISAID$VOC = GISAID$Variant
# GISAID$VOC[GISAID$covv_lineage=="A"] = "A" 
# GISAID$VOC[GISAID$VOC==""&(grepl("S",GISAID$covv_clade)|grepl("A.",GISAID$covv_lineage,fixed=T))] = "A*" # = 19B, cf https://www.news-medical.net/health/Viral-Clades-of-SARS-CoV-2.aspx
# GISAID$VOC[GISAID$covv_lineage=="B"] = "B"
# GISAID$VOC[GISAID$VOC==""&(grepl("L|O|V|G",GISAID$covv_clade)|grepl("B.",GISAID$covv_lineage,fixed=T))&(!grepl("D614G",GISAID$covsurver_prot_mutations))] = "B*" # B* without D614G
# GISAID$VOC[GISAID$VOC==""&(grepl("G",GISAID$covv_clade)|grepl("B.",GISAID$covv_lineage,fixed=T)|grepl("AE.",GISAID$covv_lineage,fixed=T)|grepl("AH.",GISAID$covv_lineage,fixed=T)|grepl("AK.",GISAID$covv_lineage,fixed=T)|grepl("AN.",GISAID$covv_lineage,fixed=T)|grepl("AS.",GISAID$covv_lineage,fixed=T)|grepl("AU.",GISAID$covv_lineage,fixed=T)|grepl("AV.",GISAID$covv_lineage,fixed=T)|grepl("AZ.",GISAID$covv_lineage,fixed=T))&(grepl("D614G",GISAID$covsurver_prot_mutations))] = "B* (+S:D614G)" # B* with D614G
# GISAID$VOC[grepl("Alpha",GISAID$VOC)] = "Alpha"
# GISAID$VOC[grepl("Beta",GISAID$VOC)] = "Beta"
# GISAID$VOC[grepl("Gamma",GISAID$VOC)|grepl("P.1.",GISAID$covv_lineage,fixed=T)] = "Gamma"
# GISAID$VOC[grepl("Delta",GISAID$VOC)|grepl("AY.",GISAID$covv_lineage,fixed=T)] = "Delta"
# GISAID$VOC[grepl("Lambda",GISAID$VOC)] = "Lambda"
# GISAID$VOC[grepl("Mu",GISAID$VOC)] = "Mu"
# table(GISAID$VOC)
# table(GISAID$covv_lineage[GISAID$VOC==""])  
# # sum(GISAID$VOC=="") # 0 remaining unassigned
# # GISAID = GISAID[GISAID$VOC!="",] # we remove these samples; instead we could also recode them as category "other"
# nrow(GISAID) # 3981577
# 
# table(GISAID$covv_lineage[GISAID$covv_location=="Asia / China / Hubei / Wuhan"])
# 
# # use main clades/VOCs as my LINEAGE variable
# levels_LINEAGE = c("A", "A*", "B", "B*", "B* (+S:D614G)", "Alpha", "Beta", "Gamma", "Delta", "Lambda", "Mu")
# GISAID$LINEAGE = factor(GISAID$VOC, 
#                              levels=levels_LINEAGE)





