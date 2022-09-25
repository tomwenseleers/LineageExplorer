# IMPORT GISAID json METADATA 
# T. Wenseleers, 20 JULY 2022

# to download latest GISAID JSON stream execute
# source(".//download_GISAID_json.R")

# PS best to run this on a workstation with a lot of memory
# takes ca half an hour to parse
library(jsonlite)
message("Parsing JSON GISAID records...")
GISAID_json = jsonlite::stream_in(gzfile(".//data//GISAID_json//provision.json.xz")) 
saveRDS(GISAID_json, ".//data/GISAID_json//GISAID_json.rds")
# GISAID_json = readRDS(".//data/GISAID_json//GISAID_json.rds")

nrow(GISAID_json) # data frame with 12025647 rows

sum(missingIDs %in% GISAID_json$covv_accession_id) # 113

# only keep human samples
GISAID_json = GISAID_json[GISAID_json$covv_host=="Human",] 
nrow(GISAID_json) # data frame with 11853650  rows
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
GISAID_json = mutate_at(GISAID_json, 
                        "covv_sampling_strategy", .funs=toupper)
GISAID_json = mutate_at(GISAID_json, 
                        "covv_add_host_info", .funs=toupper)
GISAID_json = mutate_at(GISAID_json, 
                        "covv_add_location", .funs=toupper)

# remove travel-related cases, active surveillance, surge testing, targeted S dropout sequencing, etc
filterout_actsurv = FALSE

if (filterout_actsurv) {
filter_sampling_strategy = "ACTIVE|S GENE|COVIGEN|S-GENE|SUSPECT|LONGITUDINAL|ASSAY|DROPOUT|CLUSTER|OUTBREAK|LONGITUDINAL|S GENE|SAME-PATIENT|TRAVEL|VOC|VARIANT|CONTACT|TRACING|PCR TEST|RETURNING|SGTF|TRIAL|FAMILY|DROPOUT|TAQPATH|WEEKLY|NON-RANDOM|SAME PATIENT|TIME-SERIES|TARGET|QUARANTINE|PASSAGE"
filter_add_host_info = "ACTIVE|S GENE|COVIGEN|S-GENE|SUSPECT|LONGITUDINAL|ASSAY|DROPOUT|CLUSTER|OUTBREAK|LONGITUDINAL|S GENE|SAME-PATIENT|TRAVEL|VOC|VARIANT|CONTACT|TRACING|PCR TEST|RETURNING|SGTF|TRIAL|FAMILY|DROPOUT|TAQPATH|WEEKLY|NON-RANDOM|SAME PATIENT|TIME-SERIES|TARGET|QUARANTINE|PASSAGE|PASSENGER|INFECTED IN|HOLIDAY|CONTACT|CAME FROM|RETURN|FROM|VISIT|FOREIGN|SKIING|IMPORT|RELATED WITH|ENTERED FROM|CASE FROM"
filter_add_location = "ACTIVE|S GENE|COVIGEN|S-GENE|SUSPECT|LONGITUDINAL|ASSAY|DROPOUT|CLUSTER|OUTBREAK|LONGITUDINAL|S GENE|SAME-PATIENT|TRAVEL|VOC|VARIANT|CONTACT|TRACING|PCR TEST|RETURNING|SGTF|TRIAL|FAMILY|DROPOUT|TAQPATH|WEEKLY|NON-RANDOM|SAME PATIENT|TIME-SERIES|TARGET|QUARANTINE|PASSAGE|PASSENGER|INFECTED IN|HOLIDAY|CONTACT|CAME FROM|RETURN|FROM|VISIT|FOREIGN|SKIING|IMPORT|RELATED WITH|ENTERED FROM|CASE FROM|HISTOR|SAMPLED AT|TOURIST|VISIT"

GISAID_json = GISAID_json[!grepl(filter_sampling_strategy,
                                 GISAID_json$covv_sampling_strategy),]
nrow(GISAID_json) # 11658421 rows
GISAID_json = GISAID_json[!grepl(filter_add_host_info,
                                 GISAID_json$covv_add_host_info),]
nrow(GISAID_json) # 11641947 rows
GISAID_json = GISAID_json[!grepl(filter_add_location,
                                 GISAID_json$covv_add_location),]
nrow(GISAID_json) # 11619356 rows
# 11853650 - 11619356 # 234294 records removed (2%)
}

# parse dates & remove records with invalid dates
library(stringr)
date_isvalid = (str_count(GISAID_json$covv_collection_date, 
                         pattern = "-")==2)
GISAID_json = GISAID_json[date_isvalid,]
library(lubridate)
GISAID_json$date = as.Date(fast_strptime(GISAID_json$covv_collection_date, "%Y-%m-%d"))
GISAID_json = GISAID_json[!is.na(GISAID_json$date),]
nrow(GISAID_json) # 11311172

# add numeric version of date, start of week & week & year

range(GISAID_json$date) # "2019-12-05" "2022-07-08"

GISAID_json$Week = lubridate::week(GISAID_json$date)
GISAID_json$Year = lubridate::year(GISAID_json$date)
GISAID_json$Year_Week = interaction(GISAID_json$Year,GISAID_json$Week)
GISAID_json$week_startingdate = as.Date(fast_strptime(as.character(cut(GISAID_json$date, "week")), "%Y-%m-%d"))
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

# remove sequences not assigned to pango lineage
removeunassigned = FALSE
if (removeunassigned) GISAID_json = GISAID_json[!(GISAID_json$covv_lineage=="None"|GISAID_json$covv_lineage==""|is.na(GISAID_json$covv_lineage)),]
nrow(GISAID_json)

# define variant names
sel_target_VOC = "Omicron (BA.2.75)"
sel_reference_VOC = "Omicron (BA.5)"
levels_VARIANTS = c(sel_reference_VOC, "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", "Omicron (BA.2.74)", "Omicron (BA.2.76)", sel_target_VOC)
levels_VARIANTS_plot = c("Other", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", "Omicron (BA.5)", "Omicron (BA.2.74)", "Omicron (BA.2.76)", "Omicron (BA.2.75)")

# define variant colours
n = length(levels_VARIANTS)
lineage_cols = hcl(h = seq(0, 180, length = n), l = 60, c = 180)
lineage_cols[which(levels_VARIANTS=="Alpha")] = "#0085FF"
lineage_cols[which(levels_VARIANTS=="Beta")] = "green4"
lineage_cols[which(levels_VARIANTS=="Delta")] = "mediumorchid"
# lineage_cols[which(levels_VARIANTS=="C.1.2")] = "darkorange"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.1)")] = "red" # "magenta"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2)")] = "red3"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.3)")] = "red4" 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.4)")] = "darkorange" 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.5)")] = "darkorange3" 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2.74)")] = "black"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2.76)")] = muted(muted("magenta")) 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2.75)")] = "magenta" 
lineage_cols[which(levels_VARIANTS=="Other")] = "grey65"

lineage_cols_plot = lineage_cols[match(levels_VARIANTS_plot,levels_VARIANTS)]

# code variants from pango lineages & AA substitutions present
GISAID_json$pango_lineage = GISAID_json$covv_lineage
GISAID_json$aa_substitutions = GISAID_json$covsurver_existingmutlist

library(dplyr)
# PS this is slow - optimize this, maybe use 
# multicore version multidplyr?
GISAID_json$variant = case_when(
  grepl("B.1.617.2", GISAID_json$pango_lineage, fixed=T) | grepl("AY", GISAID_json$pango_lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", GISAID_json$pango_lineage) ~ "Alpha",
  grepl("B.1.351", GISAID_json$pango_lineage, fixed=T) ~ "Beta",
  (grepl("^BA\\.1$|BA\\.1\\.", GISAID_json$pango_lineage)) ~ "Omicron (BA.1)",
  (grepl("^BA\\.3$|BA\\.3\\.", GISAID_json$pango_lineage)) ~ "Omicron (BA.3)",
  (((grepl("^BA\\.4",GISAID_json$pango_lineage)))|((grepl("BA.2",GISAID_json$pango_lineage))&
                                                ((grepl("L452R", GISAID_json$aa_substitutions)&
                                                    grepl("486V", GISAID_json$aa_substitutions)&
                                                    grepl("11F", GISAID_json$aa_substitutions)&
                                                    (!grepl("D3N",GISAID_json$aa_substitutions)) )))) ~ "Omicron (BA.4)",
  (((grepl("^BA\\.5",GISAID_json$pango_lineage)))|((grepl("BA.2",GISAID_json$pango_lineage))&
                                                ((grepl("L452R",GISAID_json$aa_substitutions)&
                                                    grepl("486V",GISAID_json$aa_substitutions)&
                                                    (!grepl("11F", GISAID_json$aa_substitutions))&
                                                    grepl("D3N",GISAID_json$aa_substitutions))))) ~ "Omicron (BA.5)",
  ((grepl("NSP3_S403L",GISAID_json$aa_substitutions)& 
    grepl("NSP8_N118S",GISAID_json$aa_substitutions))) ~ "Omicron (BA.2.75)",
  ((grepl("Spike_L452M",GISAID_json$aa_substitutions)&
    grepl("Spike_R346T",GISAID_json$aa_substitutions))) ~ "Omicron (BA.2.74)",
  ((grepl("Spike_Y248N",GISAID_json$aa_substitutions)&
   grepl("Spike_R346T",GISAID_json$aa_substitutions))) ~ "Omicron (BA.2.76)",
  (grepl("^BA\\.2",GISAID_json$pango_lineage)) ~ "Omicron (BA.2)", 
  T ~ "Other"
)

# see https://twitter.com/JosetteSchoenma/status/1545509992572272641/photo/1
# https://twitter.com/JosetteSchoenma/status/1546892094316445696/photo/1
# also: 
# BA.2.74 = S:L452M + S:R346T
# BA.2.76 = S:Y248N + S:R346T
# BA.2.77 = S:K356T + S:L452R + S:R493Q + S:D936H + S:E340K
# BA.2.78 = 
# BA.2.79 = S:N450D+S:Y248S
# or for BA.2.75 = E_T11A,NSP3_S403L+sample date since May (https://twitter.com/CorneliusRoemer/status/1547023464602734593)

sum(GISAID_json$variant=="Omicron (BA.2.75)") # 244
GISAID_json[GISAID_json$variant=="Omicron (BA.2.75)",]
df=as.data.frame(table(GISAID_json[GISAID_json$variant=="Omicron (BA.2.75)",
                  "covv_location"]))
df

df=as.data.frame(table(GISAID_json[GISAID_json$variant=="Omicron (BA.2.75)",
                  "country"]))
df=df[df$Freq!=0,]
colnames(df)=c("country","BA.2.75")
df[order(df$BA.2.75,decreasing=T),]

sum(grepl("Maharashtra",GISAID_json[GISAID_json$variant=="Omicron (BA.2.75)","covv_location"])) # 84
# most BA.2.75 seqs from Maharashtra (n=84) are from Aurangabad (n=42) or Nagpur (n=36)

sum(grepl("West Bengal",GISAID_json[GISAID_json$variant=="Omicron (BA.2.75)","covv_location"])) # 84
# most BA.2.75 seqs from West Bengal (n=30) are from Kolkata (n=8)

table(GISAID_json$continent, GISAID_json$variant)

# define X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))




# subset to desired date range
start_date = "2020-06-01"
end_date = today
GISAID_json = GISAID_json[GISAID_json$date>=as.Date(start_date)&
                            GISAID_json$date<=as.Date(end_date),]
nrow(GISAID_json) # 11346346




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
