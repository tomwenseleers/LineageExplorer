# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN SELECTED EUROPEAN COUNTRIES (GISAID records)
# T. Wenseleers
# last update 14 OCTOBER 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)

today = as.Date(Sys.time()) 
# today = as.Date("2021-10-14")
today_num = as.numeric(today)
plotdir = "global_GISAID_records"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID json metadata 
# get JSON file using 
# https://github.com/tomwenseleers/pyro-cov/blob/master/pull_gisaid.sh
# TO DO: maybe code that in R using RCurl or curl ?
# @Carl: I have a username & password (https://www.epicov.org/epi3/3p/kuleuven/export/provision.json.xz)
library(jsonlite)
GISAID_json = jsonlite::stream_in(gzfile(".//data//GISAID_json//provision.json.xz")) # takes ca half an hour to load on my laptop
nrow(GISAID_json) # data frame with 4288504  rows
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
nrow(GISAID_json) # 11800835 rows

# parse dates & remove records with invalid dates
library(stringr)
date_isvalid = sapply(GISAID_json$covv_collection_date, function (s) str_count(s, pattern = "-")==2)
GISAID_json = GISAID_json[date_isvalid,]
GISAID_json$date = as.Date(fast_strptime(GISAID_json$covv_collection_date, "%Y-%m-%d"))
GISAID_json = GISAID_json[!is.na(GISAID_json$date),]
nrow(GISAID_json) # 11492394

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


# MULLER PLOT (RAW DATA)
library(scales)
n = length(levels(GISAID_json$LINEAGE))
lineage_cols = hcl(h = seq(260, 0, length = n), l = 65, c = 200)
lineage_cols[which(levels(GISAID_json$LINEAGE)=="A")] = "red" # "grey40"
lineage_cols[which(levels(GISAID_json$LINEAGE)=="A*")] = "red4" # "grey40"
lineage_cols[which(levels(GISAID_json$LINEAGE)=="B")] = "grey50" # muted("blue",l=70,c=150)
lineage_cols[which(levels(GISAID_json$LINEAGE)=="B*")] = "grey30" # muted("blue",l=70,c=150)
lineage_cols[which(levels(GISAID_json$LINEAGE)=="B* (+S:D614G)")] = "#0085FF" # muted("blue",l=70,c=250)
lineage_cols[which(levels(GISAID_json$LINEAGE)=="Alpha")] = "blue"
lineage_cols[which(levels(GISAID_json$LINEAGE)=="Beta")] = muted("green", l=55, c=150) # "#9A9D00"
lineage_cols[which(levels(GISAID_json$LINEAGE)=="Gamma")] = "orange" # "cyan3"
lineage_cols[which(levels(GISAID_json$LINEAGE)=="Delta")] = "magenta"
lineage_cols[which(levels(GISAID_json$LINEAGE)=="Lambda")] = muted("magenta", l=25, c=110)
lineage_cols[which(levels(GISAID_json$LINEAGE)=="Mu")] = muted("magenta", l=15, c=50)

library(RColorBrewer)
lineage_cols = colorRampPalette( c("red","orange","yellow3","green4","blue2","magenta",muted("magenta", l=15, c=30)) )(n)
lineage_cols = muted(rainbow(n, start=0.3), c=150, l=50)
lineage_cols[9] = lineage_cols[10]
lineage_cols[10] = muted(rainbow(n, start=0.3), c=150, l=50)[9]
lineage_cols = colorRampPalette( c("green3","springgreen4","royalblue","magenta","red","gold"), bias=0.9, interpolate="spline" )(n)
# muted("magenta", l=15, c=30)

# Muller plot, overall by continent
muller_raw_continent = ggplot(data=data_agbyweek_bycontinent, 
                              aes(x=collection_date, y=count, group=LINEAGE)) + 
  facet_wrap(~ continent) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2019-12-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN (GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_json$Week)))
muller_raw_continent

ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_raw data_by continent.png"), width=12, height=6)

# muller plot, by country
muller_raw_country = ggplot(data=data_agbyweek_bycontinent_country, aes(x=collection_date, y=count, group=LINEAGE)) + 
  facet_wrap(~ country) +
  geom_col(aes(colour=NULL, fill=LINEAGE), width=I(7), position="fill") +
  # geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2019-12-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN (GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_json$Week)))
muller_raw_country

ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_raw data_by country.png"), width=20, height=12)


# multinomial fits

library(nnet)
library(splines)
set.seed(1)
fit1_multi = nnet::multinom(LINEAGE ~ scale(collection_date_num)*continent+country, weights=count, 
                            data=data_agbyweek_bycontinent_country, maxit=1000, MaxNWts=10000)
BIC(fit1_multi) 
# 3939555


# TO DO: try fit on data_agbyweek_bycontinent_country using mblogit & country coded 
# as random intercept in multinomial mixed model

# growth rate advantage compared to the original type A (difference in growth rate per day) 
emtr = emtrends(fit1_multi, trt.vs.ctrl ~ LINEAGE, by=c("collection_date_num","continent"),  
                   var="collection_date_num",  mode="latent",
                   at=list(collection_date_num=today_num))
delta_r = data.frame(confint(emtr, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtr$contrasts)$p.value)
delta_r
# contrast     estimate           SE df    asymp.LCL    asymp.UCL      p.value
# 1              B - A -0.002256622 8.417457e-05 NA -0.002421602 -0.002091643 6.391887e-12
# 2 (B (+S:D614G)) - A  0.007193565 6.440158e-05 NA  0.007067340  0.007319790 6.391887e-12
# 3          Alpha - A  0.033571472 7.615868e-05 NA  0.033422204  0.033720740 6.391887e-12
# 4           Beta - A  0.028832889 1.328787e-04 NA  0.028572451  0.029093326 6.391887e-12
# 5          Gamma - A  0.044544851 1.173837e-04 NA  0.044314783  0.044774918 6.391887e-12
# 6          Delta - A  0.098587921 1.121696e-04 NA  0.098368073  0.098807769 6.391887e-12
# 7         Lambda - A  0.047636271 2.569825e-04 NA  0.047132595  0.048139948 6.391887e-12
# 8             Mu - A  0.063764933 2.385987e-04 NA  0.063297289  0.064232578 6.391887e-12


# fitted prop of different LINEAGES today
multinom_preds_today_avg = data.frame(emmeans(fit1_multi, ~ LINEAGE|1,
                                              at=list(collection_date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE        prob           SE df     asymp.LCL   asymp.UCL
# 1   B.1.1.7 0.008421286 0.0018579564 NA  0.0047797585 0.012062814
# 2    A.23.1 0.001045963 0.0004335694 NA  0.0001961824 0.001895743
# 3       B.1 0.001286566 0.0002919790 NA  0.0007142979 0.001858835
# 4 B.1.1.318 0.056542496 0.0113195641 NA  0.0343565580 0.078728434
# 5   B.1.351 0.010441396 0.0014117131 NA  0.0076744890 0.013208302
# 6   B.1.525 0.003671328 0.0012003110 NA  0.0013187618 0.006023894
# 7 B.1.617.1 0.002997387 0.0020830367 NA -0.0010852897 0.007080064
# 8 B.1.617.2 0.914131630 0.0133441025 NA  0.8879776701 0.940285591
# 9     other 0.001461948 0.0003504848 NA  0.0007750100 0.002148885



# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2019-12-01"))
extrapolate=30
date.to = today_num+extrapolate

# multinomial model predictions by country (fastest, but no confidence intervals)
predgrid = expand.grid(list(collection_date_num=seq(date.from, date.to), 
                            continent=levels_continents))
fitmulti_preds = data.frame(predgrid, as.data.frame(predict(fit1_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fitmulti_preds = gather(fitmulti_preds, LINEAGE, prob, all_of(levels_LINEAGE), factor_key=TRUE)
fitmulti_preds$collection_date = as.Date(fitmulti_preds$collection_date_num, origin="1970-01-01")
fitmulti_preds$LINEAGE = factor(fitmulti_preds$LINEAGE, levels=levels_LINEAGE) 

muller_mfit_continent = ggplot(data=fitmulti_preds, 
                                   aes(x=collection_date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ continent) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  annotate("rect", xmin=max(as.numeric(GISAID_json$date))+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2019-12-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN (GISAID data, multinomial fit)")
muller_mfit_continent

ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plot_multinom fit_by continent.png"), width=10, height=6)


library(ggpubr)
ggarrange(muller_raw_continent + coord_cartesian(xlim=c(as.Date("2019-12-01"),today+extrapolate))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_mfit_continent+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\muller plots multipanel_raw data plus multinom fit.png"), width=10, height=10)




# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions with confidence intervals by continent
fitmulti_preds_bycontinent_withCI = data.frame(emmeans(fit1_multi,
                                                        ~ LINEAGE,
                                                        by=c("collection_date_num","continent"),
                                                        at=list(collection_date_num=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                        mode="prob", df=NA))
fitmulti_preds_bycontinent_withCI$collection_date = as.Date(fitmulti_preds_bycontinent_withCI$DATE_NUM, origin="1970-01-01")
fitmulti_preds_bycontinent_withCI$LINEAGE = factor(fitmulti_preds_bycontinent_withCI$LINEAGE, levels=levels_LINEAGE)
fitmulti_preds_bycontinent_withCI$continent = factor(fitmulti_preds_bycontinent_withCI$continent, levels=levels_countries)
fitmulti_preds3 = fitmulti_preds_bycontinent_withCI

# on logit scale:

ymin = 0.001
ymax = 0.999
fitmulti_preds3$asymp.LCL[fitmulti_preds3$asymp.LCL<ymin] = ymin
fitmulti_preds3$asymp.UCL[fitmulti_preds3$asymp.UCL<ymin] = ymin
fitmulti_preds3$asymp.UCL[fitmulti_preds3$asymp.UCL>ymax] = ymax
fitmulti_preds3$prob[fitmulti_preds3$prob<ymin] = ymin

plotbycontinent_mfit_logit = qplot(data=fitmulti_preds3, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~ continent) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  geom_point(data=data_agbyweek_bycontinent2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plotbycontinent_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by continent_multinom fit_logit scale.png"), width=10, height=6)


# on response scale:
plotbycontinent_mfit = qplot(data=fitmulti_preds3, x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~ continent) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-11-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  geom_point(data=data_agbyweek_bycontinent2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plotbycontinent_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by continent_multinom fit_response scale.png"), width=10, height=6)





# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot_DRC = as.data.frame(get_national_data(countries = levels_countries[grepl("Congo",levels_countries)], source="WHO"))
cases_tot_DRC = cases_tot_DRC[,c("date","country","cases_new","hosp_new","deaths_new","tested_new")]
cases_tot = as.data.frame(get_national_data(countries = levels_countries, source="Google"))
cases_tot = cases_tot[,c("date","country","cases_new","hosp_new","deaths_new","tested_new")]
cases_tot = rbind(cases_tot_DRC, cases_tot)
cases_tot = cases_tot[cases_tot$date>=as.Date("2021-01-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
cases_tot$country = factor(cases_tot$country, levels=unique(cases_tot$country),
                              labels=gsub("Congo - Kinshasa", "Democratic Republic of the Congo", unique(cases_tot$country)))
cases_tot$country = factor(cases_tot$country, levels=levels_countries)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
library(mgcv)
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=7, fx=F, by=country) + country +
                  WEEKDAY # +
                # s(tested_new, bs="cs", k=8, fx=F, by=region)
                ,
                family=poisson(log), data=cases_tot,
) 
BIC(fit_cases)

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fitmulti_preds$totcases = cases_tot$cases_new[match(interaction(fitmulti_preds$DATE_NUM,fitmulti_preds$country),
                                                                     interaction(cases_tot$DATE_NUM,cases_tot$country))]
fitmulti_preds$cases = fitmulti_preds$totcases * fitmulti_preds$prob
fitmulti_preds$cases[fitmulti_preds$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM|country, at=list(DATE_NUM=seq(date.from, date.to, by=0.5)
), type="response"))
fitmulti_preds$smoothed_totcases = cases_emmeans$rate[match(interaction(fitmulti_preds$DATE_NUM,fitmulti_preds$country),
                                                                             interaction(cases_emmeans$DATE_NUM,cases_emmeans$country))]
fitmulti_preds$smoothed_cases = fitmulti_preds$smoothed_totcases * fitmulti_preds$prob
fitmulti_preds$smoothed_cases[fitmulti_preds$smoothed_cases<=0.001] = NA

ggplot(data=fitmulti_preds, 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ country, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fitmulti_preds, 
       aes(x=collection_date, y=cases+1, group=LINEAGE)) + 
  facet_wrap(~ country, scale="free") +
  geom_col(data=fitmulti_preds[fitmulti_preds$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_col(data=fitmulti_preds[fitmulti_preds$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  geom_col(data=fitmulti_preds[fitmulti_preds$LINEAGE=="B.1.351",], aes(lwd=I(1.2), colour=NULL), fill="#9A9D00", position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fitmulti_preds,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE)) +
  facet_wrap(~ country, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fitmulti_preds, 
       aes(x=collection_date, y=smoothed_cases+1, group=LINEAGE)) + 
  facet_wrap(~ country, scale="free") +
  geom_area(data=fitmulti_preds[fitmulti_preds$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fitmulti_preds[fitmulti_preds$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  geom_area(data=fitmulti_preds[fitmulti_preds$LINEAGE=="B.1.351",], aes(lwd=I(1.2), colour=NULL), fill="#9A9D00", position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fitmulti_preds, 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE, group=LINEAGE)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fitmulti_preds,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE)) +
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE, group=LINEAGE)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),today), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit smoothed case data.png"), width=8, height=6)






# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####
# TO DO : need to finish this part

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}


# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, by="region", var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to, by=3)
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[3] = "r"
colnames(avg_r_cases)[6] = "r_LOWER"
colnames(avg_r_cases)[7] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  facet_wrap(~ region) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN THE UK AT MOMENT OF INFECTION BASED ON NEW CASES\n(data gov.uk)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants4 = do.call(rbind, lapply(levels_REGION, function(region) { do.call(rbind, 
                                                                                        lapply(seq(date.from,
                                                                                                   date.to, by=3), 
                                                                                               function (d) { 
                                                                                                 wt = as.data.frame(emmeans(fit4_cogukp2_multi, ~ LINEAGE , by="REGION", 
                                                                                                                            at=list(DATE_NUM=d, REGION=region), type="response"))$prob   # important: these should sum to 1
                                                                                                 # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) 
                                                                                                 # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
                                                                                                 cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                                                                 names(cons) = seq_along(cons)
                                                                                                 EMT = emtrends(fit4_cogukp2_multi,  ~ LINEAGE , by=c("DATE_NUM", "REGION"),
                                                                                                                var="DATE_NUM", mode="latent",
                                                                                                                at=list(DATE_NUM=d, REGION=region))
                                                                                                 out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                                                                 # sum(out$estimate*wt) # should sum to zero
                                                                                                 return(out) } )) } ))
above_avg_r_variants = above_avg_r_variants4
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE), 
                                       labels=levels(data_agbyweekregion1$LINEAGE))
above_avg_r_variants$LINEAGE = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2020-09-01" "2021-07-31"
# average growth rate of all lineages calculated from case nrs
above_avg_r_variants$avg_r = avg_r_cases$r[match(interaction(above_avg_r_variants$collection_date,above_avg_r_variants$REGION),
                                                 interaction(avg_r_cases$DATE,avg_r_cases$region))]  
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                REGION=avg_r_cases$region,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE,
                LINEAGE="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$LINEAGE = factor(above_avg_r_variants$LINEAGE, levels=c(levels_LINEAGE_plot,"avg"))
above_avg_r_variants$prob = fitmulti_preds$prob[match(interaction(round(above_avg_r_variants$DATE_NUM),
                                                                                   as.character(above_avg_r_variants$LINEAGE),
                                                                                   as.character(above_avg_r_variants$REGION)),
                                                                       interaction(round(fitmulti_preds$DATE_NUM),
                                                                                   as.character(fitmulti_preds$LINEAGE),
                                                                                   as.character(fitmulti_preds$country)))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 2
ymin = 1/2
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$LINEAGE %in% c("other"))|(above_avg_r_variants2$collection_date>=max(cases_tot$date))),], # |above_avg_r_variants2$collection_date>max(cases_tot$DATE)
      x=collection_date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE, fill=LINEAGE, alpha=I(0.5),
      group=LINEAGE, linetype=I(0)) +
  facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN THE UK\nAT MOMENT OF INFECTION\n(based on gov.uk case data & multinomial fit to\nCOG-UK lineage frequencies)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-11-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  scale_colour_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)


