# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN INDIA BASED ON ANALYSIS OF GISAID GENOMIC EPIDEMIOLOGY METADATA
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
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata (file version metadata_2021-05-21_13-00.tsv.gz)
GISAID = read_tsv(".//data//GISAID_genomic_epidemiology//metadata.tsv", col_types = cols(.default = "c"))
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
unique(GISAID$host)
# [1] "Human"                               "Environment"                         "Rhinolophus shameli"                 "Rhinolophus malayanus"              
# [5] "Rhinolophus pusillus"                "Rhinolophus sinicus"                 "Rhinolophus stheno"                  "Rhinolophus affinis"                
# [9] "Felis catus"                         "Canis lupus familiaris"              "Gorilla gorilla gorilla"             "Mesocricetus auratus"               
# [13] "Prionailurus bengalensis euptilurus" "Panthera leo"                        "Mink"                                "Mustela putorius furo"              
# [17] "Chlorocebus sabaeus"                 "Mus musculus"                        "Mus musculus (BALB/c mice)"          "Manis javanica"                     
# [21] "Manis pentadactyla"                  "Panthera tigris jacksoni" 
GISAID[GISAID$host!="Human","strain"]
GISAID = GISAID[GISAID$host=="Human",]
GISAID = GISAID[GISAID$date>=as.Date("2019-09-01"),]
range(GISAID$date) # "2019-12-24" "2021-05-18"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 10 dec Maharashtra India

# GISAID = GISAID[grepl("2021-", GISAID$date),]
GISAID = GISAID[GISAID$date>=as.Date("2020-06-01"),]
sum(is.na(GISAID$purpose_of_sequencing)) == nrow(GISAID) # field purpose_of_sequencing left blank unfortunately
range(GISAID$date) # "2021-06-01" "2021-05-18"
nrow(GISAID) # 1446226
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
library(lubridate)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
colnames(GISAID)
unique(GISAID$region)
# "Asia"          "Europe"        "Africa"        "South America" "Oceania"       "North America"
unique(GISAID$country)
unique(GISAID$division) # = city or province or region, sometimes just country
unique(GISAID$location) # = city

length(unique(GISAID$country[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 52 countries
table(GISAID$pango_lineage[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])
# B.1.617 B.1.617.1 B.1.617.2 B.1.617.3 
# 1      2225      6319        79 

GISAID$pango_lineage[grepl("B.1.177",GISAID$pango_lineage,fixed=T)] = "B.1.177+"
GISAID$pango_lineage[grepl("B.1.36\\>",GISAID$pango_lineage)] = "B.1.36+"

sel_target_VOC = "B.1.617"
GISAID$LINEAGE1 = GISAID$pango_lineage
GISAID$LINEAGE2 = GISAID$pango_lineage
GISAID[grepl(sel_target_VOC, GISAID$LINEAGE1, fixed=TRUE),"LINEAGE1"] = paste0(sel_target_VOC,"+") # in LINEAGE1 we recode B.1.617.1,2&3 all as B.1.617+

table_country_lineage = as.data.frame(table(GISAID$country, GISAID$LINEAGE1))
colnames(table_country_lineage) = c("Country","Lineage","Count")
tblB1617 = table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage, fixed=T)&table_country_lineage$Count>10,]
tblB1617
#               Country  Lineage Count
# 149628      Australia B.1.617+    95
# 149632        Bahrain B.1.617+    19
# 149636        Belgium B.1.617+    62
# 149660        Denmark B.1.617+    90
# 149670         France B.1.617+    50
# 149675        Germany B.1.617+   260
# 149688          India B.1.617+  2337
# 149692        Ireland B.1.617+    90
# 149693         Israel B.1.617+    36
# 149694          Italy B.1.617+    57
# 149696          Japan B.1.617+   141
# 149723    Netherlands B.1.617+    31
# 149724    New Zealand B.1.617+    15
# 149737         Poland B.1.617+    26
# 149750      Singapore B.1.617+   156
# 149757          Spain B.1.617+    42
# 149760         Sweden B.1.617+    13
# 149761    Switzerland B.1.617+    44
# 149772 United Kingdom B.1.617+  4001
# 149774            USA B.1.617+   923

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>10,]$Country))
sel_countries_target
# [1] "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "Israel"        
# [10] "Italy"          "Japan"          "Netherlands"    "New Zealand"    "Poland"         "Singapore"      "Spain"          "Sweden"         "Switzerland"   
# [19] "United Kingdom" "USA" 

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "Israel"        
# [10] "Italy"          "Japan"          "Netherlands"    "New Zealand"    "Poland"         "Singapore"      "Spain"          "Sweden"         "Switzerland"   
# [19] "United Kingdom" "USA"    

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "Israel"        
# [10] "Italy"          "Japan"          "Netherlands"    "New Zealand"    "Poland"         "Singapore"      "Spain"          "Sweden"         "Switzerland"   
# [19] "United Kingdom" "USA"      

sel_ref_lineage = "B.1.1.7"
tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB117
#              Country Lineage  Count
# 70287      Australia B.1.1.7    367
# 70291        Bahrain B.1.1.7     12
# 70295        Belgium B.1.1.7  11283
# 70319        Denmark B.1.1.7  42504
# 70329         France B.1.1.7  19435
# 70334        Germany B.1.1.7  84395
# 70347          India B.1.1.7    708
# 70351        Ireland B.1.1.7   9939
# 70352         Israel B.1.1.7   8085
# 70353          Italy B.1.1.7  16029
# 70355          Japan B.1.1.7   2560
# 70382    Netherlands B.1.1.7  16897
# 70383    New Zealand B.1.1.7    134
# 70396         Poland B.1.1.7   9375
# 70409      Singapore B.1.1.7    170
# 70416          Spain B.1.1.7  11171
# 70419         Sweden B.1.1.7  27164
# 70420    Switzerland B.1.1.7  15574
# 70431 United Kingdom B.1.1.7 236514
# 70433            USA B.1.1.7 122229

data.frame(Country=tblB1617$Country, Lineage="B.1.617", Perc=100*tblB1617$Count / (tblB1617$Count+tblB117$Count))
#           Country Lineage        Perc
# 1       Australia B.1.617 20.56277056
# 2         Bahrain B.1.617 61.29032258
# 3         Belgium B.1.617  0.54649625
# 4         Denmark B.1.617  0.21129737
# 5          France B.1.617  0.25660765
# 6         Germany B.1.617  0.30712894
# 7           India B.1.617 76.74876847
# 8         Ireland B.1.617  0.89739755
# 9          Israel B.1.617  0.44329516
# 10          Italy B.1.617  0.35434539
# 11          Japan B.1.617  5.22028878
# 12    Netherlands B.1.617  0.18312854
# 13    New Zealand B.1.617 10.06711409
# 14         Poland B.1.617  0.27656632
# 15      Singapore B.1.617 47.85276074
# 16          Spain B.1.617  0.37456524
# 17         Sweden B.1.617  0.04783457
# 18    Switzerland B.1.617  0.28172621
# 19 United Kingdom B.1.617  1.66351371
# 20            USA B.1.617  0.74948032


GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_sel) # 1297772

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))




# ANALYSIS VOCs IN INDIA ####

GISAID_india = GISAID_sel[GISAID_sel$country=="India",]
nrow(GISAID_india[is.na(GISAID_india$LINEAGE1),]) # 0 unknown pango clade
GISAID_india = GISAID_india[!is.na(GISAID_india$LINEAGE1),]
nrow(GISAID_india) # 11317

unique(GISAID_india$division) # best data for West Bengal, Maharashtra & Karnataka (B.1.617 most common, B.1.1.7 most common in Punjab IN, Telangana)
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.617+"])
sum(GISAID_india$LINEAGE1=="B.1.617+") # 2337
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.1.7"])
sum(GISAID_india$LINEAGE1=="B.1.1.7") # 708

table(GISAID_india$LINEAGE1)
table(GISAID_india$LINEAGE2)

main_lineages = names(table(GISAID_india$LINEAGE1))[100*table(GISAID_india$LINEAGE1)/sum(table(GISAID_india$LINEAGE1)) > 3]
main_lineages
# "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.1.326" "B.1.1.7"   "B.1.36+"   "B.1.617+" 
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526")
main_lineages = union(main_lineages, VOCs)
GISAID_india$LINEAGE1[!(GISAID_india$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_india$LINEAGE2[!(GISAID_india$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove = names(table(GISAID_india$LINEAGE1))[table(GISAID_india$LINEAGE1) < 10]
GISAID_india$LINEAGE1[(GISAID_india$LINEAGE1 %in% remove)] = "other" # minority VOCs
GISAID_india$LINEAGE2[(GISAID_india$LINEAGE2 %in% remove)] = "other" # minority VOCs
table(GISAID_india$LINEAGE1)
GISAID_india$LINEAGE1 = factor(GISAID_india$LINEAGE1)
GISAID_india$LINEAGE1 = relevel(GISAID_india$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_india$LINEAGE1)
# [1] "B.1.1.7"   "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.1.326" "B.1.351"   "B.1.36+"   "B.1.525"   "B.1.617+"  "B.1.618"   "other" 
levels_LINEAGE1 = c("B.1.1.7","B.1","B.1.1","B.1.1.216","B.1.1.306","B.1.1.326","B.1.36+",
                    "B.1.525","B.1.351","B.1.618","B.1.617+","other")
GISAID_india$LINEAGE1 = factor(GISAID_india$LINEAGE1, levels=levels_LINEAGE1)

GISAID_india$LINEAGE2 = factor(GISAID_india$LINEAGE2)
GISAID_india$LINEAGE2 = relevel(GISAID_india$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_india$LINEAGE2)
# "B.1.1.7"   "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.1.326" "B.1.351"   "B.1.36+"   "B.1.525"   "B.1.617.1" "B.1.617.2" "B.1.618"   "other" 
levels_LINEAGE2 = c("B.1.1.7","B.1","B.1.1","B.1.1.216","B.1.1.306","B.1.1.326","B.1.36+",
                    "B.1.525","B.1.351","B.1.618","B.1.617.1","B.1.617.2","other")
GISAID_india$LINEAGE2 = factor(GISAID_india$LINEAGE2, levels=levels_LINEAGE2)
firstB16172 = GISAID_india[GISAID_india$LINEAGE2=="B.1.617.2",]
firstB16172 = firstB16172[firstB16172$date==min(firstB16172$date),]
firstB16172 # Maharashtra B.1.617.2

# select states of India with a total of > 300 sequences submitted

GISAID_india = GISAID_india[GISAID_india$division!="India",]
table(GISAID_india$division)

# states with at least 1 B.1.617.2 sequence
sel_states2 = names(table(GISAID_india[GISAID_india$pango_lineage=="B.1.617.2",]$division))[table(GISAID_india[GISAID_india$pango_lineage=="B.1.617.2",]$division) > 1]
sel_states2
# "Maharashtra"    "Chhattisgarh"   "Gujarat"        "Delhi"          "Karnataka"      "Andhra Pradesh" "West Bengal"    "Telangana"      "Odisha"  

# states with at least 300 seqs over the last year
sel_states = names(table(GISAID_india$division))[table(GISAID_india$division) > 300]
sel_states = sel_states[sel_states %in% sel_states2]
sel_states
# "Maharashtra"    "Chhattisgarh"   "Gujarat"        "Delhi"          "Karnataka"      "Andhra Pradesh" "West Bengal"    "Telangana"      "Odisha"  
GISAID_india = GISAID_india[GISAID_india$division %in% sel_states,]
GISAID_india$division = factor(GISAID_india$division)
GISAID_india$division = relevel(GISAID_india$division, ref="Maharashtra")
levels_STATES = c("Maharashtra","Chhattisgarh","Gujarat","Delhi","Karnataka","Andhra Pradesh","West Bengal","Telangana","Odisha")
GISAID_india$division = factor(GISAID_india$division, levels=levels_STATES)
# levels_STATES = levels(GISAID_india$division)
levels(GISAID_india$division)
# "Maharashtra"    "Chhattisgarh"   "Gujarat"        "Delhi"          "Karnataka"      "Andhra Pradesh" "West Bengal"    "Telangana"      "Odisha"
table(GISAID_india$division)
#    Maharashtra   Chhattisgarh        Gujarat          Delhi      Karnataka Andhra Pradesh    West Bengal      Telangana         Odisha 
#           3051            347           1413            422            806            640           1224           1741            309 
sum(table(GISAID_india$division)) # 9953
table(GISAID_india$LINEAGE1)

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of India
data_agbyweek1 = as.data.frame(table(GISAID_india$floor_date, GISAID_india$LINEAGE1))
colnames(data_agbyweek1) = c("floor_date", "LINEAGE1", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_india) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$LINEAGE1 = factor(data_agbyweek1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

data_agbyweek2 = as.data.frame(table(GISAID_india$floor_date, GISAID_india$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_india) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL


# aggregated by week and state for selected variant lineages
data_agbyweekregion1 = as.data.frame(table(GISAID_india$floor_date, GISAID_india$division, GISAID_india$LINEAGE1))
colnames(data_agbyweekregion1) = c("floor_date", "division", "LINEAGE1", "count")
data_agbyweekregion1_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$floor_date,data_agbyweekregion1$division), 
                                                                  interaction(data_agbyweekregion1_sum$floor_date,data_agbyweekregion1_sum$division))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_india) # correct
data_agbyweekregion1$collection_date = as.Date(as.character(data_agbyweekregion1$floor_date))
data_agbyweekregion1$LINEAGE1 = factor(data_agbyweekregion1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweekregion1$division = factor(data_agbyweekregion1$division, levels=levels_STATES)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]
data_agbyweekregion1$floor_date = NULL

data_agbyweekregion2 = as.data.frame(table(GISAID_india$floor_date, GISAID_india$division, GISAID_india$LINEAGE2))
colnames(data_agbyweekregion2) = c("floor_date", "division", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$division), 
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$division))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_india) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$division = factor(data_agbyweekregion2$division, levels=levels_STATES)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]
data_agbyweekregion2$floor_date = NULL


# MULLER PLOT (RAW DATA)
library(scales)
n1 = length(levels(GISAID_india$LINEAGE1))
lineage_cols1 = hcl(h = seq(15, 320, length = n1), l = 65, c = 200)
lineage_cols1[which(levels(GISAID_india$LINEAGE1)=="B.1.617+")] = "magenta"
lineage_cols1[which(levels(GISAID_india$LINEAGE1)=="other")] = "grey75"

n2 = length(levels(GISAID_india$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_india$LINEAGE2)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(GISAID_india$LINEAGE2)=="B.1.617.2")] = "magenta"
lineage_cols2[which(levels(GISAID_india$LINEAGE2)=="other")] = "grey75"

muller_india_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE1)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1, group=LINEAGE1), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANT B.1.617+ IN INDIA") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_india_raw1

muller_india_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\nRaw GISAID data") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_india_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_raw data.pdf"), width=8, height=6)


muller_indiabystate_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ division, ncol=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\nRaw GISAID data") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_indiabystate_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.png"), width=12, height=5)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.pdf"), width=12, height=5)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_india_multi = nnet::multinom(LINEAGE2 ~ division + DATE_NUM, data=GISAID_india, maxit=1000)
fit2_india_multi = nnet::multinom(LINEAGE2 ~ division * DATE_NUM, data=GISAID_india, maxit=1000)
fit3_india_multi = nnet::multinom(LINEAGE2 ~ division + ns(DATE_NUM, df=2), data=GISAID_india, maxit=1000)
fit4_india_multi = nnet::multinom(LINEAGE2 ~ division * ns(DATE_NUM, df=2), data=GISAID_india, maxit=1000)
fit5_india_multi = nnet::multinom(LINEAGE2 ~ division + ns(DATE_NUM, df=3), data=GISAID_india, maxit=1000)
fit6_india_multi = nnet::multinom(LINEAGE2 ~ division * ns(DATE_NUM, df=3), data=GISAID_india, maxit=1000)
BIC(fit1_india_multi, fit2_india_multi, fit3_india_multi, fit4_india_multi, fit5_india_multi, fit6_india_multi) 
# fit4_india_multi fits best (lowest BIC) but I will use model fit5_india_multi which has almost as good a BIC and which has identical growth advantages in different states

# equivalent fit with B.1.617.1,2&3 all recoded to B.1.617+
fit5_india_multi1 = nnet::multinom(LINEAGE1 ~ division + ns(DATE_NUM, df=3), data=GISAID_india, maxit=1000)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
max(GISAID_india$date) # 2021-05-07
emtrindia = emtrends(fit5_india_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india = data.frame(confint(emtrindia, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrindia$contrasts)$p.value)
delta_r_india
#               contrast     estimate           SE df    asymp.LCL    asymp.UCL p.value
# 1        B.1 - B.1.1.7  0.036775122 0.006601208 NA  0.023836991  0.049713253 1.434866e-06
# 2      B.1.1 - B.1.1.7 -0.001804869 0.007533115 NA -0.016569503  0.012959764 9.996986e-01
# 3  B.1.1.216 - B.1.1.7 -0.069686207 0.012682336 NA -0.094543129 -0.044829286 2.051313e-06
# 4  B.1.1.306 - B.1.1.7 -0.007082194 0.008575939 NA -0.023890726  0.009726337 9.473634e-01
# 5  B.1.1.326 - B.1.1.7  0.011968347 0.015598741 NA -0.018604623  0.042541317 9.602321e-01
# 6  (B.1.36+) - B.1.1.7 -0.062552254 0.007471739 NA -0.077196594 -0.047907914 5.266898e-13
# 7    B.1.525 - B.1.1.7 -0.094681896 0.021397274 NA -0.136619783 -0.052744010 2.180655e-04
# 8    B.1.351 - B.1.1.7 -0.013715182 0.014172249 NA -0.041492280  0.014061916 9.058252e-01
# 9    B.1.618 - B.1.1.7 -0.001884644 0.015838417 NA -0.032927371  0.029158082 9.999861e-01
# 10 B.1.617.1 - B.1.1.7  0.004413148 0.007043782 NA -0.009392410  0.018218707 9.821285e-01
# 11 B.1.617.2 - B.1.1.7  0.090596891 0.008505146 NA  0.073927111  0.107266670 0.000000e+00
# 12     other - B.1.1.7  0.022263986 0.006659817 NA  0.009210986  0.035316987 1.099427e-02

# avg growth advantage of B.1.617+ over B.1.1.7 :
max(GISAID_india$date) # 2021-05-07
emtrindia1 = emtrends(fit5_india_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india1 = data.frame(confint(emtrindia1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrindia1$contrasts)$p.value)
delta_r_india1
#                contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1         B.1 - B.1.1.7  0.0384291869 0.006287973 NA  0.02610499  0.05075339 1.136125e-07
# 2       B.1.1 - B.1.1.7  0.0032446368 0.007248955 NA -0.01096305  0.01745233 9.942285e-01
# 3   B.1.1.216 - B.1.1.7 -0.0632560878 0.012478837 NA -0.08771416 -0.03879802 1.428670e-05
# 4   B.1.1.306 - B.1.1.7 -0.0002688521 0.008280548 NA -0.01649843  0.01596072 9.999999e-01
# 5   B.1.1.326 - B.1.1.7  0.0160374449 0.015357456 NA -0.01406262  0.04613751 8.635479e-01
# 6   (B.1.36+) - B.1.1.7 -0.0569778585 0.007193407 NA -0.07107668 -0.04287904 9.467871e-12
# 7     B.1.525 - B.1.1.7 -0.0904114000 0.021100513 NA -0.13176765 -0.04905515 3.675195e-04
# 8     B.1.351 - B.1.1.7 -0.0121666430 0.013661396 NA -0.03894249  0.01460920 9.205374e-01
# 9     B.1.618 - B.1.1.7  0.0066017051 0.015737047 NA -0.02424234  0.03744575 9.955562e-01
# 10 (B.1.617+) - B.1.1.7  0.0535415162 0.006181704 NA  0.04142560  0.06565743 1.838529e-13
# 11      other - B.1.1.7  0.0270051236 0.006448026 NA  0.01436722  0.03964302 5.336922e-04


# pairwise growth rate difference (differences in growth rate per day) 
max(GISAID_india$date) # 2021-05-07
emtrindia_pairw = emtrends(fit5_india_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india_pairw = data.frame(confint(emtrindia_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindia_pairw$contrasts)$p.value)
delta_r_india_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1          B.1.1.7 - B.1 -3.677512e-02 0.006601208 NA -0.049713253 -0.023836991 9.128474e-06
# 2        B.1.1.7 - B.1.1  1.804869e-03 0.007533115 NA -0.012959764  0.016569503 1.000000e+00
# 3    B.1.1.7 - B.1.1.216  6.968621e-02 0.012682336 NA  0.044829286  0.094543129 1.302283e-05
# 4    B.1.1.7 - B.1.1.306  7.082194e-03 0.008575939 NA -0.009726337  0.023890726 9.997980e-01
# 5    B.1.1.7 - B.1.1.326 -1.196835e-02 0.015598741 NA -0.042541317  0.018604623 9.999067e-01
# 6    B.1.1.7 - (B.1.36+)  6.255225e-02 0.007471739 NA  0.047907914  0.077196594 3.598233e-12
# 7      B.1.1.7 - B.1.525  9.468190e-02 0.021397274 NA  0.052744010  0.136619783 1.302602e-03
# 8      B.1.1.7 - B.1.351  1.371518e-02 0.014172249 NA -0.014061916  0.041492280 9.989934e-01
# 9      B.1.1.7 - B.1.618  1.884644e-03 0.015838417 NA -0.029158082  0.032927371 1.000000e+00
# 10   B.1.1.7 - B.1.617.1 -4.413148e-03 0.007043782 NA -0.018218707  0.009392410 9.999896e-01
# 11   B.1.1.7 - B.1.617.2 -9.059689e-02 0.008505146 NA -0.107266670 -0.073927111 2.142730e-14
# 12       B.1.1.7 - other -2.226399e-02 0.006659817 NA -0.035316987 -0.009210986 5.394813e-02
# 13           B.1 - B.1.1  3.857999e-02 0.005836111 NA  0.027141424  0.050018558 5.374847e-08
# 14       B.1 - B.1.1.216  1.064613e-01 0.011638464 NA  0.083650358  0.129272300 8.870682e-14
# 15       B.1 - B.1.1.306  4.385732e-02 0.007034182 NA  0.030070574  0.057644059 3.641419e-07
# 16       B.1 - B.1.1.326  2.480677e-02 0.014816227 NA -0.004232496  0.053846045 9.016217e-01
# 17       B.1 - (B.1.36+)  9.932738e-02 0.005636296 NA  0.088280438  0.110374314 0.000000e+00
# 18         B.1 - B.1.525  1.314570e-01 0.021061684 NA  0.090176877  0.172737160 3.522216e-07
# 19         B.1 - B.1.351  5.049030e-02 0.013621286 NA  0.023793074  0.077187534 1.754637e-02
# 20         B.1 - B.1.618  3.865977e-02 0.015156639 NA  0.008953300  0.068366233 3.502548e-01
# 21       B.1 - B.1.617.1  3.236197e-02 0.005718161 NA  0.021154584  0.043569363 6.017662e-06
# 22       B.1 - B.1.617.2 -5.382177e-02 0.007412347 NA -0.068349702 -0.039293835 1.718181e-09
# 23           B.1 - other  1.451114e-02 0.004754662 NA  0.005192169  0.023830102 1.185066e-01
# 24     B.1.1 - B.1.1.216  6.788134e-02 0.011849817 NA  0.044656123  0.091106554 4.337210e-06
# 25     B.1.1 - B.1.1.306  5.277325e-03 0.007578310 NA -0.009575889  0.020130539 9.999670e-01
# 26     B.1.1 - B.1.1.326 -1.377322e-02 0.015026148 NA -0.043223925  0.015677493 9.994133e-01
# 27     B.1.1 - (B.1.36+)  6.074738e-02 0.006197945 NA  0.048599635  0.072895135 4.807266e-14
# 28       B.1.1 - B.1.525  9.287703e-02 0.021329728 NA  0.051071529  0.134682525 1.715674e-03
# 29       B.1.1 - B.1.351  1.191031e-02 0.014071336 NA -0.015668999  0.039489624 9.997393e-01
# 30       B.1.1 - B.1.618  7.977528e-05 0.015379987 NA -0.030064446  0.030223996 1.000000e+00
# 31     B.1.1 - B.1.617.1 -6.218017e-03 0.006842512 NA -0.019629094  0.007193059 9.994621e-01
# 32     B.1.1 - B.1.617.2 -9.240176e-02 0.008485899 NA -0.109033816 -0.075769703 3.996803e-15
# 33         B.1.1 - other -2.406886e-02 0.005677476 NA -0.035196505 -0.012941206 2.663761e-03
# 34 B.1.1.216 - B.1.1.306 -6.260401e-02 0.012590126 NA -0.087280207 -0.037927819 1.355936e-04
# 35 B.1.1.216 - B.1.1.326 -8.165455e-02 0.018124148 NA -0.117177232 -0.046131877 9.475967e-04
# 36 B.1.1.216 - (B.1.36+) -7.133953e-03 0.011634530 NA -0.029937214  0.015669307 9.999918e-01
# 37   B.1.1.216 - B.1.525  2.499569e-02 0.023774344 NA -0.021601169  0.071592546 9.977575e-01
# 38   B.1.1.216 - B.1.351 -5.597103e-02 0.017488048 NA -0.090246969 -0.021695082 8.036758e-02
# 39   B.1.1.216 - B.1.618 -6.780156e-02 0.018093214 NA -0.103263612 -0.032339514 1.534199e-02
# 40 B.1.1.216 - B.1.617.1 -7.409936e-02 0.012336630 NA -0.098278706 -0.049920005 1.130480e-06
# 41 B.1.1.216 - B.1.617.2 -1.602831e-01 0.013331971 NA -0.186413282 -0.134152914 0.000000e+00
# 42     B.1.1.216 - other -9.195019e-02 0.011555924 NA -0.114599389 -0.069300998 3.714296e-11
# 43 B.1.1.306 - B.1.1.326 -1.905054e-02 0.015293751 NA -0.049025743  0.010924660 9.895800e-01
# 44 B.1.1.306 - (B.1.36+)  5.547006e-02 0.007213857 NA  0.041331160  0.069608959 1.648727e-10
# 45   B.1.1.306 - B.1.525  8.759970e-02 0.021710892 NA  0.045047135  0.130152269 5.665992e-03
# 46   B.1.1.306 - B.1.351  6.632988e-03 0.014652576 NA -0.022085533  0.035351508 9.999997e-01
# 47   B.1.1.306 - B.1.618 -5.197550e-03 0.016027991 NA -0.036611835  0.026216735 1.000000e+00
# 48 B.1.1.306 - B.1.617.1 -1.149534e-02 0.007905306 NA -0.026989458  0.003998773 9.634819e-01
# 49 B.1.1.306 - B.1.617.2 -9.767908e-02 0.009380784 NA -0.116065084 -0.079293086 3.674838e-14
# 50     B.1.1.306 - other -2.934618e-02 0.006802753 NA -0.042679331 -0.016013030 2.005535e-03
# 51 B.1.1.326 - (B.1.36+)  7.452060e-02 0.014802275 NA  0.045508675  0.103532527 1.036120e-04
# 52   B.1.1.326 - B.1.525  1.066502e-01 0.025258185 NA  0.057145110  0.156155376 2.839639e-03
# 53   B.1.1.326 - B.1.351  2.568353e-02 0.019556644 NA -0.012646789  0.064013847 9.837209e-01
# 54   B.1.1.326 - B.1.618  1.385299e-02 0.020670765 NA -0.026660963  0.054366946 9.999783e-01
# 55 B.1.1.326 - B.1.617.1  7.555199e-03 0.015267490 NA -0.022368532  0.037478930 9.999993e-01
# 56 B.1.1.326 - B.1.617.2 -7.862854e-02 0.016059699 NA -0.110104974 -0.047152113 1.883415e-04
# 57     B.1.1.326 - other -1.029564e-02 0.014717449 NA -0.039141308  0.018550030 9.999653e-01
# 58   (B.1.36+) - B.1.525  3.212964e-02 0.021377193 NA -0.009768887  0.074028171 9.532852e-01
# 59   (B.1.36+) - B.1.351 -4.883707e-02 0.014126620 NA -0.076524739 -0.021149405 3.853577e-02
# 60   (B.1.36+) - B.1.618 -6.066761e-02 0.015347521 NA -0.090748198 -0.030587021 7.585574e-03
# 61 (B.1.36+) - B.1.617.1 -6.696540e-02 0.006859359 NA -0.080409499 -0.053521306 5.229150e-14
# 62 (B.1.36+) - B.1.617.2 -1.531491e-01 0.008539817 NA -0.169886878 -0.136411411 0.000000e+00
# 63     (B.1.36+) - other -8.481624e-02 0.005393420 NA -0.095387149 -0.074245332 0.000000e+00
# 64     B.1.525 - B.1.351 -8.096671e-02 0.024059623 NA -0.128122709 -0.033810719 5.058576e-02
# 65     B.1.525 - B.1.618 -9.279725e-02 0.025675895 NA -0.143121081 -0.042473423 2.367225e-02
# 66   B.1.525 - B.1.617.1 -9.909504e-02 0.021163036 NA -0.140573832 -0.057616257 4.617132e-04
# 67   B.1.525 - B.1.617.2 -1.852788e-01 0.021805718 NA -0.228017209 -0.142540364 1.792344e-12
# 68       B.1.525 - other -1.169459e-01 0.021059780 NA -0.158222292 -0.075669473 9.926807e-06
# 69     B.1.351 - B.1.618 -1.183054e-02 0.019884115 NA -0.050802687  0.027141612 9.999942e-01
# 70   B.1.351 - B.1.617.1 -1.812833e-02 0.013794920 NA -0.045165876  0.008909215 9.836344e-01
# 71   B.1.351 - B.1.617.2 -1.043121e-01 0.014590421 NA -0.132908771 -0.075715374 3.137234e-09
# 72       B.1.351 - other -3.597917e-02 0.013636069 NA -0.062705372 -0.009252965 2.976169e-01
# 73   B.1.618 - B.1.617.1 -6.297793e-03 0.015571287 NA -0.036816954  0.024221368 9.999999e-01
# 74   B.1.618 - B.1.617.2 -9.248153e-02 0.016402520 NA -0.124629884 -0.060333186 6.653363e-06
# 75       B.1.618 - other -2.414863e-02 0.015123362 NA -0.053789875  0.005492614 9.282381e-01
# 76 B.1.617.1 - B.1.617.2 -8.618374e-02 0.007247627 NA -0.100388829 -0.071978655 0.000000e+00
# 77     B.1.617.1 - other -1.785084e-02 0.005764178 NA -0.029148419 -0.006553258 1.057061e-01
# 78     B.1.617.2 - other  6.833290e-02 0.007576246 NA  0.053483735  0.083182073 1.317835e-13



# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-06-01"))
date.to = as.numeric(as.Date("2021-06-14")) # max(GISAID_india$DATE_NUM)+extrapolate

fit_india_multi_predsbystate = data.frame(emmeans(fit5_india_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "division"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=3)), # by=3 just to speed up things a bit
                                                  mode="prob", df=NA))
fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
fit_india_multi_predsbystate$LINEAGE2 = factor(fit_india_multi_predsbystate$LINEAGE2, levels=levels_LINEAGE2) 
fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$division, levels=levels_STATES) 

fit_india_multi_preds = data.frame(emmeans(fit5_india_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=3)), # by=3 just to speed up things a bit
                                           mode="prob", df=NA))
fit_india_multi_preds$collection_date = as.Date(fit_india_multi_preds$DATE_NUM, origin="1970-01-01")
fit_india_multi_preds$LINEAGE2 = factor(fit_india_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 

muller_india_mfit = ggplot(data=fit_india_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
           xmax=as.Date("2021-06-14"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\n(multinomial fit)")
muller_india_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_multinom fit.pdf"), width=8, height=6)


library(ggpubr)
ggarrange(muller_india_raw2+coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-14")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_india_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots multipanel_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots multipanel_multinom fit.pdf"), width=8, height=6)


muller_indiabystate_mfit = ggplot(data=fit_india_multi_predsbystate, 
                                  aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ STATE, ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
           xmax=as.Date("2021-06-14"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\n(multinomial fit)")
muller_indiabystate_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.png"), width=12, height=5)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.pdf"), width=12, height=5)

ggarrange(muller_indiabystate_raw2+
            coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-14")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0)))+
            ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\nRaw GISAID data"), 
          muller_indiabystate_mfit+ggtitle("\nMultinomial fit")+
            coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-05-31"))), nrow=2)

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state multipanel_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state multipanel_multinom fit.pdf"), width=8, height=6)

# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

fit_india_multi_preds2 = fit_india_multi_preds
fit_india_multi_preds2$LINEAGE2 = factor(fit_india_multi_preds2$LINEAGE2, levels=levels_LINEAGE1)
fit_india_multi_preds2$LINEAGE1 = fit_india_multi_preds2$LINEAGE2
levels(fit_india_multi_preds2$LINEAGE1)

# on logit scale:

fit_india_multi_preds2 = fit_india_multi_preds
ymin = 0.001
ymax = 0.998
fit_india_multi_preds2$asymp.LCL[fit_india_multi_preds2$asymp.LCL<ymin] = ymin
fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL<ymin] = ymin
fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL>ymax] = ymax
fit_india_multi_preds2$prob[fit_india_multi_preds2$prob<ymin] = ymin

plot_india_mfit_logit = qplot(data=fit_india_multi_preds2, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\n(multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
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
                        range=c(0.001, 5), limits=c(10,10^3), breaks=c(10,100)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-01")), ylim=c(0.005, 0.95), expand=c(0,0))
plot_india_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_multinom fit_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_multinom fit_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_india_mfit = qplot(data=fit_india_multi_preds, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share among new infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\n(multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01","2021-05-01")),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.001, 5), limits=c(1,10^3), breaks=c(10,100)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_india_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_multinom fit_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_multinom fit_response scale.pdf"), width=8, height=6)


# project multinomial fit onto case data ####
cases_india_bystate = read.csv("https://api.covid19india.org/csv/latest/states.csv") # cumulative cases
cases_india_bystate$Date = as.Date(cases_india_bystate$Date)
cases_india_bystate = cases_india_bystate[cases_india_bystate$Date >= as.Date("2020-06-01"),]
head(cases_india_bystate)
levels_STATES
# "Maharashtra"    "Chhattisgarh"   "Gujarat"        "West Bengal"    "Telangana"      "Andhra Pradesh" "Karnataka" 

cases_india_bystate = do.call(rbind,lapply(unique(cases_india_bystate$State), function (state) { df =  cases_india_bystate[cases_india_bystate$State==state,]
                                                             df$newcases = c(NA, diff(df$Confirmed))
                                                             return(df)
                                                             } ))
cases_india_bystate = cases_india_bystate[cases_india_bystate$State!="State Unassigned",]

# plot new cases per day by state
ggplot(data=cases_india_bystate,
       aes(x=Date, y=newcases, 
           group=State)) +
  facet_wrap(~ State, scale="free", ncol=5) +
  geom_smooth(aes(lwd=I(1), colour=State), method="loess", span=0.3, se=FALSE) +
  # geom_line(aes(lwd=I(1), colour=State)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-14",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY STATE IN INDIA") +
  scale_y_log10() +
  theme(legend.position = "none") # +
#  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_cases per day by state.png"), width=12, height=12)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_cases per day by state.pdf"), width=12, height=12)


cases_india_bystate2 = cases_india_bystate[cases_india_bystate$State %in% levels_STATES,]
colnames(cases_india_bystate2)[2]="STATE"

newdat = expand.grid(DATE_NUM=seq(as.numeric(min(cases_india_bystate2$Date)),as.numeric(max(cases_india_bystate2$Date))),
                     division=unique(as.character(cases_india_bystate2$STATE)))
fit_india_multi_predsbystate = data.frame(newdat,
                        predict(fit5_india_multi, 
                                  newdata = newdat,
                                  type = "prob"), check.names=F)  
fit_india_multi_predsbystate = gather(fit_india_multi_predsbystate, LINEAGE2, prob, all_of(levels_LINEAGE2))
fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
fit_india_multi_predsbystate$LINEAGE2 = factor(fit_india_multi_predsbystate$LINEAGE2, levels=levels_LINEAGE2)
colnames(fit_india_multi_predsbystate)[2] = "STATE"
fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$STATE, levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
                                                           "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))
fit_india_multi_predsbystate$totnewcases = cases_india_bystate2$newcases[match(interaction(fit_india_multi_predsbystate$STATE,fit_india_multi_predsbystate$collection_date),
                                                                               interaction(cases_india_bystate2$STATE,cases_india_bystate2$Date))]
fit_india_multi_predsbystate$cases = fit_india_multi_predsbystate$totnewcases*fit_india_multi_predsbystate$prob
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$cases==0] = NA
fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$STATE,
                                             levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
                                                      "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))

fit_india_multi_predsbystate2 = fit_india_multi_predsbystate
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases==0] = NA
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases<=1] = NA
fit_india_multi_predsbystate2$STATE = factor(fit_india_multi_predsbystate2$STATE,
                                             levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
                                                      "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))
cases_india_bystate2$STATE = factor(cases_india_bystate2$STATE,
                                    levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
                                             "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana"))
# sorted by date of introduction of B.1.617.2
ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases)) + 
  facet_wrap(~ STATE, scale="free", ncol=3) +
  geom_smooth(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2), method="loess", span=0.3, se=FALSE) +
  geom_smooth(data=cases_india_bystate2, aes(x=Date, y=newcases, lwd=I(1.5)), method="loess", span=0.3, se=FALSE, colour=alpha("black",0.6)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-05-31",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN INDIA\n(multinomial fit)") +
  scale_colour_manual("lineage", values=lineage_cols2) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases multinomial fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases multinomial fit.pdf"), width=8, height=6)

# TO DO: group together some strains in category other

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ STATE, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
           xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA\n(multinomial fit)")

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases stacked area multinomial fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases stacked area multinomial fit.pdf"), width=8, height=6)

