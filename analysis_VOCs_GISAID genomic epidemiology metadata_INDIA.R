# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN INDIA BASED ON ANALYSIS OF GISAID GENOMIC EPIDEMIOLOGY METADATA
# T. Wenseleers
# last update 1 JUNE 2021

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
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata (file version metadata_2021-05-31_10-15.tsv.gz)
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-05-31_10-15.tsv.gz"), col_types = cols(.default = "c")) 
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
range(GISAID$date) # "2019-12-24" "2021-05-26"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 10 dec Maharashtra India


# GISAID = GISAID[grepl("2021-", GISAID$date),]
GISAID = GISAID[GISAID$date>=as.Date("2020-06-01"),]
sum(is.na(GISAID$purpose_of_sequencing)) == nrow(GISAID) # field purpose_of_sequencing left blank unfortunately
range(GISAID$date) # "2020-06-01" "2021-05-26"
nrow(GISAID) # 1589229
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

length(unique(GISAID$country[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 60 countries
table(GISAID$pango_lineage[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])
# B.1.617 B.1.617.1 B.1.617.2 B.1.617.3 
# 1      2737     13613        79 

GISAID$pango_lineage[grepl("B.1.177",GISAID$pango_lineage,fixed=T)] = "B.1.177+"
GISAID$pango_lineage[grepl("B.1.36\\>",GISAID$pango_lineage)] = "B.1.36+"

sel_target_VOC = "B.1.617"
GISAID$LINEAGE1 = GISAID$pango_lineage
GISAID$LINEAGE2 = GISAID$pango_lineage
GISAID[grepl(sel_target_VOC, GISAID$LINEAGE1, fixed=TRUE),"LINEAGE1"] = paste0(sel_target_VOC,"+") # in LINEAGE1 we recode B.1.617.1,2&3 all as B.1.617+

table_country_lineage = as.data.frame(table(GISAID[GISAID$date>=as.Date("2021-04-01"),]$country, GISAID[GISAID$date>=as.Date("2021-04-01"),]$LINEAGE1))
colnames(table_country_lineage) = c("Country","Lineage","Count")
tblB1617 = table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage, fixed=T)&table_country_lineage$Count>10,]
tblB1617
# Country  Lineage Count
# 29903      Australia B.1.617+   140
# 29905        Bahrain B.1.617+    19
# 29906     Bangladesh B.1.617+    20
# 29907        Belgium B.1.617+   116
# 29922        Denmark B.1.617+    92
# 29929         France B.1.617+    85
# 29932        Germany B.1.617+   411
# 29938          India B.1.617+  1948
# 29939      Indonesia B.1.617+    24
# 29941        Ireland B.1.617+   154
# 29942         Israel B.1.617+    36
# 29943          Italy B.1.617+    81
# 29944          Japan B.1.617+   162
# 29954         Mexico B.1.617+    24
# 29957    Netherlands B.1.617+    50
# 29958    New Zealand B.1.617+    13
# 29962         Norway B.1.617+    20
# 29968         Poland B.1.617+    34
# 29969       Portugal B.1.617+    52
# 29974      Singapore B.1.617+   146
# 29978   South Africa B.1.617+    14
# 29980          Spain B.1.617+    65
# 29982         Sweden B.1.617+    14
# 29983    Switzerland B.1.617+    61
# 29991 United Kingdom B.1.617+  9191
# 29992            USA B.1.617+  1485

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>10,]$Country))
sel_countries_target
# [1] "Australia"      "Bahrain"        "Bangladesh"     "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Indonesia"     
# [10] "Ireland"        "Israel"         "Italy"          "Japan"          "Mexico"         "Netherlands"    "New Zealand"    "Norway"         "Poland"        
# [19] "Portugal"       "Singapore"      "South Africa"   "Spain"          "Sweden"         "Switzerland"    "United Kingdom" "USA"  

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Bangladesh"     "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Indonesia"      "Ireland"       
# [10] "Israel"         "Italy"          "Japan"          "Mexico"         "Netherlands"    "New Zealand"    "Norway"         "Poland"         "Portugal"      
# [19] "Singapore"      "South Africa"   "Spain"          "Sweden"         "Switzerland"    "United Kingdom" "USA"   

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Bangladesh"     "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Indonesia"      "Ireland"       
# [10] "Israel"         "Italy"          "Japan"          "Mexico"         "Netherlands"    "New Zealand"    "Norway"         "Poland"         "Portugal"      
# [19] "Singapore"      "South Africa"   "Spain"          "Sweden"         "Switzerland"    "United Kingdom" "USA"  

sel_ref_lineage = "B.1.1.7"
tblB1617 = table_country_lineage[grepl("B.1.617",table_country_lineage$Lineage)&table_country_lineage$Country %in% sel_countries,]
tblB1617

tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Country %in% sel_countries,]
tblB117
#              Country Lineage Count
# 13619      Australia B.1.1.7    135
# 13622     Bangladesh B.1.1.7     12
# 13623        Belgium B.1.1.7   5874
# 13638        Denmark B.1.1.7  25084
# 13645         France B.1.1.7   9455
# 13648        Germany B.1.1.7  51817
# 13654          India B.1.1.7    243
# 13655      Indonesia B.1.1.7     11
# 13657        Ireland B.1.1.7   3811
# 13658         Israel B.1.1.7   1433
# 13659          Italy B.1.1.7   7514
# 13660          Japan B.1.1.7   6660
# 13670         Mexico B.1.1.7    290
# 13673    Netherlands B.1.1.7   9357
# 13674    New Zealand B.1.1.7     42
# 13678         Norway B.1.1.7   2280
# 13684         Poland B.1.1.7   6111
# 13685       Portugal B.1.1.7   2183
# 13690      Singapore B.1.1.7     61
# 13694   South Africa B.1.1.7     12
# 13696          Spain B.1.1.7   4352
# 13698         Sweden B.1.1.7   9154
# 13699    Switzerland B.1.1.7   9542
# 13707 United Kingdom B.1.1.7  42437
# 13708            USA B.1.1.7 100743

data.frame(Country=tblB1617$Country, Lineage="B.1.617", Perc=100*tblB1617$Count / (tblB1617$Count+tblB117$Count))
# Country Lineage       Perc
# 1       Australia B.1.617 50.9090909
# 2      Bangladesh B.1.617 62.5000000
# 3         Belgium B.1.617  1.9365609
# 4         Denmark B.1.617  0.3654274
# 5          France B.1.617  0.8909853
# 6         Germany B.1.617  0.7869342
# 7           India B.1.617 88.9091739
# 8       Indonesia B.1.617 68.5714286
# 9         Ireland B.1.617  3.8839849
# 10         Israel B.1.617  2.4506467
# 11          Italy B.1.617  1.0664911
# 12          Japan B.1.617  2.3746702
# 13         Mexico B.1.617  7.6433121
# 14    Netherlands B.1.617  0.5315191
# 15    New Zealand B.1.617 23.6363636
# 16         Norway B.1.617  0.8695652
# 17         Poland B.1.617  0.5532954
# 18       Portugal B.1.617  2.3266219
# 19      Singapore B.1.617 70.5314010
# 20   South Africa B.1.617 53.8461538
# 21          Spain B.1.617  1.4715871
# 22         Sweden B.1.617  0.1527051
# 23    Switzerland B.1.617  0.6352182
# 24 United Kingdom B.1.617 17.8023553
# 25            USA B.1.617  1.4526353

GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_sel) # 1457634

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))




# ANALYSIS VOCs IN INDIA ####

GISAID_india = GISAID_sel[GISAID_sel$country=="India",]
nrow(GISAID_india[is.na(GISAID_india$LINEAGE1),]) # 0 unknown pango clade
GISAID_india = GISAID_india[!is.na(GISAID_india$LINEAGE1),]
nrow(GISAID_india) # 13650

unique(GISAID_india$division) # best data for West Bengal, Maharashtra & Karnataka (B.1.617 most common, B.1.1.7 most common in Punjab IN, Telangana)
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.617+"])
sum(GISAID_india$LINEAGE1=="B.1.617+") # 3577
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.1.7"])
sum(GISAID_india$LINEAGE1=="B.1.1.7") # 1032

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
# [1] "Andhra Pradesh"    "Bihar"             "Chhattisgarh"      "Delhi"             "Gujarat"           "Himachal Pradesh"  "Jammu and Kashmir"
# [8] "Jharkhand"         "Karnataka"         "Maharashtra"       "Odisha"            "Puducherry"        "Punjab IN"         "Sikkim"           
# [15] "Tamil Nadu"        "Telangana"         "West Bengal"    

# states with at least 300 seqs over the last year
sel_states = names(table(GISAID_india$division))[table(GISAID_india$division) > 300]
sel_states = sel_states[sel_states %in% sel_states2]
sel_states
# "Maharashtra"    "Chhattisgarh"   "Gujarat"        "Delhi"          "Karnataka"      "Andhra Pradesh" "West Bengal"    "Telangana"      "Odisha" "Jharkhand"  
GISAID_india = GISAID_india[GISAID_india$division %in% sel_states,]
GISAID_india$division = factor(GISAID_india$division)
GISAID_india$division = relevel(GISAID_india$division, ref="Maharashtra")
levels_STATES = c("Maharashtra","Chhattisgarh","Gujarat","Delhi","Andhra Pradesh","Telangana","Karnataka","West Bengal","Odisha","Jharkhand")
GISAID_india$division = factor(GISAID_india$division, levels=levels_STATES)
# levels_STATES = levels(GISAID_india$division)
levels(GISAID_india$division)
# "Maharashtra"    "Chhattisgarh"   "Gujarat"        "Delhi"          "Karnataka"      "Andhra Pradesh" "West Bengal"    "Telangana"      "Odisha"
table(GISAID_india$division)
# Maharashtra   Chhattisgarh        Gujarat          Delhi Andhra Pradesh      Telangana      Karnataka    West Bengal         Odisha      Jharkhand 
# 3733            550           1572            437           1005           1673            806           1224            594            378 
sum(table(GISAID_india$division)) # 11972
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
  facet_wrap(~ division, ncol=2) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.png"), width=5, height=7)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.pdf"), width=5, height=7)



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
# fit5_india_multi has best BIC

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
# 1        B.1 - B.1.1.7  0.063229201 0.006795533 NA  0.049910201  0.076548201 4.884981e-14
# 2      B.1.1 - B.1.1.7  0.054263232 0.008132227 NA  0.038324361  0.070202104 4.958875e-09
# 3  B.1.1.216 - B.1.1.7 -0.069605734 0.016508685 NA -0.101962162 -0.037249305 4.786372e-04
# 4  B.1.1.306 - B.1.1.7  0.002916035 0.010907500 NA -0.018462272  0.024294342 9.995147e-01
# 5  B.1.1.326 - B.1.1.7  0.052494605 0.022554098 NA  0.008289386  0.096699824 1.672380e-01
# 6  (B.1.36+) - B.1.1.7 -0.069338572 0.008814103 NA -0.086613896 -0.052063248 6.840750e-12
# 7    B.1.525 - B.1.1.7  0.036057035 0.021357585 NA -0.005803061  0.077917132 5.042862e-01
# 8    B.1.351 - B.1.1.7 -0.019205299 0.018484612 NA -0.055434473  0.017023874 8.792733e-01
# 9    B.1.618 - B.1.1.7 -0.045218572 0.024534766 NA -0.093305829  0.002868685 4.059583e-01
# 10 B.1.617.1 - B.1.1.7  0.009201712 0.007250477 NA -0.005008962  0.023412385 7.682902e-01
# 11 B.1.617.2 - B.1.1.7  0.112322340 0.007970182 NA  0.096701071  0.127943609 1.221245e-15
# 12     other - B.1.1.7  0.043613413 0.006910477 NA  0.030069127  0.057157699 3.273390e-08

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
# 1         B.1 - B.1.1.7  0.06188047 0.006544554 NA  0.049053379  0.07470756 3.830269e-14
# 2       B.1.1 - B.1.1.7  0.05508620 0.007873572 NA  0.039654283  0.07051812 1.026655e-09
# 3   B.1.1.216 - B.1.1.7 -0.06537531 0.016383637 NA -0.097486644 -0.03326397 1.084230e-03
# 4   B.1.1.306 - B.1.1.7  0.01020836 0.010709352 NA -0.010781589  0.03119830 8.995092e-01
# 5   B.1.1.326 - B.1.1.7  0.05417255 0.022077210 NA  0.010902011  0.09744308 1.190887e-01
# 6   (B.1.36+) - B.1.1.7 -0.06519412 0.008669192 NA -0.082185423 -0.04820281 6.075929e-11
# 7     B.1.525 - B.1.1.7  0.03734195 0.020880986 NA -0.003584028  0.07826793 4.211301e-01
# 8     B.1.351 - B.1.1.7 -0.01414275 0.018129538 NA -0.049675989  0.02139049 9.505942e-01
# 9     B.1.618 - B.1.1.7 -0.03789139 0.024473297 NA -0.085858173  0.01007539 5.756644e-01
# 10 (B.1.617+) - B.1.1.7  0.07306817 0.006650680 NA  0.060033078  0.08610327 1.754152e-14
# 11      other - B.1.1.7  0.04774345 0.006754216 NA  0.034505431  0.06098147 6.984704e-10


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
# 1          B.1.1.7 - B.1 -0.0632292013 0.006795533 NA -0.076548201 -0.049910201 1.314504e-13
# 2        B.1.1.7 - B.1.1 -0.0542632323 0.008132227 NA -0.070202104 -0.038324361 3.207649e-08
# 3    B.1.1.7 - B.1.1.216  0.0696057335 0.016508685 NA  0.037249305  0.101962162 2.802719e-03
# 4    B.1.1.7 - B.1.1.306 -0.0029160352 0.010907500 NA -0.024294342  0.018462272 1.000000e+00
# 5    B.1.1.7 - B.1.1.326 -0.0524946048 0.022554098 NA -0.096699824 -0.008289386 5.008241e-01
# 6    B.1.1.7 - (B.1.36+)  0.0693385718 0.008814103 NA  0.052063248  0.086613896 4.423306e-11
# 7      B.1.1.7 - B.1.525 -0.0360570352 0.021357585 NA -0.077917132  0.005803061 8.965055e-01
# 8      B.1.1.7 - B.1.351  0.0192052992 0.018484612 NA -0.017023874  0.055434473 9.980139e-01
# 9      B.1.1.7 - B.1.618  0.0452185720 0.024534766 NA -0.002868685  0.093305829 8.246837e-01
# 10   B.1.1.7 - B.1.617.1 -0.0092017118 0.007250477 NA -0.023412385  0.005008962 9.878563e-01
# 11   B.1.1.7 - B.1.617.2 -0.1123223400 0.007970182 NA -0.127943609 -0.096701071 1.443290e-15
# 12       B.1.1.7 - other -0.0436134133 0.006910477 NA -0.057157699 -0.030069127 2.111267e-07
# 13           B.1 - B.1.1  0.0089659690 0.006395497 NA -0.003568974  0.021500912 9.725888e-01
# 14       B.1 - B.1.1.216  0.1328349349 0.015702093 NA  0.102059398  0.163610472 1.545541e-12
# 15       B.1 - B.1.1.306  0.0603131662 0.009619443 NA  0.041459404  0.079166929 2.608418e-07
# 16       B.1 - B.1.1.326  0.0107345965 0.021971318 NA -0.032328396  0.053797589 9.999994e-01
# 17       B.1 - (B.1.36+)  0.1325677731 0.007274891 NA  0.118309249  0.146826298 1.332268e-15
# 18         B.1 - B.1.525  0.0271721662 0.020904620 NA -0.013800135  0.068144468 9.851476e-01
# 19         B.1 - B.1.351  0.0824345006 0.017999546 NA  0.047156039  0.117712962 6.678519e-04
# 20         B.1 - B.1.618  0.1084477733 0.024054392 NA  0.061302031  0.155593516 8.925536e-04
# 21       B.1 - B.1.617.1  0.0540274895 0.005457174 NA  0.043331625  0.064723354 1.217915e-13
# 22       B.1 - B.1.617.2 -0.0490931387 0.006104432 NA -0.061057605 -0.037128672 1.623157e-11
# 23           B.1 - other  0.0196157880 0.004812901 NA  0.010182676  0.029048900 4.737454e-03
# 24     B.1.1 - B.1.1.216  0.1238689658 0.016153812 NA  0.092208077  0.155529855 1.366798e-10
# 25     B.1.1 - B.1.1.306  0.0513471971 0.010412949 NA  0.030938193  0.071756201 1.515831e-04
# 26     B.1.1 - B.1.1.326  0.0017686275 0.022310811 NA -0.041959759  0.045497014 1.000000e+00
# 27     B.1.1 - (B.1.36+)  0.1236018041 0.008259413 NA  0.107413653  0.139789956 1.332268e-15
# 28       B.1.1 - B.1.525  0.0182061971 0.021350829 NA -0.023640658  0.060053052 9.997216e-01
# 29       B.1.1 - B.1.351  0.0734685315 0.018510473 NA  0.037188672  0.109748391 6.969541e-03
# 30       B.1.1 - B.1.618  0.0994818043 0.024398829 NA  0.051660979  0.147302630 4.708848e-03
# 31     B.1.1 - B.1.617.1  0.0450615205 0.007134900 NA  0.031077373  0.059045668 2.063622e-07
# 32     B.1.1 - B.1.617.2 -0.0580591077 0.007734750 NA -0.073218939 -0.042899277 3.398379e-10
# 33         B.1.1 - other  0.0106498189 0.006402961 NA -0.001899753  0.023199391 9.060118e-01
# 34 B.1.1.216 - B.1.1.306 -0.0725217687 0.017624604 NA -0.107065358 -0.037978179 4.100930e-03
# 35 B.1.1.216 - B.1.1.326 -0.1221003383 0.026503572 NA -0.174046386 -0.070154291 5.975234e-04
# 36 B.1.1.216 - (B.1.36+) -0.0002671618 0.016134136 NA -0.031889487  0.031355164 1.000000e+00
# 37   B.1.1.216 - B.1.525 -0.1056627687 0.025947384 NA -0.156518706 -0.054806831 4.798220e-03
# 38   B.1.1.216 - B.1.351 -0.0504004343 0.023651578 NA -0.096756676 -0.004044193 6.420777e-01
# 39   B.1.1.216 - B.1.618 -0.0243871616 0.027839553 NA -0.078951684  0.030177361 9.996330e-01
# 40 B.1.1.216 - B.1.617.1 -0.0788074453 0.016181090 NA -0.110521800 -0.047093091 1.971441e-04
# 41 B.1.1.216 - B.1.617.2 -0.1819280736 0.016528524 NA -0.214323386 -0.149532761 7.027712e-14
# 42     B.1.1.216 - other -0.1132191469 0.015618764 NA -0.143831362 -0.082606932 1.421494e-09
# 43 B.1.1.306 - B.1.1.326 -0.0495785696 0.023082567 NA -0.094819569 -0.004337570 6.300762e-01
# 44 B.1.1.306 - (B.1.36+)  0.0722546069 0.010647097 NA  0.051386680  0.093122534 1.752528e-08
# 45   B.1.1.306 - B.1.525 -0.0331410000 0.022596706 NA -0.077429731  0.011147731 9.612208e-01
# 46   B.1.1.306 - B.1.351  0.0221213344 0.019909221 NA -0.016900022  0.061142691 9.962696e-01
# 47   B.1.1.306 - B.1.618  0.0481346071 0.025428898 NA -0.001705117  0.097974332 7.970312e-01
# 48 B.1.1.306 - B.1.617.1 -0.0062856766 0.010209625 NA -0.026296175  0.013724822 9.999916e-01
# 49 B.1.1.306 - B.1.617.2 -0.1094063049 0.010777533 NA -0.130529882 -0.088282728 1.153522e-13
# 50     B.1.1.306 - other -0.0406973782 0.009529943 NA -0.059375724 -0.022019032 2.279168e-03
# 51 B.1.1.326 - (B.1.36+)  0.1218331766 0.022377748 NA  0.077973597  0.165692756 1.493016e-05
# 52   B.1.1.326 - B.1.525  0.0164375696 0.029865962 NA -0.042098641  0.074973780 9.999976e-01
# 53   B.1.1.326 - B.1.351  0.0716999040 0.027895586 NA  0.017025560  0.126374248 3.373871e-01
# 54   B.1.1.326 - B.1.618  0.0977131767 0.032229926 NA  0.034543683  0.160882670 1.237826e-01
# 55 B.1.1.326 - B.1.617.1  0.0432928930 0.022232010 NA -0.000281046  0.086866832 7.646327e-01
# 56 B.1.1.326 - B.1.617.2 -0.0598277352 0.022467667 NA -0.103863553 -0.015791917 2.830882e-01
# 57     B.1.1.326 - other  0.0088811914 0.021953829 NA -0.034147523  0.051909906 9.999999e-01
# 58   (B.1.36+) - B.1.525 -0.1053956069 0.021796875 NA -0.148116697 -0.062674517 2.290913e-04
# 59   (B.1.36+) - B.1.351 -0.0501332726 0.019000622 NA -0.087373807 -0.012892738 2.968391e-01
# 60   (B.1.36+) - B.1.618 -0.0241199998 0.024488996 NA -0.072117550  0.023877551 9.988149e-01
# 61 (B.1.36+) - B.1.617.1 -0.0785402836 0.008171860 NA -0.094556835 -0.062523732 1.281197e-13
# 62 (B.1.36+) - B.1.617.2 -0.1816609118 0.008880237 NA -0.199065856 -0.164255967 1.332268e-15
# 63     (B.1.36+) - other -0.1129519851 0.007104970 NA -0.126877471 -0.099026500 1.332268e-15
# 64     B.1.525 - B.1.351  0.0552623344 0.026401015 NA  0.003517296  0.107007372 6.685467e-01
# 65     B.1.525 - B.1.618  0.0812756071 0.031945444 NA  0.018663687  0.143887527 3.536363e-01
# 66   B.1.525 - B.1.617.1  0.0268553234 0.020987312 NA -0.014279052  0.067989698 9.869800e-01
# 67   B.1.525 - B.1.617.2 -0.0762653049 0.021217040 NA -0.117849939 -0.034680671 2.475498e-02
# 68       B.1.525 - other -0.0075563782 0.020976558 NA -0.048669676  0.033556920 1.000000e+00
# 69     B.1.351 - B.1.618  0.0260132727 0.029860968 NA -0.032513148  0.084539694 9.996533e-01
# 70   B.1.351 - B.1.617.1 -0.0284070110 0.018093008 NA -0.063868656  0.007054633 9.363740e-01
# 71   B.1.351 - B.1.617.2 -0.1315276393 0.018442071 NA -0.167673433 -0.095381845 2.703544e-09
# 72       B.1.351 - other -0.0628187126 0.018062656 NA -0.098220867 -0.027416558 3.565174e-02
# 73   B.1.618 - B.1.617.1 -0.0544202837 0.024327896 NA -0.102102084 -0.006738483 5.660634e-01
# 74   B.1.618 - B.1.617.2 -0.1575409120 0.024582146 NA -0.205721032 -0.109360792 1.276890e-07
# 75       B.1.618 - other -0.0888319853 0.024009127 NA -0.135889009 -0.041774962 1.758160e-02
# 76 B.1.617.1 - B.1.617.2 -0.1031206282 0.006445088 NA -0.115752769 -0.090488487 1.332268e-15
# 77     B.1.617.1 - other -0.0344117016 0.005684718 NA -0.045553543 -0.023269860 7.804298e-07
# 78     B.1.617.2 - other  0.0687089267 0.006503390 NA  0.055962516  0.081455337 1.015854e-13



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
  facet_wrap(~ STATE, ncol=2) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.png"), width=6, height=8)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.pdf"), width=6, height=8)

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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state multipanel_multinom fit.png"), width=7, height=11)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state multipanel_multinom fit.pdf"), width=7, height=11)

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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-14")), ylim=c(0.005, 0.95), expand=c(0,0))
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
  coord_cartesian(xlim=as.Date(c("2021-01-01","2021-06-14")),
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
# [1] "Maharashtra"    "Chhattisgarh"   "Gujarat"        "Delhi"          "Andhra Pradesh" "Telangana"      "Karnataka"      "West Bengal"    "Odisha"        
# [10] "Jharkhand"


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
                                                           "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana", "Jharkhand"))
fit_india_multi_predsbystate$totnewcases = cases_india_bystate2$newcases[match(interaction(fit_india_multi_predsbystate$STATE,fit_india_multi_predsbystate$collection_date),
                                                                               interaction(cases_india_bystate2$STATE,cases_india_bystate2$Date))]
fit_india_multi_predsbystate$cases = fit_india_multi_predsbystate$totnewcases*fit_india_multi_predsbystate$prob
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$cases==0] = NA
fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$STATE,
                                             levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
                                                      "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana", "Jharkhand"))

fit_india_multi_predsbystate2 = fit_india_multi_predsbystate
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases==0] = NA
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases<=1] = NA
fit_india_multi_predsbystate2$STATE = factor(fit_india_multi_predsbystate2$STATE,
                                             levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
                                                      "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana", "Jharkhand"))
cases_india_bystate2$STATE = factor(cases_india_bystate2$STATE,
                                    levels=c("Maharashtra","Chhattisgarh","Gujarat","Delhi",
                                             "Karnataka", "West Bengal", "Odisha", "Andhra Pradesh", "Telangana", "Jharkhand"))
# sorted by date of introduction of B.1.617.2
ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases)) + 
  facet_wrap(~ STATE, scale="free", ncol=2) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases multinomial fit.png"), width=8, height=10)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases multinomial fit.pdf"), width=8, height=10)

# TO DO: group together some strains in category other

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ STATE, scale="free", ncol=2) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases stacked area multinomial fit.png"), width=8, height=10)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases stacked area multinomial fit.pdf"), width=8, height=10)

