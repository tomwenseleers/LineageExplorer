# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN INDIA BASED ON ANALYSIS OF GISAID GENOMIC EPIDEMIOLOGY METADATA
# T. Wenseleers
# last update 8 JUNE 2021

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
today = as.Date("2021-06-08")
today_num = as.numeric(today)
today # "2021-06-08"
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata (file version metadata_2021-06-04_12-59.tsv.gz)
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-06-04_12-59.tsv.gz"), col_types = cols(.default = "c")) 
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
range(GISAID$date) # "2019-12-24" "2021-06-02"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 21 nov Varanasi, Uttar Pradesh


# GISAID = GISAID[grepl("2021-", GISAID$date),]
GISAID = GISAID[GISAID$date>=as.Date("2020-06-01"),]
sum(is.na(GISAID$purpose_of_sequencing)) == nrow(GISAID) # field purpose_of_sequencing left blank unfortunately
range(GISAID$date) # "2020-06-01" "2021-06-02"
nrow(GISAID) # 1640010
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
# 1      2922     21684        84

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
nrow(GISAID_sel) # 1506684

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))




# ANALYSIS VOCs IN INDIA ####

GISAID_india = GISAID_sel[GISAID_sel$country=="India",]
nrow(GISAID_india[is.na(GISAID_india$LINEAGE1),]) # 0 unknown pango clade
GISAID_india = GISAID_india[!is.na(GISAID_india$LINEAGE1),]
nrow(GISAID_india) # 14814

unique(GISAID_india$division) # best data for West Bengal, Maharashtra & Karnataka (B.1.617 most common, B.1.1.7 most common in Punjab IN, Telangana)
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.617+"])
sum(GISAID_india$LINEAGE1=="B.1.617+") # 4321
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.1.7"])
sum(GISAID_india$LINEAGE1=="B.1.1.7") # 1077

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
# "B.1.1.7"   "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.351"   "B.1.36+"   "B.1.525"   "B.1.617+"  "B.1.618"   "other"  
levels_LINEAGE1 = c("B.1.1.7","B.1","B.1.1","B.1.1.216","B.1.1.306","B.1.36+",
                    "B.1.525","B.1.351","B.1.618","B.1.617+","other")
GISAID_india$LINEAGE1 = factor(GISAID_india$LINEAGE1, levels=levels_LINEAGE1)

GISAID_india$LINEAGE2 = factor(GISAID_india$LINEAGE2)
GISAID_india$LINEAGE2 = relevel(GISAID_india$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_india$LINEAGE2)
# "B.1.1.7"   "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.351"   "B.1.36+"   "B.1.525"   "B.1.617.1" "B.1.617.2" "B.1.618"   "other"   
levels_LINEAGE2 = c("B.1.1.7","B.1","B.1.1","B.1.1.216","B.1.1.306","B.1.36+",
                    "B.1.525","B.1.351","B.1.618","B.1.617.1","B.1.617.2","other")
GISAID_india$LINEAGE2 = factor(GISAID_india$LINEAGE2, levels=levels_LINEAGE2)
firstB16172 = GISAID_india[GISAID_india$LINEAGE2=="B.1.617.2",]
firstB16172 = firstB16172[firstB16172$date==min(firstB16172$date),]
firstB16172 # 21 nov Varanasi, Uttar Pradesh

# select states of India with a total of > 300 sequences submitted

GISAID_india = GISAID_india[GISAID_india$division!="India",]
table(GISAID_india$division)

# nrow(GISAID_india[GISAID_india$division=="Delhi"&grepl("Dhar",GISAID_india$authors),]) # 59
# nrow(GISAID_india[grepl("Dhar",GISAID_india$authors),]) # 59


# states with at least 1 B.1.617.2 sequence
sel_states2 = names(table(GISAID_india[GISAID_india$pango_lineage=="B.1.617.2",]$division))[table(GISAID_india[GISAID_india$pango_lineage=="B.1.617.2",]$division) > 1]
sel_states2
# [1] "Andhra Pradesh"    "Bihar"             "Chhattisgarh"      "Delhi"             "Gujarat"           "Himachal Pradesh"  "Jammu and Kashmir"
# [8] "Jharkhand"         "Karnataka"         "Maharashtra"       "Odisha"            "Puducherry"        "Punjab IN"         "Sikkim"           
# [15] "Tamil Nadu"        "Telangana"         "West Bengal"    

# states with at least 350 seqs over the last year
sel_states = names(table(GISAID_india$division))[table(GISAID_india$division) > 350]
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
# 3797            550           1619            437           1219           1898            922           1315            594            378 
sum(table(GISAID_india$division)) # 12729
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
# fit4_india_multi has best BIC, but I'll use the slightly simpler fit3_india_multi

# equivalent fit with B.1.617.1,2&3 all recoded to B.1.617+
fit3_india_multi1 = nnet::multinom(LINEAGE1 ~ division + ns(DATE_NUM, df=2), data=GISAID_india, maxit=1000)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
max(GISAID_india$date) # 2021-05-16
emtrindia = emtrends(fit3_india_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india = data.frame(confint(emtrindia, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrindia$contrasts)$p.value)
delta_r_india
#               contrast     estimate           SE df    asymp.LCL    asymp.UCL p.value
# 1        B.1 - B.1.1.7  0.042277631 0.003807543 NA  0.034814984  0.049740278 2.775558e-15
# 2      B.1.1 - B.1.1.7  0.011051415 0.004173099 NA  0.002872292  0.019230539 7.499933e-02
# 3  B.1.1.216 - B.1.1.7 -0.033484155 0.005843331 NA -0.044936874 -0.022031436 7.067854e-07
# 4  B.1.1.306 - B.1.1.7 -0.006155248 0.004471662 NA -0.014919545  0.002609048 6.859396e-01
# 5  (B.1.36+) - B.1.1.7 -0.017178358 0.004048167 NA -0.025112619 -0.009244097 4.313630e-04
# 6    B.1.525 - B.1.1.7 -0.051636247 0.015142153 NA -0.081314321 -0.021958172 8.359990e-03
# 7    B.1.351 - B.1.1.7 -0.018432704 0.010154006 NA -0.038334191  0.001468782 4.049730e-01
# 8    B.1.618 - B.1.1.7  0.012732586 0.008399823 NA -0.003730765  0.029195938 5.970366e-01
# 9  B.1.617.1 - B.1.1.7  0.005721352 0.004788643 NA -0.003664215  0.015106919 7.911520e-01
# 10 B.1.617.2 - B.1.1.7  0.089176885 0.005597130 NA  0.078206712  0.100147057 0.000000e+00
# 11     other - B.1.1.7  0.023236473 0.003834935 NA  0.015720138  0.030752808 1.465993e-07

# avg growth advantage of B.1.617+ over B.1.1.7 :
max(GISAID_india$date) # 2021-05-16
emtrindia1 = emtrends(fit3_india_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india1 = data.frame(confint(emtrindia1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrindia1$contrasts)$p.value)
delta_r_india1
#                contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1         B.1 - B.1.1.7  0.042653240 0.003750514 NA  0.0353023676  0.050004112 3.808065e-14
# 2       B.1.1 - B.1.1.7  0.012388836 0.004129562 NA  0.0042950446  0.020482628 2.726749e-02
# 3   B.1.1.216 - B.1.1.7 -0.031914320 0.005807350 NA -0.0432965159 -0.020532123 2.200478e-06
# 4   B.1.1.306 - B.1.1.7 -0.003587874 0.004445151 NA -0.0123002103  0.005124462 9.347585e-01
# 5   (B.1.36+) - B.1.1.7 -0.015745203 0.004004480 NA -0.0235938396 -0.007896566 1.323542e-03
# 6     B.1.525 - B.1.1.7 -0.050028978 0.015105907 NA -0.0796360116 -0.020421944 1.072388e-02
# 7     B.1.351 - B.1.1.7 -0.016160646 0.009937404 NA -0.0356375993  0.003316307 5.029251e-01
# 8     B.1.618 - B.1.1.7  0.015909515 0.008369239 NA -0.0004938925  0.032312922 3.361977e-01
# 9  (B.1.617+) - B.1.1.7  0.057198880 0.004042808 NA  0.0492751223  0.065122638 2.164935e-14
# 10      other - B.1.1.7  0.023902906 0.003800319 NA  0.0164544180  0.031351393 5.340649e-08


# pairwise growth rate difference (differences in growth rate per day) 
max(GISAID_india$date) # 2021-05-16
emtrindia_pairw = emtrends(fit3_india_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india_pairw = data.frame(confint(emtrindia_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindia_pairw$contrasts)$p.value)
delta_r_india_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1          B.1.1.7 - B.1 -0.042277631 0.003807543 NA -0.049740278 -0.0348149843 3.774758e-14
# 2        B.1.1.7 - B.1.1 -0.011051415 0.004173099 NA -0.019230539 -0.0028722920 2.644421e-01
# 3    B.1.1.7 - B.1.1.216  0.033484155 0.005843331 NA  0.022031436  0.0449368737 4.167164e-06
# 4    B.1.1.7 - B.1.1.306  0.006155248 0.004471662 NA -0.002609048  0.0149195445 9.660465e-01
# 5    B.1.1.7 - (B.1.36+)  0.017178358 0.004048167 NA  0.009244097  0.0251126187 2.345443e-03
# 6      B.1.1.7 - B.1.525  0.051636247 0.015142153 NA  0.021958172  0.0813143213 3.928512e-02
# 7      B.1.1.7 - B.1.351  0.018432704 0.010154006 NA -0.001468782  0.0383341907 8.067440e-01
# 8      B.1.1.7 - B.1.618 -0.012732586 0.008399823 NA -0.029195938  0.0037307649 9.336814e-01
# 9    B.1.1.7 - B.1.617.1 -0.005721352 0.004788643 NA -0.015106919  0.0036642150 9.885528e-01
# 10   B.1.1.7 - B.1.617.2 -0.089176885 0.005597130 NA -0.100147057 -0.0782067122 0.000000e+00
# 11       B.1.1.7 - other -0.023236473 0.003834935 NA -0.030752808 -0.0157201379 8.692612e-07
# 12           B.1 - B.1.1  0.031226216 0.002390849 NA  0.026540238  0.0359121943 0.000000e+00
# 13       B.1 - B.1.1.216  0.075761786 0.004695025 NA  0.066559706  0.0849638654 0.000000e+00
# 14       B.1 - B.1.1.306  0.048432880 0.002808971 NA  0.042927397  0.0539383625 0.000000e+00
# 15       B.1 - (B.1.36+)  0.059455989 0.002056396 NA  0.055425527  0.0634864514 0.000000e+00
# 16         B.1 - B.1.525  0.093913878 0.014897040 NA  0.064716216  0.1231115401 2.620343e-07
# 17         B.1 - B.1.351  0.060710335 0.009784962 NA  0.041532163  0.0798885077 4.282276e-07
# 18         B.1 - B.1.618  0.029545045 0.007753979 NA  0.014347525  0.0447425651 1.101173e-02
# 19       B.1 - B.1.617.1  0.036556279 0.003719128 NA  0.029266922  0.0438456363 8.826273e-14
# 20       B.1 - B.1.617.2 -0.046899253 0.004596381 NA -0.055907996 -0.0378905111 8.226753e-14
# 21           B.1 - other  0.019041158 0.001762884 NA  0.015585969  0.0224963475 5.984102e-14
# 22     B.1.1 - B.1.1.216  0.044535570 0.004845682 NA  0.035038207  0.0540329328 1.333378e-13
# 23     B.1.1 - B.1.1.306  0.017206664 0.003065968 NA  0.011197477  0.0232158501 7.227233e-06
# 24     B.1.1 - (B.1.36+)  0.028229773 0.002418563 NA  0.023489476  0.0329700701 2.331468e-15
# 25       B.1.1 - B.1.525  0.062687662 0.015002349 NA  0.033283599  0.0920917250 2.988017e-03
# 26       B.1.1 - B.1.351  0.029484119 0.009943194 NA  0.009995817  0.0489724222 1.310502e-01
# 27       B.1.1 - B.1.618 -0.001681171 0.007901345 NA -0.017167523  0.0138051805 1.000000e+00
# 28     B.1.1 - B.1.617.1  0.005330063 0.004167623 NA -0.002838327  0.0134984537 9.804170e-01
# 29     B.1.1 - B.1.617.2 -0.078125469 0.005015651 NA -0.087955964 -0.0682949747 0.000000e+00
# 30         B.1.1 - other -0.012185057 0.002205094 NA -0.016506962 -0.0078631530 1.075547e-05
# 31 B.1.1.216 - B.1.1.306 -0.027328906 0.005070084 NA -0.037266088 -0.0173917241 1.993595e-05
# 32 B.1.1.216 - (B.1.36+) -0.016305797 0.004645145 NA -0.025410113 -0.0072014810 2.901911e-02
# 33   B.1.1.216 - B.1.525  0.018152092 0.015578533 NA -0.012381272  0.0486854565 9.906643e-01
# 34   B.1.1.216 - B.1.351 -0.015051451 0.010765081 NA -0.036150621  0.0060477202 9.619979e-01
# 35   B.1.1.216 - B.1.618 -0.046216741 0.008755952 NA -0.063378092 -0.0290553901 3.292123e-05
# 36 B.1.1.216 - B.1.617.1 -0.039205507 0.005874883 NA -0.050720066 -0.0276909466 4.110579e-08
# 37 B.1.1.216 - B.1.617.2 -0.122661039 0.006492062 NA -0.135385246 -0.1099368322 0.000000e+00
# 38     B.1.1.216 - other -0.056720627 0.004605914 NA -0.065748053 -0.0476932018 0.000000e+00
# 39 B.1.1.306 - (B.1.36+)  0.011023109 0.002743064 NA  0.005646802  0.0163994167 5.343771e-03
# 40   B.1.1.306 - B.1.525  0.045480998 0.015085428 NA  0.015914102  0.0750478948 1.159609e-01
# 41   B.1.1.306 - B.1.351  0.012277456 0.010070838 NA -0.007461024  0.0320159353 9.865443e-01
# 42   B.1.1.306 - B.1.618 -0.018887835 0.008064496 NA -0.034693957 -0.0030817130 4.527322e-01
# 43 B.1.1.306 - B.1.617.1 -0.011876600 0.004469184 NA -0.020636040 -0.0031171605 2.596357e-01
# 44 B.1.1.306 - B.1.617.2 -0.095332133 0.005262016 NA -0.105645495 -0.0850187713 0.000000e+00
# 45     B.1.1.306 - other -0.029391721 0.002567881 NA -0.034424675 -0.0243587678 1.387779e-14
# 46   (B.1.36+) - B.1.525  0.034457889 0.014973289 NA  0.005110781  0.0638049967 4.810140e-01
# 47   (B.1.36+) - B.1.351  0.001254346 0.009897608 NA -0.018144609  0.0206533015 1.000000e+00
# 48   (B.1.36+) - B.1.618 -0.029910944 0.007823843 NA -0.045245394 -0.0145764939 1.054715e-02
# 49 (B.1.36+) - B.1.617.1 -0.022899710 0.004049926 NA -0.030837420 -0.0149619996 5.941757e-06
# 50 (B.1.36+) - B.1.617.2 -0.106355242 0.004916705 NA -0.115991808 -0.0967186767 0.000000e+00
# 51     (B.1.36+) - other -0.040414830 0.001808231 NA -0.043958898 -0.0368707626 0.000000e+00
# 52     B.1.525 - B.1.351 -0.033203543 0.017189993 NA -0.066895310  0.0004882245 7.373554e-01
# 53     B.1.525 - B.1.618 -0.064368833 0.016773280 NA -0.097243857 -0.0314938088 1.003936e-02
# 54   B.1.525 - B.1.617.1 -0.057357599 0.015081378 NA -0.086916556 -0.0277986409 1.127868e-02
# 55   B.1.525 - B.1.617.2 -0.140813131 0.015415243 NA -0.171026453 -0.1105998093 1.568745e-13
# 56       B.1.525 - other -0.074872719 0.014897524 NA -0.104071330 -0.0456741093 9.942397e-05
# 57     B.1.351 - B.1.618 -0.031165291 0.012339906 NA -0.055351063 -0.0069795185 3.337467e-01
# 58   B.1.351 - B.1.617.1 -0.024154056 0.010111493 NA -0.043972218 -0.0043358937 4.210456e-01
# 59   B.1.351 - B.1.617.2 -0.107609589 0.010558499 NA -0.128303867 -0.0869153103 8.026912e-14
# 60       B.1.351 - other -0.041669177 0.009790760 NA -0.060858714 -0.0224796394 2.237984e-03
# 61   B.1.618 - B.1.617.1  0.007011235 0.008452295 NA -0.009554960  0.0235774290 9.995476e-01
# 62   B.1.618 - B.1.617.2 -0.076444298 0.008943639 NA -0.093973508 -0.0589150880 1.881051e-12
# 63       B.1.618 - other -0.010503886 0.007738975 NA -0.025671999  0.0046642258 9.693626e-01
# 64 B.1.617.1 - B.1.617.2 -0.083455533 0.005095159 NA -0.093441861 -0.0734692047 0.000000e+00
# 65     B.1.617.1 - other -0.017515121 0.003802441 NA -0.024967769 -0.0100624726 5.713438e-04
# 66     B.1.617.2 - other  0.065940412 0.004710899 NA  0.056707220  0.0751736032 0.000000e+00



# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-06-01"))
date.to = as.numeric(as.Date("2021-06-30")) # max(GISAID_india$DATE_NUM)+extrapolate

fit_india_multi_predsbystate = data.frame(emmeans(fit3_india_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "division"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                                  mode="prob", df=NA))
fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
fit_india_multi_predsbystate$LINEAGE2 = factor(fit_india_multi_predsbystate$LINEAGE2, levels=levels_LINEAGE2) 
fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$division, levels=levels_STATES) 

fit_india_multi_preds = data.frame(emmeans(fit3_india_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
                                           mode="prob", df=NA))
fit_india_multi_preds$collection_date = as.Date(fit_india_multi_preds$DATE_NUM, origin="1970-01-01")
fit_india_multi_preds$LINEAGE2 = factor(fit_india_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 

muller_india_mfit = ggplot(data=fit_india_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
           xmax=as.Date("2021-06-30"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
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
ggarrange(muller_india_raw2+coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-30")))+
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
           xmax=as.Date("2021-06-30"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
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
            coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-30")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0)))+
            ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\nRaw GISAID data"), 
          muller_indiabystate_mfit+ggtitle("\nMultinomial fit")+
            coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-30"))), nrow=2)

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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-30")), ylim=c(0.005, 0.95), expand=c(0,0))
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
  coord_cartesian(xlim=as.Date(c("2021-01-01","2021-06-30")),
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

