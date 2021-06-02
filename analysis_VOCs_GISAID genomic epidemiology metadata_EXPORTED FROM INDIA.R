# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN CASES EXPORTED FROM INDIA (GISAID GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 29 MAY 2021

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
range(GISAID$date) # "2021-06-01" "2021-05-23"
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

table_country_lineage = as.data.frame(table(GISAID$country, GISAID$LINEAGE1))
colnames(table_country_lineage) = c("Country","Lineage","Count")
tblB1617 = table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage, fixed=T)&table_country_lineage$Count>10,]
tblB1617
#              Country  Lineage Count
# 151209      Australia B.1.617+   144
# 151213        Bahrain B.1.617+    19
# 151214     Bangladesh B.1.617+    20
# 151217        Belgium B.1.617+   120
# 151241        Denmark B.1.617+   105
# 151251         France B.1.617+    86
# 151256        Germany B.1.617+   420
# 151270          India B.1.617+  3577
# 151271      Indonesia B.1.617+    32
# 151274        Ireland B.1.617+   155
# 151275         Israel B.1.617+    36
# 151276          Italy B.1.617+    81
# 151278          Japan B.1.617+   165
# 151296         Mexico B.1.617+    27
# 151305    Netherlands B.1.617+    50
# 151306    New Zealand B.1.617+    15
# 151310         Norway B.1.617+    20
# 151319         Poland B.1.617+    34
# 151320       Portugal B.1.617+    52
# 151323        Romania B.1.617+    19
# 151332      Singapore B.1.617+   156
# 151337   South Africa B.1.617+    14
# 151339          Spain B.1.617+    65
# 151342         Sweden B.1.617+    17
# 151343    Switzerland B.1.617+    64
# 151354 United Kingdom B.1.617+  9273
# 151356            USA B.1.617+  1544

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>10,]$Country))
sel_countries_target
# [1] "Australia"      "Bahrain"        "Bangladesh"     "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Indonesia"     
# [10] "Ireland"        "Israel"         "Italy"          "Japan"          "Mexico"         "Netherlands"    "New Zealand"    "Norway"         "Poland"        
# [19] "Portugal"       "Romania"        "Singapore"      "South Africa"   "Spain"          "Sweden"         "Switzerland"    "United Kingdom" "USA"   

sel_ref_lineage = "B.1.1.7"

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Bahrain"        "Bangladesh"     "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Indonesia"     
# [10] "Ireland"        "Israel"         "Italy"          "Japan"          "Mexico"         "Netherlands"    "New Zealand"    "Norway"         "Poland"        
# [19] "Portugal"       "Romania"        "Singapore"      "South Africa"   "Spain"          "Sweden"         "Switzerland"    "United Kingdom" "USA"  

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Bahrain"        "Bangladesh"     "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Indonesia"     
# [10] "Ireland"        "Israel"         "Italy"          "Japan"          "Mexico"         "Netherlands"    "New Zealand"    "Norway"         "Poland"        
# [19] "Portugal"       "Romania"        "Singapore"      "South Africa"   "Spain"          "Sweden"         "Switzerland"    "United Kingdom" "USA"  


tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB117
#              Country Lineage  Count
# 70569      Australia B.1.1.7    388
# 70573        Bahrain B.1.1.7     12
# 70574     Bangladesh B.1.1.7     73
# 70577        Belgium B.1.1.7  13008
# 70601        Denmark B.1.1.7  46633
# 70611         France B.1.1.7  22532
# 70616        Germany B.1.1.7  86477
# 70630          India B.1.1.7   1032
# 70631      Indonesia B.1.1.7     23
# 70634        Ireland B.1.1.7  10655
# 70635         Israel B.1.1.7   7544
# 70636          Italy B.1.1.7  17674
# 70638          Japan B.1.1.7   9020
# 70656         Mexico B.1.1.7    369
# 70665    Netherlands B.1.1.7  20280
# 70666    New Zealand B.1.1.7    133
# 70670         Norway B.1.1.7   5751
# 70679         Poland B.1.1.7  10886
# 70680       Portugal B.1.1.7   3947
# 70683        Romania B.1.1.7    501
# 70692      Singapore B.1.1.7    170
# 70697   South Africa B.1.1.7     30
# 70699          Spain B.1.1.7  12198
# 70702         Sweden B.1.1.7  41027
# 70703    Switzerland B.1.1.7  17363
# 70714 United Kingdom B.1.1.7 242202
# 70716            USA B.1.1.7 150809

data.frame(Country=tblB1617$Country, Lineage="B.1.617", Perc=100*tblB1617$Count / (tblB1617$Count+tblB117$Count))
#           Country Lineage        Perc
# 1       Australia B.1.617 27.06766917
# 2         Bahrain B.1.617 61.29032258
# 3      Bangladesh B.1.617 21.50537634
# 4         Belgium B.1.617  0.91407678
# 5         Denmark B.1.617  0.22465660
# 6          France B.1.617  0.38022814
# 7         Germany B.1.617  0.48333084
# 8           India B.1.617 77.60902582
# 9       Indonesia B.1.617 58.18181818
# 10        Ireland B.1.617  1.43385754
# 11         Israel B.1.617  0.47493404
# 12          Italy B.1.617  0.45620952
# 13          Japan B.1.617  1.79640719
# 14         Mexico B.1.617  6.81818182
# 15    Netherlands B.1.617  0.24594196
# 16    New Zealand B.1.617 10.13513514
# 17         Norway B.1.617  0.34656039
# 18         Poland B.1.617  0.31135531
# 19       Portugal B.1.617  1.30032508
# 20        Romania B.1.617  3.65384615
# 21      Singapore B.1.617 47.85276074
# 22   South Africa B.1.617 31.81818182
# 23          Spain B.1.617  0.53004974
# 24         Sweden B.1.617  0.04141897
# 25    Switzerland B.1.617  0.36724623
# 26 United Kingdom B.1.617  3.68744408
# 27            USA B.1.617  1.01343590


GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_sel) # 1458602

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))




# ANALYSIS VOCs EXPORTED FROM INDIA ####
            
GISAID_indiaexp = GISAID_sel[GISAID_sel$country_exposure=="India"&GISAID_sel$country!="India",]
nrow(GISAID_indiaexp[is.na(GISAID_indiaexp$LINEAGE1),]) # 0 unknown pango clade
GISAID_indiaexp = GISAID_indiaexp[!is.na(GISAID_indiaexp$LINEAGE1),]
nrow(GISAID_indiaexp) # 486

unique(GISAID_indiaexp$division_exposure) # just says "India"
unique(GISAID_indiaexp$country[GISAID_indiaexp$LINEAGE1=="B.1.617+"]) # "Bangladesh"   "France"       "Italy"        "Japan"        "Singapore"    "South Africa" "Spain"        "USA" 
table(GISAID_indiaexp$country, GISAID_indiaexp$LINEAGE2)

sum(GISAID_indiaexp$LINEAGE1=="B.1.617+") # 223
unique(GISAID_indiaexp$country[GISAID_indiaexp$LINEAGE1=="B.1.1.7"]) # Italy Japan Singapore
sum(GISAID_indiaexp$LINEAGE1=="B.1.1.7") # 86

table(GISAID_indiaexp$LINEAGE1)
table(GISAID_indiaexp$LINEAGE2)

main_lineages = names(table(GISAID_indiaexp$LINEAGE1))[100*table(GISAID_indiaexp$LINEAGE1)/sum(table(GISAID_indiaexp$LINEAGE1)) > 3]
main_lineages
# "B.1.1"     "B.1.1.354" "B.1.1.7"   "B.1.36+"   "B.1.617+" 
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526")
main_lineages = union(main_lineages, VOCs)
GISAID_indiaexp$LINEAGE1[!(GISAID_indiaexp$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_indiaexp$LINEAGE2[!(GISAID_indiaexp$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove = names(table(GISAID_indiaexp$LINEAGE1))[table(GISAID_indiaexp$LINEAGE1) < 10]
GISAID_indiaexp$LINEAGE1[(GISAID_indiaexp$LINEAGE1 %in% remove)] = "other" # minority VOCs
GISAID_indiaexp$LINEAGE2[(GISAID_indiaexp$LINEAGE2 %in% remove)] = "other" # minority VOCs
table(GISAID_indiaexp$LINEAGE1)
# B.1     B.1.1 B.1.1.354   B.1.1.7   B.1.351   B.1.36+  B.1.617+     other 
# 15        43        25        86        16        43       223        35 
GISAID_indiaexp$LINEAGE1 = factor(GISAID_indiaexp$LINEAGE1)
GISAID_indiaexp$LINEAGE1 = relevel(GISAID_indiaexp$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_indiaexp$LINEAGE1)
# "B.1.1.7"   "B.1.1"     "B.1.1.354" "B.1.351"   "B.1.36+"   "B.1.617+"  "other"  
levels_LINEAGE1 = c("B.1.1.7","B.1.1","B.1.1.354","B.1.36+",
                    "B.1.351","B.1.617+","other")
GISAID_indiaexp$LINEAGE1 = factor(GISAID_indiaexp$LINEAGE1, levels=levels_LINEAGE1)

GISAID_indiaexp$LINEAGE2 = factor(GISAID_indiaexp$LINEAGE2)
GISAID_indiaexp$LINEAGE2 = relevel(GISAID_indiaexp$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_indiaexp$LINEAGE2)
# "B.1.1.7"   "B.1.1"     "B.1.1.354" "B.1.351"   "B.1.36+"   "B.1.617.1" "B.1.617.2" "other"  
levels_LINEAGE2 = c("B.1.1.7","B.1.1","B.1.1.354","B.1.36+",
                    "B.1.351","B.1.617.1","B.1.617.2","other")
GISAID_indiaexp$LINEAGE2 = factor(GISAID_indiaexp$LINEAGE2, levels=levels_LINEAGE2)

# select countries with >100 imported sequences

# GISAID_indiaexp = GISAID_indiaexp[GISAID_indiaexp$division!="India",]
table(GISAID_indiaexp$country)
# Bangladesh       France        Italy        Japan  New Zealand    Singapore South Africa        Spain          USA 
#         21            3           11          132            3          303           10            1            2 

sel_countries = names(table(GISAID_indiaexp$country))[table(GISAID_indiaexp$country) >= 100]
sel_countries # "Japan"     "Singapore"

GISAID_indiaexp = GISAID_indiaexp[GISAID_indiaexp$country %in% sel_countries,]
GISAID_indiaexp$country = factor(GISAID_indiaexp$country)
GISAID_indiaexp$country = relevel(GISAID_indiaexp$country, ref="Singapore")
levels_countries = c("Singapore","Japan")
GISAID_indiaexp$country = factor(GISAID_indiaexp$country, levels=levels_countries)
levels(GISAID_indiaexp$country)
table(GISAID_indiaexp$country)
sum(table(GISAID_indiaexp$country)) # 435
table(GISAID_indiaexp$LINEAGE2)
range(GISAID_indiaexp$date) # "2020-06-17" "2021-05-18"

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek1 = as.data.frame(table(GISAID_indiaexp$floor_date, GISAID_indiaexp$LINEAGE1))
colnames(data_agbyweek1) = c("floor_date", "LINEAGE1", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_indiaexp) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$LINEAGE1 = factor(data_agbyweek1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

data_agbyweek2 = as.data.frame(table(GISAID_indiaexp$floor_date, GISAID_indiaexp$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_indiaexp) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL


# aggregated by week and country for selected variant lineages
data_agbyweekregion1 = as.data.frame(table(GISAID_indiaexp$floor_date, GISAID_indiaexp$country, GISAID_indiaexp$LINEAGE1))
colnames(data_agbyweekregion1) = c("floor_date", "division", "LINEAGE1", "count")
data_agbyweekregion1_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$floor_date,data_agbyweekregion1$division), 
                                                                  interaction(data_agbyweekregion1_sum$floor_date,data_agbyweekregion1_sum$division))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_indiaexp) # correct
data_agbyweekregion1$collection_date = as.Date(as.character(data_agbyweekregion1$floor_date))
data_agbyweekregion1$LINEAGE1 = factor(data_agbyweekregion1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweekregion1$division = factor(data_agbyweekregion1$division, levels=levels_countries)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]
data_agbyweekregion1$floor_date = NULL

data_agbyweekregion2 = as.data.frame(table(GISAID_indiaexp$floor_date, GISAID_indiaexp$country, GISAID_indiaexp$LINEAGE2))
colnames(data_agbyweekregion2) = c("floor_date", "division", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$division), 
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$division))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_indiaexp) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$division = factor(data_agbyweekregion2$division, levels=levels_country)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]
data_agbyweekregion2$floor_date = NULL


# MULLER PLOT (RAW DATA)
library(scales)
n1 = length(levels(GISAID_indiaexp$LINEAGE1))
lineage_cols1 = hcl(h = seq(15, 320, length = n1), l = 65, c = 200)
lineage_cols1[which(levels(GISAID_indiaexp$LINEAGE1)=="B.1.617+")] = "magenta"
lineage_cols1[which(levels(GISAID_indiaexp$LINEAGE1)=="other")] = "grey75"

n2 = length(levels(GISAID_indiaexp$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_indiaexp$LINEAGE2)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(GISAID_indiaexp$LINEAGE2)=="B.1.617.2")] = "magenta"
lineage_cols2[which(levels(GISAID_indiaexp$LINEAGE2)=="other")] = "grey75"

muller_indiaexp_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE1)) + 
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANT B.1.617+ EXPORTED FROM INDIA") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_indiaexp$Week)))
muller_indiaexp_raw1

muller_indiaexp_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2\nAMONG CASES EXPORTED FROM INDIA (GISAID data Singapore & Japan)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_indiaexp$Week)))
muller_indiaexp_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots_raw data.pdf"), width=8, height=6)


muller_indiabystate_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ division, nrow=3) +
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2\nAMONG CASES EXPORTED FROM INDIA (GISAID data Singapore & Japan)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_indiaexp$Week)))
muller_indiabystate_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots by country_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots by country_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_indiaexp_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM), data=GISAID_indiaexp, maxit=1000)
fit2_indiaexp_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2), data=GISAID_indiaexp, maxit=1000)
fit3_indiaexp_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3), data=GISAID_indiaexp, maxit=1000)
BIC(fit1_indiaexp_multi, fit2_indiaexp_multi, fit3_indiaexp_multi) 
# fit1_indiaexp_multi fits best (lowest BIC)

# equivalent fit with B.1.617.1,2&3 all recoded to B.1.617+
fit1_indiaexp_multi1 = nnet::multinom(LINEAGE1 ~ scale(DATE_NUM), data=GISAID_indiaexp, maxit=1000)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
max(GISAID_indiaexp$date) # 2021-05-18
emtrindiaexp = emtrends(fit1_indiaexp_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp = data.frame(confint(emtrindiaexp, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrindiaexp$contrasts)$p.value)
delta_r_indiaexp
#              contrast      estimate         SE df   asymp.LCL  asymp.UCL   p.value
# 1     B.1.1 - B.1.1.7 -0.05097728 0.007586084 NA -0.0658457367 -0.03610883 6.067628e-05
# 2 B.1.1.354 - B.1.1.7 -0.05723849 0.007992997 NA -0.0729044757 -0.04157250 3.024051e-05
# 3 (B.1.36+) - B.1.1.7 -0.05378562 0.007685445 NA -0.0688488115 -0.03872242 3.898059e-05
# 4   B.1.351 - B.1.1.7  0.05933189 0.030197494 NA  0.0001458888  0.11851789 2.886990e-01
# 5 B.1.617.1 - B.1.1.7  0.03778126 0.012741647 NA  0.0128080945  0.06275443 5.316742e-02
# 6 B.1.617.2 - B.1.1.7  0.13846057 0.017180158 NA  0.1047880824  0.17213307 7.900114e-06
# 7     other - B.1.1.7 -0.04810091 0.007593806 NA -0.0629844988 -0.03321733 1.138353e-04

# avg growth advantage of B.1.617+ over B.1.1.7 :
max(GISAID_indiaexp$date) # 2021-05-18
emtrindiaexp1 = emtrends(fit1_indiaexp_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp1 = data.frame(confint(emtrindiaexp1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrindiaexp1$contrasts)$p.value)
delta_r_indiaexp1
#         contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1      B.1.1 - B.1.1.7 -0.05021939 0.007487918 NA -0.06489544 -0.03554334 1.142863e-04
# 2  B.1.1.354 - B.1.1.7 -0.05648026 0.007899235 NA -0.07196248 -0.04099804 6.146195e-05
# 3  (B.1.36+) - B.1.1.7 -0.05303216 0.007588416 NA -0.06790519 -0.03815914 7.682811e-05
# 4    B.1.351 - B.1.1.7  0.04990058 0.026964275 NA -0.00294843  0.10274958 3.200341e-01
# 5 (B.1.617+) - B.1.1.7  0.08607712 0.012412628 NA  0.06174882  0.11040543 8.282994e-05
# 6      other - B.1.1.7 -0.04734199 0.007496183 NA -0.06203424 -0.03264974 2.016813e-04


# pairwise growth rate difference (differences in growth rate per day) 
max(GISAID_indiaexp$date) # 2021-05-18
emtrindiaexp_pairw = emtrends(fit1_indiaexp_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp_pairw = data.frame(confint(emtrindiaexp_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindiaexp_pairw$contrasts)$p.value)
delta_r_indiaexp_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1        B.1.1.7 - B.1.1  0.050977284 0.007586084 NA  0.0361088320  0.0658457367 1.983622e-04
# 2    B.1.1.7 - B.1.1.354  0.057238490 0.007992997 NA  0.0415725049  0.0729044757 9.975035e-05
# 3    B.1.1.7 - (B.1.36+)  0.053785616 0.007685445 NA  0.0387224205  0.0688488115 1.281817e-04
# 4      B.1.1.7 - B.1.351 -0.059331890 0.030197494 NA -0.1185178916 -0.0001458888 5.342010e-01
# 5    B.1.1.7 - B.1.617.1 -0.037781263 0.012741647 NA -0.0627544316 -0.0128080945 1.312506e-01
# 6    B.1.1.7 - B.1.617.2 -0.138460574 0.017180158 NA -0.1721330651 -0.1047880824 2.642337e-05
# 7        B.1.1.7 - other  0.048100913 0.007593806 NA  0.0332173270  0.0629844988 3.686705e-04
# 8      B.1.1 - B.1.1.354  0.006261206 0.003489584 NA -0.0005782524  0.0131006642 6.332679e-01
# 9      B.1.1 - (B.1.36+)  0.002808332 0.002916199 NA -0.0029073129  0.0085239761 9.728385e-01
# 10       B.1.1 - B.1.351 -0.110309175 0.031024035 NA -0.1711151665 -0.0495031826 4.767623e-02
# 11     B.1.1 - B.1.617.1 -0.088758547 0.014391916 NA -0.1169661844 -0.0605509105 4.848098e-04
# 12     B.1.1 - B.1.617.2 -0.189437858 0.018656428 NA -0.2260037851 -0.1528719313 1.675353e-06
# 13         B.1.1 - other -0.002876371 0.003175047 NA -0.0090993489  0.0033466059 9.804590e-01
# 14 B.1.1.354 - (B.1.36+) -0.003452874 0.003478457 NA -0.0102705238  0.0033647752 9.681427e-01
# 15   B.1.1.354 - B.1.351 -0.116570380 0.031126994 NA -0.1775781675 -0.0555625934 3.403988e-02
# 16 B.1.1.354 - B.1.617.1 -0.095019753 0.014612811 NA -0.1236603360 -0.0663791707 2.806939e-04
# 17 B.1.1.354 - B.1.617.2 -0.195699064 0.018826840 NA -0.2325989922 -0.1587991359 1.254298e-06
# 18     B.1.1.354 - other -0.009137577 0.003779385 NA -0.0165450354 -0.0017301193 3.038122e-01
# 19   (B.1.36+) - B.1.351 -0.113117506 0.031049264 NA -0.1739729463 -0.0522620660 4.081829e-02
# 20 (B.1.36+) - B.1.617.1 -0.091566879 0.014446387 NA -0.1198812776 -0.0632524805 3.661774e-04
# 21 (B.1.36+) - B.1.617.2 -0.192246190 0.018698123 NA -0.2288938371 -0.1555985424 1.436137e-06
# 22     (B.1.36+) - other -0.005684703 0.003232466 NA -0.0120202195  0.0006508133 6.539721e-01
# 23   B.1.351 - B.1.617.1  0.021550627 0.030616973 NA -0.0384575381  0.0815587923 9.954356e-01
# 24   B.1.351 - B.1.617.2 -0.079128684 0.031108610 NA -0.1401004379 -0.0181569293 2.536195e-01
# 25       B.1.351 - other  0.107432803 0.031024151 NA  0.0466265842  0.1682390220 5.614076e-02
# 26 B.1.617.1 - B.1.617.2 -0.100679311 0.016943954 NA -0.1338888498 -0.0674697716 7.051872e-04
# 27     B.1.617.1 - other  0.085882176 0.014391950 NA  0.0576744730  0.1140898790 6.757439e-04
# 28     B.1.617.2 - other  0.186561487 0.018657078 NA  0.1499942850  0.2231286884 2.023200e-06



# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-06-01"))
date.to = as.numeric(as.Date("2021-06-14")) # max(GISAID_indiaexp$DATE_NUM)+extrapolate

fit_indiaexp_multi_preds = data.frame(emmeans(fit1_indiaexp_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=1)), 
                                           mode="prob", df=NA))
fit_indiaexp_multi_preds$collection_date = as.Date(fit_indiaexp_multi_preds$DATE_NUM, origin="1970-01-01")
fit_indiaexp_multi_preds$LINEAGE2 = factor(fit_indiaexp_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 

muller_indiaexp_mfit = ggplot(data=fit_indiaexp_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_indiaexp$DATE_NUM)+1, 
           xmax=as.Date("2021-06-14"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-14",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2\nAMONG CASES EXPORTED FROM INDIA (GISAID data Singapore & Japan)")
muller_indiaexp_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots_multinom fit.pdf"), width=8, height=6)


library(ggpubr)
ggarrange(muller_indiaexp_raw2 + coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-14")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_indiaexp_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots multipanel_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots multipanel_multinom fit.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

fit_indiaexp_multi_preds2 = fit_indiaexp_multi_preds
fit_indiaexp_multi_preds2$LINEAGE2 = factor(fit_indiaexp_multi_preds2$LINEAGE2, levels=levels_LINEAGE1)
fit_indiaexp_multi_preds2$LINEAGE1 = fit_indiaexp_multi_preds2$LINEAGE2
levels(fit_indiaexp_multi_preds2$LINEAGE1)

# on logit scale:

fit_indiaexp_multi_preds2 = fit_indiaexp_multi_preds
ymin = 0.001
ymax = 0.998
fit_indiaexp_multi_preds2$asymp.LCL[fit_indiaexp_multi_preds2$asymp.LCL<ymin] = ymin
fit_indiaexp_multi_preds2$asymp.UCL[fit_indiaexp_multi_preds2$asymp.UCL<ymin] = ymin
fit_indiaexp_multi_preds2$asymp.UCL[fit_indiaexp_multi_preds2$asymp.UCL>ymax] = ymax
fit_indiaexp_multi_preds2$prob[fit_indiaexp_multi_preds2$prob<ymin] = ymin

plot_indiaexp_mfit_logit = qplot(data=fit_indiaexp_multi_preds2, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2\nAMONG CASES EXPORTED FROM INDIA (multinomial fit)") +
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
                        range=c(0.001, 5), limits=c(1,10^3), breaks=c(10,100)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-14")), ylim=c(0.005, 0.95), expand=c(0,0))
plot_indiaexp_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_multinom fit_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_multinom fit_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_indiaexp_mfit = qplot(data=fit_indiaexp_multi_preds, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2\nAMONG CASES EXPORTED FROM INDIA (multinomial fit)") +
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
plot_indiaexp_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_multinom fit_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_multinom fit_response scale.pdf"), width=8, height=6)

