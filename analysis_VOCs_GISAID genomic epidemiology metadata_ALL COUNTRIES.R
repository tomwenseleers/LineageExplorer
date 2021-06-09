# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs ACROSS DIFFERENT COUNTRIES (GISAID GENOMIC EPIDEMIOLOGY METADATA)
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

# import GISAID genomic epidemiology metadata (file version metadata_2021-06-07_10-10.tsv.gz)
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-06-07_10-10.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)

GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
unique(GISAID$host)
# [1] "Human"                               "Environment"                         "Feline"                              "unknown"                            
# [5] "Rhinolophus shameli"                 "Rhinolophus malayanus"               "Rhinolophus pusillus"                "Rhinolophus sinicus"                
# [9] "Rhinolophus stheno"                  "Rhinolophus affinis"                 "Felis catus"                         "Canis lupus familiaris"             
# [13] "Gorilla gorilla gorilla"             "Mesocricetus auratus"                "Prionailurus bengalensis euptilurus" "Panthera leo"                       
# [17] "Mink"                                "Mustela putorius furo"               "Chlorocebus sabaeus"                 "Mus musculus"                       
# [21] "Mus musculus (BALB/c mice)"          "Manis javanica"                      "Manis pentadactyla"                  "Panthera tigris jacksoni" 

GISAID[GISAID$host!="Human","strain"]
GISAID = GISAID[GISAID$host=="Human",]
GISAID = GISAID[GISAID$date>=as.Date("2020-11-01"),] # we use data from Nov 1 2020 onwards
range(GISAID$date) # "2020-11-01" "2021-06-03"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 21 nov Uttar Pradesh, Varanasi India, 28 yr old male B.1.617.2

# GISAID = GISAID[grepl("2021-", GISAID$date),]
sum(is.na(GISAID$purpose_of_sequencing)) == nrow(GISAID) # field purpose_of_sequencing left blank unfortunately
nrow(GISAID) # 1743096
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

length(unique(GISAID$country[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 61 countries
table(GISAID$pango_lineage[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])
# B.1.617 B.1.617.1 B.1.617.2 B.1.617.3 
# 1      2921     21684        84

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
# 122932      Australia B.1.617+   193
# 122936        Bahrain B.1.617+    19
# 122937     Bangladesh B.1.617+    37
# 122940        Belgium B.1.617+   192
# 122963 Czech Republic B.1.617+    13
# 122965        Denmark B.1.617+   115
# 122974         France B.1.617+    93
# 122980        Germany B.1.617+   578
# 122994          India B.1.617+  4430
# 122995      Indonesia B.1.617+    32
# 122996           Iran B.1.617+    11
# 122998        Ireland B.1.617+   218
# 122999         Israel B.1.617+    36
# 123000          Italy B.1.617+   136
# 123002          Japan B.1.617+   165
# 123013     Luxembourg B.1.617+    58
# 123016       Malaysia B.1.617+    12
# 123019         Mexico B.1.617+    29
# 123028    Netherlands B.1.617+    56
# 123029    New Zealand B.1.617+    17
# 123033         Norway B.1.617+    42
# 123042         Poland B.1.617+    58
# 123043       Portugal B.1.617+    77
# 123044          Qatar B.1.617+    21
# 123046        Romania B.1.617+    19
# 123047         Russia B.1.617+    89
# 123055      Singapore B.1.617+   292
# 123059   South Africa B.1.617+    15
# 123061          Spain B.1.617+   105
# 123064         Sweden B.1.617+    42
# 123065    Switzerland B.1.617+    83
# 123067       Thailand B.1.617+    83
# 123076 United Kingdom B.1.617+ 16820
# 123078            USA B.1.617+  1880

table_country_lineage = as.data.frame(table(GISAID$country, GISAID$LINEAGE1))
colnames(table_country_lineage) = c("Country","Lineage","Count")
table_country_total = as.data.frame(table(GISAID$country))
table_country_lineage$Total = table_country_total$Freq[match(table_country_lineage$Country, table_country_total$Var1)]
table_country_lineage$Prop = table_country_lineage$Count / table_country_lineage$Total

table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage, fixed=T)&table_country_lineage$Prop>0.05,]
#                                 Country  Lineage Count  Total       Prop
# 122928                         Anguilla B.1.617+     1      2 0.50000000
# 122932                        Australia B.1.617+   193   1374 0.14046579
# 122936                          Bahrain B.1.617+    19     94 0.20212766
# 122956                            China B.1.617+     7     58 0.12068966
# 122964 Democratic Republic of the Congo B.1.617+     6     25 0.24000000
# 122979                          Georgia B.1.617+     4     44 0.09090909
# 122994                            India B.1.617+  4430  11180 0.39624329
# 122996                             Iran B.1.617+    11     59 0.18644068
# 123027                            Nepal B.1.617+     4     10 0.40000000
# 123029                      New Zealand B.1.617+    17    321 0.05295950
# 123055                        Singapore B.1.617+   292   1105 0.26425339
# 123067                         Thailand B.1.617+    83   1144 0.07255245
# 123076                   United Kingdom B.1.617+ 16820 331964 0.05066814
# 123080                          Vietnam B.1.617+     1     10 0.10000000

tblB1617 = table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage, fixed=T)&table_country_lineage$Count>10,]
tblB1617


sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>100,]$Country))
sel_countries_target
# [1] "Australia"      "Belgium"        "Denmark"        "Germany"        "India"          "Ireland"        "Italy"          "Japan"          "Singapore"     
# [10] "Spain"          "United Kingdom" "USA" 


sel_ref_lineage = "B.1.1.7"

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Belgium"        "Denmark"        "Germany"        "India"          "Ireland"        "Italy"          "Japan"          "Singapore"     
# [10] "Spain"          "United Kingdom" "USA"

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Belgium"        "Denmark"        "Germany"        "India"          "Ireland"        "Italy"          "Japan"          "Singapore"     
# [10] "Spain"          "United Kingdom" "USA"

sel_countries = sel_countries[!(sel_countries %in% c("Japan","USA","Singapore","Australia","India"))] # Japan is almost only import & for USA we do separate analysis by state

tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB117
#              Country Lineage  Count
# 77750      Australia B.1.1.7    388
# 77758        Belgium B.1.1.7  13532
# 77785        Denmark B.1.1.7  46633
# 77800        Germany B.1.1.7  86479
# 77814          India B.1.1.7   1050
# 77818        Ireland B.1.1.7  11146
# 77820          Italy B.1.1.7  18194
# 77822          Japan B.1.1.7  11349
# 77879      Singapore B.1.1.7    182
# 77902 United Kingdom B.1.1.7 244477
# 77904            USA B.1.1.7 153269

data.frame(Country=tblB1617$Country, Lineage="B.1.617", Perc=100*tblB1617$Count / (tblB1617$Count+tblB117$Count))
#           Country Lineage        Perc
# 1       Australia B.1.617 27.067669173
# 2         Bahrain B.1.617  0.140211055
# 3      Bangladesh B.1.617  0.042869698
# 4         Belgium B.1.617  0.155863948
# 5  Czech Republic B.1.617  1.222953904
# 6         Denmark B.1.617  0.933250378
# 7          France B.1.617  0.497675690
# 8         Germany B.1.617  3.576890399
# 9           India B.1.617 95.576081672
# 10      Indonesia B.1.617  0.013087453
# 11           Iran B.1.617  0.007176409
# 12        Ireland B.1.617 32.167832168
# 13         Israel B.1.617  0.265330189
# 14          Italy B.1.617  0.284394646
# 15          Japan B.1.617  0.190434421
# 16         Mexico B.1.617  2.597402597
# 17    Netherlands B.1.617  0.446588067
# 18    New Zealand B.1.617  0.082376847
# 19         Norway B.1.617  0.324960478
# 20         Poland B.1.617 18.385650224
# 21       Portugal B.1.617  0.031485889
# 22        Romania B.1.617  0.012394969
# 23      Singapore B.1.617 42.941176471
# 24   South Africa B.1.617  0.103351543
# 25          Spain B.1.617  0.184079283
# 26         Sweden B.1.617  0.020809970
# 27    Switzerland B.1.617  6.666666667
# 28 United Kingdom B.1.617 52.365485705
# 29            USA B.1.617  8.351803345





# ANALYSIS OF VOCs ACROSS THESE COUNTRIES WITH >100 B.1.617+ SEQUENCES ####
            
# GISAID_all = GISAID_all[GISAID_all$country_exposure=="India"&GISAID_all$country!="India",]
# nrow(GISAID_all[is.na(GISAID_all$LINEAGE1),]) # 0 unknown pango clade


GISAID_all = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_all) # 1316040
unique(GISAID_all$country)

rowSums(table(GISAID_all$LINEAGE1,GISAID_all$country))


GISAID_all = GISAID_all[!is.na(GISAID_all$LINEAGE1),]
nrow(GISAID_all) # 750541

GISAID_all = GISAID_all[GISAID_all$country==GISAID_all$country_exposure,] # we remove travel-related cases
nrow(GISAID_all) # 749220

sum(GISAID_all$LINEAGE1=="B.1.617+") # 17410
unique(GISAID_all$country[GISAID_all$LINEAGE1=="B.1.1.7"])
sum(GISAID_all$LINEAGE1=="B.1.1.7") # 421862

table(GISAID_all$LINEAGE1)
table(GISAID_all$LINEAGE2)

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those for all countries except India, Singapore & Australia
GISAID_all = GISAID_all[-which((grepl("B.1.617", GISAID_all$pango_lineage, fixed=TRUE)&(!GISAID_all$country %in% c("India","Singapore","Australia")))&
                                 GISAID_all$date<=as.Date("2021-04-14")),]  

main_lineages = names(table(GISAID_all$LINEAGE1))[100*table(GISAID_all$LINEAGE1)/sum(table(GISAID_all$LINEAGE1)) > 1]
main_lineages
# "B.1.1.7"  "B.1.160"  "B.1.177+" "B.1.617+"
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526")
main_lineages = union(main_lineages, VOCs)
GISAID_all$LINEAGE1[!(GISAID_all$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_all$LINEAGE2[!(GISAID_all$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove = names(table(GISAID_all$LINEAGE1))[table(GISAID_all$LINEAGE1)/sum(table(GISAID_all$LINEAGE1)) < 0.01]
remove = remove[!(remove %in% c("B.1.351","B.1.1.7","P.1","B.1.617.2"))]
GISAID_all$LINEAGE1[(GISAID_all$LINEAGE1 %in% remove)] = "other" # minority VOCs
GISAID_all$LINEAGE2[(GISAID_all$LINEAGE2 %in% remove)] = "other" # minority VOCs
table(GISAID_all$LINEAGE1)
#   B.1    B.1.1  B.1.1.7 B.1.177+  B.1.351 B.1.617+    other      P.1 
#  17781    20805   421862   119976     3867    17024   146023     1496 
GISAID_all$LINEAGE1 = factor(GISAID_all$LINEAGE1)
GISAID_all$LINEAGE1 = relevel(GISAID_all$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_all$LINEAGE1)
# "B.1.1.7"  "B.1.160"  "B.1.177+" "B.1.351"  "B.1.617+" "other"    "P.1"  
levels_LINEAGE1 = c("B.1.1.7","B.1.177+","B.1.160",
                    "B.1.351","P.1","B.1.617+","other")
GISAID_all$LINEAGE1 = factor(GISAID_all$LINEAGE1, levels=levels_LINEAGE1)

GISAID_all$LINEAGE2 = factor(GISAID_all$LINEAGE2)
GISAID_all$LINEAGE2 = relevel(GISAID_all$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_all$LINEAGE2)
# "B.1.1.7"   "B.1.160"   "B.1.177+"  "B.1.351"   "B.1.617.1" "B.1.617.2" "other"     "P.1" 
levels_LINEAGE2 = c("B.1.1.7","B.1.177+","B.1.160",
                    "B.1.351","P.1","B.1.617.1","B.1.617.2","other")
GISAID_all$LINEAGE2 = factor(GISAID_all$LINEAGE2, levels=levels_LINEAGE2)

# select countries with >100 imported sequences

# GISAID_all = GISAID_all[GISAID_all$division!="India",]
table(GISAID_all$country)
# Australia        Belgium        Denmark        Germany          India        Ireland          Italy      Singapore United Kingdom 
#     17798          23262         102013         119774          16175          14302          28407           1503         425600 

GISAID_all$country = factor(GISAID_all$country)
GISAID_all$country = relevel(GISAID_all$country, ref="United Kingdom")
levels_countries = c("United Kingdom","Ireland","Belgium","Italy","Spain","Germany","Denmark") # "India","Singapore","Australia",
GISAID_all$country = factor(GISAID_all$country, levels=levels_countries)
levels(GISAID_all$country)
table(GISAID_all$country)
sum(table(GISAID_all$country)) # 749220
table(GISAID_all$LINEAGE2)
# B.1.1.7       B.1  B.1.177+     B.1.1   B.1.351       P.1 B.1.617.1 B.1.617.2     other 
# 421862     17781    119976     20805      3867      1496      2169     14777    146101 

range(GISAID_all$date) # "2020-01-23" "2021-06-02"

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek1 = as.data.frame(table(GISAID_all$floor_date, GISAID_all$LINEAGE1))
colnames(data_agbyweek1) = c("floor_date", "LINEAGE1", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_all) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$LINEAGE1 = factor(data_agbyweek1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

data_agbyweek2 = as.data.frame(table(GISAID_all$floor_date, GISAID_all$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_all) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL


# aggregated by week and country for selected variant lineages
data_agbyweekregion1 = as.data.frame(table(GISAID_all$floor_date, GISAID_all$country, GISAID_all$LINEAGE1))
colnames(data_agbyweekregion1) = c("floor_date", "division", "LINEAGE1", "count")
data_agbyweekregion1_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$floor_date,data_agbyweekregion1$division), 
                                                                  interaction(data_agbyweekregion1_sum$floor_date,data_agbyweekregion1_sum$division))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_all) # correct
data_agbyweekregion1$collection_date = as.Date(as.character(data_agbyweekregion1$floor_date))
data_agbyweekregion1$LINEAGE1 = factor(data_agbyweekregion1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweekregion1$division = factor(data_agbyweekregion1$division, levels=levels_countries)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]
data_agbyweekregion1$floor_date = NULL

data_agbyweekregion2 = as.data.frame(table(GISAID_all$floor_date, GISAID_all$country, GISAID_all$LINEAGE2))
colnames(data_agbyweekregion2) = c("floor_date", "division", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$division), 
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$division))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_all) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$division = factor(data_agbyweekregion2$division, levels=levels_countries)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]
data_agbyweekregion2$floor_date = NULL


# MULLER PLOT (RAW DATA)
library(scales)
n1 = length(levels(GISAID_all$LINEAGE1))
lineage_cols1 = hcl(h = seq(15, 320, length = n1), l = 65, c = 200)
lineage_cols1[which(levels(GISAID_all$LINEAGE1)=="B.1.617+")] = "magenta"
lineage_cols1[which(levels(GISAID_all$LINEAGE1)=="other")] = "grey75"

n2 = length(levels(GISAID_all$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_all$LINEAGE2)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(GISAID_all$LINEAGE2)=="B.1.617.2")] = "magenta"
lineage_cols2[which(levels(GISAID_all$LINEAGE2)=="other")] = "grey75"

muller_sel_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE1)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1, group=LINEAGE1), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_all$Week)))
muller_sel_raw1

muller_sel_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_all$Week)))
muller_sel_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots_raw data.pdf"), width=8, height=6)


muller_selbycountry_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ division, nrow=3) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_all$Week)))
muller_selbycountry_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots by country_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots by country_raw data.pdf"), width=8, height=6)



# multinomial fits
data_agbyweekregion2$LINEAGE2 = relevel(data_agbyweekregion2$LINEAGE2, ref="B.1.1.7")
data_agbyweekregion2$DATE_NUM = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$COUNTRY = data_agbyweekregion2$division

library(nnet)
library(splines)
set.seed(1)
fit1_sel_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)+COUNTRY, weights=count, data=data_agbyweekregion2, maxit=1000)
fit2_sel_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)*COUNTRY, weights=count, data=data_agbyweekregion2, maxit=1000)
fit3_sel_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)+COUNTRY, weights=count, data=data_agbyweekregion2, maxit=1000)
fit4_sel_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)*COUNTRY, weights=count, data=data_agbyweekregion2, maxit=1000)
fit5_sel_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3)+COUNTRY, weights=count, data=data_agbyweekregion2, maxit=1000)
fit6_sel_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3)*COUNTRY, weights=count, data=data_agbyweekregion2, maxit=1000)
BIC(fit1_sel_multi, fit2_sel_multi, fit3_sel_multi, fit4_sel_multi, fit5_sel_multi, fit6_sel_multi) 
# fit6_sel_multi fits best (lowest BIC), but I will use the slightly simpler fit2_sel_multi


# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrsel = emtrends(fit2_sel_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_all$DATE_NUM)))
delta_r_sel = data.frame(confint(emtrsel, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrsel$contrasts)$p.value)
delta_r_sel
#              contrast      estimate         SE df   asymp.LCL  asymp.UCL   p.value
# 1 (B.1.177+) - B.1.1.7 -0.054431874 0.0003184658 NA -0.05505606 -0.053807692 2.642847e-10
# 2    B.1.160 - B.1.1.7 -0.057540768 0.0012113326 NA -0.05991494 -0.055166600 2.642847e-10
# 3    B.1.351 - B.1.1.7 -0.007881344 0.0010407316 NA -0.00992114 -0.005841547 4.072304e-10
# 4        P.1 - B.1.1.7  0.013077394 0.0014301431 NA  0.01027436  0.015880423 2.643800e-10
# 5  B.1.617.1 - B.1.1.7  0.052038335 0.0070343073 NA  0.03825135  0.065825324 5.979512e-10
# 6  B.1.617.2 - B.1.1.7  0.075807440 0.0020184120 NA  0.07185143  0.079763455 2.642847e-10
# 7      other - B.1.1.7 -0.038605673 0.0002437926 NA -0.03908350 -0.038127849 2.642847e-10


# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) per country 
emtrsel_bycountry = emtrends(fit2_sel_multi, trt.vs.ctrl ~ LINEAGE2, by=c("COUNTRY"), 
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_all$DATE_NUM)))
delta_r_sel_bycountry = data.frame(confint(emtrsel_bycountry, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrsel_bycountry$contrasts)$p.value)
delta_r_sel_bycountry
#             contrast        COUNTRY      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1  (B.1.177+) - B.1.1.7 United Kingdom -0.056303647 0.0002287789 NA -0.056752046 -0.055855249 2.642847e-10
# 2     B.1.160 - B.1.1.7 United Kingdom -0.063458042 0.0012937932 NA -0.065993830 -0.060922254 2.642847e-10
# 3     B.1.351 - B.1.1.7 United Kingdom  0.006134172 0.0010131746 NA  0.004148387  0.008119958 1.812466e-07
# 4         P.1 - B.1.1.7 United Kingdom  0.035640174 0.0027545872 NA  0.030241282  0.041039065 2.642857e-10
# 5   B.1.617.1 - B.1.1.7 United Kingdom  0.050442795 0.0026982750 NA  0.045154274  0.055731317 2.642847e-10
# 6   B.1.617.2 - B.1.1.7 United Kingdom  0.127085849 0.0010044825 NA  0.125117099  0.129054598 2.642847e-10
# 7       other - B.1.1.7 United Kingdom -0.048575012 0.0002826686 NA -0.049129033 -0.048020992 2.642847e-10
# 8  (B.1.177+) - B.1.1.7        Ireland -0.058200135 0.0016128681 NA -0.061361298 -0.055038971 2.642847e-10
# 9     B.1.160 - B.1.1.7        Ireland -0.068650907 0.0079447105 NA -0.084222254 -0.053079561 2.650662e-10
# 10    B.1.351 - B.1.1.7        Ireland -0.016087821 0.0034381584 NA -0.022826488 -0.009349155 6.316364e-05
# 11        P.1 - B.1.1.7        Ireland -0.007754429 0.0056439552 NA -0.018816378  0.003307520 5.822215e-01
# 12  B.1.617.1 - B.1.1.7        Ireland  0.077144184 0.0078180471 NA  0.061821093  0.092467274 2.643196e-10
# 13  B.1.617.2 - B.1.1.7        Ireland  0.063285495 0.0056530599 NA  0.052205701  0.074365288 2.643007e-10
# 14      other - B.1.1.7        Ireland -0.015220857 0.0011489636 NA -0.017472785 -0.012968930 2.642851e-10
# 15 (B.1.177+) - B.1.1.7        Belgium -0.045916294 0.0009857214 NA -0.047848273 -0.043984316 2.642847e-10
# 16    B.1.160 - B.1.1.7        Belgium -0.041903230 0.0010383951 NA -0.043938447 -0.039868013 2.642847e-10
# 17    B.1.351 - B.1.1.7        Belgium -0.018616122 0.0009265816 NA -0.020432189 -0.016800056 2.642847e-10
# 18        P.1 - B.1.1.7        Belgium  0.010135010 0.0009614110 NA  0.008250679  0.012019341 2.643112e-10
# 19  B.1.617.1 - B.1.1.7        Belgium  0.031212742 0.0200186501 NA -0.008023092  0.070448575 4.641095e-01
# 20  B.1.617.2 - B.1.1.7        Belgium  0.059455523 0.0043438642 NA  0.050941706  0.067969340 2.642849e-10
# 21      other - B.1.1.7        Belgium -0.033703770 0.0005970251 NA -0.034873918 -0.032533622 2.642847e-10
# 22 (B.1.177+) - B.1.1.7          Italy -0.049170768 0.0007272656 NA -0.050596182 -0.047745353 2.642847e-10
# 23    B.1.160 - B.1.1.7          Italy -0.046105066 0.0012998274 NA -0.048652681 -0.043557451 2.642847e-10
# 24    B.1.351 - B.1.1.7          Italy -0.015260929 0.0051954004 NA -0.025443727 -0.005078131 2.480503e-02
# 25        P.1 - B.1.1.7          Italy  0.013939408 0.0025232715 NA  0.008993887  0.018884930 1.889223e-06
# 26  B.1.617.1 - B.1.1.7          Italy  0.049682043 0.0245689855 NA  0.001527716  0.097836369 2.182713e-01
# 27  B.1.617.2 - B.1.1.7          Italy  0.074948247 0.0059956533 NA  0.063196983  0.086699512 2.642869e-10
# 28      other - B.1.1.7          Italy -0.024122219 0.0006113714 NA -0.025320485 -0.022923953 2.642847e-10
# 29 (B.1.177+) - B.1.1.7          Spain -0.045911048 0.0006359682 NA -0.047157523 -0.044664573 2.642847e-10
# 30    B.1.160 - B.1.1.7          Spain -0.040791920 0.0017814472 NA -0.044283492 -0.037300347 2.642847e-10
# 31    B.1.351 - B.1.1.7          Spain  0.012609054 0.0021425350 NA  0.008409763  0.016808346 3.870551e-07
# 32        P.1 - B.1.1.7          Spain  0.023571914 0.0015992723 NA  0.020437398  0.026706430 2.642847e-10
# 33  B.1.617.1 - B.1.1.7          Spain  0.082725343 0.0323840263 NA  0.019253818  0.146196868 6.801497e-02
# 34  B.1.617.2 - B.1.1.7          Spain  0.083000895 0.0074030169 NA  0.068491248  0.097510541 2.643015e-10
# 35      other - B.1.1.7          Spain -0.022251968 0.0006193729 NA -0.023465916 -0.021038019 2.642847e-10
# 36 (B.1.177+) - B.1.1.7        Germany -0.058743034 0.0004427850 NA -0.059610877 -0.057875191 2.642847e-10
# 37    B.1.160 - B.1.1.7        Germany -0.064760023 0.0007257545 NA -0.066182475 -0.063337570 2.642847e-10
# 38    B.1.351 - B.1.1.7        Germany -0.013000524 0.0007669679 NA -0.014503753 -0.011497294 2.642847e-10
# 39        P.1 - B.1.1.7        Germany  0.017871726 0.0030998543 NA  0.011796123  0.023947328 6.584282e-07
# 40  B.1.617.1 - B.1.1.7        Germany  0.043100289 0.0066231193 NA  0.030119214  0.056081364 2.298662e-08
# 41  B.1.617.2 - B.1.1.7        Germany  0.072426324 0.0032023453 NA  0.066149842  0.078702806 2.642847e-10
# 42      other - B.1.1.7        Germany -0.054266860 0.0003791359 NA -0.055009953 -0.053523767 2.642847e-10
# 43 (B.1.177+) - B.1.1.7        Denmark -0.066778189 0.0004635109 NA -0.067686654 -0.065869725 2.642847e-10
# 44    B.1.160 - B.1.1.7        Denmark -0.077116189 0.0007991162 NA -0.078682428 -0.075549950 2.642847e-10
# 45    B.1.351 - B.1.1.7        Denmark -0.010947238 0.0026825556 NA -0.016204951 -0.005689526 6.123625e-04
# 46        P.1 - B.1.1.7        Denmark -0.001862044 0.0064280788 NA -0.014460847  0.010736758 9.954055e-01
# 47  B.1.617.1 - B.1.1.7        Denmark  0.029960950 0.0161002370 NA -0.001594934  0.061516835 2.923249e-01
# 48  B.1.617.2 - B.1.1.7        Denmark  0.050449747 0.0068397811 NA  0.037044023  0.063855472 6.350708e-10
# 49      other - B.1.1.7        Denmark -0.072099026 0.0005046839 NA -0.073088188 -0.071109864 2.642847e-10

delta_r_sel_bycountry_B16172 = delta_r_sel_bycountry[delta_r_sel_bycountry$contrast=="B.1.617.2 - B.1.1.7",]
delta_r_sel_bycountry_B16172
#               contrast        COUNTRY   estimate          SE df  asymp.LCL  asymp.UCL      p.value
# 6  B.1.617.2 - B.1.1.7 United Kingdom 0.12708585 0.001004482 NA 0.12511710 0.12905460 2.642847e-10
# 13 B.1.617.2 - B.1.1.7        Ireland 0.06328549 0.005653060 NA 0.05220570 0.07436529 2.643007e-10
# 20 B.1.617.2 - B.1.1.7        Belgium 0.05945552 0.004343864 NA 0.05094171 0.06796934 2.642849e-10
# 27 B.1.617.2 - B.1.1.7          Italy 0.07494825 0.005995653 NA 0.06319698 0.08669951 2.642869e-10
# 34 B.1.617.2 - B.1.1.7          Spain 0.08300089 0.007403017 NA 0.06849125 0.09751054 2.643015e-10
# 41 B.1.617.2 - B.1.1.7        Germany 0.07242632 0.003202345 NA 0.06614984 0.07870281 2.642847e-10
# 48 B.1.617.2 - B.1.1.7        Denmark 0.05044975 0.006839781 NA 0.03704402 0.06385547 6.350708e-10


# plot of growth advantage of B.1.617.2 vs. B.1.1.7 across countries
qplot(data=delta_r_sel_bycountry_B16172, x=COUNTRY, y=estimate*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, geom="linerange", 
      colour=I("steelblue")) + 
      geom_point(colour=I("steelblue")) + xlab("") +
      ylab("Growth rate advantage of B.1.617.2 vs. B.1.1.7 (% per day)") + ggtitle("GROWTH RATE ADVANTAGE OF B.1.617.2 vs. B.1.1.7 BY COUNTRY\n(GISAID data, multinomial fit)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_growth rate advantages B16172.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_growth rate advantages B16172.pdf"), width=8, height=6)


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-11-01"))
date.to = as.numeric(as.Date("2021-06-30")) # max(GISAID_all$DATE_NUM)+extrapolate

fit_sel_multi_preds_bycountry = data.frame(emmeans(fit4_sel_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM","COUNTRY"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                           mode="prob", df=NA))
fit_sel_multi_preds_bycountry$collection_date = as.Date(fit_sel_multi_preds_bycountry$DATE_NUM, origin="1970-01-01")
fit_sel_multi_preds_bycountry$LINEAGE2 = factor(fit_sel_multi_preds_bycountry$LINEAGE2, levels=levels_LINEAGE2) 

muller_sel_mfit_bycountry = ggplot(data=fit_sel_multi_preds_bycountry, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ COUNTRY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_all$DATE_NUM)+1, 
           xmax=as.Date("2021-06-30"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\n(GISAID data, multinomial fit)")
muller_sel_mfit_bycountry

ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots_multinom fit by country.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots_multinom fit by country.pdf"), width=8, height=6)


library(ggpubr)
ggarrange(muller_selbycountry_raw2 + coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date("2021-06-30")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sel_mfit_bycountry+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots multipanel_multinom fit by country.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_muller plots multipanel_multinom fit by country.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

fit_sel_multi_preds2 = fit_sel_multi_preds_bycountry
# fit_sel_multi_preds2$LINEAGE2 = factor(fit_sel_multi_preds2$LINEAGE2, levels=levels_LINEAGE1)
# fit_sel_multi_preds2$LINEAGE1 = fit_sel_multi_preds2$LINEAGE2
# levels(fit_sel_multi_preds2$LINEAGE1)

# on logit scale:

fit_sel_multi_preds2 = fit_sel_multi_preds_bycountry
ymin = 0.001
ymax = 0.999
fit_sel_multi_preds2$asymp.LCL[fit_sel_multi_preds2$asymp.LCL<ymin] = ymin
fit_sel_multi_preds2$asymp.UCL[fit_sel_multi_preds2$asymp.UCL<ymin] = ymin
fit_sel_multi_preds2$asymp.UCL[fit_sel_multi_preds2$asymp.UCL>ymax] = ymax
fit_sel_multi_preds2$prob[fit_sel_multi_preds2$prob<ymin] = ymin

plot_sel_mfit_logit = qplot(data=fit_sel_multi_preds2, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~ COUNTRY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
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
                        range=c(0.5, 3), limits=c(1,max(data_agbyweekregion2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date("2021-06-30")), ylim=c(0.001, 0.99901), expand=c(0,0))
plot_sel_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_multinom fit by country_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_multinom fit by country_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_sel_mfit = qplot(data=fit_sel_multi_preds2, x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~ COUNTRY) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-11-01","2021-06-30")),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweekregion2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 3), limits=c(1,max(data_agbyweekregion2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_sel_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_multinom fit by country_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\all countries_multinom fit by country_response scale.pdf"), width=8, height=6)

