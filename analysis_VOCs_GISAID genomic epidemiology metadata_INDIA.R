# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN INDIA BASED ON ANALYSIS OF GISAID GENOMIC EPIDEMIOLOGY METADATA
# T. Wenseleers
# last update 16 MAY 2021

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
today = as.Date("2021-05-17")
today_num = as.numeric(today)
today # "2021-05-17"
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata (file version metadata_2021-05-14_13-02.tsv.gz)
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
range(GISAID$date) # "2019-12-24" "2021-05-08"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 12 dec Bihar India

# GISAID = GISAID[grepl("2021-", GISAID$date),]
GISAID = GISAID[GISAID$date>=as.Date("2020-06-01"),]
sum(is.na(GISAID$purpose_of_sequencing)) == nrow(GISAID) # field purpose_of_sequencing left blank unfortunately
range(GISAID$date) # "2021-06-01" "2021-05-08"
nrow(GISAID) # 1345739 | 906318
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

length(unique(GISAID$country[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 46 countries
table(GISAID$pango_lineage[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])
# B.1.617 B.1.617.1 B.1.617.2 B.1.617.3 
# 1      1968      2490        63 

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
# 144153      Australia B.1.617+    84
# 144157        Bahrain B.1.617+    22
# 144161        Belgium B.1.617+    24
# 144185        Denmark B.1.617+    39
# 144195         France B.1.617+    16
# 144200        Germany B.1.617+   103
# 144210          India B.1.617+  1824
# 144214        Ireland B.1.617+    36
# 144216          Italy B.1.617+    29
# 144218          Japan B.1.617+    29
# 144245    Netherlands B.1.617+    12
# 144246    New Zealand B.1.617+    13
# 144271      Singapore B.1.617+   156
# 144281    Switzerland B.1.617+    32
# 144292 United Kingdom B.1.617+  1532
# 144294            USA B.1.617+   454

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>10,]$Country))
sel_countries_target
# [1] "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "Italy"         
# [10] "Japan"          "Netherlands"    "New Zealand"    "Singapore"      "Switzerland"    "United Kingdom" "USA" 

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "Italy"         
# [10] "Japan"          "Netherlands"    "New Zealand"    "Singapore"      "Switzerland"    "United Kingdom" "USA"    

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "Italy"         
# [10] "Japan"          "Netherlands"    "New Zealand"    "Singapore"      "Switzerland"    "United Kingdom" "USA"     

sel_ref_lineage = "B.1.1.7"
tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB117
#              Country Lineage  Count
# 66845      Australia B.1.1.7    351
# 66849        Bahrain B.1.1.7     12
# 66853        Belgium B.1.1.7  10103
# 66877        Denmark B.1.1.7  31844
# 66887         France B.1.1.7  16708
# 66892        Germany B.1.1.7  66931
# 66902          India B.1.1.7    653
# 66906        Ireland B.1.1.7   8981
# 66908          Italy B.1.1.7  13694
# 66910          Japan B.1.1.7   2433
# 66937    Netherlands B.1.1.7  14504
# 66938    New Zealand B.1.1.7    134
# 66963      Singapore B.1.1.7    170
# 66973    Switzerland B.1.1.7  13044
# 66984 United Kingdom B.1.1.7 229141
# 66986            USA B.1.1.7  93577

data.frame(Country=tblB1617$Country, Lineage="B.1.617", Perc=100*tblB1617$Count / (tblB1617$Count+tblB117$Count))
#           Country Lineage        Perc
# 1       Australia B.1.617 19.31034483
# 2         Bahrain B.1.617 64.70588235
# 3         Belgium B.1.617  0.23699022
# 4         Denmark B.1.617  0.12232224
# 5          France B.1.617  0.09567089
# 6         Germany B.1.617  0.15365337
# 7           India B.1.617 73.63746468
# 8         Ireland B.1.617  0.39924587
# 9           Italy B.1.617  0.21132405
# 10          Japan B.1.617  1.17790414
# 11    Netherlands B.1.617  0.08266740
# 12    New Zealand B.1.617  8.84353741
# 13      Singapore B.1.617 47.85276074
# 14    Switzerland B.1.617  0.24472316
# 15 United Kingdom B.1.617  0.66414361
# 16            USA B.1.617  0.48281950


GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_sel) # 1142699 | 763637

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))




# ANALYSIS VOCs IN INDIA ####

GISAID_india = GISAID_sel[GISAID_sel$country=="India",]
nrow(GISAID_india[is.na(GISAID_india$LINEAGE1),]) # 97 unknown pango clade
GISAID_india = GISAID_india[!is.na(GISAID_india$LINEAGE1),]
nrow(GISAID_india) # 9802 | 4961

unique(GISAID_india$division) # best data for West Bengal, Maharashtra & Karnataka (B.1.617 most common, B.1.1.7 most common in Punjab IN, Telangana)
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.617+"])
sum(GISAID_india$LINEAGE1=="B.1.617+") # 1824
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.1.7"])
sum(GISAID_india$LINEAGE1=="B.1.1.7") # 653

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
firstB16172

# select states of India with a total of > 300 sequences submitted

GISAID_india = GISAID_india[GISAID_india$division!="India",]
table(GISAID_india$division)

sel_states = names(table(GISAID_india$division))[table(GISAID_india$division) > 300]
sel_states
# "Andhra Pradesh" "Chhattisgarh"   "Delhi"          "Gujarat"        "Karnataka"      "Maharashtra"    "Odisha"     "Telangana"      "West Bengal" 
GISAID_india = GISAID_india[GISAID_india$division %in% sel_states,]
GISAID_india$division = factor(GISAID_india$division)
GISAID_india$division = relevel(GISAID_india$division, ref="Maharashtra")
# levels_STATES = c("Maharashtra","Karnataka","West Bengal")
levels_STATES = c("Maharashtra","Chhattisgarh","Gujarat","Delhi","Karnataka","Andhra Pradesh","West Bengal","Telangana","Odisha")
GISAID_india$division = factor(GISAID_india$division, levels=levels_STATES)
# levels_STATES = levels(GISAID_india$division)
levels(GISAID_india$division)
# "Maharashtra"    "Chhattisgarh"   "Gujarat"        "Delhi"          "Karnataka"      "Andhra Pradesh" "West Bengal"    "Telangana"      "Odisha"   
table(GISAID_india$division)
# Maharashtra   Chhattisgarh        Gujarat          Delhi      Karnataka Andhra Pradesh    West Bengal      Telangana         Odisha 
# 2859            334           1294            323            806            609           1176           1367            309 
sum(table(GISAID_india$division)) # 9077
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
# fit5_india_multi fits best (lowest BIC)

# equivalent fit with B.1.617.1,2&3 all recoded to B.1.617+
fit5_india_multi1 = nnet::multinom(LINEAGE1 ~ division + ns(DATE_NUM, df=3), data=GISAID_india, maxit=1000)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
max(GISAID_india$date) # 2021-04-30
emtrindia = emtrends(fit5_india_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india = data.frame(confint(emtrindia, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrindia$contrasts)$p.value)
delta_r_india
#               contrast     estimate           SE df    asymp.LCL    asymp.UCL p.value
# 1        B.1 - B.1.1.7  0.022614779 0.006960121 NA  0.008973192  0.036256366 1.476186e-02
# 2      B.1.1 - B.1.1.7 -0.036957533 0.008092518 NA -0.052818577 -0.021096489 1.223022e-04
# 3  B.1.1.216 - B.1.1.7 -0.072885040 0.011216519 NA -0.094869013 -0.050901066 1.486133e-08
# 4  B.1.1.306 - B.1.1.7 -0.009651664 0.008544000 NA -0.026397595  0.007094268 8.398906e-01
# 5  B.1.1.326 - B.1.1.7 -0.009675804 0.015883056 NA -0.040806022  0.021454415 9.840456e-01
# 6  (B.1.36+) - B.1.1.7 -0.068172930 0.007594587 NA -0.083058046 -0.053287813 1.332268e-15
# 7    B.1.525 - B.1.1.7 -0.085238508 0.027597613 NA -0.139328835 -0.031148180 2.398633e-02
# 8    B.1.351 - B.1.1.7 -0.032660931 0.018019073 NA -0.067977666  0.002655804 4.251242e-01
# 9    B.1.618 - B.1.1.7 -0.012770300 0.013378016 NA -0.038990730  0.013450130 9.103203e-01
# 10 B.1.617.1 - B.1.1.7  0.018251076 0.007153860 NA  0.004229769  0.032272383 1.011091e-01
# 11 B.1.617.2 - B.1.1.7  0.109285248 0.010252277 NA  0.089191155  0.129379342 0.000000e+00
# 12     other - B.1.1.7 -0.009137929 0.007255463 NA -0.023358375  0.005082516 7.736422e-01

# avg growth advantage of B.1.617+ over B.1.1.7 :
max(GISAID_india$date) # 2021-04-30
emtrindia1 = emtrends(fit5_india_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india1 = data.frame(confint(emtrindia1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrindia1$contrasts)$p.value)
delta_r_india1
#                contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1         B.1 - B.1.1.7  0.0263968343 0.006623799 NA  0.01341443  0.039379242 1.143362e-03
# 2       B.1.1 - B.1.1.7 -0.0299932466 0.007848604 NA -0.04537623 -0.014610265 2.067804e-03
# 3   B.1.1.216 - B.1.1.7 -0.0653314095 0.011018056 NA -0.08692640 -0.043736416 2.741571e-07
# 4   B.1.1.306 - B.1.1.7 -0.0001950526 0.008306333 NA -0.01647517  0.016085061 1.000000e+00
# 5   B.1.1.326 - B.1.1.7 -0.0028847947 0.015803607 NA -0.03385930  0.028089707 9.998572e-01
# 6   (B.1.36+) - B.1.1.7 -0.0612521168 0.007325884 NA -0.07561059 -0.046893647 8.643086e-13
# 7     B.1.525 - B.1.1.7 -0.0783396359 0.027472579 NA -0.13218490 -0.024494370 4.410400e-02
# 8     B.1.351 - B.1.1.7 -0.0282647242 0.017531478 NA -0.06262579  0.006096341 5.340290e-01
# 9     B.1.618 - B.1.1.7 -0.0043940890 0.013309714 NA -0.03048065  0.021692471 9.983272e-01
# 10 (B.1.617+) - B.1.1.7  0.0559934737 0.006928861 NA  0.04241316  0.069573791 3.952505e-12
# 11      other - B.1.1.7 -0.0049344983 0.007084030 NA -0.01881894  0.008949946 9.675161e-01


# pairwise growth rate difference (differences in growth rate per day) 
max(GISAID_india$date) # 2021-04-30
emtrindia_pairw = emtrends(fit5_india_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india_pairw = data.frame(confint(emtrindia_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindia_pairw$contrasts)$p.value)
delta_r_india_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1          B.1.1.7 - B.1 -2.261478e-02 0.006960121 NA -0.036256366 -0.0089731918 7.033639e-02
# 2        B.1.1.7 - B.1.1  3.695753e-02 0.008092518 NA  0.021096489  0.0528185771 7.400042e-04
# 3    B.1.1.7 - B.1.1.216  7.288504e-02 0.011216519 NA  0.050901066  0.0948690133 9.594371e-08
# 4    B.1.1.7 - B.1.1.306  9.651664e-03 0.008544000 NA -0.007094268  0.0263975952 9.956231e-01
# 5    B.1.1.7 - B.1.1.326  9.675804e-03 0.015883056 NA -0.021454415  0.0408060223 9.999924e-01
# 6    B.1.1.7 - (B.1.36+)  6.817293e-02 0.007594587 NA  0.053287813  0.0830580462 1.560974e-13
# 7      B.1.1.7 - B.1.525  8.523851e-02 0.027597613 NA  0.031148180  0.1393288355 1.079693e-01
# 8      B.1.1.7 - B.1.351  3.266093e-02 0.018019073 NA -0.002655804  0.0679776664 8.402711e-01
# 9      B.1.1.7 - B.1.618  1.277030e-02 0.013378016 NA -0.013450130  0.0389907300 9.991207e-01
# 10   B.1.1.7 - B.1.617.1 -1.825108e-02 0.007153860 NA -0.032272383 -0.0042297694 3.499180e-01
# 11   B.1.1.7 - B.1.617.2 -1.092852e-01 0.010252277 NA -0.129379342 -0.0891911551 2.009504e-14
# 12       B.1.1.7 - other  9.137929e-03 0.007255463 NA -0.005082516  0.0233583748 9.885479e-01
# 13           B.1 - B.1.1  5.957231e-02 0.005993539 NA  0.047825191  0.0713194331 4.052314e-14
# 14       B.1 - B.1.1.216  9.549982e-02 0.009717435 NA  0.076453995  0.1145456421 4.740652e-14
# 15       B.1 - B.1.1.306  3.226644e-02 0.006577321 NA  0.019375130  0.0451577553 1.806932e-04
# 16       B.1 - B.1.1.326  3.229058e-02 0.014870543 NA  0.003144853  0.0614363121 6.133317e-01
# 17       B.1 - (B.1.36+)  9.078771e-02 0.005179539 NA  0.080635998  0.1009394194 0.000000e+00
# 18         B.1 - B.1.525  1.078533e-01 0.027365251 NA  0.054218380  0.1614881928 7.903753e-03
# 19         B.1 - B.1.351  5.527571e-02 0.017557131 NA  0.020864366  0.0896870547 9.241686e-02
# 20         B.1 - B.1.618  3.538508e-02 0.012384156 NA  0.011112580  0.0596575775 1.885213e-01
# 21       B.1 - B.1.617.1  4.363703e-03 0.005639437 NA -0.006689392  0.0154167968 9.998979e-01
# 22       B.1 - B.1.617.2 -8.667047e-02 0.009416235 NA -0.105125951 -0.0682149880 7.882583e-14
# 23           B.1 - other  3.175271e-02 0.004979081 NA  0.021993889  0.0415115277 1.776350e-07
# 24     B.1.1 - B.1.1.216  3.592751e-02 0.009984197 NA  0.016358840  0.0554961737 2.488690e-02
# 25     B.1.1 - B.1.1.306 -2.730587e-02 0.007278671 NA -0.041571802 -0.0130399366 1.513192e-02
# 26     B.1.1 - B.1.1.326 -2.728173e-02 0.015126640 NA -0.056929399  0.0023659401 8.448072e-01
# 27     B.1.1 - (B.1.36+)  3.121540e-02 0.005851006 NA  0.019747635  0.0426831581 2.712591e-05
# 28       B.1.1 - B.1.525  4.828097e-02 0.027699165 NA -0.006008391  0.1025703401 8.733320e-01
# 29       B.1.1 - B.1.351 -4.296601e-03 0.018056748 NA -0.039687176  0.0310939733 1.000000e+00
# 30       B.1.1 - B.1.618 -2.418723e-02 0.012710398 NA -0.049099155  0.0007246886 7.909908e-01
# 31     B.1.1 - B.1.617.1 -5.520861e-02 0.007053656 NA -0.069033521 -0.0413836972 7.680101e-11
# 32     B.1.1 - B.1.617.2 -1.462428e-01 0.010462014 NA -0.166747952 -0.1257376110 0.000000e+00
# 33         B.1.1 - other -2.781960e-02 0.005926402 NA -0.039435139 -0.0162040683 4.399055e-04
# 34 B.1.1.216 - B.1.1.306 -6.323338e-02 0.010569616 NA -0.083949442 -0.0425173094 1.270861e-06
# 35 B.1.1.216 - B.1.1.326 -6.320924e-02 0.017001796 NA -0.096532143 -0.0298863284 1.691853e-02
# 36 B.1.1.216 - (B.1.36+) -4.712110e-03 0.009554865 NA -0.023439301  0.0140150811 9.999993e-01
# 37   B.1.1.216 - B.1.525  1.235347e-02 0.028824999 NA -0.044142492  0.0688494278 9.999999e-01
# 38   B.1.1.216 - B.1.351 -4.022411e-02 0.019690240 NA -0.078816269 -0.0016319474 7.028754e-01
# 39   B.1.1.216 - B.1.618 -6.011474e-02 0.014596223 NA -0.088722811 -0.0315066687 4.179073e-03
# 40 B.1.1.216 - B.1.617.1 -9.113612e-02 0.010486345 NA -0.111688975 -0.0705832572 6.216139e-13
# 41 B.1.1.216 - B.1.617.2 -1.821703e-01 0.013040688 NA -0.207729566 -0.1566110095 0.000000e+00
# 42     B.1.1.216 - other -6.374711e-02 0.009684404 NA -0.082728194 -0.0447660264 6.214900e-08
# 43 B.1.1.306 - B.1.1.326  2.413997e-05 0.015152485 NA -0.029674186  0.0297224657 1.000000e+00
# 44 B.1.1.306 - (B.1.36+)  5.852127e-02 0.006440033 NA  0.045899033  0.0711434994 1.074696e-13
# 45   B.1.1.306 - B.1.525  7.558684e-02 0.027774462 NA  0.021149898  0.1300237893 2.524338e-01
# 46   B.1.1.306 - B.1.351  2.300927e-02 0.018216483 NA -0.012694384  0.0587129193 9.882622e-01
# 47   B.1.1.306 - B.1.618  3.118636e-03 0.013220494 NA -0.022793057  0.0290303290 1.000000e+00
# 48 B.1.1.306 - B.1.617.1 -2.790274e-02 0.007403035 NA -0.042412422 -0.0133930586 1.426833e-02
# 49 B.1.1.306 - B.1.617.2 -1.189369e-01 0.010696547 NA -0.139901758 -0.0979720659 0.000000e+00
# 50     B.1.1.306 - other -5.137343e-04 0.006481338 NA -0.013216924  0.0121894554 1.000000e+00
# 51 B.1.1.326 - (B.1.36+)  5.849713e-02 0.014692089 NA  0.029701161  0.0872930913 6.855270e-03
# 52   B.1.1.326 - B.1.525  7.556270e-02 0.030779487 NA  0.015236018  0.1358893898 4.124807e-01
# 53   B.1.1.326 - B.1.351  2.298513e-02 0.022565549 NA -0.021242536  0.0672127914 9.983431e-01
# 54   B.1.1.326 - B.1.618  3.094496e-03 0.018789714 NA -0.033732666  0.0399216580 1.000000e+00
# 55 B.1.1.326 - B.1.617.1 -2.792688e-02 0.015353511 NA -0.058019208  0.0021654477 8.370351e-01
# 56 B.1.1.326 - B.1.617.2 -1.189611e-01 0.017186982 NA -0.152646918 -0.0852751863 1.056308e-08
# 57     B.1.1.326 - other -5.378743e-04 0.014843185 NA -0.029629982  0.0285542332 1.000000e+00
# 58   (B.1.36+) - B.1.525  1.706558e-02 0.027560796 NA -0.036952590  0.0710837454 9.999909e-01
# 59   (B.1.36+) - B.1.351 -3.551200e-02 0.017853694 NA -0.070504595 -0.0005194015 7.382121e-01
# 60   (B.1.36+) - B.1.618 -5.540263e-02 0.012478299 NA -0.079859647 -0.0309456126 1.228098e-03
# 61 (B.1.36+) - B.1.617.1 -8.642401e-02 0.006470136 NA -0.099105240 -0.0737427718 0.000000e+00
# 62 (B.1.36+) - B.1.617.2 -1.774582e-01 0.010089034 NA -0.197232321 -0.1576840347 0.000000e+00
# 63     (B.1.36+) - other -5.903500e-02 0.004979561 NA -0.068794761 -0.0492752395 0.000000e+00
# 64     B.1.525 - B.1.351 -5.257758e-02 0.031505828 NA -0.114327865  0.0091727127 9.036850e-01
# 65     B.1.525 - B.1.618 -7.246821e-02 0.029751191 NA -0.130779470 -0.0141569451 4.254537e-01
# 66   B.1.525 - B.1.617.1 -1.034896e-01 0.027385015 NA -0.157163227 -0.0498159406 1.379965e-02
# 67   B.1.525 - B.1.617.2 -1.945238e-01 0.028441074 NA -0.250267236 -0.1387802762 1.628259e-08
# 68       B.1.525 - other -7.610058e-02 0.027465118 NA -0.129931219 -0.0222699369 2.277419e-01
# 69     B.1.351 - B.1.618 -1.989063e-02 0.020917883 NA -0.060888928  0.0211076648 9.991537e-01
# 70   B.1.351 - B.1.617.1 -5.091201e-02 0.017647954 NA -0.085501362 -0.0163226543 1.770902e-01
# 71   B.1.351 - B.1.617.2 -1.419462e-01 0.019216175 NA -0.179609191 -0.1042831688 8.685641e-10
# 72       B.1.351 - other -2.352300e-02 0.017705825 NA -0.058225782  0.0111797777 9.820966e-01
# 73   B.1.618 - B.1.617.1 -3.102138e-02 0.012811729 NA -0.056131904 -0.0059108489 4.353569e-01
# 74   B.1.618 - B.1.617.2 -1.220555e-01 0.014981955 NA -0.151419641 -0.0926914552 1.278488e-11
# 75       B.1.618 - other -3.632370e-03 0.012439960 NA -0.028014244  0.0207495030 1.000000e+00
# 76 B.1.617.1 - B.1.617.2 -9.103417e-02 0.008703469 NA -0.108092657 -0.0739756867 3.530509e-14
# 77     B.1.617.1 - other  2.738901e-02 0.006008334 NA  0.015612888  0.0391651240 7.654373e-04
# 78     B.1.617.2 - other  1.184232e-01 0.009762067 NA  0.099289877  0.1375564781 0.000000e+00



# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-06-01"))
date.to = as.numeric(as.Date("2021-06-01")) # max(GISAID_india$DATE_NUM)+extrapolate

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
           xmax=as.Date("2021-06-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
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
ggarrange(muller_india_raw2+coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-01")))+
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
           xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
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
            coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-05-31")))+
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
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
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
  geom_smooth(data=cases_india_bystate2, aes(x=Date, y=newcases, lwd=I(1.5)), method="loess", span=0.3, se=FALSE, colour="black") +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
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
           xmax=as.Date("2021-06-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
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
