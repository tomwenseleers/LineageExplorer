# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN CASES EXPORTED FROM INDIA (GISAID GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 17 MAY 2021

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




# ANALYSIS VOCs EXPORTED FROM INDIA ####
            
GISAID_indiaexp = GISAID_sel[GISAID_sel$country_exposure=="India"&GISAID_sel$country!="India",]
nrow(GISAID_indiaexp[is.na(GISAID_indiaexp$LINEAGE1),]) # 0 unknown pango clade
GISAID_indiaexp = GISAID_indiaexp[!is.na(GISAID_indiaexp$LINEAGE1),]
nrow(GISAID_indiaexp) # 388

unique(GISAID_indiaexp$division_exposure) # just says "India"
unique(GISAID_indiaexp$country[GISAID_indiaexp$LINEAGE1=="B.1.617+"]) # France Italy Japan Singapore
table(GISAID_indiaexp$country, GISAID_indiaexp$LINEAGE2)

sum(GISAID_indiaexp$LINEAGE1=="B.1.617+") # 150
unique(GISAID_indiaexp$country[GISAID_indiaexp$LINEAGE1=="B.1.1.7"])
sum(GISAID_indiaexp$LINEAGE1=="B.1.1.7") # 83

table(GISAID_indiaexp$LINEAGE1)
table(GISAID_indiaexp$LINEAGE2)

main_lineages = names(table(GISAID_indiaexp$LINEAGE1))[100*table(GISAID_indiaexp$LINEAGE1)/sum(table(GISAID_indiaexp$LINEAGE1)) > 3]
main_lineages
# "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.1.326" "B.1.1.7"   "B.1.36+"   "B.1.617+" 
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526")
main_lineages = union(main_lineages, VOCs)
GISAID_indiaexp$LINEAGE1[!(GISAID_indiaexp$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_indiaexp$LINEAGE2[!(GISAID_indiaexp$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
# remove = names(table(GISAID_indiaexp$LINEAGE1))[table(GISAID_indiaexp$LINEAGE1) < 10]
# GISAID_indiaexp$LINEAGE1[(GISAID_indiaexp$LINEAGE1 %in% remove)] = "other" # minority VOCs
# GISAID_indiaexp$LINEAGE2[(GISAID_indiaexp$LINEAGE2 %in% remove)] = "other" # minority VOCs
table(GISAID_indiaexp$LINEAGE1)
GISAID_indiaexp$LINEAGE1 = factor(GISAID_indiaexp$LINEAGE1)
GISAID_indiaexp$LINEAGE1 = relevel(GISAID_indiaexp$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_indiaexp$LINEAGE1)
# [1] "B.1.1.7"   "B.1.1"     "B.1.1.354" "B.1.351"   "B.1.36+"   "B.1.525"   "B.1.617+"  "other"
levels_LINEAGE1 = c("B.1.1.7","B.1.1","B.1.1.354","B.1.36+",
                    "B.1.525","B.1.351","B.1.617+","other")
GISAID_indiaexp$LINEAGE1 = factor(GISAID_indiaexp$LINEAGE1, levels=levels_LINEAGE1)

GISAID_indiaexp$LINEAGE2 = factor(GISAID_indiaexp$LINEAGE2)
GISAID_indiaexp$LINEAGE2 = relevel(GISAID_indiaexp$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_indiaexp$LINEAGE2)
# "B.1.1.7"   "B.1.1"     "B.1.1.354" "B.1.351"   "B.1.36+"   "B.1.525"   "B.1.617.1" "B.1.617.2" "other" 
levels_LINEAGE2 = c("B.1.1.7","B.1.1","B.1.1.354","B.1.36+",
                    "B.1.525","B.1.351","B.1.617.1","B.1.617.2","other")
GISAID_indiaexp$LINEAGE2 = factor(GISAID_indiaexp$LINEAGE2, levels=levels_LINEAGE2)
firstB16172 = GISAID_indiaexp[GISAID_indiaexp$LINEAGE2=="B.1.617.2",]
firstB16172 = firstB16172[firstB16172$date==min(firstB16172$date),]
firstB16172 # Japan 28/3/2021

# select states of India with a total of > 300 sequences submitted

# GISAID_indiaexp = GISAID_indiaexp[GISAID_indiaexp$division!="India",]
table(GISAID_indiaexp$country)
# France       Italy       Japan New Zealand   Singapore 
# 2          10          70           3         303 

sel_countries = names(table(GISAID_indiaexp$country))[table(GISAID_indiaexp$country) >= 10]
sel_countries # "Italy"     "Japan"     "Singapore"
sel_countries = sel_countries[!sel_countries %in% c("Italy")]
sel_countries # "Japan" "Singapore"

GISAID_indiaexp = GISAID_indiaexp[GISAID_indiaexp$country %in% sel_countries,]
GISAID_indiaexp$country = factor(GISAID_indiaexp$country)
GISAID_indiaexp$country = relevel(GISAID_indiaexp$country, ref="Singapore")
levels_countries = c("Singapore","Japan")
GISAID_indiaexp$country = factor(GISAID_indiaexp$country, levels=levels_countries)
levels(GISAID_indiaexp$country)
table(GISAID_indiaexp$country)
sum(table(GISAID_indiaexp$country)) # 373
table(GISAID_indiaexp$LINEAGE2)

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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots by state_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots by state_raw data.pdf"), width=8, height=6)



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
max(GISAID_indiaexp$date) # 2021-05-03
emtrindiaexp = emtrends(fit1_indiaexp_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp = data.frame(confint(emtrindiaexp, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrindiaexp$contrasts)$p.value)
delta_r_indiaexp
#              contrast      estimate         SE df   asymp.LCL  asymp.UCL   p.value
# 1     B.1.1 - B.1.1.7 -0.05782807 0.008863536 NA -0.07520028 -0.04045586 5.005985e-05
# 2 B.1.1.354 - B.1.1.7 -0.06364057 0.009207200 NA -0.08168635 -0.04559479 2.513379e-05
# 3 (B.1.36+) - B.1.1.7 -0.06009681 0.008934380 NA -0.07760787 -0.04258574 3.486596e-05
# 4   B.1.525 - B.1.1.7  0.04809387 0.051024337 NA -0.05191200  0.14809973 8.553402e-01
# 5   B.1.351 - B.1.1.7  0.06909957 0.036166396 NA -0.00178526  0.13998441 3.292978e-01
# 6 B.1.617.1 - B.1.1.7  0.05072292 0.014920396 NA  0.02147948  0.07996636 2.273220e-02
# 7 B.1.617.2 - B.1.1.7  0.13848473 0.020551992 NA  0.09820357  0.17876590 3.414351e-05
# 8     other - B.1.1.7 -0.05698502 0.008914641 NA -0.07445740 -0.03951265 6.357054e-05

# avg growth advantage of B.1.617+ over B.1.1.7 :
max(GISAID_indiaexp$date) # 2021-04-28
emtrindiaexp1 = emtrends(fit1_indiaexp_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp1 = data.frame(confint(emtrindiaexp1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrindiaexp1$contrasts)$p.value)
delta_r_indiaexp1
#         contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1      B.1.1 - B.1.1.7 -0.05746669 0.008820792 NA -0.074755126 -0.04017826 8.456525e-05
# 2  B.1.1.354 - B.1.1.7 -0.06328159 0.009165857 NA -0.081246344 -0.04531684 4.523097e-05
# 3  (B.1.36+) - B.1.1.7 -0.05973678 0.008891847 NA -0.077164477 -0.04230908 6.084119e-05
# 4    B.1.525 - B.1.1.7  0.04474306 0.048830198 NA -0.050962375  0.14044848 8.459999e-01
# 5    B.1.351 - B.1.1.7  0.06401147 0.034636145 NA -0.003874125  0.13189707 3.408010e-01
# 6 (B.1.617+) - B.1.1.7  0.08753170 0.014393410 NA  0.059321134  0.11574227 1.737915e-04
# 7      other - B.1.1.7 -0.05662413 0.008872198 NA -0.074013315 -0.03923494 1.051542e-04


# pairwise growth rate difference (differences in growth rate per day) 
max(GISAID_indiaexp$date) # 2021-04-28
emtrindiaexp_pairw = emtrends(fit1_indiaexp_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp_pairw = data.frame(confint(emtrindiaexp_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindiaexp_pairw$contrasts)$p.value)
delta_r_indiaexp_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1        B.1.1.7 - B.1.1  0.0578280686 0.008863536 NA  0.040455857  0.0752002802 1.857651e-04
# 2    B.1.1.7 - B.1.1.354  0.0636405743 0.009207200 NA  0.045594794  0.0816863548 9.417837e-05
# 3    B.1.1.7 - (B.1.36+)  0.0600968072 0.008934380 NA  0.042585745  0.0776078698 1.300674e-04
# 4      B.1.1.7 - B.1.525 -0.0480938655 0.051024337 NA -0.148099728  0.0519119970 9.863073e-01
# 5      B.1.1.7 - B.1.351 -0.0690995734 0.036166396 NA -0.139984407  0.0017852601 6.167631e-01
# 6    B.1.1.7 - B.1.617.1 -0.0507229179 0.014920396 NA -0.079966357 -0.0214794787 6.691537e-02
# 7    B.1.1.7 - B.1.617.2 -0.1384847346 0.020551992 NA -0.178765900 -0.0982035695 1.274097e-04
# 8        B.1.1.7 - other  0.0569850232 0.008914641 NA  0.039512649  0.0744573978 2.350247e-04
# 9      B.1.1 - B.1.1.354  0.0058125057 0.003548194 NA -0.001141826  0.0127668378 7.722156e-01
# 10     B.1.1 - (B.1.36+)  0.0022687386 0.002999491 NA -0.003610155  0.0081476322 9.967055e-01
# 11       B.1.1 - B.1.525 -0.1059219341 0.051714259 NA -0.207280019 -0.0045638492 5.353233e-01
# 12       B.1.1 - B.1.351 -0.1269276421 0.037201057 NA -0.199840374 -0.0540149102 6.541814e-02
# 13     B.1.1 - B.1.617.1 -0.1085509866 0.017159227 NA -0.142182454 -0.0749195196 2.646547e-04
# 14     B.1.1 - B.1.617.2 -0.1963128032 0.022379069 NA -0.240174972 -0.1524506341 4.646534e-06
# 15         B.1.1 - other -0.0008430454 0.003224352 NA -0.007162660  0.0054765689 9.999988e-01
# 16 B.1.1.354 - (B.1.36+) -0.0035437671 0.003521208 NA -0.010445208  0.0033576736 9.796301e-01
# 17   B.1.1.354 - B.1.525 -0.1117344398 0.051774835 NA -0.213211252 -0.0102576271 4.721707e-01
# 18   B.1.1.354 - B.1.351 -0.1327401478 0.037285127 NA -0.205817654 -0.0596626412 4.975139e-02
# 19 B.1.1.354 - B.1.617.1 -0.1143634923 0.017340910 NA -0.148351052 -0.0803759327 1.638957e-04
# 20 B.1.1.354 - B.1.617.2 -0.2021253089 0.022518327 NA -0.246260419 -0.1579901988 3.423613e-06
# 21     B.1.1.354 - other -0.0066555511 0.003757127 NA -0.014019384  0.0007082816 6.984242e-01
# 22   (B.1.36+) - B.1.525 -0.1081906726 0.051726843 NA -0.209573421 -0.0068079242 5.100984e-01
# 23   (B.1.36+) - B.1.351 -0.1291963806 0.037218493 NA -0.202143286 -0.0562494751 5.865683e-02
# 24 (B.1.36+) - B.1.617.1 -0.1108197251 0.017197100 NA -0.144525421 -0.0771140291 2.142475e-04
# 25 (B.1.36+) - B.1.617.2 -0.1985815418 0.022407897 NA -0.242500214 -0.1546628697 4.058376e-06
# 26     (B.1.36+) - other -0.0031117840 0.003236082 NA -0.009454387  0.0032308193 9.845219e-01
# 27     B.1.525 - B.1.351 -0.0210057080 0.061069296 NA -0.140699329  0.0986879127 9.999899e-01
# 28   B.1.525 - B.1.617.1 -0.0026290525 0.051617595 NA -0.103797680  0.0985395751 1.000000e+00
# 29   B.1.525 - B.1.617.2 -0.0903908691 0.052963138 NA -0.194196712  0.0134149736 7.350598e-01
# 30       B.1.525 - other  0.1050788887 0.051722814 NA  0.003704036  0.2064537413 5.450829e-01
# 31   B.1.351 - B.1.617.1  0.0183766555 0.036656277 NA -0.053468327  0.0902216381 9.998239e-01
# 32   B.1.351 - B.1.617.2 -0.0693851612 0.038112880 NA -0.144085034  0.0053147113 6.699569e-01
# 33       B.1.351 - other  0.1260845967 0.037212976 NA  0.053148503  0.1990206901 6.832116e-02
# 34 B.1.617.1 - B.1.617.2 -0.0877618167 0.020406806 NA -0.127758422 -0.0477652118 1.211260e-02
# 35     B.1.617.1 - other  0.1077079412 0.017185000 NA  0.074025960  0.1413899224 2.940802e-04
# 36     B.1.617.2 - other  0.1954697578 0.022398956 NA  0.151568610  0.2393709055 4.976883e-06



# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-06-01"))
date.to = as.numeric(as.Date("2021-06-01")) # max(GISAID_indiaexp$DATE_NUM)+extrapolate

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
           xmax=as.Date("2021-06-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2\nAMONG CASES EXPORTED FROM INDIA (GISAID data Singapore & Japan)")
muller_indiaexp_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_muller plots_multinom fit.pdf"), width=8, height=6)


library(ggpubr)
ggarrange(muller_indiaexp_raw2 + coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-01")))+
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
  coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-01")), ylim=c(0.005, 0.95), expand=c(0,0))
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
plot_indiaexp_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_multinom fit_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india exported_multinom fit_response scale.pdf"), width=8, height=6)
