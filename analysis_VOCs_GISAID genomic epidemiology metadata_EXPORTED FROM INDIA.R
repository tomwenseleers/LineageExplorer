# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN CASES EXPORTED FROM INDIA (GISAID GENOMIC EPIDEMIOLOGY METADATA)
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
#              Country  Lineage Count
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




# ANALYSIS VOCs EXPORTED FROM INDIA ####
            
GISAID_indiaexp = GISAID_sel[GISAID_sel$country_exposure=="India"&GISAID_sel$country!="India",]
nrow(GISAID_indiaexp[is.na(GISAID_indiaexp$LINEAGE1),]) # 0 unknown pango clade
GISAID_indiaexp = GISAID_indiaexp[!is.na(GISAID_indiaexp$LINEAGE1),]
nrow(GISAID_indiaexp) # 443

unique(GISAID_indiaexp$division_exposure) # just says "India"
unique(GISAID_indiaexp$country[GISAID_indiaexp$LINEAGE1=="B.1.617+"]) # France Italy Japan Singapore Spain
table(GISAID_indiaexp$country, GISAID_indiaexp$LINEAGE2)

sum(GISAID_indiaexp$LINEAGE1=="B.1.617+") # 197
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
# B.1.1 B.1.1.354   B.1.1.7   B.1.351   B.1.36+  B.1.617+     other 
# 45        25        86        11        43       197        36
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
firstB16172 = GISAID_indiaexp[GISAID_indiaexp$LINEAGE2=="B.1.617.2",]
firstB16172 = firstB16172[firstB16172$date==min(firstB16172$date),]
firstB16172 # Japan 28/3/2021

# select countries with >100 imported sequences

# GISAID_indiaexp = GISAID_indiaexp[GISAID_indiaexp$division!="India",]
table(GISAID_indiaexp$country)
# France       Italy       Japan New Zealand   Singapore       Spain 
# 2          11         123           3         303           1 

sel_countries = names(table(GISAID_indiaexp$country))[table(GISAID_indiaexp$country) >= 100]
sel_countries # "Japan"     "Singapore"

GISAID_indiaexp = GISAID_indiaexp[GISAID_indiaexp$country %in% sel_countries,]
GISAID_indiaexp$country = factor(GISAID_indiaexp$country)
GISAID_indiaexp$country = relevel(GISAID_indiaexp$country, ref="Singapore")
levels_countries = c("Singapore","Japan")
GISAID_indiaexp$country = factor(GISAID_indiaexp$country, levels=levels_countries)
levels(GISAID_indiaexp$country)
table(GISAID_indiaexp$country)
sum(table(GISAID_indiaexp$country)) # 426
table(GISAID_indiaexp$LINEAGE2)
range(GISAID_indiaexp$date) # "2020-06-17" "2021-05-10"

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
max(GISAID_indiaexp$date) # 2021-05-03
emtrindiaexp = emtrends(fit1_indiaexp_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp = data.frame(confint(emtrindiaexp, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrindiaexp$contrasts)$p.value)
delta_r_indiaexp
#              contrast      estimate         SE df   asymp.LCL  asymp.UCL   p.value
# 1     B.1.1 - B.1.1.7 -0.04871833 0.007221062 NA -0.0628713474 -0.03456531 5.811651e-05
# 2 B.1.1.354 - B.1.1.7 -0.05474720 0.007618034 NA -0.0696782731 -0.03981613 2.907219e-05
# 3 (B.1.36+) - B.1.1.7 -0.05142842 0.007314199 NA -0.0657639893 -0.03709285 3.701764e-05
# 4   B.1.351 - B.1.1.7  0.06078333 0.031292142 NA -0.0005481407  0.12211480 2.982309e-01
# 5 B.1.617.1 - B.1.1.7  0.03786766 0.012865118 NA  0.0126524884  0.06308282 5.533983e-02
# 6 B.1.617.2 - B.1.1.7  0.14158751 0.017865548 NA  0.1065716778  0.17660334 9.594031e-06
# 7     other - B.1.1.7 -0.04584409 0.007181736 NA -0.0599200367 -0.03176815 1.049410e-04

# avg growth advantage of B.1.617+ over B.1.1.7 :
max(GISAID_indiaexp$date) # 2021-05-10
emtrindiaexp1 = emtrends(fit1_indiaexp_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp1 = data.frame(confint(emtrindiaexp1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrindiaexp1$contrasts)$p.value)
delta_r_indiaexp1
#         contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1      B.1.1 - B.1.1.7 -0.04802674 0.007126377 NA -0.061994188 -0.03405930 1.091042e-04
# 2  B.1.1.354 - B.1.1.7 -0.05405761 0.007527852 NA -0.068811925 -0.03930329 5.891306e-05
# 3  (B.1.36+) - B.1.1.7 -0.05073738 0.007220470 NA -0.064889246 -0.03658552 7.284621e-05
# 4    B.1.351 - B.1.1.7  0.05183113 0.028334812 NA -0.003704084  0.10736634 3.295470e-01
# 5 (B.1.617+) - B.1.1.7  0.08631062 0.012665660 NA  0.061486379  0.11113485 9.806339e-05
# 6      other - B.1.1.7 -0.04514971 0.007086709 NA -0.059039400 -0.03126001 1.858295e-04


# pairwise growth rate difference (differences in growth rate per day) 
max(GISAID_indiaexp$date) # 2021-05-10
emtrindiaexp_pairw = emtrends(fit1_indiaexp_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_indiaexp$DATE_NUM)))
delta_r_indiaexp_pairw = data.frame(confint(emtrindiaexp_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindiaexp_pairw$contrasts)$p.value)
delta_r_indiaexp_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1        B.1.1.7 - B.1.1  0.048718327 0.007221062 NA  0.0345653063  0.0628713474 1.901069e-04
# 2    B.1.1.7 - B.1.1.354  0.054747201 0.007618034 NA  0.0398161279  0.0696782731 9.594132e-05
# 3    B.1.1.7 - (B.1.36+)  0.051428422 0.007314199 NA  0.0370928549  0.0657639893 1.218055e-04
# 4      B.1.1.7 - B.1.351 -0.060783331 0.031292142 NA -0.1221148019  0.0005481407 5.470466e-01
# 5    B.1.1.7 - B.1.617.1 -0.037867656 0.012865118 NA -0.0630828244 -0.0126524884 1.360129e-01
# 6    B.1.1.7 - B.1.617.2 -0.141587508 0.017865548 NA -0.1766033378 -0.1065716778 3.203204e-05
# 7        B.1.1.7 - other  0.045844093 0.007181736 NA  0.0317681500  0.0599200367 3.403049e-04
# 8      B.1.1 - B.1.1.354  0.006028874 0.003420229 NA -0.0006746521  0.0127323995 6.516053e-01
# 9      B.1.1 - (B.1.36+)  0.002710095 0.002862429 NA -0.0029001633  0.0083203539 9.751905e-01
# 10       B.1.1 - B.1.351 -0.109501657 0.031993581 NA -0.1722079244 -0.0467953905 6.024620e-02
# 11     B.1.1 - B.1.617.1 -0.086585983 0.014287093 NA -0.1145881718 -0.0585837946 5.786184e-04
# 12     B.1.1 - B.1.617.2 -0.190305835 0.019129139 NA -0.2277982577 -0.1528134115 2.154288e-06
# 13         B.1.1 - other -0.002874233 0.002988413 NA -0.0087314151  0.0029829482 9.730204e-01
# 14 B.1.1.354 - (B.1.36+) -0.003318778 0.003411074 NA -0.0100043609  0.0033668041 9.713235e-01
# 15   B.1.1.354 - B.1.351 -0.115530531 0.032086474 NA -0.1784188636 -0.0526421986 4.402389e-02
# 16 B.1.1.354 - B.1.617.1 -0.092614857 0.014494189 NA -0.1210229463 -0.0642067676 3.368058e-04
# 17 B.1.1.354 - B.1.617.2 -0.196334708 0.019283806 NA -0.2341302731 -0.1585391435 1.620806e-06
# 18     B.1.1.354 - other -0.008903107 0.003596377 NA -0.0159518769 -0.0018543375 2.799294e-01
# 19   (B.1.36+) - B.1.351 -0.112211753 0.032015511 NA -0.1749610016 -0.0494625039 5.213882e-02
# 20 (B.1.36+) - B.1.617.1 -0.089296079 0.014336293 NA -0.1173946970 -0.0611974600 4.381853e-04
# 21 (B.1.36+) - B.1.617.2 -0.193015930 0.019165576 NA -0.2305797682 -0.1554520915 1.853640e-06
# 22     (B.1.36+) - other -0.005584329 0.003045870 NA -0.0115541247  0.0003854672 6.104306e-01
# 23   B.1.351 - B.1.617.1  0.022915674 0.031762292 NA -0.0393372732  0.0851686217 9.947071e-01
# 24   B.1.351 - B.1.617.2 -0.080804177 0.032533740 NA -0.1445691354 -0.0170392189 2.766854e-01
# 25       B.1.351 - other  0.106627424 0.031982916 NA  0.0439420609  0.1693127870 7.031795e-02
# 26 B.1.617.1 - B.1.617.2 -0.103719851 0.017765102 NA -0.1385388105 -0.0689008922 8.392346e-04
# 27     B.1.617.1 - other  0.083711750 0.014263003 NA  0.0557567776  0.1116667219 7.968469e-04
# 28     B.1.617.2 - other  0.187431601 0.019111738 NA  0.1499732820  0.2248899202 2.566746e-06



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

