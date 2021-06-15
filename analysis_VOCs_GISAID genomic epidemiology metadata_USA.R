# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE USA (GISAID GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 12 JUNE 2021

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
today = as.Date("2021-06-11")
today_num = as.numeric(today)
today # "2021-06-11"
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata (file version metadata_2021-06-10_09-46.tsv.gz)
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-06-10_09-46.tsv.gz"), col_types = cols(.default = "c")) 
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
GISAID = GISAID[GISAID$date>=as.Date("2020-01-01"),]
range(GISAID$date) # "2020-01-01" "2021-06-06"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 7 sept 63r old male from Madhya Pradesh

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
# 1      2791     17357        81 

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
# 169500      Australia B.1.617+   219
# 169504        Bahrain B.1.617+    22
# 169505     Bangladesh B.1.617+    44
# 169508        Belgium B.1.617+   206
# 169516         Brazil B.1.617+    14
# 169522         Canada B.1.617+    52
# 169533 Czech Republic B.1.617+    15
# 169535        Denmark B.1.617+   121
# 169545         France B.1.617+   121
# 169550        Germany B.1.617+   634
# 169564          India B.1.617+  6250
# 169565      Indonesia B.1.617+    34
# 169566           Iran B.1.617+    11
# 169568        Ireland B.1.617+   236
# 169569         Israel B.1.617+    36
# 169570          Italy B.1.617+   150
# 169572          Japan B.1.617+   167
# 169583     Luxembourg B.1.617+    58
# 169586       Malaysia B.1.617+    12
# 169590         Mexico B.1.617+    37
# 169599    Netherlands B.1.617+    62
# 169600    New Zealand B.1.617+    17
# 169604         Norway B.1.617+    56
# 169613         Poland B.1.617+    64
# 169614       Portugal B.1.617+    74
# 169615          Qatar B.1.617+    23
# 169617        Romania B.1.617+    19
# 169618         Russia B.1.617+   168
# 169629      Singapore B.1.617+   305
# 169634   South Africa B.1.617+    21
# 169635    South Korea B.1.617+    32
# 169637          Spain B.1.617+   172
# 169640         Sweden B.1.617+    42
# 169641    Switzerland B.1.617+    95
# 169643       Thailand B.1.617+    93
# 169653 United Kingdom B.1.617+ 22792
# 169655            USA B.1.617+  2198
# 169658        Vietnam B.1.617+    54

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>100,]$Country))
sel_countries_target
# [1] "Australia"      "Belgium"        "Denmark"        "Germany"        "India"          "Ireland"        "Italy"          "Japan"          "Singapore"     
# [10] "United Kingdom" "USA"     

sel_ref_lineage = "B.1.1.7"

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Belgium"        "Denmark"        "Germany"        "India"          "Ireland"        "Italy"          "Japan"          "Singapore"     
# [10] "United Kingdom" "USA" 

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Belgium"        "Denmark"        "Germany"        "India"          "Ireland"        "Italy"          "Japan"          "Singapore"     
# [10] "United Kingdom" "USA"

# sel_countries = sel_countries[!(sel_countries %in% c("Japan","USA"))] # Japan is almost only import & for USA we do separate analysis by state

# ANALYSIS FOR THE USA
sel_countries = "USA"

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






# ANALYSIS OF VOCs ACROSS DIFFERENT STATES OF THE USA ####

GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_sel) # 1316040
unique(GISAID_sel$country)

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))
            
# GISAID_sel = GISAID_sel[GISAID_sel$country_exposure=="India"&GISAID_sel$country!="India",]
# nrow(GISAID_sel[is.na(GISAID_sel$LINEAGE1),]) # 0 unknown pango clade
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$LINEAGE1),]
nrow(GISAID_sel) # 516488

GISAID_sel = GISAID_sel[GISAID_sel$country==GISAID_sel$country_exposure,] # we remove travel-related cases
nrow(GISAID_sel) # 516470

sum(GISAID_sel$LINEAGE1=="B.1.617+") # 1656
unique(GISAID_sel$country[GISAID_sel$LINEAGE1=="B.1.1.7"])
sum(GISAID_sel$LINEAGE1=="B.1.1.7") # 153268
sum(GISAID_sel$LINEAGE1=="B.1.1.519") # 153268

table(GISAID_sel$LINEAGE1)
table(GISAID_sel$LINEAGE2)

main_lineages = names(table(GISAID_sel$LINEAGE1))[100*table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) > 3]
main_lineages
# "B.1.1.7" "B.1"     "B.1.2"   "B.1.429" "B.1.526" "other"
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526","B.1.1.519")
main_lineages = union(main_lineages, VOCs)
GISAID_sel$LINEAGE1[!(GISAID_sel$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_sel$LINEAGE2[!(GISAID_sel$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove1 = names(table(GISAID_sel$LINEAGE1))[table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) < 0.01]
remove1 = remove1[!(remove1 %in% c("B.1.351","B.1.1.7","P.1","B.1.617+","B.1.1.519"))]
remove2 = names(table(GISAID_sel$LINEAGE2))[table(GISAID_sel$LINEAGE2)/sum(table(GISAID_sel$LINEAGE2)) < 0.01]
remove2 = remove2[!(remove2 %in% c("B.1.351","B.1.1.7","P.1","B.1.617.2","B.1.617.1","B.1.1.519"))]
GISAID_sel$LINEAGE1[(GISAID_sel$LINEAGE1 %in% remove1)] = "other" # minority VOCs
GISAID_sel$LINEAGE2[(GISAID_sel$LINEAGE2 %in% remove2)] = "other" # minority VOCs
table(GISAID_sel$LINEAGE1)
#     B.1 B.1.1.519   B.1.1.7     B.1.2   B.1.351   B.1.429   B.1.526  B.1.617+     other       P.1 
#   36382     11407    153268     84349      1724     30058     19870      1656    165778     11978 
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1)
GISAID_sel$LINEAGE1 = relevel(GISAID_sel$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE1)
# "B.1.1.7"   "B.1"       "B.1.2"     "B.1.1.519" "B.1.351"   "B.1.429"   "B.1.526"   "P.1"       "B.1.617+"  "other"  
levels_LINEAGE1 = c("B.1.1.7","B.1","B.1.2","B.1.1.519",
                    "B.1.351","B.1.429","B.1.526",
                    "P.1","B.1.617+","other")
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1, levels=levels_LINEAGE1)

GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = relevel(GISAID_sel$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE2)
# "B.1.1.7"   "B.1"       "B.1.1.519" "B.1.2"     "B.1.351"   "B.1.429"   "B.1.526"   "B.1.617.1" "B.1.617.2" "other"     "P.1"   
levels_LINEAGE2 = c("B.1.1.7","B.1","B.1.2","B.1.1.519",
                    "B.1.351","B.1.429","B.1.526",
                    "P.1","B.1.617.1","B.1.617.2","other")
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2, levels=levels_LINEAGE2)

# GISAID_sel = GISAID_sel[GISAID_sel$division!="India",]
table(GISAID_sel$country)

GISAID_sel$state  = GISAID_sel$division
GISAID_sel = GISAID_sel[!grepl("Princess",GISAID_sel$state),] # remove cruise boat data
GISAID_sel = GISAID_sel[!grepl("Islands|Guam",GISAID_sel$state),] # remove Northern Mariana Islands & Virgin Islands & Guam
GISAID_sel = GISAID_sel[!grepl("USA",GISAID_sel$state),] # remove data with unspecified state
GISAID_sel$state[grepl("Washington",GISAID_sel$state)] = "Washington" # Washtington DC -> Washington
sel_states = c("Arizona", "California", "Colorado", "Florida", "Illinois", "Indiana", "Kansas", 
               "Maryland", "Massachusetts", "Minnesota", "Missouri", "Nebraska", "New Jersey", "New York", 
               "Oregon", "Texas", "Utah", "Virginia", "Washington", "Wisconsin") # 20 states with most data for B.1.617.2
GISAID_sel = GISAID_sel[GISAID_sel$state %in% sel_states,]

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
GISAID_sel = GISAID_sel[-which(grepl("B.1.617", GISAID_sel$pango_lineage, fixed=TRUE)&GISAID_sel$date<=as.Date("2021-04-14")),]

GISAID_sel$state = factor(GISAID_sel$state)
GISAID_sel$state = relevel(GISAID_sel$state, ref="Alaska")
levels(GISAID_sel$state)

levels_states = levels(GISAID_sel$state)
GISAID_sel$state = factor(GISAID_sel$state, levels=levels_states)
levels(GISAID_sel$state)
table(GISAID_sel$state)
sum(table(GISAID_sel$state)) # 
table(GISAID_sel$LINEAGE2)
# B.1.1.7       B.1     B.1.2   B.1.351   B.1.429   B.1.526       P.1 B.1.617.1 B.1.617.2     other 
#  153268     36382     84349      1724     30058     19870     11978       235      1420    177186

range(GISAID_sel$date) # "2020-01-03" "2021-06-04"

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek1 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$LINEAGE1))
colnames(data_agbyweek1) = c("floor_date", "LINEAGE1", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_sel) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$LINEAGE1 = factor(data_agbyweek1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

data_agbyweek2 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_sel) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL


# aggregated by week and state for selected variant lineages
data_agbyweekregion1 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$state, GISAID_sel$LINEAGE1))
colnames(data_agbyweekregion1) = c("floor_date", "division", "LINEAGE1", "count")
data_agbyweekregion1_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$floor_date,data_agbyweekregion1$division), 
                                                                  interaction(data_agbyweekregion1_sum$floor_date,data_agbyweekregion1_sum$division))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_sel) # correct
data_agbyweekregion1$collection_date = as.Date(as.character(data_agbyweekregion1$floor_date))
data_agbyweekregion1$LINEAGE1 = factor(data_agbyweekregion1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweekregion1$division = factor(data_agbyweekregion1$division, levels=levels_states)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]
data_agbyweekregion1$floor_date = NULL

data_agbyweekregion2 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$state, GISAID_sel$LINEAGE2))
colnames(data_agbyweekregion2) = c("floor_date", "division", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$division), 
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$division))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_sel) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$division = factor(data_agbyweekregion2$division, levels=levels_states)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]
data_agbyweekregion2$floor_date = NULL


# MULLER PLOT (RAW DATA)
library(scales)
n1 = length(levels(GISAID_sel$LINEAGE1))
lineage_cols1 = hcl(h = seq(15, 320, length = n1), l = 65, c = 200)
lineage_cols1[which(levels(GISAID_sel$LINEAGE1)=="B.1.617+")] = "magenta"
lineage_cols1[which(levels(GISAID_sel$LINEAGE1)=="other")] = "grey75"

n2 = length(levels(GISAID_sel$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.617.2")] = "magenta"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="other")] = "grey75"

muller_sel_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE1)) + 
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_sel_raw1

muller_sel_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_sel_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots_raw data.pdf"), width=8, height=6)


muller_usabystate_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ division) +
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_usabystate_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots by state_raw data.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots by state_raw data.pdf"), width=10, height=6)



# multinomial fits
data_agbyweekregion2$LINEAGE2 = relevel(data_agbyweekregion2$LINEAGE2, ref="B.1.1.7")
data_agbyweekregion2$DATE_NUM = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$STATE = data_agbyweekregion2$division

library(nnet)
library(splines)
set.seed(1)
fit1_usa_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)+STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
# fit2_usa_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)*STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
fit3_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)+STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
# fit4_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)*STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
fit5_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3)+STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
# fit6_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3)*STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
BIC(fit1_usa_multi, fit3_usa_multi, fit5_usa_multi) 
# fit5_usa_multi fits best (lowest BIC), but I will use the slightly simpler fit3_usa_multi

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrusa = emtrends(fit3_usa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_usa = data.frame(confint(emtrusa, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrusa$contrasts)$p.value)
delta_r_usa
#                         contrast      estimate         SE df   asymp.LCL  asymp.UCL   p.value
# 1        B.1 - B.1.1.7 -0.044903415 0.0004486334 NA -0.045782721 -0.044024110 4.551914e-15
# 2      B.1.2 - B.1.1.7 -0.062716997 0.0004330041 NA -0.063565669 -0.061868325 4.551914e-15
# 3  B.1.1.519 - B.1.1.7 -0.051655999 0.0010049880 NA -0.053625739 -0.049686258 4.551914e-15
# 4    B.1.351 - B.1.1.7 -0.008783653 0.0017264225 NA -0.012167379 -0.005399927 7.666043e-06
# 5    B.1.429 - B.1.1.7 -0.055160278 0.0005827699 NA -0.056302486 -0.054018070 4.551914e-15
# 6    B.1.526 - B.1.1.7 -0.015917926 0.0007399809 NA -0.017368262 -0.014467590 4.551914e-15
# 7        P.1 - B.1.1.7  0.007315964 0.0011573052 NA  0.005047687  0.009584240 1.419977e-08
# 8  B.1.617.1 - B.1.1.7  0.041004528 0.0051546644 NA  0.030901572  0.051107485 9.845458e-13
# 9  B.1.617.2 - B.1.1.7  0.080220388 0.0019608439 NA  0.076377205  0.084063572 4.551914e-15
# 10     other - B.1.1.7 -0.042480694 0.0003901294 NA -0.043245334 -0.041716054 4.551914e-15


# fitted prop of different LINEAGES in the USA today
# 76% [74%-78%] now estimated to be B.1.617.2 across all regions
multinom_preds_today_avg = data.frame(emmeans(fit3_usa_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1    B.1.1.7 0.5675888177 1.191686e-02 NA 0.544232193 0.5909454428
# 2        B.1 0.0007500372 3.323333e-05 NA 0.000684901 0.0008151733
# 3      B.1.2 0.0012017610 4.577302e-05 NA 0.001112048 0.0012914745
# 4  B.1.1.519 0.0018026213 1.208983e-04 NA 0.001565665 0.0020395777
# 5    B.1.351 0.0059424142 5.928059e-04 NA 0.004780536 0.0071042923
# 6    B.1.429 0.0021943493 9.992451e-05 NA 0.001998501 0.0023901977
# 7    B.1.526 0.0240660663 1.065506e-03 NA 0.021977713 0.0261544194
# 8        P.1 0.0935439279 4.100515e-03 NA 0.085507067 0.1015807892
# 9  B.1.617.1 0.0044971781 1.105183e-03 NA 0.002331059 0.0066632975
# 10 B.1.617.2 0.2865557253 1.416023e-02 NA 0.258802185 0.3143092654
# 11     other 0.0118571017 3.807589e-04 NA 0.011110828 0.0126033755

# 43% [39%-47%] non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.4324112 0.3898705 0.4749519 


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-11-01"))
date.to = as.numeric(as.Date("2021-07-31")) # max(GISAID_sel$DATE_NUM)+extrapolate

# multinomial model predictions (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to), STATE=levels_states))
fit_usa_multi_preds_bystate = data.frame(predgrid, as.data.frame(predict(fit3_usa_multi, newdata=predgrid, type="prob")))
fit_usa_multi_preds_bystate = gather(fit_usa_multi_preds_bystate, LINEAGE2, prob, all_of(levels_LINEAGE2), factor_key=TRUE)
fit_usa_multi_preds_bystate$collection_date = as.Date(fit_usa_multi_preds_bystate$DATE_NUM, origin="1970-01-01")
fit_usa_multi_preds_bystate$LINEAGE2 = factor(fit_usa_multi_preds_bystate$LINEAGE2, levels=levels_LINEAGE2) 

muller_usa_mfit_bystate = ggplot(data=fit_usa_multi_preds_bystate, 
                                   aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ STATE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data, multinomial fit)")
muller_usa_mfit_bystate

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots_multinom fit by state.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots_multinom fit by state.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_usabystate_raw2 + coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_usa_mfit_bystate+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots multipanel_multinom fit by state.png"), width=10, height=10)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots multipanel_multinom fit by state.pdf"), width=10, height=10)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions by state with confidence intervals (but slower)
fit_usa_multi_preds_bystate_withCI = data.frame(emmeans(fit3_usa_multi,
                                                        ~ LINEAGE2,
                                                        by=c("DATE_NUM","STATE"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                        mode="prob", df=NA))
fit_usa_multi_preds_bystate_withCI$collection_date = as.Date(fit_usa_multi_preds_bystate_withCI$DATE_NUM, origin="1970-01-01")
fit_usa_multi_preds_bystate_withCI$LINEAGE2 = factor(fit_usa_multi_preds_bystate_withCI$LINEAGE2, levels=levels_LINEAGE2)
fit_usa_multi_preds2 = fit_usa_multi_preds_bystate_withCI

fit_usa_multi_preds_bystate_withCI[fit_usa_multi_preds_bystate_withCI$collection_date==as.Date("2021-06-13")&fit_usa_multi_preds_bystate_withCI$LINEAGE2=="B.1.617.2",]
#        LINEAGE2 DATE_NUM         STATE       prob          SE df  asymp.LCL  asymp.UCL collection_date
# 362  B.1.617.2    18791       Arizona 0.2425780 0.04367862 NA 0.15696950 0.3281865      2021-06-13
# 791  B.1.617.2    18791    California 0.4198529 0.02409480 NA 0.37262801 0.4670779      2021-06-13
# 1220 B.1.617.2    18791      Colorado 0.4302180 0.02681023 NA 0.37767095 0.4827651      2021-06-13
# 1649 B.1.617.2    18791       Florida 0.1935387 0.02109506 NA 0.15219310 0.2348842      2021-06-13
# 2078 B.1.617.2    18791      Illinois 0.2177921 0.02336059 NA 0.17200621 0.2635780      2021-06-13
# 2507 B.1.617.2    18791       Indiana 0.2896677 0.03429757 NA 0.22244566 0.3568897      2021-06-13
# 2936 B.1.617.2    18791        Kansas 0.3933147 0.03587644 NA 0.32299815 0.4636312      2021-06-13
# 3365 B.1.617.2    18791      Maryland 0.2455503 0.03563444 NA 0.17570808 0.3153925      2021-06-13
# 3794 B.1.617.2    18791 Massachusetts 0.3283536 0.02410497 NA 0.28110868 0.3755984      2021-06-13
# 4223 B.1.617.2    18791     Minnesota 0.1185615 0.01774051 NA 0.08379073 0.1533322      2021-06-13
# 4652 B.1.617.2    18791      Missouri 0.5313979 0.04471123 NA 0.44376549 0.6190303      2021-06-13
# 5081 B.1.617.2    18791      Nebraska 0.1823259 0.05796164 NA 0.06872313 0.2959286      2021-06-13
# 5510 B.1.617.2    18791    New Jersey 0.4502863 0.03413220 NA 0.38338838 0.5171842      2021-06-13
# 5939 B.1.617.2    18791      New York 0.3216130 0.02618866 NA 0.27028420 0.3729419      2021-06-13
# 6368 B.1.617.2    18791        Oregon 0.1149993 0.03685883 NA 0.04275734 0.1872413      2021-06-13
# 6797 B.1.617.2    18791         Texas 0.3397006 0.02477120 NA 0.29114989 0.3882512      2021-06-13
# 7226 B.1.617.2    18791          Utah 0.6378122 0.02987422 NA 0.57925977 0.6963646      2021-06-13
# 7655 B.1.617.2    18791      Virginia 0.3498320 0.04350665 NA 0.26456054 0.4351035      2021-06-13
# 8084 B.1.617.2    18791    Washington 0.3283849 0.02540988 NA 0.27858249 0.3781874      2021-06-13
# 8513 B.1.617.2    18791     Wisconsin 0.2237085 0.04136560 NA 0.14263341 0.3047836      2021-06-13

# fit_usa_multi_preds2 = fit_usa_multi_preds_bystate # without CIs
# fit_usa_multi_preds2$asymp.LCL = NA
# fit_usa_multi_preds2$asymp.UCL = NA


# on logit scale:

ymin = 0.001
ymax = 0.999
fit_usa_multi_preds2$asymp.LCL[fit_usa_multi_preds2$asymp.LCL<ymin] = ymin
fit_usa_multi_preds2$asymp.UCL[fit_usa_multi_preds2$asymp.UCL<ymin] = ymin
fit_usa_multi_preds2$asymp.UCL[fit_usa_multi_preds2$asymp.UCL>ymax] = ymax
fit_usa_multi_preds2$prob[fit_usa_multi_preds2$prob<ymin] = ymin

plot_usa_mfit_logit = qplot(data=fit_usa_multi_preds2, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~ STATE) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
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
                        range=c(0.5/2, 3/2), limits=c(1,max(data_agbyweekregion2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_usa_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_logit scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_usa_mfit = qplot(data=fit_usa_multi_preds2, x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~ STATE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-11-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweekregion2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5/2, 3/2), limits=c(1,max(data_agbyweekregion2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_usa_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_response scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_response scale.pdf"), width=10, height=6)




# overall multinomial model predictions with confidence intervals
fit_usa_multi_preds_withCI = data.frame(emmeans(fit3_usa_multi,
                                                        ~ LINEAGE2,
                                                        by=c("DATE_NUM"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                        mode="prob", df=NA))
fit_usa_multi_preds_withCI$collection_date = as.Date(fit_usa_multi_preds_withCI$DATE_NUM, origin="1970-01-01")
fit_usa_multi_preds_withCI$LINEAGE2 = factor(fit_usa_multi_preds_withCI$LINEAGE2, levels=levels_LINEAGE2)
fit_usa_multi_preds3 = fit_usa_multi_preds_withCI

fit_usa_multi_preds_withCI[fit_usa_multi_preds_withCI$collection_date==as.Date("2021-06-13")&fit_usa_multi_preds_withCI$LINEAGE2=="B.1.617.2",]
#        LINEAGE2 DATE_NUM         STATE       prob          SE df  asymp.LCL  asymp.UCL collection_date
# 362 B.1.617.2    18791 0.3179744 0.01571049 NA 0.2871824 0.3487664      2021-06-13

# fit_usa_multi_preds2 = fit_usa_multi_preds # without CIs
# fit_usa_multi_preds2$asymp.LCL = NA
# fit_usa_multi_preds2$asymp.UCL = NA


# on logit scale:

ymin = 0.001
ymax = 0.999
fit_usa_multi_preds3$asymp.LCL[fit_usa_multi_preds3$asymp.LCL<ymin] = ymin
fit_usa_multi_preds3$asymp.UCL[fit_usa_multi_preds3$asymp.UCL<ymin] = ymin
fit_usa_multi_preds3$asymp.UCL[fit_usa_multi_preds3$asymp.UCL>ymax] = ymax
fit_usa_multi_preds3$prob[fit_usa_multi_preds3$prob<ymin] = ymin

plot_usa_avg_mfit_logit = qplot(data=fit_usa_multi_preds3, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data 20 states, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
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
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_usa_avg_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_logit scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_usa_avg_mfit = qplot(data=fit_usa_multi_preds3, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data 20 states, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-11-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_usa_avg_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_response scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_response scale.pdf"), width=10, height=6)





# plot new cases by variant & state ####
# TO DO: still need to finish this part

us_cases_by_state = read.csv("https://github.com/nytimes/covid-19-data/raw/master/us-states.csv")
us_data_by_state$date = as.Date(us_data_by_state$date)
us_data_by_state$state = factor(us_data_by_state$state, 
                                levels=c("Washington","Illinois","California",
                                         "Arizona","Massachusetts","Wisconsin",
                                         "Texas","Nebraska","Utah","Oregon",
                                         "Florida","New York","Rhode Island",
                                         "Georgia","New Hampshire","North Carolina",
                                         "New Jersey","Colorado","Maryland","Nevada",
                                         "Tennessee","Hawaii","Indiana","Kentucky","Minnesota",
                                         "Oklahoma","Pennsylvania","South Carolina","District of Columbia",
                                         "Kansas","Missouri","Vermont","Virginia","Connecticut",
                                         "Iowa","Louisiana","Ohio","Michigan","South Dakota",
                                         "Arkansas","Delaware","Mississippi","New Mexico","North Dakota",
                                         "Wyoming","Alaska","Maine","Alabama","Idaho","Montana",
                                         "Puerto Rico","Virgin Islands","Guam","West Virginia","Northern Mariana Islands"))
data_florida = us_data_by_state[us_data_by_state$state=="Florida",]
data_florida$newcases = c(0,diff(data_florida$cases))
data_florida$newcases[data_florida$newcases<0] = 0

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

