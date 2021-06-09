# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE USA (GISAID GENOMIC EPIDEMIOLOGY METADATA)
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
today = as.Date("2021-06-05")
today_num = as.numeric(today)
today # "2021-06-04"
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
GISAID = GISAID[GISAID$date>=as.Date("2020-01-01"),]
range(GISAID$date) # "2020-01-01" "2021-05-31"

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
# 167996      Australia B.1.617+   144
# 168000        Bahrain B.1.617+    19
# 168001     Bangladesh B.1.617+    20
# 168004        Belgium B.1.617+   135
# 168029 Czech Republic B.1.617+    13
# 168031        Denmark B.1.617+   105
# 168041         France B.1.617+    91
# 168046        Germany B.1.617+   421
# 168060          India B.1.617+  3932
# 168061      Indonesia B.1.617+    32
# 168062           Iran B.1.617+    11
# 168064        Ireland B.1.617+   184
# 168065         Israel B.1.617+    36
# 168066          Italy B.1.617+   133
# 168068          Japan B.1.617+   165
# 168086         Mexico B.1.617+    28
# 168095    Netherlands B.1.617+    50
# 168096    New Zealand B.1.617+    15
# 168100         Norway B.1.617+    37
# 168109         Poland B.1.617+    41
# 168110       Portugal B.1.617+    77
# 168113        Romania B.1.617+    19
# 168125      Singapore B.1.617+   292
# 168130   South Africa B.1.617+    14
# 168132          Spain B.1.617+    86
# 168135         Sweden B.1.617+    18
# 168136    Switzerland B.1.617+    75
# 168148 United Kingdom B.1.617+ 12253
# 168150            USA B.1.617+  1658

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

range(GISAID_sel$date) # "2020-01-03" "2021-06-01"

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
# 1        B.1 - B.1.1.7 -0.047397512 0.0003975381 NA -0.048176672 -0.046618351 4.551914e-15
# 2      B.1.2 - B.1.1.7 -0.065623372 0.0003772277 NA -0.066362725 -0.064884020 4.551914e-15
# 3  B.1.1.519 - B.1.1.7 -0.051627731 0.0009920338 NA -0.053572082 -0.049683381 4.551914e-15
# 4    B.1.351 - B.1.1.7 -0.010651130 0.0018366240 NA -0.014250847 -0.007051413 2.286988e-07
# 5    B.1.429 - B.1.1.7 -0.057214440 0.0005614526 NA -0.058314867 -0.056114013 4.551914e-15
# 6    B.1.526 - B.1.1.7 -0.017357306 0.0007556718 NA -0.018838395 -0.015876216 4.551914e-15
# 7        P.1 - B.1.1.7  0.005939672 0.0012129550 NA  0.003562324  0.008317021 1.856053e-05
# 8  B.1.617.1 - B.1.1.7  0.024899792 0.0063396675 NA  0.012474272  0.037325312 1.085243e-03
# 9  B.1.617.2 - B.1.1.7  0.070329433 0.0023281113 NA  0.065766419  0.074892447 4.551914e-15
# 10     other - B.1.1.7 -0.045866540 0.0003325045 NA -0.046518237 -0.045214844 4.551914e-15


# fitted prop of different LINEAGES in the USA today
# 76% [74%-78%] now estimated to be B.1.617.2 across all regions
multinom_preds_today_avg = data.frame(emmeans(fit3_usa_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1    B.1.1.7 0.665085245 9.537197e-03 NA 0.6463926827 0.683777807
# 2        B.1 0.001076100 4.081378e-05 NA 0.0009961066 0.001156094
# 3      B.1.2 0.001750906 5.527119e-05 NA 0.0016425761 0.001859235
# 4  B.1.1.519 0.003013862 1.803726e-04 NA 0.0026603378 0.003367386
# 5    B.1.351 0.005903325 5.750782e-04 NA 0.0047761926 0.007030458
# 6    B.1.429 0.003212750 1.273375e-04 NA 0.0029631730 0.003462327
# 7    B.1.526 0.028532949 1.085721e-03 NA 0.0264049747 0.030660923
# 8        P.1 0.103897931 3.941544e-03 NA 0.0961726456 0.111623216
# 9  B.1.617.1 0.002086205 5.948964e-04 NA 0.0009202299 0.003252181
# 10 B.1.617.2 0.170913478 1.079366e-02 NA 0.1497582914 0.192068664
# 11     other 0.014527250 3.663528e-04 NA 0.0138092116 0.015245288

# 33% [30%-37%] non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.3349148 0.3001037 0.3697258  


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

fit_usa_multi_preds_bystate_withCI[fit_usa_multi_preds_bystate_withCI$collection_date==as.Date("2021-06-06")&fit_usa_multi_preds_bystate_withCI$LINEAGE2=="B.1.617.2",]
#        LINEAGE2 DATE_NUM         STATE       prob          SE df  asymp.LCL  asymp.UCL collection_date
# 351  B.1.617.2    18784       Arizona 0.09958093 0.026942362 NA 0.04677487 0.15238699      2021-06-06
# 780  B.1.617.2    18784    California 0.27977698 0.020279348 NA 0.24003019 0.31952377      2021-06-06
# 1209 B.1.617.2    18784      Colorado 0.24743118 0.022232888 NA 0.20385552 0.29100684      2021-06-06
# 1638 B.1.617.2    18784       Florida 0.10825029 0.013355180 NA 0.08207462 0.13442596      2021-06-06
# 2067 B.1.617.2    18784      Illinois 0.11127194 0.014085080 NA 0.08366569 0.13887819      2021-06-06
# 2496 B.1.617.2    18784       Indiana 0.18384167 0.025889384 NA 0.13309941 0.23458393      2021-06-06
# 2925 B.1.617.2    18784        Kansas 0.21435103 0.028973733 NA 0.15756356 0.27113850      2021-06-06
# 3354 B.1.617.2    18784      Maryland 0.13447865 0.023356819 NA 0.08870013 0.18025717      2021-06-06
# 3783 B.1.617.2    18784 Massachusetts 0.18646085 0.016972203 NA 0.15319594 0.21972576      2021-06-06
# 4212 B.1.617.2    18784     Minnesota 0.05709520 0.009506131 NA 0.03846352 0.07572687      2021-06-06
# 4641 B.1.617.2    18784      Missouri 0.24650529 0.044577100 NA 0.15913578 0.33387480      2021-06-06
# 5070 B.1.617.2    18784      Nebraska 0.11007176 0.037991369 NA 0.03561005 0.18453347      2021-06-06
# 5499 B.1.617.2    18784    New Jersey 0.27974676 0.028308616 NA 0.22426290 0.33523063      2021-06-06
# 5928 B.1.617.2    18784      New York 0.18491593 0.018647379 NA 0.14836774 0.22146412      2021-06-06
# 6357 B.1.617.2    18784        Oregon 0.05574835 0.020272323 NA 0.01601533 0.09548137      2021-06-06
# 6786 B.1.617.2    18784         Texas 0.18787694 0.019606641 NA 0.14944863 0.22630525      2021-06-06
# 7215 B.1.617.2    18784          Utah 0.41954464 0.046294134 NA 0.32880981 0.51027948      2021-06-06
# 7644 B.1.617.2    18784      Virginia 0.19976433 0.031238521 NA 0.13853795 0.26099071      2021-06-06
# 8073 B.1.617.2    18784    Washington 0.18881715 0.017361974 NA 0.15478830 0.22284599      2021-06-06
# 8502 B.1.617.2    18784     Wisconsin 0.12031708 0.025105098 NA 0.07111199 0.16952217      2021-06-06

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
  coord_cartesian(xlim=as.Date(c("2020-11-01","2021-06-14")),
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

