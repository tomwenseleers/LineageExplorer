# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN VARIOUS COUNTRIES (INDIA, UK, BELGIUM) BASED ON ANALYSIS OF GISAID GENOMIC EPIDEMIOLOGY METADATA
# T. Wenseleers
# last update 8 MAY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-08")
today_num = as.numeric(today)
today # "2021-05-08"
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata (file version metadata_2021-05-07_12-55.tsv.gz)
GISAID = read_tsv(".//data//GISAID_genomic_epidemiology//metadata.tsv", col_types = cols(.default = "c"))
GISAID = as.data.frame(GISAID)
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
range(GISAID$date) # "2021-12-06" "2021-05-05"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 12 dec Bihar India

GISAID = GISAID[grepl("2021-", GISAID$date),]
unique(GISAID$host)
# [1] "Human"                               "Hunan"                               "Felis catus"                         "Canis lupus familiaris"             
# [5] "Environment"                         "Gorilla gorilla gorilla"             "Prionailurus bengalensis euptilurus" "Panthera leo"                       
# [9] "Mink"                                "Mus musculus"                        "Panthera tigris jacksoni" 
GISAID = GISAID[GISAID$host=="Human",]
sum(is.na(GISAID$purpose_of_sequencing)) == nrow(GISAID) # field purpose_of_sequencing left blank unfortunately
GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
range(GISAID$date) # "2021-01-01" "2021-05-05"
nrow(GISAID) # 866125
GISAID$Week = lubridate::week(GISAID$date)
GISAID$DATE_NUM = as.numeric(GISAID$date)
colnames(GISAID)
unique(GISAID$region)
# "Asia"          "Europe"        "Africa"        "South America" "Oceania"       "North America"
unique(GISAID$country)
unique(GISAID$division) # = city or province or region, sometimes just country
unique(GISAID$location) # = city

length(unique(GISAID$country[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 36 countries
table(GISAID$pango_lineage[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])
# B.1.617 B.1.617.1 B.1.617.2 B.1.617.3 
# 160      1547      1143        63 

GISAID$pango_lineage[grepl("B.1.177",GISAID$pango_lineage,fixed=T)] = "B.1.177+"
sel_target_VOC = "B.1.617"
GISAID$LINEAGE1 = GISAID$pango_lineage
GISAID$LINEAGE2 = GISAID$pango_lineage
GISAID[grepl(sel_target_VOC, GISAID$LINEAGE1, fixed=TRUE),"LINEAGE1"] = paste0(sel_target_VOC,"+") # in LINEAGE1 we recode B.1.617.1,2&3 all as B.1.617+

table_country_lineage = as.data.frame(table(GISAID$country, GISAID$LINEAGE1))
colnames(table_country_lineage) = c("Country","Lineage","Count")
tblB1617 = table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage, fixed=T)&table_country_lineage$Count>10,]
tblB1617
# Country  Lineage Count
# 94951      Australia B.1.617+    78
# 94955        Bahrain B.1.617+    22
# 94959        Belgium B.1.617+    14
# 94982        Denmark B.1.617+    39
# 94991         France B.1.617+    15
# 94996        Germany B.1.617+    58
# 95006          India B.1.617+  1301
# 95010        Ireland B.1.617+    20
# 95036    New Zealand B.1.617+    13
# 95057      Singapore B.1.617+   155
# 95067    Switzerland B.1.617+    17
# 95078 United Kingdom B.1.617+   833
# 95079            USA B.1.617+   270

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>10,]$Country))
sel_countries_target
# "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "New Zealand"   
# "Singapore"      "Switzerland"    "United Kingdom" "USA" 

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "New Zealand"   
# "Singapore"      "Switzerland"    "United Kingdom" "USA"  

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# "Australia"      "Bahrain"        "Belgium"        "Denmark"        "France"         "Germany"        "India"          "Ireland"        "New Zealand"   
# "Singapore"      "Switzerland"    "United Kingdom" "USA"   

sel_ref_lineage = "B.1.1.7"
tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB117
# Country Lineage  Count
# 41959      Australia B.1.1.7    330
# 41963        Bahrain B.1.1.7     12
# 41967        Belgium B.1.1.7   9436
# 41990        Denmark B.1.1.7  31635
# 41999         France B.1.1.7  15078
# 42004        Germany B.1.1.7  56372
# 42014          India B.1.1.7    400
# 42018        Ireland B.1.1.7   8342
# 42044    New Zealand B.1.1.7    115
# 42065      Singapore B.1.1.7    160
# 42075    Switzerland B.1.1.7  12148
# 42086 United Kingdom B.1.1.7 203373
# 42087            USA B.1.1.7  78739

data.frame(Country=tblB1617$Country, Lineage="B.1.617", Perc=100*tblB1617$Count / (tblB1617$Count+tblB117$Count))
#           Country Lineage        Perc
# 1       Australia B.1.617 19.11764706
# 2         Bahrain B.1.617 64.70588235
# 3         Belgium B.1.617  0.14814815
# 4         Denmark B.1.617  0.12312938
# 5          France B.1.617  0.09938382
# 6         Germany B.1.617  0.10278221
# 7           India B.1.617 76.48442093
# 8         Ireland B.1.617  0.23917723
# 9     New Zealand B.1.617 10.15625000
# 10      Singapore B.1.617 49.20634921
# 11    Switzerland B.1.617  0.13974517
# 12 United Kingdom B.1.617  0.40792141
# 13            USA B.1.617  0.34173322


GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_sel) # 684 932

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))




# ANALYSIS VOCs IN INDIA ####

GISAID_india = GISAID_sel[GISAID_sel$country=="India",]
nrow(GISAID_india) # 4097

unique(GISAID_india$division) # best data for West Bengal, Maharashtra & Karnataka (B.1.617 most common, B.1.1.7 most common in Punjab IN, Telangana)
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.617+"])
sum(GISAID_india$LINEAGE1=="B.1.617+") # 1301
unique(GISAID_india$division[GISAID_india$LINEAGE1=="B.1.1.7"])
sum(GISAID_india$LINEAGE1=="B.1.1.7") # 400

table(GISAID_india$LINEAGE1)
table(GISAID_india$LINEAGE2)

main_lineages = names(table(GISAID_india$LINEAGE1))[100*table(GISAID_india$LINEAGE1)/sum(table(GISAID_india$LINEAGE1)) > 3]
main_lineages
# "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.1.7"   "B.1.36"    "B.1.36.29" "B.1.617+"  "B.1.618"  
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526")
main_lineages = union(main_lineages, VOCs)
GISAID_india$LINEAGE1[!(GISAID_india$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_india$LINEAGE2[!(GISAID_india$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_india$LINEAGE1 = factor(GISAID_india$LINEAGE1)
GISAID_india$LINEAGE1 = relevel(GISAID_india$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_india$LINEAGE1)
# "B.1.1.7"   "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.1.318" "B.1.36"    "B.1.36.29" "B.1.525"   "B.1.351"   "B.1.618"   "P.1"       "B.1.617+" 
# "other" 
levels_LINEAGE1 = c("B.1.1.7","B.1","B.1.1","B.1.1.216","B.1.1.306","B.1.1.318","B.1.36","B.1.36.29",
                    "B.1.525","B.1.351","B.1.618","P.1","B.1.617+","other")
GISAID_india$LINEAGE1 = factor(GISAID_india$LINEAGE1, levels=levels_LINEAGE1)

GISAID_india$LINEAGE2 = factor(GISAID_india$LINEAGE2)
GISAID_india$LINEAGE2 = relevel(GISAID_india$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_india$LINEAGE2)
# "B.1.1.7"   "B.1"       "B.1.1"     "B.1.1.216" "B.1.1.306" "B.1.1.318" "B.1.351"   "B.1.36"    "B.1.36.29" "B.1.525"   "B.1.617.1" "B.1.617.2" "B.1.618"  
# "other"     "P.1"
levels_LINEAGE2 = c("B.1.1.7","B.1","B.1.1","B.1.1.216","B.1.1.306","B.1.1.318","B.1.36","B.1.36.29",
                    "B.1.525","B.1.351","B.1.618","P.1","B.1.617.1","B.1.617.2","other")
GISAID_india$LINEAGE2 = factor(GISAID_india$LINEAGE2, levels=levels_LINEAGE2)
firstB16172 = GISAID_india[GISAID_india$LINEAGE2=="B.1.617.2",]
firstB16172 = firstB16172[firstB16172$date==min(firstB16172$date),]
firstB16172

GISAID_india = GISAID_india[GISAID_india$division!="India",]
sel_states = names(table(GISAID_india$division))[table(GISAID_india$division) > 200]
# "Andhra Pradesh" "Chhattisgarh"   "Gujarat"        "Karnataka"      "Maharashtra"    "Telangana"      "West Bengal" 
GISAID_india = GISAID_india[GISAID_india$division %in% sel_states,]
GISAID_india$division = factor(GISAID_india$division)
GISAID_india$division = relevel(GISAID_india$division, ref="Maharashtra")
# levels_STATES = c("Maharashtra","Karnataka","West Bengal")
levels_STATES = c("Maharashtra","Chhattisgarh","Gujarat","West Bengal","Telangana","Andhra Pradesh","Karnataka")
GISAID_india$division = factor(GISAID_india$division, levels=levels_STATES)
# levels_STATES = levels(GISAID_india$division)
levels(GISAID_india$division)
# "Maharashtra"    "Andhra Pradesh" "Chhattisgarh"   "Gujarat"        "Karnataka"      "Telangana"      "West Bengal"  
table(GISAID_india$division)
# Maharashtra Andhra Pradesh   Chhattisgarh        Gujarat      Karnataka      Telangana    West Bengal 
#        1401            317            310            334            292            211            843 

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of India
data_agbyweek1 = as.data.frame(table(GISAID_india$Week, GISAID_india$LINEAGE1))
colnames(data_agbyweek1) = c("Week", "LINEAGE1", "count")
data_agbyweek1_sum = aggregate(count ~ Week, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$Week, data_agbyweek1_sum$Week)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_india) # correct
data_agbyweek1$Week = as.numeric(as.character(data_agbyweek1$Week))
data_agbyweek1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) + 3.5
data_agbyweek1$LINEAGE1 = factor(data_agbyweek1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

data_agbyweek2 = as.data.frame(table(GISAID_india$Week, GISAID_india$LINEAGE2))
colnames(data_agbyweek2) = c("Week", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ Week, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$Week, data_agbyweek2_sum$Week)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_india) # correct
data_agbyweek2$Week = as.numeric(as.character(data_agbyweek2$Week))
data_agbyweek2$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek2$Week - 1 ) + 3.5
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total


# aggregated by week and state for selected variant lineages
data_agbyweekregion1 = as.data.frame(table(GISAID_india$Week, GISAID_india$division, GISAID_india$LINEAGE1))
colnames(data_agbyweekregion1) = c("Week", "division", "LINEAGE1", "count")
data_agbyweekregion1_sum = aggregate(count ~ Week + division, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$Week,data_agbyweekregion1$division), 
                                                                  interaction(data_agbyweekregion1_sum$Week,data_agbyweekregion1_sum$division))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE1=="B.1.617+","total"]) == nrow(GISAID_india) # correct
data_agbyweekregion1$Week = as.numeric(as.character(data_agbyweekregion1$Week))
data_agbyweekregion1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) + 3.5
data_agbyweekregion1$LINEAGE1 = factor(data_agbyweekregion1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweekregion1$division = factor(data_agbyweekregion1$division, levels=levels_STATES)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion2 = as.data.frame(table(GISAID_india$Week, GISAID_india$division, GISAID_india$LINEAGE2))
colnames(data_agbyweekregion2) = c("Week", "division", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ Week + division, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$Week,data_agbyweekregion2$division), 
                                                                  interaction(data_agbyweekregion2_sum$Week,data_agbyweekregion2_sum$division))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_india) # correct
data_agbyweekregion2$Week = as.numeric(as.character(data_agbyweekregion2$Week))
data_agbyweekregion2$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion2$Week - 1 ) + 3.5
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$division = factor(data_agbyweekregion2$division, levels=levels_STATES)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]


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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
  facet_wrap(~ division) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.pdf"), width=8, height=6)



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
# fit3_india_multi fits best (lowest BIC), but fit1_india_multi almost as good, so we will use that simpler model

# equivalent fit with B.1.617.1,2&3 all recoded to B.1.617+
fit1_india_multi1 = nnet::multinom(LINEAGE1 ~ division + ns(DATE_NUM, df=2), data=GISAID_india, maxit=1000)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
max(GISAID_india$date) # 2021-04-24
emtrindia = emtrends(fit1_india_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india = data.frame(confint(emtrindia, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrindia$contrasts)$p.value)
delta_r_india
#               contrast     estimate           SE df    asymp.LCL    asymp.UCL     p.value
# 1        B.1 - B.1.1.7 -0.014702491 4.314671e-06 NA -0.014710947 -0.014694034 6.37268e-14
# 2      B.1.1 - B.1.1.7 -0.053740586 4.880637e-06 NA -0.053750152 -0.053731020 6.37268e-14
# 3  B.1.1.216 - B.1.1.7 -0.044519341 5.295595e-06 NA -0.044529720 -0.044508962 6.37268e-14
# 4  B.1.1.306 - B.1.1.7 -0.024028683 5.397706e-06 NA -0.024039262 -0.024018104 6.37268e-14
# 5  B.1.1.318 - B.1.1.7 -0.010225377 5.366675e-05 NA -0.010330562 -0.010120192 6.37268e-14
# 6     B.1.36 - B.1.1.7 -0.053008704 4.857640e-06 NA -0.053018225 -0.052999183 6.37268e-14
# 7  B.1.36.29 - B.1.1.7 -0.040869889 4.519973e-06 NA -0.040878748 -0.040861030 6.37268e-14
# 8    B.1.525 - B.1.1.7  0.006263769 8.531250e-06 NA  0.006247048  0.006280490 6.37268e-14
# 9    B.1.351 - B.1.1.7  0.021186758 7.934710e-06 NA  0.021171207  0.021202310 6.37268e-14
# 10   B.1.618 - B.1.1.7 -0.007212274 5.841102e-06 NA -0.007223722 -0.007200825 6.37268e-14
# 11       P.1 - B.1.1.7 -0.002328037 3.800482e-05 NA -0.002402525 -0.002253549 6.37268e-14
# 12 B.1.617.1 - B.1.1.7  0.022791282 3.638474e-06 NA  0.022784151  0.022798414 6.37268e-14
# 13 B.1.617.2 - B.1.1.7  0.094613178 5.643592e-06 NA  0.094602117  0.094624239 6.37268e-14
# 14     other - B.1.1.7 -0.031940881 3.823739e-06 NA -0.031948375 -0.031933386 6.37268e-14

# avg growth advantage of B.1.617+ over B.1.1.7 :
max(GISAID_india$date) # 2021-04-24
emtrindia1 = emtrends(fit1_india_multi1, trt.vs.ctrl ~ LINEAGE1,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india1 = data.frame(confint(emtrindia1, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrindia1$contrasts)$p.value)
delta_r_india1
#               contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1         B.1 - B.1.1.7  3.786155e-02 0.0117903330 NA  1.475292e-02  6.097018e-02 1.851633e-02
# 2       B.1.1 - B.1.1.7  3.441713e-02 0.0151494595 NA  4.724737e-03  6.410953e-02 2.007606e-01
# 3   B.1.1.216 - B.1.1.7 -9.123635e-02 0.0251496350 NA -1.405287e-01 -4.194397e-02 4.924095e-03
# 4   B.1.1.306 - B.1.1.7 -4.494566e-02 0.0206583342 NA -8.543525e-02 -4.456069e-03 2.424513e-01
# 5   B.1.1.318 - B.1.1.7 -6.431407e+01 0.0022148054 NA -6.431841e+01 -6.430973e+01 1.210143e-14
# 6      B.1.36 - B.1.1.7 -1.393216e-01 0.0294823264 NA -1.971059e-01 -8.153733e-02 8.102533e-05
# 7   B.1.36.29 - B.1.1.7 -8.261673e-02 0.0184845645 NA -1.188458e-01 -4.638765e-02 2.266531e-04
# 8     B.1.525 - B.1.1.7 -5.609454e-02 0.0300432191 NA -1.149782e-01  2.789085e-03 4.098794e-01
# 9     B.1.351 - B.1.1.7 -4.064231e-02 0.0277029340 NA -9.493906e-02  1.365445e-02 6.677445e-01
# 10    B.1.618 - B.1.1.7 -9.991769e-02 0.0263634988 NA -1.515892e-01 -4.824618e-02 2.829737e-03
# 11        P.1 - B.1.1.7 -1.332133e+04 0.0001686798 NA -1.332133e+04 -1.332133e+04 1.210143e-14
# 12 (B.1.617+) - B.1.1.7  5.219053e-02 0.0100924200 NA  3.240975e-02  7.197131e-02 1.238689e-05
# 13      other - B.1.1.7 -2.455534e-02 0.0127490678 NA -4.954305e-02  4.323758e-04 3.744123e-01


# pairwise growth rate difference (differences in growth rate per day) 
max(GISAID_india$date) # 2021-04-24
emtrindia_pairw = emtrends(fit1_india_multi, pairwise ~ LINEAGE2,  
                           var="DATE_NUM",  mode="latent",
                           at=list(DATE_NUM=max(GISAID_india$DATE_NUM)))
delta_r_india_pairw = data.frame(confint(emtrindia_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindia_pairw$contrasts)$p.value)
delta_r_india_pairw
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1           B.1.1.7 - B.1  0.0147024905 4.314671e-06 NA  0.0146940339  0.0147109471 6.861178e-14
# 2         B.1.1.7 - B.1.1  0.0537405862 4.880637e-06 NA  0.0537310203  0.0537501521 6.861178e-14
# 3     B.1.1.7 - B.1.1.216  0.0445193410 5.295595e-06 NA  0.0445089619  0.0445297202 6.861178e-14
# 4     B.1.1.7 - B.1.1.306  0.0240286829 5.397706e-06 NA  0.0240181036  0.0240392622 6.861178e-14
# 5     B.1.1.7 - B.1.1.318  0.0102253767 5.366675e-05 NA  0.0101201918  0.0103305616 6.861178e-14
# 6        B.1.1.7 - B.1.36  0.0530087038 4.857640e-06 NA  0.0529991830  0.0530182246 6.861178e-14
# 7     B.1.1.7 - B.1.36.29  0.0408698890 4.519973e-06 NA  0.0408610300  0.0408787480 6.861178e-14
# 8       B.1.1.7 - B.1.525 -0.0062637689 8.531250e-06 NA -0.0062804899 -0.0062470480 6.861178e-14
# 9       B.1.1.7 - B.1.351 -0.0211867584 7.934710e-06 NA -0.0212023102 -0.0211712067 6.861178e-14
# 10      B.1.1.7 - B.1.618  0.0072122735 5.841102e-06 NA  0.0072008252  0.0072237219 6.861178e-14
# 11          B.1.1.7 - P.1  0.0023280374 3.800482e-05 NA  0.0022535493  0.0024025254 6.861178e-14
# 12    B.1.1.7 - B.1.617.1 -0.0227912823 3.638474e-06 NA -0.0227984136 -0.0227841511 6.861178e-14
# 13    B.1.1.7 - B.1.617.2 -0.0946131780 5.643592e-06 NA -0.0946242392 -0.0946021167 6.861178e-14
# 14        B.1.1.7 - other  0.0319408806 3.823739e-06 NA  0.0319333862  0.0319483750 6.861178e-14
# 15            B.1 - B.1.1  0.0390380957 4.781477e-06 NA  0.0390287241  0.0390474672 6.861178e-14
# 16        B.1 - B.1.1.216  0.0298168505 5.247793e-06 NA  0.0298065650  0.0298271360 6.861178e-14
# 17        B.1 - B.1.1.306  0.0093261924 5.294039e-06 NA  0.0093158162  0.0093365685 6.861178e-14
# 18        B.1 - B.1.1.318 -0.0044771139 5.363133e-05 NA -0.0045822294 -0.0043719984 6.861178e-14
# 19           B.1 - B.1.36  0.0383062132 4.757294e-06 NA  0.0382968891  0.0383155373 6.861178e-14
# 20        B.1 - B.1.36.29  0.0261673985 4.393896e-06 NA  0.0261587866  0.0261760104 6.861178e-14
# 21          B.1 - B.1.525 -0.0209662595 8.486802e-06 NA -0.0209828933 -0.0209496256 6.861178e-14
# 22          B.1 - B.1.351 -0.0358892490 7.904035e-06 NA -0.0359047406 -0.0358737573 6.861178e-14
# 23          B.1 - B.1.618 -0.0074902170 5.893712e-06 NA -0.0075017685 -0.0074786656 6.861178e-14
# 24              B.1 - P.1 -0.0123744532 3.801632e-05 NA -0.0124489638 -0.0122999426 6.861178e-14
# 25        B.1 - B.1.617.1 -0.0374937729 3.516957e-06 NA -0.0375006660 -0.0374868798 6.861178e-14
# 26        B.1 - B.1.617.2 -0.1093156685 5.568855e-06 NA -0.1093265833 -0.1093047538 6.861178e-14
# 27            B.1 - other  0.0172383900 3.687719e-06 NA  0.0172311622  0.0172456178 6.861178e-14
# 28      B.1.1 - B.1.1.216 -0.0092212451 5.511828e-06 NA -0.0092320481 -0.0092104422 6.861178e-14
# 29      B.1.1 - B.1.1.306 -0.0297119033 5.783704e-06 NA -0.0297232392 -0.0297005675 6.861178e-14
# 30      B.1.1 - B.1.1.318 -0.0435152095 5.370063e-05 NA -0.0436204608 -0.0434099582 6.861178e-14
# 31         B.1.1 - B.1.36 -0.0007318824 5.143582e-06 NA -0.0007419637 -0.0007218012 6.861178e-14
# 32      B.1.1 - B.1.36.29 -0.0128706972 4.944301e-06 NA -0.0128803878 -0.0128610065 6.861178e-14
# 33        B.1.1 - B.1.525 -0.0600043551 8.778630e-06 NA -0.0600215609 -0.0599871493 6.861178e-14
# 34        B.1.1 - B.1.351 -0.0749273446 8.237437e-06 NA -0.0749434897 -0.0749111995 6.861178e-14
# 35        B.1.1 - B.1.618 -0.0465283127 6.128099e-06 NA -0.0465403235 -0.0465163018 6.861178e-14
# 36            B.1.1 - P.1 -0.0514125488 3.805267e-05 NA -0.0514871307 -0.0513379670 6.861178e-14
# 37      B.1.1 - B.1.617.1 -0.0765318685 4.301638e-06 NA -0.0765402996 -0.0765234375 6.861178e-14
# 38      B.1.1 - B.1.617.2 -0.1483537642 6.148446e-06 NA -0.1483658149 -0.1483417134 6.861178e-14
# 39          B.1.1 - other -0.0217997056 4.322423e-06 NA -0.0218081774 -0.0217912338 6.861178e-14
# 40  B.1.1.216 - B.1.1.306 -0.0204906582 6.147255e-06 NA -0.0205027066 -0.0204786098 6.861178e-14
# 41  B.1.1.216 - B.1.1.318 -0.0342939644 5.376310e-05 NA -0.0343993381 -0.0341885907 6.861178e-14
# 42     B.1.1.216 - B.1.36  0.0084893627 5.519028e-06 NA  0.0084785456  0.0085001798 6.861178e-14
# 43  B.1.1.216 - B.1.36.29 -0.0036494520 5.397276e-06 NA -0.0036600305 -0.0036388736 6.861178e-14
# 44    B.1.1.216 - B.1.525 -0.0507831100 9.045646e-06 NA -0.0508008391 -0.0507653808 6.861178e-14
# 45    B.1.1.216 - B.1.351 -0.0657060995 8.514533e-06 NA -0.0657227877 -0.0656894113 6.861178e-14
# 46    B.1.1.216 - B.1.618 -0.0373070675 6.420150e-06 NA -0.0373196508 -0.0372944843 6.861178e-14
# 47        B.1.1.216 - P.1 -0.0421913037 3.810009e-05 NA -0.0422659785 -0.0421166289 6.861178e-14
# 48  B.1.1.216 - B.1.617.1 -0.0673106234 4.805592e-06 NA -0.0673200422 -0.0673012046 6.861178e-14
# 49  B.1.1.216 - B.1.617.2 -0.1391325190 6.504294e-06 NA -0.1391452672 -0.1391197708 6.861178e-14
# 50      B.1.1.216 - other -0.0125784605 4.810510e-06 NA -0.0125878889 -0.0125690321 6.861178e-14
# 51  B.1.1.306 - B.1.1.318 -0.0138033062 5.376635e-05 NA -0.0139086863 -0.0136979261 6.861178e-14
# 52     B.1.1.306 - B.1.36  0.0289800209 5.684500e-06 NA  0.0289688795  0.0289911623 6.861178e-14
# 53  B.1.1.306 - B.1.36.29  0.0168412061 5.420353e-06 NA  0.0168305824  0.0168518298 6.861178e-14
# 54    B.1.1.306 - B.1.525 -0.0302924518 9.172578e-06 NA -0.0303104297 -0.0302744739 6.861178e-14
# 55    B.1.1.306 - B.1.351 -0.0452154413 8.628471e-06 NA -0.0452323528 -0.0451985298 6.861178e-14
# 56    B.1.1.306 - B.1.618 -0.0168164094 6.741966e-06 NA -0.0168296234 -0.0168031954 6.861178e-14
# 57        B.1.1.306 - P.1 -0.0217006455 3.815683e-05 NA -0.0217754315 -0.0216258595 6.861178e-14
# 58  B.1.1.306 - B.1.617.1 -0.0468199652 4.742874e-06 NA -0.0468292611 -0.0468106694 6.861178e-14
# 59  B.1.1.306 - B.1.617.2 -0.1186418609 6.473239e-06 NA -0.1186545482 -0.1186291736 6.861178e-14
# 60      B.1.1.306 - other  0.0079121977 4.896459e-06 NA  0.0079026008  0.0079217946 6.861178e-14
# 61     B.1.1.318 - B.1.36  0.0427833271 5.372109e-05 NA  0.0426780357  0.0428886185 6.861178e-14
# 62  B.1.1.318 - B.1.36.29  0.0306445124 5.368703e-05 NA  0.0305392877  0.0307497370 6.861178e-14
# 63    B.1.1.318 - B.1.525 -0.0164891456 5.407454e-05 NA -0.0165951297 -0.0163831614 6.861178e-14
# 64    B.1.1.318 - B.1.351 -0.0314121351 5.398931e-05 NA -0.0315179522 -0.0313063180 6.861178e-14
# 65    B.1.1.318 - B.1.618 -0.0030131032 5.382187e-05 NA -0.0031185921 -0.0029076142 6.861178e-14
# 66        B.1.1.318 - P.1 -0.0078973393 6.562871e-05 NA -0.0080259692 -0.0077687094 6.861178e-14
# 67  B.1.1.318 - B.1.617.1 -0.0330166590 5.360475e-05 NA -0.0331217224 -0.0329115956 6.861178e-14
# 68  B.1.1.318 - B.1.617.2 -0.1048385546 5.377772e-05 NA -0.1049439570 -0.1047331522 6.861178e-14
# 69      B.1.1.318 - other  0.0217155039 5.362927e-05 NA  0.0216103925  0.0218206154 6.861178e-14
# 70     B.1.36 - B.1.36.29 -0.0121388147 4.882126e-06 NA -0.0121483835 -0.0121292459 6.861178e-14
# 71       B.1.36 - B.1.525 -0.0592724727 8.848666e-06 NA -0.0592898158 -0.0592551296 6.861178e-14
# 72       B.1.36 - B.1.351 -0.0741954622 8.294824e-06 NA -0.0742117197 -0.0741792046 6.861178e-14
# 73       B.1.36 - B.1.618 -0.0457964303 6.215382e-06 NA -0.0458086122 -0.0457842483 6.861178e-14
# 74           B.1.36 - P.1 -0.0506806664 3.806753e-05 NA -0.0507552774 -0.0506060554 6.861178e-14
# 75     B.1.36 - B.1.617.1 -0.0757999861 4.285221e-06 NA -0.0758083850 -0.0757915872 6.861178e-14
# 76     B.1.36 - B.1.617.2 -0.1476218817 6.141892e-06 NA -0.1476339196 -0.1476098438 6.861178e-14
# 77         B.1.36 - other -0.0210678232 4.266772e-06 NA -0.0210761859 -0.0210594605 6.861178e-14
# 78    B.1.36.29 - B.1.525 -0.0471336580 8.696736e-06 NA -0.0471507033 -0.0471166127 6.861178e-14
# 79    B.1.36.29 - B.1.351 -0.0620566475 8.125028e-06 NA -0.0620725722 -0.0620407227 6.861178e-14
# 80    B.1.36.29 - B.1.618 -0.0336576155 6.087638e-06 NA -0.0336695471 -0.0336456840 6.861178e-14
# 81        B.1.36.29 - P.1 -0.0385418517 3.804710e-05 NA -0.0386164226 -0.0384672807 6.861178e-14
# 82  B.1.36.29 - B.1.617.1 -0.0636611714 3.807328e-06 NA -0.0636686336 -0.0636537091 6.861178e-14
# 83  B.1.36.29 - B.1.617.2 -0.1354830670 5.832804e-06 NA -0.1354944991 -0.1354716349 6.861178e-14
# 84      B.1.36.29 - other -0.0089290085 3.881151e-06 NA -0.0089366154 -0.0089214015 6.861178e-14
# 85      B.1.525 - B.1.351 -0.0149229895 1.058515e-05 NA -0.0149437360 -0.0149022430 6.861178e-14
# 86      B.1.525 - B.1.618  0.0134760424 9.266384e-06 NA  0.0134578807  0.0134942042 6.861178e-14
# 87          B.1.525 - P.1  0.0085918063 3.867599e-05 NA  0.0085160027  0.0086676098 6.861178e-14
# 88    B.1.525 - B.1.617.1 -0.0165275134 8.212121e-06 NA -0.0165436089 -0.0165114179 6.861178e-14
# 89    B.1.525 - B.1.617.2 -0.0883494090 9.271517e-06 NA -0.0883675809 -0.0883312372 6.861178e-14
# 90        B.1.525 - other  0.0382046495 8.325466e-06 NA  0.0381883319  0.0382209671 6.861178e-14
# 91      B.1.351 - B.1.618  0.0283990319 8.742268e-06 NA  0.0283818974  0.0284161665 6.861178e-14
# 92          B.1.351 - P.1  0.0235147958 3.855312e-05 NA  0.0234392330  0.0235903585 6.861178e-14
# 93    B.1.351 - B.1.617.1 -0.0016045239 7.592062e-06 NA -0.0016194041 -0.0015896437 6.861178e-14
# 94    B.1.351 - B.1.617.2 -0.0734264196 8.705504e-06 NA -0.0734434820 -0.0734093571 6.861178e-14
# 95        B.1.351 - other  0.0531276390 7.731112e-06 NA  0.0531124863  0.0531427917 6.861178e-14
# 96          B.1.618 - P.1 -0.0048842361 3.813001e-05 NA -0.0049589696 -0.0048095027 6.861178e-14
# 97    B.1.618 - B.1.617.1 -0.0300035558 5.435551e-06 NA -0.0300142093 -0.0299929024 6.861178e-14
# 98    B.1.618 - B.1.617.2 -0.1018254515 6.951590e-06 NA -0.1018390763 -0.1018118266 6.861178e-14
# 99        B.1.618 - other  0.0247286071 5.527729e-06 NA  0.0247177729  0.0247394412 6.861178e-14
# 100       P.1 - B.1.617.1 -0.0251193197 3.794625e-05 NA -0.0251936930 -0.0250449464 6.861178e-14
# 101       P.1 - B.1.617.2 -0.0969412153 3.819251e-05 NA -0.0970160713 -0.0968663594 6.861178e-14
# 102           P.1 - other  0.0296128432 3.796071e-05 NA  0.0295384416  0.0296872448 6.861178e-14
# 103 B.1.617.1 - B.1.617.2 -0.0718218956 4.894182e-06 NA -0.0718314881 -0.0718123032 6.861178e-14
# 104     B.1.617.1 - other  0.0547321629 2.971165e-06 NA  0.0547263395  0.0547379863 6.861178e-14
# 105     B.1.617.2 - other  0.1265540585 5.289363e-06 NA  0.1265436916  0.1265644255 6.861178e-14

# PS confidence intervals with this model artificially narrow - probably better to use mblogit model to take into account overdispersion

# PLOT MULTINOMIAL FIT

extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
date.to = max(GISAID_india$DATE_NUM)+extrapolate

fit_india_multi_predsbystate = data.frame(emmeans(fit1_india_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "division"),
                                                  at=list(DATE_NUM=seq(date.from, date.to)), 
                                                  mode="prob", df=NA))
fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
fit_india_multi_predsbystate$LINEAGE2 = factor(fit_india_multi_predsbystate$LINEAGE2, levels=levels_LINEAGE2) 
fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$division, levels=levels_STATES) 

fit_india_multi_preds = data.frame(emmeans(fit1_india_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to)), 
                                           mode="prob", df=NA))
fit_india_multi_preds$collection_date = as.Date(fit_india_multi_preds$DATE_NUM, origin="1970-01-01")
fit_india_multi_preds$LINEAGE2 = factor(fit_india_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 

muller_india_mfit = ggplot(data=fit_india_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
           xmax=as.Date("2021-06-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\n(multinomial fit)")
muller_india_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_multinom fit.pdf"), width=8, height=6)


ggarrange(muller_india_raw2+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_india_mfit+ggtitle("Multinomial fit"), ncol=1)


muller_indiabystate_mfit = ggplot(data=fit_india_multi_predsbystate, 
                                  aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ STATE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
           xmax=as.Date("2021-06-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\n(multinomial fit)")
muller_indiabystate_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.pdf"), width=8, height=6)

ggarrange(muller_indiabystate_raw2+
            coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0)))+
            ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\nRaw GISAID data"), 
          muller_indiabystate_mfit+ggtitle("\nMultinomial fit")+
            coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-01"))), nrow=2)



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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-01")), 
                     expand=c(0,0)) +
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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-01")), 
                     expand=c(0,0)) +
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