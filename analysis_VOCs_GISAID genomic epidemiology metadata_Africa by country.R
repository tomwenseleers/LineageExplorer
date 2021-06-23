# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN SELECTED AFRICAN COUNTRIES (GISAID METADATA+GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 22 JUNE 2021

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
today = as.Date("2021-06-22")
today_num = as.numeric(today)
today # "2021-06-20"
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID metadata 
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-06-21_05-59.tsv.gz"), col_types = cols(.default = "c")) 
# GISAID = read_tsv(gzfile(".//data//GISAID//metadata.tsv"), col_types = cols(.default = "c")) # using  metadata_tsv_2021_06_18.tar.xz
GISAID = as.data.frame(GISAID)
colnames(GISAID)
# [1] "Virus name"                      "Type"                            "Accession ID"                    "Collection date"                
# [5] "Location"                        "Additional location information" "Sequence length"                 "Host"                           
# [9] "Patient age"                     "Gender"                          "Clade"                           "Pango lineage"                  
# [13] "Pangolin version"                "Variant"                         "AA Substitutions"                "Submission date"                
# [17] "Is reference?"                   "Is complete?"                    "Is high coverage?"               "Is low coverage?"               
# [21] "N-Content"                       "GC-Content"    


## SELECTED  COUNTRIES
# sel_countries = c("Democratic Republic of the Congo","Uganda","Malawi","Kenya") # "Angola","Senegal" "Botswana"
# GISAID = GISAID[grepl(paste0(sel_countries,collapse="|"), GISAID[,"Location"]),]
# nrow(GISAID) # 1956

# library(stringr)
# date_isvalid = sapply(GISAID[,"Collection date"], function (s) str_count(s, pattern = "-")==2)
# sum(date_isvalid) # 1914

GISAID$date = as.Date(GISAID$date) # as.Date(GISAID[,"Collection date"])
GISAID = GISAID[!is.na(GISAID$date),]
# GISAID$host = GISAID$Host
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
range(GISAID$date) # "2020-01-01" "2021-06-16"

firstdetB16172 = GISAID[GISAID$pango_lineage=="B.1.617.2",]
firstdetB16172 = firstdetB16172[!is.na(firstdetB16172$date),]
firstdetB16172 = firstdetB16172[firstdetB16172$date==min(firstdetB16172$date),]
firstdetB16172 # 7 sept 63r old male from Madhya Pradesh

# GISAID = GISAID[grepl("2021-", GISAID$date),]
sum(is.na(GISAID$purpose_of_sequencing)) == nrow(GISAID) # field purpose_of_sequencing left blank unfortunately
nrow(GISAID) # 1901769
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

length(unique(GISAID$country[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 67 countries
table(GISAID$pango_lineage[grepl("B.1.617",GISAID$pango_lineage,fixed=T)])
# B.1.617 B.1.617.1 B.1.617.2 B.1.617.3 
# 2      4388     51068       147

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
# 170497      Australia B.1.617+   252
# 170501        Bahrain B.1.617+    24
# 170502     Bangladesh B.1.617+    44
# 170505        Belgium B.1.617+   247
# 170513         Brazil B.1.617+    16
# 170519         Spain B.1.617+   346
# 170530 Czech Republic B.1.617+    17
# 170532        Denmark B.1.617+   121
# 170542        Finland B.1.617+    19
# 170543         France B.1.617+   141
# 170548        Germany B.1.617+   836
# 170562          India B.1.617+  7435
# 170563      Indonesia B.1.617+    75
# 170564           Iran B.1.617+    11
# 170566        Ireland B.1.617+   299
# 170567         Israel B.1.617+    63
# 170568          Italy B.1.617+   184
# 170570          Japan B.1.617+   170
# 170581     Luxembourg B.1.617+    58
# 170583         Malawi B.1.617+    26
# 170584       Malaysia B.1.617+    12
# 170588         Mexico B.1.617+    48
# 170596          Nepal B.1.617+    34
# 170597    Netherlands B.1.617+    85
# 170598    New Zealand B.1.617+    17
# 170602         Norway B.1.617+    69
# 170611         Poland B.1.617+    71
# 170612       Portugal B.1.617+   126
# 170613          Qatar B.1.617+    23
# 170615        Romania B.1.617+    19
# 170616         Russia B.1.617+   278
# 170627      Singapore B.1.617+   762
# 170632   South Africa B.1.617+    21
# 170633    South Korea B.1.617+    32
# 170635          Spain B.1.617+   264
# 170638         Sweden B.1.617+    42
# 170639    Switzerland B.1.617+   113
# 170641       Thailand B.1.617+    94
# 170651 United Kingdom B.1.617+ 40092
# 170653            USA B.1.617+  2859
# 170656        Vietnam B.1.617+    54

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>100,]$Country))
sel_countries_target
# [1] "Australia"      "Belgium"        "Spain"         "Denmark"        "France"         "Germany"        "India"          "Ireland"       
# [9] "Italy"          "Japan"          "Portugal"       "Russia"         "Singapore"      "Spain"          "Switzerland"    "United Kingdom"
# [17] "USA"    

sel_ref_lineage = "B.1.1.7"

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Belgium"        "Spain"         "Denmark"        "France"         "Germany"        "India"          "Ireland"       
# [9] "Italy"          "Japan"          "Portugal"       "Russia"         "Singapore"      "Spain"          "Switzerland"    "United Kingdom"
# [17] "USA"

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Belgium"        "Spain"         "Denmark"        "France"         "Germany"        "India"          "Ireland"       
# [9] "Italy"          "Japan"          "Portugal"       "Russia"         "Singapore"      "Spain"          "Switzerland"    "United Kingdom"
# [17] "USA" 

# sel_countries = sel_countries[!(sel_countries %in% c("Japan","USA"))] # Japan is almost only import & for USA we do separate analysis by state





# ANALYSIS OF VOCs IN SELECTED AFRICAN COUNTRIES ####

sel_countries = c("Democratic Republic of the Congo","Uganda","Malawi","Kenya") # "Angola","Senegal" "Botswana"
levels_countries = sel_countries

tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB117

GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
nrow(GISAID_sel) # 2619
unique(GISAID_sel$country)

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))
            
# GISAID_sel = GISAID_sel[GISAID_sel$country_exposure=="India"&GISAID_sel$country!="India",]
# nrow(GISAID_sel[is.na(GISAID_sel$LINEAGE1),]) # 0 unknown pango clade
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$LINEAGE1),]
nrow(GISAID_sel) # 2619

GISAID_sel = GISAID_sel[GISAID_sel$country==GISAID_sel$country_exposure,] # we remove travel-related cases
nrow(GISAID_sel) # 2619

sum(GISAID_sel$LINEAGE1=="B.1.617+") # 81
unique(GISAID_sel$country[GISAID_sel$LINEAGE1=="B.1.1.7"])
sum(GISAID_sel$LINEAGE1=="B.1.1.7") # 392
sum(GISAID_sel$LINEAGE1=="B.1.1.519") # 0
sum(GISAID_sel$LINEAGE1=="B.1.351") # 477

table(GISAID_sel$LINEAGE1)
table(GISAID_sel$LINEAGE2)

main_lineages = names(table(GISAID_sel$LINEAGE1))[100*table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) > 3]
main_lineages
# "A.23.1"  "B.1"     "B.1.1"   "B.1.1.7" "B.1.351" "B.1.416" "C.16" 
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
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1)
GISAID_sel$LINEAGE1 = relevel(GISAID_sel$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE1)
levels_LINEAGE1 = c("B.1.1.7",levels(GISAID_sel$LINEAGE1)[!levels(GISAID_sel$LINEAGE1) %in% c("B.1.1.7","B.1.617+","B.1.617.1","B.1.617.2","other")],
                    "B.1.617+","other")
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1, levels=levels_LINEAGE1)

GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = relevel(GISAID_sel$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE2)
# "B.1.1.7"   "B.1"       "B.1.1"     "B.1.160"   "B.1.177+"  "B.1.351"   "B.1.617.1" "B.1.617.2" "B.1.91"    "other"     "P.1"  
levels_LINEAGE2 = c("B.1.1.7",levels(GISAID_sel$LINEAGE2)[!levels(GISAID_sel$LINEAGE2) %in% c("B.1.1.7","B.1.617+","B.1.617.1","B.1.617.2","other")],
                    "B.1.617.1","B.1.617.2","other")
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2, levels=levels_LINEAGE2)

# GISAID_sel = GISAID_sel[GISAID_sel$division!="India",]
table(GISAID_sel$country)

GISAID_sel$country = factor(GISAID_sel$country, levels=levels_countries)

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
# GISAID_sel = GISAID_sel[-which(grepl("B.1.617", GISAID_sel$pango_lineage, fixed=TRUE)&GISAID_sel$date<=as.Date("2021-04-14")),]

table(GISAID_sel$LINEAGE2)

range(GISAID_sel$date) # "2020-01-24" "2021-06-02"

GISAID_sel = GISAID_sel[GISAID_sel$date>="2020-11-01",]

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
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

# aggregated by week & country for selected variant lineages
data_agbyweek_bycountry2 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$country, GISAID_sel$LINEAGE2))
colnames(data_agbyweek_bycountry2) = c("floor_date", "country", "LINEAGE2", "count")
data_agbyweek_bycountry2_sum = aggregate(count ~ floor_date+country, data=data_agbyweek_bycountry2, sum)
data_agbyweek_bycountry2$total = data_agbyweek_bycountry2_sum$count[match(interaction(data_agbyweek_bycountry2$floor_date,data_agbyweek_bycountry2$country), 
                                                                interaction(data_agbyweek_bycountry2_sum$floor_date,data_agbyweek_bycountry2_sum$country))]
sum(data_agbyweek_bycountry2[data_agbyweek_bycountry2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_sel) # correct
data_agbyweek_bycountry2$collection_date = as.Date(as.character(data_agbyweek_bycountry2$floor_date))
data_agbyweek_bycountry2$LINEAGE2 = factor(data_agbyweek_bycountry2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek_bycountry2$collection_date_num = as.numeric(data_agbyweek_bycountry2$collection_date)
data_agbyweek_bycountry2$prop = data_agbyweek_bycountry2$count/data_agbyweek_bycountry2$total
data_agbyweek_bycountry2$floor_date = NULL
data_agbyweek_bycountry2$country = factor(data_agbyweek_bycountry2$country, levels=levels_countries)

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

# muller plot, overall
muller_africa_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2), width=1, position="fill") +
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_africa_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plots_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plots_raw data.pdf"), width=8, height=6)

# muller plot, by country
muller_africa_bycountry_raw2 = ggplot(data=data_agbyweek_bycountry2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ country) +
  geom_col(aes(colour=NULL, fill=LINEAGE2), width=I(7), position="fill") +
  # geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_africa_bycountry_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by country_muller plots_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by country_muller plots_raw data.pdf"), width=8, height=6)


# multinomial fits
data_agbyweek_bycountry2$LINEAGE2 = relevel(data_agbyweek_bycountry2$LINEAGE2, ref="B.1.1.7")
data_agbyweek_bycountry2$DATE_NUM = as.numeric(data_agbyweek_bycountry2$collection_date)

library(nnet)
library(splines)
set.seed(1)
fit1_africa_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)+country, weights=count, data=data_agbyweek_bycountry2, maxit=1000)
fit2_africa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)+country, weights=count, data=data_agbyweek_bycountry2, maxit=1000)
BIC(fit1_africa_multi, fit2_africa_multi) 
#             df      BIC
# fit1_africa_multi 45 3338.943
# fit2_africa_multi 54 3302.149

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrafrica = emtrends(fit1_africa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_africa = data.frame(confint(emtrafrica, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrafrica$contrasts)$p.value)
delta_r_africa
# contrast     estimate           SE df    asymp.LCL     asymp.UCL      p.value
# 1      A.23 - B.1.1.7 -0.091547729 0.010807932 NA -0.1127308859 -0.070364573 6.430304e-10
# 2    A.23.1 - B.1.1.7 -0.049757026 0.004358294 NA -0.0582991246 -0.041214927 2.949863e-13
# 3       B.1 - B.1.1.7 -0.061701373 0.004007935 NA -0.0695567816 -0.053845965 2.137179e-13
# 4     B.1.1 - B.1.1.7 -0.044901831 0.007172338 NA -0.0589593563 -0.030844306 1.128545e-06
# 5   B.1.351 - B.1.1.7 -0.021235381 0.003120084 NA -0.0273506324 -0.015120129 1.751507e-07
# 6   B.1.525 - B.1.1.7 -0.003612109 0.006925874 NA -0.0171865727  0.009962356 9.824483e-01
# 7 B.1.617.1 - B.1.1.7  0.031546400 0.016549286 NA -0.0008896036  0.063982404 3.240238e-01
# 8 B.1.617.2 - B.1.1.7  0.072844877 0.010112703 NA  0.0530243426  0.092665411 4.525167e-08
# 9     other - B.1.1.7 -0.062852738 0.004053698 NA -0.0707978403 -0.054907636 2.137179e-13


# fitted prop of different LINEAGES in the africa today
multinom_preds_today_avg = data.frame(emmeans(fit1_africa_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE2         prob           SE df     asymp.LCL    asymp.UCL
# 1    B.1.1.7 4.420254e-02 1.789822e-02 NA  9.122664e-03 7.928242e-02
# 2       A.23 3.292225e-09 7.367667e-09 NA -1.114814e-08 1.773259e-08
# 3     A.23.1 6.676990e-05 5.779578e-05 NA -4.650774e-05 1.800475e-04
# 4        B.1 1.631274e-05 9.819370e-06 NA -2.932876e-06 3.555835e-05
# 5      B.1.1 1.309549e-05 1.367024e-05 NA -1.369769e-05 3.988868e-05
# 6    B.1.351 6.010401e-03 3.201895e-03 NA -2.651989e-04 1.228600e-02
# 7    B.1.525 3.873887e-03 3.559363e-03 NA -3.102336e-03 1.085011e-02
# 8  B.1.617.1 9.706502e-03 1.210247e-02 NA -1.401391e-02 3.342691e-02
# 9  B.1.617.2 9.360993e-01 2.938448e-02 NA  8.785068e-01 9.936918e-01
# 10     other 1.118525e-05 7.030696e-06 NA -2.594662e-06 2.496516e-05

# % non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.9557975 0.8610596 1.0505353


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2020-11-01"))
date.to = as.numeric(as.Date("2021-07-31")) # max(GISAID_sel$DATE_NUM)+extrapolate

# multinomial model predictions by country (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to), country=levels(data_agbyweek_bycountry2$country)))
fit_africa_multi_preds = data.frame(predgrid, as.data.frame(predict(fit1_africa_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fit_africa_multi_preds = gather(fit_africa_multi_preds, LINEAGE2, prob, all_of(levels_LINEAGE2), factor_key=TRUE)
fit_africa_multi_preds$collection_date = as.Date(fit_africa_multi_preds$DATE_NUM, origin="1970-01-01")
fit_africa_multi_preds$LINEAGE2 = factor(fit_africa_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 

muller_africa_bycountry_mfit = ggplot(data=fit_africa_multi_preds, 
                                   aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ country) +
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)")
muller_africa_bycountry_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_by country_muller plots_multinom fit.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_by country_muller plots_multinom fit.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_africa_bycountry_raw2 + coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_africa_bycountry_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_by country_muller plots multipanel_multinom fit.png"), width=10, height=10)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_by country_muller plots multipanel_multinom fit.pdf"), width=10, height=10)


# multinomial model predictions on avg across countries
fit_africa_multi_preds_withCI = data.frame(emmeans(fit1_africa_multi,
                                                           ~ LINEAGE2,
                                                           by=c("DATE_NUM"),
                                                           at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                           mode="prob", df=NA))
fit_africa_multi_preds_withCI$collection_date = as.Date(fit_africa_multi_preds_withCI$DATE_NUM, origin="1970-01-01")
fit_africa_multi_preds_withCI$LINEAGE2 = factor(fit_africa_multi_preds_withCI$LINEAGE2, levels=levels_LINEAGE2)

muller_africa_mfit = ggplot(data=fit_africa_multi_preds_withCI, 
                                      aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ country) +
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)")
muller_africa_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plot_multinom fit avg.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plot_multinom fit avg.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_africa_raw2 + coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_africa_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plot multipanel_multinom fit avg.png"), width=10, height=10)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plot multipanel_multinom fit avg.pdf"), width=10, height=10)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# overall average multinomial model predictions over all selected countries with confidence intervals


fit_africa_multi_preds_withCI[fit_africa_multi_preds_withCI$collection_date==(as.Date("2021-06-22")-2)&fit_africa_multi_preds_withCI$LINEAGE2=="B.1.617.2",]
# LINEAGE2 DATE_NUM      prob         SE df asymp.LCL asymp.UCL collection_date
# 339 B.1.617.2    18798 0.9281428 0.03116128 NA 0.8670678 0.9892178      2021-06-20

# fit_africa_multi_preds2 = fit_africa_multi_preds_bystate # without CIs
# fit_africa_multi_preds2$asymp.LCL = NA
# fit_africa_multi_preds2$asymp.UCL = NA


# on logit scale:

fit_africa_multi_preds2 = fit_africa_multi_preds_withCI
ymin = 0.001
ymax = 0.999
fit_africa_multi_preds2$asymp.LCL[fit_africa_multi_preds2$asymp.LCL<ymin] = ymin
fit_africa_multi_preds2$asymp.UCL[fit_africa_multi_preds2$asymp.UCL<ymin] = ymin
fit_africa_multi_preds2$asymp.UCL[fit_africa_multi_preds2$asymp.UCL>ymax] = ymax
fit_africa_multi_preds2$prob[fit_africa_multi_preds2$prob<ymin] = ymin

plot_africa_mfit_logit = qplot(data=fit_africa_multi_preds2, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~ country) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)") +
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
                        range=c(0.5/2, 3/2), limits=c(1,max(data_agbyweek2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_africa_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_multinom fit avg_logit scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_multinom fit avg_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_africa_mfit = qplot(data=fit_africa_multi_preds2, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~ country) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)") +
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
                        range=c(0.5/2, 3/2), limits=c(1,max(data_agbyweek2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_africa_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_multinom fit avg_response scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_multinom fit avg_response scale.pdf"), width=10, height=6)




# multinomial model predictions with confidence intervals by country
fit_africa_multi_preds_bycountry_withCI = data.frame(emmeans(fit1_africa_multi,
                                                        ~ LINEAGE2,
                                                        by=c("DATE_NUM","country"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                        mode="prob", df=NA))
fit_africa_multi_preds_bycountry_withCI$collection_date = as.Date(fit_africa_multi_preds_bycountry_withCI$DATE_NUM, origin="1970-01-01")
fit_africa_multi_preds_bycountry_withCI$LINEAGE2 = factor(fit_africa_multi_preds_bycountry_withCI$LINEAGE2, levels=levels_LINEAGE2)
fit_africa_multi_preds_bycountry_withCI$country = factor(fit_africa_multi_preds_bycountry_withCI$country, levels=levels_countries)
fit_africa_multi_preds3 = fit_africa_multi_preds_bycountry_withCI

fit_africa_multi_preds_bycountry_withCI[fit_africa_multi_preds_bycountry_withCI$collection_date==(as.Date("2021-06-22")-2)&fit_africa_multi_preds_bycountry_withCI$LINEAGE2=="B.1.617.2",]
#       LINEAGE2 DATE_NUM                          country      prob          SE df asymp.LCL asymp.UCL collection_date
# 339  B.1.617.2    18798 Democratic Republic of the Congo 0.9949276 0.005100895 NA 0.9849300 1.0049252      2021-06-20
# 729  B.1.617.2    18798                           Uganda 0.9607608 0.038417975 NA 0.8854629 1.0360586      2021-06-20
# 1119 B.1.617.2    18798                           Malawi 0.9763441 0.020628610 NA 0.9359128 1.0167755      2021-06-20
# 1509 B.1.617.2    18798                            Kenya 0.7805385 0.077915140 NA 0.6278277 0.9332494      2021-06-20

# on logit scale:

ymin = 0.001
ymax = 0.999
fit_africa_multi_preds3$asymp.LCL[fit_africa_multi_preds3$asymp.LCL<ymin] = ymin
fit_africa_multi_preds3$asymp.UCL[fit_africa_multi_preds3$asymp.UCL<ymin] = ymin
fit_africa_multi_preds3$asymp.UCL[fit_africa_multi_preds3$asymp.UCL>ymax] = ymax
fit_africa_multi_preds3$prob[fit_africa_multi_preds3$prob<ymin] = ymin

plot_africa_bycountry_mfit_logit = qplot(data=fit_africa_multi_preds3, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~ country) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek_bycountry2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek2$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-11-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_africa_bycountry_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by country_multinom fit_logit scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by country_multinom fit_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_africa_bycountry_mfit = qplot(data=fit_africa_multi_preds3, x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~ country) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-11-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-11-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek_bycountry2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek2$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_africa_bycountry_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by country_multinom fit_response scale.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by country_multinom fit_response scale.pdf"), width=10, height=6)





# project multinomial fit onto cases ####
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

