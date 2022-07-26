# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN SWITZERLAND (GISAID GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 6 July 2021

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
today = as.Date("2021-07-06")
today_num = as.numeric(today)
plotdir = "Switzerland_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata (file version metadata_2021-07-05_12-45.tsv.gz)
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-07-05_12-45.tsv.gz"), col_types = cols(.default = "c")) 
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
range(GISAID$date) # "2020-01-01" "2021-06-27"

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
# 170519         Canada B.1.617+   346
# 170530 Czech Republic B.1.617+    17
# 170532        Switzerland B.1.617+   121
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
# 170612       Switzerland B.1.617+   126
# 170613          Qatar B.1.617+    23
# 170615        Romania B.1.617+    19
# 170616         Switzerland B.1.617+   278
# 170627      Singapore B.1.617+   762
# 170632   South Africa B.1.617+    21
# 170633    South Korea B.1.617+    32
# 170635          Switzerland B.1.617+   264
# 170638         Sweden B.1.617+    42
# 170639    Switzerland B.1.617+   113
# 170641       Thailand B.1.617+    94
# 170651 United Kingdom B.1.617+ 40092
# 170653            USA B.1.617+  2859
# 170656        Vietnam B.1.617+    54

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>100,]$Country))
sel_countries_target
# [1] "Australia"      "Belgium"        "Canada"         "Switzerland"        "France"         "Germany"        "India"          "Ireland"       
# [9] "Italy"          "Japan"          "Switzerland"       "Switzerland"         "Singapore"      "Switzerland"          "Switzerland"    "United Kingdom"
# [17] "USA"    

sel_ref_lineage = "B.1.1.7"

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# [1] "Australia"      "Belgium"        "Canada"         "Switzerland"        "France"         "Germany"        "India"          "Ireland"       
# [9] "Italy"          "Japan"          "Switzerland"       "Switzerland"         "Singapore"      "Switzerland"          "Switzerland"    "United Kingdom"
# [17] "USA"

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# [1] "Australia"      "Belgium"        "Canada"         "Switzerland"        "France"         "Germany"        "India"          "Ireland"       
# [9] "Italy"          "Japan"          "Switzerland"       "Switzerland"         "Singapore"      "Switzerland"          "Switzerland"    "United Kingdom"
# [17] "USA" 

# sel_countries = sel_countries[!(sel_countries %in% c("Japan","USA"))] # Japan is almost only import & for USA we do separate analysis by state





# ANALYSIS OF VOCs IN SWITZERLAND ####

sel_countries = "Switzerland"

tblB117 = table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
tblB117

GISAID_sel = GISAID[GISAID$country %in% sel_countries,]

# use data from Jan  1 onwards
GISAID_sel = GISAID_sel[GISAID_sel$date>=as.Date("2021-01-01"),]
# GISAID_sel = GISAID_sel[-which(GISAID_sel$LINEAGE1=="B.1.617+"&GISAID_sel$date<=as.Date("2021-05-01")),] # B.1.617+ cases before 1st of May are assumed to be all travel related and are removed
nrow(GISAID_sel) # 78813
range(GISAID_sel$date) # "2021-01-01" "2021-06-26"

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))
            
# GISAID_sel = GISAID_sel[GISAID_sel$country_exposure=="India"&GISAID_sel$country!="India",]
# nrow(GISAID_sel[is.na(GISAID_sel$LINEAGE1),]) # 0 unknown pango clade
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$LINEAGE1),]
nrow(GISAID_sel) # 25935

GISAID_sel = GISAID_sel[GISAID_sel$country==GISAID_sel$country_exposure,] # we remove travel-related cases (none indicated here as such)
nrow(GISAID_sel) # 25924

sum(GISAID_sel$LINEAGE1=="B.1.617+") # 126
unique(GISAID_sel$country[GISAID_sel$LINEAGE1=="B.1.1.7"])
sum(GISAID_sel$LINEAGE1=="B.1.1.7") # 4206
sum(GISAID_sel$LINEAGE1=="B.1.1.519") # 0

table(GISAID_sel$LINEAGE1)
table(GISAID_sel$LINEAGE2)

main_lineages = names(table(GISAID_sel$LINEAGE1))[100*table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) > 5]
main_lineages
# "B.1.1.7"  "B.1.177+"  
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
# B.1.1.519   B.1.1.7  B.1.177+   B.1.351  B.1.617+     other       P.1 
# 61     16588      4887       249       590      2897       652 
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1)
GISAID_sel$LINEAGE1 = relevel(GISAID_sel$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE1)
# "B.1.1.7"   "B.1.1.519" "B.1.177+"  "B.1.351"   "B.1.617+"  "other"     "P.1"
levels_LINEAGE1 = c("B.1.1.7","B.1.1.519","B.1.177+",
                    "B.1.351","P.1","B.1.617+","other")
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1, levels=levels_LINEAGE1)

GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = relevel(GISAID_sel$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE2)
# "B.1.1.7"   "B.1.1.519" "B.1.177+"  "B.1.351"   "B.1.617.1" "B.1.617.2" "other"     "P.1"  
levels_LINEAGE2 = c("B.1.1.7","B.1.1.519","B.1.177+",
                    "B.1.351","P.1","B.1.617.1","B.1.617.2","other")
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2, levels=levels_LINEAGE2)

# GISAID_sel = GISAID_sel[GISAID_sel$division!="India",]
table(GISAID_sel$country)


# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
# GISAID_sel = GISAID_sel[-which(grepl("B.1.617", GISAID_sel$pango_lineage, fixed=TRUE)&GISAID_sel$date<=as.Date("2021-04-14")),]

table(GISAID_sel$LINEAGE1)

range(GISAID_sel$date) # "2020-01-01" "2021-06-26"

# aggregated data to make Muller plots of raw data
# aggregate by day to identify days on which INSA performed (days with a lot of sequences)
# we subset the data to just those days to avoid sampling biases (delta infection clusters etc)
data_agbyday2 = as.data.frame(table(GISAID_sel$date, GISAID_sel$LINEAGE2))
colnames(data_agbyday2) = c("date", "LINEAGE2", "count")
data_agbyday2_sum = aggregate(count ~ date, data=data_agbyday2, sum)
data_agbyday2$total = data_agbyday2_sum$count[match(data_agbyday2$date, data_agbyday2_sum$date)]
sum(data_agbyday2[data_agbyday2$LINEAGE2=="B.1.617.1","total"]) == nrow(GISAID_sel) # correct
data_agbyday2$date = as.Date(as.character(data_agbyday2$date))
data_agbyday2$LINEAGE2 = factor(data_agbyday2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyday2$date_num = as.numeric(data_agbyday2$date)
data_agbyday2$prop = data_agbyday2$count/data_agbyday2$total
data_agbyday2$floor_date = NULL
# qplot(data=data_agbyday2, x=date, y=total, colour=total>20, fill=total>20, geom="col")
GISAID_sel$total_sequenced_on_that_day = data_agbyday2$total[match(GISAID_sel$date, data_agbyday2$date)]
# GISAID_sel = GISAID_sel[GISAID_sel$total_sequenced_on_that_day>20,] # dates on which was performed
# nrow(GISAID_sel) # 5710

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


# MULLER PLOT (RAW DATA)
library(scales)
n1 = length(levels(GISAID_sel$LINEAGE1))
lineage_cols1 = hcl(h = seq(15, 320, length = n1), l = 65, c = 200)
lineage_cols1[which(levels(GISAID_sel$LINEAGE1)=="B.1.1.7")] = "#0085FF"
lineage_cols1[which(levels(GISAID_sel$LINEAGE1)=="B.1.351")] = "#9A9D00"
lineage_cols1[which(levels(GISAID_sel$LINEAGE1)=="P.1")] = "cyan3"
lineage_cols1[which(levels(GISAID_sel$LINEAGE1)=="B.1.617+")] = "magenta"
lineage_cols1[which(levels(GISAID_sel$LINEAGE1)=="other")] = "grey75"

n2 = length(levels(GISAID_sel$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.1.7")] = "#0085FF"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.351")] = "#9A9D00"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="P.1")] = "cyan3"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="B.1.617.2")] = "magenta"
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="other")] = "grey75"

muller_switzerland_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SWITZERLAND\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_switzerland_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_muller plots_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_muller plots_raw data.pdf"), width=8, height=6)


# multinomial fits
data_agbyweek2$LINEAGE2 = relevel(data_agbyweek2$LINEAGE2, ref="B.1.1.7")
data_agbyweek2$DATE_NUM = as.numeric(data_agbyweek2$collection_date)

library(nnet)
library(splines)
set.seed(1)
fit1_switzerland_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM), weights=count, data=data_agbyweek2, maxit=1000)
fit2_switzerland_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2), weights=count, data=data_agbyweek2, maxit=1000)
BIC(fit1_switzerland_multi, fit2_switzerland_multi) 
# df      BIC
# fit1_switzerland_multi 14 40394.98
# fit2_switzerland_multi 21 38962.15

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrswitzerland = emtrends(fit2_switzerland_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_switzerland = data.frame(confint(emtrswitzerland, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrswitzerland$contrasts)$p.value)
delta_r_switzerland
# contrast    estimate          SE df   asymp.LCL   asymp.UCL      p.value
# 1  B.1.1.519 - B.1.1.7  0.021297834 0.012803000 NA -0.003795586  0.046391253 4.218317e-01
# 2 (B.1.177+) - B.1.1.7 -0.039363841 0.010198544 NA -0.059352621 -0.019375061 5.467585e-03
# 3    B.1.351 - B.1.1.7  0.006570061 0.006166011 NA -0.005515098  0.018655220 7.724506e-01
# 4        P.1 - B.1.1.7  0.037231319 0.005264798 NA  0.026912505  0.047550133 3.706214e-06
# 5  B.1.617.1 - B.1.1.7 -0.135382286 0.071230517 NA -0.274991533  0.004226962 3.004059e-01
# 6  B.1.617.2 - B.1.1.7  0.089978290 0.005580738 NA  0.079040243  0.100916336 1.897038e-12
# 7      other - B.1.1.7  0.028773558 0.001624913 NA  0.025588787  0.031958328 3.226308e-13


# fitted prop of different LINEAGES in the switzerland today
# 53% [43%-63%] now estimated to be B.1.617.2 across all regions
multinom_preds_today_avg = data.frame(emmeans(fit2_switzerland_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE2         prob           SE df     asymp.LCL    asymp.UCL
# 1   B.1.1.7 3.563482e-01 3.914811e-02 NA  2.796193e-01 4.330771e-01
# 2 B.1.1.519 5.820049e-04 5.695709e-04 NA -5.343335e-04 1.698343e-03
# 3  B.1.177+ 1.765097e-05 1.625588e-05 NA -1.420997e-05 4.951192e-05
# 4   B.1.351 2.427778e-03 1.042974e-03 NA  3.835856e-04 4.471970e-03
# 5       P.1 2.742067e-02 7.949625e-03 NA  1.183969e-02 4.300165e-02
# 6 B.1.617.1 7.136964e-08 3.154441e-07 NA -5.468894e-07 6.896287e-07
# 7 B.1.617.2 5.269709e-01 5.111224e-02 NA  4.267928e-01 6.271491e-01
# 8     other 8.623274e-02 1.229949e-02 NA  6.212618e-02 1.103393e-01

# % non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.6436518 0.5005931 0.7867105 


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
date.to = as.numeric(as.Date("2021-07-31")) # max(GISAID_sel$DATE_NUM)+extrapolate

# multinomial model predictions (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to)))
fit_switzerland_multi_preds = data.frame(predgrid, as.data.frame(predict(fit2_switzerland_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fit_switzerland_multi_preds = gather(fit_switzerland_multi_preds, LINEAGE2, prob, all_of(levels_LINEAGE2), factor_key=TRUE)
fit_switzerland_multi_preds$collection_date = as.Date(fit_switzerland_multi_preds$DATE_NUM, origin="1970-01-01")
fit_switzerland_multi_preds$LINEAGE2 = factor(fit_switzerland_multi_preds$LINEAGE2, levels=levels_LINEAGE2) 

muller_switzerland_mfit = ggplot(data=fit_switzerland_multi_preds, 
                                   aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SWITZERLAND\n(GISAID data, multinomial fit)")
muller_switzerland_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_muller plots_multinom fit.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_muller plots_multinom fit.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_switzerland_raw2 + coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_switzerland_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_muller plots multipanel_multinom fit.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_muller plots multipanel_multinom fit.pdf"), width=10, height=10)





# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions by state with confidence intervals (but slower)
fit_switzerland_multi_preds_withCI = data.frame(emmeans(fit2_switzerland_multi,
                                                        ~ LINEAGE2,
                                                        by=c("DATE_NUM"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=1)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_switzerland_multi_preds_withCI$collection_date = as.Date(fit_switzerland_multi_preds_withCI$DATE_NUM, origin="1970-01-01")
fit_switzerland_multi_preds_withCI$LINEAGE2 = factor(fit_switzerland_multi_preds_withCI$LINEAGE2, levels=levels_LINEAGE2)
fit_switzerland_multi_preds2 = fit_switzerland_multi_preds_withCI


# on logit scale:

ymin = 0.001
ymax = 0.999
fit_switzerland_multi_preds2$asymp.LCL[fit_switzerland_multi_preds2$asymp.LCL<ymin] = ymin
fit_switzerland_multi_preds2$asymp.UCL[fit_switzerland_multi_preds2$asymp.UCL<ymin] = ymin
fit_switzerland_multi_preds2$asymp.UCL[fit_switzerland_multi_preds2$asymp.UCL>ymax] = ymax
fit_switzerland_multi_preds2$prob[fit_switzerland_multi_preds2$prob<ymin] = ymin

plot_switzerland_mfit_logit = qplot(data=fit_switzerland_multi_preds2, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SWITZERLAND\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
                        range=c(0.5, 4), limits=c(1,max(data_agbyweek2$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_switzerland_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_multinom fit_logit scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_multinom fit_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_switzerland_mfit = qplot(data=fit_switzerland_multi_preds2, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SWITZERLAND\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek2,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 5), limits=c(1,max(data_agbyweek2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_switzerland_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_multinom fit_response scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\switzerland_multinom fit_response scale.pdf"), width=10, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot = as.data.frame(get_national_data(countries = "Switzerland"))
cases_tot = cases_tot[cases_tot$date>=as.Date("2020-08-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
# cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
range(cases_tot$date)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
library(mgcv)
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=30, fx=F) + 
                  WEEKDAY, # + 
                # BANKHOLIDAY,
                # s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=poisson(log), data=cases_tot,
) 
BIC(fit_cases)

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_switzerland_multi_preds_withCI$totcases = cases_tot$cases_new[match(round(fit_switzerland_multi_preds_withCI$DATE_NUM),cases_tot$DATE_NUM)]
fit_switzerland_multi_preds_withCI$cases = fit_switzerland_multi_preds_withCI$totcases * fit_switzerland_multi_preds_withCI$prob
fit_switzerland_multi_preds_withCI$cases[fit_switzerland_multi_preds_withCI$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM, at=list(DATE_NUM=seq(date.from, date.to, by=0.5), BANHOLIDAY="no"), type="response"))
fit_switzerland_multi_preds_withCI$smoothed_totcases = cases_emmeans$rate[match(fit_switzerland_multi_preds_withCI$DATE_NUM,cases_emmeans$DATE_NUM)]
fit_switzerland_multi_preds_withCI$smoothed_cases = fit_switzerland_multi_preds_withCI$smoothed_totcases * fit_switzerland_multi_preds_withCI$prob
fit_switzerland_multi_preds_withCI$smoothed_cases[fit_switzerland_multi_preds_withCI$smoothed_cases<=0.001] = NA

ggplot(data=fit_switzerland_multi_preds_withCI[fit_switzerland_multi_preds_withCI$collection_date>=as.Date("2021-02-18"),], 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),max(cases_tot$date)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SWITZERLAND\n(case data & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-02-18"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_switzerland_multi_preds_withCI[fit_switzerland_multi_preds_withCI$collection_date>=as.Date("2021-02-18")&
                                              fit_switzerland_multi_preds_withCI$collection_date<=max(cases_tot$date),], 
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),today), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SWITZERLAND\n(case data & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-02-18"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit case data.png"), width=8, height=6)



# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}


# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to)#,
                                             # BANKHOLIDAY="no"
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[2] = "r"
colnames(avg_r_cases)[5] = "r_LOWER"
colnames(avg_r_cases)[6] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN SWITZERLAND AT MOMENT OF INFECTION BASED ON NEW CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit2_switzerland_multi, ~ LINEAGE2 , at=list(DATE_NUM=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(fit2_switzerland_multi,  ~ LINEAGE2 , by=c("DATE_NUM"),
                                                               var="DATE_NUM", mode="latent",
                                                               at=list(DATE_NUM=d))
                                                out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                # sum(out$estimate*wt) # should sum to zero
                                                return(out) } ))
above_avg_r_variants = above_avg_r_variants0
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE2), 
                                       labels=levels_LINEAGE2)
above_avg_r_variants$variant = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2021-01-04" "2021-07-30"
above_avg_r_variants$avg_r = avg_r_cases$r[match(above_avg_r_variants$collection_date,
                                                 avg_r_cases$DATE)]  # average growth rate of all lineages calculated from case nrs
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                # REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE,
                variant="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$variant = factor(above_avg_r_variants$variant, levels=c(levels_LINEAGE2,"avg"))
above_avg_r_variants$prob = fit_switzerland_multi_preds_withCI$prob[match(interaction(above_avg_r_variants$DATE_NUM,
                                                                      above_avg_r_variants$variant),
                                                          interaction(fit_switzerland_multi_preds_withCI$DATE_NUM,
                                                                      fit_switzerland_multi_preds_withCI$LINEAGE2))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 2
ymin = 1/2
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$variant %in% c("other"))|above_avg_r_variants2$collection_date>max(cases_tot$DATE)),], 
      x=collection_date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
      group=variant, linetype=I(0)) +
  # facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN SWITZERLAND\nAT MOMENT OF INFECTION\n(based on case data & multinomial fit to GISAID data)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  scale_colour_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)
        
