# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN VARIOUS COUNTRIES (INDIA, UK, BELGIUM) BASED ON GISAID DATA

# last update 4 MAY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-04")
today_num = as.numeric(today)
today # "2021-05-04"
plotdir = "VOCs_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID metadata (file version metadata_tsv_2021_05_04.tar.xz)
GISAID = read_tsv(".//data//GISAID//metadata.tsv", col_types = cols(.default = "c"))
GISAID = as.data.frame(GISAID)
GISAID = GISAID[grepl("2021-", GISAID[,"Collection date"]),]
nrow(GISAID) # 838048
GISAID[,"Collection date"] = as.Date(GISAID[,"Collection date"])
GISAID = GISAID[!is.na(GISAID[,"Collection date"]),]
range(GISAID[,"Collection date"]) # "2021-01-01" "2021-04-30"
nrow(GISAID) # 818009
GISAID$Week = lubridate::week(GISAID[,"Collection date"])
GISAID$DATE_NUM = as.numeric(GISAID[,"Collection date"])
colnames(GISAID)
GISAID$Continent = sapply(unlist(GISAID[,"Location"]), function(loc) trimws(strsplit(loc,"/", fixed = TRUE)[[1]][1]))
unique(GISAID$Continent)
# "Oceania"       "Europe"        "North America" "Asia"          "South America" "Africa"  
GISAID$Country = sapply(unlist(GISAID[,"Location"]), function(loc) trimws(strsplit(loc,"/", fixed = TRUE)[[1]][2]))
unique(GISAID$Country)
GISAID$Region = sapply(unlist(GISAID[,"Location"]), function(loc) trimws(strsplit(loc,"/", fixed = TRUE)[[1]][3]))
unique(GISAID$Region)


sel_target_VOC = "B.1.617"
GISAID[grepl(sel_target_VOC, GISAID[,"Pango lineage"], fixed=TRUE),"Pango lineage"] = sel_target_VOC

table_country_lineage = as.data.frame(table(GISAID$Country, GISAID[,"Pango lineage"]))
colnames(table_country_lineage) = c("Country","Lineage","Count")
table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>10,]
#               Country Lineage Count
# 120440      Australia B.1.617    47
# 120444        Bahrain B.1.617    22
# 120491        Germany B.1.617    26
# 120503          India B.1.617  1235
# 120538    New Zealand B.1.617    11
# 120561      Singapore B.1.617   121
# 120572    Switzerland B.1.617    11
# 120582 United Kingdom B.1.617   728
# 120584            USA B.1.617   256

sel_countries_target = unique(as.character(table_country_lineage[grepl(sel_target_VOC, table_country_lineage$Lineage)&table_country_lineage$Count>10,]$Country))
sel_countries_target
# "Australia"      "Bahrain"        "Germany"        "India"          "New Zealand"    "Singapore"      "United Kingdom" "USA"

sel_countries_ref = as.character(table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries_target,]$Country)
sel_countries_ref
# "Australia"      "Bahrain"        "Germany"        "India"          "New Zealand"    "Singapore"      "United Kingdom" "USA"  

sel_countries = intersect(sel_countries_target, sel_countries_ref)
sel_countries
# "Australia"      "Bahrain"        "Germany"        "India"          "New Zealand"    "Singapore"      "United Kingdom" "USA"  

sel_ref_lineage = "B.1.1.7"
table_country_lineage[table_country_lineage$Lineage==sel_ref_lineage&table_country_lineage$Count>10&table_country_lineage$Country %in% sel_countries,]
#              Country Lineage  Count
# 47265      Australia B.1.1.7    304
# 47269        Bahrain B.1.1.7     12
# 47316        Germany B.1.1.7  46469
# 47328          India B.1.1.7    299
# 47363    New Zealand B.1.1.7    112
# 47386      Singapore B.1.1.7    156
# 47408 United Kingdom B.1.1.7 198903
# 47410            USA B.1.1.7  70693

GISAID_sel = GISAID[GISAID$Country %in% sel_countries,]
nrow(GISAID_sel) # 533 832
sum(grepl("hCoV-19", GISAID_sel[,"Virus name"])) # 533 832

table(GISAID_sel[,"Pango lineage"],GISAID_sel$Country)

# PS some fields are missing from GISAID data - check how these can be obtained
# Additional location information - ZIP code
# Sampling strategy : Baseline surveillance / Active surveillance / Longitudinal sampling on same patient / Non-sentinel-surveillance (Hospital)
# Additional host information


# ANALYSIS VOCs ALL COUNTRIES COMBINED ####

GISAID_sel$E484Q = grepl("484Q",GISAID_sel[,"AA Substitutions"])
GISAID_sel$LINEAGE = GISAID_sel[,"Pango lineage"]
GISAID_sel$LINEAGE[grepl("B.1.177",GISAID_sel$LINEAGE, fixed=TRUE)] = "B.1.177+" 
GISAID_sel$B1617 = grepl("B.1.617",GISAID_sel$LINEAGE, fixed=TRUE)
GISAID_sel$LINEAGE1 = GISAID_sel$LINEAGE # with extra category "other" for minority or non-VOC strains
GISAID_sel$LINEAGE2 = GISAID_sel$LINEAGE # LINEAGE2 with B.1.617 with & without E484Q
GISAID_sel[GISAID_sel$B1617 & (!GISAID_sel$E484Q),"LINEAGE2"] = "B.1.617_without_E484Q"
GISAID_sel[GISAID_sel$B1617 & GISAID_sel$E484Q,"LINEAGE2"] = "B.1.617_with_E484Q"
GISAID_sel$LINEAGE1[GISAID_sel$B1617] = "B.1.617+"
main_lineages = names(table(GISAID_sel$LINEAGE1))[(100*table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) > 1)&(table(GISAID_sel$LINEAGE1)>10)]
VOCs = c("B.1.617","B.1.617+","B.1.618","B.1.617_with_E484Q","B.1.617_without_E484Q","B.1.1.7","B.1.1.177+","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526")
main_lineages = union(main_lineages, VOCs)
GISAID_sel$LINEAGE1[!(GISAID_sel$LINEAGE %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_sel$LINEAGE2[!(GISAID_sel$LINEAGE %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_sel$LINEAGE = factor(GISAID_sel$LINEAGE)
GISAID_sel$LINEAGE = relevel(GISAID_sel$LINEAGE, ref="B.1.1.7") # we code UK strain as the reference level
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1)
GISAID_sel$LINEAGE1 = relevel(GISAID_sel$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = relevel(GISAID_sel$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
GISAID_sel$Country = factor(GISAID_sel$Country)
levels(GISAID_sel$LINEAGE1)
levels(GISAID_sel$LINEAGE2)
levels(GISAID_sel$Country)

# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_all_multi1 = nnet::multinom(LINEAGE1 ~ Country + DATE_NUM, data=GISAID_sel, maxit=1000)
fit2_all_multi1 = nnet::multinom(LINEAGE1 ~ Country * DATE_NUM, data=GISAID_sel, maxit=1000)
fit3_all_multi1 = nnet::multinom(LINEAGE1 ~ Country + ns(DATE_NUM, df=2), data=GISAID_sel, maxit=1000)
fit4_all_multi1 = nnet::multinom(LINEAGE1 ~ Country * ns(DATE_NUM, df=2), data=GISAID_sel, maxit=1000)
BIC(fit1_all_multi1, fit2_all_multi1, fit3_all_multi1, fit4_all_multi1)
# fit3_all_multi1 fits best (lowest BIC)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
max(as.Date(unlist(GISAID_sel[,"Collection date"]))) # 2021-04-28
emtr_all = emtrends(fit3_all_multi1, trt.vs.ctrl ~ LINEAGE1,  
                    var="DATE_NUM",  mode="latent",
                    at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)), 
                    adjust="none", df=NA)
delta_r_all = data.frame(confint(emtr_all)$contrasts)
delta_r_all
#                contrast     estimate           SE df    asymp.LCL    asymp.UCL
# 1        A.26 - B.1.1.7 -0.003117571 4.272723e-05 NA -0.003201315 -0.003033827
# 2   B.1.1.117 - B.1.1.7 -0.006477870 3.121176e-05 NA -0.006539044 -0.006416696
# 3   B.1.1.207 - B.1.1.7 -0.065885737 2.492884e-06 NA -0.065890623 -0.065880851
# 4   B.1.1.318 - B.1.1.7  0.012456003 2.514590e-06 NA  0.012451075  0.012460932
# 5  B.1.160.32 - B.1.1.7 -0.044399356 1.529046e-05 NA -0.044429325 -0.044369387
# 6  B.1.177.83 - B.1.1.7 -0.056335359 5.609591e-06 NA -0.056346354 -0.056324365
# 7     B.1.351 - B.1.1.7 -0.002652430 9.794316e-07 NA -0.002654350 -0.002650511
# 8     B.1.429 - B.1.1.7 -0.045860152 4.000542e-07 NA -0.045860936 -0.045859367
# 9     B.1.474 - B.1.1.7 -0.008519707 3.922997e-05 NA -0.008596597 -0.008442818
# 10    B.1.525 - B.1.1.7  0.003521109 1.483996e-06 NA  0.003518200  0.003524017
# 11    B.1.526 - B.1.1.7 -0.003530566 5.117686e-07 NA -0.003531570 -0.003529563
# 12    B.1.581 - B.1.1.7 -0.009531769 3.698265e-05 NA -0.009604254 -0.009459284
# 13    B.1.617 - B.1.1.7  0.067185947 1.519963e-06 NA  0.067182967  0.067188926
# 14    B.1.618 - B.1.1.7  0.010978936 4.639514e-06 NA  0.010969843  0.010988029
# 15      other - B.1.1.7 -0.056984522 2.280423e-07 NA -0.056984968 -0.056984075
# 16        P.1 - B.1.1.7  0.037237685 8.513062e-07 NA  0.037236016  0.037239353
# 17        W.2 - B.1.1.7 -0.007708989 3.901124e-05 NA -0.007785450 -0.007632529

# tests for differences in growth rate with UK variant B.1.1.7
contr_all_againstB117 = data.frame(emtr_all$contrasts)
contr_all_againstB117
#                contrast     estimate           SE df       z.ratio p.value
# 1        A.26 - B.1.1.7 -0.003117571 4.272723e-05 NA     -72.96451       0
# 2   B.1.1.117 - B.1.1.7 -0.006477870 3.121176e-05 NA    -207.54583       0
# 3   B.1.1.207 - B.1.1.7 -0.065885737 2.492884e-06 NA  -26429.51950       0
# 4   B.1.1.318 - B.1.1.7  0.012456003 2.514590e-06 NA    4953.49233       0
# 5  B.1.160.32 - B.1.1.7 -0.044399356 1.529046e-05 NA   -2903.72913       0
# 6  B.1.177.83 - B.1.1.7 -0.056335359 5.609591e-06 NA  -10042.68541       0
# 7     B.1.351 - B.1.1.7 -0.002652430 9.794316e-07 NA   -2708.13230       0
# 8     B.1.429 - B.1.1.7 -0.045860152 4.000542e-07 NA -114634.85849       0
# 9     B.1.474 - B.1.1.7 -0.008519707 3.922997e-05 NA    -217.17343       0
# 10    B.1.525 - B.1.1.7  0.003521109 1.483996e-06 NA    2372.72077       0
# 11    B.1.526 - B.1.1.7 -0.003530566 5.117686e-07 NA   -6898.75555       0
# 12    B.1.581 - B.1.1.7 -0.009531769 3.698265e-05 NA    -257.73626       0
# 13    B.1.617 - B.1.1.7  0.067185947 1.519963e-06 NA   44202.35754       0
# 14    B.1.618 - B.1.1.7  0.010978936 4.639514e-06 NA    2366.39777       0
# 15      other - B.1.1.7 -0.056984522 2.280423e-07 NA -249885.80440       0
# 16        P.1 - B.1.1.7  0.037237685 8.513062e-07 NA   43741.82261       0
# 17        W.2 - B.1.1.7 -0.007708989 3.901124e-05 NA    -197.60943       0









# ANALYSIS VOCs IN INDIA ####

GISAID_bengal = GISAID[grepl("India / West Bengal",GISAID$Location),]
GISAID_bengal$COUNTRY = "India"
GISAID_bengal$STATE = "West Bengal"
GISAID_maharashtra = GISAID[grepl("India / Maharashtra",GISAID$Location),]
GISAID_maharashtra$COUNTRY = "India"
GISAID_maharashtra$STATE = "Maharashtra"
GISAID_karnataka = GISAID[grepl("India / Karnataka",GISAID$Location),]
GISAID_karnataka$COUNTRY = "India"
GISAID_karnataka$STATE = "Karnataka"
# GISAID_punjab = GISAID[grepl("India / Punjab",GISAID$Location),]
# GISAID_punjab$COUNTRY = "India"
# GISAID_punjab$STATE = "Punjab"
# GISAID_telangana = GISAID[grepl("India / Telangana",GISAID$Location),]
# GISAID_telangana$COUNTRY = "India"
# GISAID_telangana$STATE = "Telangana"
# GISAID_delhi = GISAID[grepl("India / Delhi",GISAID$Location),] # no data for Delhi in 2021
# GISAID_india_bystate = rbind(GISAID_bengal, GISAID_maharashtra, GISAID_karnataka, GISAID_punjab, GISAID_telangana)
GISAID_india_bystate = rbind(GISAID_bengal, GISAID_maharashtra, GISAID_karnataka)
unique(GISAID_india_bystate$STATE)

GISAID_india_bystate$B617 = unlist(GISAID_india_bystate[,"Pango lineage"])=="B.1.617"
sum(GISAID_india_bystate$B617) # 340
GISAID_india_bystate$B617_1 = grepl("484Q",unlist(GISAID_india_bystate[,"AA Substitutions"]))&(unlist(GISAID_india_bystate[,"Pango lineage"])=="B.1.617") 
sum(GISAID_india_bystate$B617_1) # 307
GISAID_india_bystate$B617_2_3 = (!grepl("484Q",unlist(GISAID_india_bystate[,"AA Substitutions"])))&(unlist(GISAID_india_bystate[,"Pango lineage"])=="B.1.617")
sum(GISAID_india_bystate$B617_2_3) # 33
length(unique(unlist(GISAID_india_bystate[,"Pango lineage"]))) # 83
table(unlist(GISAID_india_bystate[,"Pango lineage"]))

colnames(GISAID_india_bystate)
GISAID_india_bystate = as.data.frame(GISAID_india_bystate)
GISAID_india_bystate$E484Q = grepl("484Q",GISAID_india_bystate[,"AA Substitutions"])
GISAID_india_bystate$B617 = (GISAID_india_bystate[,"Pango lineage"]=="B.1.617") 
GISAID_india_bystate$LINEAGE = GISAID_india_bystate[,"Pango lineage"]
GISAID_india_bystate$LINEAGE1 = GISAID_india_bystate$LINEAGE # with extra category "other" for minority or non-VOC strains
GISAID_india_bystate$LINEAGE2 = GISAID_india_bystate$LINEAGE # LINEAGE2 with B.1.617 with & without E484Q
GISAID_india_bystate[GISAID_india_bystate$B617 & (!GISAID_india_bystate$E484Q),"LINEAGE2"] = "B.1.617_without_E484Q"
GISAID_india_bystate[GISAID_india_bystate$B617 & GISAID_india_bystate$E484Q,"LINEAGE2"] = "B.1.617_with_E484Q"
main_lineages = levels(GISAID_india_bystate$LINEAGE)[100*table(GISAID_india_bystate$LINEAGE)/sum(table(GISAID_india_bystate$LINEAGE)) > 1]
VOCs = c("B.1.617","B.1.618","B.1.617_with_E484Q","B.1.617_without_E484Q","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526")
main_lineages = union(main_lineages, VOCs)
GISAID_india_bystate$LINEAGE1[!(GISAID_india_bystate$LINEAGE %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_india_bystate$LINEAGE2[!(GISAID_india_bystate$LINEAGE %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_india_bystate$LINEAGE = factor(GISAID_india_bystate$LINEAGE)
GISAID_india_bystate$LINEAGE = relevel(GISAID_india_bystate$LINEAGE, ref="B.1.1.7") # we code UK strain as the reference level
GISAID_india_bystate$LINEAGE1 = factor(GISAID_india_bystate$LINEAGE1)
GISAID_india_bystate$LINEAGE1 = relevel(GISAID_india_bystate$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_india_bystate$LINEAGE1)
# "B.1.1.7"   "B.1.351"   "B.1.525"   "B.1.617"   "B.1.618"   "other"     "P.1"    
levels_LINEAGE1 = c("B.1.1.7","B.1.525","B.1.351","B.1.618","P.1","B.1.617","other")
GISAID_india_bystate$LINEAGE1 = factor(GISAID_india_bystate$LINEAGE1, levels=levels_LINEAGE1)
GISAID_india_bystate$LINEAGE2 = factor(GISAID_india_bystate$LINEAGE2)
GISAID_india_bystate$LINEAGE2 = relevel(GISAID_india_bystate$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_india_bystate$LINEAGE2)
# "B.1.1.7"               "B.1.351"               "B.1.525"               "B.1.617_with_E484Q"    "B.1.617_without_E484Q"      "B.1.618"
# "other"                 "P.1"                   
levels_LINEAGE2 = c("B.1.1.7","B.1.351","B.1.525",
                    "B.1.618","P.1","B.1.617_with_E484Q","B.1.617_without_E484Q","other")
GISAID_india_bystate$LINEAGE2 = factor(GISAID_india_bystate$LINEAGE2, levels=levels_LINEAGE2)
GISAID_india_bystate$count = 1
levels_STATES = c("West Bengal","Maharashtra","Karnataka")
GISAID_india_bystate$STATE = factor(GISAID_india_bystate$STATE, levels=levels_STATES)
table(GISAID_india_bystate$STATE)


# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of India
data_agbyweek1 = as.data.frame(table(GISAID_india_bystate$Week, GISAID_india_bystate$LINEAGE1))
colnames(data_agbyweek1) = c("Week", "LINEAGE1", "count")
data_agbyweek1_sum = aggregate(count ~ Week, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$Week, data_agbyweek1_sum$Week)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE1=="B.1.617","total"]) == nrow(GISAID_india_bystate) # correct
data_agbyweek1$Week = as.numeric(as.character(data_agbyweek1$Week))
data_agbyweek1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) + 3.5
data_agbyweek1$LINEAGE1 = factor(data_agbyweek1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

data_agbyweek2 = as.data.frame(table(GISAID_india_bystate$Week, GISAID_india_bystate$LINEAGE2))
colnames(data_agbyweek2) = c("Week", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ Week, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$Week, data_agbyweek2_sum$Week)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.617_with_E484Q","total"]) == nrow(GISAID_india_bystate) # correct
data_agbyweek2$Week = as.numeric(as.character(data_agbyweek2$Week))
data_agbyweek2$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek2$Week - 1 ) + 3.5
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total


# aggregated by week and state for selected variant lineages
data_agbyweekregion1 = as.data.frame(table(GISAID_india_bystate$Week, GISAID_india_bystate$STATE, GISAID_india_bystate$LINEAGE1))
colnames(data_agbyweekregion1) = c("Week", "STATE", "LINEAGE1", "count")
data_agbyweekregion1_sum = aggregate(count ~ Week + STATE, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$Week,data_agbyweekregion1$STATE), 
                                                                  interaction(data_agbyweekregion1_sum$Week,data_agbyweekregion1_sum$STATE))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE1=="B.1.617","total"]) == nrow(GISAID_india_bystate) # correct
data_agbyweekregion1$Week = as.numeric(as.character(data_agbyweekregion1$Week))
data_agbyweekregion1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) + 3.5
data_agbyweekregion1$LINEAGE1 = factor(data_agbyweekregion1$LINEAGE1, levels=levels_LINEAGE1)
data_agbyweekregion1$STATE = factor(data_agbyweekregion1$STATE, levels=levels_STATES)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion2 = as.data.frame(table(GISAID_india_bystate$Week, GISAID_india_bystate$STATE, GISAID_india_bystate$LINEAGE2))
colnames(data_agbyweekregion2) = c("Week", "STATE", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ Week + STATE, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$Week,data_agbyweekregion2$STATE), 
                                                                  interaction(data_agbyweekregion2_sum$Week,data_agbyweekregion2_sum$STATE))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.617_with_E484Q","total"]) == nrow(GISAID_india_bystate) # correct
data_agbyweekregion2$Week = as.numeric(as.character(data_agbyweekregion2$Week))
data_agbyweekregion2$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion2$Week - 1 ) + 3.5
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$STATE = factor(data_agbyweekregion2$STATE, levels=levels_STATES)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]


# MULLER PLOT (RAW DATA)
library(scales)
n1 = length(levels(GISAID_india_bystate$LINEAGE1))
lineage_cols1 = hcl(h = seq(15, 370, length = n1), l = 65, c = 200)
lineage_cols1[which(levels(GISAID_india_bystate$LINEAGE1)=="B.1.617")] = "magenta"
lineage_cols1[which(levels(GISAID_india_bystate$LINEAGE1)=="other")] = "grey75"

n2 = length(levels(GISAID_india_bystate$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 420, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_india_bystate$LINEAGE2)=="B.1.617_without_E484Q")] = muted("magenta")
lineage_cols2[which(levels(GISAID_india_bystate$LINEAGE2)=="B.1.617_with_E484Q")] = "magenta"
lineage_cols2[which(levels(GISAID_india_bystate$LINEAGE2)=="other")] = "grey75"

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
  labs(title = "SPREAD OF SARS-CoV2 VARIANT B.1.617 IN INDIA\n(West Bengal, Maharashtra & Karnataka)") 
  # +
  # coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANT B.1.617 IN INDIA\n(West Bengal, Maharashtra & Karnataka)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_india_raw2



muller_indiabystate_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ STATE, ncol=1) +
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANT B.1.617 IN INDIA\n(West Bengal, Maharashtra & Karnataka)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_indiabystate_raw2


# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_india_multi = nnet::multinom(LINEAGE2 ~ STATE + DATE_NUM, data=GISAID_india_bystate, maxit=1000)
fit2_india_multi = nnet::multinom(LINEAGE2 ~ STATE * DATE_NUM, data=GISAID_india_bystate, maxit=1000)
fit3_india_multi = nnet::multinom(LINEAGE2 ~ STATE + ns(DATE_NUM, df=2), data=GISAID_india_bystate, maxit=1000)
fit4_india_multi = nnet::multinom(LINEAGE2 ~ STATE * ns(DATE_NUM, df=2), data=GISAID_india_bystate, maxit=1000)
fit5_india_multi = nnet::multinom(LINEAGE2 ~ STATE + ns(DATE_NUM, df=3), data=GISAID_india_bystate, maxit=1000)
fit6_india_multi = nnet::multinom(LINEAGE2 ~ STATE * ns(DATE_NUM, df=3), data=GISAID_india_bystate, maxit=1000)
BIC(fit1_india_multi2, fit2_india_multi2, fit3_india_multi2, fit4_india_multi2, fit5_india_multi2, fit6_india_multi2) 
# fit3_india_multi2 fits best (lowest BIC)

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
max(as.Date(unlist(GISAID_india_bystate[,"Collection date"]))) # 2021-04-03
delta_r_india = data.frame(confint(emtrends(fit3_india_multi, trt.vs.ctrl ~ LINEAGE2,  
                                                          var="DATE_NUM",  mode="latent",
                                                          at=list(DATE_NUM=max(GISAID_india_bystate$DATE_NUM))), 
                                                 adjust="none", df=NA)$contrasts)
delta_r_india
# contrast     estimate         SE df   asymp.LCL    asymp.UCL
# 1               B.1.351 - B.1.1.7 -0.009593587 0.04008937 NA -0.08816731  0.068980134
# 2               B.1.525 - B.1.1.7 -0.079001200 0.04301190 NA -0.16330298  0.005300584
# 3               B.1.618 - B.1.1.7 -0.088143156 0.02071142 NA -0.12873680 -0.047549517
# 4                   P.1 - B.1.1.7 -2.606706291 0.53946112 NA -3.66403065 -1.549381933
# 5    B.1.617_with_E484Q - B.1.1.7  0.062464086 0.01566046 NA  0.03177015  0.093158024
# 6 B.1.617_without_E484Q - B.1.1.7  0.113475163 0.03026871 NA  0.05414958  0.172800747
# 7                 other - B.1.1.7 -0.073941556 0.01403207 NA -0.10144390 -0.046439209

# PS this is growth advantage of B.1.617.1+.2+.3 over B.1.1.7 :
fit3_india_multi1 = nnet::multinom(LINEAGE1 ~ STATE + ns(DATE_NUM, df=2), data=GISAID_india_bystate, maxit=1000)
delta_r_india1 = data.frame(confint(emtrends(fit3_india_multi1, trt.vs.ctrl ~ LINEAGE1,  
                                            var="DATE_NUM",  mode="latent",
                                            at=list(DATE_NUM=max(GISAID_india_bystate$DATE_NUM))), 
                                   adjust="none", df=NA)$contrasts)
delta_r_india1
# contrast     estimate         SE df   asymp.LCL    asymp.UCL
# 1 B.1.525 - B.1.1.7 -0.076883758 0.04260682 NA -0.16039160  0.006624083
# 2 B.1.351 - B.1.1.7 -0.008581747 0.03977493 NA -0.08653918  0.069375683
# 3 B.1.618 - B.1.1.7 -0.087714328 0.02070782 NA -0.12830092 -0.047127740
# 4     P.1 - B.1.1.7  0.117365007 0.15958744 NA -0.19542063  0.430150640
# 5 B.1.617 - B.1.1.  7  0.069604861 0.01543212 NA  0.03935847  0.099851256
# 6   other - B.1.1.7 -0.072765435 0.01401637 NA -0.10023701 -0.045293858



# growth advantage compared to category "other"
GISAID_india_bystate_refother = GISAID_india_bystate
GISAID_india_bystate_refother$LINEAGE1 = relevel(GISAID_india_bystate_refother$LINEAGE1, ref="other")
GISAID_india_bystate_refother$LINEAGE2 = relevel(GISAID_india_bystate_refother$LINEAGE2, ref="other")
fit3_india_multi_other = nnet::multinom(LINEAGE2 ~ STATE + ns(DATE_NUM, df=2), data=GISAID_india_bystate_refother, maxit=1000)
delta_r_india_other = data.frame(confint(emtrends(fit3_india_multi_other, trt.vs.ctrl ~ LINEAGE2,  
                                             var="DATE_NUM",  mode="latent",
                                             at=list(DATE_NUM=max(GISAID_india_bystate$DATE_NUM))), 
                                    adjust="none", df=NA)$contrasts)
delta_r_india_other
#                        contrast    estimate         SE df   asymp.LCL   asymp.UCL
# 1               B.1.1.7 - other  0.073246899 0.01404114 NA  0.04572677 0.10076703
# 2               B.1.351 - other  0.064311463 0.03950891 NA -0.01312458 0.14174751
# 3               B.1.525 - other -0.006657234 0.04283937 NA -0.09062086 0.07730639
# 4               B.1.618 - other -0.014595332 0.01879819 NA -0.05143910 0.02224843
# 5                   P.1 - other  0.156777020 0.16906537 NA -0.17458502 0.48813906
# 6    B.1.617_with_E484Q - other  0.136341607 0.01277712 NA  0.11129890 0.16138431
# 7 B.1.617_without_E484Q - other  0.189660026 0.02871247 NA  0.13338461 0.24593544

# pairwise tests for differences in growth rates among lineages
contr_india = data.frame(emtrends(fit3_india_multi, pairwise ~ LINEAGE2|1, 
                                            var="DATE_NUM",  mode="latent",
                                            at=list(DATE_NUM=max(GISAID_india_bystate$DATE_NUM)),
                                            adjust="none", df=NA)$contrasts)
contr_india

# tests for differences in growth rate with UK variant B.1.1.7
contr_india_againstB117 = data.frame(emtrends(fit3_india_multi, trt.vs.ctrl ~ LINEAGE2|1, 
                                  var="DATE_NUM",  mode="latent",
                                  at=list(DATE_NUM=max(GISAID_india_bystate$DATE_NUM)),
                                  adjust="none", df=NA)$contrasts)
contr_india_againstB117
#                          contrast     estimate         SE df   z.ratio      p.value
# 1               B.1.351 - B.1.1.7 -0.009593587 0.04008937 NA -0.239305 8.108691e-01
# 2               B.1.525 - B.1.1.7 -0.079001200 0.04301190 NA -1.836729 6.624995e-02
# 3               B.1.618 - B.1.1.7 -0.088143156 0.02071142 NA -4.255775 2.083258e-05
# 4                   P.1 - B.1.1.7 -2.606706291 0.53946112 NA -4.832056 1.351303e-06
# 5    B.1.617_with_E484Q - B.1.1.7  0.062464086 0.01566046 NA  3.988650 6.645051e-05
# 6 B.1.617_without_E484Q - B.1.1.7  0.113475163 0.03026871 NA  3.748926 1.775934e-04
# 7                 other - B.1.1.7 -0.073941556 0.01403207 NA -5.269470 1.368184e-07

contr_india_againstB117_B1617 = data.frame(emtrends(fit3_india_multi1, trt.vs.ctrl ~ LINEAGE1|1, 
                                              var="DATE_NUM",  mode="latent",
                                              at=list(DATE_NUM=max(GISAID_india_bystate$DATE_NUM)),
                                              adjust="none", df=NA)$contrasts)
contr_india_againstB117_B1617




# PLOT MULTINOMIAL FIT

extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
date.to = max(GISAID_india_bystate$DATE_NUM)+extrapolate

fit_india_multi_predsbystate = data.frame(emmeans(fit3_india_multi, 
                                            ~ LINEAGE2,
                                            by=c("DATE_NUM", "STATE"),
                                            at=list(DATE_NUM=seq(date.from, date.to)), 
                                            mode="prob", df=NA))
fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
fit_india_multi_predsbystate$LINEAGE2 = factor(fit_india_multi_predsbystate$LINEAGE2, levels=levels_LINEAGE2) 
fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$STATE, levels=levels_STATES) 

fit_india_multi_preds = data.frame(emmeans(fit3_india_multi, 
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
  annotate("rect", xmin=max(GISAID_india_bystate$DATE_NUM)+1, 
           xmax=as.Date("2021-05-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANT B.1.617 IN INDIA\n(West Bengal, Maharashtra & Karnataka, multinomial fit)")
muller_india_mfit

ggarrange(muller_india_raw2+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_india_mfit+ggtitle("Multinomial fit"), ncol=1)


muller_indiabystate_mfit = ggplot(data=fit_india_multi_predsbystate, 
                             aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india_bystate$DATE_NUM)+1, 
           xmax=as.Date("2021-05-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANT B.1.617 IN INDIA\n(West Bengal, Maharashtra & Karnataka, multinomial fit)")
muller_indiabystate_mfit

ggarrange(muller_indiabystate_raw2+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_indiabystate_mfit+ggtitle("Multinomial fit"), ncol=1)





# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

fit_india_multi_preds2 = fit_india_multi_preds
fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_with_E484Q","prob"] = fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_with_E484Q","prob"] + fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_without_E484Q","prob"]
fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_with_E484Q","asymp.LCL"] = fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_with_E484Q","asymp.LCL"] + fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_without_E484Q","asymp.LCL"]
fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_with_E484Q","asymp.UCL"] = fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_with_E484Q","asymp.UCL"] + fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2=="B.1.617_without_E484Q","asymp.UCL"]
fit_india_multi_preds2 = fit_india_multi_preds2[fit_india_multi_preds2$LINEAGE2!="B.1.617_without_E484Q",]
fit_india_multi_preds2$LINEAGE2 = as.character(fit_india_multi_preds2$LINEAGE2)
fit_india_multi_preds2$LINEAGE2[fit_india_multi_preds2$LINEAGE2=="B.1.617_with_E484Q"] = "B.1.617"
fit_india_multi_preds2$LINEAGE2 = factor(fit_india_multi_preds2$LINEAGE2, levels=levels_LINEAGE1)
fit_india_multi_preds2$LINEAGE1 = fit_india_multi_preds2$LINEAGE2
levels(fit_india_multi_preds2$LINEAGE1)


# on response scale:
plot_india_mfit = qplot(data=fit_india_multi_preds2, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE1
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE1
  ), alpha=I(1)) +
  ylab("Share among new infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANT B.1.617, B.1.17 & B.1.351 IN INDIA\n(West Bengal, Maharashtra & Karnataka)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-01")), 
                     expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01","2021-05-01")),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  geom_point(data=data_agbyweek1,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE1
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(1, 6), limits=c(1,10^3), breaks=c(10,100)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_india_mfit

ggsave(file=paste0(".\\plots\\",dat,"\\XXX.png"), width=8, height=6)


# on logit scale:

fit_india_multi_preds3 = fit_india_multi_preds2
ymin = 0.005
ymax = 0.95
fit_india_multi_preds3$asymp.LCL[fit_india_multi_preds3$asymp.LCL<ymin] = ymin
fit_india_multi_preds3$asymp.UCL[fit_india_multi_preds3$asymp.UCL<ymin] = ymin
fit_india_multi_preds3$asymp.UCL[fit_india_multi_preds3$asymp.UCL>ymax] = ymax
fit_india_multi_preds3$prob[fit_india_multi_preds3$prob<ymin] = ymin

plot_india_mfit_logit = qplot(data=fit_india_multi_preds3, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE1
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE1
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANT B.1.617, B.1.17 & B.1.351 IN INDIA\n(West Bengal, Maharashtra & Karnataka)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-01")), 
                     expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  geom_point(data=data_agbyweek1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE1
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(1, 6), limits=c(10,10^3), breaks=c(10,100)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-01")), ylim=c(0.005, 0.95), expand=c(0,0))
plot_india_mfit_logit


ggsave(file=paste0(".\\plots\\",dat,"\\XXX.png"), width=8, height=6)



# project multinomial fit onto case data
cases_india_bystate = read.csv("https://api.covid19india.org/csv/latest/states.csv")
cases_india_bystate$Date = as.Date(cases_india_bystate$Date)
cases_india_bystate = cases_india_bystate[cases_india_bystate$Date >= as.Date("2021-01-01"),]
head(cases_india_bystate)
cases_india_bystate_WB = cases_india_bystate[cases_india_bystate$State=="West Bengal",]
cases_india_bystate_M = cases_india_bystate[cases_india_bystate$State=="Maharashtra",]
cases_india_bystate_K = cases_india_bystate[cases_india_bystate$State=="Karnataka",]
cases_india_bystate_WB$newcases = c(0,diff(cases_india_bystate_WB$Confirmed))
cases_india_bystate_M$newcases = c(0,diff(cases_india_bystate_M$Confirmed))
cases_india_bystate_K$newcases = c(0,diff(cases_india_bystate_K$Confirmed))

cases_india_bystate2 = rbind(cases_india_bystate_WB, cases_india_bystate_M, cases_india_bystate_K)

fit_india_multi_predsbystate2$totnewcases = cases_india_bystate2$newcases[match(interaction(fit_india_multi_predsbystate2$STATE,fit_india_multi_predsbystate2$collection_date),
                                                                            interaction(cases_india_bystate2$State,cases_india_bystate2$Date))]
fit_india_multi_predsbystate2$cases = fit_india_multi_predsbystate2$totnewcases*fit_india_multi_predsbystate2$prob

ggplot(data=fit_india_multi_predsbystate2, 
                                  aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ STATE, ncol=1, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_india_bystate$DATE_NUM)+1, 
           xmax=as.Date("2021-05-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA\n(West Bengal, Maharashtra & Karnataka, multinomial fit)")










# IMPORT GISAID DATA FROM UK
GISAID_uk = GISAID[grepl("United Kingdom",GISAID$Location),]
nrow(GISAID_uk) # 206653
GISAID_uk$B617_1 = grepl("484Q",unlist(GISAID_uk[,"AA Substitutions"]))&(unlist(GISAID_uk[,"Pango lineage"])=="B.1.617") 
sum(GISAID_uk$B617_1) # 146
GISAID_uk$B617_2_3 = (!grepl("484Q",unlist(GISAID_uk[,"AA Substitutions"])))&(unlist(GISAID_uk[,"Pango lineage"])=="B.1.617")
sum(GISAID_uk$B617_2_3) # 122





# IMPORT GISAID DATA FROM BELGIUM
GISAID_be = GISAID[grepl("Belgium",GISAID$Location),]
nrow(GISAID_be) # 12869
GISAID_be$B617_1 = grepl("484Q",unlist(GISAID_be[,"AA Substitutions"]))&(unlist(GISAID_be[,"Pango lineage"])=="B.1.617") 
sum(GISAID_be$B617_1) # 3
GISAID_be$B617_2_3 = (!grepl("484Q",unlist(GISAID_be[,"AA Substitutions"])))&(unlist(GISAID_be[,"Pango lineage"])=="B.1.617")
sum(GISAID_be$B617_2_3) # 1
length(unique(unlist(GISAID_be[,"Pango lineage"]))) # 132
table(unlist(GISAID_be[,"Pango lineage"]))


