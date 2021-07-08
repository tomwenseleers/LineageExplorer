# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN SELECTED AFRICAN COUNTRIES (GISAID records)
# T. Wenseleers
# last update 8 JULY 2021

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
today = as.Date("2021-07-08")
today_num = as.numeric(today)
plotdir = "Africa_GISAID_records"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID metadata 
d1 = read_tsv(".//data//GISAID//Africa//gisaid_hcov-19_2021_06_28_08_africa_subm_jan_apr_2021.tsv", col_types = cols(.default = "c")) 
d2 = read_tsv(".//data//GISAID//Africa//gisaid_hcov-19_2021_07_04_15_africa_subm_may_2021.tsv", col_types = cols(.default = "c")) 
d3 = read_tsv(".//data//GISAID//Africa//gisaid_hcov-19_2021_07_08_21_africa_subm_june_july_8_2021.tsv", col_types = cols(.default = "c")) 

GISAID = as.data.frame(rbind(d1,d2,d3))
colnames(GISAID)
# [1] "Virus name"                      "Accession ID"                    "Collection date"                 "Location"                       
# [5] "Host"                            "Additional location information" "Sampling strategy"               "Gender"                         
# [9] "Patient age"                     "Patient status"                  "Last vaccinated"                 "Passage"                        
# [13] "Specimen"                        "Additional host information"     "Lineage"                         "Clade"                          
# [17] "AA Substitutions"     

GISAID = GISAID[!GISAID$`Sampling strategy` %in% c("Travel Surveillance", "S gene dropout",
                                                   "Outbreak investigation","Same-patient sampling strategy","S-gene dropout","Active surveillance",
                                                   "Quarantine testing","Contact Tracing","Longitudinal sampling on same patient(s)"),]
GISAID = GISAID[!grepl("Suspect|ontact|ravel",GISAID$`Additional host information`),] # remove travel related cases
GISAID = GISAID[!grepl("ravel",GISAID$`Additional location information`),] # remove travel related cases
GISAID = GISAID[GISAID$Host=="Human",]

library(stringr)
date_isvalid = sapply(GISAID[,"Collection date"], function (s) str_count(s, pattern = "-")==2)
sum(date_isvalid) # 19008
GISAID = GISAID[date_isvalid,]
GISAID$date = as.Date(GISAID[,"Collection date"]) 
GISAID = GISAID[!is.na(GISAID$date),]

loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID[,"Location"]), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
unique(loc[,1]) # continent
unique(loc[,2]) # countries
unique(loc[,3]) # province

GISAID$Continent = loc[,1]
GISAID$Country = loc[,2]
unique(GISAID$Country)
# [1] "South Africa"                     "Botswana"                         "Egypt"                            "Mayotte"                         
# [5] "Zimbabwe"                         "Nigeria"                          "Malawi"                           "Togo"                            
# [9] "Rwanda"                           "Reunion"                          "Cote d'Ivoire"                    "Democratic Republic of the Congo"
# [13] "Zambia"                           "Morocco"                          "Cameroon"                         "Senegal"                         
# [17] "Equatorial Guinea"                "Tunisia"                          "Mozambique"                       "Gambia"                          
# [21] "Algeria"                          "Lesotho"                          "Kenya"                            "Angola"                          
# [25] "Republic of the Congo"            "Madagascar"                       "Mali"                             "Ghana"                           
# [29] "Gabon"                            "Uganda"                           "Guinea"                           "Mauritius"                       
# [33] "Burkina Faso"                     "Comoros"                          "Ethiopia"                         "Eswatini"                        
# [37] "Canary Islands"                   "Guinea Bissau"                    "Somalia"                          "Sierra Leone"                    
# [41] "Niger"                            "Djibouti"                         "South Sudan"                      "Libya"                           
# [45] "Central African Republic"
table(GISAID[GISAID$Lineage=="B.1.617.2",]$Country)
# Angola                         Botswana Democratic Republic of the Congo                           Gambia 
# 4                              143                                6                                2 
# Ghana                            Kenya                           Malawi                          Nigeria 
# 9                               37                               14                                1 
# Reunion                           Rwanda                     South Africa                           Uganda 
# 2                               41                              612                               38 
# Zambia 
# 82 
GISAID$Region = loc[,3]
unique(GISAID$Region)
table(GISAID$Region)

GISAID = GISAID[GISAID$Lineage!="None",] # remove invalid lineages


# ANALYSIS OF VOCs IN SELECTED AFRICAN COUNTRIES ####

sel_countries = c("South Africa", "Botswana", "Zambia", "Namibia", "Malawi", "Kenya", "Uganda", "Rwanda", "Democratic Republic of the Congo", "Ghana") # "Angola","Senegal" "Botswana"
levels_countries = sel_countries

GISAID = GISAID[GISAID$date>=as.Date("2021-01-01"),]
range(GISAID$date) # "2021-01-01" "2021-06-27"

GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
library(lubridate)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
colnames(GISAID)

# conflate lineages
GISAID$Lineage[grepl("B.1.177",GISAID$Lineage,fixed=T)] = "B.1.177+"
GISAID$Lineage[grepl("B.1.36\\>",GISAID$Lineage)] = "B.1.36+"

sel_target_VOC = "B.1.617"
GISAID$LINEAGE1 = GISAID$Lineage
GISAID$LINEAGE2 = GISAID$Lineage
GISAID[grepl(sel_target_VOC, GISAID$LINEAGE1, fixed=TRUE),"LINEAGE1"] = paste0(sel_target_VOC,"+") # in LINEAGE1 we recode B.1.617.1,2&3 all as B.1.617+

# subset to selected countries
GISAID_sel = GISAID[GISAID$Country %in% sel_countries,]
nrow(GISAID_sel) # 7691

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$Country))

main_lineages = names(table(GISAID_sel$LINEAGE1))[100*table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) > 1]
main_lineages
# "A.23.1"    "B.1"       "B.1.1.318" "B.1.1.7"   "B.1.351"   "B.1.525"   "B.1.617+"
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
levels_LINEAGE2 = c("B.1.1.7",levels(GISAID_sel$LINEAGE2)[!levels(GISAID_sel$LINEAGE2) %in% c("B.1.1.7","B.1.617+","B.1.617.1","B.1.617.2","other")],
                    "B.1.617.1","B.1.617.2","other")
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2, levels=levels_LINEAGE2)

GISAID_sel$country = factor(GISAID_sel$Country, levels=levels_countries)

table(GISAID_sel$Country)


# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
GISAID_sel = GISAID_sel[-which(grepl("B.1.617", GISAID_sel$Lineage, fixed=TRUE)&GISAID_sel$date<=as.Date("2021-04-14")),]
GISAID_sel_namibia = data.frame("Virus name"=NA, # 17 B.1.617.2 Namibian samples manually added, https://twitter.com/Lizzy_Lang7/status/1412127100535029766
                                "Accession ID"=NA,
                                "Collection date"=as.Date("2021-05-30"),
                                "Location"="Africa / Namibia",
                                "Host"=NA,
                                "Additional location information"=NA,
                                "Sampling strategy"=NA,
                                "Gender"=NA,
                                "Patient age"=NA,
                                "Patient status"=NA,
                                "Last vaccinated"=NA,
                                "Passage"=NA,
                                "Specimen"=NA,
                                "Additional host information"=NA,
                                "Lineage"="B.1.617.2",
                                "Clade"=NA,
                                "AA Substitutions"=NA,
                                "date"=as.Date("2021-05-30"),
                                "Continent"="Africa",
                                "Country"="Namibia",                        
                                "Region"=NA,
                                "Week"=NA,
                                "Year"=2021,
                                "Year_Week"=NA,
                                "floor_date"=as.Date("2021-05-30"),
                                "DATE_NUM"=as.numeric(as.Date("2021-05-30")),
                                "LINEAGE1"="B.1.617+",
                                "LINEAGE2"="B.1.617.2",                       
                                "country"="Namibia", check.names=F)
GISAID_sel_namibia2 = data.frame("Virus name"=NA, # 1 B.1.351 Namibian sample manually added, https://twitter.com/Lizzy_Lang7/status/1412127100535029766
                                "Accession ID"=NA,
                                "Collection date"=as.Date("2021-05-30"),
                                "Location"="Africa / Namibia",
                                "Host"=NA,
                                "Additional location information"=NA,
                                "Sampling strategy"=NA,
                                "Gender"=NA,
                                "Patient age"=NA,
                                "Patient status"=NA,
                                "Last vaccinated"=NA,
                                "Passage"=NA,
                                "Specimen"=NA,
                                "Additional host information"=NA,
                                "Lineage"="B.1.351",
                                "Clade"=NA,
                                "AA Substitutions"=NA,
                                "date"=as.Date("2021-05-30"),
                                "Continent"="Africa",
                                "Country"="Namibia",                        
                                "Region"=NA,
                                "Week"=NA,
                                "Year"=2021,
                                "Year_Week"=NA,
                                "floor_date"=as.Date("2021-05-30"),
                                "DATE_NUM"=as.numeric(as.Date("2021-05-30")),
                                "LINEAGE1"="B.1.351",
                                "LINEAGE2"="B.1.351",                       
                                "country"="Namibia", check.names=F)
GISAID_sel = rbind(GISAID_sel,
                   GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia,
                   GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia,
                   GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia,
                   GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia, GISAID_sel_namibia,
                   GISAID_sel_namibia, GISAID_sel_namibia2)

table(GISAID_sel$LINEAGE2)

range(GISAID_sel$date) # "2021-01-01" "2021-06-27"



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

# muller plot, overall
muller_africa_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data South Africa+Botswana+Zambia+Uganda+Kenya+Malawi+\nGhana+Nigeria+Rwanda+DRC)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_africa_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plots_raw data.png"), width=8, height=6)

# muller plot, by country
muller_africa_bycountry_raw2 = ggplot(data=data_agbyweek_bycountry2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ country) +
  geom_col(aes(colour=NULL, fill=LINEAGE2), width=I(7), position="fill") +
  # geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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


# multinomial fits
data_agbyweek_bycountry2$LINEAGE2 = relevel(data_agbyweek_bycountry2$LINEAGE2, ref="B.1.1.7")
data_agbyweek_bycountry2$DATE_NUM = as.numeric(data_agbyweek_bycountry2$collection_date)

library(nnet)
library(splines)
set.seed(1)
fit1_africa_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)+country, weights=count, data=data_agbyweek_bycountry2, maxit=1000)
BIC(fit1_africa_multi) 
# 11738.57

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrafrica = emtrends(fit1_africa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_africa = data.frame(confint(emtrafrica, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrafrica$contrasts)$p.value)
delta_r_africa
# contrast     estimate          SE df     asymp.LCL   asymp.UCL      p.value
# 1    A.23.1 - B.1.1.7 -0.029682606 0.002806205 NA -0.0351826668 -0.02418254 3.116813e-10
# 2       B.1 - B.1.1.7 -0.017290403 0.001954987 NA -0.0211221062 -0.01345870 3.123779e-10
# 3 B.1.1.318 - B.1.1.7  0.043569881 0.004069417 NA  0.0355939708  0.05154579 3.116787e-10
# 4   B.1.351 - B.1.1.7 -0.022961511 0.001500486 NA -0.0259024092 -0.02002061 3.116536e-10
# 5   B.1.525 - B.1.1.7  0.004588359 0.002788310 NA -0.0008766279  0.01005335 4.413733e-01
# 6 B.1.617.1 - B.1.1.7  0.036089815 0.009688977 NA  0.0170997697  0.05507986 2.549492e-03
# 7 B.1.617.2 - B.1.1.7  0.082207480 0.004087268 NA  0.0741965826  0.09021838 3.116536e-10
# 8     other - B.1.1.7 -0.024712586 0.001830991 NA -0.0283012624 -0.02112391 3.116540e-10


# fitted prop of different LINEAGES in the africa today
multinom_preds_today_avg = data.frame(emmeans(fit1_africa_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE2        prob           SE df     asymp.LCL   asymp.UCL
# 1   B.1.1.7 0.008421286 0.0018579564 NA  0.0047797585 0.012062814
# 2    A.23.1 0.001045963 0.0004335694 NA  0.0001961824 0.001895743
# 3       B.1 0.001286566 0.0002919790 NA  0.0007142979 0.001858835
# 4 B.1.1.318 0.056542496 0.0113195641 NA  0.0343565580 0.078728434
# 5   B.1.351 0.010441396 0.0014117131 NA  0.0076744890 0.013208302
# 6   B.1.525 0.003671328 0.0012003110 NA  0.0013187618 0.006023894
# 7 B.1.617.1 0.002997387 0.0020830367 NA -0.0010852897 0.007080064
# 8 B.1.617.2 0.914131630 0.0133441025 NA  0.8879776701 0.940285591
# 9     other 0.001461948 0.0003504848 NA  0.0007750100 0.002148885

# % non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.9915787 0.9319277 1.0512297 


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data, multinomial fit)")
muller_africa_bycountry_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_by country_muller plots_multinom fit.png"), width=10, height=6)


library(ggpubr)
ggarrange(muller_africa_bycountry_raw2 + coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_africa_bycountry_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_by country_muller plots multipanel_multinom fit.png"), width=10, height=10)


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
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data South Africa+Botswana+Zambia+Namibia+\nMalawi+Kenya+Uganda+Rwanda+DRC+Ghana,\nmultinomial fit)")
muller_africa_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plot_multinom fit avg.png"), width=10, height=6)


library(ggpubr)
ggarrange(muller_africa_raw2 + coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_africa_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_muller plot multipanel_multinom fit avg.png"), width=10, height=10)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# overall average multinomial model predictions over all selected countries with confidence intervals

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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data South Africa+Botswana+Zambia+Namibia+\nMalawi+Kenya+Uganda+Rwanda+DRC+Ghana,\nmultinomial fit)") +
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
plot_africa_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_multinom fit avg_logit scale.png"), width=10, height=6)


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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data South Africa+Botswana+Zambia+Namibia+\nMalawi+Kenya+Uganda+Rwanda+DRC+Ghana, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
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
                        range=c(0.5, 4), limits=c(1,max(data_agbyweek2$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_africa_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa_multinom fit avg_response scale.png"), width=10, height=6)




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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_africa_bycountry_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\africa by country_multinom fit_logit scale.png"), width=10, height=6)


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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data, multinomial fit)") +
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





# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot_DRC = as.data.frame(get_national_data(countries = levels_countries[grepl("Congo",levels_countries)], source="WHO"))
cases_tot_DRC = cases_tot_DRC[,c("date","country","cases_new","hosp_new","deaths_new","tested_new")]
cases_tot = as.data.frame(get_national_data(countries = levels_countries, source="Google"))
cases_tot = cases_tot[,c("date","country","cases_new","hosp_new","deaths_new","tested_new")]
cases_tot = rbind(cases_tot_DRC, cases_tot)
cases_tot = cases_tot[cases_tot$date>=as.Date("2021-01-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
cases_tot$country = factor(cases_tot$country, levels=unique(cases_tot$country),
                              labels=gsub("Congo - Kinshasa", "Democratic Republic of the Congo", unique(cases_tot$country)))
cases_tot$country = factor(cases_tot$country, levels=levels_countries)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
library(mgcv)
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=7, fx=F, by=country) + country +
                  WEEKDAY # +
                # s(tested_new, bs="cs", k=8, fx=F, by=region)
                ,
                family=poisson(log), data=cases_tot,
) 
BIC(fit_cases)

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_africa_multi_preds$totcases = cases_tot$cases_new[match(interaction(fit_africa_multi_preds$DATE_NUM,fit_africa_multi_preds$country),
                                                                     interaction(cases_tot$DATE_NUM,cases_tot$country))]
fit_africa_multi_preds$cases = fit_africa_multi_preds$totcases * fit_africa_multi_preds$prob
fit_africa_multi_preds$cases[fit_africa_multi_preds$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM|country, at=list(DATE_NUM=seq(date.from, date.to, by=0.5)
), type="response"))
fit_africa_multi_preds$smoothed_totcases = cases_emmeans$rate[match(interaction(fit_africa_multi_preds$DATE_NUM,fit_africa_multi_preds$country),
                                                                             interaction(cases_emmeans$DATE_NUM,cases_emmeans$country))]
fit_africa_multi_preds$smoothed_cases = fit_africa_multi_preds$smoothed_totcases * fit_africa_multi_preds$prob
fit_africa_multi_preds$smoothed_cases[fit_africa_multi_preds$smoothed_cases<=0.001] = NA

ggplot(data=fit_africa_multi_preds, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ country, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_africa_multi_preds, 
       aes(x=collection_date, y=cases+1, group=LINEAGE2)) + 
  facet_wrap(~ country, scale="free") +
  geom_col(data=fit_africa_multi_preds[fit_africa_multi_preds$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_col(data=fit_africa_multi_preds[fit_africa_multi_preds$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  geom_col(data=fit_africa_multi_preds[fit_africa_multi_preds$LINEAGE2=="B.1.351",], aes(lwd=I(1.2), colour=NULL), fill="#9A9D00", position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_africa_multi_preds,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) +
  facet_wrap(~ country, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_africa_multi_preds, 
       aes(x=collection_date, y=smoothed_cases+1, group=LINEAGE2)) + 
  facet_wrap(~ country, scale="free") +
  geom_area(data=fit_africa_multi_preds[fit_africa_multi_preds$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_africa_multi_preds[fit_africa_multi_preds$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_africa_multi_preds[fit_africa_multi_preds$LINEAGE2=="B.1.351",], aes(lwd=I(1.2), colour=NULL), fill="#9A9D00", position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_africa_multi_preds, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_africa_multi_preds,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) +
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),today), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit smoothed case data.png"), width=8, height=6)






# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####
# TO DO : need to finish this part

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
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, by="region", var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to, by=3)
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[3] = "r"
colnames(avg_r_cases)[6] = "r_LOWER"
colnames(avg_r_cases)[7] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  facet_wrap(~ region) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN THE UK AT MOMENT OF INFECTION BASED ON NEW CASES\n(data gov.uk)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants4 = do.call(rbind, lapply(levels_REGION, function(region) { do.call(rbind, 
                                                                                        lapply(seq(date.from,
                                                                                                   date.to, by=3), 
                                                                                               function (d) { 
                                                                                                 wt = as.data.frame(emmeans(fit4_cogukp2_multi, ~ LINEAGE2 , by="REGION", 
                                                                                                                            at=list(DATE_NUM=d, REGION=region), type="response"))$prob   # important: these should sum to 1
                                                                                                 # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) 
                                                                                                 # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                                                                                 cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                                                                 names(cons) = seq_along(cons)
                                                                                                 EMT = emtrends(fit4_cogukp2_multi,  ~ LINEAGE2 , by=c("DATE_NUM", "REGION"),
                                                                                                                var="DATE_NUM", mode="latent",
                                                                                                                at=list(DATE_NUM=d, REGION=region))
                                                                                                 out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                                                                 # sum(out$estimate*wt) # should sum to zero
                                                                                                 return(out) } )) } ))
above_avg_r_variants = above_avg_r_variants4
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE2), 
                                       labels=levels(data_agbyweekregion1$LINEAGE2))
above_avg_r_variants$LINEAGE2 = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2020-09-01" "2021-07-31"
# average growth rate of all lineages calculated from case nrs
above_avg_r_variants$avg_r = avg_r_cases$r[match(interaction(above_avg_r_variants$collection_date,above_avg_r_variants$REGION),
                                                 interaction(avg_r_cases$DATE,avg_r_cases$region))]  
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                REGION=avg_r_cases$region,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE,
                LINEAGE2="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$LINEAGE2 = factor(above_avg_r_variants$LINEAGE2, levels=c(levels_LINEAGE2_plot,"avg"))
above_avg_r_variants$prob = fit_africa_multi_preds$prob[match(interaction(round(above_avg_r_variants$DATE_NUM),
                                                                                   as.character(above_avg_r_variants$LINEAGE2),
                                                                                   as.character(above_avg_r_variants$REGION)),
                                                                       interaction(round(fit_africa_multi_preds$DATE_NUM),
                                                                                   as.character(fit_africa_multi_preds$LINEAGE2),
                                                                                   as.character(fit_africa_multi_preds$country)))]
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
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$LINEAGE2 %in% c("other"))|(above_avg_r_variants2$collection_date>=max(cases_tot$date))),], # |above_avg_r_variants2$collection_date>max(cases_tot$DATE)
      x=collection_date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE2, fill=LINEAGE2, alpha=I(0.5),
      group=LINEAGE2, linetype=I(0)) +
  facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE2), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN THE UK\nAT MOMENT OF INFECTION\n(based on gov.uk case data & multinomial fit to\nCOG-UK lineage frequencies)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-11-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  scale_colour_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)


