# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN SELECTED AFRICAN COUNTRIES (GISAID records)
# T. Wenseleers
# last update 4 JULY 2021

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
today = as.Date("2021-07-04")
today_num = as.numeric(today)
plotdir = "Africa_GISAID_records"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID metadata 
d1 = read_tsv(".//data//GISAID//Africa//gisaid_hcov-19_2021_06_28_08_africa_subm_jan_apr_2021.tsv", col_types = cols(.default = "c")) 
d2 = read_tsv(".//data//GISAID//Africa//gisaid_hcov-19_2021_07_04_15_africa_subm_may_2021.tsv", col_types = cols(.default = "c")) 
d3 = read_tsv(".//data//GISAID//Africa//gisaid_hcov-19_2021_07_04_15_africa_subm_june_july_4_2021.tsv", col_types = cols(.default = "c")) 

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
# [41] "Niger"                            "Djibouti"                         "South Sudan"                      "Central African Republic" 
table(GISAID[GISAID$Lineage=="B.1.617.2",]$Country)
# Angola                         Botswana Democratic Republic of the Congo                            Ghana 
# 4                               14                                6                                2 
# Kenya                           Malawi                          Reunion                     South Africa 
# 37                               14                                2                              448 
# Uganda                           Zambia 
# 38                               82 
GISAID$Region = loc[,3]
unique(GISAID$Region)
table(GISAID$Region)

GISAID = GISAID[GISAID$Lineage!="None",] # remove invalid lineages


# ANALYSIS OF VOCs IN SELECTED AFRICAN COUNTRIES ####

sel_countries = c("South Africa", "Zambia", "Uganda", "Botswana", "Kenya", "Malawi", "Democratic Republic of the Congo") # "Angola","Senegal" "Botswana"
levels_countries = sel_countries

GISAID = GISAID[GISAID$date>=as.Date("2021-01-01"),]
range(GISAID$date) # "2021-01-01" "2021-06-23"

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
nrow(GISAID_sel) # 6392

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$Country))

main_lineages = names(table(GISAID_sel$LINEAGE1))[100*table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) > 2]
main_lineages
# "B.1"      "B.1.1.7"  "B.1.351"  "B.1.617+"
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

table(GISAID_sel$Country)

GISAID_sel$country = factor(GISAID_sel$Country, levels=levels_countries)

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
GISAID_sel = GISAID_sel[-which(grepl("B.1.617", GISAID_sel$Lineage, fixed=TRUE)&GISAID_sel$date<=as.Date("2021-04-14")),]

table(GISAID_sel$LINEAGE2)

range(GISAID_sel$date) # "2021-01-01" "2021-06-23"



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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi)") 
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
fit2_africa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)+country, weights=count, data=data_agbyweek_bycountry2, maxit=1000)
BIC(fit1_africa_multi, fit2_africa_multi) 
# df      BIC
# fit1_africa_multi 48 8398.625
# fit2_africa_multi 54 8338.649

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrafrica = emtrends(fit1_africa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_africa = data.frame(confint(emtrafrica, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrafrica$contrasts)$p.value)
delta_r_africa
# contrast     estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1       B.1 - B.1.1.7 -0.022887137 0.002285126 NA -0.027365900 -0.018408373 2.405853e-13
# 2   B.1.351 - B.1.1.7 -0.028240031 0.001807628 NA -0.031782916 -0.024697146 0.000000e+00
# 3   B.1.525 - B.1.1.7 -0.004277424 0.004029495 NA -0.012175090  0.003620242 7.412186e-01
# 4 B.1.617.1 - B.1.1.7  0.036923291 0.014001860 NA  0.009480151  0.064366432 5.506474e-02
# 5 B.1.617.2 - B.1.1.7  0.071212756 0.004548427 NA  0.062298004  0.080127508 0.000000e+00
# 6     other - B.1.1.7 -0.031605417 0.002135395 NA -0.035790714 -0.027420120 0.000000e+00


# fitted prop of different LINEAGES in the africa today
multinom_preds_today_avg = data.frame(emmeans(fit1_africa_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE2        prob          SE df     asymp.LCL   asymp.UCL
# 1   B.1.1.7 0.018137588 0.004624153 NA  0.0090744147 0.027200760
# 2       B.1 0.004076623 0.001252059 NA  0.0016226326 0.006530613
# 3   B.1.351 0.037862791 0.007636248 NA  0.0228960199 0.052829563
# 4   B.1.525 0.003085489 0.001284387 NA  0.0005681375 0.005602841
# 5 B.1.617.1 0.007436145 0.006251546 NA -0.0048166599 0.019688949
# 6 B.1.617.2 0.925685641 0.014551780 NA  0.8971646754 0.954206607
# 7     other 0.003715723 0.001058946 NA  0.0016402278 0.005791218

# % non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.9818624 0.9190750 1.0446498 


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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)")
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)")
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)") +
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)") +
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN AFRICA\n(GISAID data DRC+Uganda+Kenya+Malawi, multinomial fit)") +
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

