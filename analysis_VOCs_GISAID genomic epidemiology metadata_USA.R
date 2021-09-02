# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE USA (GISAID GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 3 AUGUST 2021

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
today = as.Date("2021-08-03")
today_num = as.numeric(today)
plotdir = "USA_GISAID_genomic_epidemiology"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-07-26_11-18.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)

GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID[GISAID$host!="Human","strain"]
GISAID = GISAID[GISAID$host=="Human",]
GISAID = GISAID[GISAID$date>=as.Date("2020-01-01"),]
range(GISAID$date) # "2020-01-01" "2021-07-21"
nrow(GISAID) #  2339555
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
library(lubridate)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
colnames(GISAID)

length(unique(GISAID$country[grepl("AY",GISAID$pango_lineage,fixed=T)|grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 103 countries
table(GISAID$pango_lineage[grepl("AY",GISAID$pango_lineage,fixed=T)|grepl("B.1.617",GISAID$pango_lineage,fixed=T)])
# AY.1      AY.2      AY.3   B.1.617 B.1.617.1 B.1.617.2 B.1.617.3 
# 371       680      2344         2      5772    241814       279 

GISAID$pango_lineage[grepl("B.1.177",GISAID$pango_lineage,fixed=T)] = "B.1.177+"
GISAID$pango_lineage[grepl("B.1.36\\>",GISAID$pango_lineage)] = "B.1.36+"

sel_target_VOC = "B.1.617"
GISAID$LINEAGE1 = GISAID$pango_lineage
GISAID$LINEAGE2 = GISAID$pango_lineage
GISAID[grepl(sel_target_VOC, GISAID$LINEAGE1, fixed=TRUE)|grepl("AY",GISAID$pango_lineage,fixed=T),"LINEAGE1"] = paste0(sel_target_VOC,"+") # in LINEAGE1 we recode B.1.617.1,2&3 & AY.1&2&3 all as B.1.617+


# ANALYSIS VOCs IN THE USA ####
sel_countries = "USA"
GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
GISAID_sel = GISAID_sel[GISAID_sel$date>=as.Date("2021-01-01"),]
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$LINEAGE1),]
GISAID_sel = GISAID_sel[GISAID_sel$country==GISAID_sel$country_exposure,] # we remove travel-related cases
range(GISAID_sel$date) # "2021-01-01" "2021-07-20"
nrow(GISAID_sel) # 501941

rowSums(table(GISAID_sel$LINEAGE1,GISAID_sel$country))
            
main_lineages = names(table(GISAID_sel$LINEAGE1))[100*table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) > 3]
main_lineages
# "B.1.1.7"  "B.1.2"    "B.1.427"  "B.1.429"  "B.1.526"  "B.1.617+" "P.1" 
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526","B.1.1.519","AY.1","AY.2","AY.3")
main_lineages = union(main_lineages, VOCs)
GISAID_sel$LINEAGE1[!(GISAID_sel$LINEAGE1 %in% main_lineages)] = "other" # minority lineages & non-VOCs
GISAID_sel$LINEAGE2[!(GISAID_sel$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove1 = names(table(GISAID_sel$LINEAGE1))[table(GISAID_sel$LINEAGE1)/sum(table(GISAID_sel$LINEAGE1)) < 0.01]
remove1 = remove1[!(remove1 %in% c("B.1.351","B.1.1.7","P.1","B.1.617+","B.1.1.519","AY.3"))]
remove2 = names(table(GISAID_sel$LINEAGE2))[table(GISAID_sel$LINEAGE2)/sum(table(GISAID_sel$LINEAGE2)) < 0.01]
remove2 = remove2[!(remove2 %in% c("B.1.351","B.1.1.7","P.1","B.1.617.2","B.1.1.519","AY.3"))]
GISAID_sel$LINEAGE1[(GISAID_sel$LINEAGE1 %in% remove1)] = "other" # minority VOCs
GISAID_sel$LINEAGE2[(GISAID_sel$LINEAGE2 %in% remove2)] = "other" # minority VOCs
table(GISAID_sel$LINEAGE1)
# B.1.1.519   B.1.1.7     B.1.2   B.1.351   B.1.427   B.1.429   B.1.526  B.1.617+     other       P.1 
# 13334    205676     60912      2449     17207     35145     46535     23760     75860     21063 
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1)
GISAID_sel$LINEAGE1 = relevel(GISAID_sel$LINEAGE1, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE1)
# "B.1.1.7"   "B.1.1.519" "B.1.2"     "B.1.351"   "B.1.427"   "B.1.429"   "B.1.526"   "B.1.617+"  "other"     "P.1"  
levels_LINEAGE1 = c("B.1.1.7","B.1.2","B.1.1.519",
                    "B.1.351","B.1.427","B.1.429","B.1.526",
                    "P.1","B.1.617+","other")
GISAID_sel$LINEAGE1 = factor(GISAID_sel$LINEAGE1, levels=levels_LINEAGE1)

GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = relevel(GISAID_sel$LINEAGE2, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE2)
# "B.1.1.7"   "AY.1"      "B.1.1.519" "B.1.2"     "B.1.351"   "B.1.427"   "B.1.429"   "B.1.526"   "B.1.617.2" "other"     "P.1"    
levels_LINEAGE2 = c("B.1.1.7","B.1.2","B.1.1.519",
                    "B.1.351","B.1.427","B.1.429","B.1.526",
                    "P.1","B.1.617.2","AY.3","other")
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2, levels=levels_LINEAGE2)

GISAID_sel$state  = GISAID_sel$division
GISAID_sel = GISAID_sel[!grepl("Princess",GISAID_sel$state),] # remove cruise boat data
GISAID_sel = GISAID_sel[!grepl("Islands|Guam",GISAID_sel$state),] # remove Northern Mariana Islands & Virgin Islands & Guam
GISAID_sel = GISAID_sel[!grepl("USA",GISAID_sel$state),] # remove data with unspecified state
GISAID_sel$state[grepl("Washington",GISAID_sel$state)] = "Washington" # Washtington DC -> Washington
sel_states <- levels_states <- c("Arkansas", "Arizona", "California", "Colorado", "Connecticut", "Florida", "Indiana", "Kansas", 
               "Massachusetts", "Missouri", "Mississippi", "Louisiana", "Nebraska", "Nevada", "New Jersey", "New York", "Oklahoma",
               "Texas", "Utah", "Washington") # 20 states with most data for B.1.617.2 or increasing cases
length(sel_states)
GISAID_sel = GISAID_sel[GISAID_sel$state %in% sel_states,]

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
GISAID_sel = GISAID_sel[-which(grepl("B.1.617", GISAID_sel$pango_lineage, fixed=TRUE)&GISAID_sel$date<=as.Date("2021-04-14")),]

GISAID_sel$state = factor(GISAID_sel$state, levels=levels_states)
levels(GISAID_sel$state)

table(GISAID_sel$LINEAGE2)
# B.1.1.7     B.1.2 B.1.1.519   B.1.351   B.1.427   B.1.429   B.1.526       P.1 B.1.617.2      AY.3     other 
# 111607     33246      7882      1030     13875     27622     30094     13542     16829      1831     50591 

range(GISAID_sel$date) # "2020-01-01" "2021-07-20"

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
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="AY.3")] = muted("magenta")
lineage_cols2[which(levels(GISAID_sel$LINEAGE2)=="other")] = "grey75"

muller_sel_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots_raw data.pdf"), width=8, height=6)


muller_usabystate_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ division) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots by state_raw data.pdf"), width=10, height=6)



# multinomial fits
data_agbyweekregion2$LINEAGE2 = relevel(data_agbyweekregion2$LINEAGE2, ref="B.1.617.2")
data_agbyweekregion2$DATE_NUM = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$STATE = data_agbyweekregion2$division

library(nnet)
library(splines)
set.seed(1)
fit1_usa_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)+STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
fit2_usa_multi = nnet::multinom(LINEAGE2 ~ scale(DATE_NUM)*STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
fit3_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)+STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
fit4_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=2)*STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
fit5_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3)+STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
fit6_usa_multi = nnet::multinom(LINEAGE2 ~ ns(DATE_NUM, df=3)*STATE, weights=count, data=data_agbyweekregion2, maxit=1000)
BIC(fit1_usa_multi, fit2_usa_multi, fit3_usa_multi, fit4_usa_multi, fit5_usa_multi, fit6_usa_multi) 
# fit4_usa_multi fits best (lowest BIC)

# growth rate advantage compared to Delta / B.1.617.2 (difference in growth rate per day) 
emtrusa = emtrends(fit4_usa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_usa = data.frame(confint(emtrusa, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrusa$contrasts)$p.value)
delta_r_usa
#                 contrast     estimate          SE df   asymp.LCL   asymp.UCL      p.value
# 1    B.1.1.7 - B.1.617.2 -0.096669504 0.002064927 NA -0.10071669 -0.09262232 4.028523e-10
# 2      B.1.2 - B.1.617.2 -0.183434898 0.005675645 NA -0.19455896 -0.17231084 4.028523e-10
# 3  B.1.1.519 - B.1.617.2 -0.157022778 0.005699977 NA -0.16819453 -0.14585103 4.028523e-10
# 4    B.1.351 - B.1.617.2 -0.147301754 0.010753115 NA -0.16837747 -0.12622604 4.028523e-10
# 5    B.1.427 - B.1.617.2 -0.147876149 0.010472476 NA -0.16840182 -0.12735047 4.028523e-10
# 6    B.1.429 - B.1.617.2 -0.190899176 0.005453771 NA -0.20158837 -0.18020998 4.028523e-10
# 7    B.1.526 - B.1.617.2 -0.113493862 0.003449500 NA -0.12025476 -0.10673297 4.028523e-10
# 8        P.1 - B.1.617.2 -0.084216954 0.002779842 NA -0.08966535 -0.07876856 4.028523e-10
# 9       AY.3 - B.1.617.2  0.005911306 0.011670138 NA -0.01696175  0.02878436 9.880398e-01
# 10     other - B.1.617.2 -0.072329362 0.003119414 NA -0.07844330 -0.06621542 4.028523e-10

# growth rate advantage compared to Delta / B.1.617.2 (difference in growth rate per day) with simple model fit1_usa_multi
emtrusa1 = emtrends(fit1_usa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_usa1 = data.frame(confint(emtrusa1, 
                                 adjust="none")$contrasts, 
                         p.value=as.data.frame(emtrusa1$contrasts, adjust="none")$p.value)
delta_r_usa1
#                 contrast    estimate           SE  df    lower.CL    upper.CL       p.value
# 1    B.1.1.7 - B.1.617.2 -0.08607682 0.0006007248 210 -0.08726104 -0.08489260 2.028816e-211
# 2      B.1.2 - B.1.617.2 -0.15120206 0.0006649525 210 -0.15251289 -0.14989122 2.939653e-253
# 3  B.1.1.519 - B.1.617.2 -0.11735341 0.0007002093 210 -0.11873375 -0.11597307 1.377341e-225
# 4    B.1.351 - B.1.617.2 -0.09032868 0.0011942114 210 -0.09268286 -0.08797451 2.517323e-154
# 5    B.1.427 - B.1.617.2 -0.12948136 0.0006811168 210 -0.13082406 -0.12813866 5.289104e-237
# 6    B.1.429 - B.1.617.2 -0.13015868 0.0006545719 210 -0.13144905 -0.12886830 4.411794e-241
# 7    B.1.526 - B.1.617.2 -0.09226210 0.0006382735 210 -0.09352034 -0.09100385 3.280951e-212
# 8        P.1 - B.1.617.2 -0.06701099 0.0006435433 210 -0.06827962 -0.06574236 1.031271e-182
# 9       AY.3 - B.1.617.2  0.02237076 0.0019395857 210  0.01854720  0.02619431  3.683296e-24
# 10     other - B.1.617.2 -0.13573484 0.0006411076 210 -0.13699867 -0.13447101 8.943659e-247

# pairwise growth rates advantages for simple multinomial model fit1_usa_multi
emtrusa_pairw1 = emtrends(fit1_usa_multi, pairwise ~ LINEAGE2,  
                    var="DATE_NUM",  mode="latent",
                    at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_pairw_usa1 = data.frame(confint(emtrusa_pairw1, 
                                  adjust="none")$contrasts, 
                          p.value=as.data.frame(emtrusa_pairw1$contrasts, adjust="none")$p.value)
delta_r_pairw_usa1[grepl("\\<B.1.1.7\\>|\\<B.1.617.2\\>|\\<B.1.351\\>|\\<P.1\\>|^AY.3\\>", delta_r_pairw_usa1$contrast)&
                   !grepl("\\<B.1.2\\>|\\<B.1.1.519\\>|\\<B.1.351\\>|\\<B.1.427\\>|^B.1.429\\>|^B.1.526\\>", delta_r_pairw_usa1$contrast)
                   ,]
#               contrast     estimate           SE  df     lower.CL     upper.CL       p.value
# 1  B.1.617.2 - B.1.1.7  0.086076821 0.0006007248 210  0.084892598  0.087261045 2.028816e-211
# 6  B.1.617.2 - B.1.429  0.130158675 0.0006545719 210  0.128868302  0.131449049 4.411794e-241
# 7  B.1.617.2 - B.1.526  0.092262095 0.0006382735 210  0.091003851  0.093520340 3.280951e-212
# 8      B.1.617.2 - P.1  0.067010991 0.0006435433 210  0.065742358  0.068279624 1.031271e-182
# 9     B.1.617.2 - AY.3 -0.022370757 0.0019395857 210 -0.026194310 -0.018547204  3.683296e-24
# 10   B.1.617.2 - other  0.135734841 0.0006411076 210  0.134471010  0.136998672 8.943659e-247
# 15   B.1.1.7 - B.1.429  0.044081854 0.0002615028 210  0.043566347  0.044597361 4.117339e-226
# 16   B.1.1.7 - B.1.526  0.006185274 0.0002413189 210  0.005709556  0.006660992  1.390999e-66
# 17       B.1.1.7 - P.1 -0.019065830 0.0003306657 210 -0.019717680 -0.018413981 1.026625e-130
# 18      B.1.1.7 - AY.3 -0.108447578 0.0019642643 210 -0.112319781 -0.104575375 5.313215e-127
# 19     B.1.1.7 - other  0.049658019 0.0002268189 210  0.049210886  0.050105153 8.124731e-250
# 53          P.1 - AY.3 -0.089381748 0.0019779236 210 -0.093280878 -0.085482618 3.736995e-110
# 54         P.1 - other  0.068723850 0.0003889606 210  0.067957082  0.069490617 2.273356e-230
# 55        AY.3 - other  0.158105598 0.0019770268 210  0.154208236  0.162002960 3.095347e-159



# fitted prop of different LINEAGES in the USA today
multinom_preds_today_avg = data.frame(emmeans(fit4_usa_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE2         prob           SE df     asymp.LCL    asymp.UCL
# 1  B.1.617.2 7.311699e-01 8.627455e-02 NA  5.620749e-01 9.002649e-01
# 2    B.1.1.7 1.475671e-02 3.482241e-03 NA  7.931643e-03 2.158178e-02
# 3      B.1.2 5.724127e-07 5.459968e-07 NA -4.977213e-07 1.642547e-06
# 4  B.1.1.519 1.981487e-05 6.433251e-06 NA  7.205935e-06 3.242381e-05
# 5    B.1.351 4.114541e-03 6.733128e-03 NA -9.082148e-03 1.731123e-02
# 6    B.1.427 1.313007e-05 8.797265e-06 NA -4.112258e-06 3.037239e-05
# 7    B.1.429 9.436433e-06 2.651604e-06 NA  4.239385e-06 1.463348e-05
# 8    B.1.526 8.475557e-04 2.652055e-04 NA  3.277625e-04 1.367349e-03
# 9        P.1 9.428002e-03 2.178736e-03 NA  5.157757e-03 1.369825e-02
# 10      AY.3 2.349791e-01 9.012996e-02 NA  5.832759e-02 4.116305e-01
# 11     other 4.661319e-03 6.839249e-04 NA  3.320851e-03 6.001788e-03

# % non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.26883014 0.06599029 0.47166999 

# fitted prop of delta by state
lineages_bystate = data.frame(emmeans(fit4_usa_multi,
                   ~ LINEAGE2,
                   by=c("DATE_NUM","STATE"),
                   at=list(DATE_NUM=today_num),  
                   mode="prob", df=NA))
lineages_bystate$DATE = as.Date(lineages_bystate$DATE_NUM, origin="1970-01-01")
delta_bystate = lineages_bystate[lineages_bystate$LINEAGE2=="B.1.617.2",]
deltaplusAY3_bystate = delta_bystate
delta_bystate = delta_bystate[order(delta_bystate$prob, decreasing=T),]
delta_bystate
# LINEAGE2 DATE_NUM         STATE      prob          SE df   asymp.LCL asymp.UCL       DATE
# 1   B.1.617.2    18837      Arkansas 0.9725283 0.006355744 NA  0.96007131 0.9849854 2021-07-29
# 177 B.1.617.2    18837      Oklahoma 0.9617354 0.036209243 NA  0.89076663 1.0327043 2021-07-29
# 56  B.1.617.2    18837       Florida 0.9615188 0.010672231 NA  0.94060161 0.9824360 2021-07-29
# 188 B.1.617.2    18837         Texas 0.9550654 0.006794461 NA  0.94174849 0.9683823 2021-07-29
# 210 B.1.617.2    18837    Washington 0.9366521 0.007997447 NA  0.92097743 0.9523268 2021-07-29
# 12  B.1.617.2    18837       Arizona 0.9121968 0.017661853 NA  0.87758021 0.9468134 2021-07-29
# 23  B.1.617.2    18837    California 0.8920910 0.016742726 NA  0.85927582 0.9249061 2021-07-29
# 34  B.1.617.2    18837      Colorado 0.8783566 0.079948557 NA  0.72166036 1.0350529 2021-07-29
# 155 B.1.617.2    18837    New Jersey 0.8714939 0.065518314 NA  0.74308036 0.9999074 2021-07-29
# 78  B.1.617.2    18837        Kansas 0.7884555 0.028982193 NA  0.73165146 0.8452596 2021-07-29
# 67  B.1.617.2    18837       Indiana 0.7639072 0.076466516 NA  0.61403563 0.9137789 2021-07-29
# 122 B.1.617.2    18837     Louisiana 0.7605992 0.100768666 NA  0.56309628 0.9581022 2021-07-29
# 100 B.1.617.2    18837      Missouri 0.6231711 0.069413129 NA  0.48712390 0.7592184 2021-07-29
# 133 B.1.617.2    18837      Nebraska 0.5921500 0.052928320 NA  0.48841238 0.6958876 2021-07-29
# 45  B.1.617.2    18837   Connecticut 0.5866333 0.233630625 NA  0.12872571 1.0445409 2021-07-29
# 166 B.1.617.2    18837      New York 0.5584603 0.102746525 NA  0.35708084 0.7598398 2021-07-29
# 199 B.1.617.2    18837          Utah 0.4890154 1.622069819 NA -2.69018307 3.6682138 2021-07-29
# 144 B.1.617.2    18837        Nevada 0.4777433 0.278182943 NA -0.06748527 1.0229718 2021-07-29
# 111 B.1.617.2    18837   Mississippi 0.3795860 0.074545092 NA  0.23348034 0.5256917 2021-07-29
# 89  B.1.617.2    18837 Massachusetts 0.2620373 0.376413619 NA -0.47571980 0.9997945 2021-07-29

AY3_bystate = lineages_bystate[lineages_bystate$LINEAGE2=="AY.3",]
deltaplusAY3_bystate$prob = deltaplusAY3_bystate$prob + AY3_bystate$prob
AY3_bystate = AY3_bystate[order(AY3_bystate$prob, decreasing=T),]
AY3_bystate
# LINEAGE2 DATE_NUM         STATE        prob          SE df    asymp.LCL  asymp.UCL       DATE
# 120     AY.3    18837   Mississippi 0.611325715 0.074754318 NA  0.464809945 0.75784149 2021-07-29
# 98      AY.3    18837 Massachusetts 0.594533640 0.573260388 NA -0.529036074 1.71810335 2021-07-29
# 153     AY.3    18837        Nevada 0.513241264 0.283411066 NA -0.042234218 1.06871675 2021-07-29
# 208     AY.3    18837          Utah 0.504525390 1.643488779 NA -2.716653426 3.72570421 2021-07-29
# 175     AY.3    18837      New York 0.410816471 0.107676179 NA  0.199775039 0.62185790 2021-07-29
# 142     AY.3    18837      Nebraska 0.389637202 0.053581391 NA  0.284619605 0.49465480 2021-07-29
# 109     AY.3    18837      Missouri 0.374474353 0.069748563 NA  0.237769681 0.51117903 2021-07-29
# 54      AY.3    18837   Connecticut 0.357711192 0.254175476 NA -0.140463587 0.85588597 2021-07-29
# 76      AY.3    18837       Indiana 0.226272666 0.077304130 NA  0.074759355 0.37778598 2021-07-29
# 131     AY.3    18837     Louisiana 0.223848427 0.101849978 NA  0.024226138 0.42347072 2021-07-29
# 87      AY.3    18837        Kansas 0.200869317 0.029091290 NA  0.143851435 0.25788720 2021-07-29
# 43      AY.3    18837      Colorado 0.106321455 0.081172299 NA -0.052773328 0.26541624 2021-07-29
# 164     AY.3    18837    New Jersey 0.095827797 0.067010260 NA -0.035509899 0.22716549 2021-07-29
# 32      AY.3    18837    California 0.043512990 0.016862789 NA  0.010462530 0.07656345 2021-07-29
# 65      AY.3    18837       Florida 0.015969891 0.010047128 NA -0.003722117 0.03566190 2021-07-29
# 10      AY.3    18837      Arkansas 0.010384319 0.004470036 NA  0.001623209 0.01914543 2021-07-29
# 197     AY.3    18837         Texas 0.006199843 0.002625712 NA  0.001053542 0.01134615 2021-07-29
# 186     AY.3    18837      Oklahoma 0.005692452 0.006932587 NA -0.007895168 0.01928007 2021-07-29
# 219     AY.3    18837    Washington 0.005550881 0.003521063 NA -0.001350276 0.01245204 2021-07-29
# 21      AY.3    18837       Arizona 0.002865894 0.005214746 NA -0.007354820 0.01308661 2021-07-29

deltaplusAY3_bystate = deltaplusAY3_bystate[order(deltaplusAY3_bystate$prob, decreasing=T),]
deltaplusAY3_bystate
# LINEAGE2 DATE_NUM         STATE      prob          SE df   asymp.LCL asymp.UCL       DATE
# 1   B.1.617.2    18837      Arkansas 0.9725283 0.006355744 NA  0.96007131 0.9849854 2021-07-29
# 12  B.1.617.2    18837       Arizona 0.9121968 0.017661853 NA  0.87758021 0.9468134 2021-07-29
# 23  B.1.617.2    18837    California 0.8920910 0.016742726 NA  0.85927582 0.9249061 2021-07-29
# 34  B.1.617.2    18837      Colorado 0.8783566 0.079948557 NA  0.72166036 1.0350529 2021-07-29
# 45  B.1.617.2    18837   Connecticut 0.5866333 0.233630625 NA  0.12872571 1.0445409 2021-07-29
# 56  B.1.617.2    18837       Florida 0.9615188 0.010672231 NA  0.94060161 0.9824360 2021-07-29
# 67  B.1.617.2    18837       Indiana 0.7639072 0.076466516 NA  0.61403563 0.9137789 2021-07-29
# 78  B.1.617.2    18837        Kansas 0.7884555 0.028982193 NA  0.73165146 0.8452596 2021-07-29
# 89  B.1.617.2    18837 Massachusetts 0.2620373 0.376413619 NA -0.47571980 0.9997945 2021-07-29
# 100 B.1.617.2    18837      Missouri 0.6231711 0.069413129 NA  0.48712390 0.7592184 2021-07-29
# 111 B.1.617.2    18837   Mississippi 0.3795860 0.074545092 NA  0.23348034 0.5256917 2021-07-29
# 122 B.1.617.2    18837     Louisiana 0.7605992 0.100768666 NA  0.56309628 0.9581022 2021-07-29
# 133 B.1.617.2    18837      Nebraska 0.5921500 0.052928320 NA  0.48841238 0.6958876 2021-07-29
# 144 B.1.617.2    18837        Nevada 0.4777433 0.278182943 NA -0.06748527 1.0229718 2021-07-29
# 155 B.1.617.2    18837    New Jersey 0.8714939 0.065518314 NA  0.74308036 0.9999074 2021-07-29
# 166 B.1.617.2    18837      New York 0.5584603 0.102746525 NA  0.35708084 0.7598398 2021-07-29
# 177 B.1.617.2    18837      Oklahoma 0.9617354 0.036209243 NA  0.89076663 1.0327043 2021-07-29
# 188 B.1.617.2    18837         Texas 0.9550654 0.006794461 NA  0.94174849 0.9683823 2021-07-29
# 199 B.1.617.2    18837          Utah 0.4890154 1.622069819 NA -2.69018307 3.6682138 2021-07-29
# 210 B.1.617.2    18837    Washington 0.9366521 0.007997447 NA  0.92097743 0.9523268 2021-07-29

# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
date.to = as.numeric(as.Date("2021-08-14")) # max(GISAID_sel$DATE_NUM)+extrapolate

# multinomial model predictions (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to), STATE=levels_states))
fit_usa_multi_preds_bystate4 = data.frame(predgrid, as.data.frame(predict(fit4_usa_multi, newdata=predgrid, type="prob")), check.names=F)
fit_usa_multi_preds_bystate1 = data.frame(predgrid, as.data.frame(predict(fit1_usa_multi, newdata=predgrid, type="prob")), check.names=F)
fit_usa_multi_preds_bystate = fit_usa_multi_preds_bystate1
library(tidyr)
fit_usa_multi_preds_bystate = gather(fit_usa_multi_preds_bystate, LINEAGE2, prob, all_of(levels_LINEAGE2), factor_key=TRUE)
fit_usa_multi_preds_bystate$collection_date = as.Date(fit_usa_multi_preds_bystate$DATE_NUM, origin="1970-01-01")
fit_usa_multi_preds_bystate$LINEAGE2 = factor(fit_usa_multi_preds_bystate$LINEAGE2, levels=levels_LINEAGE2) 
# order states from highest to lowest incidence of delta+AY.3
levels_states2 = levels_states[match(deltaplusAY3_bystate$STATE,levels_states)]
fit_usa_multi_preds_bystate$STATE = factor(fit_usa_multi_preds_bystate$STATE, levels=levels_states2)
data_agbyweekregion2$division = factor(data_agbyweekregion2$division, levels=levels_states2)
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)

muller_usa_mfit_bystate = ggplot(data=fit_usa_multi_preds_bystate, 
                                   aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ STATE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE USA\n(GISAID data, multinomial fit)")
muller_usa_mfit_bystate

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots_multinom fit by state_fit1.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots_multinom fit by state.pdf"), width=10, height=6)

# redo plot of raw data using this order
muller_usabystate_raw2 = ggplot(data=data_agbyweekregion2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ division) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots by state_raw data_fit1.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots by state_raw data.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_usabystate_raw2 + coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_usa_mfit_bystate+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots multipanel_multinom fit by state_fit1.png"), width=10, height=15)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_muller plots multipanel_multinom fit by state.pdf"), width=10, height=15)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions by state with confidence intervals (slower)
fit_usa_multi_preds_bystate_withCI4 = data.frame(emmeans(fit4_usa_multi,
                                                        ~ LINEAGE2,
                                                        by=c("DATE_NUM","STATE"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                        mode="prob", df=NA))
fit_usa_multi_preds_bystate_withCI1 = data.frame(emmeans(fit1_usa_multi,
                                                         ~ LINEAGE2,
                                                         by=c("DATE_NUM","STATE"),
                                                         at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                         mode="prob", df=NA))
fit_usa_multi_preds_bystate_withCI = fit_usa_multi_preds_bystate_withCI1
fit_usa_multi_preds_bystate_withCI$collection_date = as.Date(fit_usa_multi_preds_bystate_withCI$DATE_NUM, origin="1970-01-01")
fit_usa_multi_preds_bystate_withCI$LINEAGE2 = factor(fit_usa_multi_preds_bystate_withCI$LINEAGE2, levels=levels_LINEAGE2)
fit_usa_multi_preds_bystate_withCI$STATE = factor(fit_usa_multi_preds_bystate_withCI$STATE, levels=levels_states2)
fit_usa_multi_preds2 = fit_usa_multi_preds_bystate_withCI

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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_usa_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_logit scale_fit1.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_logit scale.pdf"), width=10, height=6)


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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",NA)),
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_response scale_fit1.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit by state_response scale.pdf"), width=10, height=6)




# overall multinomial model predictions with confidence intervals
fit_usa_multi_preds_withCI4 = data.frame(emmeans(fit4_usa_multi,
                                                        ~ LINEAGE2,
                                                        by=c("DATE_NUM"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                        mode="prob", df=NA))
fit_usa_multi_preds_withCI1 = data.frame(emmeans(fit1_usa_multi,
                                                 ~ LINEAGE2,
                                                 by=c("DATE_NUM"),
                                                 at=list(DATE_NUM=seq(date.from, date.to, by=7)),  # by=7 to speed up things a bit
                                                 mode="prob", df=NA))
fit_usa_multi_preds_withCI = fit_usa_multi_preds_withCI1
fit_usa_multi_preds_withCI$collection_date = as.Date(fit_usa_multi_preds_withCI$DATE_NUM, origin="1970-01-01")
fit_usa_multi_preds_withCI$LINEAGE2 = factor(fit_usa_multi_preds_withCI$LINEAGE2, levels=levels_LINEAGE2)
fit_usa_multi_preds3 = fit_usa_multi_preds_withCI

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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
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
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_usa_avg_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_logit scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_logit scale.pdf"), width=10, height=6)


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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
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
                        range=c(0.5, 3), limits=c(1,max(data_agbyweek2$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_usa_avg_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_response scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\usa_multinom fit_response scale.pdf"), width=10, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot = as.data.frame(get_regional_data(country = "USA", level="1"))
cases_tot = cases_tot[,c("date","state","cases_new","hosp_new","deaths_new","tested_new")]
colnames(cases_tot)[2] = "STATE"
cases_tot = cases_tot[cases_tot$date>=as.Date("2021-01-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
# cases_tot = cases_tot[cases_tot$STATE %in% levels_states2,]
cases_tot$STATE = factor(cases_tot$STATE) # , levels=levels_states2

vaccinations = read.csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv")



# smooth out weekday effects in case nrs using GAM (& potentially correct for unequal testing intensity)
library(mgcv)
k=7
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F, by=STATE) + STATE +
                  WEEKDAY, # + 
                # s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=poisson(log), data=cases_tot,
                method = "REML",
                knots = list(DATE_NUM = c(min(cases_tot$DATE_NUM)-14,
                                          seq(min(cases_tot$DATE_NUM)+1*diff(range(cases_tot$DATE_NUM))/(k-2), 
                                              max(cases_tot$DATE_NUM)-1*diff(range(cases_tot$DATE_NUM))/(k-2), length.out=k-2),
                                          max(cases_tot$DATE_NUM)+14))
) 
BIC(fit_cases)


# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_usa_multi_preds_bystate_withCI$totcases = cases_tot$cases_new[match(interaction(fit_usa_multi_preds_bystate_withCI$DATE_NUM,fit_usa_multi_preds_bystate_withCI$STATE),
                                                            interaction(cases_tot$DATE_NUM,cases_tot$STATE))]
fit_usa_multi_preds_bystate_withCI$cases = fit_usa_multi_preds_bystate_withCI$totcases * fit_usa_multi_preds_bystate_withCI$prob
fit_usa_multi_preds_bystate_withCI$cases[fit_usa_multi_preds_bystate_withCI$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM|STATE, at=list(DATE_NUM=seq(date.from, date.to, by=0.5)
), type="response"))
fit_usa_multi_preds_bystate_withCI$smoothed_totcases = cases_emmeans$rate[match(interaction(fit_usa_multi_preds_bystate_withCI$DATE_NUM,fit_usa_multi_preds_bystate_withCI$STATE),
                                                                    interaction(cases_emmeans$DATE_NUM,cases_emmeans$STATE))]
fit_usa_multi_preds_bystate_withCI$smoothed_cases = fit_usa_multi_preds_bystate_withCI$smoothed_totcases * fit_usa_multi_preds_bystate_withCI$prob
fit_usa_multi_preds_bystate_withCI$smoothed_cases[fit_usa_multi_preds_bystate_withCI$smoothed_cases<=0.001] = NA



fit_usa_multi_preds_bystate$totcases = cases_tot$cases_new[match(interaction(fit_usa_multi_preds_bystate$DATE_NUM,fit_usa_multi_preds_bystate$STATE),
                                                                        interaction(cases_tot$DATE_NUM,cases_tot$STATE))]
fit_usa_multi_preds_bystate$cases = fit_usa_multi_preds_bystate$totcases * fit_usa_multi_preds_bystate$prob
fit_usa_multi_preds_bystate$cases[fit_usa_multi_preds_bystate$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM|STATE, at=list(DATE_NUM=seq(date.from, date.to, by=0.5)
), type="response"))
fit_usa_multi_preds_bystate$smoothed_totcases = cases_emmeans$rate[match(interaction(fit_usa_multi_preds_bystate$DATE_NUM,fit_usa_multi_preds_bystate$STATE),
                                                                                interaction(cases_emmeans$DATE_NUM,cases_emmeans$STATE))]
fit_usa_multi_preds_bystate$smoothed_cases = fit_usa_multi_preds_bystate$smoothed_totcases * fit_usa_multi_preds_bystate$prob
fit_usa_multi_preds_bystate$smoothed_cases[fit_usa_multi_preds_bystate$smoothed_cases<=0.001] = NA


ggplot(data=fit_usa_multi_preds_bystate, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ STATE, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE USA\n(case data NYT & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_usa_multi_preds_bystate, 
       aes(x=collection_date, y=cases+1, group=LINEAGE2)) + 
  facet_wrap(~ STATE, scale="free") +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="AY.3",], aes(lwd=I(1.2), colour=NULL), fill=muted(I("magenta")), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="B.1.351",], aes(lwd=I(1.2), colour=NULL), fill="#9A9D00", position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE USA\n(case data NYT & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_usa_multi_preds_bystate,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) +
  facet_wrap(~ STATE, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE USA\n(case data NYT & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_usa_multi_preds_bystate, 
       aes(x=collection_date, y=smoothed_cases+1, group=LINEAGE2)) + 
  facet_wrap(~ STATE, scale="free") +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="AY.3",], aes(lwd=I(1.2), colour=NULL), fill=muted(I("magenta")), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_usa_multi_preds_bystate[fit_usa_multi_preds_bystate$LINEAGE2=="B.1.351",], aes(lwd=I(1.2), colour=NULL), fill="#9A9D00", position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE USA\n(case data NYT & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_usa_multi_preds_bystate, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ STATE, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE USA\n(case data NYT & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10()

# ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_usa_multi_preds_bystate,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) +
  facet_wrap(~ STATE, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),today), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE USA\n(case data NYT & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit smoothed case data.png"), width=8, height=6)






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
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, by="STATE", var="DATE_NUM", 
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
  facet_wrap(~ STATE) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES IN SELECTED STATES IN THE USA\n(case data NYT)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for simple multinomial model fit1_usa_multi
above_avg_r_variants1 = do.call(rbind, lapply(levels_states2, function(STATE) { do.call(rbind, 
                                                                                            lapply(seq(date.from,
                                                                                                       date.to, by=3), 
                                                                                                   function (d) { 
                                                                                                     wt = as.data.frame(emmeans(fit1_usa_multi, ~ LINEAGE2 , by="STATE", 
                                                                                                                                at=list(DATE_NUM=d, STATE=STATE), type="response"))$prob   # important: these should sum to 1
                                                                                                     # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) 
                                                                                                     # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                                                                                     cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                                                                     names(cons) = seq_along(cons)
                                                                                                     EMT = emtrends(fit1_usa_multi,  ~ LINEAGE2 , by=c("DATE_NUM", "STATE"),
                                                                                                                    var="DATE_NUM", mode="latent",
                                                                                                                    at=list(DATE_NUM=d, STATE=STATE))
                                                                                                     out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                                                                     # sum(out$estimate*wt) # should sum to zero
                                                                                                     return(out) } )) } ))
above_avg_r_variants = above_avg_r_variants1
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE2), 
                                       labels=levels_LINEAGE2)
above_avg_r_variants$LINEAGE2 = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2021-01-01" "2021-07-30"
# average growth rate of all lineages calculated from case nrs
above_avg_r_variants$avg_r = avg_r_cases$r[match(interaction(above_avg_r_variants$collection_date,above_avg_r_variants$STATE),
                                                 interaction(avg_r_cases$DATE,avg_r_cases$STATE))]  
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                country=avg_r_cases$STATE,
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
above_avg_r_variants$LINEAGE2 = factor(above_avg_r_variants$LINEAGE2, levels=c(levels_LINEAGE2,"avg"))
above_avg_r_variants$prob = fit_usa_multi_preds_bystate_withCI$prob[match(interaction(round(above_avg_r_variants$DATE_NUM),
                                                                          as.character(above_avg_r_variants$LINEAGE2),
                                                                          as.character(above_avg_r_variants$STATE)),
                                                              interaction(round(fit_usa_multi_preds_bystate_withCI$DATE_NUM),
                                                                          as.character(fit_usa_multi_preds_bystate_withCI$LINEAGE2),
                                                                          as.character(fit_usa_multi_preds_bystate_withCI$STATE)))]
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
  facet_wrap(~ STATE) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE2), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES BY VARIANT IN SELECTED AFRICAN COUNTRIES\n(case data WHO & Google)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-11-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  scale_colour_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)






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

