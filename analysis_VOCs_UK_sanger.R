# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN ENGLAND BASED ON SANGER INSTITUTE BASELINE SURVEILLANCE SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & most surge testing/active surveillance)

# last update 15 JUNE 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-16")
today_num = as.numeric(today)
today # "2021-06-11"
plotdir = "VOCs_SANGER"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))


sanger = read_tsv("https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv")
sanger = as.data.frame(sanger)
# code ONS regions
LTLAs_regions = read.csv(".//data/ONS/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")
sanger$REGION = LTLAs_regions$RGN20NM[match(sanger$LTLA, LTLAs_regions$LAD20CD)]
sanger = sanger[!is.na(sanger$REGION),]
sanger$REGION = factor(sanger$REGION)
levels(sanger$REGION)
levels_REGION = c("London","North West","South West","South East","East of England","East Midlands","West Midlands",
                  "North East","Yorkshire and The Humber")
sanger$REGION = factor(sanger$REGION, levels=levels_REGION)
head(sanger)
# code NHS England regions
LTLAs_NHSregions = read.csv(".//data/ONS/LTLA_to_NHSregion_lookup.csv")
sanger$NHSREGION = LTLAs_NHSregions$NHSREGION[match(sanger$LTLA, LTLAs_NHSregions$LTLA)]
sanger = sanger[!is.na(sanger$REGION),]
sanger$NHSREGION = factor(sanger$NHSREGION)
levels(sanger$NHSREGION)
levels_NHSREGION = c("London","North West","South West","South East","East of England","Midlands",
                     "North East and Yorkshire")
sanger$NHSREGION = factor(sanger$NHSREGION, levels=levels_NHSREGION)
head(sanger)


# sanger = sanger[grepl("2021-", sanger[,"WeekEndDate"]),]
sanger = sanger[!(sanger$Lineage=="None"|sanger$Lineage=="Lineage data suppressed"),]
range(sanger$WeekEndDate) # "2021-09-05" "2021-06-05"
# sanger$Week = lubridate::week(sanger$WeekEndDate)
sanger$DATE_NUM = as.numeric(sanger$WeekEndDate)-3.5 # using week midpoint
colnames(sanger)

sanger = sanger[rep(seq_len(nrow(sanger)), sanger$Count),] # convert to long format
sanger$Count = NULL
nrow(sanger) # 241356

nrow(sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 27378 (last 2 weeks)
nrow(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 23055 (last 2 weeks)
nrow(sanger[grepl("B.1.1.7",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 4223

length(unique(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),"LTLA"])) # 299 LTLAs

unique(sanger$Lineage)
sel_target_VOC = "B.1.617"
table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
sanger[sanger$Lineage=="B.1.617",]
sum(table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])) # 27493 B.1.617

unique(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
# "B.1.617"   "B.1.617.1" "B.1.617.3" "B.1.617.2"
sum(grepl("B.1.617",sanger$Lineage, fixed=TRUE)) # 27493 B.1.617+
sanger$LINEAGE1 = sanger$Lineage
sanger$LINEAGE2 = sanger$Lineage
sanger[grepl(sel_target_VOC, sanger$LINEAGE1, fixed=T),"LINEAGE1"] = paste0(sel_target_VOC,"+")
sel_target_VOC = paste0(sel_target_VOC, "+")
sanger[grepl("B.1.177", sanger$LINEAGE1, fixed=T),"LINEAGE1"] = "B.1.177+"
sanger[grepl("B.1.177", sanger$LINEAGE2, fixed=T),"LINEAGE2"] = "B.1.177+"

sum("B.1.1.7"==sanger$LINEAGE1) # 152892 B.1.1.7


table_lineage = as.data.frame(table(sanger$LINEAGE1))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger$LINEAGE1))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# Lineage        Prop
# 13        B.1 0.003712359
# 15    B.1.1.1 0.003604634
# 110  B.1.1.29 0.001379705
# 121 B.1.1.307 0.005485673
# 126 B.1.1.311 0.006355757
# 129 B.1.1.315 0.009550208
# 133  B.1.1.37 0.009285040
# 148   B.1.1.7 0.633470889
# 160   B.1.160 0.005191501
# 163  B.1.177+ 0.160861963
# 179   B.1.235 0.001649016
# 184   B.1.258 0.005141782
# 190 B.1.258.3 0.001027528
# 202   B.1.351 0.001201545
# 203    B.1.36 0.002100631
# 208 B.1.36.17 0.004938763
# 222   B.1.389 0.001847893
# 263  B.1.617+ 0.113910572

# just for the last week worth of data
table_lineage = as.data.frame(table(sanger[sanger$WeekEndDate==max(sanger$WeekEndDate),]$LINEAGE2))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger[sanger$WeekEndDate==max(sanger$WeekEndDate),]$LINEAGE2))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
#      Lineage        Prop
# 7    B.1.1.7 0.07432432
# 11 B.1.617.2 0.92347769


sel_ref_lineage = "B.1.1.7"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"Lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.1.7","B.1.177+","B.1.351", "B.1.525","B.1.617+","B.1.617.1","B.1.617.2","P.1")

sanger$LINEAGE1[!(sanger$LINEAGE1 %in% sel_lineages)] = "other"
sanger$LINEAGE2[!(sanger$LINEAGE2 %in% sel_lineages)] = "other"
# sanger = sanger[sanger$LINEAGE1 %in% sel_lineages, ]
sum(table(sanger$LINEAGE1))
table(sanger$LINEAGE1)
# B.1.1.7 B.1.177+  B.1.351  B.1.525 B.1.617+    other      P.1 
# 152892    38825      290      142    27493    21649       65
table(sanger$LINEAGE2)
# B.1.1.7  B.1.177+   B.1.351   B.1.525 B.1.617.1 B.1.617.2     other       P.1 
# 152892     38825       290       142       105     27341     21696        65 

levels_LINEAGE1 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617+","P.1")
levels_LINEAGE2 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617.1","B.1.617.2","P.1")
sanger$LINEAGE1 = factor(sanger$LINEAGE1, levels=levels_LINEAGE1)
sanger$LINEAGE2 = factor(sanger$LINEAGE2, levels=levels_LINEAGE2)

sanger_lastmonth = sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-30),]
table(sanger_lastmonth$Lineage, sanger_lastmonth$REGION)

str(sanger)

# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of England
data_agbyweek1 = as.data.frame(table(sanger$WeekEndDate, sanger$LINEAGE2))
colnames(data_agbyweek1) = c("WeekEndDate", "LINEAGE2", "count")
data_agbyweek1_sum = aggregate(count ~ WeekEndDate, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$WeekEndDate, data_agbyweek1_sum$WeekEndDate)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE2=="B.1.617.1","total"]) == nrow(sanger) # TRUE
data_agbyweek1$WeekEndDate = as.Date(as.character(data_agbyweek1$WeekEndDate))
data_agbyweek1$collection_date = data_agbyweek1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and ONS region
data_agbyweekregion1 = as.data.frame(table(sanger$WeekEndDate, sanger$REGION, sanger$LINEAGE2))
colnames(data_agbyweekregion1) = c("WeekEndDate", "REGION", "LINEAGE2", "count")
data_agbyweekregion1_sum = aggregate(count ~ WeekEndDate + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$WeekEndDate,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$WeekEndDate,data_agbyweekregion1_sum$REGION))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE2=="B.1.617.1","total"]) == nrow(sanger) # TRUE
data_agbyweekregion1$WeekEndDate = as.Date(as.character(data_agbyweekregion1$WeekEndDate))
data_agbyweekregion1$collection_date = data_agbyweekregion1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion1[data_agbyweekregion1$collection_date==max(data_agbyweekregion1$collection_date),]

# aggregated by week and NHS region
data_agbyweeknhsregion1 = as.data.frame(table(sanger$WeekEndDate, sanger$NHSREGION, sanger$LINEAGE2))
colnames(data_agbyweeknhsregion1) = c("WeekEndDate", "NHSREGION", "LINEAGE2", "count")
data_agbyweeknhsregion1_sum = aggregate(count ~ WeekEndDate + NHSREGION, data=data_agbyweeknhsregion1, sum)
data_agbyweeknhsregion1$total = data_agbyweeknhsregion1_sum$count[match(interaction(data_agbyweeknhsregion1$WeekEndDate,data_agbyweeknhsregion1$NHSREGION), 
                                                                  interaction(data_agbyweeknhsregion1_sum$WeekEndDate,data_agbyweeknhsregion1_sum$NHSREGION))]
sum(data_agbyweeknhsregion1[data_agbyweeknhsregion1$LINEAGE2=="B.1.617.1","total"]) == nrow(sanger) # TRUE
data_agbyweeknhsregion1$WeekEndDate = as.Date(as.character(data_agbyweeknhsregion1$WeekEndDate))
data_agbyweeknhsregion1$collection_date = data_agbyweeknhsregion1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweeknhsregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweeknhsregion1$LINEAGE2 = factor(data_agbyweeknhsregion1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweeknhsregion1$NHSREGION = factor(data_agbyweeknhsregion1$NHSREGION, levels=levels_NHSREGION)
data_agbyweeknhsregion1$collection_date_num = as.numeric(data_agbyweeknhsregion1$collection_date)
data_agbyweeknhsregion1$prop = data_agbyweeknhsregion1$count/data_agbyweeknhsregion1$total
data_agbyweeknhsregion1 = data_agbyweeknhsregion1[data_agbyweeknhsregion1$total!=0,]

data_agbyweeknhsregion1[data_agbyweeknhsregion1$collection_date==max(data_agbyweeknhsregion1$collection_date),]



# MULLER PLOT (RAW DATA)
unique(sanger$LINEAGE2)
# levels_LINEAGE2_plot = rev(c("B.1.1.7","B.1.617.2","B.1.617.1","P.1","B.1.351","B.1.525","B.1.177+","other")) # "B.1.617.1","B.1.617.2","B.1.617.3"
levels_LINEAGE2_plot = c("other","B.1.177+","B.1.525","B.1.351","P.1","B.1.1.7","B.1.617.1","B.1.617.2")

library(scales)
n1 = length(levels_LINEAGE2_plot)
lineage_cols1 = c(hcl(h = seq(0, 260, length = n1-1), l = 60, c = 200), muted("dodgerblue", l=55, c=120)) # grey70
# lineage_cols1 = c(hcl(h = seq(0, 265, length = n1), l = 60, c = 200)) # grey70
# lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617+")] = "magenta" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.1")] = muted("magenta") # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=50,c=100)

# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.1")] = "magenta" # muted("magenta",l=50,c=100)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=80,c=255)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.3")] = "magenta" # muted("magenta",l=100,c=200)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")] = "#0085FF"  
lineage_cols1[which(levels_LINEAGE2_plot=="other")] = "grey70"  
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.177+")] = "grey55"  
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.351")] = muted("cyan")  
lineage_cols1[which(levels_LINEAGE2_plot=="P.1")] = "cyan3"  

data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2_plot)
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2_plot)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweeknhsregion1$LINEAGE2 = factor(data_agbyweeknhsregion1$LINEAGE2, levels=levels_LINEAGE2_plot)
data_agbyweeknhsregion1$NHSREGION = factor(data_agbyweeknhsregion1$NHSREGION, levels=levels_NHSREGION)


library(ggplot2)
library(ggthemes)
muller_sanger_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + #  
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN ENGLAND\n(Sanger Institute baseline surveillance data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sanger_raw1


muller_sangerbyregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + #  
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN ENGLAND\n(Sanger Institute baseline surveillance data)") +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbyregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_raw data.pdf"), width=8, height=6)


muller_sangerbynhsregion_raw1 = ggplot(data=data_agbyweeknhsregion1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ NHSREGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + #  
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN ENGLAND\n(Sanger Institute baseline surveillance data)") +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbynhsregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
sanger$LINEAGE1 = relevel(sanger$LINEAGE1, ref="B.1.1.7") # we take B.1.1.7 as baseline / reference level
sanger$LINEAGE2 = relevel(sanger$LINEAGE2, ref="B.1.1.7")
# sanger$DATE_NUM_WEEK = sanger$DATE_NUM/7

data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
data_agbyweekregion1$LINEAGE2 = relevel(data_agbyweekregion1$LINEAGE2, ref="B.1.1.7")
data_agbyweeknhsregion1$DATE_NUM = data_agbyweeknhsregion1$collection_date_num
data_agbyweeknhsregion1$LINEAGE2 = relevel(data_agbyweeknhsregion1$LINEAGE2, ref="B.1.1.7")


# by ONS region
set.seed(1)
fit1_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION + DATE_NUM, data=data_agbyweekregion1, weights=count, maxit=1000)
fit2_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * DATE_NUM, data=data_agbyweekregion1, weights=count, maxit=1000)
fit3_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION + ns(DATE_NUM, df=2), data=data_agbyweekregion1, weights=count, maxit=1000)
fit4_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=data_agbyweekregion1, weights=count, maxit=1000)
BIC(fit1_sanger_multi, fit2_sanger_multi, fit3_sanger_multi, fit4_sanger_multi) 
# fit4_sanger_multi fits best (lowest BIC) but I will use the slightly simpler fit3_sanger_multi below

# by NHS region
set.seed(1)
fit1_sanger_nhs_multi = nnet::multinom(LINEAGE2 ~ NHSREGION + DATE_NUM, data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit2_sanger_nhs_multi = nnet::multinom(LINEAGE2 ~ NHSREGION * DATE_NUM, data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit3_sanger_nhs_multi = nnet::multinom(LINEAGE2 ~ NHSREGION + ns(DATE_NUM, df=2), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit4_sanger_nhs_multi = nnet::multinom(LINEAGE2 ~ NHSREGION * ns(DATE_NUM, df=2), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
BIC(fit1_sanger_nhs_multi, fit2_sanger_nhs_multi, fit3_sanger_nhs_multi, fit4_sanger_nhs_multi) 
# fit4_sanger_multi fits best (lowest BIC) but I will use the slightly simpler fit3_sanger_nhs_multi below


# growth rate advantages of different VOCs compared to UK type B.1.1.7 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrsanger = emtrends(fit3_sanger_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger = data.frame(confint(emtrsanger, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger$contrasts)$p.value)
delta_r_sanger
#               contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1      other - B.1.1.7  0.01668160 0.001008907 NA  0.0147041741  0.01865902 0.0000000000
# 2 (B.1.177+) - B.1.1.7 -0.04938597 0.001623129 NA -0.0525672449 -0.04620470 0.0000000000
# 3    B.1.525 - B.1.1.7  0.00952153 0.004523594 NA  0.0006554494  0.01838761 0.1880617462
# 4    B.1.351 - B.1.1.7  0.01400645 0.003059078 NA  0.0080107718  0.02000214 0.0001199795
# 5        P.1 - B.1.1.7  0.01989233 0.008813869 NA  0.0026174637  0.03716720 0.1376627109
# 6  B.1.617.1 - B.1.1.7 -0.07906623 0.007383808 NA -0.0935382312 -0.06459423 0.0000000000
# 7  B.1.617.2 - B.1.1.7  0.11757769 0.001845234 NA  0.1139611026  0.12119429 0.0000000000

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 5.5 days, Ferretti et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*5.5)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
# contrast  estimate  asymp.LCL   asymp.UCL
# 1      other - B.1.1.7   9.608942   8.4233144  10.80754
# 2 (B.1.177+) - B.1.1.7 -31.209203 -33.5251745 -28.93340
# 3    B.1.525 - B.1.1.7   5.376389   0.3611478  10.64225
# 4    B.1.351 - B.1.1.7   8.008041   4.5044267  11.62912
# 5        P.1 - B.1.1.7  11.561723   1.4500173  22.68128
# 6  B.1.617.1 - B.1.1.7 -54.475339 -67.2735436 -42.65633
# 7  B.1.617.2 - B.1.1.7  90.918666  87.1585747  94.75430

# so this would estimate that B.1.617.2 was 91% more infectious than B.1.1.7 [87%-95%] 95% CLs

# or with generation time of 4.7 days (Nishiura et al. 2020)
transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*4.7)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
# contrast  estimate  asymp.LCL   asymp.UCL
# 1      other - B.1.1.7   8.155898   7.1553664   9.165771
# 2 (B.1.177+) - B.1.1.7 -26.126358 -28.0263673 -24.254546
# 3    B.1.525 - B.1.1.7   4.576763   0.3085362   9.026607
# 4    B.1.351 - B.1.1.7   6.804549   3.8368392   9.857077
# 5        P.1 - B.1.1.7   9.800396   1.2378061  19.087202
# 6  B.1.617.1 - B.1.1.7 -45.006922 -55.2132337 -35.471744
# 7  B.1.617.2 - B.1.1.7  73.779169  70.8502375  76.758312

# pairwise growth rate advantages for all strain comparisons (i.e. pairwise differences in growth rate per day among the different lineages)
emtrsanger_pairw = emtrends(fit3_sanger_multi, pairwise ~ LINEAGE2,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_pairw = data.frame(confint(emtrsanger_pairw, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger_pairw$contrasts)$p.value)
delta_r_sanger_pairw
#                  contrast      estimate          SE df     asymp.LCL    asymp.UCL      p.value
# 1         B.1.1.7 - other -0.016681595 0.001008907 NA -0.018659017 -0.0147041741 0.000000e+00
# 2    B.1.1.7 - (B.1.177+)  0.049385970 0.001623129 NA  0.046204695  0.0525672449 0.000000e+00
# 3       B.1.1.7 - B.1.525 -0.009521530 0.004523594 NA -0.018387611 -0.0006554494 4.211893e-01
# 4       B.1.1.7 - B.1.351 -0.014006454 0.003059078 NA -0.020002136 -0.0080107718 4.529264e-04
# 5           B.1.1.7 - P.1 -0.019892330 0.008813869 NA -0.037167197 -0.0026174637 3.306958e-01
# 6     B.1.1.7 - B.1.617.1  0.079066233 0.007383808 NA  0.064594234  0.0935382312 0.000000e+00
# 7     B.1.1.7 - B.1.617.2 -0.117577694 0.001845234 NA -0.121194286 -0.1139611026 0.000000e+00
# 8      other - (B.1.177+)  0.066067566 0.001394475 NA  0.063334444  0.0688006870 0.000000e+00
# 9         other - B.1.525  0.007160065 0.004633482 NA -0.001921393  0.0162415232 7.804032e-01
# 10        other - B.1.351  0.002675142 0.003219824 NA -0.003635597  0.0089858806 9.907702e-01
# 11            other - P.1 -0.003210735 0.008870805 NA -0.020597194  0.0141757237 9.999578e-01
# 12      other - B.1.617.1  0.095747828 0.007451604 NA  0.081142953  0.1103527034 0.000000e+00
# 13      other - B.1.617.2 -0.100896099 0.002096220 NA -0.105004614 -0.0967875836 0.000000e+00
# 14   (B.1.177+) - B.1.525 -0.058907500 0.004806048 NA -0.068327181 -0.0494878200 0.000000e+00
# 15   (B.1.177+) - B.1.351 -0.063392424 0.003463037 NA -0.070179852 -0.0566049956 0.000000e+00
# 16       (B.1.177+) - P.1 -0.069278300 0.008961717 NA -0.086842943 -0.0517136575 8.540557e-10
# 17 (B.1.177+) - B.1.617.1  0.029680263 0.007559699 NA  0.014863525  0.0444970001 4.418915e-03
# 18 (B.1.177+) - B.1.617.2 -0.166963664 0.002453625 NA -0.171772681 -0.1621546477 0.000000e+00
# 19      B.1.525 - B.1.351 -0.004484923 0.005449386 NA -0.015165524  0.0061956770 9.912775e-01
# 20          B.1.525 - P.1 -0.010370800 0.009896920 NA -0.029768407  0.0090268065 9.652464e-01
# 21    B.1.525 - B.1.617.1  0.088587763 0.008651768 NA  0.071630609  0.1055449173 0.000000e+00
# 22    B.1.525 - B.1.617.2 -0.108056164 0.004869099 NA -0.117599423 -0.0985129050 0.000000e+00
# 23          B.1.351 - P.1 -0.005885877 0.009309681 NA -0.024132516  0.0123607625 9.983010e-01
# 24    B.1.351 - B.1.617.1  0.093072687 0.007978274 NA  0.077435558  0.1087098153 0.000000e+00
# 25    B.1.351 - B.1.617.2 -0.103571240 0.003540747 NA -0.110510976 -0.0966315045 0.000000e+00
# 26        P.1 - B.1.617.1  0.098958563 0.011477357 NA  0.076463357  0.1214537692 0.000000e+00
# 27        P.1 - B.1.617.2 -0.097685364 0.008957957 NA -0.115242637 -0.0801280908 0.000000e+00
# 28  B.1.617.1 - B.1.617.2 -0.196643927 0.007616913 NA -0.211572802 -0.1817150523 0.000000e+00


# predicted incidences on average over all ONS regions from multinomial fit
# fitted prop of different LINEAGES in England today
# 94.7% [94.3%-95.1%] now estimated to be B.1.617.2 across all regions
multinom_preds_today_avg = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1   B.1.1.7 2.983606e-02 1.476131e-03 NA 2.694290e-02 3.272922e-02
# 2     other 4.600333e-04 5.338259e-05 NA 3.554053e-04 5.646612e-04
# 3  B.1.177+ 7.494838e-07 1.448432e-07 NA 4.655964e-07 1.033371e-06
# 4   B.1.525 6.346875e-05 1.997884e-05 NA 2.431094e-05 1.026266e-04
# 5   B.1.351 2.085682e-04 4.202209e-05 NA 1.262065e-04 2.909300e-04
# 6       P.1 1.198245e-04 4.689932e-05 NA 2.790349e-05 2.117455e-04
# 7 B.1.617.1 8.651984e-06 3.561596e-06 NA 1.671385e-06 1.563258e-05
# 8 B.1.617.2 9.693026e-01 1.515620e-03 NA 9.663321e-01 9.722732e-01

# 94.8% [94.4%-95.3%] non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.9337270 0.9190977 0.9483562 

# here given by region:
multinom_preds_today_byregion = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|DATE_NUM, by=c("REGION"),
                                                   at=list(DATE_NUM=today_num), 
                                                   mode="prob", df=NA))
multinom_preds_today_byregion
#     LINEAGE2 DATE_NUM                   REGION         prob           SE df     asymp.LCL    asymp.UCL
# 1    B.1.1.7    18794                   London 1.221728e-02 7.510665e-04 NA  1.074521e-02 1.368934e-02
# 2      other    18794                   London 3.597603e-05 4.591137e-06 NA  2.697756e-05 4.497449e-05
# 3   B.1.177+    18794                   London 3.898404e-08 7.841445e-09 NA  2.361509e-08 5.435298e-08
# 4    B.1.525    18794                   London 1.051592e-04 3.503556e-05 NA  3.649077e-05 1.738276e-04
# 5    B.1.351    18794                   London 6.831805e-04 1.375837e-04 NA  4.135214e-04 9.528396e-04
# 6        P.1    18794                   London 4.741142e-04 1.818493e-04 NA  1.176961e-04 8.305322e-04
# 7  B.1.617.1    18794                   London 1.519318e-05 6.233698e-06 NA  2.975353e-06 2.741100e-05
# 8  B.1.617.2    18794                   London 9.864691e-01 8.356850e-04 NA  9.848311e-01 9.881070e-01
# 9    B.1.1.7    18794               North West 7.403741e-03 4.101228e-04 NA  6.599916e-03 8.207567e-03
# 10     other    18794               North West 8.375314e-05 1.025096e-05 NA  6.366163e-05 1.038446e-04
# 11  B.1.177+    18794               North West 2.659150e-07 5.205595e-08 NA  1.638872e-07 3.679428e-07
# 12   B.1.525    18794               North West 3.084107e-05 1.015606e-05 NA  1.093556e-05 5.074658e-05
# 13   B.1.351    18794               North West 3.910634e-05 9.788233e-06 NA  1.992175e-05 5.829092e-05
# 14       P.1    18794               North West 1.498698e-05 8.184145e-06 NA -1.053652e-06 3.102761e-05
# 15 B.1.617.1    18794               North West 8.712793e-08 9.291743e-08 NA -9.498688e-08 2.692427e-07
# 16 B.1.617.2    18794               North West 9.924272e-01 4.190542e-04 NA  9.916059e-01 9.932485e-01
# 17   B.1.1.7    18794               South West 2.561645e-02 2.399804e-03 NA  2.091292e-02 3.031998e-02
# 18     other    18794               South West 2.189299e-04 3.487928e-05 NA  1.505678e-04 2.872920e-04
# 19  B.1.177+    18794               South West 5.022637e-07 1.100402e-07 NA  2.865888e-07 7.179385e-07
# 20   B.1.525    18794               South West 9.124128e-05 5.934176e-05 NA -2.506643e-05 2.075490e-04
# 21   B.1.351    18794               South West 1.888387e-04 1.013602e-04 NA -9.823716e-06 3.875010e-04
# 22       P.1    18794               South West 2.362108e-04 1.829100e-04 NA -1.222863e-04 5.947078e-04
# 23 B.1.617.1    18794               South West 2.390202e-05 1.489152e-05 NA -5.284820e-06 5.308885e-05
# 24 B.1.617.2    18794               South West 9.736239e-01 2.466688e-03 NA  9.687893e-01 9.784585e-01
# 25   B.1.1.7    18794               South East 1.700245e-02 1.143400e-03 NA  1.476142e-02 1.924347e-02
# 26     other    18794               South East 2.349058e-05 3.125338e-06 NA  1.736503e-05 2.961613e-05
# 27  B.1.177+    18794               South East 3.775570e-08 7.667687e-09 NA  2.272731e-08 5.278409e-08
# 28   B.1.525    18794               South East 1.153438e-04 4.193274e-05 NA  3.315712e-05 1.975305e-04
# 29   B.1.351    18794               South East 2.157172e-04 5.927478e-05 NA  9.954080e-05 3.318937e-04
# 30       P.1    18794               South East 2.382114e-04 1.121240e-04 NA  1.845236e-05 4.579704e-04
# 31 B.1.617.1    18794               South East 7.559904e-06 3.787592e-06 NA  1.363614e-07 1.498345e-05
# 32 B.1.617.2    18794               South East 9.823972e-01 1.182078e-03 NA  9.800804e-01 9.847140e-01
# 33   B.1.1.7    18794          East of England 1.671053e-02 1.047595e-03 NA  1.465728e-02 1.876378e-02
# 34     other    18794          East of England 6.286114e-05 8.220456e-06 NA  4.674935e-05 7.897294e-05
# 35  B.1.177+    18794          East of England 7.993897e-08 1.621552e-08 NA  4.815714e-08 1.117208e-07
# 36   B.1.525    18794          East of England 2.506478e-05 1.244107e-05 NA  6.807330e-07 4.944882e-05
# 37   B.1.351    18794          East of England 1.116800e-04 3.340823e-05 NA  4.620106e-05 1.771589e-04
# 38       P.1    18794          East of England 8.198542e-05 4.476341e-05 NA -5.749250e-06 1.697201e-04
# 39 B.1.617.1    18794          East of England 3.578136e-06 1.831544e-06 NA -1.162350e-08 7.167896e-06
# 40 B.1.617.2    18794          East of England 9.830042e-01 1.064439e-03 NA  9.809180e-01 9.850905e-01
# 41   B.1.1.7    18794            East Midlands 3.044034e-02 1.785288e-03 NA  2.694124e-02 3.393944e-02
# 42     other    18794            East Midlands 4.006995e-04 5.085901e-05 NA  3.010176e-04 5.003813e-04
# 43  B.1.177+    18794            East Midlands 8.531977e-07 1.702229e-07 NA  5.195670e-07 1.186828e-06
# 44   B.1.525    18794            East Midlands 1.653046e-05 1.058860e-05 NA -4.222814e-06 3.728373e-05
# 45   B.1.351    18794            East Midlands 1.311173e-04 3.974928e-05 NA  5.321015e-05 2.090245e-04
# 46       P.1    18794            East Midlands 1.531544e-05 1.633394e-05 NA -1.669850e-05 4.732938e-05
# 47 B.1.617.1    18794            East Midlands 6.512770e-06 3.042276e-06 NA  5.500193e-07 1.247552e-05
# 48 B.1.617.2    18794            East Midlands 9.689886e-01 1.817783e-03 NA  9.654258e-01 9.725514e-01
# 49   B.1.1.7    18794            West Midlands 2.734092e-02 1.613048e-03 NA  2.417941e-02 3.050244e-02
# 50     other    18794            West Midlands 3.512833e-04 4.452467e-05 NA  2.640166e-04 4.385501e-04
# 51  B.1.177+    18794            West Midlands 6.545783e-07 1.309863e-07 NA  3.978499e-07 9.113066e-07
# 52   B.1.525    18794            West Midlands 9.855233e-05 3.649137e-05 NA  2.703057e-05 1.700741e-04
# 53   B.1.351    18794            West Midlands 2.181837e-04 5.836516e-05 NA  1.037901e-04 3.325773e-04
# 54       P.1    18794            West Midlands 1.759585e-05 1.872488e-05 NA -1.910423e-05 5.429594e-05
# 55 B.1.617.1    18794            West Midlands 1.690362e-05 7.135790e-06 NA  2.917726e-06 3.088951e-05
# 56 B.1.617.2    18794            West Midlands 9.719559e-01 1.652929e-03 NA  9.687162e-01 9.751956e-01
# 57   B.1.1.7    18794               North East 4.805563e-02 4.956143e-03 NA  3.834177e-02 5.776949e-02
# 58     other    18794               North East 9.825066e-04 1.500632e-04 NA  6.883880e-04 1.276625e-03
# 59  B.1.177+    18794               North East 1.064622e-06 2.320277e-07 NA  6.098559e-07 1.519388e-06
# 60   B.1.525    18794               North East 2.155727e-15 6.115035e-16 NA  9.572025e-16 3.354252e-15
# 61   B.1.351    18794               North East 1.264158e-04 5.936771e-05 NA  1.005721e-05 2.427744e-04
# 62       P.1    18794               North East 3.293275e-10 4.701577e-10 NA -5.921647e-10 1.250820e-09
# 63 B.1.617.1    18794               North East 3.624720e-34 1.581171e-34 NA  5.256812e-35 6.723759e-34
# 64 B.1.617.2    18794               North East 9.508344e-01 5.069733e-03 NA  9.408979e-01 9.607709e-01
# 65   B.1.1.7    18794 Yorkshire and The Humber 8.373722e-02 5.644085e-03 NA  7.267502e-02 9.479943e-02
# 66     other    18794 Yorkshire and The Humber 1.980799e-03 2.436207e-04 NA  1.503312e-03 2.458287e-03
# 67  B.1.177+    18794 Yorkshire and The Humber 3.248099e-06 6.417114e-07 NA  1.990368e-06 4.505830e-06
# 68   B.1.525    18794 Yorkshire and The Humber 8.848586e-05 3.761316e-05 NA  1.476543e-05 1.622063e-04
# 69   B.1.351    18794 Yorkshire and The Humber 1.628746e-04 5.411313e-05 NA  5.681485e-05 2.689344e-04
# 70       P.1    18794 Yorkshire and The Humber 1.686905e-56 2.129616e-55 NA -4.005280e-55 4.342661e-55
# 71 B.1.617.1    18794 Yorkshire and The Humber 4.131104e-06 2.551584e-06 NA -8.699095e-07 9.132117e-06
# 72 B.1.617.2    18794 Yorkshire and The Humber 9.140232e-01 5.792867e-03 NA  9.026694e-01 9.253770e-01


# predicted incidences on average over all ONS regions from multinomial fit
# fitted prop of different LINEAGES in England on the 8th of May, when the UK had 
# similar vaccination coverage than Belgium now (on the 15th of June)
# 94.7% [94.3%-95.1%] now estimated to be B.1.617.2 across all regions
multinom_preds_8may_avg = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=as.numeric(as.Date("2021-05-08"))), 
                                              mode="prob", df=NA))
multinom_preds_8may_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1   B.1.1.7 6.939269e-01 4.411317e-03 NA 6.852809e-01 0.7025729185
# 2     other 4.380954e-03 3.005299e-04 NA 3.791926e-03 0.0049699821
# 3  B.1.177+ 9.900226e-05 1.241198e-05 NA 7.467523e-05 0.0001233293
# 4   B.1.525 1.414500e-03 2.372646e-04 NA 9.494696e-04 0.0018795296
# 5   B.1.351 4.101970e-03 4.062202e-04 NA 3.305793e-03 0.0048981465
# 6       P.1 2.080631e-03 3.550046e-04 NA 1.384835e-03 0.0027764271
# 7 B.1.617.1 4.366731e-03 7.722712e-04 NA 2.853108e-03 0.0058803551
# 8 B.1.617.2 2.896293e-01 4.422646e-03 NA 2.809611e-01 0.2982975433




# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME ####

gentime = 4.7 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantages of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit4_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=sanger, maxit=1000) here
emtrsanger3 = emtrends(fit3_sanger_multi, trt.vs.ctrl1 ~ LINEAGE2, by=c("DATE_NUM"), 
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=seq(as.numeric(as.Date("2021-04-01")),as.numeric(as.Date("2021-07-01")))))
delta_r_sanger3 = data.frame(confint(emtrsanger3, 
                                    adjust="none", df=NA)$contrasts)
delta_r_sanger3 = delta_r_sanger3[delta_r_sanger3$contrast=="B.1.617.2 - B.1.1.7",]
delta_r_sanger3

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger3 = delta_r_sanger3
transmadv_sanger3$estimate =  100*(exp(transmadv_sanger3$estimate*gentime)-1)
transmadv_sanger3$asymp.LCL =  100*(exp(transmadv_sanger3$asymp.LCL*gentime)-1)
transmadv_sanger3$asymp.UCL =  100*(exp(transmadv_sanger3$asymp.UCL*gentime)-1)
transmadv_sanger3$collection_date = as.Date(transmadv_sanger3$DATE_NUM, origin="1970-01-01")
transmadv_sanger3$REGION = factor(transmadv_sanger3$REGION, levels=levels_REGION)

plot_sanger_mfit3_transmadv = qplot(data=transmadv_sanger3, 
                               x=collection_date, y=estimate, ymin=asymp.LCL, ymax=asymp.UCL, # colour=REGION, fill=REGION, 
                               geom="blank") +
  # facet_wrap(~ REGION) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-04-01",NA)), expand=c(0,0)) + #  
  scale_y_continuous(limits=c(0,max(transmadv_sanger3$asymp.UCL))) +
  geom_ribbon(aes(fill=I("steelblue"), colour=NULL), alpha=I(0.3)) +
  geom_line(aes(colour=I("steelblue"), fill=NULL)) +
  ylab("Transmission advantage of B.1.617.2 over B.1.1.7 (%)") +
  theme_hc() + xlab("") +
  ggtitle("TRANSMISSION ADVANTAGE OF B.1.617.2 OVER B.1.1.7 BY REGION\n(Sanger Institute baseline surveillance data, multinomial fit)") +
#  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                     limits=as.Date(c("2021-01-01","2021-05-31")), 
#                     expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  # scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  # scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "none") +
  xlab("") # +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) # +
  # coord_cartesian(xlim=c(as.Date("2021-05-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit3_transmadv

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit3_transm advantage B16172.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit3_transm advantage B16172.pdf"), width=8, height=6)




# PLOT MULTINOMIAL FIT ####

# extrapolate = 60
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-06-30")) # max(sanger$DATE_NUM)+extrapolate

# predictions by ONS region ####

fit_sanger_multi_predsbyregion = data.frame(emmeans(fit3_sanger_multi, 
                                                    ~ LINEAGE2,
                                                    by=c("DATE_NUM", "REGION"),
                                                    at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                                    mode="prob", df=NA))
fit_sanger_multi_predsbyregion$collection_date = as.Date(fit_sanger_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion$LINEAGE2 = factor(fit_sanger_multi_predsbyregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_sanger_multi_predsbyregion$REGION = factor(fit_sanger_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of England today

muller_sangerbyregion_mfit = ggplot(data=fit_sanger_multi_predsbyregion, 
                                    aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + #  
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbyregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit.pdf"), width=8, height=6)

library(ggpubr)
ggarrange(muller_sangerbyregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbyregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.pdf"), width=8, height=6)


# predictions by NHS region

fit_sanger_multi_predsbynhsregion = data.frame(emmeans(fit3_sanger_nhs_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "NHSREGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                                  mode="prob", df=NA))
fit_sanger_multi_predsbynhsregion$collection_date = as.Date(fit_sanger_multi_predsbynhsregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbynhsregion$LINEAGE2 = factor(fit_sanger_multi_predsbynhsregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_sanger_multi_predsbynhsregion$NHSREGION = factor(fit_sanger_multi_predsbynhsregion$NHSREGION, levels=levels_NHSREGION) 
# predicted incidence in different parts of England today
# fit_sanger_multi_predsbynhsregion[fit_sanger_multi_predsbynhsregion$collection_date==today,]

fit_sanger_multi_preds = data.frame(emmeans(fit3_sanger_nhs_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                           mode="prob", df=NA))
fit_sanger_multi_preds$collection_date = as.Date(fit_sanger_multi_preds$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_preds$LINEAGE2 = factor(fit_sanger_multi_preds$LINEAGE2, levels=levels_LINEAGE2_plot) 

muller_sanger_mfit = ggplot(data=fit_sanger_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + #  
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN ENGLAND\n(Sanger Institute baseline surveillance data)")
muller_sanger_mfit

library(ggpubr)
ggarrange(muller_sanger_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sanger_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots_multinom fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots_multinom fit multipanel.pdf"), width=8, height=6)


# predictions by NHS region
muller_sangerbynhsregion_mfit = ggplot(data=fit_sanger_multi_predsbynhsregion, 
                                  aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ NHSREGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + #  
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbynhsregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit.pdf"), width=8, height=6)

library(ggpubr)
ggarrange(muller_sangerbynhsregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbynhsregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit multipanel.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & 95% CONFIDENCE INTERVALS

# PHE SGTF data to superimpose on model predictions as independent estimate of prop of B.1.1.7
sgtf_UK = read.csv(".\\data\\SGTF_UK\\sgtf_2021-05-31_by adm region.csv") # sgtf_2021-05-11_by NHS region, TO DO : ask Nick for update
colnames(sgtf_UK) = c("collection_date", "ereg_code", "sgtf", "non_sgtf", "REGION")
sgtf_UK$ereg_code = NULL
sgtf_UK$collection_date = as.Date(sgtf_UK$collection_date)
sgtf_UK$DATE_NUM = as.numeric(sgtf_UK$collection_date)
sgtf_UK$REGION = factor(sgtf_UK$REGION, levels=levels_REGION)
sgtf_UK = sgtf_UK[sgtf_UK$collection_date>="2021-01-01",]
sgtf_UK$total = sgtf_UK$sgtf + sgtf_UK$non_sgtf
sgtf_UK = sgtf_UK[sgtf_UK$total!=0,]
sgtf_UK$Week = lubridate::week(sgtf_UK$collection_date)
sgtf_UK$LINEAGE2 = "B.1.1.7 (S dropout)"
sgtf_UK$LINEAGE2 = factor(sgtf_UK$LINEAGE2, levels=c(levels_LINEAGE2_plot, "B.1.1.7 (S dropout)"))
sgtf_UK2 = sgtf_UK
sgtf_UK2$count = sgtf_UK2$sgtf 
sgtf_UK2$collection_date_num = sgtf_UK2$DATE_NUM
sgtf_UK2$DATE_NUM = NULL
sgtf_UK2$sgtf = NULL
sgtf_UK2$non_sgtf = NULL
sgtf_UK2$prop = sgtf_UK2$count / sgtf_UK2$total
head(sgtf_UK2)  
range(sgtf_UK$collection_date) # "2021-01-01" "2021-05-30"
library(tidyr)
sgtf_UK_long = gather(sgtf_UK, Sdropout, count, sgtf:non_sgtf, factor_key=TRUE)
sgtf_UK_long$collection_date


# plot of SGTF data by region
ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-01-01"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
  facet_wrap(~REGION, scales="free_y") +
  geom_col(position="fill") +
  scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                 lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                   lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  ylab("Share") +
  ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN ENGLAND\n(data PHE)") +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) + # 
  xlab("Collection date")
  
ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by ONS region_stacked as proportion.png"), width=10, height=8)
ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by ONS region_stacked as proportion.pdf"), width=10, height=8)

ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-01-01"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
  # facet_wrap(~REGION, scales="free_y") +
  geom_col(position="fill") +
  scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                 lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                   lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  ylab("Share") +
  ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN ENGLAND\n(data PHE)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     expand=c(0,0)) +
  xlab("Collection date")

ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.png"), width=10, height=8)
ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.pdf"), width=10, height=8)



# logistic spline fit to SGTF data

fit_SGTF = glm(cbind(sgtf,non_sgtf) ~ REGION * ns(DATE_NUM, df=3), family=binomial, data=sgtf_UK) # at level of ONS regions
BIC(fit_SGTF)

SGTF_predsbyregion = data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                                                      at=list(DATE_NUM=seq(date.from, date.to)), 
                                                      type="response"))
SGTF_predsbyregion$collection_date = as.Date(SGTF_predsbyregion$DATE_NUM, origin="1970-01-01")
SGTF_predsbyregion$LINEAGE2 = "B.1.1.7 (S dropout)"
SGTF_predsbyregion$LINEAGE2 = factor(SGTF_predsbyregion$LINEAGE2, levels=c(levels_LINEAGE2_plot, "B.1.1.7 (S dropout)"))

# fitted S dropout in different parts of England today
# 98% [97%-98%] now estimated to be S positive across all regions
1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM,
        at=list(DATE_NUM=today_num), 
        type="response"))[,c(2,6,5)]
#        prob asymp.UCL asymp.LCL
# 1 0.9689404 0.9641571 0.9731032

# here given by NHS region:
data.frame(REGION=data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                   at=list(DATE_NUM=today_num), 
                   type="response"))[,2], 
           1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                        at=list(DATE_NUM=today_num), 
                        type="response"))[,c(3,7,6)])
#                    REGION      prob  asymp.UCL asymp.LCL
# 1                   London 0.9173624 0.8995004 0.9322887
# 2               North West 0.9947801 0.9939643 0.9954861
# 3               South West 0.9824263 0.9485692 0.9941330
# 4               South East 0.9786745 0.9713814 0.9841395
# 5          East of England 0.9706332 0.9629772 0.9767442
# 6            East Midlands 0.9523325 0.9412206 0.9614298
# 7            West Midlands 0.9615704 0.9523972 0.9690333
# 8               North East 0.9467920 0.9217294 0.9641416
# 9 Yorkshire and The Humber 0.9460530 0.9316387 0.9575665




# plot predictions multinomial fit by ONS region
# on logit scale:

fit_sanger_multi_preds2 = fit_sanger_multi_preds
ymin = 0.001
ymax = 0.998
fit_sanger_multi_preds2$asymp.LCL[fit_sanger_multi_preds2$asymp.LCL<ymin] = ymin
fit_sanger_multi_preds2$asymp.UCL[fit_sanger_multi_preds2$asymp.UCL<ymin] = ymin
fit_sanger_multi_preds2$asymp.UCL[fit_sanger_multi_preds2$asymp.UCL>ymax] = ymax
fit_sanger_multi_preds2$prob[fit_sanger_multi_preds2$prob<ymin] = ymin

fit_sanger_multi_predsbyregion2 = fit_sanger_multi_predsbyregion
fit_sanger_multi_predsbyregion2$asymp.LCL[fit_sanger_multi_predsbyregion2$asymp.LCL<ymin] = ymin
fit_sanger_multi_predsbyregion2$asymp.UCL[fit_sanger_multi_predsbyregion2$asymp.UCL<ymin] = ymin
fit_sanger_multi_predsbyregion2$asymp.UCL[fit_sanger_multi_predsbyregion2$asymp.UCL>ymax] = ymax
fit_sanger_multi_predsbyregion2$prob[fit_sanger_multi_predsbyregion2$prob<ymin] = ymin

# plots of multinomial fit to Sanger Inst data by ONS region

plot_sanger_mfit_logit = qplot(data=fit_sanger_multi_predsbyregion2, 
                               x=collection_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + # 
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(1), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.pdf"), width=8, height=6)
library(svglite)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.svg"), width=8, height=6)

# on response scale:
plot_sanger_mfit = qplot(data=fit_sanger_multi_predsbyregion2, 
                         x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + # 
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-09-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(1), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_sanger_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_response scale.pdf"), width=8, height=6)


# plots of multinomial fit to Sanger Inst data by ONS region with PHE S dropout data overlaid

plot_sanger_mfit_logit = qplot(data=rbind(fit_sanger_multi_predsbyregion2[,1:9], SGTF_predsbyregion[,1:9]), 
                               x=collection_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + # 
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(1), pch=I(16)) +
  geom_point(data=sgtf_UK2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(0.5), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.pdf"), width=8, height=6)
library(svglite)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.svg"), width=8, height=6)
write.csv(fit_sanger_multi_predsbyregion2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_multinomial fit sanger inst data.csv"), row.names=F)
write.csv(SGTF_predsbyregion, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logistic fit S dropout data.csv"), row.names=F)
write.csv(data_agbyweekregion1, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_sanger inst data aggregated by week.csv"), row.names=F)
write.csv(sgtf_UK2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_S dropout data.csv"), row.names=F)

# on response scale:
plot_sanger_mfit = qplot(data=rbind(fit_sanger_multi_predsbyregion2[,1:9], SGTF_predsbyregion[,1:9]), 
                         x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + # 
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-09-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(1), pch=I(16)) +
  geom_point(data=sgtf_UK2,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(0.5), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_sanger_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_response scale.pdf"), width=8, height=6)


# PROJECT MULTINOMIAL FIT ONTO CASE DATA ####

# case data from https://coronavirus.data.gov.uk/details/download

# cases_uk_ltla = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesByPublishDateRollingRate&metric=newCasesBySpecimenDate&format=csv")
# cases_uk_ltla$date = as.Date(cases_uk_ltla$date)
# cases_uk_ltla$DATE_NUM = as.numeric(cases_uk_ltla$date)
# cases_uk_ltla$LTLA = factor(cases_uk_ltla$areaCode)
# cases_uk_ltla$areaCode = NULL
# cases_uk_ltla$REGION = LTLAs_regions$RGN20NM[match(cases_uk_ltla$LTLA, LTLAs_regions$LAD20CD)]
# cases_uk_ltla = cases_uk_ltla[!is.na(cases_uk_ltla$REGION),]
# cases_uk_ltla$REGION = factor(cases_uk_ltla$REGION, levels=levels_REGION)

# daily new cases by region :
cases_uk_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDate&format=csv")
cases_uk_region$date = as.Date(cases_uk_region$date)
cases_uk_region$DATE_NUM = as.numeric(cases_uk_region$date)
cases_uk_region$REGION = factor(cases_uk_region$areaName, levels=levels_REGION)
cases_uk_region$areaName = NULL

ggplot(data=cases_uk_region, 
       aes(x=date, y=newCasesBySpecimenDate, group=REGION, colour=REGION, fill=REGION)) +
  facet_wrap(~REGION) +
  geom_point(cex=I(0.2)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  scale_y_log10() + ylab("New cases (per day)") + ggtitle("CONFIRMED SARS-CoV2 CASES PER DAY IN ENGLAND BY ONS REGION")

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases England by ONS region.png"), width=8, height=6)

# daily  new cases by ONS region & age group :
cases_uk_age_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateAgeDemographics&format=csv")
cases_uk_age_region$date = as.Date(cases_uk_age_region$date)
cases_uk_age_region$DATE_NUM = as.numeric(cases_uk_age_region$date)
cases_uk_age_region$REGION = factor(cases_uk_age_region$areaName, levels=levels_REGION)
cases_uk_age_region$areaName = NULL

ggplot(data=cases_uk_age_region[cases_uk_age_region$age %in% c("00_59","60+"),], 
        aes(x=date, y=cases, group=age, colour=age, fill=age)) +
          geom_point(cex=I(0.2)) +
          geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
          facet_wrap(~REGION) + scale_y_log10()

ggplot(data=cases_uk_age_region[cases_uk_age_region$age %in% c("00_59","60+"),], 
       aes(x=date, y=rollingSum/7, group=age, colour=age, fill=age)) +
  geom_point(cex=I(0.2)) +
  # geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  facet_wrap(~REGION) + scale_y_log10()

# daily new cases by age group for England :
cases_ENG_age = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDateAgeDemographics&format=csv")
cases_ENG_age$date = as.Date(cases_ENG_age$date)
cases_ENG_age$DATE_NUM = as.numeric(cases_ENG_age$date)
cases_ENG_age$areaName = NULL

ggplot(data=cases_ENG_age[cases_ENG_age$age %in% c("00_59","60+"),], 
       aes(x=date, y=cases, group=age, colour=age, fill=age)) +
  geom_point(cex=I(0.2)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  scale_y_log10()

# daily hospital admissions by NHS region
hosps_uk_nhsregion = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=newAdmissions&format=csv") # &metric=newCasesBySpecimenDate gives NAs
hosps_uk_nhsregion$date = as.Date(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$DATE_NUM = as.numeric(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$REGION = factor(hosps_uk_nhsregion$areaName, levels=levels_NHSREGION)

ggplot(data=hosps_uk_nhsregion, 
       aes(x=date, y=newAdmissions, group=REGION, colour=REGION, fill=REGION)) +
  facet_wrap(~REGION) +
  geom_point(cex=I(0.2)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  scale_y_log10() + 
  ylab("New hospital admissions (per day)") + 
  ggtitle("HOSPITAL ADMISSIONS PER DAY IN ENGLAND BY NHS REGION")

ggsave(file=paste0(".\\plots\\",plotdir,"\\hospital admissions England by NHS region.png"), width=8, height=6)


# daily hospital admission also see 
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/COVID-19-daily-admissions-and-beds-20210611.xlsx

# daily hospital admission by age for England
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/Covid-Publication-10-06-2021-Supplementary-Data.xlsx
library(rio)
dat = rio::import(file = "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/Covid-Publication-10-06-2021-Supplementary-Data.xlsx",
            which = 1)
dat[,1] = NULL
dates = as.Date(as.vector(unlist(dat[11,-1])), origin="1899-12-30")
ages = gsub("Total reported admissions and diagnoses  ", "", dat[14:(nrow(dat)-3),1])
dat = dat[14:(nrow(dat)-3),-1]
colnames(dat) = dates
dat = data.frame(age=ages, dat, check.names=F)
library(tidyr)
hosps_ENG_age = gather(dat, date, newAdmissions, all_of(as.character(dates)), factor_key=TRUE)
hosps_ENG_age = hosps_ENG_age[hosps_ENG_age$age!="Unknown age",]
hosps_ENG_age$age = factor(hosps_ENG_age$age, levels=c("0-5","6-17","18-54","55-64","65-74","75-84","85+"))
hosps_ENG_age$date = as.Date(as.character(hosps_ENG_age$date))

# combined with daily new cases by age (matched to closest age category)
cases_hosps_ENG_age = rbind(data.frame(type="hospitalisations",hosps_ENG_age), 
                            data.frame(type="cases",hosps_ENG_age))
colnames(cases_hosps_ENG_age)[4] = "count"
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="00_04"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                    cases_ENG_age[cases_ENG_age$age=="00_04","date"])]  
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="6-17"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="05_09"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                    cases_ENG_age[cases_ENG_age$age=="05_09","date"])] +
                                                                                               cases_ENG_age$cases[cases_ENG_age$age=="10_14"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                     cases_ENG_age[cases_ENG_age$age=="10_14","date"])] +
                                                                                               cases_ENG_age$cases[cases_ENG_age$age=="15_19"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                     cases_ENG_age[cases_ENG_age$age=="15_19","date"])]
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="18-54"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="20_24"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="20_24","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="25_29"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="25_29","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="30_34"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="30_34","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="35_39"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="35_39","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="40_44"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="40_44","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="45_49"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="45_49","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="50_54"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="50_54","date"])] 
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="55-64"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="55_59"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="55_59","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="60_64"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="60_64","date"])]
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="65-74"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="65_69"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="65_69","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="70_74"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="70_74","date"])] 
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="75-84"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="75_79"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="75_79","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="80_84"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="80_84","date"])] 
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="85+"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="85_89"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="85_89","date"])] +
                                                                                              cases_ENG_age$cases[cases_ENG_age$age=="90+"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                    cases_ENG_age[cases_ENG_age$age=="90+","date"])]
cases_hosps_ENG_age$age = factor(cases_hosps_ENG_age$age, levels=c("0-5","6-17","18-54","55-64","65-74","75-84","85+"))
ggplot(data=cases_hosps_ENG_age, 
       aes(x=date, y=count, group=age, colour=age, fill=age)) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) + # 
  geom_point(cex=I(0.8)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 35)) + 
  facet_wrap(~ type) + scale_y_log10() + ylab("Daily new confirmed cases & hospital admissions") + 
  ggtitle("Daily new SARS-CoV2 cases & hospitalisations in England\n(data gov.uk & NHS)")

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.pdf"), width=8, height=6)

write.csv(cases_hosps_ENG_age, file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.csv"), row.names=F)


# MAP MULTINOMIAL FIT ONTO CASE NUMBERS ####

newdat = expand.grid(REGION=levels_REGION,
                     DATE_NUM=seq(date.from, date.to))
fit_sanger_multi_predsbyregion_day = data.frame(newdat,
                                          predict(fit3_sanger_multi, 
                                                  newdata = newdat,
                                                  type = "prob"), check.names=F)  
fit_sanger_multi_predsbyregion_day = gather(fit_sanger_multi_predsbyregion_day, LINEAGE2, prob, all_of(levels_LINEAGE2_plot))
fit_sanger_multi_predsbyregion_day$collection_date = as.Date(fit_sanger_multi_predsbyregion_day$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion_day$LINEAGE2 = factor(fit_sanger_multi_predsbyregion_day$LINEAGE2, levels=levels_LINEAGE2_plot)

fit_sanger_multi_predsbyregion_day$totcases_rollingmean = cases_uk_region$newCasesBySpecimenDateRollingSum[match(interaction(fit_sanger_multi_predsbyregion_day$collection_date,
                                                                                                                         fit_sanger_multi_predsbyregion_day$REGION),
                                                                                        interaction(cases_uk_region$date,
                                                                                                    cases_uk_region$REGION))]/7
fit_sanger_multi_predsbyregion_day$totcases = cases_uk_region$newCasesBySpecimenDate[match(interaction(fit_sanger_multi_predsbyregion_day$collection_date,
                                                                                                       fit_sanger_multi_predsbyregion_day$REGION),
                                                                                                 interaction(cases_uk_region$date,
                                                                                                             cases_uk_region$REGION))]
fit_sanger_multi_predsbyregion_day$WEEKDAY = factor(weekdays(fit_sanger_multi_predsbyregion_day$collection_date))
library(mgcv)
gamfittotcases = gam(totcases ~ s(DATE_NUM, bs="cs", k=7, by=REGION) + WEEKDAY, family=poisson(log), 
                     data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$collection_date<=max(sanger$WeekEndDate),])
BIC(gamfittotcases) 
gamfitemmeans = as.data.frame(emmeans(gamfittotcases, ~ DATE_NUM, by=c("REGION","DATE_NUM"),
                                                                         at=list(REGION=levels(fit_sanger_multi_predsbyregion_day$REGION),
                                                                           DATE_NUM=seq(min(fit_sanger_multi_predsbyregion_day$DATE_NUM),
                                                         max(fit_sanger_multi_predsbyregion_day$DATE_NUM)+14)), type="response"))
fit_sanger_multi_predsbyregion_day$totcases_smoothed = gamfitemmeans$rate[match(interaction(fit_sanger_multi_predsbyregion_day$REGION,
                                                                                            fit_sanger_multi_predsbyregion_day$DATE_NUM),
                                                                            interaction(gamfitemmeans$REGION,gamfitemmeans$DATE_NUM))]
fit_sanger_multi_predsbyregion_day$cases = fit_sanger_multi_predsbyregion_day$totcases * fit_sanger_multi_predsbyregion_day$prob
fit_sanger_multi_predsbyregion_day$cases_rollingmean = fit_sanger_multi_predsbyregion_day$totcases_rollingmean * fit_sanger_multi_predsbyregion_day$prob
fit_sanger_multi_predsbyregion_day$cases_smoothed = fit_sanger_multi_predsbyregion_day$totcases_smoothed * fit_sanger_multi_predsbyregion_day$prob
fit_sanger_multi_predsbyregion_day$REGION = factor(fit_sanger_multi_predsbyregion_day$REGION, levels=levels_REGION)
fit_sanger_multi_predsbyregion_day$cases[fit_sanger_multi_predsbyregion_day$cases<=0.001] = NA
fit_sanger_multi_predsbyregion_day$cases_smoothed[fit_sanger_multi_predsbyregion_day$cases_smoothed<=0.001] = NA

cases_uk_region$totcases_smoothed = gamfitemmeans$rate[match(interaction(cases_uk_region$REGION,cases_uk_region$DATE_NUM),
                                                             interaction(gamfitemmeans$REGION,gamfitemmeans$DATE_NUM))]

# plot new cases per day by region
ggplot(data=fit_sanger_multi_predsbyregion_day,
       aes(x=collection_date, y=cases_smoothed, 
           group=LINEAGE2)) +
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_line(aes(lwd=I(1.3), colour=LINEAGE2)) +
#  geom_smooth(aes(lwd=I(1.3), colour=LINEAGE2), method = "gam", family = "poisson", formula = y ~ s(x, k=4), se=FALSE) +
# geom_smooth(aes(lwd=I(1), colour=LINEAGE2), method="loess", span=0.8, se=FALSE) +
  geom_line(data=data.frame(collection_date=cases_uk_region$date,
                            cases_smoothed=cases_uk_region$totcases_smoothed, 
                              REGION=cases_uk_region$REGION,
                              LINEAGE2="total"), aes(lwd=I(1.9), alpha=I(0.5)), 
              colour=alpha(I("black"),0.7)) +
#  geom_smooth(data=data.frame(collection_date=cases_uk_region$date,
#                              cases=cases_uk_region$newCasesBySpecimenDate, 
#                              REGION=cases_uk_region$REGION,
#                              LINEAGE2="total"), aes(lwd=I(1.9), alpha=I(0.5)), 
#              method = "gam", family = "poisson", formula = y ~ s(x, k=4), se=FALSE, colour=alpha(I("black"),0.7)) +
  #  geom_smooth(data=data.frame(collection_date=cases_uk_region$date,
#                              cases=cases_uk_region$newCasesBySpecimenDate, 
#                              REGION=cases_uk_region$REGION,
#                              LINEAGE2="total"), aes(lwd=I(1.5)), method="loess", span=0.8, se=FALSE, colour=I("black")) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_uk_region$date)+14), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_y_log10() +
  scale_fill_manual("variant", values=c(lineage_cols1,I("black"))) +
  scale_colour_manual("variant", values=c(lineage_cols1,I("black"))) 

# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day by lineage multinomial fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day by lineage multinomial fit.pdf"), width=8, height=6)


ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-03-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases_rollingmean, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-03-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data.pdf"), width=8, height=6)

# ggplot(data=fit_sanger_multi_predsbyregion_day, 
#        aes(x=collection_date, y=cases_smoothed, group=LINEAGE2)) + 
#   facet_wrap(~ REGION, scale="free", ncol=3) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                      limits=as.Date(c("2021-03-01","2021-06-30")), expand=c(0,0)) +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("New confirmed cases per day") +
#   ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
#   scale_fill_manual("variant", values=lineage_cols1) +
#   scale_colour_manual("variant", values=lineage_cols1) 
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.pdf"), width=8, height=6)



# CALCULATE Re VALUES PER VARIANT ####

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}
library(mgcv)

qplot(data=fit_sanger_multi_predsbyregion_day, x=collection_date, y=totcases, geom="line") +
  facet_wrap(~ REGION) + scale_y_log10()

fit_sanger_multi_predsbyregion_day$weekday = factor(weekdays(fit_sanger_multi_predsbyregion_day$collection_date))
fit_cases = gam(totcases ~ s(DATE_NUM, bs="cs", k=8, by=REGION) + weekday,
                family=poisson(log), data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE2==levels_LINEAGE2[[1]],],
)
BIC(fit_cases)
# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", by="REGION",
                              at=list(DATE_NUM=seq(date.from,
                                                   date.to),
                                      REGION=levels_REGION
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
avg_r_cases$REGION = factor(avg_r_cases$REGION, levels=levels_REGION)
qplot(data=avg_r_cases, x=DATE, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue"), group=REGION) +
  facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J")) +
  scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re VALUES IN ENGLAND BY REGION BASED ON NEW CONFIRMED CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
  # coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants3 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
       function (dat) { 
  wt = as.data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
  # wt = rep(1/length(levels_LINEAGE2), length(levels_LINEAGE2)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
  cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
  names(cons) = seq_along(cons)
  EMT = emtrends(fit3_sanger_multi,  ~ LINEAGE2, by=c("DATE_NUM","REGION"),
               var="DATE_NUM", mode="latent",
               at=list(DATE_NUM=dat,
                       REGION=levels(SGTF_predsbyregion$REGION)))
  out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
  # sum(out$estimate*wt) # should sum to zero
  return(out) } ))

# for simpler model fit1_sanger_multi
above_avg_r_variants1 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
                                              function (dat) { 
                                                wt = as.data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_LINEAGE2), length(levels_LINEAGE2)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(fit1_sanger_multi,  ~ LINEAGE2, by=c("DATE_NUM","REGION"),
                                                               var="DATE_NUM", mode="latent",
                                                               at=list(DATE_NUM=dat,
                                                                       REGION=levels(SGTF_predsbyregion$REGION)))
                                                out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                # sum(out$estimate*wt) # should sum to zero
                                                return(out) } ))
above_avg_r_variants = above_avg_r_variants3
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE2), 
                                       labels=levels(data_agbyweeknhsregion1$LINEAGE2))
above_avg_r_variants$LINEAGE2 = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$DATE_NUM = above_avg_r_variants$DATE_NUM-0 # the -7 to calculate back to time of infection
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2020-12-22" "2021-05-25"
above_avg_r_variants$avg_r = avg_r_cases$r[match(interaction(above_avg_r_variants$REGION,above_avg_r_variants$collection_date),
                                                 interaction(avg_r_cases$REGION,avg_r_cases$DATE))]  # average growth rate of all lineages calculated from case nrs
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM-0, # -7 to calculate back to time of infection
                REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE-0,
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
above_avg_r_variants$prob = fit_sanger_multi_predsbyregion_day$prob[match(interaction(above_avg_r_variants$DATE_NUM,
                                                                                      above_avg_r_variants$LINEAGE2,
                                                                                      above_avg_r_variants$REGION),
                                                                          interaction(fit_sanger_multi_predsbyregion_day$DATE_NUM,
                                                                                      fit_sanger_multi_predsbyregion_day$LINEAGE2,
                                                                                      fit_sanger_multi_predsbyregion_day$REGION))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 4
ymin = 1/ymax
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
qplot(data=above_avg_r_variants2[!(above_avg_r_variants2$LINEAGE2 %in% c("other")),], 
      x=collection_date, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE2, fill=LINEAGE2, alpha=I(0.5),
      group=LINEAGE2, linetype=I(0)) +
  facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=LINEAGE2, colour=LINEAGE2), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE2), lwd=I(0.72)) + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J")) +
  scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN ENGLAND BY REGION") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) +
  scale_fill_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  scale_colour_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  theme(legend.position="right",  
        axis.title.x=element_blank())

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant_with clipping.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant_with clipping.pdf"), width=8, height=6)


