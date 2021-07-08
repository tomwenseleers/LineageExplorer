# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN ENGLAND BASED ON SANGER INSTITUTE BASELINE SURVEILLANCE SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & most surge testing/active surveillance)

# last update 5 JULY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-07-05")
today_num = as.numeric(today)
plotdir = "UK_SANGER"
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
range(sanger$WeekEndDate) # "2021-09-05" "2021-06-26"
# sanger$Week = lubridate::week(sanger$WeekEndDate)
sanger$DATE_NUM = as.numeric(sanger$WeekEndDate)-3.5 # using week midpoint
colnames(sanger)

sanger = sanger[rep(seq_len(nrow(sanger)), sanger$Count),] # convert to long format
sanger$Count = NULL
nrow(sanger) # 280752

nrow(sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 36754 (last 2 weeks)
nrow(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 35791 (last 2 weeks)
nrow(sanger[grepl("B.1.1.7",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 907

length(unique(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),"LTLA"])) # 311 LTLAs

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

table_lineage = as.data.frame(table(sanger$LINEAGE1))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger$LINEAGE1))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# just for the last week worth of data
table_lineage = as.data.frame(table(sanger[sanger$WeekEndDate==max(sanger$WeekEndDate),]$LINEAGE2))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger[sanger$WeekEndDate==max(sanger$WeekEndDate),]$LINEAGE2))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# Lineage        Prop
# 3   B.1.1.7 0.009735627
# 5 B.1.617.2 0.989341167


sel_ref_lineage = "B.1.1.7"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"Lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.1.7","B.1.177+","B.1.351", "B.1.525","B.1.617+","B.1.617.1","B.1.617.2","P.1")

sanger$LINEAGE1[!(sanger$LINEAGE1 %in% sel_lineages)] = "other"
sanger$LINEAGE2[!(sanger$LINEAGE2 %in% sel_lineages)] = "other"
# sanger = sanger[sanger$LINEAGE1 %in% sel_lineages, ]
sum(table(sanger$LINEAGE1))
table(sanger$LINEAGE1)
table(sanger$LINEAGE2)

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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_raw data.pdf"), width=8, height=6)


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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_raw data.pdf"), width=8, height=6)



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
# fit4_sanger_multi fits best (lowest BIC)


# growth rate advantages of different VOCs compared to UK type B.1.1.7 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrsanger = emtrends(fit4_sanger_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger = data.frame(confint(emtrsanger, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger$contrasts)$p.value)
delta_r_sanger
# contrast    estimate           SE df     asymp.LCL   asymp.UCL      p.value
# 1      other - B.1.1.7  0.03881435 0.001544900 NA  0.03578640  0.04184230 0.0000000000
# 2 (B.1.177+) - B.1.1.7 -0.03319271 0.002629929 NA -0.03834728 -0.02803815 0.0000000000
# 3    B.1.525 - B.1.1.7 -0.04780911 0.013747618 NA -0.07475395 -0.02086428 0.0040781597
# 4    B.1.351 - B.1.1.7  0.01920245 0.004481360 NA  0.01041914  0.02798575 0.0001980897
# 5        P.1 - B.1.1.7 -0.11705918 0.026790560 NA -0.16956771 -0.06455064 0.0001401820
# 6  B.1.617.1 - B.1.1.7 -0.12864249 0.011682278 NA -0.15153933 -0.10574564 0.0000000000
# 7  B.1.617.2 - B.1.1.7  0.11785610 0.001835523 NA  0.11425854  0.12145366 0.0000000000

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*5.5)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
# contrast   estimate   asymp.LCL asymp.UCL
# 1      other - B.1.1.7   23.79774   21.753122  25.87670
# 2 (B.1.177+) - B.1.1.7  -20.02861  -23.480126 -16.67356
# 3    B.1.525 - B.1.1.7  -30.07618  -50.854669 -12.15970
# 4    B.1.351 - B.1.1.7   11.13922    5.897904  16.63995
# 5        P.1 - B.1.1.7  -90.37497 -154.116440 -42.62214
# 6  B.1.617.1 - B.1.1.7 -102.89810 -130.128194 -78.89002
# 7  B.1.617.2 - B.1.1.7   91.21123   87.464996  95.03232

# so this would estimate that B.1.617.2 was 91% more infectious than B.1.1.7 [87%-95%] 95% CLs

# or with generation time of 4.7 days (Nishiura et al. 2020)
transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*4.7)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
# contrast   estimate   asymp.LCL  asymp.UCL
# 1      other - B.1.1.7  20.01271   18.316860  21.73286
# 2 (B.1.177+) - B.1.1.7 -16.88329  -19.749541 -14.08565
# 3    B.1.525 - B.1.1.7 -25.19506  -42.097544 -10.30313
# 4    B.1.351 - B.1.1.7   9.44495    5.018881  14.05756
# 5        P.1 - B.1.1.7 -73.35618 -121.880383 -35.44399
# 6  B.1.617.1 - B.1.1.7 -83.05559 -103.854204 -64.37900
# 7  B.1.617.2 - B.1.1.7  74.00671   71.089243  76.97392

# pairwise growth rate advantages for all strain comparisons (i.e. pairwise differences in growth rate per day among the different lineages)
emtrsanger_pairw = emtrends(fit4_sanger_multi, pairwise ~ LINEAGE2,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_pairw = data.frame(confint(emtrsanger_pairw, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger_pairw$contrasts)$p.value)
delta_r_sanger_pairw
# contrast    estimate          SE df   asymp.LCL   asymp.UCL      p.value
# 1         B.1.1.7 - other -0.03881435 0.001544900 NA -0.04184230 -0.03578640 0.000000e+00
# 2    B.1.1.7 - (B.1.177+)  0.03319271 0.002629929 NA  0.02803815  0.03834728 0.000000e+00
# 3       B.1.1.7 - B.1.525  0.04780911 0.013747618 NA  0.02086428  0.07475395 1.425967e-02
# 4       B.1.1.7 - B.1.351 -0.01920245 0.004481360 NA -0.02798575 -0.01041914 7.528841e-04
# 5           B.1.1.7 - P.1  0.11705918 0.026790560 NA  0.06455064  0.16956771 5.355400e-04
# 6     B.1.1.7 - B.1.617.1  0.12864249 0.011682278 NA  0.10574564  0.15153933 0.000000e+00
# 7     B.1.1.7 - B.1.617.2 -0.11785610 0.001835523 NA -0.12145366 -0.11425854 0.000000e+00
# 8      other - (B.1.177+)  0.07200707 0.002303176 NA  0.06749292  0.07652121 0.000000e+00
# 9         other - B.1.525  0.08662347 0.013832191 NA  0.05951287  0.11373406 6.943774e-08
# 10        other - B.1.351  0.01961191 0.004727180 NA  0.01034680  0.02887701 1.288021e-03
# 11            other - P.1  0.15587353 0.026832358 NA  0.10328307  0.20846398 7.266475e-07
# 12      other - B.1.617.1  0.16745684 0.011782163 NA  0.14436423  0.19054946 0.000000e+00
# 13      other - B.1.617.2 -0.07904174 0.002359554 NA -0.08366639 -0.07441710 0.000000e+00
# 14   (B.1.177+) - B.1.525  0.01461640 0.013995547 NA -0.01281437  0.04204717 9.668974e-01
# 15   (B.1.177+) - B.1.351 -0.05239516 0.005183890 NA -0.06255540 -0.04223492 1.854072e-14
# 16       (B.1.177+) - P.1  0.08386646 0.026916540 NA  0.03111101  0.13662191 4.330189e-02
# 17 (B.1.177+) - B.1.617.1  0.09544978 0.011973591 NA  0.07198197  0.11891758 4.071965e-12
# 18 (B.1.177+) - B.1.617.2 -0.15104881 0.003177206 NA -0.15727602 -0.14482160 0.000000e+00
# 19      B.1.525 - B.1.351 -0.06701156 0.014453849 NA -0.09534058 -0.03868254 1.759214e-04
# 20          B.1.525 - P.1  0.06925006 0.030107712 NA  0.01024003  0.12826009 2.990893e-01
# 21    B.1.525 - B.1.617.1  0.08083337 0.018033000 NA  0.04548934  0.11617740 3.363488e-04
# 22    B.1.525 - B.1.617.2 -0.16566521 0.013855222 NA -0.19282095 -0.13850948 0.000000e+00
# 23          B.1.351 - P.1  0.13626162 0.027157130 NA  0.08303463  0.18948862 3.268526e-05
# 24    B.1.351 - B.1.617.1  0.14784493 0.012502068 NA  0.12334133  0.17234854 0.000000e+00
# 25    B.1.351 - B.1.617.2 -0.09865365 0.004787434 NA -0.10803685 -0.08927045 0.000000e+00
# 26        P.1 - B.1.617.1  0.01158331 0.029237597 NA -0.04572132  0.06888795 9.999266e-01
# 27        P.1 - B.1.617.2 -0.23491527 0.026819609 NA -0.28748074 -0.18234981 4.785061e-14
# 28  B.1.617.1 - B.1.617.2 -0.24649859 0.011811466 NA -0.26964863 -0.22334854 0.000000e+00


# predicted incidences on average over all ONS regions from multinomial fit
# fitted prop of different LINEAGES in England today
# 99.63% [99.57%-999.69%] now estimated to be B.1.617.2 across all regions
multinom_preds_today_avg = data.frame(emmeans(fit4_sanger_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# LINEAGE2         prob           SE df     asymp.LCL    asymp.UCL
# 1   B.1.1.7 2.954502e-03 2.265854e-04 NA  2.510403e-03 3.398601e-03
# 2     other 5.183601e-04 6.227155e-05 NA  3.963101e-04 6.404101e-04
# 3  B.1.177+ 1.036337e-07 4.334769e-08 NA  1.867378e-08 1.885936e-07
# 4   B.1.525 5.036763e-05 4.831555e-05 NA -4.432911e-05 1.450644e-04
# 5   B.1.351 1.415020e-04 7.154772e-05 NA  1.271045e-06 2.817329e-04
# 6       P.1 1.541517e-05 8.470995e-06 NA -1.187677e-06 3.201801e-05
# 7 B.1.617.1 1.377736e-05 3.378513e-05 NA -5.244028e-05 7.999501e-05
# 8 B.1.617.2 9.963060e-01 2.863106e-04 NA  9.957448e-01 9.968671e-01

# 99.7% [99.6%-99.8%] non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.9970455 0.9960445 0.9980465 

# here given by region:
multinom_preds_today_byregion = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|DATE_NUM, by=c("REGION"),
                                                   at=list(DATE_NUM=today_num), 
                                                   mode="prob", df=NA))
multinom_preds_today_byregion
# LINEAGE2 DATE_NUM                   REGION         prob           SE df     asymp.LCL    asymp.UCL
# 1    B.1.1.7    18813                   London 1.666193e-03 1.193941e-04 NA  1.432185e-03 1.900201e-03
# 2      other    18813                   London 2.800717e-05 3.567232e-06 NA  2.101553e-05 3.499882e-05
# 3   B.1.177+    18813                   London 2.445328e-09 6.092829e-10 NA  1.251155e-09 3.639500e-09
# 4    B.1.525    18813                   London 1.845143e-05 6.807011e-06 NA  5.109930e-06 3.179292e-05
# 5    B.1.351    18813                   London 1.863911e-04 3.925581e-05 NA  1.094511e-04 2.633311e-04
# 6        P.1    18813                   London 9.333754e-05 4.250577e-05 NA  1.002776e-05 1.766473e-04
# 7  B.1.617.1    18813                   London 1.807454e-07 8.419723e-08 NA  1.572188e-08 3.457690e-07
# 8  B.1.617.2    18813                   London 9.980074e-01 1.468807e-04 NA  9.977196e-01 9.982953e-01
# 9    B.1.1.7    18813               North West 1.017427e-03 6.879610e-05 NA  8.825892e-04 1.152265e-03
# 10     other    18813               North West 6.265656e-05 7.708334e-06 NA  4.754851e-05 7.776462e-05
# 11  B.1.177+    18813               North West 1.626935e-08 3.984869e-09 NA  8.459150e-09 2.407955e-08
# 12   B.1.525    18813               North West 4.793657e-06 1.802408e-06 NA  1.261002e-06 8.326312e-06
# 13   B.1.351    18813               North West 9.779989e-06 2.554857e-06 NA  4.772561e-06 1.478742e-05
# 14       P.1    18813               North West 2.693128e-06 1.578714e-06 NA -4.010939e-07 5.787350e-06
# 15 B.1.617.1    18813               North West 1.490135e-09 1.243997e-09 NA -9.480546e-10 3.928324e-09
# 16 B.1.617.2    18813               North West 9.989026e-01 7.404598e-05 NA  9.987575e-01 9.990478e-01
# 17   B.1.1.7    18813               South West 2.782273e-03 2.313716e-04 NA  2.328793e-03 3.235753e-03
# 18     other    18813               South West 1.299041e-04 1.917492e-05 NA  9.232195e-05 1.674862e-04
# 19  B.1.177+    18813               South West 2.429684e-08 6.253294e-09 NA  1.204060e-08 3.655307e-08
# 20   B.1.525    18813               South West 1.362977e-05 8.138257e-06 NA -2.320923e-06 2.958046e-05
# 21   B.1.351    18813               South West 5.467179e-05 2.291603e-05 NA  9.757198e-06 9.958638e-05
# 22       P.1    18813               South West 2.256014e-05 1.785026e-05 NA -1.242573e-05 5.754601e-05
# 23 B.1.617.1    18813               South West 3.591016e-07 1.976281e-07 NA -2.824240e-08 7.464456e-07
# 24 B.1.617.2    18813               South West 9.969966e-01 2.489434e-04 NA  9.965087e-01 9.974845e-01
# 25   B.1.1.7    18813               South East 2.278296e-03 1.731680e-04 NA  1.938893e-03 2.617699e-03
# 26     other    18813               South East 1.844802e-05 2.443870e-06 NA  1.365813e-05 2.323792e-05
# 27  B.1.177+    18813               South East 2.364301e-09 5.919424e-10 NA  1.204116e-09 3.524487e-09
# 28   B.1.525    18813               South East 2.016280e-05 7.980177e-06 NA  4.521937e-06 3.580366e-05
# 29   B.1.351    18813               South East 7.634639e-05 1.986231e-05 NA  3.741698e-05 1.152758e-04
# 30       P.1    18813               South East 3.488857e-05 1.869638e-05 NA -1.755652e-06 7.153279e-05
# 31 B.1.617.1    18813               South East 8.660939e-08 4.593787e-08 NA -3.427186e-09 1.766460e-07
# 32 B.1.617.2    18813               South East 9.975718e-01 1.842510e-04 NA  9.972106e-01 9.979329e-01
# 33   B.1.1.7    18813          East of England 2.408574e-03 1.771156e-04 NA  2.061434e-03 2.755714e-03
# 34     other    18813          East of England 5.078026e-05 6.644017e-06 NA  3.775823e-05 6.380230e-05
# 35  B.1.177+    18813          East of England 5.237053e-09 1.313501e-09 NA  2.662639e-09 7.811467e-09
# 36   B.1.525    18813          East of England 6.631263e-06 3.063060e-06 NA  6.277766e-07 1.263475e-05
# 37   B.1.351    18813          East of England 4.520119e-05 1.238967e-05 NA  2.091788e-05 6.948450e-05
# 38       P.1    18813          East of England 1.314101e-05 7.990231e-06 NA -2.519551e-06 2.880158e-05
# 39 B.1.617.1    18813          East of England 3.169998e-08 1.796072e-08 NA -3.502377e-09 6.690233e-08
# 40 B.1.617.2    18813          East of England 9.974756e-01 1.851675e-04 NA  9.971127e-01 9.978386e-01
# 41   B.1.1.7    18813            East Midlands 4.151149e-03 2.899840e-04 NA  3.582791e-03 4.719507e-03
# 42     other    18813            East Midlands 2.918748e-04 3.693768e-05 NA  2.194783e-04 3.642714e-04
# 43  B.1.177+    18813            East Midlands 5.126980e-08 1.269717e-08 NA  2.638380e-08 7.615579e-08
# 44   B.1.525    18813            East Midlands 2.465403e-06 1.646718e-06 NA -7.621053e-07 5.692911e-06
# 45   B.1.351    18813            East Midlands 3.279779e-05 1.017072e-05 NA  1.286354e-05 5.273203e-05
# 46       P.1    18813            East Midlands 2.360292e-06 2.593646e-06 NA -2.723161e-06 7.443745e-06
# 47 B.1.617.1    18813            East Midlands 6.470771e-08 3.347685e-08 NA -9.057159e-10 1.303211e-07
# 48 B.1.617.2    18813            East Midlands 9.955192e-01 3.127572e-04 NA  9.949062e-01 9.961322e-01
# 49   B.1.1.7    18813            West Midlands 3.373648e-03 2.363506e-04 NA  2.910410e-03 3.836887e-03
# 50     other    18813            West Midlands 2.384492e-04 3.014389e-05 NA  1.793683e-04 2.975301e-04
# 51  B.1.177+    18813            West Midlands 3.624275e-08 9.003455e-09 NA  1.859631e-08 5.388920e-08
# 52   B.1.525    18813            West Midlands 1.572839e-05 6.342010e-06 NA  3.298278e-06 2.815850e-05
# 53   B.1.351    18813            West Midlands 4.680206e-05 1.311040e-05 NA  2.110615e-05 7.249798e-05
# 54       P.1    18813            West Midlands 2.441384e-06 2.676050e-06 NA -2.803577e-06 7.686344e-06
# 55 B.1.617.1    18813            West Midlands 1.689266e-07 8.120541e-08 NA  9.766947e-09 3.280863e-07
# 56 B.1.617.2    18813            West Midlands 9.963227e-01 2.572636e-04 NA  9.958185e-01 9.968270e-01
# 57   B.1.1.7    18813               North East 5.995176e-03 5.871899e-04 NA  4.844305e-03 7.146047e-03
# 58     other    18813               North East 6.660327e-04 9.633013e-05 NA  4.772291e-04 8.548363e-04
# 59  B.1.177+    18813               North East 5.907007e-08 1.524674e-08 NA  2.918701e-08 8.895312e-08
# 60   B.1.525    18813               North East 3.749864e-06 3.016991e-06 NA -2.163331e-06 9.663058e-06
# 61   B.1.351    18813               North East 3.974176e-05 1.696362e-05 NA  6.493674e-06 7.298985e-05
# 62       P.1    18813               North East 9.932360e-06 1.102504e-05 NA -1.167632e-05 3.154104e-05
# 63 B.1.617.1    18813               North East 2.279342e-08 2.547463e-08 NA -2.713594e-08 7.272279e-08
# 64 B.1.617.2    18813               North East 9.932853e-01 6.579160e-04 NA  9.919958e-01 9.945748e-01
# 65   B.1.1.7    18813 Yorkshire and The Humber 9.592254e-03 7.023372e-04 NA  8.215698e-03 1.096881e-02
# 66     other    18813 Yorkshire and The Humber 1.205835e-03 1.468447e-04 NA  9.180245e-04 1.493645e-03
# 67  B.1.177+    18813 Yorkshire and The Humber 1.630458e-07 3.995088e-08 NA  8.474353e-08 2.413481e-07
# 68   B.1.525    18813 Yorkshire and The Humber 1.125562e-05 5.260396e-06 NA  9.454334e-07 2.156581e-05
# 69   B.1.351    18813 Yorkshire and The Humber 3.694422e-05 1.245305e-05 NA  1.253669e-05 6.135175e-05
# 70       P.1    18813 Yorkshire and The Humber 3.160175e-13 1.451496e-13 NA  3.152959e-14 6.005054e-13
# 71 B.1.617.1    18813 Yorkshire and The Humber 3.184170e-08 2.119718e-08 NA -9.704013e-09 7.338742e-08
# 72 B.1.617.2    18813 Yorkshire and The Humber 9.891535e-01 7.955072e-04 NA  9.875944e-01 9.907127e-01

multinom_preds_delta_today_byregion = multinom_preds_today_byregion[multinom_preds_today_byregion$LINEAGE2=="B.1.617.2",]
multinom_preds_delta_today_byregion
# LINEAGE2 DATE_NUM                   REGION      prob           SE df asymp.LCL asymp.UCL
# 8  B.1.617.2    18813                   London 0.9980074 1.468807e-04 NA 0.9977196 0.9982953
# 16 B.1.617.2    18813               North West 0.9989026 7.404598e-05 NA 0.9987575 0.9990478
# 24 B.1.617.2    18813               South West 0.9969966 2.489434e-04 NA 0.9965087 0.9974845
# 32 B.1.617.2    18813               South East 0.9975718 1.842510e-04 NA 0.9972106 0.9979329
# 40 B.1.617.2    18813          East of England 0.9974756 1.851675e-04 NA 0.9971127 0.9978386
# 48 B.1.617.2    18813            East Midlands 0.9955192 3.127572e-04 NA 0.9949062 0.9961322
# 56 B.1.617.2    18813            West Midlands 0.9963227 2.572636e-04 NA 0.9958185 0.9968270
# 64 B.1.617.2    18813               North East 0.9932853 6.579160e-04 NA 0.9919958 0.9945748
# 72 B.1.617.2    18813 Yorkshire and The Humber 0.9891535 7.955072e-04 NA 0.9875944 0.9907127

# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME ####

gentime = 4.7 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantages of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit4_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=sanger, maxit=1000) here
emtrsanger4 = emtrends(fit4_sanger_multi, trt.vs.ctrl1 ~ LINEAGE2, by=c("DATE_NUM","REGION"), 
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=seq(as.numeric(as.Date("2021-04-01")),as.numeric(as.Date("2021-07-01")))))
delta_r_sanger4 = data.frame(confint(emtrsanger4, 
                                    adjust="none", df=NA)$contrasts)
delta_r_sanger4 = delta_r_sanger4[delta_r_sanger4$contrast=="B.1.617.2 - B.1.1.7",]
delta_r_sanger4

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger4 = delta_r_sanger4
transmadv_sanger4$estimate =  100*(exp(transmadv_sanger4$estimate*gentime)-1)
transmadv_sanger4$asymp.LCL =  100*(exp(transmadv_sanger4$asymp.LCL*gentime)-1)
transmadv_sanger4$asymp.UCL =  100*(exp(transmadv_sanger4$asymp.UCL*gentime)-1)
transmadv_sanger4$collection_date = as.Date(transmadv_sanger4$DATE_NUM, origin="1970-01-01")
transmadv_sanger4$REGION = factor(transmadv_sanger4$REGION, levels=levels_REGION)

plot_sanger_mfit4_transmadv = qplot(data=transmadv_sanger4, 
                               x=collection_date, y=estimate, ymin=asymp.LCL, ymax=asymp.UCL, colour=REGION, fill=REGION, group=REGION,
                               geom="blank") +
  facet_wrap(~ REGION) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2021-04-01",NA)), expand=c(0,0)) + #  
  scale_y_continuous(limits=c(0,max(transmadv_sanger4$asymp.UCL))) +
  geom_ribbon(aes(fill=I("steelblue"), colour=NULL), alpha=I(0.3)) +
  geom_line(aes(colour=I("steelblue"), fill=NULL)) +
  ylab("Transmission advantage of B.1.617.2 over B.1.1.7 (%)") +
  theme_hc() + xlab("") +
  ggtitle("TRANSMISSION ADVANTAGE OF B.1.617.2 OVER B.1.1.7\n(Sanger Institute baseline surveillance data, multinomial fit)") +
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
plot_sanger_mfit4_transmadv

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit3_transm advantage B16172.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit3_transm advantage B16172.pdf"), width=8, height=6)




# PLOT MULTINOMIAL FIT ####

# extrapolate = 60
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-07-31")) # max(sanger$DATE_NUM)+extrapolate

# predictions by ONS region ####

fit_sanger_multi_predsbyregion = data.frame(emmeans(fit4_sanger_multi, 
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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit.pdf"), width=8, height=6)

library(ggpubr)
ggarrange(muller_sangerbyregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbyregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.png"), width=8, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.pdf"), width=8, height=6)


# predictions by NHS region

fit_sanger_multi_predsbynhsregion = data.frame(emmeans(fit4_sanger_nhs_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "NHSREGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                                  mode="prob", df=NA))
fit_sanger_multi_predsbynhsregion$collection_date = as.Date(fit_sanger_multi_predsbynhsregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbynhsregion$LINEAGE2 = factor(fit_sanger_multi_predsbynhsregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_sanger_multi_predsbynhsregion$NHSREGION = factor(fit_sanger_multi_predsbynhsregion$NHSREGION, levels=levels_NHSREGION) 
# predicted incidence in different parts of England today
# fit_sanger_multi_predsbynhsregion[fit_sanger_multi_predsbynhsregion$collection_date==today,]

fit_sanger_multi_preds = data.frame(emmeans(fit4_sanger_nhs_multi, 
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots_multinom fit multipanel.png"), width=8, height=8)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots_multinom fit multipanel.pdf"), width=8, height=6)


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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit.pdf"), width=8, height=6)

library(ggpubr)
ggarrange(muller_sangerbynhsregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbynhsregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit multipanel.png"), width=8, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit multipanel.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & 95% CONFIDENCE INTERVALS

# PHE SGTF data to superimpose on model predictions as independent estimate of prop of B.1.1.7 (these data are not up to date anymore)
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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by ONS region_stacked as proportion.pdf"), width=10, height=8)

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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.pdf"), width=10, height=8)



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
# 1 0.9973022 0.9966437 0.9978318

# here given by NHS region:
data.frame(REGION=data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                   at=list(DATE_NUM=today_num), 
                   type="response"))[,2], 
           1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                        at=list(DATE_NUM=today_num), 
                        type="response"))[,c(3,7,6)])
#                     REGION      prob asymp.UCL asymp.LCL
# 1                   London 0.9744834 0.9649967 0.9814484
# 2               North West 0.9997466 0.9996810 0.9997987
# 3               South West 0.9981229 0.9909365 0.9996135
# 4               South East 0.9980492 0.9969066 0.9987702
# 5          East of England 0.9968837 0.9954483 0.9978675
# 6            East Midlands 0.9954652 0.9935775 0.9967998
# 7            West Midlands 0.9967837 0.9954200 0.9977423
# 8               North East 0.9970730 0.9942716 0.9985064
# 9 Yorkshire and The Humber 0.9976482 0.9964750 0.9984316




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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.pdf"), width=8, height=6)
library(svglite)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.svg"), width=8, height=6)

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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_response scale.pdf"), width=8, height=6)


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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.pdf"), width=8, height=6)
library(svglite)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.svg"), width=8, height=6)
# write.csv(fit_sanger_multi_predsbyregion2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_multinomial fit sanger inst data.csv"), row.names=F)
# write.csv(SGTF_predsbyregion, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logistic fit S dropout data.csv"), row.names=F)
# write.csv(data_agbyweekregion1, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_sanger inst data aggregated by week.csv"), row.names=F)
# write.csv(sgtf_UK2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_S dropout data.csv"), row.names=F)

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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_response scale.pdf"), width=8, height=6)


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
cases_uk_region = cases_uk_region[cases_uk_region$date<=(max(cases_uk_region$date)-3),] # cut off data from last 3 days (incomplete)

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
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/07/Covid-Publication-08-07-2021-Supplementary-Data.xlsx
# PS Next publication: Thursday 8 July 2021
library(rio)
dat = rio::import(file = "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/07/Covid-Publication-08-07-2021-Supplementary-Data.xlsx",
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
# ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.pdf"), width=8, height=6)
# write.csv(cases_hosps_ENG_age, file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.csv"), row.names=F)


# MAP MULTINOMIAL FIT ONTO CASE NUMBERS ####

newdat = expand.grid(REGION=levels_REGION,
                     DATE_NUM=seq(date.from, date.to))
fit_sanger_multi_predsbyregion_day = data.frame(newdat,
                                          predict(fit4_sanger_multi, 
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
gamfittotcases = gam(totcases ~ s(DATE_NUM, bs="cs", k=7, by=REGION) + REGION + WEEKDAY, family=poisson(log), 
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

fit_sanger_multi_predsbyregion_day$cases[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.617.2"&
                                           fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA
fit_sanger_multi_predsbyregion_day$cases_smoothed[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.617.2"&
                                           fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA
fit_sanger_multi_predsbyregion_day$cases_rollingmean[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.617.2"&
                                                    fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA

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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases+1, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases_smoothed+1, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit smoothed case data_log10 Y scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases_rollingmean, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases_rollingmean+1, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND\n(blue=B.1.1.7 (alpha), magenta=B.1.617.2 (delta))") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data_log10 Y scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

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
fit_cases = gam(totcases ~ s(DATE_NUM, bs="cs", k=8, by=REGION) + REGION + weekday,
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
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re VALUES IN ENGLAND BY REGION BASED ON NEW CONFIRMED CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
  # coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants4 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
       function (dat) { 
  wt = as.data.frame(emmeans(fit4_sanger_multi, ~ LINEAGE2, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
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

# # for simpler model fit1_sanger_multi
# above_avg_r_variants1 = do.call(rbind, lapply(seq(date.from,
#                                                   date.to), 
#                                               function (dat) { 
#                                                 wt = as.data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
#                                                 # wt = rep(1/length(levels_LINEAGE2), length(levels_LINEAGE2)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE2)
#                                                 cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
#                                                 names(cons) = seq_along(cons)
#                                                 EMT = emtrends(fit1_sanger_multi,  ~ LINEAGE2, by=c("DATE_NUM","REGION"),
#                                                                var="DATE_NUM", mode="latent",
#                                                                at=list(DATE_NUM=dat,
#                                                                        REGION=levels(SGTF_predsbyregion$REGION)))
#                                                 out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
#                                                 # sum(out$estimate*wt) # should sum to zero
#                                                 return(out) } ))
above_avg_r_variants = above_avg_r_variants4
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE2), 
                                       labels=levels(data_agbyweeknhsregion1$LINEAGE2))
above_avg_r_variants$LINEAGE2 = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$DATE_NUM = above_avg_r_variants$DATE_NUM-0 # the -7 to calculate back to time of infection
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2020-09-01" "2021-07-31"
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
ymax = 2
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
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$LINEAGE2 %in% c("other"))|(above_avg_r_variants2$collection_date>max(cases_uk_region$date))),], 
      x=collection_date-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE2, fill=LINEAGE2, alpha=I(0.5),
      group=LINEAGE2, linetype=I(0)) +
  facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=LINEAGE2, colour=LINEAGE2), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE2), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN ENGLAND BY REGION\n(case data gov.uk & multinomial fit to Sanger Institute sequence data)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),max(cases_uk_region$date))) +
  scale_fill_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  scale_colour_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  theme(legend.position="right") # ,  
# axis.title.x=element_blank()

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant_with clipping.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant_with clipping.pdf"), width=8, height=6)


