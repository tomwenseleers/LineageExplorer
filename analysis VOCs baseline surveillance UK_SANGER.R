# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN ENGLAND BASED ON SANGER INSTITUTE BASELINE SURVEILLANCE SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & most surge testing/active surveillance)

# last update 18 OCTOBER 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
# today = as.Date("2021-08-26")
today # "2021-10-18"
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
range(sanger$WeekEndDate) # "2020-09-05" "2021-10-09"
# sanger$Week = lubridate::week(sanger$WeekEndDate)
sanger$DATE_NUM = as.numeric(sanger$WeekEndDate)-3.5 # using week midpoint
colnames(sanger)

sanger = sanger[rep(seq_len(nrow(sanger)), sanger$Count),] # convert to long format
sanger$Count = NULL
nrow(sanger) # 601759

nrow(sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 77504 (last 2 weeks)
nrow(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 11186 (last 2 weeks)
nrow(sanger[grepl("AY.3",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 708 (last 2 weeks)
nrow(sanger[grepl("AY.4",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 60059 (last 2 weeks)
nrow(sanger[grepl("AY.4.2",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 6319 (last 2 weeks)
nrow(sanger[grepl("B.1.1.7",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 0

length(unique(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),"LTLA"])) # 310 LTLAs
length(unique(sanger[grepl("AY.4.2",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),"LTLA"])) # 297 LTLAs

unique(sanger$Lineage)
sel_target_VOC = "B.1.617.2"
table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
table(sanger$Lineage[grepl("AY.",sanger$Lineage, fixed=TRUE)]) 
# main ones: AY.10, AY.11, AY.3, AY.34, AY.36, AY.4, AY.4.1, AY.4.2, AY.4.3, AY.4.5, AY.5, AY.6, AY.7, AY.8, AY.9

unique(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
# "B.1.617.1" "B.1.617.3" "B.1.617.2"
sanger$LINEAGE = sanger$Lineage
# sanger[grepl(sel_target_VOC, sanger$LINEAGE1, fixed=T),"LINEAGE"] = paste0(sel_target_VOC,"+")
sanger[grepl("B.1.177", sanger$LINEAGE, fixed=T),"LINEAGE"] = "B.1.177+"

table_lineage = as.data.frame(table(sanger$LINEAGE))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger$LINEAGE))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# just for the last week worth of data
table_lineage = as.data.frame(table(sanger[sanger$WeekEndDate==max(sanger$WeekEndDate),]$LINEAGE))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger[sanger$WeekEndDate==max(sanger$WeekEndDate),]$LINEAGE))
colnames(table_lineage) = c("Lineage","Prop")
tbl = table_lineage[table_lineage$Prop>0.00001,]
tbl = tbl[order(tbl$Prop, decreasing=T),]
tbl
#      Lineage         Prop
# 7       AY.4 6.768393e-01
# 20 B.1.617.2 1.475225e-01
# 9     AY.4.2 9.241742e-02
# 14      AY.5 3.460961e-02
# 15      AY.6 2.271021e-02
# 19      AY.9 7.732733e-03
# 4      AY.34 5.442943e-03
# 5      AY.36 4.354354e-03
# 16      AY.7 3.078078e-03
# 2       AY.3 1.238739e-03
# 1      AY.25 9.384384e-04
# 8     AY.4.1 8.258258e-04
# 10    AY.4.3 7.882883e-04
# 12    AY.4.5 5.630631e-04
# 3      AY.33 4.129129e-04
# 11    AY.4.4 1.876877e-04
# 13     AY.41 1.501502e-04
# 17    AY.7.1 1.126126e-04
# 6      AY.39 3.753754e-05
# 18    AY.7.2 3.753754e-05

table(sanger$Lineage)[order(table(sanger$Lineage), decreasing=T)]

sel_ref_lineage = "B.1.617.2"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"Lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.177+","B.1.1.7","B.1.617.2", # "B.1.351", "B.1.525","P.1","B.1.621","C.37"
                 rev(c("AY.4.2","AY.4","AY.5","AY.6","AY.9","AY.34","AY.36","AY.7","AY.3")))

sanger$LINEAGE[!(sanger$LINEAGE %in% sel_lineages)] = "other"
# sanger = sanger[sanger$LINEAGE1 %in% sel_lineages, ]
sum(table(sanger$LINEAGE))
table(sanger$LINEAGE)

levels_LINEAGE = c("other",sel_lineages)
sanger$LINEAGE = factor(sanger$LINEAGE, levels=levels_LINEAGE)

sanger_lastmonth = sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-30),]
table(sanger_lastmonth$Lineage, sanger_lastmonth$REGION)

str(sanger)

# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of England
data_agbyweek1 = as.data.frame(table(sanger$WeekEndDate, sanger$LINEAGE))
colnames(data_agbyweek1) = c("WeekEndDate", "LINEAGE", "count")
data_agbyweek1_sum = aggregate(count ~ WeekEndDate, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$WeekEndDate, data_agbyweek1_sum$WeekEndDate)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE=="B.1.617.2","total"]) == nrow(sanger) # TRUE
data_agbyweek1$WeekEndDate = as.Date(as.character(data_agbyweek1$WeekEndDate))
data_agbyweek1$collection_date = data_agbyweek1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweek1$LINEAGE = factor(data_agbyweek1$LINEAGE, levels=levels_LINEAGE)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and ONS region
data_agbyweekregion1 = as.data.frame(table(sanger$WeekEndDate, sanger$REGION, sanger$LINEAGE))
colnames(data_agbyweekregion1) = c("WeekEndDate", "REGION", "LINEAGE", "count")
data_agbyweekregion1_sum = aggregate(count ~ WeekEndDate + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$WeekEndDate,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$WeekEndDate,data_agbyweekregion1_sum$REGION))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE=="B.1.617.2","total"]) == nrow(sanger) # TRUE
data_agbyweekregion1$WeekEndDate = as.Date(as.character(data_agbyweekregion1$WeekEndDate))
data_agbyweekregion1$collection_date = data_agbyweekregion1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweekregion1$LINEAGE = factor(data_agbyweekregion1$LINEAGE, levels=levels_LINEAGE)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion1[data_agbyweekregion1$collection_date==max(data_agbyweekregion1$collection_date),]

# aggregated by week and NHS region
data_agbyweeknhsregion1 = as.data.frame(table(sanger$WeekEndDate, sanger$NHSREGION, sanger$LINEAGE))
colnames(data_agbyweeknhsregion1) = c("WeekEndDate", "NHSREGION", "LINEAGE", "count")
data_agbyweeknhsregion1_sum = aggregate(count ~ WeekEndDate + NHSREGION, data=data_agbyweeknhsregion1, sum)
data_agbyweeknhsregion1$total = data_agbyweeknhsregion1_sum$count[match(interaction(data_agbyweeknhsregion1$WeekEndDate,data_agbyweeknhsregion1$NHSREGION), 
                                                                  interaction(data_agbyweeknhsregion1_sum$WeekEndDate,data_agbyweeknhsregion1_sum$NHSREGION))]
sum(data_agbyweeknhsregion1[data_agbyweeknhsregion1$LINEAGE=="B.1.617.2","total"]) == nrow(sanger) # TRUE
data_agbyweeknhsregion1$WeekEndDate = as.Date(as.character(data_agbyweeknhsregion1$WeekEndDate))
data_agbyweeknhsregion1$collection_date = data_agbyweeknhsregion1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweeknhsregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweeknhsregion1$LINEAGE = factor(data_agbyweeknhsregion1$LINEAGE, levels=levels_LINEAGE)
data_agbyweeknhsregion1$NHSREGION = factor(data_agbyweeknhsregion1$NHSREGION, levels=levels_NHSREGION)
data_agbyweeknhsregion1$collection_date_num = as.numeric(data_agbyweeknhsregion1$collection_date)
data_agbyweeknhsregion1$prop = data_agbyweeknhsregion1$count/data_agbyweeknhsregion1$total
data_agbyweeknhsregion1 = data_agbyweeknhsregion1[data_agbyweeknhsregion1$total!=0,]

data_agbyweeknhsregion1[data_agbyweeknhsregion1$collection_date==max(data_agbyweeknhsregion1$collection_date),]



# MULLER PLOT (RAW DATA)
unique(sanger$LINEAGE)
levels_LINEAGE_plot = levels_LINEAGE # c("other","B.1.177+","B.1.525","B.1.351","P.1","B.1.1.7","B.1.617.1","B.1.617.2","AY.3","B.1.621")

library(scales)
n1 = length(levels_LINEAGE_plot)
lineage_cols1 = c(hcl(h = seq(0, 260, length = n1-1), l = 60, c = 200), muted("dodgerblue", l=55, c=120)) # grey70
# lineage_cols1 = c(hcl(h = seq(0, 265, length = n1), l = 60, c = 200)) # grey70
# lineage_cols1[which(levels_LINEAGE_plot=="B.1.617+")] = "magenta" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.1")] = muted("magenta") # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE_plot=="B.1.621")] = "limegreen" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")] = "#0085FF"  
lineage_cols1[which(levels_LINEAGE_plot=="other")] = "grey70"  
lineage_cols1[which(levels_LINEAGE_plot=="B.1.177+")] = "grey55"  
lineage_cols1[which(levels_LINEAGE_plot=="B.1.351")] = muted("cyan")  
lineage_cols1[which(levels_LINEAGE_plot=="P.1")] = "cyan3"  
library(colorspace)
AY_cols = colorRampPalette(c("blue", "orange", "red3"))(9) # rainbow_hcl(n=11, c=200, l=45) # 
lineage_cols1[which(levels_LINEAGE_plot=="AY.3")] = AY_cols[1]
lineage_cols1[which(levels_LINEAGE_plot=="AY.7")] = AY_cols[2]  
lineage_cols1[which(levels_LINEAGE_plot=="AY.36")] = AY_cols[3]
lineage_cols1[which(levels_LINEAGE_plot=="AY.34")] = AY_cols[4]  
lineage_cols1[which(levels_LINEAGE_plot=="AY.9")] = AY_cols[5] 
lineage_cols1[which(levels_LINEAGE_plot=="AY.6")] = AY_cols[6] 
lineage_cols1[which(levels_LINEAGE_plot=="AY.5")] = AY_cols[7]  
lineage_cols1[which(levels_LINEAGE_plot=="AY.4")] = AY_cols[8]  
lineage_cols1[which(levels_LINEAGE_plot=="AY.4.2")] = muted(AY_cols[9])  

data_agbyweek1$LINEAGE = factor(data_agbyweek1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweekregion1$LINEAGE = factor(data_agbyweekregion1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweeknhsregion1$LINEAGE = factor(data_agbyweeknhsregion1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweeknhsregion1$NHSREGION = factor(data_agbyweeknhsregion1$NHSREGION, levels=levels_NHSREGION)


library(ggplot2)
library(ggthemes)
muller_sanger_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sanger_raw1


muller_sangerbyregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)") +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbyregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_raw data.pdf"), width=8, height=6)


muller_sangerbynhsregion_raw1 = ggplot(data=data_agbyweeknhsregion1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  facet_wrap(~ NHSREGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)") +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbynhsregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
sanger$LINEAGE = relevel(sanger$LINEAGE, ref="B.1.617.2")
# sanger$DATE_NUM_WEEK = sanger$DATE_NUM/7

data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
data_agbyweekregion1$LINEAGE = relevel(data_agbyweekregion1$LINEAGE, ref="B.1.617.2")
data_agbyweeknhsregion1$DATE_NUM = data_agbyweeknhsregion1$collection_date_num
data_agbyweeknhsregion1$LINEAGE = relevel(data_agbyweeknhsregion1$LINEAGE, ref="B.1.617.2")

# by ONS region
set.seed(1)
fit1_sanger_multi = nnet::multinom(LINEAGE ~ REGION + DATE_NUM, data=data_agbyweekregion1, weights=count, maxit=1000)
fit2_sanger_multi = nnet::multinom(LINEAGE ~ REGION * DATE_NUM, data=data_agbyweekregion1, weights=count, maxit=1000)
fit3_sanger_multi = nnet::multinom(LINEAGE ~ REGION + ns(DATE_NUM, df=2), data=data_agbyweekregion1, weights=count, maxit=1000)
fit4_sanger_multi = nnet::multinom(LINEAGE ~ REGION * ns(DATE_NUM, df=2), data=data_agbyweekregion1, weights=count, maxit=1000)
fit5_sanger_multi = nnet::multinom(LINEAGE ~ REGION + ns(DATE_NUM, df=3), data=data_agbyweekregion1, weights=count, maxit=1000)
fit6_sanger_multi = nnet::multinom(LINEAGE ~ REGION * ns(DATE_NUM, df=3), data=data_agbyweekregion1, weights=count, maxit=1000)
BIC(fit1_sanger_multi, fit2_sanger_multi, fit3_sanger_multi, fit4_sanger_multi, fit5_sanger_multi, fit6_sanger_multi) 
# fit6_sanger_multi fits best (lowest BIC)

# by NHS region
set.seed(1)
fit1_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION + DATE_NUM, data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit2_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION * DATE_NUM, data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit3_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION + ns(DATE_NUM, df=2), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit4_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION * ns(DATE_NUM, df=2), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit5_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION + ns(DATE_NUM, df=3), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit6_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION * ns(DATE_NUM, df=3), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
BIC(fit1_sanger_nhs_multi, fit2_sanger_nhs_multi, fit3_sanger_nhs_multi, fit4_sanger_nhs_multi, fit5_sanger_nhs_multi, fit6_sanger_nhs_multi) 
# fit6_sanger_multi fits best (lowest BIC)

# growth rate advantage of different lineages over B.1.617.2 based on model fit6_sanger_multi
emtrsanger6 = emtrends(fit6_sanger_multi, trt.vs.ctrl1 ~ LINEAGE,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger6 = data.frame(confint(emtrsanger6, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger6$contrasts)$p.value)
delta_r_sanger6
#                  contrast      estimate           SE df    asymp.LCL    asymp.UCL      p.value
# 1       other - B.1.617.2 -0.0057820585 0.0010565826 NA -0.007852922 -0.003711195 8.994271e-07
# 2  (B.1.177+) - B.1.617.2 -1.0421383337 0.0537443058 NA -1.147475238 -0.936801430 0.000000e+00
# 3     B.1.1.7 - B.1.617.2 -0.0976570294 0.0038168732 NA -0.105137963 -0.090176095 0.000000e+00
# 4        AY.3 - B.1.617.2 -0.1691843193 0.0148133316 NA -0.198217916 -0.140150723 0.000000e+00
# 5        AY.7 - B.1.617.2 -0.0393355084 0.0018962777 NA -0.043052144 -0.035618873 0.000000e+00
# 6       AY.36 - B.1.617.2 -0.0267934139 0.0072477481 NA -0.040998739 -0.012588089 2.731606e-03
# 7       AY.34 - B.1.617.2 -0.0365927216 0.0078062399 NA -0.051892671 -0.021292773 4.367935e-05
# 8        AY.9 - B.1.617.2 -0.0371259913 0.0014329111 NA -0.039934445 -0.034317537 0.000000e+00
# 9        AY.6 - B.1.617.2 -0.0312333548 0.0010349504 NA -0.033261820 -0.029204889 0.000000e+00
# 10       AY.5 - B.1.617.2 -0.0209473437 0.0007740163 NA -0.022464388 -0.019430300 0.000000e+00
# 11       AY.4 - B.1.617.2 -0.0213013040 0.0003148839 NA -0.021918465 -0.020684143 0.000000e+00
# 12     AY.4.2 - B.1.617.2  0.0006794459 0.0009690568 NA -0.001219871  0.002578762 9.722237e-01

# pairwise growth rate advantages for all strain comparisons (i.e. pairwise differences in growth rate per day among the different lineages)
emtrsanger_pairw6 = emtrends(fit6_sanger_multi, pairwise ~ LINEAGE,  
                            var="DATE_NUM",  mode="latent",
                            at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_pairw6 = data.frame(confint(emtrsanger_pairw6, 
                                          adjust="none", df=NA)$contrasts, 
                                  p.value=as.data.frame(emtrsanger_pairw6$contrasts)$p.value)
delta_r_sanger_pairw6

delta_r_sanger_pairw6[delta_r_sanger_pairw6$contrast %in% c("AY.4.2 - B.1.617.2",
                                                            "AY.4 - B.1.617.2")]
#                  contrast      estimate           SE df    asymp.LCL     asymp.UCL      p.value
# 1       B.1.617.2 - other  0.0057820585 0.0010565826 NA  0.003711195  7.852922e-03 5.772407e-06
# 2  B.1.617.2 - (B.1.177+)  1.0421383337 0.0537443058 NA  0.936801430  1.147475e+00 0.000000e+00
# 3     B.1.617.2 - B.1.1.7  0.0976570294 0.0038168732 NA  0.090176095  1.051380e-01 0.000000e+00
# 4        B.1.617.2 - AY.3  0.1691843193 0.0148133316 NA  0.140150723  1.982179e-01 0.000000e+00
# 5        B.1.617.2 - AY.7  0.0393355084 0.0018962777 NA  0.035618873  4.305214e-02 0.000000e+00
# 6       B.1.617.2 - AY.36  0.0267934139 0.0072477481 NA  0.012588089  4.099874e-02 1.508839e-02
# 7       B.1.617.2 - AY.34  0.0365927216 0.0078062399 NA  0.021292773  5.189267e-02 2.720425e-04
# 8        B.1.617.2 - AY.9  0.0371259913 0.0014329111 NA  0.034317537  3.993445e-02 0.000000e+00
# 9        B.1.617.2 - AY.6  0.0312333548 0.0010349504 NA  0.029204889  3.326182e-02 0.000000e+00
# 10       B.1.617.2 - AY.5  0.0209473437 0.0007740163 NA  0.019430300  2.246439e-02 0.000000e+00
# 11       B.1.617.2 - AY.4  0.0213013040 0.0003148839 NA  0.020684143  2.191847e-02 0.000000e+00
# 12     B.1.617.2 - AY.4.2 -0.0006794459 0.0009690568 NA -0.002578762  1.219871e-03 9.999680e-01
# 13     other - (B.1.177+)  1.0363562752 0.0537475271 NA  0.931013058  1.141699e+00 0.000000e+00
# 14        other - B.1.1.7  0.0918749709 0.0033401328 NA  0.085328431  9.842151e-02 0.000000e+00
# 15           other - AY.3  0.1634022608 0.0148487134 NA  0.134299317  1.925052e-01 0.000000e+00
# 16           other - AY.7  0.0335534499 0.0021508405 NA  0.029337880  3.776902e-02 0.000000e+00
# 17          other - AY.36  0.0210113554 0.0073198915 NA  0.006664632  3.535808e-02 1.749716e-01
# 18          other - AY.34  0.0308106631 0.0078732287 NA  0.015379418  4.624191e-02 6.867266e-03
# 19           other - AY.9  0.0313439328 0.0017558380 NA  0.027902554  3.478531e-02 0.000000e+00
# 20           other - AY.6  0.0254512963 0.0014549950 NA  0.022599558  2.830303e-02 0.000000e+00
# 21           other - AY.5  0.0151652852 0.0012801492 NA  0.012656239  1.767433e-02 0.000000e+00
# 22           other - AY.4  0.0155192455 0.0010689795 NA  0.013424084  1.761441e-02 0.000000e+00
# 23         other - AY.4.2 -0.0064615044 0.0014087103 NA -0.009222526 -3.700483e-03 4.285390e-04
# 24   (B.1.177+) - B.1.1.7 -0.9444813043 0.0543773728 NA -1.051058997 -8.379036e-01 0.000000e+00
# 25      (B.1.177+) - AY.3 -0.8729540144 0.0557478223 NA -0.982217738 -7.636903e-01 0.000000e+00
# 26      (B.1.177+) - AY.7 -1.0028028253 0.0537764528 NA -1.108202736 -8.974029e-01 0.000000e+00
# 27     (B.1.177+) - AY.36 -1.0153449198 0.0542297818 NA -1.121633339 -9.090565e-01 0.000000e+00
# 28     (B.1.177+) - AY.34 -1.0055456122 0.0543071626 NA -1.111985695 -8.991055e-01 0.000000e+00
# 29      (B.1.177+) - AY.9 -1.0050123424 0.0537622853 NA -1.110384485 -8.996402e-01 0.000000e+00
# 30      (B.1.177+) - AY.6 -1.0109049790 0.0537530915 NA -1.116259102 -9.055509e-01 0.000000e+00
# 31      (B.1.177+) - AY.5 -1.0211909900 0.0537485603 NA -1.126536232 -9.158457e-01 0.000000e+00
# 32      (B.1.177+) - AY.4 -1.0208370297 0.0537439432 NA -1.126173223 -9.155008e-01 0.000000e+00
# 33    (B.1.177+) - AY.4.2 -1.0428177796 0.0537519309 NA -1.148169628 -9.374659e-01 0.000000e+00
# 34         B.1.1.7 - AY.3  0.0715272899 0.0152967607 NA  0.041546190  1.015084e-01 2.868510e-04
# 35         B.1.1.7 - AY.7 -0.0583215210 0.0042621727 NA -0.066675226 -4.996782e-02 0.000000e+00
# 36        B.1.1.7 - AY.36 -0.0708636155 0.0081900775 NA -0.086915872 -5.481136e-02 0.000000e+00
# 37        B.1.1.7 - AY.34 -0.0610643078 0.0086885930 NA -0.078093637 -4.403498e-02 6.140310e-10
# 38         B.1.1.7 - AY.9 -0.0605310381 0.0040781837 NA -0.068524131 -5.253794e-02 0.000000e+00
# 39         B.1.1.7 - AY.6 -0.0664236746 0.0039534148 NA -0.074172225 -5.867512e-02 0.000000e+00
# 40         B.1.1.7 - AY.5 -0.0767096856 0.0038937681 NA -0.084341331 -6.907804e-02 0.000000e+00
# 41         B.1.1.7 - AY.4 -0.0763557254 0.0038301867 NA -0.083862753 -6.884870e-02 0.000000e+00
# 42       B.1.1.7 - AY.4.2 -0.0983364753 0.0039353097 NA -0.106049541 -9.062341e-02 0.000000e+00
# 43            AY.3 - AY.7 -0.1298488109 0.0149293858 NA -0.159109869 -1.005878e-01 0.000000e+00
# 44           AY.3 - AY.36 -0.1423909054 0.0164866577 NA -0.174704161 -1.100777e-01 0.000000e+00
# 45           AY.3 - AY.34 -0.1325915978 0.0167386386 NA -0.165398727 -9.978447e-02 0.000000e+00
# 46            AY.3 - AY.9 -0.1320583280 0.0148779509 NA -0.161218576 -1.028981e-01 0.000000e+00
# 47            AY.3 - AY.6 -0.1379509646 0.0148439566 NA -0.167044585 -1.088573e-01 0.000000e+00
# 48            AY.3 - AY.5 -0.1482369756 0.0148283535 NA -0.177300014 -1.191739e-01 0.000000e+00
# 49            AY.3 - AY.4 -0.1478830153 0.0148114074 NA -0.176912840 -1.188532e-01 0.000000e+00
# 50          AY.3 - AY.4.2 -0.1698637652 0.0148394476 NA -0.198948548 -1.407790e-01 0.000000e+00
# 51           AY.7 - AY.36 -0.0125420945 0.0074831910 NA -0.027208879  2.124690e-03 9.033021e-01
# 52           AY.7 - AY.34 -0.0027427869 0.0080246919 NA -0.018470894  1.298532e-02 1.000000e+00
# 53            AY.7 - AY.9 -0.0022095171 0.0023422739 NA -0.006800290  2.381255e-03 9.992825e-01
# 54            AY.7 - AY.6 -0.0081021537 0.0021228082 NA -0.012262781 -3.941526e-03 9.823218e-03
# 55            AY.7 - AY.5 -0.0183881647 0.0020093159 NA -0.022326351 -1.444998e-02 0.000000e+00
# 56            AY.7 - AY.4 -0.0180342044 0.0018809675 NA -0.021720833 -1.434758e-02 0.000000e+00
# 57          AY.7 - AY.4.2 -0.0400149543 0.0020961252 NA -0.044123284 -3.590662e-02 0.000000e+00
# 58          AY.36 - AY.34  0.0097993077 0.0106431372 NA -0.011060858  3.065947e-02 9.994379e-01
# 59           AY.36 - AY.9  0.0103325774 0.0073794200 NA -0.004130820  2.479597e-02 9.739243e-01
# 60           AY.36 - AY.6  0.0044399409 0.0073112931 NA -0.009889930  1.876981e-02 9.999934e-01
# 61           AY.36 - AY.5 -0.0058460702 0.0072794203 NA -0.020113472  8.421331e-03 9.998631e-01
# 62           AY.36 - AY.4 -0.0054921099 0.0072445620 NA -0.019691191  8.706971e-03 9.999258e-01
# 63         AY.36 - AY.4.2 -0.0274728598 0.0073002410 NA -0.041781069 -1.316465e-02 1.191824e-02
# 64           AY.34 - AY.9  0.0005332698 0.0079280748 NA -0.015005471  1.607201e-02 1.000000e+00
# 65           AY.34 - AY.6 -0.0053593668 0.0078644875 NA -0.020773479  1.005475e-02 9.999766e-01
# 66           AY.34 - AY.5 -0.0156453778 0.0078348466 NA -0.031001395 -2.893606e-04 7.345376e-01
# 67           AY.34 - AY.4 -0.0152914176 0.0078026107 NA -0.030584253  1.418325e-06 7.584076e-01
# 68         AY.34 - AY.4.2 -0.0372721674 0.0078536779 NA -0.052665093 -2.187924e-02 2.083492e-04
# 69            AY.9 - AY.6 -0.0058926366 0.0017225663 NA -0.009268804 -2.516469e-03 3.789442e-02
# 70            AY.9 - AY.5 -0.0161786476 0.0015801027 NA -0.019275592 -1.308170e-02 0.000000e+00
# 71            AY.9 - AY.4 -0.0158246873 0.0014133654 NA -0.018594833 -1.305454e-02 0.000000e+00
# 72          AY.9 - AY.4.2 -0.0378054372 0.0016896983 NA -0.041117185 -3.449369e-02 0.000000e+00
# 73            AY.6 - AY.5 -0.0102860110 0.0012262821 NA -0.012689480 -7.882542e-03 0.000000e+00
# 74            AY.6 - AY.4 -0.0099320508 0.0010022506 NA -0.011896426 -7.967676e-03 0.000000e+00
# 75          AY.6 - AY.4.2 -0.0319128006 0.0013593580 NA -0.034577093 -2.924851e-02 0.000000e+00
# 76            AY.5 - AY.4  0.0003539603 0.0007329955 NA -0.001082685  1.790605e-03 9.999995e-01
# 77          AY.5 - AY.4.2 -0.0216267896 0.0011739813 NA -0.023927751 -1.932583e-02 0.000000e+00
# 78          AY.4 - AY.4.2 -0.0219807499 0.0009380731 NA -0.023819339 -2.014216e-02 0.000000e+00


# growth rate advantages of different VOCs compared to B.1.617.2 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrsanger = emtrends(fit4_sanger_multi, trt.vs.ctrl1 ~ LINEAGE,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger = data.frame(confint(emtrsanger, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger$contrasts)$p.value)
delta_r_sanger
#                  contrast      estimate           SE df     asymp.LCL    asymp.UCL      p.value
# 1       other - B.1.617.2  0.0001912504 0.0007388665 NA -0.0012569014  0.001639402 9.995818e-01
# 2  (B.1.177+) - B.1.617.2 -0.1945106995 0.0048361343 NA -0.2039893485 -0.185032051 7.169820e-13
# 3     B.1.1.7 - B.1.617.2 -0.1361957569 0.0015829803 NA -0.1392983413 -0.133093172 7.169820e-13
# 4        AY.3 - B.1.617.2 -0.0292885882 0.0060468603 NA -0.0411402167 -0.017436960 2.335912e-05
# 5        AY.7 - B.1.617.2 -0.0388530323 0.0015351588 NA -0.0418618882 -0.035844176 7.169820e-13
# 6       AY.36 - B.1.617.2 -0.0023610228 0.0049589784 NA -0.0120804419  0.007358396 9.942810e-01
# 7       AY.34 - B.1.617.2  0.0016676927 0.0036592347 NA -0.0055042754  0.008839661 9.952415e-01
# 8        AY.9 - B.1.617.2 -0.0368181517 0.0010563291 NA -0.0388885188 -0.034747785 7.169820e-13
# 9        AY.6 - B.1.617.2 -0.0307251749 0.0007964834 NA -0.0322862537 -0.029164096 7.169820e-13
# 10       AY.5 - B.1.617.2 -0.0200640261 0.0006155078 NA -0.0212703992 -0.018857653 7.169820e-13
# 11       AY.4 - B.1.617.2 -0.0203320179 0.0003070369 NA -0.0209337992 -0.019730237 7.169820e-13
# 12     AY.4.2 - B.1.617.2  0.0019324240 0.0008503946 NA  0.0002656812  0.003599167 1.844430e-01

# If we take the exponent of the product of these growth rate advantages/disadvantages 
# and the generation time (e.g. 4.7 days, Nishiura et al 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*5.5)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
#                  contrast     estimate    asymp.LCL    asymp.UCL
# 1       other - B.1.617.2    0.1052430   -0.6936907    0.9057485
# 2  (B.1.177+) - B.1.617.2 -191.4822271 -207.0810143 -176.6758111
# 3     B.1.1.7 - B.1.617.2 -111.5046214 -115.1447512 -107.9260806
# 4        AY.3 - B.1.617.2  -17.4787447  -25.3915671  -10.0652602
# 5        AY.7 - B.1.617.2  -23.8240818  -25.8902595  -21.7918154
# 6       AY.36 - B.1.617.2   -1.3070305   -6.8699438    4.1301298
# 7       AY.34 - B.1.617.2    0.9214504   -3.0736417    4.9819385
# 8        AY.9 - B.1.617.2  -22.4459882  -23.8482517  -21.0596018
# 9        AY.6 - B.1.617.2  -18.4106476  -19.4316906  -17.3983338
# 10       AY.5 - B.1.617.2  -11.6671230  -12.4105036  -10.9286583
# 11       AY.4 - B.1.617.2  -11.8318367  -12.2025906  -11.4623078
# 12     AY.4.2 - B.1.617.2    1.0685013    0.1462315    1.9992646


# or with generation time of 4.7 days (Nishiura et al. 2020)
transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*4.7)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
#                  contrast      estimate    asymp.LCL    asymp.UCL
# 1       other - B.1.617.2    0.08992808   -0.5924920    0.7734951
# 2  (B.1.177+) - B.1.617.2 -149.47793484 -160.8433731 -138.6077102
# 3     B.1.1.7 - B.1.617.2  -89.67085793  -92.4569303  -86.9251177
# 4        AY.3 - B.1.617.2  -14.75811332  -21.3318319   -8.5405566
# 5        AY.7 - B.1.617.2  -20.03452835  -21.7440707  -18.3489915
# 6       AY.36 - B.1.617.2   -1.11586053   -5.8420896    3.5189459
# 7       AY.34 - B.1.617.2    0.78689545   -2.6207630    4.2421535
# 8        AY.9 - B.1.617.2  -18.89199775  -20.0545502  -17.7407029
# 9        AY.6 - B.1.617.2  -15.53557694  -16.3863870  -14.6909865
# 10       AY.5 - B.1.617.2   -9.88903770  -10.5138732   -9.2677350
# 11       AY.4 - B.1.617.2  -10.02753693  -10.3391762   -9.7167778
# 12     AY.4.2 - B.1.617.2    0.91237627    0.1249481    1.7059971

emtrsanger_byregion = emtrends(fit4_sanger_multi, trt.vs.ctrl1 ~ LINEAGE|REGION,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_byregion = data.frame(confint(emtrsanger_byregion, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger_byregion$contrasts)$p.value)
delta_r_sanger_byregion
# contrast                   REGION      estimate           SE df     asymp.LCL     asymp.UCL      p.value
# 1        other - B.1.617.2                   London -0.0079562470 0.0013334977 NA -1.056985e-02 -0.0053426395 7.601798e-08
# 2   (B.1.177+) - B.1.617.2                   London -0.2810692187 0.0162661534 NA -3.129503e-01 -0.2491881438 7.169820e-13
# 3      B.1.1.7 - B.1.617.2                   London -0.1416040969 0.0027971714 NA -1.470865e-01 -0.1361217417 7.169820e-13
# 4         AY.3 - B.1.617.2                   London -0.0201752719 0.0090597114 NA -3.793198e-02 -0.0024185638 2.027544e-01
# 5         AY.7 - B.1.617.2                   London -0.0613832418 0.0021065919 NA -6.551209e-02 -0.0572543975 7.169820e-13
# 6        AY.36 - B.1.617.2                   London -0.0135626485 0.0058242170 NA -2.497790e-02 -0.0021473929 1.634446e-01
# 7        AY.34 - B.1.617.2                   London -0.0305173265 0.0041601477 NA -3.867107e-02 -0.0223635869 2.216149e-11
# 8         AY.9 - B.1.617.2                   London -0.0540887092 0.0017924024 NA -5.760175e-02 -0.0505756650 7.169820e-13
# 9         AY.6 - B.1.617.2                   London -0.0281531932 0.0012069926 NA -3.051886e-02 -0.0257875313 7.169820e-13
# 10        AY.5 - B.1.617.2                   London -0.0310999508 0.0011172376 NA -3.328970e-02 -0.0289102053 7.169820e-13
# 11        AY.4 - B.1.617.2                   London -0.0265147492 0.0006876756 NA -2.786257e-02 -0.0251669298 7.169820e-13
# 12      AY.4.2 - B.1.617.2                   London -0.0146976164 0.0012884037 NA -1.722284e-02 -0.0121723915 7.288614e-13
# 13       other - B.1.617.2               North West  0.0283981526 0.0017203319 NA  2.502636e-02  0.0317699411 7.169820e-13
# 14  (B.1.177+) - B.1.617.2               North West -0.1533670766 0.0081663716 NA -1.693729e-01 -0.1373612824 7.169820e-13
# 15     B.1.1.7 - B.1.617.2               North West -0.0773594411 0.0039232598 NA -8.504889e-02 -0.0696699931 7.169820e-13
# 16        AY.3 - B.1.617.2               North West -0.0143637772 0.0117775812 NA -3.744741e-02  0.0087198577 7.952507e-01
# 17        AY.7 - B.1.617.2               North West -0.0411298136 0.0030331894 NA -4.707476e-02 -0.0351848715 7.169820e-13
# 18       AY.36 - B.1.617.2               North West -0.0100317509 0.0068846112 NA -2.352534e-02  0.0034618390 6.545798e-01
# 19       AY.34 - B.1.617.2               North West  0.0025143544 0.0082097352 NA -1.357643e-02  0.0186051396 9.991317e-01
# 20        AY.9 - B.1.617.2               North West -0.0490831474 0.0033495845 NA -5.564821e-02 -0.0425180824 7.169820e-13
# 21        AY.6 - B.1.617.2               North West -0.0418021156 0.0021155974 NA -4.594861e-02 -0.0376556209 7.169820e-13
# 22        AY.5 - B.1.617.2               North West -0.0387221904 0.0020004256 NA -4.264295e-02 -0.0348014283 7.169820e-13
# 23        AY.4 - B.1.617.2               North West -0.0393400430 0.0006299745 NA -4.057477e-02 -0.0381053157 7.169820e-13
# 24      AY.4.2 - B.1.617.2               North West -0.0112755008 0.0019282432 NA -1.505479e-02 -0.0074962135 1.457023e-07
# 25       other - B.1.617.2               South West -0.0295667053 0.0021979835 NA -3.387467e-02 -0.0252587367 7.169820e-13
# 26  (B.1.177+) - B.1.617.2               South West -0.2328378143 0.0220798755 NA -2.761136e-01 -0.1895620535 7.586154e-13
# 27     B.1.1.7 - B.1.617.2               South West -0.1705004137 0.0049509144 NA -1.802040e-01 -0.1607967998 7.169820e-13
# 28        AY.3 - B.1.617.2               South West -0.0338195331 0.0074840958 NA -4.848809e-02 -0.0191509749 1.019483e-04
# 29        AY.7 - B.1.617.2               South West -0.0305476613 0.0031759019 NA -3.677231e-02 -0.0243230079 7.650547e-13
# 30       AY.36 - B.1.617.2               South West -0.0053243847 0.0234215922 NA -5.122986e-02  0.0405810925 9.997627e-01
# 31       AY.34 - B.1.617.2               South West  0.0160702841 0.0082404841 NA -8.076802e-05  0.0322213363 3.399593e-01
# 32        AY.9 - B.1.617.2               South West -0.0636745546 0.0033960185 NA -7.033063e-02 -0.0570184807 7.169820e-13
# 33        AY.6 - B.1.617.2               South West -0.0245140255 0.0019494625 NA -2.833490e-02 -0.0206931493 7.170931e-13
# 34        AY.5 - B.1.617.2               South West -0.0200347289 0.0016219546 NA -2.321370e-02 -0.0168557563 7.172041e-13
# 35        AY.4 - B.1.617.2               South West -0.0115340633 0.0008881847 NA -1.327487e-02 -0.0097932534 7.169820e-13
# 36      AY.4.2 - B.1.617.2               South West  0.0305717649 0.0014633916 NA  2.770357e-02  0.0334399598 7.169820e-13
# 37       other - B.1.617.2               South East -0.0016642287 0.0023027743 NA -6.177583e-03  0.0028491260 9.686688e-01
# 38  (B.1.177+) - B.1.617.2               South East -0.1920228066 0.0169100279 NA -2.251659e-01 -0.1588797609 7.310819e-13
# 39     B.1.1.7 - B.1.617.2               South East -0.1253026398 0.0041746525 NA -1.334848e-01 -0.1171204712 7.169820e-13
# 40        AY.3 - B.1.617.2               South East -0.0247985632 0.0090537580 NA -4.254360e-02 -0.0070535237 6.016675e-02
# 41        AY.7 - B.1.617.2               South East -0.0373510073 0.0025014565 NA -4.225377e-02 -0.0324482427 7.169820e-13
# 42       AY.36 - B.1.617.2               South East  0.0168416078 0.0063619064 NA  4.372500e-03  0.0293107153 7.655929e-02
# 43       AY.34 - B.1.617.2               South East -0.0263884299 0.0053349075 NA -3.684466e-02 -0.0159322032 1.439484e-05
# 44        AY.9 - B.1.617.2               South East -0.0454166384 0.0018428138 NA -4.902849e-02 -0.0418047898 7.169820e-13
# 45        AY.6 - B.1.617.2               South East -0.0357498361 0.0011732276 NA -3.804932e-02 -0.0334503522 7.169820e-13
# 46        AY.5 - B.1.617.2               South East -0.0265023508 0.0011949492 NA -2.884441e-02 -0.0241602935 7.169820e-13
# 47        AY.4 - B.1.617.2               South East -0.0208630307 0.0007227855 NA -2.227966e-02 -0.0194463971 7.169820e-13
# 48      AY.4.2 - B.1.617.2               South East  0.0092744304 0.0012349849 NA  6.853904e-03  0.0116949563 7.714940e-12
# 49       other - B.1.617.2          East of England  0.0031886615 0.0016970351 NA -1.374661e-04  0.0065147891 3.817560e-01
# 50  (B.1.177+) - B.1.617.2          East of England -0.1931807070 0.0169856741 NA -2.264720e-01 -0.1598893974 7.309708e-13
# 51     B.1.1.7 - B.1.617.2          East of England -0.1313533277 0.0040287358 NA -1.392495e-01 -0.1234571506 7.169820e-13
# 52        AY.3 - B.1.617.2          East of England -0.0415279160 0.0224320833 NA -8.549399e-02  0.0024381593 3.985945e-01
# 53        AY.7 - B.1.617.2          East of England -0.0648558795 0.0033677841 NA -7.145662e-02 -0.0582551440 7.169820e-13
# 54       AY.36 - B.1.617.2          East of England -0.0057995766 0.0059026645 NA -1.736859e-02  0.0057694332 9.009509e-01
# 55       AY.34 - B.1.617.2          East of England  0.0164666627 0.0068822004 NA  2.977798e-03  0.0299555276 1.417666e-01
# 56        AY.9 - B.1.617.2          East of England -0.0287316486 0.0017026901 NA -3.206886e-02 -0.0253944374 7.169820e-13
# 57        AY.6 - B.1.617.2          East of England -0.0469873347 0.0017511326 NA -5.041949e-02 -0.0435551780 7.169820e-13
# 58        AY.5 - B.1.617.2          East of England -0.0241567747 0.0011257134 NA -2.636313e-02 -0.0219504170 7.169820e-13
# 59        AY.4 - B.1.617.2          East of England -0.0296087977 0.0007863647 NA -3.115004e-02 -0.0280675513 7.169820e-13
# 60      AY.4.2 - B.1.617.2          East of England -0.0063689573 0.0015549207 NA -9.416546e-03 -0.0033213687 6.078001e-04
# 61       other - B.1.617.2            East Midlands -0.0321045317 0.0014547154 NA -3.495572e-02 -0.0292533419 7.169820e-13
# 62  (B.1.177+) - B.1.617.2            East Midlands -0.2463194461 0.0109509428 NA -2.677829e-01 -0.2248559926 7.169820e-13
# 63     B.1.1.7 - B.1.617.2            East Midlands -0.1820829571 0.0031527452 NA -1.882622e-01 -0.1759036900 7.169820e-13
# 64        AY.3 - B.1.617.2            East Midlands -0.0978488517 0.0272925186 NA -1.513412e-01 -0.0443564982 4.243533e-03
# 65        AY.7 - B.1.617.2            East Midlands -0.0503619797 0.0035523941 NA -5.732454e-02 -0.0433994152 7.169820e-13
# 66       AY.36 - B.1.617.2            East Midlands  0.0449776812 0.0106873404 NA  2.403088e-02  0.0659244835 3.837510e-04
# 67       AY.34 - B.1.617.2            East Midlands  0.0384437996 0.0104496170 NA  1.796293e-02  0.0589246727 3.024973e-03
# 68        AY.9 - B.1.617.2            East Midlands -0.0292133004 0.0040834986 NA -3.721681e-02 -0.0212097902 6.858147e-11
# 69        AY.6 - B.1.617.2            East Midlands -0.0493320895 0.0028315674 NA -5.488186e-02 -0.0437823193 7.169820e-13
# 70        AY.5 - B.1.617.2            East Midlands -0.0221915789 0.0017786886 NA -2.567774e-02 -0.0187054133 7.172041e-13
# 71        AY.4 - B.1.617.2            East Midlands -0.0232872950 0.0009665654 NA -2.518173e-02 -0.0213928617 7.169820e-13
# 72      AY.4.2 - B.1.617.2            East Midlands  0.0046204318 0.0023575981 NA -3.756412e-07  0.0092412392 3.344804e-01
# 73       other - B.1.617.2            West Midlands  0.0003263261 0.0022342522 NA -4.052728e-03  0.0047053799 9.999660e-01
# 74  (B.1.177+) - B.1.617.2            West Midlands -0.2223721247 0.0127396260 NA -2.473413e-01 -0.1974029166 7.169820e-13
# 75     B.1.1.7 - B.1.617.2            West Midlands -0.1520270266 0.0045436509 NA -1.609324e-01 -0.1431216344 7.169820e-13
# 76        AY.3 - B.1.617.2            West Midlands -0.0541810351 0.0204622108 NA -9.428623e-02 -0.0140758388 7.644116e-02
# 77        AY.7 - B.1.617.2            West Midlands -0.0299665172 0.0052256308 NA -4.020857e-02 -0.0197244691 2.679163e-07
# 78       AY.36 - B.1.617.2            West Midlands  0.0430039873 0.0132327493 NA  1.706828e-02  0.0689396993 1.326923e-02
# 79       AY.34 - B.1.617.2            West Midlands  0.0616155313 0.0075875178 NA  4.674427e-02  0.0764867928 8.885115e-13
# 80        AY.9 - B.1.617.2            West Midlands -0.0131572892 0.0029219412 NA -1.888419e-02 -0.0074303897 1.093422e-04
# 81        AY.6 - B.1.617.2            West Midlands -0.0163947731 0.0024356981 NA -2.116865e-02 -0.0116208925 9.196803e-10
# 82        AY.5 - B.1.617.2            West Midlands -0.0127778164 0.0018639955 NA -1.643118e-02 -0.0091244524 4.336832e-10
# 83        AY.4 - B.1.617.2            West Midlands -0.0078180833 0.0007504209 NA -9.288881e-03 -0.0063472853 7.619461e-13
# 84      AY.4.2 - B.1.617.2            West Midlands  0.0034972522 0.0014150007 NA  7.239019e-04  0.0062706026 1.181138e-01
# 85       other - B.1.617.2               North East  0.0373906529 0.0032194331 NA  3.108068e-02  0.0437006258 7.237544e-13
# 86  (B.1.177+) - B.1.617.2               North East -0.0646697290 0.0127015307 NA -8.956427e-02 -0.0397751863 7.163430e-06
# 87     B.1.1.7 - B.1.617.2               North East -0.1050547674 0.0076105525 NA -1.199712e-01 -0.0901383585 7.169820e-13
# 88        AY.3 - B.1.617.2               North East -0.0253139182 0.0277910279 NA -7.978333e-02  0.0291554956 9.245983e-01
# 89        AY.7 - B.1.617.2               North East -0.0215367189 0.0048917837 NA -3.112444e-02 -0.0119489991 1.691023e-04
# 90       AY.36 - B.1.617.2               North East -0.0596249217 0.0299634157 NA -1.183521e-01 -0.0008977062 3.176636e-01
# 91       AY.34 - B.1.617.2               North East -0.0342159608 0.0256526953 NA -8.449432e-02  0.0160623980 7.309680e-01
# 92        AY.9 - B.1.617.2               North East -0.0362411647 0.0047731575 NA -4.559638e-02 -0.0268859478 4.803491e-12
# 93        AY.6 - B.1.617.2               North East -0.0228177343 0.0039395371 NA -3.053909e-02 -0.0150963835 1.967971e-07
# 94        AY.5 - B.1.617.2               North East  0.0084543236 0.0028429740 NA  2.882197e-03  0.0140264503 3.110964e-02
# 95        AY.4 - B.1.617.2               North East -0.0117391301 0.0015262859 NA -1.473060e-02 -0.0087476646 2.873146e-12
# 96      AY.4.2 - B.1.617.2               North East  0.0057247809 0.0044876392 NA -3.070830e-03  0.0145203921 7.646085e-01
# 97       other - B.1.617.2 Yorkshire and The Humber  0.0037091726 0.0029949349 NA -2.160792e-03  0.0095791371 7.851415e-01
# 98  (B.1.177+) - B.1.617.2 Yorkshire and The Humber -0.1647573728 0.0080852539 NA -1.806042e-01 -0.1489105663 7.169820e-13
# 99     B.1.1.7 - B.1.617.2 Yorkshire and The Humber -0.1404771416 0.0056987956 NA -1.516466e-01 -0.1293077074 7.169820e-13
# 100       AY.3 - B.1.617.2 Yorkshire and The Humber  0.0484315723 0.0128210827 NA  2.330271e-02  0.0735604326 2.100681e-03
# 101       AY.7 - B.1.617.2 Yorkshire and The Humber -0.0125444711 0.0037333812 NA -1.986176e-02 -0.0052271783 9.234974e-03
# 102      AY.36 - B.1.617.2 Yorkshire and The Humber -0.0317291995 0.0099861749 NA -5.130174e-02 -0.0121566565 1.672627e-02
# 103      AY.34 - B.1.617.2 Yorkshire and The Humber -0.0289796806 0.0055828559 NA -3.992188e-02 -0.0180374842 4.400699e-06
# 104       AY.9 - B.1.617.2 Yorkshire and The Humber -0.0117569130 0.0025562291 NA -1.676703e-02 -0.0067467960 7.134077e-05
# 105       AY.6 - B.1.617.2 Yorkshire and The Humber -0.0107754720 0.0017355694 NA -1.417713e-02 -0.0073738185 1.959586e-08
# 106       AY.5 - B.1.617.2 Yorkshire and The Humber -0.0135451671 0.0017135122 NA -1.690359e-02 -0.0101867450 1.275979e-12
# 107       AY.4 - B.1.617.2 Yorkshire and The Humber -0.0122829689 0.0008807864 NA -1.400928e-02 -0.0105566592 7.169820e-13
# 108     AY.4.2 - B.1.617.2 Yorkshire and The Humber -0.0039547701 0.0017794471 NA -7.442422e-03 -0.0004671179 2.046124e-01




# results using simpler model fit2_sanger_multi
emtrsanger2 = emtrends(fit2_sanger_multi, trt.vs.ctrl1 ~ LINEAGE,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger2 = data.frame(confint(emtrsanger2, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger2$contrasts)$p.value)
delta_r_sanger2
#                  contrast     estimate           SE df    asymp.LCL    asymp.UCL p.value
# 1       other - B.1.617.2 -0.117463378 8.212644e-07 NA -0.117464988 -0.117461768       0
# 2  (B.1.177+) - B.1.617.2 -0.126114580 8.324083e-07 NA -0.126116212 -0.126112949       0
# 3     B.1.1.7 - B.1.617.2 -0.084749905 6.697615e-07 NA -0.084751218 -0.084748592       0
# 4        AY.3 - B.1.617.2  0.027924500 6.227197e-06 NA  0.027912294  0.027936705       0
# 5        AY.7 - B.1.617.2 -0.011485677 1.508384e-06 NA -0.011488634 -0.011482721       0
# 6       AY.36 - B.1.617.2  0.041237519 4.985606e-06 NA  0.041227747  0.041247290       0
# 7       AY.34 - B.1.617.2  0.040621623 3.709901e-06 NA  0.040614352  0.040628895       0
# 8        AY.9 - B.1.617.2 -0.012533279 1.138186e-06 NA -0.012535510 -0.012531048       0
# 9        AY.6 - B.1.617.2  0.003686990 7.974531e-07 NA  0.003685427  0.003688553       0
# 10       AY.5 - B.1.617.2  0.002183768 6.527144e-07 NA  0.002182489  0.002185047       0
# 11       AY.4 - B.1.617.2  0.002230561 2.963084e-07 NA  0.002229980  0.002231142       0
# 12     AY.4.2 - B.1.617.2  0.036269047 6.519412e-07 NA  0.036267769  0.036270325       0

# results using simpler model fit1_sanger_multi
emtrsanger1 = emtrends(fit1_sanger_multi, trt.vs.ctrl1 ~ LINEAGE,  
                       var="DATE_NUM",  mode="latent",
                       at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger1 = data.frame(confint(emtrsanger1, 
                                     adjust="none", df=NA)$contrasts, 
                             p.value=as.data.frame(emtrsanger2$contrasts)$p.value)
delta_r_sanger1
#                  contrast      estimate           SE df     asymp.LCL     asymp.UCL p.value
# 1       other - B.1.617.2 -0.1160508461 6.969309e-07 NA -0.1160522121 -0.1160494802       0
# 2  (B.1.177+) - B.1.617.2 -0.1243448479 6.947756e-07 NA -0.1243462096 -0.1243434861       0
# 3     B.1.1.7 - B.1.617.2 -0.0844836053 5.526911e-07 NA -0.0844846886 -0.0844825221       0
# 4        AY.3 - B.1.617.2  0.0291081670 3.799426e-06 NA  0.0291007203  0.0291156138       0
# 5        AY.7 - B.1.617.2 -0.0160410526 9.904481e-07 NA -0.0160429939 -0.0160391114       0
# 6       AY.36 - B.1.617.2  0.0445657644 2.438578e-06 NA  0.0445609848  0.0445705439       0
# 7       AY.34 - B.1.617.2  0.0366963869 2.103749e-06 NA  0.0366922637  0.0367005102       0
# 8        AY.9 - B.1.617.2 -0.0195890602 6.468192e-07 NA -0.0195903279 -0.0195877924       0
# 9        AY.6 - B.1.617.2 -0.0005600488 5.262820e-07 NA -0.0005610803 -0.0005590173       0
# 10       AY.5 - B.1.617.2 -0.0030056426 4.697117e-07 NA -0.0030065632 -0.0030047220       0
# 11       AY.4 - B.1.617.2  0.0023933945 2.670646e-07 NA  0.0023928710  0.0023939179       0
# 12     AY.4.2 - B.1.617.2  0.0336399137 5.406211e-07 NA  0.0336388541  0.0336409733       0


# predicted incidences on average over all ONS regions from multinomial fit
# fitted prop of different LINEAGES in England today
multinom_preds_today_avg = data.frame(emmeans(fit4_sanger_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#      LINEAGE         prob           SE df     asymp.LCL    asymp.UCL
# 1  B.1.617.2 2.078192e-01 2.244213e-03 NA  2.034206e-01 2.122178e-01
# 2      other 1.084557e-02 6.960358e-04 NA  9.481365e-03 1.220978e-02
# 3   B.1.177+ 3.272671e-10 7.340559e-10 NA -1.111456e-09 1.765990e-09
# 4    B.1.1.7 2.000803e-07 8.775228e-08 NA  2.808903e-08 3.720716e-07
# 5       AY.3 1.288851e-03 2.819744e-04 NA  7.361911e-04 1.841510e-03
# 6       AY.7 1.473056e-03 1.466811e-04 NA  1.185566e-03 1.760546e-03
# 7      AY.36 5.613567e-03 6.841538e-04 NA  4.272650e-03 6.954483e-03
# 8      AY.34 7.950355e-03 8.672091e-04 NA  6.250657e-03 9.650054e-03
# 9       AY.9 3.874581e-03 2.118446e-04 NA  3.459373e-03 4.289789e-03
# 10      AY.6 1.642100e-02 5.161566e-04 NA  1.540935e-02 1.743265e-02
# 11      AY.5 3.144465e-02 8.984985e-04 NA  2.968363e-02 3.320567e-02
# 12      AY.4 6.045102e-01 2.943547e-03 NA  5.987409e-01 6.102794e-01
# 13    AY.4.2 1.087588e-01 2.351650e-03 NA  1.041497e-01 1.133680e-01

# % non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.7921808 0.7733694 0.8109922 

# here given by region:
multinom_preds_today_byregion = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE|DATE_NUM, by=c("REGION"),
                                                   at=list(DATE_NUM=today_num), 
                                                   mode="prob", df=NA))
multinom_preds_today_byregion

multinom_preds_delta_today_byregion = multinom_preds_today_byregion[multinom_preds_today_byregion$LINEAGE=="AY.4.2",]
multinom_preds_delta_today_byregion
#     LINEAGE DATE_NUM                   REGION       prob          SE df  asymp.LCL  asymp.UCL
# 13   AY.4.2    18918                   London 0.15118941 0.004255598 NA 0.14284859 0.15953023
# 26   AY.4.2    18918               North West 0.09420867 0.002796490 NA 0.08872765 0.09968969
# 39   AY.4.2    18918               South West 0.13807667 0.004206891 NA 0.12983131 0.14632202
# 52   AY.4.2    18918               South East 0.15441578 0.003792914 NA 0.14698181 0.16184976
# 65   AY.4.2    18918          East of England 0.09752246 0.003197061 NA 0.09125634 0.10378859
# 78   AY.4.2    18918            East Midlands 0.09456926 0.003606529 NA 0.08750059 0.10163793
# 91   AY.4.2    18918            West Midlands 0.14088294 0.003626969 NA 0.13377422 0.14799167
# 104  AY.4.2    18918               North East 0.06838891 0.004302484 NA 0.05995619 0.07682162
# 117  AY.4.2    18918 Yorkshire and The Humber 0.09552138 0.003110142 NA 0.08942562 0.10161715


# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME ####

gentime = 4.7 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantages of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit4_sanger_multi = nnet::multinom(LINEAGE ~ REGION * ns(DATE_NUM, df=2), data=sanger, maxit=1000) here
emtrsanger4 = emtrends(fit4_sanger_multi, trt.vs.ctrl1 ~ LINEAGE, by=c("DATE_NUM","REGION"), 
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
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
date.to = today_num # max(sanger$DATE_NUM)+extrapolate

# predictions by ONS region (without CIs) ####
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to), 
                            REGION=levels_REGION))
fit_multi_preds = data.frame(predgrid, as.data.frame(predict(fit6_sanger_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fit_multi_preds = gather(fit_multi_preds, LINEAGE, prob, all_of(levels_LINEAGE), factor_key=TRUE)
fit_multi_preds$collection_date = as.Date(fit_multi_preds$DATE_NUM, origin="1970-01-01")
fit_multi_preds$LINEAGE = factor(fit_multi_preds$LINEAGE, levels=levels_LINEAGE) 

# library(glm.predict)
# fit_multi_preds_withCIs = basepredict.multinom(model=fit6_sanger_multi, values=predgrid, type="simulation",
#                                            sim.count=10, conf.int=0.95, sigma=NULL, set.seed=NULL)

muller_sangerbyregion_mfit = ggplot(data=fit_multi_preds, 
                                      aes(x=collection_date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbyregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit.png"), width=8, height=6)


library(ggpubr)
ggarrange(muller_sangerbyregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbyregion_mfit+ggtitle("Multinomial fit (LINEAGE~REGION*ns(DATE,df=3))"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.png"), width=8, height=10)





# predictions by ONS region (with CIs) ####

fit_sanger_multi_predsbyregion = data.frame(emmeans(fit4_sanger_multi, 
                                                    ~ LINEAGE,
                                                    by=c("DATE_NUM", "REGION"),
                                                    at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                                    mode="prob", df=NA))
fit_sanger_multi_predsbyregion$collection_date = as.Date(fit_sanger_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion$LINEAGE = factor(fit_sanger_multi_predsbyregion$LINEAGE, levels=levels_LINEAGE_plot) 
fit_sanger_multi_predsbyregion$REGION = factor(fit_sanger_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of England today

muller_sangerbyregion_mfit = ggplot(data=fit_sanger_multi_predsbyregion, 
                                    aes(x=collection_date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
                                                  ~ LINEAGE,
                                                  by=c("DATE_NUM", "NHSREGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                                  mode="prob", df=NA))
fit_sanger_multi_predsbynhsregion$collection_date = as.Date(fit_sanger_multi_predsbynhsregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbynhsregion$LINEAGE = factor(fit_sanger_multi_predsbynhsregion$LINEAGE, levels=levels_LINEAGE_plot) 
fit_sanger_multi_predsbynhsregion$NHSREGION = factor(fit_sanger_multi_predsbynhsregion$NHSREGION, levels=levels_NHSREGION) 
# predicted incidence in different parts of England today
# fit_sanger_multi_predsbynhsregion[fit_sanger_multi_predsbynhsregion$collection_date==today,]

fit_sanger_multi_preds = data.frame(emmeans(fit4_sanger_nhs_multi, 
                                           ~ LINEAGE,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                           mode="prob", df=NA))
fit_sanger_multi_preds$collection_date = as.Date(fit_sanger_multi_preds$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_preds$LINEAGE = factor(fit_sanger_multi_preds$LINEAGE, levels=levels_LINEAGE_plot) 

muller_sanger_mfit = ggplot(data=fit_sanger_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)")
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
                                  aes(x=collection_date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ NHSREGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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

# # PHE SGTF data to superimpose on model predictions as independent estimate of prop of B.1.1.7 (these data are not up to date anymore)
# sgtf_UK = read.csv(".\\data\\SGTF_UK\\sgtf_2021-05-31_by adm region.csv") # sgtf_2021-05-11_by NHS region, TO DO : ask Nick for update
# colnames(sgtf_UK) = c("collection_date", "ereg_code", "sgtf", "non_sgtf", "REGION")
# sgtf_UK$ereg_code = NULL
# sgtf_UK$collection_date = as.Date(sgtf_UK$collection_date)
# sgtf_UK$DATE_NUM = as.numeric(sgtf_UK$collection_date)
# sgtf_UK$REGION = factor(sgtf_UK$REGION, levels=levels_REGION)
# sgtf_UK = sgtf_UK[sgtf_UK$collection_date>="2021-01-01",]
# sgtf_UK$total = sgtf_UK$sgtf + sgtf_UK$non_sgtf
# sgtf_UK = sgtf_UK[sgtf_UK$total!=0,]
# sgtf_UK$Week = lubridate::week(sgtf_UK$collection_date)
# sgtf_UK$LINEAGE = "B.1.1.7 (S dropout)"
# sgtf_UK$LINEAGE = factor(sgtf_UK$LINEAGE, levels=c(levels_LINEAGE_plot, "B.1.1.7 (S dropout)"))
# sgtf_UK2 = sgtf_UK
# sgtf_UK2$count = sgtf_UK2$sgtf 
# sgtf_UK2$collection_date_num = sgtf_UK2$DATE_NUM
# sgtf_UK2$DATE_NUM = NULL
# sgtf_UK2$sgtf = NULL
# sgtf_UK2$non_sgtf = NULL
# sgtf_UK2$prop = sgtf_UK2$count / sgtf_UK2$total
# head(sgtf_UK2)  
# range(sgtf_UK$collection_date) # "2021-01-01" "2021-05-30"
# library(tidyr)
# sgtf_UK_long = gather(sgtf_UK, Sdropout, count, sgtf:non_sgtf, factor_key=TRUE)
# sgtf_UK_long$collection_date
# 
# 
# # plot of SGTF data by region
# ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-01-01"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
#   facet_wrap(~REGION, scales="free_y") +
#   geom_col(position="fill") +
#   scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                  lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                    lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   ylab("Share") +
#   ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN ENGLAND\n(data PHE)") +
#   scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
#                      labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
#                      limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
#   xlab("Collection date")
#   
# ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by ONS region_stacked as proportion.png"), width=10, height=8)
# # ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by ONS region_stacked as proportion.pdf"), width=10, height=8)
# 
# ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-01-01"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
#   # facet_wrap(~REGION, scales="free_y") +
#   geom_col(position="fill") +
#   scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                  lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                    lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   ylab("Share") +
#   ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN ENGLAND\n(data PHE)") +
#   scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
#                      labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
#                      limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
#   xlab("Collection date")
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.png"), width=10, height=8)
# # ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.pdf"), width=10, height=8)
# 
# 
# 
# # logistic spline fit to SGTF data
# 
# fit_SGTF = glm(cbind(sgtf,non_sgtf) ~ REGION * ns(DATE_NUM, df=3), family=binomial, data=sgtf_UK) # at level of ONS regions
# BIC(fit_SGTF)
# 
# SGTF_predsbyregion = data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
#                                                       at=list(DATE_NUM=seq(date.from, date.to)), 
#                                                       type="response"))
# SGTF_predsbyregion$collection_date = as.Date(SGTF_predsbyregion$DATE_NUM, origin="1970-01-01")
# SGTF_predsbyregion$LINEAGE = "B.1.1.7 (S dropout)"
# SGTF_predsbyregion$LINEAGE = factor(SGTF_predsbyregion$LINEAGE, levels=c(levels_LINEAGE_plot, "B.1.1.7 (S dropout)"))
# 
# # fitted S dropout in different parts of England today
# 1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM,
#         at=list(DATE_NUM=today_num), 
#         type="response"))[,c(2,6,5)]
# #        prob asymp.UCL asymp.LCL
# # 1 0.9997718 0.9996945 0.9998295
# 
# # here given by NHS region:
# data.frame(REGION=data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
#                    at=list(DATE_NUM=today_num), 
#                    type="response"))[,2], 
#            1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
#                         at=list(DATE_NUM=today_num), 
#                         type="response"))[,c(3,7,6)])
# # REGION      prob asymp.UCL asymp.LCL
# # 1                   London 0.9924462 0.9883269 0.9951190
# # 2               North West 0.9999878 0.9999832 0.9999911
# # 3               South West 0.9998023 0.9984491 0.9999748
# # 4               South East 0.9998247 0.9996719 0.9999063
# # 5          East of England 0.9996771 0.9994547 0.9998088
# # 6            East Midlands 0.9995856 0.9993293 0.9997439
# # 7            West Midlands 0.9997396 0.9995751 0.9998404
# # 8               North East 0.9998467 0.9996065 0.9999403
# # 9 Yorkshire and The Humber 0.9999026 0.9998289 0.9999445
# 
# 

# plot predictions multinomial fit by ONS region
# on logit scale:

fit_sanger_multi_preds2 = fit_multi_preds # fit_sanger_multi_preds
ymin = 0.001
ymax = 0.9998
fit_sanger_multi_preds2$asymp.LCL[fit_sanger_multi_preds2$asymp.LCL<ymin] = ymin
fit_sanger_multi_preds2$asymp.UCL[fit_sanger_multi_preds2$asymp.UCL<ymin] = ymin
fit_sanger_multi_preds2$asymp.UCL[fit_sanger_multi_preds2$asymp.UCL>ymax] = ymax
fit_sanger_multi_preds2$prob[fit_sanger_multi_preds2$prob<ymin] = ymin

fit_sanger_multi_predsbyregion2 = fit_multi_preds # fit_sanger_multi_predsbyregion
fit_sanger_multi_predsbyregion2$asymp.LCL[fit_sanger_multi_predsbyregion2$asymp.LCL<ymin] = ymin
fit_sanger_multi_predsbyregion2$asymp.UCL[fit_sanger_multi_predsbyregion2$asymp.UCL<ymin] = ymin
fit_sanger_multi_predsbyregion2$asymp.UCL[fit_sanger_multi_predsbyregion2$asymp.UCL>ymax] = ymax
fit_sanger_multi_predsbyregion2$prob[fit_sanger_multi_predsbyregion2$prob<ymin] = ymin

# plots of multinomial fit to Sanger Inst data by ONS region

plot_sanger_mfit_logit = qplot(data=fit_sanger_multi_predsbyregion2, 
                               x=collection_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=LINEAGE
  # ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data,\nmultinomial fit LINEAGE~REGION*ns(DATE,df=3))") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(0.4), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.001, 3), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), ylim=c(0.001, 0.9991), expand=c(0,0))
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
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-09-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE, 
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
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(1), pch=I(16)) +
  geom_point(data=sgtf_UK2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(0.5), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), ylim=c(0.001, 0.9991), expand=c(0,0))
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
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-09-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(1), pch=I(16)) +
  geom_point(data=sgtf_UK2,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE, 
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

# daily new cases by ONS region :
cases_uk_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDate&format=csv")
cases_uk_region$date = as.Date(cases_uk_region$date)
cases_uk_region$DATE_NUM = as.numeric(cases_uk_region$date)
cases_uk_region$REGION = factor(cases_uk_region$areaName, levels=levels_REGION)
cases_uk_region$areaName = NULL
cases_uk_region = cases_uk_region[cases_uk_region$date<=(max(cases_uk_region$date)-3),] # cut off data from last 3 days (incomplete)

ggplot(data=cases_uk_region, 
       aes(x=date, y=newCasesBySpecimenDate, group=REGION, colour=REGION, fill=REGION)) +
  facet_wrap(~REGION, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 60, fx=T)) + 
  # scale_y_log10() + 
  ylab("New cases (per day)") + 
  ggtitle("CONFIRMED SARS-CoV2 CASES PER DAY IN ENGLAND BY ONS REGION", "data gov.uk, https://coronavirus.data.gov.uk/details/download\n(note: low case ascertainment during 1st wave)") +
  theme_hc() +
  theme(legend.position="none")

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases England by ONS region.png"), width=8, height=6)

# daily  new cases by ONS region & age group :
cases_uk_age_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateAgeDemographics&format=csv")
cases_uk_age_region$date = as.Date(cases_uk_age_region$date)
cases_uk_age_region$DATE_NUM = as.numeric(cases_uk_age_region$date)
cases_uk_age_region$REGION = factor(cases_uk_age_region$areaName, levels=levels_REGION)
cases_uk_age_region$areaName = NULL
cases_uk_age_region=cases_uk_age_region[cases_uk_age_region$age!="unassigned",]

ggplot(data=cases_uk_age_region[!cases_uk_age_region$age %in% c("00_59","60+"),], 
       aes(x=date, y=cases, group=age, colour=age, fill=age)) +
  facet_wrap(~REGION, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  # geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 60, fx=T)) + 
  scale_y_log10() + 
  ylab("diagnosed cases (per day)") + 
  ggtitle("DIAGNOSED SARS-CoV2 CASES BY AGE IN ENGLAND BY ONS REGION", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="right")

ggsave(file=paste0(".\\plots\\",plotdir,"\\diagnosed cases per day by age by ONS region.png"), width=8, height=6)


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
hosps_uk_nhsregion = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=newAdmissions&metric=hospitalCases&format=csv") # &metric=newCasesBySpecimenDate gives NAs
hosps_uk_nhsregion$date = as.Date(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$DATE_NUM = as.numeric(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$REGION = factor(hosps_uk_nhsregion$areaName, levels=levels_NHSREGION)

ggplot(data=hosps_uk_nhsregion, 
       aes(x=date, y=newAdmissions, group=REGION, colour=REGION, fill=REGION)) +
  facet_wrap(~REGION, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 60, fx=T)) + 
  # scale_y_log10() + 
  ylab("New hospital admissions (per day)") + 
  ggtitle("COVID HOSPITAL ADMISSIONS PER DAY IN ENGLAND BY NHS REGION", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="none")

ggsave(file=paste0(".\\plots\\",plotdir,"\\hospital admissions England by NHS region.png"), width=8, height=6)

ggplot(data=hosps_uk_nhsregion, 
       aes(x=date, y=hospitalCases, group=REGION, colour=REGION, fill=REGION)) +
  facet_wrap(~REGION, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL)) +
  geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 70, fx=T)) + 
  # scale_y_log10() + 
  ylab("Covid patients currently admitted") + 
  ggtitle("COVID PATIENTS CURRENTLY IN HOSPITAL IN ENGLAND BY NHS REGION", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="none")

ggsave(file=paste0(".\\plots\\",plotdir,"\\covid patients by NHS region.png"), width=8, height=6)


# hospitalisations by age :
hosps_eng_age = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=cumAdmissionsByAge&format=csv")
hosps_eng_age$date = as.Date(hosps_eng_age$date)
hosps_eng_age$DATE_NUM = as.numeric(hosps_eng_age$date)

hosps_eng_age = data.frame(hosps_eng_age %>% group_by(age) %>% mutate(Diff = value - lag(value,1,order_by=date))) # see https://stackoverflow.com/questions/26291988/how-to-create-a-lag-variable-within-each-group/26292059
hosps_eng_age$age = factor(hosps_eng_age$age, levels=c("0_to_5", "6_to_17", "18_to_64", "65_to_84", "85+"))
hosps_eng_age$week = cut(hosps_eng_age$date, breaks="week")
hosps_eng_age_byweek = hosps_eng_age %>% group_by(week, age) %>% summarise(hosps=sum(Diff))  
hosps_eng_age_byweek$hosps[is.na(hosps_eng_age_byweek$hosps)] = 0

ggplot(data=hosps_eng_age, 
       aes(x=date, y=Diff, group=age, colour=age, fill=age)) +
  facet_wrap(~ age, scale="free_y") +
  # geom_point(cex=I(0.2)) +
  geom_col(aes(alpha=I(0.4), width=I(1), colour=NULL), position="identity") +
  # geom_line() +
  # geom_area(position="stack") +
  stat_smooth(se=FALSE, geom="line", lwd=I(1.1),
              method = 'gam', formula = y ~ s(x, k = 30, fx=T),
              alpha=1, aes(fill=age, colour=age), position="identity") +
  # geom_smooth(method="gam", se=F, formula = y ~ s(x, k = 40, fx=T)) + 
  # scale_y_log10() + 
  ylab("new hospitalisations (per day)") + 
  ggtitle("DAILY COVID HOSPITALISATIONS BY AGE IN ENGLAND", "data gov.uk, https://coronavirus.data.gov.uk/details/download") +
  theme_hc() +
  theme(legend.position="none") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\hospitalisation in england by age.png"), width=8, height=6)



# daily hospital admission also see 
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/COVID-19-daily-admissions-and-beds-20210611.xlsx

# daily hospital admission by age for England
  # https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/Covid-Publication-10-06-2021-Supplementary-Data.xlsx
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/07/Covid-Publication-08-07-2021-Supplementary-Data.xlsx
# PS Next publication: Thursday 8 July 2021
library(rio)
# dat = rio::import(file = "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/07/Covid-Publication-08-07-2021-Supplementary-Data.xlsx",
#            which = 1)
dat = rio::import(file = "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/08/Covid-Publication-12-08-2021-Supplementary-Data.xlsx",
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
hosps_ENG_age$age = factor(hosps_ENG_age$age, levels=c("0-5","6-17","18-24","25-34", "35-44", "45-54","55-64","65-74","75-84","85+"))
hosps_ENG_age$date = as.Date(as.character(hosps_ENG_age$date))

ggplot(data=hosps_ENG_age, 
       aes(x=date, y=newAdmissions, group=age, colour=age, fill=age)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  geom_point(cex=I(0.8)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 35)) + 
  # facet_wrap(~ type) + 
  scale_y_log10() + ylab("Daily new confirmed hospital admissions") + 
  ggtitle("Daily new SARS-CoV2 hospitalisations in England\n(data NHS)") +
  scale_colour_discrete(guide = guide_legend(reverse = TRUE) ) +
  scale_fill_discrete(guide = guide_legend(reverse = TRUE) )

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
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  geom_point(cex=I(0.8)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 35)) + 
  facet_wrap(~ type) + scale_y_log10() + ylab("Daily new confirmed cases & hospital admissions") + 
  ggtitle("Daily new SARS-CoV2 cases & hospitalisations in England\n(data gov.uk & NHS)")

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.pdf"), width=8, height=6)
# write.csv(cases_hosps_ENG_age, file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.csv"), row.names=F)

cases_hosps_ENG_age2 = cases_hosps_ENG_age
lag = 7 # lag = 7 from cases to hosps gives best BIC in Poisson GLM
cases_hosps_ENG_age2$date[cases_hosps_ENG_age2$type=="cases"] = cases_hosps_ENG_age2$date[cases_hosps_ENG_age2$type=="cases"]-lag
cases_hosps_ENG_age_wide = spread(cases_hosps_ENG_age2, type, count)
cases_hosps_ENG_age_wide$DATE_NUM = as.numeric(cases_hosps_ENG_age_wide$date)
cases_hosps_ENG_age_wide = cases_hosps_ENG_age_wide[complete.cases(cases_hosps_ENG_age_wide),]
cases_hosps_ENG_age_wide$CHR = cases_hosps_ENG_age_wide$hospitalisations / cases_hosps_ENG_age_wide$cases
fit = glm(hospitalisations~cases*age*DATE_NUM, data=cases_hosps_ENG_age_wide, family=poisson)
BIC(fit)
qplot(data=cases_hosps_ENG_age_wide, x=date, y=CHR*100, group=age, colour=age, fill=age, geom=c("blank")) + 
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 10)) + 
  # facet_wrap(~age) +
  # scale_x_log10() + 
  scale_y_log10() +
  ylab("Case (-7 days) to hospitalisation ratio (%)") +
  ggtitle("CASE TO HOSPITALISATION RATIO IN ENGLAND\n(data gov.uk & NHS)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\case to hospitalisation ratio England.png"), width=8, height=6)
cases_hosps_ENG_age_wide[cases_hosps_ENG_age_wide$date==max(cases_hosps_ENG_age_wide$date),]


# MAP MULTINOMIAL FIT ONTO CASE NUMBERS ####

newdat = expand.grid(REGION=levels_REGION,
                     DATE_NUM=seq(date.from, date.to))
fit_sanger_multi_predsbyregion_day = data.frame(newdat,
                                          predict(fit4_sanger_multi, 
                                                  newdata = newdat,
                                                  type = "prob"), check.names=F)  
fit_sanger_multi_predsbyregion_day = gather(fit_sanger_multi_predsbyregion_day, LINEAGE, prob, all_of(levels_LINEAGE_plot))
fit_sanger_multi_predsbyregion_day$collection_date = as.Date(fit_sanger_multi_predsbyregion_day$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion_day$LINEAGE = factor(fit_sanger_multi_predsbyregion_day$LINEAGE, levels=levels_LINEAGE_plot)

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
gamfittotcases = gam(totcases ~ s(DATE_NUM, bs="cs", k=25, by=REGION) + REGION + WEEKDAY, family=poisson(log), 
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

fit_sanger_multi_predsbyregion_day$cases[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2"&
                                           fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA
fit_sanger_multi_predsbyregion_day$cases_smoothed[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2"&
                                           fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA
fit_sanger_multi_predsbyregion_day$cases_rollingmean[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2"&
                                                    fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA

# plot new cases per day by region
ggplot(data=fit_sanger_multi_predsbyregion_day,
       aes(x=collection_date, y=cases_smoothed, 
           group=LINEAGE)) +
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_line(aes(lwd=I(1.3), colour=LINEAGE)) +
#  geom_smooth(aes(lwd=I(1.3), colour=LINEAGE), method = "gam", family = "poisson", formula = y ~ s(x, k=4), se=FALSE) +
# geom_smooth(aes(lwd=I(1), colour=LINEAGE), method="loess", span=0.8, se=FALSE) +
  geom_line(data=data.frame(collection_date=cases_uk_region$date,
                            cases_smoothed=cases_uk_region$totcases_smoothed, 
                              REGION=cases_uk_region$REGION,
                              LINEAGE="total"), aes(lwd=I(1.9), alpha=I(0.5)), 
              colour=alpha(I("black"),0.7)) +
#  geom_smooth(data=data.frame(collection_date=cases_uk_region$date,
#                              cases=cases_uk_region$newCasesBySpecimenDate, 
#                              REGION=cases_uk_region$REGION,
#                              LINEAGE="total"), aes(lwd=I(1.9), alpha=I(0.5)), 
#              method = "gam", family = "poisson", formula = y ~ s(x, k=4), se=FALSE, colour=alpha(I("black"),0.7)) +
  #  geom_smooth(data=data.frame(collection_date=cases_uk_region$date,
#                              cases=cases_uk_region$newCasesBySpecimenDate, 
#                              REGION=cases_uk_region$REGION,
#                              LINEAGE="total"), aes(lwd=I(1.5)), method="loess", span=0.8, se=FALSE, colour=I("black")) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
       aes(x=collection_date, y=cases+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
       aes(x=collection_date, y=cases_smoothed+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
       aes(x=collection_date, y=cases_rollingmean, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
       aes(x=collection_date, y=cases_rollingmean+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
#        aes(x=collection_date, y=cases_smoothed, group=LINEAGE)) + 
#   facet_wrap(~ REGION, scale="free", ncol=3) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
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
k=15
fit_cases = gam(totcases ~ s(DATE_NUM, bs="cs", k=k, m=c(2), by=REGION) + REGION + weekday,
                method = "REML",
                knots = list(DATE_NUM = c(min(fit_sanger_multi_predsbyregion_day$DATE_NUM)-14,
                                          seq(min(fit_sanger_multi_predsbyregion_day$DATE_NUM)+1*diff(range(fit_sanger_multi_predsbyregion_day$DATE_NUM))/(k-2), 
                                              max(fit_sanger_multi_predsbyregion_day$DATE_NUM)-1*diff(range(fit_sanger_multi_predsbyregion_day$DATE_NUM))/(k-2), length.out=k-2),
                                          max(fit_sanger_multi_predsbyregion_day$DATE_NUM)+14)),
                family=poisson(log), data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE==levels_LINEAGE[[1]],],
) # TO DO: add testing
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
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
  wt = as.data.frame(emmeans(fit4_sanger_multi, ~ LINEAGE, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
  # wt = rep(1/length(levels_LINEAGE), length(levels_LINEAGE)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
  cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
  names(cons) = seq_along(cons)
  EMT = emtrends(fit3_sanger_multi,  ~ LINEAGE, by=c("DATE_NUM","REGION"),
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
#                                                 wt = as.data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
#                                                 # wt = rep(1/length(levels_LINEAGE), length(levels_LINEAGE)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
#                                                 cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
#                                                 names(cons) = seq_along(cons)
#                                                 EMT = emtrends(fit1_sanger_multi,  ~ LINEAGE, by=c("DATE_NUM","REGION"),
#                                                                var="DATE_NUM", mode="latent",
#                                                                at=list(DATE_NUM=dat,
#                                                                        REGION=levels(SGTF_predsbyregion$REGION)))
#                                                 out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
#                                                 # sum(out$estimate*wt) # should sum to zero
#                                                 return(out) } ))
above_avg_r_variants = above_avg_r_variants4
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE), 
                                       labels=levels(data_agbyweeknhsregion1$LINEAGE))
above_avg_r_variants$LINEAGE = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
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
                LINEAGE="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$LINEAGE = factor(above_avg_r_variants$LINEAGE, levels=c(levels_LINEAGE_plot,"avg"))
above_avg_r_variants$prob = fit_sanger_multi_predsbyregion_day$prob[match(interaction(above_avg_r_variants$DATE_NUM,
                                                                                      above_avg_r_variants$LINEAGE,
                                                                                      above_avg_r_variants$REGION),
                                                                          interaction(fit_sanger_multi_predsbyregion_day$DATE_NUM,
                                                                                      fit_sanger_multi_predsbyregion_day$LINEAGE,
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
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$LINEAGE %in% c("other"))|(above_avg_r_variants2$collection_date>max(cases_uk_region$date))),], 
      x=collection_date-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE, fill=LINEAGE, alpha=I(0.5),
      group=LINEAGE, linetype=I(0)) +
  facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=LINEAGE, colour=LINEAGE), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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


