# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN UK BASED ON SANGER INSTITUTE BASELINE SURVEILLANCE SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & surge testing/active surveillance)

# last update 24 MAY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-24")
today_num = as.numeric(today)
today # "2021-05-24"
plotdir = "VOCs_SANGER"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))


sanger = read_tsv("https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv")
sanger = as.data.frame(sanger)
LTLAs_regions = read.csv(".//data/ONS/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")
sanger$REGION = LTLAs_regions$RGN20NM[match(sanger$LTLA, LTLAs_regions$LAD20CD)]
sanger = sanger[!is.na(sanger$REGION),]
sanger$REGION = factor(sanger$REGION)
levels(sanger$REGION)
levels_REGION = c("London","North West","South West","South East","East of England","East Midlands","West Midlands",
                  "North East","Yorkshire and The Humber")
sanger$REGION = factor(sanger$REGION, levels=levels_REGION)
head(sanger)

sanger = sanger[grepl("2021-", sanger[,"WeekEndDate"]),]
sanger = sanger[!(sanger$Lineage=="None"|sanger$Lineage=="Lineage data suppressed"),]
nrow(sanger) # 9434
range(sanger$WeekEndDate) # "2021-01-02" "2021-05-15"
sanger$Week = lubridate::week(sanger$WeekEndDate)
sanger$DATE_NUM = as.numeric(sanger$WeekEndDate)-3.5 # using week midpoint
colnames(sanger)

sanger = sanger[rep(seq_len(nrow(sanger)), sanger$Count),] # convert to long format
sanger$Count = NULL
nrow(sanger) # 148683

nrow(sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 12980 (last 2 weeks)
nrow(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 3750 (last 2 weeks)
nrow(sanger[grepl("B.1.1.7",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 8919

length(unique(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),"LTLA"])) # 192 LTLAs

unique(sanger$Lineage)
sel_target_VOC = "B.1.617"
table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
sanger[sanger$Lineage=="B.1.617",]
sum(table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])) # 4090 B.1.617

unique(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
# "B.1.617"   "B.1.617.1" "B.1.617.3" "B.1.617.2"
sum(grepl("B.1.617",sanger$Lineage, fixed=TRUE)) # 4090 B.1.617+
sanger$LINEAGE1 = sanger$Lineage
sanger$LINEAGE2 = sanger$Lineage
sanger[grepl(sel_target_VOC, sanger$LINEAGE1, fixed=T),"LINEAGE1"] = paste0(sel_target_VOC,"+")
sel_target_VOC = paste0(sel_target_VOC, "+")
sanger[grepl("B.1.177", sanger$LINEAGE1, fixed=T),"LINEAGE1"] = "B.1.177+"
sanger[grepl("B.1.177", sanger$LINEAGE2, fixed=T),"LINEAGE2"] = "B.1.177+"

sum("B.1.1.7"==sanger$LINEAGE1) # 137157 B.1.1.7


table_lineage = as.data.frame(table(sanger$LINEAGE1))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger$LINEAGE1))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# Lineage        Prop
# 8        B.1 0.001156824
# 71   B.1.1.7 0.922479369
# 80  B.1.177+ 0.036937646
# 100  B.1.351 0.001748687
# 138 B.1.617+ 0.027508189

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
# 137157     5492      260      131     4090     1498       55 
table(sanger$LINEAGE2)
# B.1.1.7  B.1.177+   B.1.351   B.1.525 B.1.617.1 B.1.617.2     other       P.1 
# 137157      5492       260       131       125      3956      1507        55 

levels_LINEAGE1 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617+","P.1")
levels_LINEAGE2 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617.1","B.1.617.2","P.1")
sanger$LINEAGE1 = factor(sanger$LINEAGE1, levels=levels_LINEAGE1)
sanger$LINEAGE2 = factor(sanger$LINEAGE2, levels=levels_LINEAGE2)

sanger_lastmonth = sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-30),]
table(sanger_lastmonth$Lineage, sanger_lastmonth$REGION)

str(sanger)

# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of UK
data_agbyweek1 = as.data.frame(table(sanger$Week, sanger$LINEAGE2))
colnames(data_agbyweek1) = c("Week", "LINEAGE2", "count")
data_agbyweek1_sum = aggregate(count ~ Week, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$Week, data_agbyweek1_sum$Week)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE2=="B.1.617.1","total"]) == nrow(sanger) # TRUE
data_agbyweek1$Week = as.numeric(as.character(data_agbyweek1$Week))
data_agbyweek1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and region
data_agbyweekregion1 = as.data.frame(table(sanger$Week, sanger$REGION, sanger$LINEAGE2))
colnames(data_agbyweekregion1) = c("Week", "REGION", "LINEAGE2", "count")
data_agbyweekregion1_sum = aggregate(count ~ Week + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$Week,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$Week,data_agbyweekregion1_sum$REGION))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE2=="B.1.617.1","total"]) == nrow(sanger) # TRUE
data_agbyweekregion1$Week = as.numeric(as.character(data_agbyweekregion1$Week))
data_agbyweekregion1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion1[data_agbyweekregion1$collection_date==max(data_agbyweekregion1$collection_date),]

# MULLER PLOT (RAW DATA)
unique(sanger$LINEAGE2)
levels_LINEAGE2_plot = rev(c("B.1.1.7","B.1.617.2","B.1.617.1","P.1","B.1.351","B.1.525","B.1.177+","other")) # "B.1.617.1","B.1.617.2","B.1.617.3"

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
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.1.7")] = "grey60"  
lineage_cols1[which(levels_LINEAGE2_plot=="other")] = "grey70"  
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.177+")] = "grey55"  
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.351")] = muted("cyan")  
lineage_cols1[which(levels_LINEAGE2_plot=="P.1")] = "cyan3"  

data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2_plot)
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2_plot)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)

library(ggplot2)
library(ggthemes)
muller_sanger_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(Sanger Institute baseline surveillance data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sanger_raw1


muller_sangerbyregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(Sanger Institute baseline surveillance data)") +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbyregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
sanger$LINEAGE1 = relevel(sanger$LINEAGE1, ref="B.1.1.7") # we take B.1.1.7 as baseline / reference level
sanger$LINEAGE2 = relevel(sanger$LINEAGE2, ref="B.1.1.7")
set.seed(1)
fit1_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION + DATE_NUM, data=sanger, maxit=1000)
fit2_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * DATE_NUM, data=sanger, maxit=1000)
fit3_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION + ns(DATE_NUM, df=2), data=sanger, maxit=1000)
fit4_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=sanger, maxit=1000)
BIC(fit1_sanger_multi, fit2_sanger_multi, fit3_sanger_multi, fit4_sanger_multi) 
# fit3_sanger_multi fits best (lowest BIC)

# growth rate advantages of different VOCs compared to UK type B.1.1.7 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-05-11
emtrsanger = emtrends(fit3_sanger_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger = data.frame(confint(emtrsanger, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger$contrasts)$p.value)
delta_r_sanger
#               contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1      other - B.1.1.7  0.05062302 0.002155993 NA  0.04639735  0.05484869 0.000000e+00
# 2 (B.1.177+) - B.1.1.7 -0.09142217 0.007611266 NA -0.10633998 -0.07650437 0.000000e+00
# 3    B.1.351 - B.1.1.7  0.02697188 0.004254006 NA  0.01863418  0.03530958 9.954887e-08
# 4    B.1.525 - B.1.1.7  0.02728603 0.006061028 NA  0.01540663  0.03916543 1.594661e-04
# 5  B.1.617.1 - B.1.1.7 -0.07382178 0.017175123 NA -0.10748441 -0.04015916 3.348412e-04
# 6  B.1.617.2 - B.1.1.7  0.10702011 0.004325031 NA  0.09854321  0.11549702 0.000000e+00
# 7        P.1 - B.1.1.7  0.04542234 0.010439079 NA  0.02496212  0.06588256 2.765853e-04

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*4.7)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
# contrast  estimate  asymp.LCL   asymp.UCL
# 1      other - B.1.1.7  26.86181  24.367108  29.40655
# 2 (B.1.177+) - B.1.1.7 -53.67722 -64.838812 -43.27140
# 3    B.1.351 - B.1.1.7  13.51535   9.153031  18.05200
# 4    B.1.525 - B.1.1.7  13.68308   7.509732  20.21090
# 5  B.1.617.1 - B.1.1.7 -41.47635 -65.727835 -20.77366
# 6  B.1.617.2 - B.1.1.7  65.36658  58.907658  72.08804
# 7        P.1 - B.1.1.7  23.79849  12.448143  36.29453

# so this would estimate that B.1.617.2 had a 65% transmission advantage over B.1.1.7 [59%-72%] 95% CLs

# pairwise growth rate advantages for all strain comparisons (i.e. pairwise differences in growth rate per day among the different lineages)
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-05-11
emtrsanger_pairw = emtrends(fit3_sanger_multi, pairwise ~ LINEAGE2,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_pairw = data.frame(confint(emtrsanger_pairw, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger_pairw$contrasts)$p.value)
delta_r_sanger_pairw
#                  contrast      estimate          SE df     asymp.LCL    asymp.UCL      p.value
# 1         B.1.1.7 - other -0.0506230219 0.002155993 NA -0.05484869 -0.046397352 0.000000e+00
# 2    B.1.1.7 - (B.1.177+)  0.0914221741 0.007611266 NA  0.07650437  0.106339982 0.000000e+00
# 3       B.1.1.7 - B.1.351 -0.0269718823 0.004254006 NA -0.03530958 -0.018634184 3.937746e-07
# 4       B.1.1.7 - B.1.525 -0.0272860316 0.006061028 NA -0.03916543 -0.015406635 5.993595e-04
# 5     B.1.1.7 - B.1.617.1  0.0738217821 0.017175123 NA  0.04015916  0.107484406 1.241841e-03
# 6     B.1.1.7 - B.1.617.2 -0.1070201144 0.004325031 NA -0.11549702 -0.098543208 0.000000e+00
# 7           B.1.1.7 - P.1 -0.0454223413 0.010439079 NA -0.06588256 -0.024962123 1.029616e-03
# 8      other - (B.1.177+)  0.1420451960 0.007884556 NA  0.12659175  0.157498642 0.000000e+00
# 9         other - B.1.351  0.0236511396 0.004754926 NA  0.01433166  0.032970623 1.024499e-04
# 10        other - B.1.525  0.0233369902 0.006421547 NA  0.01075099  0.035922990 1.123165e-02
# 11      other - B.1.617.1  0.1244448040 0.017302143 NA  0.09053323  0.158356381 9.870355e-09
# 12      other - B.1.617.2 -0.0563970925 0.004782688 NA -0.06577099 -0.047023196 0.000000e+00
# 13            other - P.1  0.0052006806 0.010651047 NA -0.01567499  0.026076348 9.996846e-01
# 14   (B.1.177+) - B.1.351 -0.1183940564 0.008718809 NA -0.13548261 -0.101305505 0.000000e+00
# 15   (B.1.177+) - B.1.525 -0.1187082058 0.009727462 NA -0.13777368 -0.099642731 0.000000e+00
# 16 (B.1.177+) - B.1.617.1 -0.0176003920 0.018786012 NA -0.05442030  0.019219515 9.813731e-01
# 17 (B.1.177+) - B.1.617.2 -0.1984422885 0.008753680 NA -0.21559919 -0.181285391 0.000000e+00
# 18       (B.1.177+) - P.1 -0.1368445154 0.012919616 NA -0.16216650 -0.111522534 0.000000e+00
# 19      B.1.351 - B.1.525 -0.0003141494 0.007387100 NA -0.01479260  0.014164301 1.000000e+00
# 20    B.1.351 - B.1.617.1  0.1007936644 0.017673222 NA  0.06615479  0.135432543 5.634364e-06
# 21    B.1.351 - B.1.617.2 -0.0800482321 0.006025500 NA -0.09185799 -0.068238469 0.000000e+00
# 22          B.1.351 - P.1 -0.0184504590 0.011240571 NA -0.04048157  0.003580656 7.238705e-01
# 23    B.1.525 - B.1.617.1  0.1011078137 0.018200906 NA  0.06543469  0.136780933 1.030024e-05
# 24    B.1.525 - B.1.617.2 -0.0797340827 0.007415384 NA -0.09426797 -0.065200198 0.000000e+00
# 25          B.1.525 - P.1 -0.0181363096 0.012052359 NA -0.04175850  0.005485879 8.025672e-01
# 26  B.1.617.1 - B.1.617.2 -0.1808418965 0.017666089 NA -0.21546679 -0.146216999 0.000000e+00
# 27        B.1.617.1 - P.1 -0.1192441234 0.020048572 NA -0.15853860 -0.079949644 2.052665e-06
# 28        B.1.617.2 - P.1  0.0615977731 0.011238843 NA  0.03957005  0.083625501 1.390676e-05


# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME BY REGION ####

gentime = 4.7 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantages of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit4_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=sanger, maxit=1000) here
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-05-11
emtrsanger4 = emtrends(fit4_sanger_multi, trt.vs.ctrl1 ~ LINEAGE2, by=c("DATE_NUM", "REGION"), 
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=seq(as.numeric(as.Date("2021-05-01")),as.numeric(as.Date("2021-06-14")))))
delta_r_sanger4 = data.frame(confint(emtrsanger4, 
                                    adjust="none", df=NA)$contrasts)
delta_r_sanger4 = delta_r_sanger4[delta_r_sanger4$contrast=="B.1.617.2 - B.1.1.7",]
delta_r_sanger4

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger4 = delta_r_sanger4
transmadv_sanger4$estimate =  100*(exp(transmadv_sanger4$estimate*4.7)-1)
transmadv_sanger4$asymp.LCL =  100*(exp(transmadv_sanger4$asymp.LCL*4.7)-1)
transmadv_sanger4$asymp.UCL =  100*(exp(transmadv_sanger4$asymp.UCL*4.7)-1)
transmadv_sanger4$collection_date = as.Date(transmadv_sanger4$DATE_NUM, origin="1970-01-01")
transmadv_sanger4$REGION = factor(transmadv_sanger4$REGION, levels=levels_REGION)

plot_sanger_mfit4_transmadv = qplot(data=transmadv_sanger4, 
                               x=collection_date, y=estimate, ymin=asymp.LCL, ymax=asymp.UCL, colour=REGION, fill=REGION, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(fill=REGION, colour=NULL), alpha=I(0.3)) +
  geom_line(aes(colour=REGION, fill=NULL)) +
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
  theme(legend.position = "right") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) # +
  # coord_cartesian(xlim=c(as.Date("2021-05-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit4_transmadv

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit4_transm advantage.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit4_transm advantage.pdf"), width=8, height=6)




# PLOT MULTINOMIAL FIT ####

# extrapolate = 60
date.from = as.numeric(as.Date("2021-01-01"))
date.to = as.numeric(as.Date("2021-06-14")) # max(sanger$DATE_NUM)+extrapolate

fit_sanger_multi_predsbyregion = data.frame(emmeans(fit3_sanger_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "REGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to)), 
                                                  mode="prob", df=NA))
fit_sanger_multi_predsbyregion$collection_date = as.Date(fit_sanger_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion$LINEAGE2 = factor(fit_sanger_multi_predsbyregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_sanger_multi_predsbyregion$REGION = factor(fit_sanger_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of the UK today
fit_sanger_multi_predsbyregion[fit_sanger_multi_predsbyregion$collection_date==today,]

fit_sanger_multi_preds = data.frame(emmeans(fit3_sanger_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to)), 
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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",date.to)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(Sanger Institute baseline surveillance data)")
muller_sanger_mfit

library(ggpubr)
ggarrange(muller_sanger_raw1+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sanger_mfit+ggtitle("Multinomial fit"), ncol=1)


muller_sangerbyregion_mfit = ggplot(data=fit_sanger_multi_predsbyregion, 
                                  aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(c("2021-06-14")), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-14")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbyregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_multinom fit.pdf"), width=8, height=6)


ggarrange(muller_sangerbyregion_raw1+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbyregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_multinom fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_multinom fit multipanel.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & 95% CONFIDENCE INTERVALS

# SGTF data UK to superimpose on model predictions as independent estimate of prop of B.1.1.7
sgtf_UK = read.csv(".\\data\\SGTF_UK\\sgtf_2021-05-17_by adm region.csv")
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
range(sgtf_UK$collection_date) # "2021-01-01" "2021-05-11"
library(tidyr)
sgtf_UK_long = gather(sgtf_UK, Sdropout, count, sgtf:non_sgtf, factor_key=TRUE)
sgtf_UK_long$collection_date

# plot of SGTF data by region
ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-02-20"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
  facet_wrap(~REGION, scales="free_y") +
  geom_col(position="fill") +
  scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                 lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                   lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  ylab("Share") +
  ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN THE UK\n(data PHE)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     expand=c(0,0)) +
  xlab("Collection date")
  
ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by region_stacked as proportion.png"), width=10, height=8)
ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by region_stacked as proportion.pdf"), width=10, height=8)

ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-02-20"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
  # facet_wrap(~REGION, scales="free_y") +
  geom_col(position="fill") +
  scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                 lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")],
                                   lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
  ylab("Share") +
  ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN THE UK\n(data PHE)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     expand=c(0,0)) +
  xlab("Collection date")

ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.png"), width=10, height=8)
ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.pdf"), width=10, height=8)



# logistic spline fit to SGTF data

fit_SGTF = glm(cbind(sgtf,non_sgtf) ~ REGION * ns(DATE_NUM, df=3), family=binomial, data=sgtf_UK) # at level of regions
BIC(fit_SGTF)

SGTF_predsbyregion = data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                                                      at=list(DATE_NUM=seq(date.from, date.to)), 
                                                      type="response"))
SGTF_predsbyregion$collection_date = as.Date(SGTF_predsbyregion$DATE_NUM, origin="1970-01-01")
SGTF_predsbyregion$LINEAGE2 = "B.1.1.7 (S dropout)"
SGTF_predsbyregion$LINEAGE2 = factor(SGTF_predsbyregion$LINEAGE2, levels=c(levels_LINEAGE2_plot, "B.1.1.7 (S dropout)"))

# fitted S dropout in different parts of the UK today
# 65% [61%-68%] now estimated to be S positive across all regions
1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM,
        at=list(DATE_NUM=today_num), 
        type="response"))[,c(2,6,5)]
#        prob asymp.UCL asymp.LCL
# 1 0.8109347  0.781617 0.8371366

# here given by region:
data.frame(REGION=data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                   at=list(DATE_NUM=today_num), 
                   type="response"))[,2], 
           1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                        at=list(DATE_NUM=today_num), 
                        type="response"))[,c(3,7,6)])
# REGION      prob  asymp.UCL asymp.LCL
# 1                   London 0.8069710 0.7599349 0.8466509
# 2               North West 0.9554706 0.9474107 0.9623443
# 3               South West 0.9746051 0.9308922 0.9909374
# 4               South East 0.6782867 0.5945573 0.7519390
# 5          East of England 0.8242979 0.7739807 0.8653620
# 6            East Midlands 0.7518362 0.7006162 0.7968343
# 7            West Midlands 0.5192187 0.4374016 0.6000178
# 8               North East 0.7729366 0.5795222 0.8937018
# 9 Yorkshire and The Humber 0.5643678 0.4461324 0.6757100

# fitted prop of different LINEAGES in the UK today
# 62% [58%-65%] now estimated to be B.1.617.2 across all regions
multinom_preds_today_avg = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|1,
                        at=list(DATE_NUM=today_num), 
                        mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1   B.1.1.7 3.562940e-01 1.656303e-02 NA  3.238311e-01 3.887570e-01
# 2     other 1.842815e-02 1.994328e-03 NA  1.451934e-02 2.233696e-02
# 3  B.1.177+ 1.224589e-06 6.731637e-07 NA -9.478783e-08 2.543966e-06
# 4   B.1.351 3.486894e-03 7.008081e-04 NA  2.113335e-03 4.860452e-03
# 5   B.1.525 1.534008e-03 4.507621e-04 NA  6.505304e-04 2.417485e-03
# 6 B.1.617.1 3.739651e-04 1.788931e-04 NA  2.334099e-05 7.245892e-04
# 7 B.1.617.2 6.169492e-01 1.774651e-02 NA  5.821667e-01 6.517317e-01
# 8       P.1 2.932535e-03 1.127743e-03 NA  7.221988e-04 5.142872e-03

# 64% [60%-69%] non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.6437060 0.6001953 0.6872166

# here given by region:
multinom_preds_today_byregion = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|DATE_NUM, by=c("REGION"),
                        at=list(DATE_NUM=today_num), 
                        mode="prob", df=NA))
multinom_preds_today_byregion
#     LINEAGE2 DATE_NUM                   REGION          prob           SE df      asymp.LCL     asymp.UCL
# 1    B.1.1.7    18771                   London  1.366028e-01 1.147808e-02 NA   1.141062e-01  1.590995e-01
# 2      other    18771                   London  7.614357e-03 1.082641e-03 NA   5.492419e-03  9.736295e-03
# 3   B.1.177+    18771                   London  5.413599e-08 3.010680e-08 NA  -4.872255e-09  1.131442e-07
# 4    B.1.351    18771                   London  1.001842e-02 2.123047e-03 NA   5.857320e-03  1.417951e-02
# 5    B.1.525    18771                   London  2.150457e-03 7.009381e-04 NA   7.766435e-04  3.524270e-03
# 6  B.1.617.1    18771                   London  8.316339e-04 4.069042e-04 NA   3.411637e-05  1.629151e-03
# 7  B.1.617.2    18771                   London  8.344130e-01 1.390087e-02 NA   8.071678e-01  8.616582e-01
# 8        P.1    18771                   London  8.369218e-03 3.090278e-03 NA   2.312384e-03  1.442605e-02
# 9    B.1.1.7    18771               North West  1.309848e-01 9.211402e-03 NA   1.129308e-01  1.490388e-01
# 10     other    18771               North West  9.692227e-03 1.214685e-03 NA   7.311488e-03  1.207297e-02
# 11  B.1.177+    18771               North West  5.672465e-07 3.159961e-07 NA  -5.209454e-08  1.186588e-06
# 12   B.1.351    18771               North West  8.672885e-04 2.204839e-04 NA   4.351480e-04  1.299429e-03
# 13   B.1.525    18771               North West  8.596248e-04 2.675127e-04 NA   3.353095e-04  1.383940e-03
# 14 B.1.617.1    18771               North West  1.072490e-05 8.966761e-06 NA  -6.849627e-06  2.829943e-05
# 15 B.1.617.2    18771               North West  8.574394e-01 1.000155e-02 NA   8.378367e-01  8.770420e-01
# 16       P.1    18771               North West  1.453819e-04 1.121337e-04 NA  -7.439613e-05  3.651600e-04
# 17   B.1.1.7    18771               South West  2.661043e-01 4.264466e-02 NA   1.825223e-01  3.496863e-01
# 18     other    18771               South West  1.280385e-02 3.330009e-03 NA   6.277156e-03  1.933055e-02
# 19  B.1.177+    18771               South West  7.203856e-07 4.162658e-07 NA  -9.548032e-08  1.536252e-06
# 20   B.1.351    18771               South West  2.469060e-03 1.547055e-03 NA  -5.631124e-04  5.501233e-03
# 21   B.1.525    18771               South West  1.971739e-03 1.301663e-03 NA  -5.794724e-04  4.522951e-03
# 22 B.1.617.1    18771               South West  4.608288e-04 3.970111e-04 NA  -3.172986e-04  1.238956e-03
# 23 B.1.617.2    18771               South West  7.098933e-01 4.636696e-02 NA   6.190157e-01  8.007709e-01
# 24       P.1    18771               South West  6.296220e-03 4.992188e-03 NA  -3.488289e-03  1.608073e-02
# 25   B.1.1.7    18771               South East  2.804334e-01 2.289085e-02 NA   2.355681e-01  3.252986e-01
# 26     other    18771               South East  1.645976e-02 2.437952e-03 NA   1.168146e-02  2.123806e-02
# 27  B.1.177+    18771               South East  1.595561e-07 8.903579e-08 NA  -1.495080e-08  3.340631e-07
# 28   B.1.351    18771               South East  4.623015e-03 1.319692e-03 NA   2.036467e-03  7.209564e-03
# 29   B.1.525    18771               South East  3.657200e-03 1.276171e-03 NA   1.155951e-03  6.158448e-03
# 30 B.1.617.1    18771               South East  5.103280e-04 2.801018e-04 NA  -3.866143e-05  1.059317e-03
# 31 B.1.617.2    18771               South East  6.848631e-01 2.559532e-02 NA   6.346972e-01  7.350290e-01
# 32       P.1    18771               South East  9.453049e-03 3.940919e-03 NA   1.728989e-03  1.717711e-02
# 33   B.1.1.7    18771          East of England  2.427291e-01 1.721326e-02 NA   2.089917e-01  2.764665e-01
# 34     other    18771          East of England  1.567742e-02 2.169165e-03 NA   1.142594e-02  1.992891e-02
# 35  B.1.177+    18771          East of England  1.833925e-07 1.021706e-07 NA  -1.685813e-08  3.836432e-07
# 36   B.1.351    18771          East of England  2.324570e-03 6.885246e-04 NA   9.750871e-04  3.674054e-03
# 37   B.1.525    18771          East of England  5.723371e-04 2.945396e-04 NA  -4.949832e-06  1.149624e-03
# 38 B.1.617.1    18771          East of England  1.814426e-04 1.042493e-04 NA  -2.288230e-05  3.857675e-04
# 39 B.1.617.2    18771          East of England  7.371803e-01 1.857925e-02 NA   7.007656e-01  7.735949e-01
# 40       P.1    18771          East of England  1.334689e-03 7.785383e-04 NA  -1.912177e-04  2.860596e-03
# 41   B.1.1.7    18771            East Midlands  3.522568e-01 2.084701e-02 NA   3.113974e-01  3.931161e-01
# 42     other    18771            East Midlands  1.363863e-02 2.083867e-03 NA   9.554322e-03  1.772293e-02
# 43  B.1.177+    18771            East Midlands  1.279050e-06 7.107241e-07 NA  -1.139436e-07  2.672044e-06
# 44   B.1.351    18771            East Midlands  2.114196e-03 6.301117e-04 NA   8.792001e-04  3.349192e-03
# 45   B.1.525    18771            East Midlands  3.474801e-04 2.193344e-04 NA  -8.240738e-05  7.773676e-04
# 46 B.1.617.1    18771            East Midlands  4.713902e-04 2.333154e-04 NA   1.410047e-05  9.286798e-04
# 47 B.1.617.2    18771            East Midlands  6.308817e-01 2.178300e-02 NA   5.881878e-01  6.735756e-01
# 48       P.1    18771            East Midlands  2.885596e-04 3.022265e-04 NA  -3.037934e-04  8.809127e-04
# 49   B.1.1.7    18771            West Midlands  4.806421e-01 2.556019e-02 NA   4.305451e-01  5.307392e-01
# 50     other    18771            West Midlands  2.875843e-02 3.696365e-03 NA   2.151369e-02  3.600317e-02
# 51  B.1.177+    18771            West Midlands  1.433428e-06 7.933745e-07 NA  -1.215579e-07  2.988413e-06
# 52   B.1.351    18771            West Midlands  5.256496e-03 1.370193e-03 NA   2.570968e-03  7.942024e-03
# 53   B.1.525    18771            West Midlands  2.927827e-03 1.029577e-03 NA   9.098924e-04  4.945761e-03
# 54 B.1.617.1    18771            West Midlands  8.197769e-04 4.074620e-04 NA   2.116602e-05  1.618388e-03
# 55 B.1.617.2    18771            West Midlands  4.810882e-01 2.740381e-02 NA   4.273777e-01  5.347987e-01
# 56       P.1    18771            West Midlands  5.057011e-04 5.273958e-04 NA  -5.279757e-04  1.539378e-03
# 57   B.1.1.7    18771               North East  5.789721e-01 5.835761e-02 NA   4.645933e-01  6.933510e-01
# 58     other    18771               North East  2.014383e-02 4.055528e-03 NA   1.219514e-02  2.809252e-02
# 59  B.1.177+    18771               North East  2.169893e-06 1.214558e-06 NA  -2.105978e-07  4.550383e-06
# 60   B.1.351    18771               North East  1.784620e-03 8.992503e-04 NA   2.212199e-05  3.547118e-03
# 61   B.1.525    18771               North East  1.318641e-16 3.964742e-17 NA   5.415662e-17  2.095717e-16
# 62 B.1.617.1    18771               North East  7.391342e-13 1.052796e-12 NA  -1.324309e-12  2.802577e-12
# 63 B.1.617.2    18771               North East  3.990972e-01 6.049291e-02 NA   2.805333e-01  5.176612e-01
# 64       P.1    18771               North East  3.208544e-16 1.265211e-16 NA   7.287769e-17  5.688311e-16
# 65   B.1.1.7    18771 Yorkshire and The Humber  7.379208e-01 2.764448e-02 NA   6.837386e-01  7.921030e-01
# 66     other    18771 Yorkshire and The Humber  4.106485e-02 5.052968e-03 NA   3.116121e-02  5.096848e-02
# 67  B.1.177+    18771 Yorkshire and The Humber  4.454213e-06 2.443762e-06 NA  -3.354732e-07  9.243899e-06
# 68   B.1.351    18771 Yorkshire and The Humber  1.924381e-03 6.229977e-04 NA   7.033283e-04  3.145435e-03
# 69   B.1.525    18771 Yorkshire and The Humber  1.319406e-03 5.443911e-04 NA   2.524192e-04  2.386393e-03
# 70 B.1.617.1    18771 Yorkshire and The Humber  7.956060e-05 5.553480e-05 NA  -2.928561e-05  1.884068e-04
# 71 B.1.617.2    18771 Yorkshire and The Humber  2.176866e-01 2.893627e-02 NA   1.609725e-01  2.744006e-01
# 72       P.1    18771 Yorkshire and The Humber 2.713376e-290 0.000000e+00 NA  2.713376e-290 2.713376e-290



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

# plots of multinomial fit to Sanger Inst data 

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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-31")), 
                     expand=c(0,0)) +
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit_logit scale.pdf"), width=8, height=6)
library(svglite)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit_logit scale.svg"), width=8, height=6)

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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-30")), 
                     expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",date.to-1)),
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit_response scale.pdf"), width=8, height=6)


# plots of multinomial fit to Sanger Inst data with PHE S dropout data overlaid

plot_sanger_mfit_logit = qplot(data=rbind(fit_sanger_multi_predsbyregion2, SGTF_predsbyregion), 
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-31")), 
                     expand=c(0,0)) +
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_logit scale.pdf"), width=8, height=6)
library(svglite)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_logit scale.svg"), width=8, height=6)
write.csv(fit_sanger_multi_predsbyregion2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_multinomial fit sanger inst data.csv"), row.names=F)
write.csv(SGTF_predsbyregion, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_logistic fit S dropout data.csv"), row.names=F)
write.csv(data_agbyweekregion1, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_sanger inst data aggregated by week.csv"), row.names=F)
write.csv(sgtf_UK2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_S dropout data.csv"), row.names=F)

# on response scale:
plot_sanger_mfit = qplot(data=rbind(fit_sanger_multi_predsbyregion2, SGTF_predsbyregion), 
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-30")), 
                     expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",date.to-1)),
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


# project multinomial fit onto case data ####
# case data from https://coronavirus.data.gov.uk/details/download

# cases_uk_ltla = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesByPublishDateRollingRate&metric=newCasesBySpecimenDate&format=csv")
# cases_uk_ltla$date = as.Date(cases_uk_ltla$date)
# cases_uk_ltla$DATE_NUM = as.numeric(cases_uk_ltla$date)
# cases_uk_ltla$LTLA = factor(cases_uk_ltla$areaCode)
# cases_uk_ltla$areaCode = NULL
# cases_uk_ltla$REGION = LTLAs_regions$RGN20NM[match(cases_uk_ltla$LTLA, LTLAs_regions$LAD20CD)]
# cases_uk_ltla = cases_uk_ltla[!is.na(cases_uk_ltla$REGION),]
# cases_uk_ltla$REGION = factor(cases_uk_ltla$REGION, levels=levels_REGION)
cases_uk_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesByPublishDateRollingRate&metric=newCasesBySpecimenDate&format=csv")
cases_uk_region$date = as.Date(cases_uk_region$date)
cases_uk_region$DATE_NUM = as.numeric(cases_uk_region$date)
cases_uk_region$REGION = factor(cases_uk_region$areaName, levels=levels_REGION)
cases_uk_region$areaName = NULL

fit_sanger_multi_predsbyregion$totcases = cases_uk_region$newCasesBySpecimenDate[match(interaction(fit_sanger_multi_predsbyregion$collection_date,
                                                                                                    fit_sanger_multi_predsbyregion$REGION),
                                                                                        interaction(cases_uk_region$date,
                                                                                                    cases_uk_region$REGION))]
fit_sanger_multi_predsbyregion$WEEKDAY = factor(weekdays(fit_sanger_multi_predsbyregion$collection_date))
gamfittotcases = gam(totcases ~ s(DATE_NUM, k=5, by=REGION) + WEEKDAY, family=poisson(log), data=fit_sanger_multi_predsbyregion)
BIC(gamfittotcases) 
gamfitemmeans = as.data.frame(emmeans(gamfittotcases, ~ DATE_NUM, by=c("REGION","DATE_NUM"),
                                                                         at=list(REGION=levels(fit_sanger_multi_predsbyregion$REGION),
                                                                           DATE_NUM=seq(min(fit_sanger_multi_predsbyregion$DATE_NUM),
                                                         max(fit_sanger_multi_predsbyregion$DATE_NUM))), type="response"))
fit_sanger_multi_predsbyregion$totcases_smoothed = gamfitemmeans$rate[match(interaction(fit_sanger_multi_predsbyregion$REGION,fit_sanger_multi_predsbyregion$DATE_NUM),
                                                                            interaction(gamfitemmeans$REGION,gamfitemmeans$DATE_NUM))]
fit_sanger_multi_predsbyregion$cases = fit_sanger_multi_predsbyregion$totcases * fit_sanger_multi_predsbyregion$prob
fit_sanger_multi_predsbyregion$cases_smoothed = fit_sanger_multi_predsbyregion$totcases_smoothed * fit_sanger_multi_predsbyregion$prob
fit_sanger_multi_predsbyregion$REGION = factor(fit_sanger_multi_predsbyregion$REGION, levels=levels_REGION)
fit_sanger_multi_predsbyregion$cases[fit_sanger_multi_predsbyregion$cases<=1] = NA
fit_sanger_multi_predsbyregion$cases_smoothed[fit_sanger_multi_predsbyregion$cases_smoothed<=1] = NA
cases_uk_region$totcases_smoothed = gamfitemmeans$rate[match(interaction(cases_uk_region$REGION,cases_uk_region$DATE_NUM),
                                                             interaction(gamfitemmeans$REGION,gamfitemmeans$DATE_NUM))]

# plot new cases per day by region
ggplot(data=fit_sanger_multi_predsbyregion,
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
                     limits=c(as.Date("2021-01-01"),max(cases_uk_region$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE UK") +
  scale_y_log10() +
  scale_fill_manual("variant", values=c(lineage_cols1,I("black"))) +
  scale_colour_manual("variant", values=c(lineage_cols1,I("black"))) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day by lineage multinomial fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day by lineage multinomial fit.pdf"), width=8, height=6)


ggplot(data=fit_sanger_multi_predsbyregion, 
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE UK") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion, 
       aes(x=collection_date, y=cases_smoothed, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-03-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN THE UK") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.pdf"), width=8, height=6)



# some further tests

SGTF_predsbyregion$newCasesBySpecimenDate = cases_uk_region$newCasesBySpecimenDate[match(interaction(SGTF_predsbyregion$collection_date, SGTF_predsbyregion$REGION),
                                                                                         interaction(cases_uk_region$date, cases_uk_region$REGION))]
SGTF_predsbyregion$newCasesByPublishDateRollingRate = cases_uk_region$newCasesByPublishDateRollingRate[match(interaction(SGTF_predsbyregion$collection_date, SGTF_predsbyregion$REGION),
                                                                                                   interaction(cases_uk_region$date, cases_uk_region$REGION))]
SGTF_predsbyregion$weekday = factor(weekdays(SGTF_predsbyregion$collection_date))
SGTF_predsbyregion$newPCRTestsByPublishDate = NULL
SGTF_predsbyregion3 = SGTF_predsbyregion[complete.cases(SGTF_predsbyregion),]  
SGTF_predsbyregion3 = SGTF_predsbyregion3[SGTF_predsbyregion3$collection_date==max(SGTF_predsbyregion3$collection_date),]

# Calculate Re values from intrinsic growth rate
# BETTER: from epiforecasts package growth_to_R
# https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (better if distribution is known to be gamma, based on the full integral)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) {
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}
library(mgcv)

qplot(data=SGTF_predsbyregion, x=collection_date, y=newCasesBySpecimenDate, geom="line") +
  facet_wrap(~ REGION) + scale_y_log10()

fit_cases = gam(newCasesBySpecimenDate ~ s(DATE_NUM, bs="cs", k=7, fx=F, by=REGION) + weekday,
                family=poisson(log), data=SGTF_predsbyregion[SGTF_predsbyregion$collection_date>=as.Date("2021-01-01"),],
)
BIC(fit_cases)
# calculate instantaneous growth rates & 95% CLs using emtrends 
# based on the slope of the GAM fit on a log link scale
extrapolate = 0
df_r = as.data.frame(emtrends(fit_cases, ~DATE_NUM|REGION, var="DATE_NUM", 
                              at=list(DATE_NUM=seq(min(SGTF_predsbyregion$DATE_NUM),
                                                   max(SGTF_predsbyregion$DATE_NUM)+extrapolate),
                                      REGION=levels(SGTF_predsbyregion$REGION)
                              ),
                              type="link"))
colnames(df_r)[3] = "r"
colnames(df_r)[6] = "r_LOWER"
colnames(df_r)[7] = "r_UPPER"
df_r$DATE = as.Date(df_r$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
df_r$Re = Re.from.r(df_r$r)
df_r$Re_LOWER = Re.from.r(df_r$r_LOWER)
df_r$Re_UPPER = Re.from.r(df_r$r_UPPER)
df_r = df_r[complete.cases(df_r),]
df_r$REGION = factor(df_r$REGION, levels=levels_REGION)
qplot(data=df_r, x=DATE, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue"), group=REGION) +
  facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J")) +
  scale_y_continuous(limits=c(1/3,3), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re IN THE UK BY REGION BASED ON NEW CONFIRMED CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))



