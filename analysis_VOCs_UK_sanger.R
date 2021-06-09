# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN ENGLAND BASED ON SANGER INSTITUTE BASELINE SURVEILLANCE SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & most surge testing/active surveillance)

# last update 8 JUNE 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-31")
today_num = as.numeric(today)
today # "2021-05-31"
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
range(sanger$WeekEndDate) # "2021-01-02" "2021-05-29"
sanger$Week = lubridate::week(sanger$WeekEndDate)
sanger$DATE_NUM = as.numeric(sanger$WeekEndDate)-3.5 # using week midpoint
colnames(sanger)

sanger = sanger[rep(seq_len(nrow(sanger)), sanger$Count),] # convert to long format
sanger$Count = NULL
nrow(sanger) # 163504

nrow(sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 18462 (last 2 weeks)
nrow(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 12530 (last 2 weeks)
nrow(sanger[grepl("B.1.1.7",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 5476

length(unique(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),"LTLA"])) # 269 LTLAs

unique(sanger$Lineage)
sel_target_VOC = "B.1.617"
table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
sanger[sanger$Lineage=="B.1.617",]
sum(table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])) # 14726 B.1.617

unique(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
# "B.1.617"   "B.1.617.1" "B.1.617.3" "B.1.617.2"
sum(grepl("B.1.617",sanger$Lineage, fixed=TRUE)) # 14726 B.1.617+
sanger$LINEAGE1 = sanger$Lineage
sanger$LINEAGE2 = sanger$Lineage
sanger[grepl(sel_target_VOC, sanger$LINEAGE1, fixed=T),"LINEAGE1"] = paste0(sel_target_VOC,"+")
sel_target_VOC = paste0(sel_target_VOC, "+")
sanger[grepl("B.1.177", sanger$LINEAGE1, fixed=T),"LINEAGE1"] = "B.1.177+"
sanger[grepl("B.1.177", sanger$LINEAGE2, fixed=T),"LINEAGE2"] = "B.1.177+"

sum("B.1.1.7"==sanger$LINEAGE1) # 140961 B.1.1.7


table_lineage = as.data.frame(table(sanger$LINEAGE1))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger$LINEAGE1))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# Lineage        Prop
# 8         B.1 0.001646760
# 73    B.1.1.7 0.896664189
# 83   B.1.177+ 0.035446996
# 105   B.1.351 0.001890250
# 151 B.1.617.1 0.001044443
# 152 B.1.617.2 0.052625846

sel_ref_lineage = "B.1.1.7"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"Lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.1.7","B.1.177+","B.1.351", "B.1.525","B.1.617+","B.1.617.1","B.1.617.2","P.1")

sanger$LINEAGE1[!(sanger$LINEAGE1 %in% sel_lineages)] = "other"
sanger$LINEAGE2[!(sanger$LINEAGE2 %in% sel_lineages)] = "other"
# sanger = sanger[sanger$LINEAGE1 %in% sel_lineages, ]
sum(table(sanger$LINEAGE1))
table(sanger$LINEAGE1)
# B.1.1.7  B.1.177+   B.1.351   B.1.525 B.1.617.1 B.1.617.2     other       P.1 
# 139937      5532       295       150       163      8213      1706        68 
table(sanger$LINEAGE2)
# B.1.1.7  B.1.177+   B.1.351   B.1.525 B.1.617.1 B.1.617.2     other       P.1 
# 139937      5532       295       150       163      8213      1706        68 

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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-05-25
emtrsanger = emtrends(fit3_sanger_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger = data.frame(confint(emtrsanger, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger$contrasts)$p.value)
delta_r_sanger
#               contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1      other - B.1.1.7  0.06896361 0.001681382 NA  0.065668165  0.07225906 0.000000e+00
# 2 (B.1.177+) - B.1.1.7 -0.08595069 0.008784854 NA -0.103168684 -0.06873269 3.774758e-15
# 3    B.1.351 - B.1.1.7  0.01309689 0.004038305 NA  0.005181960  0.02101182 1.088318e-02
# 4    B.1.525 - B.1.1.7  0.01605618 0.005793807 NA  0.004700530  0.02741184 4.060376e-02
# 5  B.1.617.1 - B.1.1.7 -0.10121464 0.016384817 NA -0.133328288 -0.06910098 1.988111e-07
# 6  B.1.617.2 - B.1.1.7  0.10807769 0.002257249 NA  0.103653561  0.11250181 0.000000e+00
# 7        P.1 - B.1.1.7  0.01443965 0.010019772 NA -0.005198745  0.03407804 5.400375e-01

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 5.5 days, Ferretti et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger =  sign(delta_r_sanger[,c(2,5,6)])*100*(exp(abs(delta_r_sanger[,c(2,5,6)])*5.5)-1)
transmadv_sanger =  data.frame(contrast=delta_r_sanger$contrast, transmadv_sanger)
transmadv_sanger
# contrast  estimate  asymp.LCL   asymp.UCL
# 1      other - B.1.1.7  46.126117   43.501444  48.79880
# 2 (B.1.177+) - B.1.1.7 -60.436618  -76.372449 -45.94064
# 3    B.1.351 - B.1.1.7   7.469071    2.891081  12.25075
# 4    B.1.525 - B.1.1.7   9.232560    2.619000  16.27235
# 5  B.1.617.1 - B.1.1.7 -74.487078 -108.195131 -46.23656
# 6  B.1.617.2 - B.1.1.7  81.199288   76.843431  85.66244
# 7        P.1 - B.1.1.7   8.265685   -2.900580  20.61449

# so this would estimate that B.1.617.2 had a 81% transmission advantage over B.1.1.7 [77%-86%] 95% CLs

# pairwise growth rate advantages for all strain comparisons (i.e. pairwise differences in growth rate per day among the different lineages)
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-05-25
emtrsanger_pairw = emtrends(fit3_sanger_multi, pairwise ~ LINEAGE2,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_pairw = data.frame(confint(emtrsanger_pairw, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger_pairw$contrasts)$p.value)
delta_r_sanger_pairw
#                  contrast      estimate          SE df     asymp.LCL    asymp.UCL      p.value
# 1         B.1.1.7 - other -0.068963614 0.001681382 NA -0.07225906 -0.065668165 0.000000e+00
# 2    B.1.1.7 - (B.1.177+)  0.085950687 0.008784854 NA  0.06873269  0.103168684 0.000000e+00
# 3       B.1.1.7 - B.1.351 -0.013096892 0.004038305 NA -0.02101182 -0.005181960 3.539163e-02
# 4       B.1.1.7 - B.1.525 -0.016056183 0.005793807 NA -0.02741184 -0.004700530 1.177104e-01
# 5     B.1.1.7 - B.1.617.1  0.101214636 0.016384817 NA  0.06910098  0.133328288 7.848433e-07
# 6     B.1.1.7 - B.1.617.2 -0.108077687 0.002257249 NA -0.11250181 -0.103653561 0.000000e+00
# 7           B.1.1.7 - P.1 -0.014439648 0.010019772 NA -0.03407804  0.005198745 8.351285e-01
# 8      other - (B.1.177+)  0.154914301 0.008925255 NA  0.13742112  0.172407480 0.000000e+00
# 9         other - B.1.351  0.055866722 0.004358568 NA  0.04732409  0.064409358 0.000000e+00
# 10        other - B.1.525  0.052907431 0.006020512 NA  0.04110744  0.064707418 0.000000e+00
# 11      other - B.1.617.1  0.170178250 0.016463035 NA  0.13791129  0.202445206 0.000000e+00
# 12      other - B.1.617.2 -0.039114074 0.002689847 NA -0.04438608 -0.033842070 0.000000e+00
# 13            other - P.1  0.054523965 0.010148621 NA  0.03463303  0.074414897 2.146755e-05
# 14   (B.1.177+) - B.1.351 -0.099047579 0.009668707 NA -0.11799790 -0.080097261 0.000000e+00
# 15   (B.1.177+) - B.1.525 -0.102006870 0.010521851 NA -0.12262932 -0.081384420 0.000000e+00
# 16 (B.1.177+) - B.1.617.1  0.015263950 0.018591242 NA -0.02117422  0.051702115 9.914028e-01
# 17 (B.1.177+) - B.1.617.2 -0.194028374 0.009068771 NA -0.21180284 -0.176253910 0.000000e+00
# 18       (B.1.177+) - P.1 -0.100390335 0.013325960 NA -0.12650874 -0.074271933 2.147324e-09
# 19      B.1.351 - B.1.525 -0.002959291 0.007046560 NA -0.01677029  0.010851712 9.998847e-01
# 20    B.1.351 - B.1.617.1  0.114311528 0.016859256 NA  0.08126799  0.147355063 5.965809e-08
# 21    B.1.351 - B.1.617.2 -0.094980795 0.004594325 NA -0.10398551 -0.085976083 0.000000e+00
# 22          B.1.351 - P.1 -0.001342757 0.010777243 NA -0.02246576  0.019780251 1.000000e+00
# 23    B.1.525 - B.1.617.1  0.117270819 0.017369657 NA  0.08322692  0.151314721 6.759204e-08
# 24    B.1.525 - B.1.617.2 -0.092021505 0.006194380 NA -0.10416227 -0.079880742 0.000000e+00
# 25          B.1.525 - P.1  0.001616534 0.011559508 NA -0.02103969  0.024272755 9.999999e-01
# 26  B.1.617.1 - B.1.617.2 -0.209292324 0.016532799 NA -0.24169601 -0.176888633 0.000000e+00
# 27        B.1.617.1 - P.1 -0.115654285 0.019165769 NA -0.15321850 -0.078090069 1.430031e-06
# 28        B.1.617.2 - P.1  0.093638039 0.010232744 NA  0.07358223  0.113693849 0.000000e+00


# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME BY REGION ####

gentime = 5.5 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantages of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit4_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=sanger, maxit=1000) here
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-05-25
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
transmadv_sanger4$estimate =  100*(exp(transmadv_sanger4$estimate*5.5)-1)
transmadv_sanger4$asymp.LCL =  100*(exp(transmadv_sanger4$asymp.LCL*5.5)-1)
transmadv_sanger4$asymp.UCL =  100*(exp(transmadv_sanger4$asymp.UCL*5.5)-1)
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
date.to = as.numeric(as.Date("2021-06-30")) # max(sanger$DATE_NUM)+extrapolate

fit_sanger_multi_predsbyregion = data.frame(emmeans(fit3_sanger_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "REGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to)), 
                                                  mode="prob", df=NA))
fit_sanger_multi_predsbyregion$collection_date = as.Date(fit_sanger_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion$LINEAGE2 = factor(fit_sanger_multi_predsbyregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_sanger_multi_predsbyregion$REGION = factor(fit_sanger_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of England today
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
                     limits=as.Date(c("2021-01-01",date.to), origin="1970-01-01"), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN ENGLAND\n(Sanger Institute baseline surveillance data)")
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
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",date.to), origin="1970-01-01"), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbyregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_multinom fit.pdf"), width=8, height=6)

library(ggpubr)
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

# PHE SGTF data to superimpose on model predictions as independent estimate of prop of B.1.1.7
sgtf_UK = read.csv(".\\data\\SGTF_UK\\sgtf_2021-05-31_by adm region.csv")
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
ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-02-20"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
  facet_wrap(~REGION, scales="free_y") +
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
  ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN ENGLAND\n(data PHE)") +
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

# fitted S dropout in different parts of England today
# 91% [89%-93%] now estimated to be S positive across all regions
1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM,
        at=list(DATE_NUM=today_num), 
        type="response"))[,c(2,6,5)]
#        prob asymp.UCL asymp.LCL
# 1 0.8817344 0.8698647 0.8926551

# here given by region:
data.frame(REGION=data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                   at=list(DATE_NUM=today_num), 
                   type="response"))[,2], 
           1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
                        at=list(DATE_NUM=today_num), 
                        type="response"))[,c(3,7,6)])
#                    REGION      prob  asymp.UCL asymp.LCL
# 1                   London 0.8444509 0.8229577 0.8637667
# 2               North West 0.9705632 0.9675806 0.9732791
# 3               South West 0.9381692 0.8672523 0.9724061
# 4               South East 0.9191101 0.9016632 0.9336892
# 5          East of England 0.8987968 0.8829933 0.9126769
# 6            East Midlands 0.8330039 0.8108668 0.8530197
# 7            West Midlands 0.8535572 0.8333523 0.8716896
# 8               North East 0.7631161 0.7103600 0.8088487
# 9 Yorkshire and The Humber 0.7349352 0.7007993 0.7664742

# fitted prop of different LINEAGES in England today
# 76% [74%-78%] now estimated to be B.1.617.2 across all regions
multinom_preds_today_avg = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|1,
                        at=list(DATE_NUM=today_num), 
                        mode="prob", df=NA))
multinom_preds_today_avg
#    LINEAGE2         prob           SE df    asymp.LCL    asymp.UCL
# 1   B.1.1.7 1.896879e-01 5.890849e-03 NA  1.781420e-01 2.012337e-01
# 2     other 2.632919e-02 1.807114e-03 NA  2.278732e-02 2.987107e-02
# 3  B.1.177+ 6.016692e-07 3.830262e-07 NA -1.490483e-07 1.352387e-06
# 4   B.1.351 1.094098e-03 2.095095e-04 NA  6.834668e-04 1.504729e-03
# 5   B.1.525 4.816908e-04 1.377317e-04 NA  2.117417e-04 7.516399e-04
# 6 B.1.617.1 9.322959e-05 4.417962e-05 NA  6.639126e-06 1.798201e-04
# 7 B.1.617.2 7.816938e-01 6.535848e-03 NA  7.688838e-01 7.945038e-01
# 8       P.1 6.195226e-04 2.122040e-04 NA  2.036105e-04 1.035435e-03

# 81% [79%-83%] non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.8103121 0.7927764 0.8278479 

# here given by region:
multinom_preds_today_byregion = data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE2|DATE_NUM, by=c("REGION"),
                        at=list(DATE_NUM=today_num), 
                        mode="prob", df=NA))
multinom_preds_today_byregion
#     LINEAGE2 DATE_NUM                   REGION          prob           SE df      asymp.LCL     asymp.UCL
# 1    B.1.1.7    18778                   London  8.410155e-02 4.190842e-03 NA   7.588765e-02  9.231545e-02
# 2      other    18778                   London  1.430649e-02 1.317448e-03 NA   1.172434e-02  1.688864e-02
# 3   B.1.177+    18778                   London  2.971694e-08 1.897985e-08 NA  -7.482884e-09  6.691677e-08
# 4    B.1.351    18778                   London  3.614924e-03 6.976901e-04 NA   2.247476e-03  4.982371e-03
# 5    B.1.525    18778                   London  7.999155e-04 2.465730e-04 NA   3.166412e-04  1.283190e-03
# 6  B.1.617.1    18778                   London  2.259873e-04 1.073639e-04 NA   1.555789e-05  4.364167e-04
# 7  B.1.617.2    18778                   London  8.947837e-01 5.161075e-03 NA   8.846681e-01  9.048992e-01
# 8        P.1    18778                   London  2.167456e-03 7.235347e-04 NA   7.493535e-04  3.585558e-03
# 9    B.1.1.7    18778               North West  5.698448e-02 2.206176e-03 NA   5.266046e-02  6.130851e-02
# 10     other    18778               North West  1.209923e-02 8.965126e-04 NA   1.034210e-02  1.385636e-02
# 11  B.1.177+    18778               North West  2.203474e-07 1.413213e-07 NA  -5.663718e-08  4.973320e-07
# 12   B.1.351    18778               North West  2.250128e-04 5.444072e-05 NA   1.183109e-04  3.317146e-04
# 13   B.1.525    18778               North West  2.418783e-04 7.233452e-05 NA   1.001053e-04  3.836514e-04
# 14 B.1.617.1    18778               North West  1.899050e-06 1.588132e-06 NA  -1.213632e-06  5.011732e-06
# 15 B.1.617.2    18778               North West  9.304223e-01 2.653661e-03 NA   9.252213e-01  9.356234e-01
# 16       P.1    18778               North West  2.492690e-05 1.907605e-05 NA  -1.246146e-05  6.231526e-05
# 17   B.1.1.7    18778               South West  1.751696e-01 1.557978e-02 NA   1.446338e-01  2.057054e-01
# 18     other    18778               South West  2.165801e-02 4.117362e-03 NA   1.358813e-02  2.972789e-02
# 19  B.1.177+    18778               South West  4.245597e-07 2.751738e-07 NA  -1.147709e-07  9.638904e-07
# 20   B.1.351    18778               South West  1.079083e-03 5.763761e-04 NA  -5.059366e-05  2.208759e-03
# 21   B.1.525    18778               South West  6.996568e-04 4.455997e-04 NA  -1.737026e-04  1.573016e-03
# 22 B.1.617.1    18778               South West  1.778873e-04 1.297766e-04 NA  -7.647019e-05  4.322449e-04
# 23 B.1.617.2    18778               South West  7.999958e-01 1.763740e-02 NA   7.654271e-01  8.345644e-01
# 24       P.1    18778               South West  1.219588e-03 9.237674e-04 NA  -5.909633e-04  3.030138e-03
# 25   B.1.1.7    18778               South East  1.297852e-01 7.104323e-03 NA   1.158610e-01  1.437094e-01
# 26     other    18778               South East  2.052055e-02 2.156073e-03 NA   1.629473e-02  2.474638e-02
# 27  B.1.177+    18778               South East  6.563479e-08 4.208464e-08 NA  -1.684959e-08  1.481192e-07
# 28   B.1.351    18778               South East  1.189187e-03 3.281668e-04 NA   5.459923e-04  1.832382e-03
# 29   B.1.525    18778               South East  1.006117e-03 3.405169e-04 NA   3.387161e-04  1.673518e-03
# 30 B.1.617.1    18778               South East  1.153141e-04 6.130332e-05 NA  -4.838199e-06  2.354664e-04
# 31 B.1.617.2    18778               South East  8.457396e-01 8.350073e-03 NA   8.293738e-01  8.621054e-01
# 32       P.1    18778               South East  1.643983e-03 6.509042e-04 NA   3.682337e-04  2.919731e-03
# 33   B.1.1.7    18778          East of England  1.160844e-01 5.402149e-03 NA   1.054964e-01  1.266724e-01
# 34     other    18778          East of England  1.933329e-02 1.879256e-03 NA   1.565002e-02  2.301657e-02
# 35  B.1.177+    18778          East of England  7.781970e-08 4.989098e-08 NA  -1.996483e-08  1.756042e-07
# 36   B.1.351    18778          East of England  6.275441e-04 1.824350e-04 NA   2.699780e-04  9.851102e-04
# 37   B.1.525    18778          East of England  1.644843e-04 8.415464e-05 NA  -4.557311e-07  3.294244e-04
# 38 B.1.617.1    18778          East of England  3.531970e-05 2.025639e-05 NA  -4.382093e-06  7.502150e-05
# 39 B.1.617.2    18778          East of England  8.633811e-01 6.300121e-03 NA   8.510331e-01  8.757291e-01
# 40       P.1    18778          East of England  3.737726e-04 1.872898e-04 NA   6.691296e-06  7.408540e-04
# 41   B.1.1.7    18778            East Midlands  1.931174e-01 7.967119e-03 NA   1.775021e-01  2.087326e-01
# 42     other    18778            East Midlands  1.915905e-02 2.187149e-03 NA   1.487232e-02  2.344578e-02
# 43  B.1.177+    18778            East Midlands  6.232810e-07 3.995550e-07 NA  -1.598325e-07  1.406394e-06
# 44   B.1.351    18778            East Midlands  6.642249e-04 1.967860e-04 NA   2.785315e-04  1.049918e-03
# 45   B.1.525    18778            East Midlands  1.187951e-04 7.458693e-05 NA  -2.739260e-05  2.649828e-04
# 46 B.1.617.1    18778            East Midlands  1.064290e-04 5.271770e-05 NA   3.104234e-06  2.097538e-04
# 47 B.1.617.2    18778            East Midlands  7.867704e-01 8.714325e-03 NA   7.696906e-01  8.038501e-01
# 48       P.1    18778            East Midlands  6.314003e-05 6.597850e-05 NA  -6.617545e-05  1.924555e-04
# 49   B.1.1.7    18778            West Midlands  1.977154e-01 8.338853e-03 NA   1.813715e-01  2.140592e-01
# 50     other    18778            West Midlands  2.986110e-02 2.797755e-03 NA   2.437761e-02  3.534460e-02
# 51  B.1.177+    18778            West Midlands  5.244830e-07 3.355034e-07 NA  -1.330915e-07  1.182058e-06
# 52   B.1.351    18778            West Midlands  1.255797e-03 3.247068e-04 NA   6.193834e-04  1.892211e-03
# 53   B.1.525    18778            West Midlands  7.880455e-04 2.719290e-04 NA   2.550744e-04  1.321017e-03
# 54 B.1.617.1    18778            West Midlands  1.431916e-04 7.133079e-05 NA   3.385779e-06  2.829973e-04
# 55 B.1.617.2    18778            West Midlands  7.701531e-01 9.543437e-03 NA   7.514484e-01  7.888579e-01
# 56       P.1    18778            West Midlands  8.283800e-05 8.644471e-05 NA  -8.659053e-05  2.522665e-04
# 57   B.1.1.7    18778               North East  2.899730e-01 2.779144e-02 NA   2.355028e-01  3.444433e-01
# 58     other    18778               North East  2.713337e-02 4.831818e-03 NA   1.766318e-02  3.660356e-02
# 59  B.1.177+    18778               North East  9.672003e-07 6.241194e-07 NA  -2.560512e-07  2.190452e-06
# 60   B.1.351    18778               North East  4.942796e-04 2.473526e-04 NA   9.477502e-06  9.790817e-04
# 61   B.1.525    18778               North East  1.657826e-94 5.315412e-95 NA   6.160248e-95  2.699628e-94
# 62 B.1.617.1    18778               North East  4.683585e-95 2.267252e-95 NA   2.398520e-96  9.127318e-95
# 63 B.1.617.2    18778               North East  6.823983e-01 3.032323e-02 NA   6.229659e-01  7.418308e-01
# 64       P.1    18778               North East  4.492056e-94 1.561886e-94 NA   1.430816e-94  7.553296e-94
# 65   B.1.1.7    18778 Yorkshire and The Humber  4.642598e-01 1.946216e-02 NA   4.261146e-01  5.024049e-01
# 66     other    18778 Yorkshire and The Humber  7.289164e-02 6.920316e-03 NA   5.932807e-02  8.645521e-02
# 67  B.1.177+    18778 Yorkshire and The Humber  2.481979e-06 1.579912e-06 NA  -6.145914e-07  5.578550e-06
# 68   B.1.351    18778 Yorkshire and The Humber  6.968286e-04 2.265461e-04 NA   2.528063e-04  1.140851e-03
# 69   B.1.525    18778 Yorkshire and The Humber  5.163247e-04 2.146118e-04 NA   9.569326e-05  9.369562e-04
# 70 B.1.617.1    18778 Yorkshire and The Humber  3.303828e-05 2.096394e-05 NA  -8.050284e-06  7.412684e-05
# 71 B.1.617.2    18778 Yorkshire and The Humber  4.615999e-01 2.189414e-02 NA   4.186882e-01  5.045116e-01
# 72       P.1    18778 Yorkshire and The Humber 1.138040e-245 0.000000e+00 NA  1.138040e-245 1.138040e-245


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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-14")), 
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-30")), ylim=c(0.001, 0.998), expand=c(0,0))
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-30")), 
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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-14")), 
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-30")), ylim=c(0.001, 0.998), expand=c(0,0))
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
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-30")), 
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
cases_uk_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDate&format=csv")
cases_uk_region$date = as.Date(cases_uk_region$date)
cases_uk_region$DATE_NUM = as.numeric(cases_uk_region$date)
cases_uk_region$REGION = factor(cases_uk_region$areaName, levels=levels_REGION)
cases_uk_region$areaName = NULL

fit_sanger_multi_predsbyregion$totcases_rollingmean = cases_uk_region$newCasesBySpecimenDateRollingSum[match(interaction(fit_sanger_multi_predsbyregion$collection_date,
                                                                                                    fit_sanger_multi_predsbyregion$REGION),
                                                                                        interaction(cases_uk_region$date,
                                                                                                    cases_uk_region$REGION))]/7
fit_sanger_multi_predsbyregion$totcases = cases_uk_region$newCasesBySpecimenDate[match(interaction(fit_sanger_multi_predsbyregion$collection_date,
                                                                                                             fit_sanger_multi_predsbyregion$REGION),
                                                                                                 interaction(cases_uk_region$date,
                                                                                                             cases_uk_region$REGION))]
fit_sanger_multi_predsbyregion$WEEKDAY = factor(weekdays(fit_sanger_multi_predsbyregion$collection_date))
library(mgcv)
gamfittotcases = gam(totcases ~ s(DATE_NUM, bs="cs", k=7, by=REGION) + WEEKDAY, family=poisson(log), 
                     data=fit_sanger_multi_predsbyregion[fit_sanger_multi_predsbyregion$collection_date<=max(sanger$WeekEndDate),])
BIC(gamfittotcases) 
gamfitemmeans = as.data.frame(emmeans(gamfittotcases, ~ DATE_NUM, by=c("REGION","DATE_NUM"),
                                                                         at=list(REGION=levels(fit_sanger_multi_predsbyregion$REGION),
                                                                           DATE_NUM=seq(min(fit_sanger_multi_predsbyregion$DATE_NUM),
                                                         max(fit_sanger_multi_predsbyregion$DATE_NUM)+14)), type="response"))
fit_sanger_multi_predsbyregion$totcases_smoothed = gamfitemmeans$rate[match(interaction(fit_sanger_multi_predsbyregion$REGION,fit_sanger_multi_predsbyregion$DATE_NUM),
                                                                            interaction(gamfitemmeans$REGION,gamfitemmeans$DATE_NUM))]
fit_sanger_multi_predsbyregion$cases = fit_sanger_multi_predsbyregion$totcases * fit_sanger_multi_predsbyregion$prob
fit_sanger_multi_predsbyregion$cases_rollingmean = fit_sanger_multi_predsbyregion$totcases_rollingmean * fit_sanger_multi_predsbyregion$prob
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
                     limits=c(as.Date("2021-01-01"),max(cases_uk_region$date)+14), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion, 
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


ggplot(data=fit_sanger_multi_predsbyregion, 
       aes(x=collection_date, y=cases_smoothed, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-03-01","2021-06-30")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.pdf"), width=8, height=6)



# some further tests

SGTF_predsbyregion$newCasesBySpecimenDate = cases_uk_region$newCasesBySpecimenDate[match(interaction(SGTF_predsbyregion$collection_date, SGTF_predsbyregion$REGION),
                                                                                         interaction(cases_uk_region$date, cases_uk_region$REGION))]
SGTF_predsbyregion$newCasesBySpecimenDateRollingAvg = cases_uk_region$newCasesBySpecimenDateRollingSum[match(interaction(SGTF_predsbyregion$collection_date, SGTF_predsbyregion$REGION),
                                                                                                   interaction(cases_uk_region$date, cases_uk_region$REGION))]/7
SGTF_predsbyregion$weekday = factor(weekdays(SGTF_predsbyregion$collection_date))
SGTF_predsbyregion$newPCRTestsByPublishDate = NULL
SGTF_predsbyregion3 = SGTF_predsbyregion[complete.cases(SGTF_predsbyregion),]  
SGTF_predsbyregion3 = SGTF_predsbyregion3[SGTF_predsbyregion3$collection_date==max(SGTF_predsbyregion3$collection_date),]

# Calculate Re values from intrinsic growth rate
# BETTER: from epiforecasts package growth_to_R
# https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (better if distribution is known to be gamma, based on the full integral)
Re.from.r <- function(r, gamma_mean=5.5, gamma_sd=1.8) { # from Ferretti et al. 2020
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}
library(mgcv)

qplot(data=SGTF_predsbyregion, x=collection_date, y=newCasesBySpecimenDate, geom="line") +
  facet_wrap(~ REGION) + scale_y_log10()


fit_cases = gam(newCasesBySpecimenDateRollingAvg ~ s(DATE_NUM, bs="cs", k=7, by=REGION) + weekday,
                family=poisson(log), data=SGTF_predsbyregion[SGTF_predsbyregion$collection_date>=as.Date("2021-01-01"),],
)
BIC(fit_cases)
# calculate instantaneous growth rates & 95% CLs using emtrends 
# based on the slope of the GAM fit on a log link scale
extrapolate = 0
avg_r_cases = as.data.frame(emtrends(fit_cases, ~DATE_NUM|REGION, var="DATE_NUM", 
                              at=list(DATE_NUM=seq(min(SGTF_predsbyregion$DATE_NUM),
                                                   max(SGTF_predsbyregion$DATE_NUM)+extrapolate),
                                      REGION=levels(SGTF_predsbyregion$REGION)
                              ),
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
  scale_y_continuous(limits=c(1/3,3), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re VALUES IN ENGLAND BY REGION BASED ON NEW CONFIRMED CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends effect contrasts
extrapolate = 0
emtrsanger_eff = emtrends(fit3_sanger_multi, eff ~ LINEAGE2, by=c("DATE_NUM","REGION"),   
                            var="DATE_NUM",  mode="latent",
                            at=list(DATE_NUM=seq(min(sanger$DATE_NUM),
                                                 max(sanger$DATE_NUM)+extrapolate),
                                    REGION=levels(SGTF_predsbyregion$REGION)
                            ))
above_avg_r_variants = data.frame(confint(emtrsanger_eff, 
                                          adjust="none", df=NA)$contrasts, 
                                  p.value=as.data.frame(emtrsanger_eff$contrasts)$p.value)
above_avg_r_variants$LINEAGE2 = gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$DATE_NUM = above_avg_r_variants$DATE_NUM-7 # the -7 to calculate back to time of infection
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
                DATE_NUM=avg_r_cases$DATE_NUM-7, # -7 to calculate back to time of infection
                REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                p.value=NA,
                collection_date=avg_r_cases$DATE-7,
                LINEAGE2="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$LINEAGE2 = factor(above_avg_r_variants$LINEAGE2, levels=c(levels_LINEAGE2_plot,"avg"))

above_avg_r_variants2 = above_avg_r_variants
ymax = 4
ymin = 1/ymax
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_fill_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  scale_colour_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  theme(legend.position="right",  
        axis.title.x=element_blank())

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant.pdf"), width=8, height=6)


