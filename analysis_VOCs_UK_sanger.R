# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN UK BASED ON SANGER INSTITUTE BASELINE SURVEILLANCE SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & surge testing/active surveillance)

# last update 10 MAY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-05-10")
today_num = as.numeric(today)
today # "2021-05-10"
plotdir = "VOCs_SANGER"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))


sanger = read_tsv("https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv")
sanger = as.data.frame(sanger)
LTLAs_regions = read.csv(".//data/ONS/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")
sanger$REGION = LTLAs_regions$RGN20NM[match(sanger$LTLA, LTLAs_regions$LAD20CD)]
sanger = sanger[!is.na(sanger$REGION),]
sanger$REGION = factor(sanger$REGION)
levels(sanger$REGION)
levels_REGION = c("London","North West","East Midlands","East of England","South East","South West","West Midlands",
                  "North East","Yorkshire and The Humber")
sanger$REGION = factor(sanger$REGION, levels=levels_REGION)
head(sanger)

sanger = sanger[grepl("2021-", sanger[,"WeekEndDate"]),]
sanger = sanger[sanger$Lineage!="None",]
nrow(sanger) # 8465
range(sanger$WeekEndDate) # "2021-01-02" "2021-05-01"
sanger$Week = lubridate::week(sanger$WeekEndDate)
sanger$DATE_NUM = as.numeric(sanger$WeekEndDate)-3.5 # using week midpoint
colnames(sanger)

sanger = sanger[rep(seq_len(nrow(sanger)), sanger$Count),] # convert to long format
sanger$Count = NULL
nrow(sanger) # 137 820

nrow(sanger[sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 11905 (last 2 weeks)
nrow(sanger[grepl("B.1.617",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 562 (last 2 weeks)
nrow(sanger[grepl("B.1.1.7",sanger$Lineage, fixed=TRUE)&sanger$WeekEndDate>=(max(sanger$WeekEndDate)-14),]) # 11 138

unique(sanger$Lineage)
sel_target_VOC = "B.1.617"
table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
sanger[sanger$Lineage=="B.1.617",]
sum(table(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])) # 594 B.1.617

unique(sanger$Lineage[grepl("B.1.617",sanger$Lineage, fixed=TRUE)])
# "B.1.617"   "B.1.617.1" "B.1.617.3" "B.1.617.2"
sum(grepl("B.1.617",sanger$Lineage, fixed=TRUE)) # 594 B.1.617+
sanger$LINEAGE1 = sanger$Lineage
sanger$LINEAGE2 = sanger$Lineage
sanger[grepl(sel_target_VOC, sanger$LINEAGE1, fixed=T),"LINEAGE1"] = paste0(sel_target_VOC,"+")
sel_target_VOC = paste0(sel_target_VOC, "+")
sanger[grepl("B.1.177", sanger$LINEAGE1, fixed=T),"LINEAGE1"] = "B.1.177+"
sanger[grepl("B.1.177", sanger$LINEAGE2, fixed=T),"LINEAGE2"] = "B.1.177+"

sum("B.1.1.7"==sanger$LINEAGE1) # 130015 B.1.1.7


table_lineage = as.data.frame(table(sanger$LINEAGE1))
table_lineage$Freq = table_lineage$Freq / sum(table(sanger$LINEAGE1))
colnames(table_lineage) = c("Lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# Lineage        Prop
# 68   B.1.1.7 0.943368161
# 77  B.1.177+ 0.039820055
# 88   B.1.258 0.001044841
# 97   B.1.351 0.001545494
# 132 B.1.617+ 0.004309970

sel_ref_lineage = "B.1.1.7"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"Lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.1.7","B.1.177+","B.1.351", "B.1.525","B.1.617+","B.1.617.1","B.1.617.2","P.1")

sanger$LINEAGE1[!(sanger$LINEAGE1 %in% sel_lineages)] = "other"
sanger$LINEAGE2[!(sanger$LINEAGE2 %in% sel_lineages)] = "other"
sanger = sanger[sanger$LINEAGE1 %in% sel_lineages, ]
sum(table(sanger$LINEAGE1))
table(sanger$LINEAGE1)
# B.1.1.7 B.1.177+  B.1.351  B.1.525 B.1.617+      P.1 
# 130015     5488      213      103      594       34 
table(sanger$LINEAGE2)
# B.1.1.7  B.1.177+   B.1.351   B.1.525 B.1.617.1 B.1.617.2     other       P.1 
# 130015      5488       213       103        90       495         9        34

levels_LINEAGE1 = c("B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617+","P.1")
levels_LINEAGE2 = c("B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617.1","B.1.617.2","P.1")
sanger$LINEAGE1 = factor(sanger$LINEAGE1, levels=levels_LINEAGE1)
sanger$LINEAGE2 = factor(sanger$LINEAGE2, levels=levels_LINEAGE2)


str(sanger)

# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of UK
data_agbyweek1 = as.data.frame(table(sanger$Week, sanger$LINEAGE2))
colnames(data_agbyweek1) = c("Week", "LINEAGE2", "count")
data_agbyweek1_sum = aggregate(count ~ Week, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$Week, data_agbyweek1_sum$Week)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE2=="B.1.617.1","total"]) == nrow(sanger) # should evaluate to TRUE
data_agbyweek1$Week = as.numeric(as.character(data_agbyweek1$Week))
data_agbyweek1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) + 3.5
data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and region
data_agbyweekregion1 = as.data.frame(table(sanger$Week, sanger$REGION, sanger$LINEAGE2))
colnames(data_agbyweekregion1) = c("Week", "REGION", "LINEAGE2", "count")
data_agbyweekregion1_sum = aggregate(count ~ Week + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$Week,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$Week,data_agbyweekregion1_sum$REGION))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE2=="B.1.617.1","total"]) == nrow(sanger) # should evaluate to TRUE
data_agbyweekregion1$Week = as.numeric(as.character(data_agbyweekregion1$Week))
data_agbyweekregion1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) + 3.5
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]


# MULLER PLOT (RAW DATA)
unique(sanger$LINEAGE2)
levels_LINEAGE2_plot = rev(c("B.1.1.7","B.1.617.2","B.1.617.1","P.1","B.1.351","B.1.525","B.1.177+")) # "B.1.617.1","B.1.617.2","B.1.617.3"

library(scales)
n1 = length(levels_LINEAGE2_plot)
lineage_cols1 = c(hcl(h = seq(0, 260, length = n1-1), l = 60, c = 200), muted("dodgerblue", l=60, c=120)) # grey70
lineage_cols1 = c(hcl(h = seq(0, 265, length = n1), l = 60, c = 200)) # grey70
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617+")] = "magenta" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.1")] = muted("magenta") # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=50,c=100)

# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.1")] = "magenta" # muted("magenta",l=50,c=100)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=80,c=255)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.3")] = "magenta" # muted("magenta",l=100,c=200)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.1.7")] = "grey60"  
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.177+")] = "grey70"  
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
  ylab("Share among major+VOC lineages") + 
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
  ylab("Share among major+VOC lineages") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(Sanger Institute baseline surveillance data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbyregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
fit1_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION + DATE_NUM, data=sanger, maxit=1000)
fit2_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * DATE_NUM, data=sanger, maxit=1000)
fit3_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION + ns(DATE_NUM, df=2), data=sanger, maxit=1000)
fit4_sanger_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=sanger, maxit=1000)
BIC(fit1_sanger_multi, fit2_sanger_multi, fit3_sanger_multi, fit4_sanger_multi) 
# fit1_sanger_multi fits best (lowest BIC), but I will use the slightly more flexible model fit3_sanger_multi

# growth rate advantages of different VOCs compared to UK type B.1.1.7 (difference in growth rate per day) 
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-04-27
emtrsanger = emtrends(fit3_sanger_multi, trt.vs.ctrl ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger = data.frame(confint(emtrsanger, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger$contrasts)$p.value)
delta_r_sanger
#               contrast    estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1 (B.1.177+) - B.1.1.7 -0.10103244 0.007362119 NA -0.115461927 -0.08660295 0.000000e+00
# 2    B.1.351 - B.1.1.7  0.03687716 0.005601636 NA  0.025898157  0.04785617 5.286363e-08
# 3    B.1.525 - B.1.1.7  0.01022481 0.009134869 NA -0.007679201  0.02812883 7.066696e-01
# 4  B.1.617.1 - B.1.1.7  0.03688804 0.021418094 NA -0.005090650  0.07886674 3.389213e-01
# 5  B.1.617.2 - B.1.1.7  0.15275070 0.012694119 NA  0.127870681  0.17763071 0.000000e+00
# 6        P.1 - B.1.1.7  0.07852189 0.014860068 NA  0.049396689  0.10764709 9.010158e-06

# pairwise growth rate advantages for all strain comparisons (i.e. pairwise differences in growth rate per day among the different lineages)
max(as.Date(sanger$WeekEndDate)-3.5) # 2021-04-47
emtrsanger_pairw = emtrends(fit3_sanger_multi, pairwise ~ LINEAGE2,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_pairw = data.frame(confint(emtrsanger_pairw, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger_pairw$contrasts)$p.value)
delta_r_sanger_pairw
#                  contrast      estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1    B.1.1.7 - (B.1.177+)  0.1009772963 0.007359498 NA  0.086552945  0.115401648 0.000000e+00
# 2       B.1.1.7 - B.1.351 -0.0368669523 0.005600951 NA -0.047844615 -0.025889290 1.838407e-07
# 3       B.1.1.7 - B.1.525 -0.0101598424 0.009138546 NA -0.028071063  0.007751378 9.223125e-01
# 4     B.1.1.7 - B.1.617.1 -0.0364712286 0.021452872 NA -0.078518084  0.005575627 6.184330e-01
# 5     B.1.1.7 - B.1.617.2 -0.1506055170 0.012824419 NA -0.175740915 -0.125470119 0.000000e+00
# 6           B.1.1.7 - P.1 -0.0784581860 0.014855118 NA -0.107573682 -0.049342690 3.090438e-05
# 7    (B.1.177+) - B.1.351 -0.1378442486 0.009246883 NA -0.155967806 -0.119720691 0.000000e+00
# 8    (B.1.177+) - B.1.525 -0.1111371387 0.011729986 NA -0.134127488 -0.088146789 0.000000e+00
# 9  (B.1.177+) - B.1.617.1 -0.1374485249 0.022680309 NA -0.181901114 -0.092995936 1.493176e-06
# 10 (B.1.177+) - B.1.617.2 -0.2515828133 0.014786390 NA -0.280563605 -0.222602021 0.000000e+00
# 11       (B.1.177+) - P.1 -0.1794354822 0.016578393 NA -0.211928536 -0.146942429 0.000000e+00
# 12      B.1.351 - B.1.525  0.0267071098 0.010703757 NA  0.005728131  0.047686089 1.779490e-01
# 13    B.1.351 - B.1.617.1  0.0003957236 0.022140492 NA -0.042998843  0.043790290 1.000000e+00
# 14    B.1.351 - B.1.617.2 -0.1137385647 0.013954876 NA -0.141089618 -0.086387511 2.887882e-10
# 15          B.1.351 - P.1 -0.0415912337 0.015843669 NA -0.072644255 -0.010538213 1.351623e-01
# 16    B.1.525 - B.1.617.1 -0.0263113862 0.023309503 NA -0.071997172  0.019374400 9.168940e-01
# 17    B.1.525 - B.1.617.2 -0.1404456746 0.015735824 NA -0.171287323 -0.109604026 3.722467e-12
# 18          B.1.525 - P.1 -0.0682983435 0.017428126 NA -0.102456842 -0.034139845 3.845576e-03
# 19  B.1.617.1 - B.1.617.2 -0.1141342884 0.024736285 NA -0.162616516 -0.065652060 3.617981e-04
# 20        B.1.617.1 - P.1 -0.0419869573 0.026006083 NA -0.092957943  0.008984028 6.735838e-01
# 21        B.1.617.2 - P.1  0.0721473310 0.019529344 NA  0.033870519  0.110424143 7.800988e-03



# PLOT MULTINOMIAL FIT

# extrapolate = 60
date.from = as.numeric(as.Date("2021-01-01"))
date.to = as.numeric(as.Date("2021-06-01")) # max(sanger$DATE_NUM)+extrapolate

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
  ylab("Share among major+VOC lineages") +
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
           xmax=as.Date(c("2021-05-31")), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-31")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share among major+VOC lineages") +
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





# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

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

plot_sanger_mfit_logit = qplot(data=fit_sanger_multi_predsbyregion2, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share among major+VOC lineages (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-31")), 
                     expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 4), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit_logit scale.pdf"), width=8, height=6)


# on response scale:
plot_sanger_mfit = qplot(data=fit_sanger_multi_predsbyregion2, x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share among major+VOC lineages (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(Sanger Institute baseline surveillance data)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-05-01")), 
                     expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",date.to-1)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE2
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.01, 4), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_sanger_mfit




