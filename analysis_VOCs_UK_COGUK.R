# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE UK BASED ON COG-UK SEQUENCING DATA ####

# last update 12 JULY 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-07-12")
today_num = as.numeric(today)
plotdir = "UK_COGUK"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))


coguk = read.csv("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv")
levels_REGION = c("England","Scotland","Wales","Northern Ireland")
coguk$REGION = factor(coguk$adm1, levels=c("UK-ENG","UK-SCT","UK-WLS","UK-NIR"), labels=levels_REGION)
coguk = coguk[!is.na(coguk$REGION),]
levels(coguk$REGION)
coguk$date = as.Date(coguk$sample_date)
head(coguk)

# use data from sept 2020 onwards & Pillar 2 (community as opposed to hospitalised) only
cogukp1 = coguk[coguk$date>=as.Date("2020-09-01")&coguk$is_pillar_2=="N",]
cogukp2 = coguk[coguk$date>=as.Date("2020-09-01")&coguk$is_pillar_2=="Y",]

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
cogukp2 = cogukp2[-which(grepl("B.1.617", cogukp2$lineage, fixed=TRUE)&cogukp2$date<=as.Date("2021-04-14")),]  

nrow(cogukp2) # 
range(cogukp2$date) # "2020-09-01" "2021-06-29"
library(lubridate)
cogukp2$WeekEndDate = floor_date(cogukp2$date,unit="week")+6 # lubridate::week(cogukp2$date)
cogukp2$DATE_NUM = as.numeric(cogukp2$date) 
colnames(cogukp2)

unique(cogukp2$lineage)
sel_target_VOC = "B.1.617"
table(cogukp2$lineage[grepl("B.1.617",cogukp2$lineage, fixed=TRUE)])

unique(cogukp2$lineage[grepl("B.1.617",cogukp2$lineage, fixed=TRUE)])
# "B.1.617.2" "B.1.617.1" "B.1.617.3"
sum(grepl("B.1.617",cogukp2$lineage, fixed=TRUE)) #  B.1.617+
cogukp2$LINEAGE1 = cogukp2$lineage
cogukp2$LINEAGE2 = cogukp2$lineage
cogukp2[grepl(sel_target_VOC, cogukp2$LINEAGE1, fixed=T),"LINEAGE1"] = paste0(sel_target_VOC,"+")
sel_target_VOC = paste0(sel_target_VOC, "+")
cogukp2[grepl("B.1.177", cogukp2$LINEAGE1, fixed=T),"LINEAGE1"] = "B.1.177+"
cogukp2[grepl("B.1.177", cogukp2$LINEAGE2, fixed=T),"LINEAGE2"] = "B.1.177+"

sum("B.1.1.7"==cogukp2$LINEAGE1) #  B.1.1.7


table_lineage = as.data.frame(table(cogukp2$LINEAGE1))
table_lineage$Freq = table_lineage$Freq / sum(table(cogukp2$LINEAGE1))
colnames(table_lineage) = c("lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]

sel_ref_lineage = "B.1.1.7"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.1.7","B.1.177+","B.1.351", "B.1.525","B.1.617+","B.1.617.1","B.1.617.2","P.1")

cogukp2$LINEAGE1[!(cogukp2$LINEAGE1 %in% sel_lineages)] = "other"
cogukp2$LINEAGE2[!(cogukp2$LINEAGE2 %in% sel_lineages)] = "other"
# cogukp2 = cogukp2[cogukp2$LINEAGE1 %in% sel_lineages, ]
sum(table(cogukp2$LINEAGE1))
table(cogukp2$LINEAGE1)
table(cogukp2$LINEAGE2)

levels_LINEAGE1 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617+","P.1")
levels_LINEAGE2 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617.1","B.1.617.2","P.1")
cogukp2$LINEAGE1 = factor(cogukp2$LINEAGE1, levels=levels_LINEAGE1)
cogukp2$LINEAGE2 = factor(cogukp2$LINEAGE2, levels=levels_LINEAGE2)

cogukp2_lastmonth = cogukp2[cogukp2$date>=(max(cogukp2$date)-30),]
table(cogukp2_lastmonth$lineage, cogukp2_lastmonth$REGION)

str(cogukp2)

# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of UK
data_agbyweek1 = as.data.frame(table(cogukp2$WeekEndDate, cogukp2$LINEAGE2))
colnames(data_agbyweek1) = c("WeekEndDate", "LINEAGE2", "count")
data_agbyweek1_sum = aggregate(count ~ WeekEndDate, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$WeekEndDate, data_agbyweek1_sum$WeekEndDate)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE2=="B.1.617.1","total"]) == nrow(cogukp2) # TRUE
data_agbyweek1$WeekEndDate = as.Date(data_agbyweek1$WeekEndDate)
data_agbyweek1$collection_date = data_agbyweek1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and region
data_agbyweekregion1 = as.data.frame(table(cogukp2$WeekEndDate, cogukp2$REGION, cogukp2$LINEAGE2))
colnames(data_agbyweekregion1) = c("WeekEndDate", "REGION", "LINEAGE2", "count")
data_agbyweekregion1_sum = aggregate(count ~ WeekEndDate + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$WeekEndDate,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$WeekEndDate,data_agbyweekregion1_sum$REGION))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE2=="B.1.617.1","total"]) == nrow(cogukp2) # TRUE
data_agbyweekregion1$WeekEndDate = as.Date(data_agbyweekregion1$WeekEndDate) # as.numeric(as.character(data_agbyweekregion1$WeekEndDate))
data_agbyweekregion1$collection_date = data_agbyweekregion1$WeekEndDate -3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$WeekEndDate - 1 ) - 3.5 # we use the week midpoint
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion1[data_agbyweekregion1$collection_date==max(data_agbyweekregion1$collection_date),]

# MULLER PLOT (RAW DATA)
unique(cogukp2$LINEAGE2)
# levels_LINEAGE2_plot = rev(c("B.1.1.7","B.1.617.2","B.1.617.1","P.1","B.1.351","B.1.525","B.1.177+","other")) # "B.1.617.1","B.1.617.2","B.1.617.3"
levels_LINEAGE2_plot = c("other","B.1.177+","B.1.525","B.1.351","P.1","B.1.1.7","B.1.617.1","B.1.617.2")

library(scales)
n1 = length(levels_LINEAGE2_plot)
lineage_cols1 = c(hcl(h = seq(0, 260, length = n1-1), l = 60, c = 200), muted("dodgerblue", l=55, c=120)) # grey70
# lineage_cols1 = c(hcl(h = seq(0, 265, length = n1), l = 60, c = 200)) # grey70
# lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617+")] = "magenta" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.1")] = muted("magenta") # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.1.7")] = "#0085FF"

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
muller_cogukp2_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(COG-UK data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_cogukp2_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plot avg_raw data.png"), width=8, height=6)


muller_cogukp2byregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(COG-UK data)") +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_cogukp2byregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
data_agbyweekregion1$LINEAGE2 = relevel(data_agbyweekregion1$LINEAGE2, ref="B.1.1.7") # we take B.1.1.7 as baseline / reference level
data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
set.seed(1)
fit1_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION + DATE_NUM, weights=count, data=data_agbyweekregion1, maxit=1000)
fit2_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION * DATE_NUM, weights=count, data=data_agbyweekregion1, maxit=1000)
fit3_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION + ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion1, maxit=1000)
fit4_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion1, maxit=1000)
BIC(fit1_cogukp2_multi, fit2_cogukp2_multi, fit3_cogukp2_multi, fit4_cogukp2_multi) 
# fit4_cogukp2_multi fits best (lowest BIC)

# growth rate advantages of different VOCs compared to UK type B.1.1.7 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2 = emtrends(fit4_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2 = data.frame(confint(emtrcogukp2, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrcogukp2$contrasts)$p.value)
delta_r_cogukp2
# contrast    estimate          SE df   asymp.LCL    asymp.UCL      p.value
# 1      other - B.1.1.7  0.02861600 0.002937140 NA  0.02285932  0.034372693 1.356497e-10
# 2 (B.1.177+) - B.1.1.7 -0.04206022 0.007281690 NA -0.05633207 -0.027788370 8.715882e-07
# 3    B.1.525 - B.1.1.7 -0.03369598 0.019907178 NA -0.07271333  0.005321373 3.846499e-01
# 4    B.1.351 - B.1.1.7  0.16029831 0.005336575 NA  0.14983882  0.170757809 1.356012e-10
# 5        P.1 - B.1.1.7  0.07724817 0.010624484 NA  0.05642456  0.098071775 1.329832e-09
# 6  B.1.617.1 - B.1.1.7 -0.27783485 0.012413103 NA -0.30216408 -0.253505614 1.356012e-10
# 7  B.1.617.2 - B.1.1.7  0.10303213 0.003709768 NA  0.09576112  0.110303146 1.356012e-10

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_cogukp2 =  sign(delta_r_cogukp2[,c(2,5,6)])*100*(exp(abs(delta_r_cogukp2[,c(2,5,6)])*4.7)-1)
transmadv_cogukp2 =  data.frame(contrast=delta_r_cogukp2$contrast, transmadv_cogukp2)
transmadv_cogukp2
# contrast   estimate  asymp.LCL   asymp.UCL
# 1      other - B.1.1.7   14.39592   11.34227   17.533317
# 2 (B.1.177+) - B.1.1.7  -21.85761  -30.31191  -13.951797
# 3    B.1.525 - B.1.1.7  -17.16009  -40.74121    2.532584
# 4    B.1.351 - B.1.1.7  112.42145  102.23141  123.124941
# 5        P.1 - B.1.1.7   43.77313   30.36858   58.555950
# 6  B.1.617.1 - B.1.1.7 -269.07282 -313.78288 -229.193779
# 7  B.1.617.2 - B.1.1.7   62.29590   56.84334   67.938018
# so this would estimate that B.1.617.2 had a 63% transmission advantage over B.1.1.7 [57%-68%] 95% CLs


# growth rate & transmission advantages of different VOCs compared to UK type B.1.1.7 by REGION (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2_region = emtrends(fit4_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  by=c("REGION"), mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2_region = data.frame(confint(emtrcogukp2_region, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrcogukp2_region$contrasts)$p.value)
delta_r_cogukp2_region
# contrast           REGION     estimate           SE df    asymp.LCL    asymp.UCL      p.value
# 1       other - B.1.1.7          England  0.042557602 0.0006423038 NA  0.041298710  0.043816495 1.356012e-10
# 2  (B.1.177+) - B.1.1.7          England -0.039894104 0.0016446881 NA -0.043117634 -0.036670575 1.356012e-10
# 3     B.1.525 - B.1.1.7          England -0.018747630 0.0045826450 NA -0.027729449 -0.009765811 6.543315e-04
# 4     B.1.351 - B.1.1.7          England  0.011521909 0.0025846888 NA  0.006456012  0.016587806 1.721620e-04
# 5         P.1 - B.1.1.7          England  0.009895817 0.0069124131 NA -0.003652264  0.023443897 5.457347e-01
# 6   B.1.617.1 - B.1.1.7          England -0.219928131 0.0070621926 NA -0.233769774 -0.206086487 1.356012e-10
# 7   B.1.617.2 - B.1.1.7          England  0.108271116 0.0015722066 NA  0.105189648  0.111352585 1.356012e-10
# 8       other - B.1.1.7         Scotland  0.016753216 0.0037408009 NA  0.009421381  0.024085051 1.593056e-04
# 9  (B.1.177+) - B.1.1.7         Scotland -0.049161655 0.0059805839 NA -0.060883384 -0.037439926 1.509911e-10
# 10    B.1.525 - B.1.1.7         Scotland -0.073617404 0.0265452281 NA -0.125645095 -0.021589712 3.978603e-02
# 11    B.1.351 - B.1.1.7         Scotland  0.031428226 0.0073515796 NA  0.017019395  0.045837057 3.377928e-04
# 12        P.1 - B.1.1.7         Scotland -0.007446851 0.0144350098 NA -0.035738950  0.020845249 9.715500e-01
# 13  B.1.617.1 - B.1.1.7         Scotland -0.298995180 0.0290586335 NA -0.355949055 -0.242041305 1.356316e-10
# 14  B.1.617.2 - B.1.1.7         Scotland  0.086613896 0.0033022365 NA  0.080141631  0.093086160 1.356012e-10
# 15      other - B.1.1.7            Wales  0.019028814 0.0037400991 NA  0.011698354  0.026359273 1.496533e-05
# 16 (B.1.177+) - B.1.1.7            Wales -0.009905621 0.0042707131 NA -0.018276065 -0.001535177 1.193962e-01
# 17    B.1.525 - B.1.1.7            Wales -0.044595731 0.0477544143 NA -0.138192663  0.049001201 8.379491e-01
# 18    B.1.351 - B.1.1.7            Wales  0.016569634 0.0193622399 NA -0.021379659  0.054518927 8.726191e-01
# 19        P.1 - B.1.1.7            Wales -0.021287443 0.0317091772 NA -0.083436289  0.040861402 9.368184e-01
# 20  B.1.617.1 - B.1.1.7            Wales -0.205610366 0.0377622050 NA -0.279622928 -0.131597805 3.492597e-06
# 21  B.1.617.2 - B.1.1.7            Wales  0.146982526 0.0062483002 NA  0.134736082  0.159228969 1.356012e-10
# 22      other - B.1.1.7 Northern Ireland  0.036124388 0.0104706006 NA  0.015602388  0.056646388 5.605878e-03
# 23 (B.1.177+) - B.1.1.7 Northern Ireland -0.069279502 0.0281363863 NA -0.124425805 -0.014133198 8.625913e-02
# 24    B.1.525 - B.1.1.7 Northern Ireland  0.002176846 0.0576596215 NA -0.110833935  0.115187628 9.999940e-01
# 25    B.1.351 - B.1.1.7 Northern Ireland  0.581673485 0.0025846888 NA  0.576607588  0.586739382 1.356012e-10
# 26        P.1 - B.1.1.7 Northern Ireland  0.327831157 0.0069124131 NA  0.314283076  0.341379237 1.356012e-10
# 27  B.1.617.1 - B.1.1.7 Northern Ireland -0.386805719 0.0070621959 NA -0.400647368 -0.372964069 1.356012e-10
# 28  B.1.617.2 - B.1.1.7 Northern Ireland  0.070260996 0.0130484516 NA  0.044686501  0.095835491 4.478597e-06


# Just for B.1.617.2 / delta
delta_r_cogukp2_region[delta_r_cogukp2_region$contrast=="B.1.617.2 - B.1.1.7",]
# contrast           REGION  estimate          SE df  asymp.LCL  asymp.UCL      p.value
# 7  B.1.617.2 - B.1.1.7          England 0.1082711 0.001572207 NA 0.10518965 0.11135258 1.356012e-10
# 14 B.1.617.2 - B.1.1.7         Scotland 0.0866139 0.003302236 NA 0.08014163 0.09308616 1.356012e-10
# 21 B.1.617.2 - B.1.1.7            Wales 0.1469825 0.006248300 NA 0.13473608 0.15922897 1.356012e-10
# 28 B.1.617.2 - B.1.1.7 Northern Ireland 0.0702610 0.013048452 NA 0.04468650 0.09583549 4.478597e-06



# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME BY REGION ####

gentime = 4.7 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantages of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit4_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=cogukp2, maxit=1000) here
emtrcogukp24 = emtrends(fit4_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2, by=c("DATE_NUM", "REGION"), 
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=seq(as.numeric(as.Date("2021-05-01")),as.numeric(as.Date("2021-06-14")))))
delta_r_cogukp24 = data.frame(confint(emtrcogukp24, 
                                    adjust="none", df=NA)$contrasts)
delta_r_cogukp24 = delta_r_cogukp24[delta_r_cogukp24$contrast=="B.1.617.2 - B.1.1.7",]
delta_r_cogukp24

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_cogukp24 = delta_r_cogukp24
transmadv_cogukp24$estimate =  100*(exp(transmadv_cogukp24$estimate*4.7)-1)
transmadv_cogukp24$asymp.LCL =  100*(exp(transmadv_cogukp24$asymp.LCL*4.7)-1)
transmadv_cogukp24$asymp.UCL =  100*(exp(transmadv_cogukp24$asymp.UCL*4.7)-1)
transmadv_cogukp24$collection_date = as.Date(transmadv_cogukp24$DATE_NUM, origin="1970-01-01")
transmadv_cogukp24$REGION = factor(transmadv_cogukp24$REGION, levels=levels_REGION)

plot_cogukp2_mfit4_transmadv = qplot(data=transmadv_cogukp24, 
                               x=collection_date, y=estimate, ymin=asymp.LCL, ymax=asymp.UCL, colour=REGION, fill=REGION, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(fill=REGION, colour=NULL), alpha=I(0.3)) +
  geom_line(aes(colour=REGION, fill=NULL)) +
  ylab("Transmission advantage of B.1.617.2 over B.1.1.7 (%)") +
  theme_hc() + xlab("") +
  ggtitle("TRANSMISSION ADVANTAGE OF B.1.617.2 OVER B.1.1.7 BY REGION\n(COG-UK data, multinomial fit)") +
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
plot_cogukp2_mfit4_transmadv

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit4_transm advantage.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit4_transm advantage.pdf"), width=8, height=6)




# PLOT MULTINOMIAL FIT ####

# extrapolate = 60
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-07-31")) # max(cogukp2$DATE_NUM)+extrapolate

fit_cogukp2_multi_predsbyregion = data.frame(emmeans(fit4_cogukp2_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "REGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=1)), 
                                                  mode="prob", df=NA))
fit_cogukp2_multi_predsbyregion$collection_date = as.Date(fit_cogukp2_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_cogukp2_multi_predsbyregion$LINEAGE2 = factor(fit_cogukp2_multi_predsbyregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_cogukp2_multi_predsbyregion$REGION = factor(fit_cogukp2_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of the UK today
fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$collection_date==today-1,]


fit_cogukp2_multi_preds = data.frame(emmeans(fit4_cogukp2_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=3)), 
                                           mode="prob", df=NA))
fit_cogukp2_multi_preds$collection_date = as.Date(fit_cogukp2_multi_preds$DATE_NUM, origin="1970-01-01")
fit_cogukp2_multi_preds$LINEAGE2 = factor(fit_cogukp2_multi_preds$LINEAGE2, levels=levels_LINEAGE2_plot) 

muller_cogukp2_mfit = ggplot(data=fit_cogukp2_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(cogukp2$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(COG-UK data)")
muller_cogukp2_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plot avg_multinom fit.png"), width=8, height=6)


library(ggpubr)
ggarrange(muller_cogukp2_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_cogukp2_mfit+ggtitle("Multinomial fit"), ncol=1)


muller_cogukp2byregion_mfit = ggplot(data=fit_cogukp2_multi_predsbyregion, 
                                  aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(cogukp2$DATE_NUM)+1, 
           xmax=as.Date(c("2021-07-31")), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit)")
muller_cogukp2byregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit.pdf"), width=8, height=6)


ggarrange(muller_cogukp2byregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_cogukp2byregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.png"), width=8, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & 95% CONFIDENCE INTERVALS


# on logit scale:

fit_cogukp2_multi_preds2 = fit_cogukp2_multi_preds
ymin = 0.001
ymax = 0.998
fit_cogukp2_multi_preds2$asymp.LCL[fit_cogukp2_multi_preds2$asymp.LCL<ymin] = ymin
fit_cogukp2_multi_preds2$asymp.UCL[fit_cogukp2_multi_preds2$asymp.UCL<ymin] = ymin
fit_cogukp2_multi_preds2$asymp.UCL[fit_cogukp2_multi_preds2$asymp.UCL>ymax] = ymax
fit_cogukp2_multi_preds2$prob[fit_cogukp2_multi_preds2$prob<ymin] = ymin

fit_cogukp2_multi_predsbyregion2 = fit_cogukp2_multi_predsbyregion
fit_cogukp2_multi_predsbyregion2$asymp.LCL[fit_cogukp2_multi_predsbyregion2$asymp.LCL<ymin] = ymin
fit_cogukp2_multi_predsbyregion2$asymp.UCL[fit_cogukp2_multi_predsbyregion2$asymp.UCL<ymin] = ymin
fit_cogukp2_multi_predsbyregion2$asymp.UCL[fit_cogukp2_multi_predsbyregion2$asymp.UCL>ymax] = ymax
fit_cogukp2_multi_predsbyregion2$prob[fit_cogukp2_multi_predsbyregion2$prob<ymin] = ymin

# plots of multinomial fit to cogukp2 Inst data 

plot_cogukp2_mfit_logit = qplot(data=fit_cogukp2_multi_predsbyregion2, 
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
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
                        range=c(1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-14")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_cogukp2_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.pdf"), width=8, height=6)
library(svglite)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.svg"), width=8, height=6)

# on response scale:
plot_cogukp2_mfit = qplot(data=fit_cogukp2_multi_predsbyregion2, 
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
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
                        range=c(1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_cogukp2_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_response scale.pdf"), width=8, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot = as.data.frame(get_regional_data(country = "United Kingdom", regions=c("England","Scotland","Wales","Northern Ireland")))
cases_tot = cases_tot[cases_tot$date>=as.Date("2020-08-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot$region = factor(cases_tot$region, levels=c("England","Scotland","Wales","Northern Ireland"))
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=20, fx=F, by=region) + region +
                  WEEKDAY # +
                  # s(tested_new, bs="cs", k=8, fx=F, by=region)
                ,
                family=poisson(log), data=cases_tot,
) 
BIC(fit_cases)

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_cogukp2_multi_predsbyregion$totcases = cases_tot$cases_new[match(interaction(fit_cogukp2_multi_predsbyregion$DATE_NUM,fit_cogukp2_multi_predsbyregion$REGION),
                                                                     interaction(cases_tot$DATE_NUM,cases_tot$region))]
fit_cogukp2_multi_predsbyregion$cases = fit_cogukp2_multi_predsbyregion$totcases * fit_cogukp2_multi_predsbyregion$prob
fit_cogukp2_multi_predsbyregion$cases[fit_cogukp2_multi_predsbyregion$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM|region, at=list(DATE_NUM=seq(date.from, date.to, by=0.5)
                                                                     ), type="response"))
fit_cogukp2_multi_predsbyregion$smoothed_totcases = cases_emmeans$rate[match(interaction(fit_cogukp2_multi_predsbyregion$DATE_NUM,fit_cogukp2_multi_predsbyregion$REGION),
                                                                             interaction(cases_emmeans$DATE_NUM,cases_emmeans$region))]
fit_cogukp2_multi_predsbyregion$smoothed_cases = fit_cogukp2_multi_predsbyregion$smoothed_totcases * fit_cogukp2_multi_predsbyregion$prob
fit_cogukp2_multi_predsbyregion$smoothed_cases[fit_cogukp2_multi_predsbyregion$smoothed_cases<=0.001] = NA

ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=cases+1, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_cogukp2_multi_predsbyregion,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) +
  facet_wrap(~ REGION, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=smoothed_cases+1, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE2=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE2=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=cases, group=LINEAGE2)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_cogukp2_multi_predsbyregion,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE2)) +
  facet_wrap(~ REGION, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),today), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
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
  above_avg_r_variants$prob = fit_cogukp2_multi_predsbyregion$prob[match(interaction(round(above_avg_r_variants$DATE_NUM),
                                                                        as.character(above_avg_r_variants$LINEAGE2),
                                                                        as.character(above_avg_r_variants$REGION)),
                                                            interaction(round(fit_cogukp2_multi_predsbyregion$DATE_NUM),
                                                                        as.character(fit_cogukp2_multi_predsbyregion$LINEAGE2),
                                                                        as.character(fit_cogukp2_multi_predsbyregion$REGION)))]
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


