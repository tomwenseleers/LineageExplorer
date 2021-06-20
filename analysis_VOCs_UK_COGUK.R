# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE UK BASED ON COG-UK SEQUENCING DATA ####

# last update 19 JUNE 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-20")
today_num = as.numeric(today)
today # "2021-06-20"
plotdir = "VOCs_COGUK"
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
range(cogukp2$date) # "2020-09-01" "2021-05-30"
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
# lineage        Prop
# 74   B.1.1.7 0.882927846
# 85  B.1.177+ 0.027505266
# 111  B.1.351 0.002023195
# 160  B.1.525 0.001810227
# 175 B.1.617+ 0.072756314

sel_ref_lineage = "B.1.1.7"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.1.7","B.1.177+","B.1.351", "B.1.525","B.1.617+","B.1.617.1","B.1.617.2","P.1")

cogukp2$LINEAGE1[!(cogukp2$LINEAGE1 %in% sel_lineages)] = "other"
cogukp2$LINEAGE2[!(cogukp2$LINEAGE2 %in% sel_lineages)] = "other"
# cogukp2 = cogukp2[cogukp2$LINEAGE1 %in% sel_lineages, ]
sum(table(cogukp2$LINEAGE1))
table(cogukp2$LINEAGE1)
# B.1.1.7 B.1.177+  B.1.351  B.1.525 B.1.617+    other      P.1 
# 190708     5941      437      391    15715     2686      117 
table(cogukp2$LINEAGE2)
# B.1.1.7  B.1.177+   B.1.351   B.1.525 B.1.617.1 B.1.617.2     other       P.1 
# 190708      5941       437       391       345     15363      2693       117 


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
muller_cogukp2_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
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


muller_cogukp2byregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
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
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_raw data.pdf"), width=8, height=6)



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
# fit4_cogukp2_multi fits best (lowest BIC) but I will use the slightly simpler fit3_cogukp2_multi

# growth rate advantages of different VOCs compared to UK type B.1.1.7 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2 = emtrends(fit3_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2 = data.frame(confint(emtrcogukp2, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrcogukp2$contrasts)$p.value)
delta_r_cogukp2
#               contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1      other - B.1.1.7  0.031940742 0.0006068275 NA  0.030751382  0.033130102 9.059420e-13
# 2 (B.1.177+) - B.1.1.7 -0.041368883 0.0012839060 NA -0.043885292 -0.038852473 9.059420e-13
# 3    B.1.525 - B.1.1.7 -0.016379003 0.0039642947 NA -0.024148878 -0.008609128 1.092326e-03
# 4    B.1.351 - B.1.1.7  0.012243393 0.0023067364 NA  0.007722273  0.016764514 2.641193e-05
# 5        P.1 - B.1.1.7  0.009996597 0.0067542034 NA -0.003241398  0.023234593 5.195633e-01
# 6  B.1.617.1 - B.1.1.7 -0.197641943 0.0068245268 NA -0.211017770 -0.184266116 9.059420e-13
# 7  B.1.617.2 - B.1.1.7  0.100082204 0.0015455350 NA  0.097053011  0.103111397 9.059420e-13

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_cogukp2 =  sign(delta_r_cogukp2[,c(2,5,6)])*100*(exp(abs(delta_r_cogukp2[,c(2,5,6)])*4.7)-1)
transmadv_cogukp2 =  data.frame(contrast=delta_r_cogukp2$contrast, transmadv_cogukp2)
transmadv_cogukp2
#            contrast       estimate      asymp.LCL      asymp.UCL
# 1      other - B.1.1.7   16.197540   15.549809   16.848902
# 2 (B.1.177+) - B.1.1.7  -21.462301  -22.907380  -20.034213
# 3    B.1.525 - B.1.1.7   -8.002190  -12.019158   -4.129268
# 4    B.1.351 - B.1.1.7    5.923182    3.696138    8.198056
# 5        P.1 - B.1.1.7    4.810525   -1.535121   11.538829
# 6  B.1.617.1 - B.1.1.7 -153.176613 -169.603862 -137.750293
# 7  B.1.617.2 - B.1.1.7   60.061248   57.798568   62.356373
# so this would estimate that B.1.617.2 had a 60% transmission advantage over B.1.1.7 [58%-62%] 95% CLs


# growth rate & transmission advantages of different VOCs compared to UK type B.1.1.7 by REGION (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2_region = emtrends(fit4_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  by=c("REGION"), mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2_region = data.frame(confint(emtrcogukp2_region, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrcogukp2_region$contrasts)$p.value)
delta_r_cogukp2_region
# contrast           REGION     estimate           SE df    asymp.LCL     asymp.UCL      p.value
# 1       other - B.1.1.7          England  0.035701674 0.0006294629 NA  0.034467950  0.0369353989 1.356012e-10
# 2  (B.1.177+) - B.1.1.7          England -0.039554582 0.0014736813 NA -0.042442944 -0.0366662200 1.356012e-10
# 3     B.1.525 - B.1.1.7          England -0.015098729 0.0041862349 NA -0.023303599 -0.0068938594 3.396157e-03
# 4     B.1.351 - B.1.1.7          England  0.010783332 0.0024685859 NA  0.005944993  0.0156216716 2.400369e-04
# 5         P.1 - B.1.1.7          England  0.017275700 0.0064707392 NA  0.004593285  0.0299581160 5.196613e-02
# 6   B.1.617.1 - B.1.1.7          England -0.114606254 0.0057083749 NA -0.125794463 -0.1034180446 1.356012e-10
# 7   B.1.617.2 - B.1.1.7          England  0.104758127 0.0016580507 NA  0.101508408  0.1080078468 1.356012e-10
# 8       other - B.1.1.7         Scotland  0.009298424 0.0037240633 NA  0.001999394  0.0165974543 7.949595e-02
# 9  (B.1.177+) - B.1.1.7         Scotland -0.050008216 0.0054455265 NA -0.060681252 -0.0393351801 1.358152e-10
# 10    B.1.525 - B.1.1.7         Scotland -0.034548389 0.0171916285 NA -0.068243361 -0.0008534161 2.249197e-01
# 11    B.1.351 - B.1.1.7         Scotland  0.031164779 0.0069938165 NA  0.017457151  0.0448724075 1.732631e-04
# 12        P.1 - B.1.1.7         Scotland -0.023606593 0.0148298135 NA -0.052672493  0.0054593074 4.447900e-01
# 13  B.1.617.1 - B.1.1.7         Scotland -0.216980546 0.0256971422 NA -0.267346019 -0.1666150730 1.410992e-10
# 14  B.1.617.2 - B.1.1.7         Scotland  0.057942643 0.0022723291 NA  0.053488960  0.0623963262 1.356012e-10
# 15      other - B.1.1.7            Wales  0.014736014 0.0037609376 NA  0.007364712  0.0221073164 1.196866e-03
# 16 (B.1.177+) - B.1.1.7            Wales -0.011657629 0.0039839956 NA -0.019466117 -0.0038491407 2.637231e-02
# 17    B.1.525 - B.1.1.7            Wales -0.084497661 0.0295647086 NA -0.142443425 -0.0265518969 3.174798e-02
# 18    B.1.351 - B.1.1.7            Wales  0.021652179 0.0189151393 NA -0.015420813  0.0587251703 7.253400e-01
# 19        P.1 - B.1.1.7            Wales  0.006267008 0.0295142182 NA -0.051579797  0.0641138124 9.983100e-01
# 20  B.1.617.1 - B.1.1.7            Wales -0.143468299 0.0337353750 NA -0.209588419 -0.0773481789 3.662905e-04
# 21  B.1.617.2 - B.1.1.7            Wales  0.120541627 0.0067223368 NA  0.107366089  0.1337171648 1.356012e-10
# 22      other - B.1.1.7 Northern Ireland  0.035800784 0.0116205653 NA  0.013024895  0.0585766735 1.705995e-02
# 23 (B.1.177+) - B.1.1.7 Northern Ireland -0.048374856 0.0286926900 NA -0.104611495  0.0078617834 3.885274e-01
# 24    B.1.525 - B.1.1.7 Northern Ireland  0.051775334 0.0634328097 NA -0.072550689  0.1761013560 8.885164e-01
# 25    B.1.351 - B.1.1.7 Northern Ireland -0.113453975 0.0024686168 NA -0.118292376 -0.1086155754 1.356012e-10
# 26        P.1 - B.1.1.7 Northern Ireland -0.401972717 0.0064704659 NA -0.414654597 -0.3892908365 1.356012e-10
# 27  B.1.617.1 - B.1.1.7 Northern Ireland -0.790845103 0.0057083749 NA -0.802033312 -0.7796568937 1.356012e-10
# 28  B.1.617.2 - B.1.1.7 Northern Ireland  0.034558836 0.0149964332 NA  0.005166367  0.0639513048 1.233872e-01







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
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit4_transm advantage.pdf"), width=8, height=6)




# PLOT MULTINOMIAL FIT ####

# extrapolate = 60
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-07-31")) # max(cogukp2$DATE_NUM)+extrapolate

fit_cogukp2_multi_predsbyregion = data.frame(emmeans(fit3_cogukp2_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "REGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=3)), 
                                                  mode="prob", df=NA))
fit_cogukp2_multi_predsbyregion$collection_date = as.Date(fit_cogukp2_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_cogukp2_multi_predsbyregion$LINEAGE2 = factor(fit_cogukp2_multi_predsbyregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_cogukp2_multi_predsbyregion$REGION = factor(fit_cogukp2_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of the UK today
fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$collection_date==today-1,]


fit_cogukp2_multi_preds = data.frame(emmeans(fit3_cogukp2_multi, 
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
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(COG-UK data)")
muller_cogukp2_mfit

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
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit)")
muller_cogukp2byregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit.pdf"), width=8, height=6)


ggarrange(muller_cogukp2byregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_cogukp2byregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.pdf"), width=8, height=6)



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
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
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
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.pdf"), width=8, height=6)
library(svglite)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.svg"), width=8, height=6)

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
  scale_x_continuous(breaks=as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
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
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_response scale.pdf"), width=8, height=6)



