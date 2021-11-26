# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE USA (GISAID GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 10 OCT 2021

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
today = as.Date("2021-10-10")
today_num = as.numeric(today)
plotdir = "USA_GISAID_genomic_epidemiology"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID genomic epidemiology metadata
GISAID = read_tsv(gzfile(".//data//GISAID_genomic_epidemiology//metadata_2021-10-01_13-45.tsv.gz"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(GISAID)

GISAID$date = as.Date(GISAID$date)
GISAID = GISAID[!is.na(GISAID$date),]
GISAID[GISAID$host!="Human","strain"]
GISAID = GISAID[GISAID$host=="Human",]
GISAID = GISAID[GISAID$date>=as.Date("2020-01-01"),]
range(GISAID$date) # "2020-01-01" "2021-10-04"
nrow(GISAID) #  3851990
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
library(lubridate)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
colnames(GISAID)

length(unique(GISAID$country[grepl("AY",GISAID$pango_lineage,fixed=T)|grepl("B.1.617",GISAID$pango_lineage,fixed=T)])) # B.1.617+ now found in 103 countries
table(GISAID$pango_lineage[grepl("AY",GISAID$pango_lineage,fixed=T)|grepl("B.1.617",GISAID$pango_lineage,fixed=T)])

GISAID$pango_lineage[grepl("B.1.177",GISAID$pango_lineage,fixed=T)] = "B.1.177+"
GISAID$pango_lineage[grepl("B.1.36\\>",GISAID$pango_lineage)] = "B.1.36+"
# GISAID$pango_lineage[grepl("B.1.617.2|AY",GISAID$pango_lineage)] = "Delta (B.1.617.2 & AY.X)"

# sel_target_VOC = "Delta (B.1.617.2 & AY.X)"
sel_target_VOC = "B.1.617.2"

GISAID$LINEAGE2 = GISAID$pango_lineage


# ANALYSIS VOCs IN THE USA ####
sel_countries = "USA"
GISAID_sel = GISAID[GISAID$country %in% sel_countries,]
GISAID_sel = GISAID_sel[GISAID_sel$date>=as.Date("2021-01-01"),]
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$LINEAGE2),]
GISAID_sel = GISAID_sel[!GISAID_sel$LINEAGE2=="None",]
GISAID_sel = GISAID_sel[GISAID_sel$country==GISAID_sel$country_exposure,] # we remove travel-related cases
range(GISAID_sel$date) # "2021-01-01" "2021-09-24"
nrow(GISAID_sel) # 967414

main_lineages = names(table(GISAID_sel$LINEAGE2))[100*table(GISAID_sel$LINEAGE2)/sum(table(GISAID_sel$LINEAGE2)) > 3]
main_lineages
# "AY.25"     "AY.3"      "AY.4"      "B.1.1.7"   "B.1.2"     "B.1.429"   "B.1.526"   "B.1.617.2"
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526","B.1.1.519","AY.1","AY.2","AY.3","AY.4","AY.4","AY.25","AY.26","AY.16","AY.34","B.1.621")
main_lineages = union(main_lineages, VOCs)
GISAID_sel$LINEAGE2[!(GISAID_sel$LINEAGE2 %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove2 = names(table(GISAID_sel$LINEAGE2))[table(GISAID_sel$LINEAGE2)/sum(table(GISAID_sel$LINEAGE2)) < 0.01]
remove2 = remove2[!(remove2 %in% c("B.1.351","B.1.1.7","P.1","B.1.617.2","B.1.1.519","AY.3","B.1.621","AY.16","B.1.621"))]
GISAID_sel$LINEAGE2[(GISAID_sel$LINEAGE2 %in% remove2)] = "other" # minority VOCs
table(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2)
GISAID_sel$LINEAGE2 = relevel(GISAID_sel$LINEAGE2, ref=sel_target_VOC) # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE2)
# [1] "B.1.617.2" "AY.16"     "AY.25"     "AY.26"     "AY.3"      "AY.4"      "B.1.1.519" "B.1.1.7"   "B.1.2"     "B.1.351"   "B.1.429"  
# [12] "B.1.526"   "B.1.621"   "other"     "P.1"   
levels_LINEAGE2 = c("B.1.1.7","B.1.1.519","B.1.2",
                    "B.1.351","B.1.429","B.1.526",
                    "P.1","B.1.621",sel_target_VOC,
                    "AY.3","AY.4","AY.16","AY.25","AY.26","other")
GISAID_sel$LINEAGE2 = factor(GISAID_sel$LINEAGE2, levels=levels_LINEAGE2)

GISAID_sel$state  = GISAID_sel$division
GISAID_sel = GISAID_sel[!grepl("Princess",GISAID_sel$state),] # remove cruise boat data
GISAID_sel = GISAID_sel[!grepl("Islands|Guam",GISAID_sel$state),] # remove Northern Mariana Islands & Virgin Islands & Guam
GISAID_sel = GISAID_sel[!grepl("USA",GISAID_sel$state),] # remove data with unspecified state
GISAID_sel$state[grepl("Washington",GISAID_sel$state)] = "Washington" # Washtington DC -> Washington
sel_states <- levels_states <- c("Arkansas", "Arizona", "California", "Colorado", "Connecticut", "Florida", "Indiana", "Kansas", 
               "Massachusetts", "Missouri", "Mississippi", "Louisiana", "Nebraska", "Nevada", "New Jersey", "New York", "Oklahoma",
               "Texas", "Utah", "Puerto Rico") # 20 states with most data for Delta or increasing cases
length(sel_states)
GISAID_sel = GISAID_sel[GISAID_sel$state %in% sel_states,]

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
GISAID_sel = GISAID_sel[-which(grepl("B.1.617", GISAID_sel$pango_lineage, fixed=TRUE)&GISAID_sel$date<=as.Date("2021-04-14")),]

GISAID_sel$state = factor(GISAID_sel$state, levels=levels_states)
levels(GISAID_sel$state)

table(GISAID_sel$LINEAGE2)
# B.1.1.7 B.1.1.519     B.1.2   B.1.351   B.1.429   B.1.526       P.1   B.1.621 B.1.617.2      AY.3      AY.4     AY.16     AY.25     AY.26 
# 107860      7543     33335       845     26047     25090     12100      2378    128020     26155     38131       385     39529     12216 
# other 
# 113834 

range(GISAID_sel$date) # "2021-01-01" "2021-09-24"

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek2 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$LINEAGE2))
colnames(data_agbyweek2) = c("floor_date", "LINEAGE2", "count")
data_agbyweek2_sum = aggregate(count ~ floor_date, data=data_agbyweek2, sum)
data_agbyweek2$total = data_agbyweek2_sum$count[match(data_agbyweek2$floor_date, data_agbyweek2_sum$floor_date)]
sum(data_agbyweek2[data_agbyweek2$LINEAGE2=="B.1.1.7","total"]) == nrow(GISAID_sel) # correct
data_agbyweek2$collection_date = as.Date(as.character(data_agbyweek2$floor_date))
data_agbyweek2$LINEAGE2 = factor(data_agbyweek2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek2$collection_date_num = as.numeric(data_agbyweek2$collection_date)
data_agbyweek2$prop = data_agbyweek2$count/data_agbyweek2$total
data_agbyweek2$floor_date = NULL


# aggregated by week and state for selected variant lineages
data_agbyweekregion2 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$state, GISAID_sel$LINEAGE2))
colnames(data_agbyweekregion2) = c("floor_date", "division", "LINEAGE2", "count")
data_agbyweekregion2_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion2, sum)
data_agbyweekregion2$total = data_agbyweekregion2_sum$count[match(interaction(data_agbyweekregion2$floor_date,data_agbyweekregion2$division), 
                                                                  interaction(data_agbyweekregion2_sum$floor_date,data_agbyweekregion2_sum$division))]
sum(data_agbyweekregion2[data_agbyweekregion2$LINEAGE2=="B.1.1.7","total"]) == nrow(GISAID_sel) # correct
data_agbyweekregion2$collection_date = as.Date(as.character(data_agbyweekregion2$floor_date))
data_agbyweekregion2$LINEAGE2 = factor(data_agbyweekregion2$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion2$division = factor(data_agbyweekregion2$division, levels=levels_states)
data_agbyweekregion2$collection_date_num = as.numeric(data_agbyweekregion2$collection_date)
data_agbyweekregion2$prop = data_agbyweekregion2$count/data_agbyweekregion2$total
data_agbyweekregion2 = data_agbyweekregion2[data_agbyweekregion2$total!=0,]
data_agbyweekregion2$floor_date = NULL


# MULLER PLOT (RAW DATA)
library(scales)
n2 = length(levels(GISAID_sel$LINEAGE2))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.1.7")] = "#0085FF"
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.351")] = "#9A9D00"
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.177+")] = "grey55"
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="P.1")] = "cyan3"
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(GISAID_sel$LINEAGE)==sel_target_VOC)] = "magenta"
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.3")] = muted("red",c=150,l=80)
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.4")] = muted("red",c=150,l=70)
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.16")] = muted("red",c=150,l=30)
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.25")] = muted("red",c=150,l=20)
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.26")] = muted("red",c=150,l=10)
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.621")] = "red"
lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="other")] = "grey75"



muller_sel_raw2 = ggplot(data=data_agbyweek2, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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
data_agbyweekregion2$LINEAGE2 = relevel(data_agbyweekregion2$LINEAGE2, ref=sel_target_VOC)
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
BIC(fit1_usa_multi, fit2_usa_multi, fit3_usa_multi, fit4_usa_multi, fit5_usa_multi) 
# fit4_usa_multi fits best (lowest BIC)

# growth rate advantage compared to Delta (difference in growth rate per day) 
emtrusa = emtrends(fit4_usa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_usa = data.frame(confint(emtrusa, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrusa$contrasts)$p.value)
delta_r_usa
# contrast     estimate           SE df    asymp.LCL    asymp.UCL      p.value
# 1    B.1.1.7 - B.1.617.2 -0.105052459 0.0014492501 NA -0.107892937 -0.102211981 4.495293e-13
# 2  B.1.1.519 - B.1.617.2 -0.201983188 0.0108313274 NA -0.223212200 -0.180754176 4.495293e-13
# 3      B.1.2 - B.1.617.2 -0.224519339 0.0108573300 NA -0.245799315 -0.203239363 4.495293e-13
# 4    B.1.351 - B.1.617.2 -0.370345839 0.0143871389 NA -0.398544113 -0.342147565 4.495293e-13
# 5    B.1.429 - B.1.617.2 -0.284598137 0.0091023563 NA -0.302438427 -0.266757846 4.495293e-13
# 6    B.1.526 - B.1.617.2 -0.136166827 0.0050362193 NA -0.146037635 -0.126296018 4.495293e-13
# 7        P.1 - B.1.617.2 -0.088006914 0.0024821280 NA -0.092871796 -0.083142033 4.495293e-13
# 8    B.1.621 - B.1.617.2 -0.106636797 0.0062915864 NA -0.118968079 -0.094305514 4.495293e-13
# 9       AY.3 - B.1.617.2 -0.010592749 0.0009811432 NA -0.012515754 -0.008669743 4.938272e-13
# 10      AY.4 - B.1.617.2 -0.021813806 0.0009394341 NA -0.023655063 -0.019972549 4.495293e-13
# 11     AY.16 - B.1.617.2 -0.077107484 0.0104408128 NA -0.097571101 -0.056643867 5.623058e-12
# 12     AY.25 - B.1.617.2 -0.009387186 0.0009900273 NA -0.011327603 -0.007446768 5.013767e-13
# 13     AY.26 - B.1.617.2 -0.018243682 0.0014347446 NA -0.021055730 -0.015431634 4.495293e-13
# 14     other - B.1.617.2  0.003247914 0.0007811895 NA  0.001716811  0.004779018 4.749604e-04

# growth rate advantage compared to Delta (difference in growth rate per day) with simpler model fit3_usa_multi
emtrusa3 = emtrends(fit3_usa_multi, trt.vs.ctrl ~ LINEAGE2,  
                    var="DATE_NUM",  mode="latent",
                    at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_usa3 = data.frame(confint(emtrusa3, 
                                  adjust="none")$contrasts, 
                          p.value=as.data.frame(emtrusa1$contrasts, adjust="none")$p.value)
delta_r_usa3
#                 contrast     estimate           SE  df     lower.CL     upper.CL       p.value
# 1    B.1.1.7 - B.1.617.2 -0.102402936 0.0006579456 308 -0.103697573 -0.101108299  0.000000e+00
# 2  B.1.1.519 - B.1.617.2 -0.145413383 0.0032809434 308 -0.151869282 -0.138957484  0.000000e+00
# 3      B.1.2 - B.1.617.2 -0.202905031 0.0043508781 308 -0.211466236 -0.194343825  0.000000e+00
# 4    B.1.351 - B.1.617.2 -0.118745078 0.0050107190 308 -0.128604649 -0.108885506 7.123708e-206
# 5    B.1.429 - B.1.617.2 -0.176262739 0.0029939944 308 -0.182154010 -0.170371468  0.000000e+00
# 6    B.1.526 - B.1.617.2 -0.112785147 0.0012868563 308 -0.115317289 -0.110253005  0.000000e+00
# 7        P.1 - B.1.617.2 -0.090434626 0.0009887015 308 -0.092380090 -0.088489162 5.730892e-314
# 8    B.1.621 - B.1.617.2 -0.072593441 0.0018707865 308 -0.076274580 -0.068912302 4.056907e-186
# 9       AY.3 - B.1.617.2 -0.006631926 0.0006210513 308 -0.007853966 -0.005409886  1.263945e-32
# 10      AY.4 - B.1.617.2 -0.017614749 0.0004963625 308 -0.018591439 -0.016638058  3.911132e-65
# 11     AY.16 - B.1.617.2 -0.015749962 0.0034847134 308 -0.022606819 -0.008893105  5.161124e-16
# 12     AY.25 - B.1.617.2 -0.007620395 0.0005053321 308 -0.008614735 -0.006626055  1.104530e-15
# 13     AY.26 - B.1.617.2 -0.012822020 0.0007988113 308 -0.014393838 -0.011250202  4.289186e-30
# 14     other - B.1.617.2 -0.009201354 0.0004000061 308 -0.009988444 -0.008414263  0.000000e+00



# growth rate advantage compared to Delta (difference in growth rate per day) with simple model fit1_usa_multi
emtrusa1 = emtrends(fit1_usa_multi, trt.vs.ctrl ~ LINEAGE2,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
delta_r_usa1 = data.frame(confint(emtrusa1, 
                                 adjust="none")$contrasts, 
                         p.value=as.data.frame(emtrusa1$contrasts, adjust="none")$p.value)
delta_r_usa1
# contrast     estimate           SE  df     lower.CL     upper.CL       p.value
# 1    B.1.1.7 - B.1.617.2 -0.054294577 0.0001865859 294 -0.054661790 -0.053927364  0.000000e+00
# 2  B.1.1.519 - B.1.617.2 -0.068947991 0.0003149778 294 -0.069567888 -0.068328094  0.000000e+00
# 3      B.1.2 - B.1.617.2 -0.094231625 0.0002864774 294 -0.094795432 -0.093667819  0.000000e+00
# 4    B.1.351 - B.1.617.2 -0.054949316 0.0006631556 294 -0.056254450 -0.053644182 7.123708e-206
# 5    B.1.429 - B.1.617.2 -0.076720007 0.0002538672 294 -0.077219635 -0.076220380  0.000000e+00
# 6    B.1.526 - B.1.617.2 -0.057038325 0.0002294013 294 -0.057489801 -0.056586848  0.000000e+00
# 7        P.1 - B.1.617.2 -0.046347872 0.0002358205 294 -0.046811982 -0.045883762 5.730892e-314
# 8    B.1.621 - B.1.617.2 -0.031811490 0.0004516679 294 -0.032700403 -0.030922578 4.056907e-186
# 9       AY.3 - B.1.617.2  0.003965319 0.0002938635 294  0.003386976  0.004543662  1.263945e-32
# 10      AY.4 - B.1.617.2 -0.004975790 0.0002232218 294 -0.005415105 -0.004536475  3.911132e-65
# 11     AY.16 - B.1.617.2 -0.013680224 0.0015924531 294 -0.016814276 -0.010546171  5.161124e-16
# 12     AY.25 - B.1.617.2  0.002063584 0.0002433216 294  0.001584711  0.002542457  1.104530e-15
# 13     AY.26 - B.1.617.2 -0.004783933 0.0003739260 294 -0.005519844 -0.004048022  4.289186e-30
# 14     other - B.1.617.2 -0.055619411 0.0001879145 294 -0.055989239 -0.055249583  0.000000e+00


# # pairwise growth rates advantages for simple multinomial model fit1_usa_multi
# emtrusa_pairw1 = emtrends(fit1_usa_multi, pairwise ~ LINEAGE2,  
#                     var="DATE_NUM",  mode="latent",
#                     at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)))
# delta_r_pairw_usa1 = data.frame(confint(emtrusa_pairw1, 
#                                   adjust="none")$contrasts, 
#                           p.value=as.data.frame(emtrusa_pairw1$contrasts, adjust="none")$p.value)
# delta_r_pairw_usa1[grepl("\\<B.1.1.7\\>|\\<B.1.617.2\\>|\\<B.1.351\\>|\\<P.1\\>|^AY.3\\>", delta_r_pairw_usa1$contrast)&
#                    !grepl("\\<B.1.2\\>|\\<B.1.1.519\\>|\\<B.1.351\\>|\\<B.1.427\\>|^B.1.429\\>|^B.1.526\\>", delta_r_pairw_usa1$contrast)
#                    ,]




# fitted prop of different LINEAGES in the USA today
multinom_preds_today_avg = data.frame(emmeans(fit4_usa_multi, ~ LINEAGE2|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#                    LINEAGE2         prob           SE df     asymp.LCL    asymp.UCL
# 1  Delta (B.1.617.2 & AY.X) 9.954898e-01 1.030055e-03 NA  9.934709e-01 9.975086e-01
# 2                   B.1.1.7 2.072669e-04 3.641782e-05 NA  1.358893e-04 2.786445e-04
# 3                     B.1.2 2.840249e-08 6.489041e-08 NA -9.878037e-08 1.555854e-07
# 4                 B.1.1.519 4.146986e-07 1.666438e-06 NA -2.851459e-06 3.680857e-06
# 5                   B.1.351 2.347053e-06 1.807046e-06 NA -1.194693e-06 5.888799e-06
# 6                   B.1.429 4.118771e-07 1.248523e-06 NA -2.035183e-06 2.858937e-06
# 7                   B.1.526 5.829119e-06 4.650586e-06 NA -3.285862e-06 1.494410e-05
# 8                       P.1 4.209269e-04 1.889681e-04 NA  5.055617e-05 7.912977e-04
# 9                   B.1.621 9.550352e-04 3.309310e-04 NA  3.064223e-04 1.603648e-03
# 10                    other 2.917963e-03 9.039799e-04 NA  1.146195e-03 4.689731e-03

# % non-Delta
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#            prob asymp.LCL asymp.UCL 
# 0.004510223 0.001629596 0.007390849 

# fitted prop of delta by state
lineages_bystate = data.frame(emmeans(fit4_usa_multi,
                   ~ LINEAGE2,
                   by=c("DATE_NUM","STATE"),
                   at=list(DATE_NUM=today_num),  
                   mode="prob", df=NA))
lineages_bystate$DATE = as.Date(lineages_bystate$DATE_NUM, origin="1970-01-01")
delta_bystate = lineages_bystate[lineages_bystate$LINEAGE2==sel_target_VOC,]
deltaplusAY3_bystate = delta_bystate
delta_bystate = delta_bystate[order(delta_bystate$prob, decreasing=T),]
delta_bystate
#                     LINEAGE2 DATE_NUM         STATE      prob           SE df asymp.LCL asymp.UCL       DATE
# 111 Delta (B.1.617.2 & AY.X)    18877     Louisiana 0.9998712 0.0000676684 NA 0.9997385 1.0000038 2021-09-07
# 121 Delta (B.1.617.2 & AY.X)    18877      Nebraska 0.9995922 0.0002316194 NA 0.9991383 1.0000462 2021-09-07
# 71  Delta (B.1.617.2 & AY.X)    18877        Kansas 0.9991467 0.0004674791 NA 0.9982304 1.0000629 2021-09-07
# 51  Delta (B.1.617.2 & AY.X)    18877       Florida 0.9987323 0.0002528331 NA 0.9982368 0.9992279 2021-09-07
# 21  Delta (B.1.617.2 & AY.X)    18877    California 0.9986610 0.0001743451 NA 0.9983193 0.9990027 2021-09-07
# 81  Delta (B.1.617.2 & AY.X)    18877 Massachusetts 0.9986015 0.0002660218 NA 0.9980801 0.9991229 2021-09-07
# 1   Delta (B.1.617.2 & AY.X)    18877      Arkansas 0.9985756 0.0005031932 NA 0.9975893 0.9995618 2021-09-07
# 131 Delta (B.1.617.2 & AY.X)    18877        Nevada 0.9982874 0.0006419429 NA 0.9970292 0.9995456 2021-09-07
# 61  Delta (B.1.617.2 & AY.X)    18877       Indiana 0.9981395 0.0008425738 NA 0.9964881 0.9997909 2021-09-07
# 41  Delta (B.1.617.2 & AY.X)    18877   Connecticut 0.9979875 0.0008114996 NA 0.9963970 0.9995781 2021-09-07
# 151 Delta (B.1.617.2 & AY.X)    18877      New York 0.9974086 0.0004838172 NA 0.9964603 0.9983568 2021-09-07
# 101 Delta (B.1.617.2 & AY.X)    18877   Mississippi 0.9971460 0.0014881653 NA 0.9942292 1.0000627 2021-09-07
# 171 Delta (B.1.617.2 & AY.X)    18877         Texas 0.9966140 0.0005301810 NA 0.9955749 0.9976531 2021-09-07
# 141 Delta (B.1.617.2 & AY.X)    18877    New Jersey 0.9960966 0.0012753263 NA 0.9935970 0.9985962 2021-09-07
# 181 Delta (B.1.617.2 & AY.X)    18877          Utah 0.9957981 0.0011178364 NA 0.9936072 0.9979890 2021-09-07
# 11  Delta (B.1.617.2 & AY.X)    18877       Arizona 0.9955032 0.0012323561 NA 0.9930878 0.9979186 2021-09-07
# 91  Delta (B.1.617.2 & AY.X)    18877      Missouri 0.9917820 0.0045478568 NA 0.9828683 1.0006956 2021-09-07
# 31  Delta (B.1.617.2 & AY.X)    18877      Colorado 0.9892757 0.0063308020 NA 0.9768675 1.0016838 2021-09-07
# 191 Delta (B.1.617.2 & AY.X)    18877   Puerto Rico 0.9823949 0.0118560008 NA 0.9591576 1.0056323 2021-09-07
# 161 Delta (B.1.617.2 & AY.X)    18877      Oklahoma 0.9801815 0.0143701588 NA 0.9520165 1.0083465 2021-09-07

# AY3_bystate = lineages_bystate[lineages_bystate$LINEAGE2=="AY.3",]
# deltaplusAY3_bystate$prob = deltaplusAY3_bystate$prob + AY3_bystate$prob
# AY3_bystate = AY3_bystate[order(AY3_bystate$prob, decreasing=T),]
# AY3_bystate
# 
# 
# deltaplusAY3_bystate = deltaplusAY3_bystate[order(deltaplusAY3_bystate$prob, decreasing=T),]
# deltaplusAY3_bystate


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
date.to = as.numeric(as.Date("2021-09-21")) # max(GISAID_sel$DATE_NUM)+extrapolate

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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01"))),1,1),
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

