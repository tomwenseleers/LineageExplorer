# ANALYSIS OF GROWTH ADVANTAGE OF OMICRON (B.1.1.529) IN SOUTH AFRICA 
# (GISAID RECORDS & 
# SGTF data (now proxy for Omicron / B.1.1.529), traced from graph by Alex Selby from briefing by Tulio de Oliveira & Richard Lessells,
# original data courtesy of Lesley Scott & NHLS team
# https://www.youtube.com/watch?v=Vh4XMueP1zQ, https://twitter.com/chrischirp/status/1463885565221384202/photo/1)

# T. Wenseleers
# last update 6 DECEMBER 2021
    
  library(nnet)
  # devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
  library(mclogit)
  # remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
  library(emmeans)
  library(readr)
  library(ggplot2)
  library(ggthemes)
  library(scales)
  library(stringr)
  library(lubridate)
  library(dplyr)  
  library(splines)
  library(tidyr)
  library(tidyselect)
  library(effects)
  library(MASS)
  library(nlme)
  library(lme4)
    
  today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
  # today = as.Date("2021-12-06")
  today_num = as.numeric(today)
  plotdir = "South Africa_GISAID"
  suppressWarnings(dir.create(paste0(".//plots//",plotdir)))
  
  # X axis for plots
  firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
  xaxis = scale_x_continuous(breaks=firststofmonth,
                             labels=substring(months(firststofmonth),1,1),
                             expand=c(0,0))

# import GISAID records for South Africa
d1 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_11_25_20_subm_2020.tsv"), col_types = cols(.default = "c")) 
d2 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_11_25_20_subm_jan_may_2021.tsv"), col_types = cols(.default = "c")) 
d3 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_11_25_20_subm_june_aug_2021.tsv"), col_types = cols(.default = "c")) 
d4 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_12_06_11_subm_sept_dec_2021.tsv"), col_types = cols(.default = "c")) 



# import SGTF data (now proxy for B.1.1.529 / Omicron)
# sgtf = read.csv(".//data//GISAID//South Africa//South Africa SGTF.csv") # SGTF data, from graph shown at press conference, now proxy for B.1.1.529 (traced by John Burn Murdoch)
sgtf = read.csv(".//data//GISAID//South Africa//SA_sgtf_Alex Selby.csv") # SGTF data, from graph shown at press conference, now proxy for B.1.1.529 (traced by Alex Selby, https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/EarlyOmicronEstimate/SA_sgtf)
# note : technically, columns sgtf and non_sgtf should be counts and they are not - this is due to rounding errors caused by tracing the data from the graph
# this should not affect estimates much though
sgtf$date = as.Date(sgtf$date)
sgtf$date_num = as.numeric(sgtf$date)
sgtf$prop = sgtf$sgtf/sgtf$tests
sgtf$baseline = 0.02 # 2% baseline subtracted
sgtf$prop = sgtf$prop - sgtf$baseline
sgtf$prop[sgtf$prop<0] = 0
sgtf$sgtf = sgtf$prop*sgtf$tests
sgtf$non_sgtf = (1-sgtf$prop)*sgtf$tests
sgtf$obs = as.factor(1:nrow(sgtf)) # for observation-level random effect to take into account overdispersion
sgtf$random = factor(1) # fake random effect to be able to run glmmPQL
names(sgtf)


# ANALYSIS OF SGTF DATA USING LOGISTIC REGRESSION ####

sgtf_subs = sgtf[sgtf$date>=as.Date("2021-10-24"),] # SGTF data sinc Oct 25 when B.1.1.529 started to spread
# B_1_1_529_nov = data.frame(variant="B.1.1.529", date=sgtf_nov$date, "count"=sgtf_nov$SGTF-2) # we subtract a 2% baseline, "count" is actually a percentage here

# fit binomial GLM to SGTF data
fit_sgtf = glm(cbind(sgtf, non_sgtf) ~ date_num, family=binomial(logit), data=sgtf_subs)
summary(fit_sgtf)

# fit_sgtf = glmer(cbind(sgtf, non_sgtf) ~ (1|obs)+date_num, family=binomial(logit), data=sgtf_nov) # taking into account overdispersion via random effect
# summary(fit_sgtf)

# fit_sgtf = glmmPQL(cbind(sgtf, non_sgtf) ~ date_num, family=binomial(logit), correlation=corAR1(), random=~1|obs, data=sgtf_nov) # taking into account lag-1 autocorrelation in residuals
# summary(fit_sgtf)

# plot(allEffects(fit_sgtf))

dateseq = seq(as.Date("2021-10-24"), today, by=1)
emmeans_sgtf = as.data.frame(emmeans(fit_sgtf, ~ date_num, at=list(date_num=dateseq)), type="response")
colnames(emmeans_sgtf) = c("date_num","prob","SE","df","asymp.LCL","asymp.UCL")
emmeans_sgtf$date = as.Date(emmeans_sgtf$date_num, origin="1970-01-01")

deltar_sgtf = as.data.frame(emtrends(fit_sgtf, ~ date_num, var="date_num", at=list(date_num=max(dateseq))))[,c(2,5,6)] # growth rate advantage of Omicron over Delta
deltar_sgtf_char = sapply(deltar_sgtf, function (x) sprintf(as.character(round(x,2)), 2) )
deltar_sgtf_char = paste0(deltar_sgtf_char[1], " [", deltar_sgtf_char[2], "-", deltar_sgtf_char[3],"] 95% CLs")
transmadv_sgtf = exp(deltar_sgtf*4.7) # transmission advantage with gen time of 4.7 days
transmadv_sgtf_char = sapply(transmadv_sgtf, function (x) sprintf(as.character(round(x,1)), 1) )
transmadv_sgtf_char = paste0(transmadv_sgtf_char[1], " [", transmadv_sgtf_char[2], "-", transmadv_sgtf_char[3],"] 95% CLs")

qplot(data=sgtf_subs, x=date, y=prop, geom="point", colour=I("steelblue"), size=I(2)) + 
  geom_ribbon(data=emmeans_sgtf, aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), alpha=I(0.5), fill=I("steelblue")) +
  geom_line(data=emmeans_sgtf, aes(y=prob), colour=I("steelblue"), alpha=I(0.5), size=1) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  xlim(as.Date("2021-10-14"), today) +
  coord_cartesian(ylim=c(0.001,0.99)) +
  ggtitle("Spread of Omicron in South Africa inferred from SGTF data", "Data courtesy of Lesley Scott & NHLS team\ntraced from original graph by Alex Selby\n(2% baseline subtracted)\n\nLogistic fit by Tom Wenseleers") +
  ylab("Share of Omicron among confirmed cases (%)") +
  annotate("text", x = c(as.Date("2021-10-14")), y = c(0.94), 
                      label = paste0("Growth rate advantage of Omicron over Delta:\n", deltar_sgtf_char,
                                     " per day\n\nTransmission advantage Omicron over Delta\n(with generation time of 4.7 days):\n", transmadv_sgtf_char),
                      color="black", hjust=0, size=3) +
  theme_hc() +
  xlab("")

ggsave(file=paste0(".\\plots\\",plotdir,"\\sgtf data_spread omicron logistic fit.png"), width=8, height=6)






# parse GISAID & SGTF data
GISAID = as.data.frame(rbind(d1,d2,d3,d4))
colnames(GISAID) = c("Virus name","Accession ID","date","Location","host",
                     "Additional location information","Sampling strategy","Gender",                         
                     "Patient age","Patient status","Last vaccinated","Passage","Specimen",
                     "Additional host information","pango_lineage","Clade","aa_substitutions")
date_isvalid = sapply(GISAID$date, function (s) str_count(s, pattern = "-")==2)
GISAID = GISAID[date_isvalid,]
GISAID$date = as.Date(GISAID$date) 
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID = GISAID[GISAID$date>=as.Date("2020-01-01"),]
range(GISAID$date) # "2020-03-06" "2021-11-27"
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$date_num = as.numeric(GISAID$date)
# GISAID = GISAID[GISAID$pango_lineage!="None",]
attach(GISAID)

GISAID$variant = case_when(
  # grepl("N679K", aa_substitutions) & grepl("H655Y", aa_substitutions) & grepl("P681H", aa_substitutions) ~ "B.1.1.529",
  grepl("B.1.1.529", pango_lineage, fixed=T) ~ "Omicron",
  grepl("B.1.617.2", pango_lineage, fixed=T) | grepl("AY", pango_lineage)  ~ "Delta",
  grepl("B.1.1.7", pango_lineage, fixed=T) ~ "Alpha",
  grepl("B.1.351", pango_lineage, fixed=T) ~ "Beta",
  grepl("C.1.2", pango_lineage, fixed=T) ~ "C.1.2",
  T ~ "Other"
)

table(GISAID[GISAID$variant=="Other"&GISAID$date>as.Date("2021-10-01"),]$pango_lineage)

GISAID = GISAID[!(GISAID$variant=="Other"&GISAID$pango_lineage=="None"),]

table(GISAID$variant)

# ANALYSIS OF VOCs IN SOUTH AFRICA ####

sel_target_VOC = "Omicron"
sel_reference_VOC = "Delta"

GISAID_sel = GISAID
nrow(GISAID_sel) # 22881
sum(GISAID_sel$variant==sel_target_VOC) # 227
table(GISAID_sel$variant)
# Alpha    Beta   C.1.2   Delta Omicron   Other 
# 224    6862     264   10907     227    4397 
range(GISAID_sel$date) # "2020-03-06" "2021-11-27"

GISAID_sel$variant = factor(GISAID_sel$variant)
GISAID_sel$variant = relevel(GISAID_sel$variant, ref=sel_reference_VOC) # we code Delta as the reference
levels_variant = c(sel_reference_VOC, "Beta", "Alpha", "C.1.2", "Other", sel_target_VOC)
GISAID_sel$variant = factor(GISAID_sel$variant, levels=levels_variant)
table(GISAID_sel$variant, GISAID_sel$`Sampling strategy`)

# age distribution of patients infected with Delta or Omicron (GISAID data, since Oct 1 2021)
ages_omicron = as.numeric(GISAID_sel[GISAID_sel$variant=="Omicron","Patient age"])
ages_delta = as.numeric(GISAID_sel[GISAID_sel$variant=="Delta"&GISAID_sel$date>=as.Date("2021-10-01"),"Patient age"])  
ages = rbind(data.frame(variant="Omicron", age=ages_omicron), 
             data.frame(variant="Delta", age=ages_delta))
library(ggplot2)
library(ggthemes)
ggplot(data=ages, aes(x=age, fill=variant)) + facet_wrap(~ variant, ncol=1) + geom_histogram() +
  scale_fill_manual(values=c("blue","red")) + theme(legend.position="none") + ggtitle("Age distribution of patients with sequenced Delta & Omicron infections\nin South Africa (GISAID data since Oct 1 2021)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\age distribution patients with delta omicron.png"), width=8, height=6)



# # for november I extrapolate a multinomial fit to the GISAID & SGTF data and use that to fill in the missing values for the frequencies of the other non-B.1.1.529 lineages
# # this would be similar to using an EM algorithm to fill in missing data
# set.seed(1)
# fit1_southafrica_multi0 = nnet::multinom(variant ~ scale(date_num), data=GISAID_sel[GISAID_sel$variant!="Omicron",], maxit=1000)
# fit2_southafrica_multi0 = nnet::multinom(variant ~ ns(date_num, df=2), data=GISAID_sel[GISAID_sel$variant!="Omicron",], maxit=1000)
# BIC(fit1_southafrica_multi0, fit2_southafrica_multi0) 
# #                         df      BIC
# # fit1_southafrica_multi0  8 22379.62
# # fit2_southafrica_multi0 12 19112.50 # best

# # multinomial model predictions for november, merged with SGTF data for november as proxy for B.1.1.529
# predgrid = expand.grid(list(date_num=seq(as.Date("2021-11-01"), max(sgtf$date), by=1)))
# data_nov = data.frame(predgrid, as.data.frame(predict(fit2_southafrica_multi0, newdata=predgrid, type="prob")),check.names=F)
# data_nov = gather(data_nov, variant, prob, all_of(levels_variant[levels_variant!="Omicron"]), factor_key=TRUE)
# data_nov$date = as.Date(data_nov$date_num, origin="1970-01-01")
# data_nov$date_num = NULL
# data_nov$count = 100*data_nov$prob*(1-B_1_1_529_nov$count/100)
# data_nov$prob = NULL
# head(data_nov)
# data_nov = rbind(data_nov, B_1_1_529_nov)
# data_nov$variant = factor(data_nov$variant, levels=levels_variant) 
# data_nov$total = 100


# aggregated data to make Muller plots of raw data
# aggregate by day to identify days on which INSA performed (days with a lot of sequences)
# we subset the data to just those days to avoid sampling biases (delta infection clusters etc)
data_agbyday = as.data.frame(table(GISAID_sel$date, GISAID_sel$variant))
colnames(data_agbyday) = c("date", "variant", "count")
data_agbyday_sum = aggregate(count ~ date, data=data_agbyday, sum)
data_agbyday$total = data_agbyday_sum$count[match(data_agbyday$date, data_agbyday_sum$date)]
data_agbyday$date = as.Date(as.character(data_agbyday$date))
# data_agbyday = rbind(data_agbyday, data_nov) # merge with SGTF+extrapolated GISAID for november
data_agbyday$variant = factor(data_agbyday$variant, levels=levels_variant)
data_agbyday$date_num = as.numeric(data_agbyday$date)
data_agbyday$prop = data_agbyday$count/data_agbyday$total
data_agbyday$floor_date = NULL


# # GISAID data merged with SGTF data
# GISAID_sel2 = GISAID_sel[GISAID_sel$date<as.Date("2021-11-01"), c("date","variant")]
# data_agbyday = as.data.frame(table(GISAID_sel2$date, GISAID_sel2$variant))
# colnames(data_agbyday) = c("date", "variant", "count")
# data_agbyday_sum = aggregate(count ~ date, data=data_agbyday, sum)
# data_agbyday$total = data_agbyday_sum$count[match(data_agbyday$date, data_agbyday_sum$date)]
# data_agbyday$date = as.Date(as.character(data_agbyday$date))
# data_agbyday = rbind(data_agbyday, data_nov) # merge with SGTF+extrapolated GISAID for november
# data_agbyday$variant = factor(data_agbyday$variant, levels=levels_variant)
# data_agbyday$date_num = as.numeric(data_agbyday$date)
# data_agbyday$prop = data_agbyday$count/data_agbyday$total
# data_agbyday$floor_date = NULL
# # GISAID_sel$total_sequenced_on_that_day = data_agbyday$total[match(GISAID_sel$date, data_agbyday$date)]

# # aggregated by week for selected variant lineages
data_agbyweek = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$variant))
colnames(data_agbyweek) = c("floor_date", "variant", "count")
data_agbyweek_sum = aggregate(count ~ floor_date, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$floor_date, data_agbyweek_sum$floor_date)]
sum(data_agbyweek[data_agbyweek$variant=="Beta","total"]) == nrow(GISAID_sel) # TRUE
data_agbyweek$date = as.Date(as.character(data_agbyweek$floor_date))
data_agbyweek$variant = factor(data_agbyweek$variant, levels=levels_variant)
data_agbyweek$date_num = as.numeric(data_agbyweek$date)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
data_agbyweek$floor_date = NULL
data_agbyweek$date_num = as.numeric(data_agbyweek$date)


# MULLER PLOT (RAW DATA)
  n2 = length(levels(GISAID_sel$variant))
  lineage_cols2 = hcl(h = seq(0, 180, length = n2), l = 60, c = 180)
  lineage_cols2[which(levels(GISAID_sel$variant)=="Alpha")] = "#0085FF"
  lineage_cols2[which(levels(GISAID_sel$variant)=="Beta")] = "green4"
  lineage_cols2[which(levels(GISAID_sel$variant)=="Delta")] = "mediumorchid"
  lineage_cols2[which(levels(GISAID_sel$variant)=="C.1.2")] = "darkorange"
  lineage_cols2[which(levels(GISAID_sel$variant)==sel_target_VOC)] = "red2" # "magenta"
  lineage_cols2[which(levels(GISAID_sel$variant)=="Other")] = "grey65"


  
muller_southafrica_raw2 = ggplot(data=data_agbyweek, aes(x=date, y=count, group=variant)) + 
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  xaxis +
  scale_y_continuous(expand=c(0,0)) +
  theme_hc() +
  # labs(title = "SARS-CoV2 VARIANTS IN SOUTH AFRICA") +
  ylab("Share") + 
  theme_hc() +
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_southafrica_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_raw data.pdf"), width=8, height=6)


# multinomial fits
library(nnet)
library(splines)
set.seed(1)
# fit1_southafrica_multi = nnet::multinom(variant ~ scale(date_num), weights=count, data=data_agbyday[data_agbyday$variant!="Other",], maxit=1000)
# fit2_southafrica_multi = nnet::multinom(variant ~ ns(date_num, df=2), weights=count, data=data_agbyday[data_agbyday$variant!="Other",], maxit=1000)
# fit3_southafrica_multi = nnet::multinom(variant ~ ns(date_num, df=3), weights=count, data=data_agbyday[data_agbyday$variant!="Other",], maxit=1000)
fit1_southafrica_multi = nnet::multinom(variant ~ scale(date_num), weights=count, data=data_agbyday, maxit=1000)
fit2_southafrica_multi = nnet::multinom(variant ~ ns(date_num, df=2), weights=count, data=data_agbyday, maxit=1000)
BIC(fit1_southafrica_multi, fit2_southafrica_multi) 
#                        df      BIC
# fit1_southafrica_multi 10 22655.05
# fit2_southafrica_multi 15 19290.92

# fit1_southafrica_multi = nnet::multinom(variant ~ scale(date_num), weights=count, data=data_agbyday[(data_agbyday$variant!="Other"),], maxit=1000)
# fit2_southafrica_multi = nnet::multinom(variant ~ ns(date_num, df=2), weights=count, data=data_agbyday[(data_agbyday$variant!="Other"),], maxit=1000)
# BIC(fit1_southafrica_multi, fit2_southafrica_multi) 


# growth rate advantage of Omicron over Delta (difference in growth rate per day) 
emtrsouthafrica = emtrends(fit2_southafrica_multi, trt.vs.ctrl ~ variant,  
                   var="date_num",  mode="latent",
                   at=list(date_num=today_num),
                   adjust="none", df=NA)
delta_r_southafrica = data.frame(confint(emtrsouthafrica, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrsouthafrica$contrasts,
                                               adjust="none", df=NA)$p.value)
delta_r_southafrica
# contrast     estimate          SE df   asymp.LCL    asymp.UCL      p.value
# 1    Beta - Delta -0.030074041 0.003542634 NA -0.037017476 -0.023130605 2.081088e-17
# 2   Alpha - Delta -0.044163917 0.007392891 NA -0.058653718 -0.029674116 2.317395e-09
# 3   C.1.2 - Delta -0.001982514 0.005068634 NA -0.011916854  0.007951825 6.956983e-01
# 4   Other - Delta  0.015626469 0.003099793 NA  0.009550986  0.021701951 4.627837e-07
# 5 Omicron - Delta  0.245066559 0.022525988 NA  0.200916433  0.289216684 1.446997e-27

# based on GISAID data alone: 3.2x [2.6-3.9]x higher effective R value of Omicron compared to Delta
exp(0.200916433*4.7) # 2.6x
exp(0.245066559*4.7) # Omicron 3.2x higher transmission than Delta
exp(0.289216684*4.7) # 3.9

# # based on GISAID+SGTF data
# exp(0.33*4.7) # 4.7x
# exp(0.38*4.7) # 6x higher effective R
# exp(0.43*4.7) # 7.5x



# fitted prop of different variantS today
multinom_preds_today_avg = data.frame(emmeans(fit2_southafrica_multi, ~ variant|1,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# variant         prob           SE df     asymp.LCL    asymp.UCL
# 1   Delta 1.965539e-03 1.028162e-03 NA -4.962175e-05 3.980699e-03
# 2    Beta 3.033600e-07 1.999241e-07 NA -8.848407e-08 6.952041e-07
# 3   Alpha 1.089642e-08 1.227779e-08 NA -1.316761e-08 3.496046e-08
# 4   C.1.2 5.110366e-05 3.383546e-05 NA -1.521262e-05 1.174199e-04
# 5   Other 8.198092e-05 4.996481e-05 NA -1.594832e-05 1.799102e-04
# 6 Omicron 9.979011e-01 1.097481e-03 NA  9.957500e-01 1.000052e+00
# (note that there is a bit of sampling bias towards Gauteng Province though)


# PLOT MULTINOMIAL FIT

extrapolate = 30
date.from = as.numeric(min(GISAID_sel$date_num))
date.to = max(GISAID_sel$date_num)+extrapolate

# multinomial model predictions (fastest, but no confidence intervals)
predgrid = expand.grid(list(date_num=seq(date.from, date.to)))
fit_southafrica_multi_preds = data.frame(predgrid, as.data.frame(predict(fit2_southafrica_multi, newdata=predgrid, type="prob")),check.names=F)
fit_southafrica_multi_preds = gather(fit_southafrica_multi_preds, variant, prob, all_of(levels_variant), factor_key=TRUE)
fit_southafrica_multi_preds$date = as.Date(fit_southafrica_multi_preds$date_num, origin="1970-01-01")
fit_southafrica_multi_preds$variant = factor(fit_southafrica_multi_preds$variant, levels=levels_variant) 

muller_southafrica_mfit = ggplot(data=fit_southafrica_multi_preds, 
                                   aes(x=date, y=prob, group=variant)) + 
  scale_y_continuous(expand=c(0,0)) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$date_num)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.2, fill="white") + # extrapolated part
  xaxis +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID data, multinomial fit)")
muller_southafrica_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_multinom fit.png"), width=10, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_multinom fit.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_southafrica_raw2 + coord_cartesian(xlim=c(min(GISAID_sel$date),max(GISAID_sel$date)+extrapolate))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_southafrica_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots multipanel_multinom fit.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots multipanel_multinom fit.pdf"), width=10, height=10)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions with confidence intervals (but slower)
fit_southafrica_multi_preds_withCI = data.frame(emmeans(fit2_southafrica_multi,
                                                        ~ variant,
                                                        by=c("date_num"),
                                                        at=list(date_num=seq(date.from, date.to, by=1)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_southafrica_multi_preds_withCI$date = as.Date(fit_southafrica_multi_preds_withCI$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_withCI$variant = factor(fit_southafrica_multi_preds_withCI$variant, levels=levels_variant)
fit_southafrica_multi_preds2 = fit_southafrica_multi_preds_withCI


# on logit scale:

ymin = 0.001
ymax = 0.999
fit_southafrica_multi_preds2$asymp.LCL[fit_southafrica_multi_preds2$asymp.LCL<ymin] = ymin
fit_southafrica_multi_preds2$asymp.UCL[fit_southafrica_multi_preds2$asymp.UCL<ymin] = ymin
fit_southafrica_multi_preds2$asymp.UCL[fit_southafrica_multi_preds2$asymp.UCL>ymax] = ymax
fit_southafrica_multi_preds2$prob[fit_southafrica_multi_preds2$prob<ymin] = ymin

plot_southafrica_mfit_logit = qplot(data=fit_southafrica_multi_preds2, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID data, multinomial fit)") +
  xaxis +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek,
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 4), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_southafrica_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_logit scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_southafrica_mfit = qplot(data=fit_southafrica_multi_preds2, x=date, y=100*prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID data, multinomial fit)") +
  xaxis +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek,
             aes(x=date, y=100*prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\ngenotyped", trans="sqrt",
                        range=c(0.5, 4), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_southafrica_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_response scale_with points.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_response scale.pdf"), width=10, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_regional = get_regional_data(country = "South Africa")
cases_regional = as.data.frame(get_regional_data(country = "South Africa"))[,c("date","province","cases_new","cases_total","deaths_new","deaths_total")]
cases_regional$date_num = as.numeric(cases_regional$date)
cases_regional$WEEKDAY = weekdays(cases_regional$date)
# cases_regional = cases_regional[!is.na(cases_regional$date),]
cases_regional = cases_regional[complete.cases(cases_regional),]

cases_gauteng = cases_regional[cases_regional$province=="Gauteng",]
tail(cases_gauteng)

# case & testing data for SA overall
cases_cum = read.csv("https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=cumulative")
cases_cum = cases_cum[cases_cum$date!="",]
head(cases_cum)
cases_cum$date = as.Date(cases_cum$date, "%d-%m-%Y")
qplot(data=cases_cum, x=date, y=cases_daily, geom="col")
names(cases_cum)
qplot(data=cases_cum, x=date, y=tests_daily, geom="col")

# for more data see
# https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=provincial-tests
# https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=deaths-age
# https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=deaths-gender

# hospitalisation data for SA overall (NICD) (no use: only runs until sept)
hosp = read.csv("https://raw.githubusercontent.com/dsfsi/covid19za/master/data/nicd_hospital_surveillance_data.csv")
hosp$date = as.Date(hosp$date, "%d-%m-%Y")
qplot(data=hosp, x=date, y=total_admissions, geom="col")
hosp$new_admissions = c(0,diff(hosp$total_admissions))
qplot(data=hosp, x=date, y=new_admissions, geom="col") + xaxis
# qplot(data=hosp, x=date, y=current_num_in_hospital, geom="col") + xaxis
# qplot(data=hosp, x=date, y=WC, geom="col") = xaxis

# case data for SA per province
cases_prov = read.csv("https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=provinces")
cases_prov = cases_prov[cases_prov$date!="",]
cases_prov = cases_prov[cases_prov$province!="Unknown",]
cases_prov$province = factor(cases_prov$province, levels=c("Gauteng", "North West", "Mpumalanga", "Limpopo", "Western Cape", "Free State", "KwaZulu Natal", "Eastern Cape", "Northern Cape"))
head(cases_prov)
cases_prov$date = as.Date(cases_prov$date)
qplot(data=cases_prov, x=date, y=daily_cases, geom="col") + facet_wrap(~ province, scale="free_y") # + scale_y_log10()
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases by province.png"), width=10, height=6)
qplot(data=cases_prov, x=date, y=daily_deaths, geom="col") + facet_wrap(~ province, scale="free_y") # + scale_y_log10()
ggsave(file=paste0(".\\plots\\",plotdir,"\\deaths by province.png"), width=10, height=6)
names(cases_cum)
qplot(data=cases_cum, x=date, y=tests_daily, geom="col")



# cases_tot = as.data.frame(get_national_data(countries = "South Africa"))
# # cases_tot = cases_tot[cases_tot$date>=as.Date("2020-01-01"),]
# cases_tot$date_num = as.numeric(cases_tot$date)
# # cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
# cases_tot$WEEKDAY = weekdays(cases_tot$date)
# cases_tot[cases_tot$date==as.Date("2021-11-24"),"cases_new"] = 868 # https://twitter.com/nicd_sa/status/1463200722615513093
# # cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
# range(cases_tot$date) # "2020-01-03" "2021-12-03"

cases_tot = cases_cum
cases_tot$cases_new = cases_tot$cases_daily
cases_tot$tests_new = abs(cases_tot$tests_daily)
cases_tot$date_num = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot$testspercase = (cases_tot$tests_new+1)/(cases_tot$cases_new+1)
qplot(data=cases_tot, x=date, y=cases_new, geom="col")
qplot(data=cases_tot, x=date, y=tests_new, geom="col")
qplot(data=cases_tot, x=date, y=testspercase, geom="col")
qplot(data=cases_tot, x=date, y=1/testspercase, geom="col")


# CALCULATE Re VALUES THROUGH TIME

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

# smooth out weekday effects in case nrs using negative binomial GAM (possibility to also correct for variable testing intensity)
library(mgcv)
k=40
fit_cases = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) + 
                  WEEKDAY, # + offset(log(testspercase)),
                # BANKHOLIDAY,
                # + s(tests_new, bs="cs", k=8, fx=F),
                family=nb(), data=cases_tot,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+1*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-1*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 
BIC(fit_cases)

# calculate average instantaneous growth rates & 95% CLs & Re values using emtrends ####
# based on the slope of the GAM fit on a log link scale
date.from = as.numeric(as.Date("2020-03-14")) # as.numeric(min(GISAID_sel$date_num))
date.to = today_num+7 # max(GISAID_sel$date_num)+extrapolate

avg_r_cases = as.data.frame(emtrends(fit_cases, ~ date_num, var="date_num", 
                                     at=list(date_num = seq(date.from,
                                                          date.to, by=1),
                                             tests_new = max(cases_tot$tests_daily,na.rm=T),
                                             testspercase = 20),
                                     type="link"))
colnames(avg_r_cases)[2] = "r"
colnames(avg_r_cases)[5] = "r_LOWER"
colnames(avg_r_cases)[6] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$date_num, origin="1970-01-01") # date of diagnosis
avg_r_cases$DATE_OF_INFECTION = avg_r_cases$DATE-7 # date of infection
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE_OF_INFECTION, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) + 
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  # xaxis +
  # scale_y_continuous(limits=c(1/2, 4), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES IN SOUTH AFRICA","Based on negative binomial GAM fit to case data NICD,\nassuming gamma distr generation time of 4.7d with SD of 2.9d\nand time from infection to diagnosis of 7d") +
  labs(tag = "Tom Wenseleers\n4 Dec 2021") +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(0.25,NA))
  # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

head(avg_r_cases,40)
avg_r_cases[avg_r_cases$DATE_OF_INFECTION==today,"Re"]/avg_r_cases[avg_r_cases$DATE_OF_INFECTION==as.Date("2021-09-23"),"Re"] # Re now = 3.34 / Re on 23 Sept = 0.74 = x4.5
avg_r_cases[avg_r_cases$DATE_OF_INFECTION==today,"Re"]/avg_r_cases[avg_r_cases$DATE_OF_INFECTION==as.Date("2021-10-07"),"Re"] # Re now (3.34) / Re on 7 Oct (0.79) = x4.2

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values South Africa.png"), width=8, height=6)


# FOR GAUTENG PROVINCE
k=40
fit_cases_gauteng = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) + 
                  WEEKDAY, # + 
                # BANKHOLIDAY,
                # s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=nb(), data=cases_gauteng,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+1*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-1*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 
BIC(fit_cases_gauteng)

avg_r_cases_gauteng = as.data.frame(emtrends(fit_cases_gauteng, ~ date_num, var="date_num", 
                                             at=list(date_num=seq(date.from,
                                                                  date.to)#,
                                                     # BANKHOLIDAY="no"
                                             ), # weekday="Wednesday",
                                             type="link"))
colnames(avg_r_cases_gauteng)[2] = "r"
colnames(avg_r_cases_gauteng)[5] = "r_LOWER"
colnames(avg_r_cases_gauteng)[6] = "r_UPPER"
avg_r_cases_gauteng$DATE = as.Date(avg_r_cases_gauteng$date_num, origin="1970-01-01") 
avg_r_cases_gauteng$DATE_OF_INFECTION = avg_r_cases_gauteng$DATE-7
avg_r_cases_gauteng$Re = Re.from.r(avg_r_cases_gauteng$r)
avg_r_cases_gauteng$Re_LOWER = Re.from.r(avg_r_cases_gauteng$r_LOWER)
avg_r_cases_gauteng$Re_UPPER = Re.from.r(avg_r_cases_gauteng$r_UPPER)
avg_r_cases_gauteng = avg_r_cases_gauteng[complete.cases(avg_r_cases_gauteng),]
qplot(data=avg_r_cases_gauteng, x=DATE_OF_INFECTION, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) + # -7 TO CALCULATE BACK TO INFECTION DATE
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  # xaxis +
  scale_y_continuous(limits=c(1/4, 4), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES IN GAUTENG PROVINCE, SA","Based on negative binomial GAM fit to case data NICD,\nassuming gamma distr generation time of 4.7d with SD of 2.9d\nand time from infection to diagnosis of 7d") +
  labs(tag = "Tom Wenseleers\n2 Dec 2021") +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

head(avg_r_cases_gauteng,40)
avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==today,"Re"]/avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==as.Date("2021-09-23"),"Re"] # Re now / Re on 23 Sept = x2.8
avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==today,"Re"]/avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==as.Date("2021-10-07"),"Re"] # Re now / Re on 14 Oct = x2.8

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values Gauteng Province.png"), width=8, height=6)
tail(avg_r_cases_gauteng)




# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_southafrica_multi_preds_withCI$totcases = cases_tot$cases_new[match(round(fit_southafrica_multi_preds_withCI$date_num),cases_tot$date_num)]
fit_southafrica_multi_preds_withCI$cases = fit_southafrica_multi_preds_withCI$totcases * fit_southafrica_multi_preds_withCI$prob
fit_southafrica_multi_preds_withCI$cases[fit_southafrica_multi_preds_withCI$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ date_num, at=list(date_num=seq(date.from, date.to, by=1),
                                                                     testspercase=20,
                                                                     tests_new=max(cases_tot$tests_new, na.rm=T)), 
                                      type="response"))
fit_southafrica_multi_preds_withCI$smoothed_totcases = cases_emmeans$response[match(fit_southafrica_multi_preds_withCI$date_num,cases_emmeans$date_num)]
fit_southafrica_multi_preds_withCI$smoothed_cases = fit_southafrica_multi_preds_withCI$smoothed_totcases * fit_southafrica_multi_preds_withCI$prob
fit_southafrica_multi_preds_withCI$smoothed_cases[fit_southafrica_multi_preds_withCI$smoothed_cases<=0.001] = NA
fit_southafrica_multi_preds_withCI$variant = factor(fit_southafrica_multi_preds_withCI$variant, levels=levels_variant)

fit_southafrica_multi_preds_withCI[fit_southafrica_multi_preds_withCI$date==as.Date("2021-12-03"),] # this date is not plotted in plot below FIX

ggplot(data=fit_southafrica_multi_preds_withCI, 
       aes(x=date, y=cases, group=variant)) + 
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant, width=I(1.1)), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA\n(case data NICD & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) # +
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)
write.csv(fit_southafrica_multi_preds_withCI, file=paste0(".\\plots\\",plotdir,"\\cases per day by variant South Africa 6 dec 2021.csv"), row.names=F)

ggplot(data=fit_southafrica_multi_preds_withCI[fit_southafrica_multi_preds_withCI$date<=today,], 
       aes(x=date, y=smoothed_cases, group=variant)) + 
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA\n(case data NICD & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit case data.png"), width=8, height=6)







# DIDN'T TRY TO RUN / UPDATE THE PART BELOW YET

# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for simplest model fit1_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from+16,
                                                  date.to-extrapolate), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit1_southafrica_multi, ~ variant , at=list(date_num=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_variant)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(fit1_southafrica_multi,  ~ variant , by=c("date_num"),
                                                               var="date_num", mode="latent",
                                                               at=list(date_num=d))
                                                out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                # sum(out$estimate*wt) # should sum to zero
                                                return(out) } ))
above_avg_r_variants = above_avg_r_variants0
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_variant), 
                                       labels=levels_variant)
above_avg_r_variants$variant = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$date = as.Date(above_avg_r_variants$date_num, origin="1970-01-01")
range(above_avg_r_variants$date) # "2021-01-04" "2021-07-30"
above_avg_r_variants$avg_r = avg_r_cases$r[match(above_avg_r_variants$date,
                                                 avg_r_cases$DATE)]  # average growth rate of all lineages calculated from case nrs
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                date_num=avg_r_cases$date_num, # -7 to calculate back to time of infection
                # REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                date=avg_r_cases$DATE,
                variant="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$date_num<=max(above_avg_r_variants$date_num)&df$date_num>=(min(above_avg_r_variants$date_num)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$variant = factor(above_avg_r_variants$variant, levels=c(levels_variant,"avg"))
above_avg_r_variants$prob = fit_southafrica_multi_preds_withCI$prob[match(interaction(above_avg_r_variants$date_num,
                                                                      above_avg_r_variants$variant),
                                                          interaction(fit_southafrica_multi_preds_withCI$date_num,
                                                                      fit_southafrica_multi_preds_withCI$variant))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 3
ymin = 1/3
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$variant %in% c("Other"))|above_avg_r_variants2$date>max(cases_tot$DATE)),], 
      x=date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
      group=variant, linetype=I(0)) +
  # facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  xaxis +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN SOUTH AFRICA\nAT MOMENT OF INFECTION\n(based on case data & multinomial fit to GISAID & SGTF data)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  scale_colour_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)
      
