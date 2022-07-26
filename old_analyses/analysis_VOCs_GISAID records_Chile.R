# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN CHILE (GISAID RECORDS)
# T. Wenseleers
# last update 25 OCTOBER 2021
    
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
  
  
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
# today = as.Date("2021-10-10")
today # 2021-10-25
today_num = as.numeric(today)
plotdir = "Chile_GISAID"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID records for Chile
d1 = read_tsv(file(".//data//GISAID//Chile//gisaid_hcov-19_2021_10_08_21_subm_2020.tsv"), col_types = cols(.default = "c")) 
d2 = read_tsv(file(".//data//GISAID//Chile//gisaid_hcov-19_2021_10_08_21_subm_jan_2021_june_2021.tsv"), col_types = cols(.default = "c")) 
d3 = read_tsv(file(".//data//GISAID//Chile//gisaid_hcov-19_2021_10_25_19_subm_july_oct_2021.tsv"), col_types = cols(.default = "c")) 
GISAID = as.data.frame(rbind(d1, d2, d3))
colnames(GISAID) = c("Virus name","Accession ID","date","Location","host",
                     "Additional location information","Sampling strategy","Gender",                         
                     "Patient age","Patient status","Last vaccinated","Passage","Specimen",
                     "Additional host information","pango_lineage","Clade","AA Substitutions")
date_isvalid = sapply(GISAID$date, function (s) str_count(s, pattern = "-")==2)
GISAID = GISAID[date_isvalid,]
GISAID$date = as.Date(GISAID$date) 
GISAID = GISAID[!is.na(GISAID$date),]
GISAID = GISAID[GISAID$host=="Human",]
GISAID = GISAID[GISAID$date>=as.Date("2021-01-01"),]
range(GISAID$date) # "2021-01-01" "2021-10-16"
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
GISAID = GISAID[GISAID$pango_lineage!="None",]

GISAID$pango_lineage[grepl("B.1.177",GISAID$pango_lineage,fixed=T)] = "B.1.177+"
GISAID$pango_lineage[grepl("B.1.621",GISAID$pango_lineage,fixed=T)] = "B.1.621+"
GISAID$pango_lineage[grepl("B.1.36\\>",GISAID$pango_lineage)] = "B.1.36+"
sum(sum(grepl("AY.34",GISAID$pango_lineage,fixed=T))) # 104
sum(grepl("AY.34",GISAID$pango_lineage,fixed=T)& # 103
      (grepl("Spike_Q677H",GISAID[,"AA Substitutions"]))&
      (grepl("Spike_P681R",GISAID[,"AA Substitutions"])))
GISAID$pango_lineage[grepl("AY.34",GISAID$pango_lineage,fixed=T)&
                       (grepl("Spike_Q677H",GISAID[,"AA Substitutions"]))&
                       (grepl("Spike_P681R",GISAID[,"AA Substitutions"]))
                       ] = "AY.34_with_S:Q677H_and_S:P681R"
GISAID$pango_lineage[grepl("AY.34",GISAID$pango_lineage,fixed=T)&!(
                       (grepl("Spike_Q677H",GISAID[,"AA Substitutions"]))&
                       (grepl("Spike_P681R",GISAID[,"AA Substitutions"])))
] = "AY.34_other"

GISAID[GISAID$pango_lineage=="AY.34_with_S:Q677H_and_S:P681R",]

# GISAID$pango_lineage[grepl("B.1.617.2|AY",GISAID$pango_lineage)] = "Delta (B.1.617.2 & AY.X)"

# sel_target_VOC = "Delta (B.1.617.2 & AY.X)"
sel_target_VOC = "B.1.617.2"

GISAID$LINEAGE = GISAID$pango_lineage
nrow(GISAID) # 10491

# ANALYSIS OF VOCs IN CHILE ####

# sel_countries = "Chile"
# GISAID[GISAID$country %in% sel_countries,]

GISAID_sel = GISAID

sum(GISAID_sel$LINEAGE==sel_target_VOC) # 1197

table(GISAID_sel$LINEAGE)

main_lineages = names(table(GISAID_sel$LINEAGE))[100*table(GISAID_sel$LINEAGE)/sum(table(GISAID_sel$LINEAGE)) > 1]
main_lineages
# "AY.25"     "AY.4", "B.1.1"     "B.1.1.348" "B.1.1.7"   "B.1.617.2" "B.1.621+"  "C.37"      "P.1"       "P.1.12"  "P.1.15"
VOCs = c("B.1.617.1","B.1.617.2","B.1.617+","B.1.618","B.1.1.7","B.1.351","P.1","B.1.1.318","B.1.1.207","B.1.429",
         "B.1.525","B.1.526","B.1.1.519","B.1.1.318","B.1.621+",
         "AY.1","AY.2","AY.3","AY.4","AY.5","AY.6","AY.7","AY.8","AY.9","AY.10","AY.11",
         "AY.12","AY.13","AY.14","AY.15","AY.16","AY.17","AY.18","AY.19","AY.20","AY.21","AY.22",
         "AY.23","AY.24","AY.25","AY.26","AY.27","AY.28","AY.29","AY.30","AY.31","AY.32",
         "AY.33","AY.34",sel_target_VOC,"AY.34_with_S:Q677H_and_S:P681R","AY.34_other")
main_lineages = union(main_lineages, VOCs)
GISAID_sel$LINEAGE[!(GISAID_sel$LINEAGE %in% main_lineages)] = "other" # minority lineages & non-VOCs
remove = unique(c(names(table(GISAID_sel$LINEAGE))[table(GISAID_sel$LINEAGE)/sum(table(GISAID_sel$LINEAGE)) < 0.0005],
                  names(table(GISAID_sel$LINEAGE))[table(GISAID_sel$LINEAGE) < 20]))
remove = remove[!(remove %in% c("B.1.351","B.1.1.7","P.1","B.1.617.2","B.1.617.1","B.1.621+",sel_target_VOC,"B.1.351","P.1","AY.30","AY.4"))]
GISAID_sel$LINEAGE[(GISAID_sel$LINEAGE %in% remove)] = "other" # minority VOCs
table(GISAID_sel$LINEAGE)

GISAID_sel$LINEAGE = factor(GISAID_sel$LINEAGE)
GISAID_sel$LINEAGE = relevel(GISAID_sel$LINEAGE, ref="B.1.1.7") # we code UK strain as the reference level
levels(GISAID_sel$LINEAGE)
# [1] "B.1.1.7"                        "AY.20"                          "AY.25"                          "AY.26"                         
# [5] "AY.3"                           "AY.33"                          "AY.34_with_S:Q677H_and_S:P681R" "AY.4"                          
# [9] "B.1.1"                          "B.1.1.348"                      "B.1.351"                        "B.1.429"                       
# [13] "B.1.617.2"                      "B.1.621+"                       "C.37"                           "other"                         
# [17] "P.1"                            "P.1.12"                         "P.1.15"   

levels_LINEAGE = c("B.1.1.7","B.1.1","B.1.1.348","B.1.351","B.1.429","B.1.621+","C.37","P.1","P.1.12","P.1.15",
                    sel_target_VOC,"AY.3","AY.4","AY.20","AY.25","AY.26","AY.33","AY.34_with_S:Q677H_and_S:P681R","other")
GISAID_sel$LINEAGE = factor(GISAID_sel$LINEAGE, levels=levels_LINEAGE)

# # comparison of mutations in AY.34 & B.1.617.2
# # count of mutations present in AY.34
# t = table(strsplit(paste0(GISAID_sel[GISAID_sel$pango_lineage=="AY.34","AA Substitutions"],collapse=","),",")[[1]])
# t = data.frame(lineage="AY.34", data.frame(sort(t/sum(GISAID_sel$pango_lineage=="AY.34"), decreasing=T)))
# t[t$Freq>0.5,]
# # count of mutations present in B.1.617.2
# t2 = table(strsplit(paste0(GISAID_sel[GISAID_sel$pango_lineage=="B.1.617.2","AA Substitutions"],collapse=","),",")[[1]])
# t2 = data.frame(lineage="B.1.617.2",data.frame(sort(t2/sum(GISAID_sel$pango_lineage=="B.1.617.2"), decreasing=T)))
# t2[t2$Freq>0.5,]
# 
# tab_comparison = merge(t,t2, by="Var1")
# tab_comparison = tab_comparison[order(tab_comparison$Freq.x, decreasing=T),]
# tab_comparison[(tab_comparison$Freq.x>0.5|tab_comparison$Freq.y>0.5)&(abs(tab_comparison$Freq.x-tab_comparison$Freq.y)>0.1),]
# # Spike_T95I, N_D377Y and N_G215C more common in AY.34 than in B.1.617.2
# # Spike_D950N more common in AY.34 than in B.1.617.2
# 
# # Var1            lineage.x     Freq.x lineage.y     Freq.y
# # 3       N_D377Y     AY.34 0.97530864 B.1.617.2 0.64344942
# # 26   NSP4_V167L     AY.34 0.97530864 B.1.617.2 0.84079602
# # 5       N_G215C     AY.34 0.92592593 B.1.617.2 0.05970149
# # 39 Spike_L452R)     AY.34 0.91358025 B.1.617.2 0.67330017
# # 45   Spike_T95I     AY.34 0.87654321 B.1.617.2 0.46766169
# # 33  Spike_D950N     AY.34 0.65432099 B.1.617.2 0.11774461
# # 32  Spike_D950B     AY.34 0.30864198 B.1.617.2 0.85737977
# # 1      (N_G215C     AY.34 0.04938272 B.1.617.2 0.90713101


# aggregated data to make Muller plots of raw data
# aggregate by day to identify days on which INSA performed (days with a lot of sequences)
# we subset the data to just those days to avoid sampling biases (delta infection clusters etc)
data_agbyday = as.data.frame(table(GISAID_sel$date, GISAID_sel$LINEAGE))
colnames(data_agbyday) = c("date", "LINEAGE", "count")
data_agbyday_sum = aggregate(count ~ date, data=data_agbyday, sum)
data_agbyday$total = data_agbyday_sum$count[match(data_agbyday$date, data_agbyday_sum$date)]
sum(data_agbyday[data_agbyday$LINEAGE=="P.1","total"]) == nrow(GISAID_sel) # correct
data_agbyday$date = as.Date(as.character(data_agbyday$date))
data_agbyday$LINEAGE = factor(data_agbyday$LINEAGE, levels=levels_LINEAGE)
data_agbyday$date_num = as.numeric(data_agbyday$date)
data_agbyday$prop = data_agbyday$count/data_agbyday$total
data_agbyday$floor_date = NULL
# qplot(data=data_agbyday, x=date, y=total, colour=total>20, fill=total>20, geom="col")
GISAID_sel$total_sequenced_on_that_day = data_agbyday$total[match(GISAID_sel$date, data_agbyday$date)]
# GISAID_sel = GISAID_sel[GISAID_sel$total_sequenced_on_that_day>20,] # dates on which was performed
# nrow(GISAID_sel) # 5710

# aggregated by week for selected variant lineages
data_agbyweek = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$LINEAGE))
colnames(data_agbyweek) = c("floor_date", "LINEAGE", "count")
data_agbyweek_sum = aggregate(count ~ floor_date, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$floor_date, data_agbyweek_sum$floor_date)]
sum(data_agbyweek[data_agbyweek$LINEAGE=="P.1","total"]) == nrow(GISAID_sel) # correct
data_agbyweek$collection_date = as.Date(as.character(data_agbyweek$floor_date))
data_agbyweek$LINEAGE = factor(data_agbyweek$LINEAGE, levels=levels_LINEAGE)
data_agbyweek$collection_date_num = as.numeric(data_agbyweek$collection_date)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
data_agbyweek$floor_date = NULL


# MULLER PLOT (RAW DATA)
  n2 = length(levels(GISAID_sel$LINEAGE))
  lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.1.7")] = "#0085FF"
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.351")] = "#9A9D00"
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.177+")] = "grey55"
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="P.1")] = "cyan3"
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="B.1.617.1")] = muted("magenta")
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)==sel_target_VOC)] = "magenta"
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.3")] = muted("red",c=150,l=65)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.4")] = muted("red",c=150,l=60)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.20")] = muted("red",c=150,l=55)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.23")] = muted("red",c=150,l=50)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.25")] = muted("red",c=150,l=40)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.26")] = muted("red",c=150,l=30)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.30")] = muted("red",c=150,l=20)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.33")] = muted("red",c=150,l=15)
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.34")] = "blue"
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="AY.34_with_S:Q677H_and_S:P681R")] = "blue"
  lineage_cols2[which(levels(GISAID_sel$LINEAGE)=="other")] = "grey75"

muller_chile_raw2 = ggplot(data=data_agbyweek, aes(x=collection_date, y=count, group=LINEAGE)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN CHILE\n(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_chile_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_muller plots_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_muller plots_raw data.pdf"), width=8, height=6)


# multinomial fits
data_agbyweek$LINEAGE = relevel(data_agbyweek$LINEAGE, ref=sel_target_VOC)
data_agbyweek$DATE_NUM = as.numeric(data_agbyweek$collection_date)

library(nnet)
library(splines)
set.seed(1)
fit1_chile_multi = nnet::multinom(LINEAGE ~ scale(DATE_NUM), weights=count, data=data_agbyweek, maxit=1000)
fit2_chile_multi = nnet::multinom(LINEAGE ~ ns(DATE_NUM, df=2), weights=count, data=data_agbyweek, maxit=1000)
fit3_chile_multi = nnet::multinom(LINEAGE ~ ns(DATE_NUM, df=3), weights=count, data=data_agbyweek, maxit=1000)
BIC(fit1_chile_multi, fit2_chile_multi, fit3_chile_multi) 
# df      BIC
# fit1_chile_multi 36 38052.99
# fit2_chile_multi 54 36840.15
# fit3_chile_multi 72 36883.20

# growth rate advantage compared to Delta (difference in growth rate per day) 
emtrchile = emtrends(fit2_chile_multi, trt.vs.ctrl ~ LINEAGE,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(GISAID_sel$DATE_NUM)),
                   adjust="none", df=NA)
delta_r_chile = data.frame(confint(emtrchile, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrchile$contrasts,
                                               adjust="none", df=NA)$p.value)
delta_r_chile
#                                      contrast      estimate          SE df   asymp.LCL     asymp.UCL       p.value
# 1                         B.1.1.7 - B.1.617.2 -0.1149014620 0.005203140 NA -0.12509943 -0.1047034955 4.594517e-108
# 2                           B.1.1 - B.1.617.2 -0.0991964038 0.004944137 NA -0.10888673 -0.0895060734  1.540436e-89
# 3                       B.1.1.348 - B.1.617.2 -0.1212460621 0.006965389 NA -0.13489797 -0.1075941503  7.309625e-68
# 4                         B.1.351 - B.1.617.2 -0.0912797202 0.019270451 NA -0.12904911 -0.0535103297  2.171500e-06
# 5                         B.1.429 - B.1.617.2 -0.1144340213 0.012132941 NA -0.13821415 -0.0906538942  4.035705e-21
# 6                      (B.1.621+) - B.1.617.2 -0.0673657791 0.004092144 NA -0.07538623 -0.0593453249  6.853717e-61
# 7                            C.37 - B.1.617.2 -0.1145674796 0.003309171 NA -0.12105334 -0.1080816239 1.212047e-262
# 8                             P.1 - B.1.617.2 -0.1125397558 0.003191180 NA -0.11879435 -0.1062851577 1.959531e-272
# 9                          P.1.12 - B.1.617.2 -0.1382642218 0.008029469 NA -0.15400169 -0.1225267523  1.892992e-66
# 10                         P.1.15 - B.1.617.2 -0.1296083123 0.003494270 NA -0.13645696 -0.1227596684 3.827856e-301
# 11                           AY.3 - B.1.617.2 -0.0911844738 0.017658760 NA -0.12579501 -0.0565739397  2.421205e-07
# 12                           AY.4 - B.1.617.2  0.0005515706 0.006473467 NA -0.01213619  0.0132393319  9.320986e-01
# 13                          AY.20 - B.1.617.2 -0.0709979029 0.016550168 NA -0.10343564 -0.0385601689  1.787859e-05
# 14                          AY.25 - B.1.617.2 -0.0069549706 0.003143875 NA -0.01311685 -0.0007930891  2.695086e-02
# 15                          AY.26 - B.1.617.2 -0.0421809408 0.007818553 NA -0.05750502 -0.0268568590  6.853094e-08
# 16                          AY.33 - B.1.617.2 -0.0501971017 0.017140330 NA -0.08379153 -0.0166026723  3.404967e-03
# 17 AY.34_with_S:Q677H_and_S:P681R - B.1.617.2  0.0025444564 0.010386123 NA -0.01781197  0.0229008830  8.064671e-01
# 18                          other - B.1.617.2 -0.0678492413 0.003248796 NA -0.07421676 -0.0614817190  7.418564e-97
  
# AY.34 1% more infectious than B.1.617.2 with GT of 4.7 days
exp(0.0025444564*4.7) # = 1

# fitted prop of different LINEAGES in the chile today
multinom_preds_today_avg = data.frame(emmeans(fit2_chile_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg



# PLOT MULTINOMIAL FIT

extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
date.to = today_num+extrapolate

# multinomial model predictions (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to)))
fit_chile_multi_preds = data.frame(predgrid, as.data.frame(predict(fit2_chile_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fit_chile_multi_preds = gather(fit_chile_multi_preds, LINEAGE, prob, all_of(levels_LINEAGE), factor_key=TRUE)
fit_chile_multi_preds$collection_date = as.Date(fit_chile_multi_preds$DATE_NUM, origin="1970-01-01")
fit_chile_multi_preds$LINEAGE = factor(fit_chile_multi_preds$LINEAGE, levels=levels_LINEAGE) 

muller_chile_mfit = ggplot(data=fit_chile_multi_preds, 
                                   aes(x=collection_date, y=prob, group=LINEAGE)) + 
  # facet_wrap(~ STATE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN CHILE\n(GISAID data, multinomial fit)")
muller_chile_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_muller plots_multinom fit.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_muller plots_multinom fit.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_chile_raw2 + coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_chile_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_muller plots multipanel_multinom fit.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_muller plots multipanel_multinom fit.pdf"), width=10, height=10)





# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions with confidence intervals (but slower)
fit_chile_multi_preds_withCI = data.frame(emmeans(fit2_chile_multi,
                                                        ~ LINEAGE,
                                                        by=c("DATE_NUM"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=2)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_chile_multi_preds_withCI$collection_date = as.Date(fit_chile_multi_preds_withCI$DATE_NUM, origin="1970-01-01")
fit_chile_multi_preds_withCI$LINEAGE = factor(fit_chile_multi_preds_withCI$LINEAGE, levels=levels_LINEAGE)
fit_chile_multi_preds2 = fit_chile_multi_preds_withCI


# on logit scale:

ymin = 0.001
ymax = 0.999
fit_chile_multi_preds2$asymp.LCL[fit_chile_multi_preds2$asymp.LCL<ymin] = ymin
fit_chile_multi_preds2$asymp.UCL[fit_chile_multi_preds2$asymp.UCL<ymin] = ymin
fit_chile_multi_preds2$asymp.UCL[fit_chile_multi_preds2$asymp.UCL>ymax] = ymax
fit_chile_multi_preds2$prob[fit_chile_multi_preds2$prob<ymin] = ymin

plot_chile_mfit_logit = qplot(data=fit_chile_multi_preds2, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN CHILE\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="identity",
                        range=c(0.1, 4), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_chile_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_multinom fit_logit scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_multinom fit_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_chile_mfit = qplot(data=fit_chile_multi_preds2, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN CHILE\n(GISAID data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyweek,
             aes(x=collection_date, y=100*prop, size=total,
                 colour=LINEAGE
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="identity",
                        range=c(0.5, 5), limits=c(1,max(data_agbyweek$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_chile_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_multinom fit_response scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\chile_multinom fit_response scale.pdf"), width=10, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot = as.data.frame(get_national_data(countries = "Chile"))
cases_tot = cases_tot[cases_tot$date>=as.Date("2020-08-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
# cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
range(cases_tot$date)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
library(mgcv)
k=25
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F) + 
                  WEEKDAY, # + 
                # BANKHOLIDAY,
                # s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=poisson(log), data=cases_tot,
                method = "REML",
                knots = list(DATE_NUM = c(min(cases_tot$DATE_NUM)-14,
                                          seq(min(cases_tot$DATE_NUM)+0.7*diff(range(cases_tot$DATE_NUM))/(k-2), 
                                              max(cases_tot$DATE_NUM)-0.7*diff(range(cases_tot$DATE_NUM))/(k-2), length.out=k-2),
                                          max(cases_tot$DATE_NUM)+14))
) 
BIC(fit_cases)

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_chile_multi_preds_withCI$totcases = cases_tot$cases_new[match(round(fit_chile_multi_preds_withCI$DATE_NUM),cases_tot$DATE_NUM)]
fit_chile_multi_preds_withCI$cases = fit_chile_multi_preds_withCI$totcases * fit_chile_multi_preds_withCI$prob
fit_chile_multi_preds_withCI$cases[fit_chile_multi_preds_withCI$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM, at=list(DATE_NUM=seq(date.from, date.to, by=0.5), BANHOLIDAY="no"), type="response"))
fit_chile_multi_preds_withCI$smoothed_totcases = cases_emmeans$rate[match(fit_chile_multi_preds_withCI$DATE_NUM,cases_emmeans$DATE_NUM)]
fit_chile_multi_preds_withCI$smoothed_cases = fit_chile_multi_preds_withCI$smoothed_totcases * fit_chile_multi_preds_withCI$prob
fit_chile_multi_preds_withCI$smoothed_cases[fit_chile_multi_preds_withCI$smoothed_cases<=0.001] = NA
fit_chile_multi_preds_withCI$LINEAGE = factor(fit_chile_multi_preds_withCI$LINEAGE, levels=levels_LINEAGE)

ggplot(data=fit_chile_multi_preds_withCI[fit_chile_multi_preds_withCI$collection_date>=as.Date("2021-01-01"),], 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),max(cases_tot$date)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN CHILE\n(case data & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_chile_multi_preds_withCI[fit_chile_multi_preds_withCI$collection_date>=as.Date("2021-01-01")&
                                              fit_chile_multi_preds_withCI$collection_date<=max(cases_tot$date),], 
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01","2021-11-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),today), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN CHILE\n(case data & multinomial fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit case data.png"), width=8, height=6)



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
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to)#,
                                             # BANKHOLIDAY="no"
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[2] = "r"
colnames(avg_r_cases)[5] = "r_LOWER"
colnames(avg_r_cases)[6] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$DATE_NUM, origin="1970-01-01") # -7 TO CALCULATE BACK TO INFECTION DATE
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN CHILE AT MOMENT OF INFECTION BASED ON NEW CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit2_chile_multi, ~ LINEAGE , at=list(DATE_NUM=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(fit2_chile_multi,  ~ LINEAGE , by=c("DATE_NUM"),
                                                               var="DATE_NUM", mode="latent",
                                                               at=list(DATE_NUM=d))
                                                out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                # sum(out$estimate*wt) # should sum to zero
                                                return(out) } ))
above_avg_r_variants = above_avg_r_variants0
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE), 
                                       labels=levels_LINEAGE)
above_avg_r_variants$variant = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2021-01-04" "2021-07-30"
above_avg_r_variants$avg_r = avg_r_cases$r[match(above_avg_r_variants$collection_date,
                                                 avg_r_cases$DATE)]  # average growth rate of all lineages calculated from case nrs
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                # REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE,
                variant="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$variant = factor(above_avg_r_variants$variant, levels=c(levels_LINEAGE,"avg"))
above_avg_r_variants$prob = fit_chile_multi_preds_withCI$prob[match(interaction(above_avg_r_variants$DATE_NUM,
                                                                      above_avg_r_variants$variant),
                                                          interaction(fit_chile_multi_preds_withCI$DATE_NUM,
                                                                      fit_chile_multi_preds_withCI$LINEAGE))]
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
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$variant %in% c("other"))|above_avg_r_variants2$collection_date>max(cases_tot$DATE)),], 
      x=collection_date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
      group=variant, linetype=I(0)) +
  # facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01","2021-08-01","2021-09-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN CHILE\nAT MOMENT OF INFECTION\n(based on case data & multinomial fit to GISAID data)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  scale_colour_manual("variant", values=c(head(lineage_cols2,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)
      
