# ANALYSIS OF GROWTH ADVANTAGE OF OMICRON (B.1.1.529) IN SOUTH AFRICA BASED ON GISAID SEQUENCE DATA

# T. Wenseleers
# last update 30 JANUARY 2022
    
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
  # today = as.Date("2021-12-14")
  today_num = as.numeric(today)
  tag = paste("@TWenseleers\n",today)
  plotdir = "South Africa_GISAID"
  suppressWarnings(dir.create(paste0(".//plots//",plotdir)))
  levels_PROVINCES = c("Gauteng", "North West", "Mpumalanga", "Limpopo", "Western Cape", "Free State", "KwaZulu Natal", "Eastern Cape", "Northern Cape")
  sel_target_VOC = "Omicron (BA.2)"
  sel_reference_VOC = "Omicron (BA.1)"
  levels_VARIANTS = c(sel_reference_VOC, "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1.1)", "Omicron (BA.3)", sel_target_VOC)
  n = length(levels_VARIANTS)
  lineage_cols = hcl(h = seq(0, 180, length = n), l = 60, c = 180)
  lineage_cols[which(levels_VARIANTS=="Alpha")] = "#0085FF"
  lineage_cols[which(levels_VARIANTS=="Beta")] = "green4"
  lineage_cols[which(levels_VARIANTS=="Delta")] = "mediumorchid"
  # lineage_cols[which(levels_VARIANTS=="C.1.2")] = "darkorange"
  lineage_cols[which(levels_VARIANTS==sel_target_VOC)] = "red2" # "magenta"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.3)")] = "darkorange" 
  lineage_cols[which(levels_VARIANTS=="Other")] = "grey65"
  
  # X axis for plots
  firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
  xaxis = scale_x_continuous(breaks=firststofmonth,
                             labels=substring(months(firststofmonth),1,1),
                             expand=c(0,0))
  
  
  # import GISAID records for South Africa
d1 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_01_30_19_subm_2020.tsv"), col_types = cols(.default = "c")) 
d2 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_01_30_19_subm_jan_may_2021.tsv"), col_types = cols(.default = "c")) 
d3 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_01_30_19_subm_june_aug_2021.tsv"), col_types = cols(.default = "c")) 
d4 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_01_30_19_subm_sept_2021_dec_2021.tsv"), col_types = cols(.default = "c")) 
d5 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_01_30_19_subm_jan_2022.tsv"), col_types = cols(.default = "c")) 

# parse GISAID & SGTF data
GISAID = as.data.frame(rbind(d1,d2,d3,d4,d5))
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
range(GISAID$date) # "2020-03-06" "2022-01-18"
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Month = lubridate::month(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
GISAID$Year_Month = interaction(GISAID$Year,GISAID$Month)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$date_num = as.numeric(GISAID$date)
# GISAID = GISAID[GISAID$pango_lineage!="None",]

# parse location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID$Location), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
unique(loc[,1]) # continent
unique(loc[,2]) # country
unique(loc[,3]) # city or province
table(GISAID[GISAID$date>as.Date("2021-10-01"),"pango_lineage"])
GISAID$country = loc[,2]
GISAID$province = factor(toupper(loc[,3]))
GISAID$province = factor(GISAID$province, levels=levels(GISAID$province),
                                          labels=c("Western Cape", "Eastern Cape", "Free State",
                                                   "Gauteng", "KwaZulu Natal", "KwaZulu Natal", "KwaZulu Natal",
                                                   "Limpopo", "Mpumalanga", "North West",
                                                   "Northern Cape", "Northern Cape", NA,
                                                   "Western Cape", "Western Cape"))
GISAID$province = factor(GISAID$province, levels=levels_PROVINCES)

GISAID$variant = case_when(
  # grepl("N679K", aa_substitutions) & grepl("H655Y", aa_substitutions) & grepl("P681H", aa_substitutions) ~ "B.1.1.529",
  (GISAID$pango_lineage=="BA.1") ~ "Omicron (BA.1)",
  (GISAID$pango_lineage=="BA.1.1") ~ "Omicron (BA.1.1)",
  (GISAID$pango_lineage=="BA.2") ~ "Omicron (BA.2)",
  (GISAID$pango_lineage=="BA.3") ~ "Omicron (BA.3)",
  grepl("B.1.617.2", GISAID$pango_lineage, fixed=T) | grepl("AY", GISAID$pango_lineage)  ~ "Delta",
  grepl("B.1.1.7", GISAID$pango_lineage, fixed=T) ~ "Alpha",
  grepl("B.1.351", GISAID$pango_lineage, fixed=T) ~ "Beta",
  # grepl("C.1.2", GISAID$pango_lineage, fixed=T) ~ "C.1.2",
  T ~ "Other"
)

GISAID[GISAID$variant == "Omicron" & GISAID$date<as.Date("2021-10-01"),"variant"] = "Other" # miscoded date??

GISAID$sampling_strategy = case_when(grepl("Pneumonia", GISAID$`Sampling strategy`) ~ "Pneumonia Surveillance",
                                     grepl("Baseline|Community|Random|routine|Routine|Sentinal|Non-sentinel|Sentinel|ILI|surveillance", GISAID$`Sampling strategy`) ~ "Baseline Surveillance",
                                     grepl("Clinic|CLINIC|Med unit|CASUALTY|CASAULTY|Casualty|Emergency|CLINIC|Hospitalised|Postmortem", GISAID$`Sampling strategy`) ~ "Hospitalised or died",
                                     T ~ "Other"
)

GISAID$age = as.numeric(GISAID$`Patient age`)
GISAID$age[GISAID$age>110] = NA
range(GISAID$age, na.rm=T)
table(GISAID$variant)
table(GISAID$sampling_strategy)
table(GISAID$variant, GISAID$sampling_strategy)
table(GISAID$variant, GISAID$sampling_strategy)/rowSums(table(GISAID$variant, GISAID$sampling_strategy))

GISAID = GISAID[!(GISAID$variant=="Other"&GISAID$pango_lineage=="None"),]
GISAID = GISAID[!is.na(GISAID$province),]

table(GISAID$variant)

# ANALYSIS OF VOCs IN SOUTH AFRICA ####
GISAID_sel = GISAID
nrow(GISAID_sel) # 27360
sum(GISAID_sel$variant==sel_target_VOC) # 779
table(GISAID_sel$variant)
range(GISAID_sel$date)

GISAID_sel$variant = factor(GISAID_sel$variant)
GISAID_sel$variant = relevel(GISAID_sel$variant, ref=sel_reference_VOC) # we code Omicron (BA.1) as the reference

GISAID_sel$variant = factor(GISAID_sel$variant, levels=levels_VARIANTS)
table(GISAID_sel$variant, GISAID_sel$sampling_strategy)

# age distribution of patients infected with Delta or Omicron (GISAID data)
library(ggplot2)
library(ggthemes)
GISAID_sel2 = GISAID_sel[(GISAID_sel$variant %in% c("Delta", "Omicron (BA.1)")) & (GISAID_sel$floor_date>as.Date("2021-05-01")),]
GISAID_sel2$variant = factor(GISAID_sel2$variant, levels=c("Omicron (BA.1)", "Delta"))
GISAID_sel2$Year_Month = factor(GISAID_sel2$Year_Month)
GISAID_sel2$Year_Month = droplevels(GISAID_sel2$Year_Month)
GISAID_sel3 = GISAID_sel
GISAID_sel3 = GISAID_sel3[(GISAID_sel3$variant %in% c("Other","Beta","Delta","Omicron (BA.1)")),]
GISAID_sel3$variant = factor(GISAID_sel3$variant, levels=c("Other","Beta","Delta","Omicron (BA.1)"))
ggplot(data=GISAID_sel3, aes(x=age, fill=variant)) + 
  facet_wrap(~ variant, scale="free_y", nrow=6, drop=FALSE) + geom_histogram() +
  scale_fill_manual(values=c(lineage_cols[5],lineage_cols[2],lineage_cols[6],lineage_cols[1])) + 
  theme(legend.position="none") + 
  ggtitle("Age distribution of individuals with sequenced\ninfections of different variants in South Africa","(GISAID data)") +
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(fill="white", colour="white",size=1)) +
  theme(strip.text.x = element_text(margin = margin(0, 0, 0, 0, "cm")))
ggsave(file=paste0(".\\plots\\",plotdir,"\\age distribution GISAID by variant.png"), width=5, height=7)

ggplot(data=GISAID_sel2, aes(x=age, fill=variant)) + 
  facet_wrap(~ Year_Month+variant, scale="free_y", ncol=2,drop=FALSE) + geom_histogram() +
  scale_fill_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1)","Delta"))]) + theme(legend.position="none") + 
  ggtitle("Age distribution of individuals with sequenced\nOmicron & Delta infections in South Africa","(GISAID data)") +
  theme(strip.text = element_text(size=8),
        strip.background = element_rect(fill="white", colour="white",size=1)) +
  theme(strip.text.x = element_text(margin = margin(0, 0, 0, 0, "cm")))
ggsave(file=paste0(".\\plots\\",plotdir,"\\age distribution GISAID delta omicron by week.png"), width=5, height=7)

ggplot(data=GISAID_sel2, aes(x=date, y=age, fill=variant, colour=variant, group=variant)) + 
  # facet_wrap(~ Year_Month+variant, scale="free_y", ncol=2,drop=FALSE) + 
  geom_point(aes(alpha=I(0.1)), pch=I(16)) +
  geom_smooth(method="gam") +
  scale_fill_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1)","Delta"))]) + 
  scale_colour_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1)","Delta"))]) + 
  theme_hc() +
  xlab("") +
  theme(legend.position="none") + 
  ggtitle("Age distribution of individuals with sequenced\nOmicron & Delta infections in South Africa","(GISAID data)") +
  theme(strip.text = element_text(size=8),
        strip.background = element_rect(fill="white", colour="white",size=1)) +
  theme(strip.text.x = element_text(margin = margin(0, 0, 0, 0, "cm"))) +
  xaxis 

# proportion of sequenced cases over 60
GISAID_sel2$over60 = (GISAID_sel2$age > 60)*1
ggplot(data=GISAID_sel2, aes(x=date, y=over60, fill=variant, colour=variant, group=variant)) + 
  # facet_wrap(~ Year_Month+variant, scale="free_y", ncol=2,drop=FALSE) + 
  # geom_point(aes(alpha=I(0.1)), pch=I(16)) +
  geom_smooth(method="glm", 
              method.args=list(family="binomial"),
              formula = y ~ ns(x, df=2)) +
  scale_fill_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1)","Delta"))]) + 
  scale_colour_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1)","Delta"))]) + 
  ylab("Proportion over 60") +
  xlab("") +
  theme_hc() +
  ggtitle("Proportion of individuals with sequenced\nOmicron & Delta infections that are\nover 60 in South Africa","(GISAID data, logistic 2 df spline fit)") +
  theme(strip.text = element_text(size=8),
        strip.background = element_rect(fill="white", colour="white",size=1)) +
  theme(strip.text.x = element_text(margin = margin(0, 0, 0, 0, "cm"))) +
  xaxis +
  theme(legend.position="none") 
ggsave(file=paste0(".\\plots\\",plotdir,"\\age distribution GISAID delta omicron proportion over 60.png"), width=5, height=7)


# aggregated data to make Muller plots of raw data

# aggregated for the whole of South Africa combined
# aggregate by day 
data_agbyday = as.data.frame(table(GISAID_sel$date, GISAID_sel$variant))
colnames(data_agbyday) = c("date", "variant", "count")
data_agbyday_sum = aggregate(count ~ date, data=data_agbyday, sum)
data_agbyday$total = data_agbyday_sum$count[match(data_agbyday$date, data_agbyday_sum$date)]
data_agbyday$date = as.Date(as.character(data_agbyday$date))
# data_agbyday = rbind(data_agbyday, data_nov) # merge with SGTF+extrapolated GISAID for november
data_agbyday$variant = factor(data_agbyday$variant, levels=levels_VARIANTS)
data_agbyday$date_num = as.numeric(data_agbyday$date)
data_agbyday$prop = data_agbyday$count/data_agbyday$total
data_agbyday$floor_date = NULL

# aggregated by week
data_agbyweek = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$variant))
colnames(data_agbyweek) = c("floor_date", "variant", "count")
data_agbyweek_sum = aggregate(count ~ floor_date, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$floor_date, data_agbyweek_sum$floor_date)]
sum(data_agbyweek[data_agbyweek$variant=="Beta","total"]) == nrow(GISAID_sel) # TRUE
data_agbyweek$date = as.Date(as.character(data_agbyweek$floor_date))
data_agbyweek$variant = factor(data_agbyweek$variant, levels=levels_VARIANTS)
data_agbyweek$date_num = as.numeric(data_agbyweek$date)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
data_agbyweek$floor_date = NULL
data_agbyweek$date_num = as.numeric(data_agbyweek$date)

# aggregated data by province
# aggregate by day 
data_agbyday_province = as.data.frame(table(GISAID_sel$date, GISAID_sel$variant, GISAID_sel$province))
colnames(data_agbyday_province) = c("date", "variant", "province", "count")
data_agbyday_province_sum = aggregate(count ~ date+province, data=data_agbyday_province, sum)
data_agbyday_province$total = data_agbyday_province_sum$count[match(interaction(data_agbyday_province$date, data_agbyday_province$province),
                                                  interaction(data_agbyday_province_sum$date, data_agbyday_province_sum$province))]
data_agbyday_province$date = as.Date(as.character(data_agbyday_province$date))
# data_agbyday = rbind(data_agbyday, data_nov) # merge with SGTF+extrapolated GISAID for november
data_agbyday_province$variant = factor(data_agbyday_province$variant, levels=levels_VARIANTS)
data_agbyday_province$date_num = as.numeric(data_agbyday_province$date)
data_agbyday_province$prop = data_agbyday_province$count/data_agbyday_province$total
data_agbyday_province$floor_date = NULL

# aggregated by week
data_agbyweek_province = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$variant, GISAID_sel$province))
colnames(data_agbyweek_province) = c("floor_date", "variant", "province", "count")
data_agbyweek_province_sum = aggregate(count ~ floor_date, data=data_agbyweek_province, sum)
data_agbyweek_province$total = data_agbyweek_province_sum$count[match(interaction(data_agbyweek_province$floor_date, data_agbyweek_province$province),
                                                                      interaction(data_agbyweek_province_sum$floor_date, data_agbyweek_province_sum$province))]
data_agbyweek_province$date = as.Date(as.character(data_agbyweek_province$floor_date))
data_agbyweek_province$variant = factor(data_agbyweek_province$variant, levels=levels_VARIANTS)
data_agbyweek_province$date_num = as.numeric(data_agbyweek_province$date)
data_agbyweek_province$prop = data_agbyweek_province$count/data_agbyweek_province$total
data_agbyweek_province$floor_date = NULL
data_agbyweek_province$date_num = as.numeric(data_agbyweek_province$date)


# MULLER PLOT (RAW DATA)
# for the whole of South Africa combined
muller_southafrica_raw = ggplot(data=data_agbyweek, aes(x=date, y=count, group=variant)) + 
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
  xaxis +
  scale_y_continuous(expand=c(0,0)) +
  theme_hc() +
  # labs(title = "SARS-CoV2 VARIANTS IN SOUTH AFRICA") +
  ylab("Share") + 
  theme_hc() +
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA","(GISAID data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_sel$Week)))
muller_southafrica_raw

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_raw data.pdf"), width=8, height=6)

# by province
muller_southafrica_byprovince_raw = ggplot(data=data_agbyweek_province, aes(x=date, y=count, group=variant)) + 
  facet_wrap(~ province) +
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill", width=I(7)) + # or geom_area
  scale_fill_manual("", values=lineage_cols) +
  xaxis +
  scale_y_continuous(expand=c(0,0)) +
  theme_hc() +
  ylab("Share") + 
  theme_hc() +
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA","(GISAID data)") 
muller_southafrica_byprovince_raw

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_byprovince_raw data.png"), width=10, height=6)



# 2. ESTIMATE OF GROWTH RATE & TRANSMISSION ADVANTAGE OF OMICRON OVER DELTA BASED ON GISAID SEQUENCING DATA ####

# multinomial fits per province ####
library(nnet)
library(splines)
set.seed(1)
# because dataset is small I do not use splines here & only add province as a main effect (allowing for different dates of introduction, but not
# differences in competitive advantage, as more data comes in this could be changed)
# one could also add age in the model if desired
fit1_southafrica_province_multi = nnet::multinom(variant ~ scale(date_num)+province, weights=count, data=data_agbyday_province, maxit=1000)
# fit2_southafrica_province_multi = nnet::multinom(variant ~ ns(date_num, df=2)+province, weights=count, data=data_agbyday_province, maxit=1000)
BIC(fit1_southafrica_province_multi) 
#                                  df      BIC
# fit1_southafrica_province_multi  50 22891.81

# and also a multinomial 2 df spline fit to the data aggregated over all provinces
fit1_southafrica_multi = nnet::multinom(variant ~ scale(date_num), weights=count, data=data_agbyday, maxit=1000)
fit2_southafrica_multi = nnet::multinom(variant ~ ns(date_num, df=2), weights=count, data=data_agbyday, maxit=1000)
BIC(fit1_southafrica_multi, fit2_southafrica_multi) 
# fit2 best

# avg growth rate advantage of Omicron (BA.2) over Omicron BA.1 (difference in growth rate per day) (today)
# based on multinomial model fit2_southafrica_multi without province included as factor
emtrsouthafrica_avg = emtrends(fit2_southafrica_multi, trt.vs.ctrl ~ variant,  
                                        var="date_num",  mode="latent",
                                        at=list(date_num=seq(today_num, today_num, by=1)),
                                        adjust="none", df=NA)
delta_r_southafrica_avg = data.frame(confint(emtrsouthafrica_avg, 
                                                      adjust="none", df=NA)$contrasts, 
                                              p.value=as.data.frame(emtrsouthafrica_avg$contrasts,
                                                                    adjust="none", df=NA)$p.value)
delta_r_southafrica_avg
# contrast   estimate          SE df   asymp.LCL    asymp.UCL       p.value
# 1             Beta - Omicron (BA.1) -0.2370734 0.008608625 NA -0.25394595 -0.220200757 5.985722e-167
# 2            Alpha - Omicron (BA.1) -0.2304367 0.009849771 NA -0.24974187 -0.211131473 4.790508e-121
# 3            Other - Omicron (BA.1) -0.1741009 0.008416516 NA -0.19059700 -0.157604868  4.666373e-95
# 4            Delta - Omicron (BA.1) -0.2142903 0.008404795 NA -0.23076337 -0.197817176 2.172567e-143
# 5 Omicron (BA.1.1) - Omicron (BA.1)  0.0151324 0.004158423 NA  0.00698204  0.023282759  2.737251e-04
# 6   Omicron (BA.3) - Omicron (BA.1) -0.0184686 0.012208664 NA -0.04239714  0.005459937  1.303443e-01
# 7   Omicron (BA.2) - Omicron (BA.1)  0.1035992 0.006256452 NA  0.09133683  0.115861671  1.383597e-61

# corresponding transmission advantage of Omicron BA.2 over BA.1 (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_avg[7,5]*2.2) # 1.22x
exp(delta_r_southafrica_avg[7,2]*2.2) # Omicron 1.25x transmission advantage over Delta
exp(delta_r_southafrica_avg[7,6]*2.2) # 1.29x


# avg growth rate advantage of Omicron BA.2 over BA.1 (difference in growth rate per day) (today)
# based on multinomial model with province included as factor
emtrsouthafrica_province_avg = emtrends(fit1_southafrica_province_multi, trt.vs.ctrl ~ variant,  
                   var="date_num",  mode="latent",
                   at=list(date_num=seq(today_num, today_num, by=1)),
                   adjust="none", df=NA)
delta_r_southafrica_province_avg = data.frame(confint(emtrsouthafrica_province_avg, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrsouthafrica_province_avg$contrasts,
                                               adjust="none", df=NA)$p.value)
delta_r_southafrica_province_avg
# contrast     estimate          SE df   asymp.LCL   asymp.UCL       p.value
# 1             Beta - Omicron (BA.1) -0.276510996 0.011553442 NA -0.29915533 -0.25386667 1.382057e-126
# 2            Alpha - Omicron (BA.1) -0.265390172 0.011587467 NA -0.28810119 -0.24267915 4.316479e-116
# 3            Other - Omicron (BA.1) -0.287322550 0.011556050 NA -0.30997199 -0.26467311 1.852902e-136
# 4            Delta - Omicron (BA.1) -0.233471703 0.011527623 NA -0.25606543 -0.21087798  3.326322e-91
# 5 Omicron (BA.1.1) - Omicron (BA.1)  0.024519103 0.004442325 NA  0.01581231  0.03322590  3.401012e-08
# 6   Omicron (BA.3) - Omicron (BA.1) -0.005913976 0.015280129 NA -0.03586248  0.02403453  6.987288e-01
# 7   Omicron (BA.2) - Omicron (BA.1)  0.112914554 0.006752651 NA  0.09967960  0.12614951  9.136586e-63

# corresponding transmission advantage of Omicron over Delta (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_province_avg[7,5]*2.2) # 1.24x
exp(delta_r_southafrica_province_avg[7,2]*2.2) # Omicron 1.28x transmission advantage over Delta
exp(delta_r_southafrica_province_avg[7,6]*2.2) # 1.32x

# fitted prop of different variantS today (on average across provinces)
# based on model with province included as a factor
multinom_preds_today_avg_province = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|1,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg_province
# variant         prob           SE df     asymp.LCL    asymp.UCL
# 1   Omicron (BA.1) 2.952433e-01 4.283087e-02 NA  2.112963e-01 3.791903e-01
# 2             Beta 2.771189e-13 2.256821e-13 NA -1.652099e-13 7.194477e-13
# 3            Alpha 1.157152e-13 9.818370e-14 NA -7.672128e-14 3.081517e-13
# 4            Other 2.301611e-15 1.908902e-15 NA -1.439769e-15 6.042991e-15
# 5            Delta 1.044378e-08 8.316977e-09 NA -5.857201e-09 2.674475e-08
# 6 Omicron (BA.1.1) 6.357452e-02 1.920182e-02 NA  2.593965e-02 1.012094e-01
# 7   Omicron (BA.3) 3.224436e-03 3.045244e-03 NA -2.744134e-03 9.193005e-03
# 8   Omicron (BA.2) 6.379577e-01 5.409274e-02 NA  5.319379e-01 7.439776e-01

# fitted prop of different variantS today (by province)
multinom_preds_today_byprovince = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|province,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_byprovince
# variant      province         prob           SE df     asymp.LCL    asymp.UCL
# 1    Omicron (BA.1)       Gauteng 7.792932e-02 1.876653e-02 NA  4.114759e-02 1.147111e-01
# 2              Beta       Gauteng 6.210921e-15 5.834212e-15 NA -5.223924e-15 1.764577e-14
# 3             Alpha       Gauteng 1.095224e-14 1.060120e-14 NA -9.825725e-15 3.173020e-14
# 4             Other       Gauteng 4.164347e-17 3.936772e-17 NA -3.551584e-17 1.188028e-16
# 5             Delta       Gauteng 2.254017e-10 2.078507e-10 NA -1.819783e-10 6.327816e-10
# 6  Omicron (BA.1.1)       Gauteng 3.249137e-02 1.052442e-02 NA  1.186388e-02 5.311886e-02
# 7    Omicron (BA.3)       Gauteng 4.366979e-04 4.486639e-04 NA -4.426671e-04 1.316063e-03
# 8    Omicron (BA.2)       Gauteng 8.891426e-01 2.636545e-02 NA  8.374673e-01 9.408179e-01
# 9    Omicron (BA.1)    North West 4.905155e-01 2.072202e-01 NA  8.437134e-02 8.966596e-01
# 10             Beta    North West 1.880803e-13 1.824340e-13 NA -1.694838e-13 5.456445e-13
# 11            Alpha    North West 1.541452e-13 1.611107e-13 NA -1.616261e-13 4.699164e-13
# 12            Other    North West 5.728380e-15 5.601371e-15 NA -5.250105e-15 1.670687e-14
# 13            Delta    North West 9.632598e-09 9.157465e-09 NA -8.315704e-09 2.758090e-08
# 14 Omicron (BA.1.1)    North West 9.907717e-02 5.196057e-02 NA -2.763678e-03 2.009180e-01
# 15   Omicron (BA.3)    North West 6.094574e-03 6.417469e-03 NA -6.483434e-03 1.867258e-02
# 16   Omicron (BA.2)    North West 4.043128e-01 2.498592e-01 NA -8.540221e-02 8.940277e-01
# 17   Omicron (BA.1)    Mpumalanga 3.599685e-01 1.793178e-01 NA  8.512118e-03 7.114249e-01
# 18             Beta    Mpumalanga 3.359807e-13 3.419806e-13 NA -3.342890e-13 1.006250e-12
# 19            Alpha    Mpumalanga 5.023548e-14 5.955161e-14 NA -6.648353e-14 1.669545e-13
# 20            Other    Mpumalanga 3.446643e-15 3.536311e-15 NA -3.484399e-15 1.037768e-14
# 21            Delta    Mpumalanga 4.688939e-09 4.690247e-09 NA -4.503777e-09 1.388165e-08
# 22 Omicron (BA.1.1)    Mpumalanga 1.541016e-01 8.520044e-02 NA -1.288820e-02 3.210914e-01
# 23   Omicron (BA.3)    Mpumalanga 1.991368e-02 2.033320e-02 NA -1.993867e-02 5.976602e-02
# 24   Omicron (BA.2)    Mpumalanga 4.660162e-01 2.631177e-01 NA -4.968494e-02 9.817173e-01
# 25   Omicron (BA.1)       Limpopo 3.536005e-02 1.938458e-02 NA -2.633033e-03 7.335314e-02
# 26             Beta       Limpopo 6.284509e-17 2.295934e-16 NA -3.871496e-16 5.128398e-16
# 27            Alpha       Limpopo 3.234481e-17 1.191666e-16 NA -2.012174e-16 2.659070e-16
# 28            Other       Limpopo 1.129070e-18 4.129294e-18 NA -6.964198e-18 9.222337e-18
# 29            Delta       Limpopo 4.229781e-12 1.543039e-11 NA -2.601323e-11 3.447279e-11
# 30 Omicron (BA.1.1)       Limpopo 1.757606e-02 1.135447e-02 NA -4.678301e-03 3.983042e-02
# 31   Omicron (BA.3)       Limpopo 2.566196e-04 3.765737e-04 NA -4.814512e-04 9.946905e-04
# 32   Omicron (BA.2)       Limpopo 9.468073e-01 2.919669e-02 NA  8.895828e-01 1.004032e+00
# 33   Omicron (BA.1)  Western Cape 3.667600e-01 5.454928e-02 NA  2.598454e-01 4.736746e-01
# 34             Beta  Western Cape 3.415944e-13 2.804667e-13 NA -2.081101e-13 8.912990e-13
# 35            Alpha  Western Cape 4.011550e-13 3.456762e-13 NA -2.763580e-13 1.078668e-12
# 36            Other  Western Cape 1.731170e-15 1.435126e-15 NA -1.081625e-15 4.543966e-15
# 37            Delta  Western Cape 1.496382e-08 1.198472e-08 NA -8.525792e-09 3.845344e-08
# 38 Omicron (BA.1.1)  Western Cape 7.929519e-02 1.894714e-02 NA  4.215948e-02 1.164309e-01
# 39   Omicron (BA.3)  Western Cape 1.326055e-03 1.182525e-03 NA -9.916519e-04 3.643763e-03
# 40   Omicron (BA.2)  Western Cape 5.526187e-01 6.491374e-02 NA  4.253901e-01 6.798473e-01
# 41   Omicron (BA.1)    Free State 1.437739e-01 5.014527e-02 NA  4.549100e-02 2.420568e-01
# 42             Beta    Free State 2.183627e-13 2.237218e-13 NA -2.201239e-13 6.568493e-13
# 43            Alpha    Free State 7.382957e-14 8.197887e-14 NA -8.684608e-14 2.345052e-13
# 44            Other    Free State 1.256480e-15 1.298418e-15 NA -1.288373e-15 3.801332e-15
# 45            Delta    Free State 3.780011e-09 3.804153e-09 NA -3.675993e-09 1.123601e-08
# 46 Omicron (BA.1.1)    Free State 4.497409e-02 2.036651e-02 NA  5.056468e-03 8.489171e-02
# 47   Omicron (BA.3)    Free State 1.441592e-05 1.053945e-04 NA -1.921536e-04 2.209854e-04
# 48   Omicron (BA.2)    Free State 8.112376e-01 6.532297e-02 NA  6.832069e-01 9.392682e-01
# 49   Omicron (BA.1) KwaZulu Natal 1.119508e-01 2.552424e-02 NA  6.192422e-02 1.619774e-01
# 50             Beta KwaZulu Natal 9.724966e-14 8.045239e-14 NA -6.043413e-14 2.549335e-13
# 51            Alpha KwaZulu Natal 2.909301e-14 2.701478e-14 NA -2.385498e-14 8.204100e-14
# 52            Other KwaZulu Natal 9.831563e-16 8.204157e-16 NA -6.248289e-16 2.591142e-15
# 53            Delta KwaZulu Natal 5.199389e-09 4.173114e-09 NA -2.979763e-09 1.337854e-08
# 54 Omicron (BA.1.1) KwaZulu Natal 1.260842e-02 4.366451e-03 NA  4.050333e-03 2.116651e-02
# 55   Omicron (BA.3) KwaZulu Natal 2.443458e-04 2.598398e-04 NA -2.649309e-04 7.536224e-04
# 56   Omicron (BA.2) KwaZulu Natal 8.751964e-01 2.828031e-02 NA  8.197680e-01 9.306248e-01
# 57   Omicron (BA.1)  Eastern Cape 9.077089e-01 6.363489e-02 NA  7.829868e-01 1.032431e+00
# 58             Beta  Eastern Cape 6.726394e-13 5.968914e-13 NA -4.972463e-13 1.842525e-12
# 59            Alpha  Eastern Cape 2.006583e-13 2.016206e-13 NA -1.945108e-13 5.958275e-13
# 60            Other  Eastern Cape 3.788481e-15 3.392806e-15 NA -2.861296e-15 1.043826e-14
# 61            Delta  Eastern Cape 4.534457e-08 3.920810e-08 NA -3.150189e-08 1.221910e-07
# 62 Omicron (BA.1.1)  Eastern Cape 9.216500e-02 6.363432e-02 NA -3.255597e-02 2.168860e-01
# 63   Omicron (BA.3)  Eastern Cape 1.260856e-04 9.933631e-04 NA -1.820870e-03 2.073042e-03
# 64   Omicron (BA.2)  Eastern Cape 4.365522e-11 1.198528e-11 NA  2.016450e-11 6.714594e-11
# 65   Omicron (BA.1) Northern Cape 1.632227e-01 7.792055e-02 NA  1.050122e-02 3.159442e-01
# 66             Beta Northern Cape 6.338893e-13 5.988363e-13 NA -5.398083e-13 1.807587e-12
# 67            Alpha Northern Cape 1.213359e-13 1.323889e-13 NA -1.381415e-13 3.808133e-13
# 68            Other Northern Cape 3.737413e-15 3.564156e-15 NA -3.248204e-15 1.072303e-14
# 69            Delta Northern Cape 1.015501e-08 9.362621e-09 NA -8.195388e-09 2.850541e-08
# 70 Omicron (BA.1.1) Northern Cape 3.988176e-02 2.292956e-02 NA -5.059362e-03 8.482287e-02
# 71   Omicron (BA.3) Northern Cape 6.074485e-04 8.377887e-04 NA -1.034587e-03 2.249484e-03
# 72   Omicron (BA.2) Northern Cape 7.962881e-01 9.696549e-02 NA  6.062392e-01 9.863370e-01





# PLOT MULTINOMIAL FIT

extrapolate = 30
date.from = as.numeric(min(GISAID_sel$date_num))
date.to = max(GISAID_sel$date_num)+extrapolate

# multinomial model predictions by province (fastest, but no confidence intervals)
predgrid = expand.grid(list(date_num=seq(date.from, date.to),
                            province=levels_PROVINCES))
fit_southafrica_multi_preds_byprovince = data.frame(predgrid, as.data.frame(predict(fit1_southafrica_province_multi, newdata=predgrid, type="prob")),check.names=F)
fit_southafrica_multi_preds_byprovince = gather(fit_southafrica_multi_preds_byprovince, variant, prob, all_of(levels_VARIANTS), factor_key=TRUE)
fit_southafrica_multi_preds_byprovince$date = as.Date(fit_southafrica_multi_preds_byprovince$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_byprovince$variant = factor(fit_southafrica_multi_preds_byprovince$variant, levels=levels_VARIANTS)
fit_southafrica_multi_preds_byprovince$province = factor(fit_southafrica_multi_preds_byprovince$province, levels=levels_PROVINCES)

# multinomial model predictions overall for South Africa (with model without province) (fastest, but no confidence intervals)
predgrid = expand.grid(list(date_num=seq(date.from, date.to)))
fit_southafrica_multi_preds = data.frame(predgrid, as.data.frame(predict(fit2_southafrica_multi, newdata=predgrid, type="prob")),check.names=F)
fit_southafrica_multi_preds = gather(fit_southafrica_multi_preds, variant, prob, all_of(levels_VARIANTS), factor_key=TRUE)
fit_southafrica_multi_preds$date = as.Date(fit_southafrica_multi_preds$date_num, origin="1970-01-01")
fit_southafrica_multi_preds$variant = factor(fit_southafrica_multi_preds$variant, levels=levels_VARIANTS)


# Muller plot for South Africa by province
muller_southafrica_province_mfit = ggplot(data=fit_southafrica_multi_preds_byprovince, 
                                 aes(x=date, y=prob, group=variant)) + 
  facet_wrap(~ province) +
  scale_y_continuous(expand=c(0,0)) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  annotate("rect", xmin=max(GISAID_sel$date_num)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA","(GISAID data, multinomial fit)")
muller_southafrica_province_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_byprovince_multinom fit.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_byprovince_multinom fit.pdf"), width=10, height=6)

library(ggpubr)
ggarrange(muller_southafrica_byprovince_raw + coord_cartesian(xlim=c(min(GISAID_sel$date),max(GISAID_sel$date)+extrapolate)) +
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_southafrica_province_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots multipanel_byprovince_multinom fit.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots multipanel_byprovince_multinom fit.pdf"), width=10, height=10)


# multinomial model predictions with confidence intervals for all of South Africa combined (slower)
# using fit without province fit2_southafrica_multi here
fit_southafrica_multi_preds_withCI = data.frame(emmeans(fit2_southafrica_multi,
                                                        ~ variant,
                                                        by=c("date_num"),
                                                        at=list(date_num=seq(date.from, date.to, by=7)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_southafrica_multi_preds_withCI$date = as.Date(fit_southafrica_multi_preds_withCI$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_withCI$variant = factor(fit_southafrica_multi_preds_withCI$variant, levels=levels_VARIANTS)

# multinomial model predictions with confidence intervals by province (slower)
fit_southafrica_multi_preds_byprovince_withCI = data.frame(emmeans(fit1_southafrica_province_multi,
                                                        ~ variant,
                                                        by=c("date_num","province"),
                                                        at=list(date_num=seq(date.from, date.to, by=7)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_southafrica_multi_preds_byprovince_withCI$date = as.Date(fit_southafrica_multi_preds_byprovince_withCI$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_byprovince_withCI$variant = factor(fit_southafrica_multi_preds_byprovince_withCI$variant, levels=levels_VARIANTS)

# Muller plot overall for South Africa
muller_southafrica_mfit = ggplot(data=fit_southafrica_multi_preds_withCI, 
                                   aes(x=date, y=prob, group=variant)) + 
  scale_y_continuous(expand=c(0,0)) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  annotate("rect", xmin=max(GISAID_sel$date_num)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.2, fill="white") + # extrapolated part
  xaxis +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID data, multinomial fit)")
muller_southafrica_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_multinom fit.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots_multinom fit.pdf"), width=10, height=6)

library(ggpubr)
ggarrange(muller_southafrica_raw + coord_cartesian(xlim=c(min(GISAID_sel$date),max(GISAID_sel$date)+extrapolate))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_southafrica_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots multipanel_multinom fit.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_muller plots multipanel_multinom fit.pdf"), width=10, height=10)



# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# BY PROVINCE
# on logit scale:

fit_southafrica_multi_preds_byprovince2 = fit_southafrica_multi_preds_byprovince_withCI
ymin = 0.001
ymax = 0.999
fit_southafrica_multi_preds_byprovince2$asymp.LCL[fit_southafrica_multi_preds_byprovince2$asymp.LCL<ymin] = ymin
fit_southafrica_multi_preds_byprovince2$asymp.UCL[fit_southafrica_multi_preds_byprovince2$asymp.UCL<ymin] = ymin
fit_southafrica_multi_preds_byprovince2$asymp.UCL[fit_southafrica_multi_preds_byprovince2$asymp.UCL>ymax] = ymax
fit_southafrica_multi_preds_byprovince2$prob[fit_southafrica_multi_preds_byprovince2$prob<ymin] = ymin

plot_southafrica_byprovince_mfit_logit = qplot(data=fit_southafrica_multi_preds_byprovince2, x=date, y=prob, geom="blank") +
  facet_wrap(~ province) +
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
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  geom_point(data=data_agbyweek,
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.1, 2), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_southafrica_byprovince_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_byprovince_logit scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_byprovince_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_southafrica_byprovince_mfit = qplot(data=fit_southafrica_multi_preds_byprovince2, x=date, y=100*prob, geom="blank") +
  facet_wrap(~ province) +
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
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  geom_point(data=data_agbyweek,
             aes(x=date, y=100*prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\ngenotyped", trans="sqrt",
                        range=c(0.1, 2), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_southafrica_byprovince_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_byprovince_response scale_with points.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_byprovince_response scale.pdf"), width=10, height=6)



# OVERALL FOR SOUTH AFRICA
# using fit without province fit2_southafrica_multi here
# on logit scale:

fit_southafrica_multi_preds2 = fit_southafrica_multi_preds_withCI
ymin = 0.001
ymax = 0.9995
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
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.995), expand=c(0,0))
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
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
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
cases_cum$tests_daily = abs(cases_cum$tests_daily)
head(cases_cum)
cases_cum$date = as.Date(cases_cum$date, "%d-%m-%Y")
qplot(data=cases_cum, x=date, y=cases_daily, geom="col")
names(cases_cum)

qplot(data=cases_cum, x=date, y=tests_daily, geom="col")

# for more data see
# https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=provincial-tests
# https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=deaths-age
# https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=deaths-gender
# https://raw.githubusercontent.com/alex1770/Covid-19/master/VOCgrowth/EarlyOmicronEstimate/extracthospdata/SouthAfricaHospData.json
# https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/EarlyOmicronEstimate/SAcasecounts.csv
# https://github.com/ggilestro/playground/tree/main/south_africa_data

# case data for SA per province
cases_prov = read.csv("https://mediahack.co.za/datastories/coronavirus/data/csv.php?table=provinces")
cases_prov = cases_prov[cases_prov$date!="",]
cases_prov = cases_prov[cases_prov$province!="Unknown",]
cases_prov$province = factor(cases_prov$province, levels=c("Gauteng", "North West", "Mpumalanga", "Limpopo", "Western Cape", "Free State", "KwaZulu Natal", "Eastern Cape", "Northern Cape"))
head(cases_prov)
cases_prov$date = as.Date(cases_prov$date)
names(cases_cum)
# qplot(data=cases_cum, x=date, y=tests_daily, geom="col")

# hospitalisation data for SA per province
hosp_prov = read.csv("https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_provincial_raw_hospitalization.csv")
colnames(hosp_prov)[1] = "province"
colnames(hosp_prov)[3] = "date"
hosp_prov = hosp_prov[hosp_prov$province!="Total",]
hosp_prov = hosp_prov[hosp_prov$Owner=="Total",]
hosp_prov$province = factor(hosp_prov$province, 
                            levels=c("Gauteng", "NorthWest", "Mpumalanga", "Limpopo", "WesternCape", "FreeState", "KwaZulu-Natal", "EasternCape", "NorthernCape"),
                            labels=c("Gauteng", "North West", "Mpumalanga", "Limpopo", "Western Cape", "Free State", "KwaZulu Natal", "Eastern Cape", "Northern Cape"))
head(hosp_prov)
hosp_prov$date = as.Date(hosp_prov$date)
unique(hosp_prov$variable) 
# "FacilitiesReporting"     "AdmissionstoDate"        "DiedtoDate"              "Dischargedtodate"        "CurrentlyAdmitted"       "CurrentlyinICU"         
# "CurrentlyVentilated"     "CurrentlyOxygenated"     "AdmissionsinPreviousDay"
colnames(hosp_prov)
# "province" "Owner"    "date"     "variable" "value"   

# combined dataset with cases, hospitalisation & deaths per province
data_prov = hosp_prov[hosp_prov$Owner=="Total",]
data_prov$Owner = NULL
data_prov$variable = as.character(data_prov$variable)
data_prov = rbind(data_prov,
                  data.frame(province=cases_prov$province, date=cases_prov$date, variable="daily_cases", value=cases_prov$daily_cases),
                  data.frame(province=cases_prov$province, date=cases_prov$date, variable="daily_deaths", value=cases_prov$daily_deaths))
data_prov$variable = as.factor(data_prov$variable)
str(data_prov)

# convert to wide format
library(tidyr)
data_prov_wide = spread(data_prov, variable, value)
# fill in missing dates with NAs
data_prov_wide2 = expand.grid(date=seq(min(data_prov_wide$date), max(data_prov_wide$date), by=1), province=levels(data_prov_wide$province)) 
hosp_prov_wide = merge(data_prov_wide2, data_prov_wide, all.x=TRUE)

# hosp_prov_wide = as.data.frame(hosp_prov_wide %>% # test to see if diff(AdmissionstoDate) has less reporting lag than AdmissionsinPreviousDay
#                       group_by(province) %>%
#                       mutate(NewAdmissions = AdmissionstoDate - lag(AdmissionstoDate)))

# plots of cases, hospitalisations & deaths by province
qplot(data=data_prov_wide, x=date, y=daily_cases, geom="col") + facet_wrap(~ province, scale="free_y")  + 
  ylab("new cases per day") + xlab("reporting date") + ggtitle("New confirmed Covid cases per day\nin South Africa") + xaxis # + scale_y_log10()
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases by province.png"), width=10, height=6)

qplot(data=data_prov_wide, x=date, y=daily_cases, geom="col", fill=I("blue"), width=I(1.2)) + facet_wrap(~ province, scale="free_y")  + 
  geom_col(aes(x=date-6, y=-AdmissionsinPreviousDay*7), fill=I("red"), width=I(1.2)) +
  xlab("case reporting date") +
  ylab("cases per day (top) & hospitalisations per week 7 days later (bottom)") + 
  ggtitle("Confirmed Covid cases & hospitalisations in South Africa",
          "New confirmed Covid cases per day (blue)\n& hospitalised patients testing positive per week one week later (red), data NICD") + xaxis + # + scale_y_log10()
  coord_cartesian(xlim=c(as.Date("2020-11-01"), NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day and hosps per week mirrored by province.png"), width=10, height=6)

qplot(data=hosp_prov_wide, x=date, y=AdmissionsinPreviousDay, geom="col", width=I(1.3)) + facet_wrap(~ province, scale="free_y") + # + # + scale_y_log10()
  ylab("new admissions per day") + ggtitle("Covid patients admitted in hospital per day\nin South Africa") + xaxis +
  coord_cartesian(xlim=c(as.Date("2020-11-01"), NA))
  # facet_wrap(~ variable)
ggsave(file=paste0(".\\plots\\",plotdir,"\\new hospital admissions by province.png"), width=10, height=6)

# qplot(data=hosp_prov_wide, x=date, y=NewAdmissions, geom="col") + facet_wrap(~ province+Owner, scale="free_y") + # + # + scale_y_log10()
#   ylab("new admissions per day") + ggtitle("Covid patients admitted in hospital per day\nin South Africa") + xaxis
# # facet_wrap(~ variable)

# qplot(data=hosp_prov_wide, x=date, y=AdmissionstoDate, geom="col") + facet_wrap(~ province+Owner, scale="free_y") + # + # + scale_y_log10()
#   ylab("new admissions per day") + ggtitle("Total nr of Covid patients admitted to date\nin South Africa") + xaxis
# facet_wrap(~ variable)



qplot(data=hosp_prov[hosp_prov$variable=="CurrentlyAdmitted",], x=date, y=value, geom="col", width=I(1.3)) + facet_wrap(~ province, scale="free_y") +
  ylab("currently admitted") + ggtitle("Covid patients currently admitted in hospital\nin South Africa") + xaxis # + # + scale_y_log10()
# facet_wrap(~ variable)
ggsave(file=paste0(".\\plots\\",plotdir,"\\total hospital admissions by province.png"), width=10, height=6)

qplot(data=cases_prov, x=date, y=daily_deaths, geom="col", width=I(1.3)) + facet_wrap(~ province, scale="free_y") + 
  ylab("new deaths per day") + ggtitle("Covid mortality in South Africa") + xaxis # + scale_y_log10()
ggsave(file=paste0(".\\plots\\",plotdir,"\\deaths by province.png"), width=10, height=6)


# cases_tot = as.data.frame(get_national_data(countries = "South Africa"))
# # cases_tot = cases_tot[cases_tot$date>=as.Date("2020-01-01"),]
# cases_tot$date_num = as.numeric(cases_tot$date)
# # cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
# cases_tot$WEEKDAY = weekdays(cases_tot$date)
# cases_tot[cases_tot$date==as.Date("2021-11-24"),"cases_new"] = 868 # https://twitter.com/nicd_sa/status/1463200722615513093
# # cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
# range(cases_tot$date) # "2020-01-03" "2021-12-03"

cases_tot = cases_cum # TO DO: replace by sum of per province data / data_prov_wide (appears more up to date)
cases_tot$cases_new = cases_tot$cases_daily
cases_tot$tests_new = abs(cases_tot$tests_daily)
cases_tot$date_num = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot$testspercase = (cases_tot$tests_new+1)/(cases_tot$cases_new+1)
qplot(data=cases_tot, x=date, y=cases_new, geom="col")
qplot(data=cases_tot, x=date, y=tests_new, geom="col")
qplot(data=cases_tot, x=date, y=testspercase, geom="col")
qplot(data=cases_tot, x=date, y=1/testspercase, geom="col") + ylab("positivity rate") + ylim(c(0,1))


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
k=45
fit_cases = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) +
                  WEEKDAY + # + offset(log(testspercase)),
                  # BANKHOLIDAY,
                + s(tests_new, bs="cs", k=8, fx=F),
                family=nb(), data=cases_tot,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+1*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-1*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 
BIC(fit_cases) # 9953.436

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
  geom_line() + theme_hc() + xlab("Date of infection") + ylab("Rt") +
  # xaxis +
  # scale_y_continuous(limits=c(1/2, 4), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Rt VALUES IN SOUTH AFRICA","Based on negative binomial GAM fit to case data NICD,\ncorrecting for variable testing intensity &\nassuming gamma distr generation time of 4.7d with SD of 2.9d\nand time from infection to diagnosis of 7d") +
  labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(0.5,NA))
  # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

head(avg_r_cases,40)
avg_r_cases[avg_r_cases$DATE_OF_INFECTION==today,"Re"]/avg_r_cases[avg_r_cases$DATE_OF_INFECTION==as.Date("2021-09-23"),"Re"] # Re now = 3.0 / Re on 23 Sept = 0.73 = x4.1
avg_r_cases[avg_r_cases$DATE_OF_INFECTION==today,"Re"]/avg_r_cases[avg_r_cases$DATE_OF_INFECTION==as.Date("2021-10-07"),"Re"] # Re now (3.02) / Re on 7 Oct (0.81) = x3.7

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values South Africa.png"), width=8, height=6)
# PS there were reporting problems/delays lately so not sure about stagnation of Re value


qplot(data=avg_r_cases, x=DATE_OF_INFECTION, y=r*700, ymin=r_LOWER*700, ymax=r_UPPER*700, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) + 
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") + ylab("Continuous growth in infections per week (%)") + 
  xaxis +
  # xaxis +
  # scale_y_continuous(limits=c(1/2, 4), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("WEEKLY GROWTH/DECLINE IN SARS-CoV2\nINFECTIONS IN SOUTH AFRICA","Based on negative binomial GAM fit to case data NICD,\ncorrecting for variable testing intensity &\nassuming time from infection to diagnosis of 7d") +
  labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(-100,NA))
# +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

head(avg_r_cases,40)
avg_r_cases[avg_r_cases$DATE_OF_INFECTION==today,"Re"]/avg_r_cases[avg_r_cases$DATE_OF_INFECTION==as.Date("2021-09-23"),"Re"] # Re now = 3.0 / Re on 23 Sept = 0.73 = x4.1
avg_r_cases[avg_r_cases$DATE_OF_INFECTION==today,"Re"]/avg_r_cases[avg_r_cases$DATE_OF_INFECTION==as.Date("2021-10-07"),"Re"] # Re now (3.02) / Re on 7 Oct (0.81) = x3.7

ggsave(file=paste0(".\\plots\\",plotdir,"\\growth infections cases South Africa.png"), width=8, height=6)
# PS there were reporting problems/delays lately so not sure about stagnation of Re value


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
  geom_line() + theme_hc() + xlab("Date of infection") + ylab("Rt") +
  # xaxis +
  scale_y_continuous(limits=c(1/4, 4), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Rt VALUES IN GAUTENG PROVINCE, SA","Based on negative binomial GAM fit to case data NICD,\nassuming gamma distr generation time of 4.7d with SD of 2.9d\nand time from infection to diagnosis of 7d") +
  labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

head(avg_r_cases_gauteng,40)
avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==today,"Re"]/avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==as.Date("2021-09-23"),"Re"] # Re now / Re on 23 Sept = x2.2
avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==today,"Re"]/avg_r_cases_gauteng[avg_r_cases_gauteng$DATE_OF_INFECTION==as.Date("2021-10-07"),"Re"] # Re now / Re on 14 Oct = x2.2

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values Gauteng Province.png"), width=8, height=6)
tail(avg_r_cases_gauteng)

qplot(data=avg_r_cases_gauteng, x=DATE_OF_INFECTION, y=r*700, ymin=r_LOWER*700, ymax=r_UPPER*700, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) + 
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") + ylab("Continuous growth in infections per week (%)") + 
  xaxis +
  # xaxis +
  # scale_y_continuous(limits=c(1/2, 4), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("WEEKLY GROWTH/DECLINE IN SARS-CoV2\nINFECTIONS IN GAUTENG PROVINCE, SA","Based on negative binomial GAM fit to case data NICD,\nassuming time from infection to diagnosis of 7d") +
  labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(-100,NA))
# +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\growth infections cases Gauteng Province.png"), width=8, height=6)






# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_southafrica_multi_preds$totcases = cases_tot$cases_new[match(round(fit_southafrica_multi_preds$date_num),cases_tot$date_num)]
fit_southafrica_multi_preds$cases = fit_southafrica_multi_preds$totcases * fit_southafrica_multi_preds$prob
fit_southafrica_multi_preds$cases[fit_southafrica_multi_preds$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ date_num, at=list(date_num=seq(date.from, date.to, by=1),
                                                                     testspercase=20,
                                                                     tests_new=max(cases_tot$tests_new, na.rm=T)), 
                                      type="response"))
fit_southafrica_multi_preds$smoothed_totcases = cases_emmeans$response[match(fit_southafrica_multi_preds$date_num,cases_emmeans$date_num)]
fit_southafrica_multi_preds$smoothed_cases = fit_southafrica_multi_preds$smoothed_totcases * fit_southafrica_multi_preds$prob
fit_southafrica_multi_preds$smoothed_cases[fit_southafrica_multi_preds$smoothed_cases<=0.001] = NA
fit_southafrica_multi_preds$variant = factor(fit_southafrica_multi_preds$variant, levels=levels_VARIANTS)

fit_southafrica_multi_preds[fit_southafrica_multi_preds$date==as.Date("2021-12-03"),] # this date is not plotted in plot below FIX

qplot(data=fit_southafrica_multi_preds[fit_southafrica_multi_preds$variant=="Omicron",], x=date, y=cases, geom="line") +
  xlim(c(as.Date("2021-10-01"), today)) + scale_y_log10() + 
  geom_line(aes(y=smoothed_cases), colour=I("red"), lwd=I(2)) +
  ylab("new Omicron cases per day") +
  ggtitle("Inferred nr. of new Omicron cases per day in SA", "(black=observed,\nred=GAM negative binomial spline smooth with correction for\nvariable testing intensity & weekday effects & marginal means\ncalculated at uniform, maximal testing effort)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\omicron cases observed and fitted.png"), width=8, height=6)

fit_southafrica_multi_preds[fit_southafrica_multi_preds$variant=="Omicron","smoothed_cases"][fit_southafrica_multi_preds$date==today]
# 36106.74
# matches latest figure of 37

ggplot(data=fit_southafrica_multi_preds, 
       aes(x=date, y=cases, group=variant)) + 
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant, width=I(1.1)), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA","(case data NICD plus multinomial spline fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)
write.csv(fit_southafrica_multi_preds_withCI, file=paste0(".\\plots\\",plotdir,"\\cases per day by variant South Africa 5 Jan 2022.csv"), row.names=F)

ggplot(data=fit_southafrica_multi_preds[fit_southafrica_multi_preds$date<=today,],  # or (today-1) or (today-2) to be on the safe side
       aes(x=date, y=smoothed_cases, group=variant)) + 
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA","(negative binomial fit to NICD case data, with correction for weekday effects &\nvariable testing intensity; marginal means calculated at constant, maximal testing effort;\nvariant frequencies based on multinomial spline fit to GISAID data)") +
  scale_fill_manual("variant", values=lineage_cols) +
  scale_colour_manual("variant", values=lineage_cols) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit case data.png"), width=8, height=6)







# DIDN'T TRY TO RUN / UPDATE THE PART BELOW YET

# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for simplest model fit1_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from+16,
                                                  date.to-extrapolate), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit1_southafrica_multi, ~ variant , at=list(date_num=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_VARIANTS)
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
                                       levels=1:length(levels_VARIANTS), 
                                       labels=levels_VARIANTS)
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
above_avg_r_variants$variant = factor(above_avg_r_variants$variant, levels=c(levels_VARIANTS,"avg"))
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
  scale_fill_manual("variant", values=c(head(lineage_cols,-1),"black")) +
  scale_colour_manual("variant", values=c(head(lineage_cols,-1),"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)
      
