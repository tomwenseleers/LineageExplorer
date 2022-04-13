# ANALYSIS OF GROWTH ADVANTAGE OF OMICRON (B.1.1.529/BA.XX) IN SOUTH AFRICA BASED ON GISAID SEQUENCE DATA

# T. Wenseleers
# last update 12 APRIL 2022
    
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
  
  sel_target_VOC = "Omicron (BA.5)"
  sel_reference_VOC = "Omicron (BA.1*)"
  levels_VARIANTS = c(sel_reference_VOC, "Beta", "Alpha", "Other", "Delta", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", sel_target_VOC)
  levels_VARIANTS_plot = c("Other", "Beta", "Alpha", "Delta", "Omicron (BA.1*)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", "Omicron (BA.5)")
  
  n = length(levels_VARIANTS)
  lineage_cols = hcl(h = seq(0, 180, length = n), l = 60, c = 180)
  lineage_cols[which(levels_VARIANTS=="Alpha")] = "#0085FF"
  lineage_cols[which(levels_VARIANTS=="Beta")] = "green4"
  lineage_cols[which(levels_VARIANTS=="Delta")] = "mediumorchid"
  # lineage_cols[which(levels_VARIANTS=="C.1.2")] = "darkorange"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.1*)")] = "red" # "magenta"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.2)")] = "red3"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.3)")] = "red4" 
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.4)")] = "darkorange" 
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.5)")] = "darkorange3" 
  lineage_cols[which(levels_VARIANTS=="Other")] = "grey65"
  
  lineage_cols_plot = lineage_cols[match(levels_VARIANTS_plot,levels_VARIANTS)]
  
  # X axis for plots
  firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
  xaxis = scale_x_continuous(breaks=firststofmonth,
                             labels=substring(months(firststofmonth),1,1),
                             expand=c(0,0))
  
  
  # import GISAID records for South Africa
d1 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_15_subm_2020.tsv"), col_types = cols(.default = "c")) 
d2 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_16_subm_jan_may_2021.tsv"), col_types = cols(.default = "c")) 
d3 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_12_19_subm_june_aug_2021.tsv"), col_types = cols(.default = "c")) 
d4 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_16_subm_sept_2021_dec_2021.tsv"), col_types = cols(.default = "c")) 
d5 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_12_19_subm_jan_feb_2022.tsv"), col_types = cols(.default = "c")) 
d6 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_16_subm_mar_apr_2022.tsv"), col_types = cols(.default = "c")) 

# parse GISAID & SGTF data
GISAID = as.data.frame(rbind(d1,d2,d3,d4,d5,d6))
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
range(GISAID$date) # "2020-03-06" "2022-04-03"
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
                                                   "Gauteng", "Western Cape", "KwaZulu Natal", "KwaZulu Natal", "KwaZulu Natal",
                                                   "Limpopo", "Mpumalanga", "North West",
                                                   "Northern Cape", "Northern Cape", # NA,
                                                   "Western Cape", "Western Cape"))
GISAID$province = factor(GISAID$province, levels=levels_PROVINCES)

GISAID$variant = case_when(
  # grepl("N679K", aa_substitutions) & grepl("H655Y", aa_substitutions) & grepl("P681H", aa_substitutions) ~ "B.1.1.529",
  (grepl("BA.1", GISAID$pango_lineage) ~ "Omicron (BA.1*)"),
  (GISAID$pango_lineage=="BA.3") ~ "Omicron (BA.3)",
  ((grepl("BA.2",GISAID$pango_lineage))&
    ((grepl("L452R", GISAID$aa_substitutions)&
            grepl("486V", GISAID$aa_substitutions)&
            grepl("11F", GISAID$aa_substitutions)&
            (!grepl("D3N",GISAID$aa_substitutions)) ))) ~ "Omicron (BA.4)",
  ((grepl("BA.2",GISAID$pango_lineage))&
    ((grepl("L452R",GISAID$aa_substitutions)&
       grepl("486V",GISAID$aa_substitutions)&
     (!grepl("11F", GISAID$aa_substitutions))&
       grepl("D3N",GISAID$aa_substitutions)))) ~ "Omicron (BA.5)",
  (grepl("BA.2",GISAID$pango_lineage)) ~ "Omicron (BA.2)",
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
nrow(GISAID_sel) # 33180
sum(GISAID_sel$variant==sel_target_VOC) # 47
table(GISAID_sel$variant)
range(GISAID_sel$date)

GISAID_sel$variant = factor(GISAID_sel$variant)
GISAID_sel$variant = relevel(GISAID_sel$variant, ref=sel_reference_VOC) # we code Omicron (BA.2) as the reference

GISAID_sel$variant = factor(GISAID_sel$variant, levels=levels_VARIANTS)
table(GISAID_sel$variant, GISAID_sel$sampling_strategy)

# age distribution of patients infected with Delta or Omicron (GISAID data)
library(ggplot2)
library(ggthemes)
GISAID_sel2 = GISAID_sel[(GISAID_sel$variant %in% c("Delta", "Omicron (BA.1*)")) & (GISAID_sel$floor_date>as.Date("2021-05-01")),]
GISAID_sel2$variant = factor(GISAID_sel2$variant, levels=c("Omicron (BA.1*)", "Delta"))
GISAID_sel2$Year_Month = factor(GISAID_sel2$Year_Month)
GISAID_sel2$Year_Month = droplevels(GISAID_sel2$Year_Month)
GISAID_sel3 = GISAID_sel
GISAID_sel3 = GISAID_sel3[(GISAID_sel3$variant %in% c("Other","Beta","Delta","Omicron (BA.1*)")),]
GISAID_sel3$variant = factor(GISAID_sel3$variant, levels=c("Other","Beta","Delta","Omicron (BA.1*)"))
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
  scale_fill_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1*)","Delta"))]) + theme(legend.position="none") + 
  ggtitle("Age distribution of individuals with sequenced\nOmicron & Delta infections in South Africa","(GISAID data)") +
  theme(strip.text = element_text(size=8),
        strip.background = element_rect(fill="white", colour="white",size=1)) +
  theme(strip.text.x = element_text(margin = margin(0, 0, 0, 0, "cm")))
ggsave(file=paste0(".\\plots\\",plotdir,"\\age distribution GISAID delta omicron by week.png"), width=5, height=7)

ggplot(data=GISAID_sel2, aes(x=date, y=age, fill=variant, colour=variant, group=variant)) + 
  # facet_wrap(~ Year_Month+variant, scale="free_y", ncol=2,drop=FALSE) + 
  geom_point(aes(alpha=I(0.1)), pch=I(16)) +
  geom_smooth(method="gam") +
  scale_fill_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1*)","Delta"))]) + 
  scale_colour_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1*)","Delta"))]) + 
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
  scale_fill_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1*)","Delta"))]) + 
  scale_colour_manual(values=lineage_cols[which(levels_VARIANTS %in% c("Omicron (BA.1*)","Delta"))]) + 
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
data_agbyweek_province_sum = aggregate(count ~ floor_date+province, data=data_agbyweek_province, sum)
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
data_agbyweek$variant = factor(data_agbyweek$variant, levels=levels_VARIANTS_plot)
muller_southafrica_raw = ggplot(data=data_agbyweek, aes(x=date, y=count, group=variant)) + 
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
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
data_agbyweek_province$variant = factor(data_agbyweek_province$variant, levels=levels_VARIANTS_plot)
muller_southafrica_byprovince_raw = ggplot(data=data_agbyweek_province, aes(x=date, y=count, group=variant)) + 
  facet_wrap(~ province) +
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill", width=I(7)) + # or geom_area
  scale_fill_manual("", values=lineage_cols_plot) +
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

# avg growth rate advantage of Omicron (BA.2) over Omicron BA.1* (difference in growth rate per day) (today)
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
# contrast    estimate          SE df    asymp.LCL     asymp.UCL       p.value
# 1           Beta - (Omicron (BA.1*))  0.05547020 0.002342672 NA  0.050878649  0.0600617559 6.052717e-124
# 2          Alpha - (Omicron (BA.1*))  0.01922592 0.009566829 NA  0.000475276  0.0379765546  4.446893e-02
# 3          Other - (Omicron (BA.1*))  0.12819447 0.001793306 NA  0.124679651  0.1317092816  0.000000e+00
# 4          Delta - (Omicron (BA.1*))  0.06822738 0.002940012 NA  0.062465061  0.0739896950 3.914732e-119
# 5 Omicron (BA.2) - (Omicron (BA.1*))  0.06096365 0.001655560 NA  0.057718813  0.0642084888 7.743800e-297
# 6 Omicron (BA.3) - (Omicron (BA.1*)) -0.01097291 0.005523485 NA -0.021798737 -0.0001470742  4.696773e-02
# 7 Omicron (BA.4) - (Omicron (BA.1*))  0.16036423 0.009722571 NA  0.141308343  0.1794201200  4.051082e-61
# 8 Omicron (BA.5) - (Omicron (BA.1*))  0.24383638 0.020866420 NA  0.202938948  0.2847338128  1.510316e-31

# corresponding transmission advantage of Omicron BA.5 over BA.1* (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_avg[8,5]*2.2) # 1.56x
exp(delta_r_southafrica_avg[8,2]*2.2) # Omicron 1.71x transmission advantage over BA.1*
exp(delta_r_southafrica_avg[8,6]*2.2) # 1.87x


# avg growth rate advantage of Omicron BA.2 over BA.1* (difference in growth rate per day) (today)
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
# contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1           Beta - (Omicron (BA.1*)) -0.122995565 0.002069553 NA -0.127051813 -0.11893932 0.000000e+00
# 2          Alpha - (Omicron (BA.1*)) -0.113599746 0.002202509 NA -0.117916585 -0.10928291 0.000000e+00
# 3          Other - (Omicron (BA.1*)) -0.129629272 0.002075849 NA -0.133697861 -0.12556068 0.000000e+00
# 4          Delta - (Omicron (BA.1*)) -0.092355453 0.002009159 NA -0.096293333 -0.08841757 0.000000e+00
# 5 Omicron (BA.2) - (Omicron (BA.1*))  0.079367380 0.001578636 NA  0.076273311  0.08246145 0.000000e+00
# 6 Omicron (BA.3) - (Omicron (BA.1*))  0.003183872 0.004740990 NA -0.006108298  0.01247604 5.018622e-01
# 7 Omicron (BA.4) - (Omicron (BA.1*))  0.192727698 0.012595247 NA  0.168041467  0.21741393 7.457598e-53
# 8 Omicron (BA.5) - (Omicron (BA.1*))  0.206645397 0.020869395 NA  0.165742134  0.24754866 4.086879e-23

# corresponding transmission advantage of Omicron BA.5 over BA.1* (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_province_avg[8,5]*2.2) # 1.44x
exp(delta_r_southafrica_province_avg[8,2]*2.2) # Omicron 1.57x transmission advantage over BA.1*
exp(delta_r_southafrica_province_avg[8,6]*2.2) # 1.73x


emtrsouthafrica_province_avg_pairw = emtrends(fit1_southafrica_province_multi, pairwise ~ variant,  
                                        var="date_num",  mode="latent",
                                        at=list(date_num=seq(today_num, today_num, by=1)),
                                        adjust="none", df=NA)
delta_r_southafrica_province_avg_pairw = data.frame(confint(emtrsouthafrica_province_avg_pairw, 
                                                      adjust="none", df=NA)$contrasts, 
                                              p.value=as.data.frame(emtrsouthafrica_province_avg_pairw$contrasts,
                                                                    adjust="none", df=NA)$p.value)
delta_r_southafrica_province_avg_pairw
#                              contrast     estimate           SE df    asymp.LCL    asymp.UCL       p.value
# 1            (Omicron (BA.1*)) - Beta  0.122995565 0.0020695525 NA  0.118939316  0.127051813  0.000000e+00
# 2           (Omicron (BA.1*)) - Alpha  0.113599746 0.0022025089 NA  0.109282908  0.117916585  0.000000e+00
# 3           (Omicron (BA.1*)) - Other  0.129629272 0.0020758489 NA  0.125560683  0.133697861  0.000000e+00
# 4           (Omicron (BA.1*)) - Delta  0.092355453 0.0020091595 NA  0.088417573  0.096293333  0.000000e+00
# 5  (Omicron (BA.1*)) - Omicron (BA.2) -0.079367380 0.0015786357 NA -0.082461449 -0.076273311  0.000000e+00
# 6  (Omicron (BA.1*)) - Omicron (BA.3) -0.003183872 0.0047409902 NA -0.012476042  0.006108298  5.018622e-01
# 7  (Omicron (BA.1*)) - Omicron (BA.4) -0.192727698 0.0125952472 NA -0.217413929 -0.168041467  7.457598e-53
# 8  (Omicron (BA.1*)) - Omicron (BA.5) -0.206645397 0.0208693948 NA -0.247548659 -0.165742134  4.086879e-23
# 9                        Beta - Alpha -0.009395818 0.0008053043 NA -0.010974186 -0.007817451  1.870243e-31
# 10                       Beta - Other  0.006633708 0.0001752496 NA  0.006290225  0.006977191  0.000000e+00
# 11                       Beta - Delta -0.030640112 0.0004838418 NA -0.031588424 -0.029691799  0.000000e+00
# 12              Beta - Omicron (BA.2) -0.202362945 0.0026151665 NA -0.207488577 -0.197237313  0.000000e+00
# 13              Beta - Omicron (BA.3) -0.126179437 0.0051351771 NA -0.136244199 -0.116114674 2.543176e-133
# 14              Beta - Omicron (BA.4) -0.315723262 0.0127667025 NA -0.340745539 -0.290700985 5.062407e-135
# 15              Beta - Omicron (BA.5) -0.329640961 0.0209733132 NA -0.370747900 -0.288534023  1.153773e-55
# 16                      Alpha - Other  0.016029526 0.0008171611 NA  0.014427920  0.017631132  1.126411e-85
# 17                      Alpha - Delta -0.021244293 0.0008948590 NA -0.022998185 -0.019490402 1.381238e-124
# 18             Alpha - Omicron (BA.2) -0.192967127 0.0027215972 NA -0.198301359 -0.187632894  0.000000e+00
# 19             Alpha - Omicron (BA.3) -0.116783618 0.0051901820 NA -0.126956188 -0.106611049 4.069862e-112
# 20             Alpha - Omicron (BA.4) -0.306327444 0.0127889283 NA -0.331393283 -0.281261605 8.691820e-127
# 21             Alpha - Omicron (BA.5) -0.320245143 0.0209868497 NA -0.361378613 -0.279111673  1.427189e-52
# 22                      Other - Delta -0.037273819 0.0005101697 NA -0.038273734 -0.036273905  0.000000e+00
# 23             Other - Omicron (BA.2) -0.208996653 0.0026201520 NA -0.214132056 -0.203861249  0.000000e+00
# 24             Other - Omicron (BA.3) -0.132813144 0.0051377181 NA -0.142882887 -0.122743402 2.394994e-147
# 25             Other - Omicron (BA.4) -0.322356970 0.0127677246 NA -0.347381251 -0.297332690 1.197331e-140
# 26             Other - Omicron (BA.5) -0.336274669 0.0209739354 NA -0.377382827 -0.295166511  7.519002e-58
# 27             Delta - Omicron (BA.2) -0.171722833 0.0025676373 NA -0.176755310 -0.166690356  0.000000e+00
# 28             Delta - Omicron (BA.3) -0.095539325 0.0051111279 NA -0.105556952 -0.085521698  5.707016e-78
# 29             Delta - Omicron (BA.4) -0.285083151 0.0127570513 NA -0.310086512 -0.260079789 1.288339e-110
# 30             Delta - Omicron (BA.5) -0.299000850 0.0209674399 NA -0.340096277 -0.257905423  3.870878e-46
# 31    Omicron (BA.2) - Omicron (BA.3)  0.076183508 0.0049147406 NA  0.066550794  0.085816223  3.414109e-54
# 32    Omicron (BA.2) - Omicron (BA.4) -0.113360317 0.0124832506 NA -0.137827039 -0.088893596  1.075882e-19
# 33    Omicron (BA.2) - Omicron (BA.5) -0.127278016 0.0208020478 NA -0.168049281 -0.086506752  9.444071e-10
# 34    Omicron (BA.3) - Omicron (BA.4) -0.189543826 0.0134269440 NA -0.215860152 -0.163227499  2.998001e-45
# 35    Omicron (BA.3) - Omicron (BA.5) -0.203461525 0.0213813333 NA -0.245368168 -0.161554882  1.802362e-21
# 36    Omicron (BA.4) - Omicron (BA.5) -0.013917699 0.0237463261 NA -0.060459643  0.032624245  5.578090e-01


# fitted prop of different variantS today (on average across provinces)
# based on model with province included as a factor
multinom_preds_today_avg_province = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|1,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg_province
#           variant         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1*) 1.251600e-03 2.290513e-04 NA  8.026680e-04 1.700532e-03
# 2            Beta 7.158945e-14 3.096418e-14 NA  1.090077e-14 1.322781e-13
# 3           Alpha 3.505428e-14 1.619606e-14 NA  3.310591e-15 6.679797e-14
# 4           Other 2.224430e-15 9.233796e-16 NA  4.146389e-16 4.034220e-15
# 5           Delta 9.217807e-10 3.389141e-10 NA  2.575213e-10 1.586040e-09
# 6  Omicron (BA.2) 5.747336e-01 6.792880e-02 NA  4.415956e-01 7.078716e-01
# 7  Omicron (BA.3) 2.242774e-05 1.285545e-05 NA -2.768475e-06 4.762395e-05
# 8  Omicron (BA.4) 3.115148e-01 4.691516e-02 NA  2.195628e-01 4.034668e-01
# 9  Omicron (BA.5) 1.124775e-01 5.613106e-02 NA  2.462644e-03 2.224924e-01

# fitted prop of different variantS today (by province)
multinom_preds_today_byprovince = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|province,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_byprovince
#            variant      province         prob           SE df     asymp.LCL    asymp.UCL
# 1  Omicron (BA.1*)       Gauteng 1.216431e-04 4.564192e-05 NA  3.218653e-05 2.110996e-04
# 2             Beta       Gauteng 1.732816e-15 8.859440e-16 NA -3.602540e-18 3.469234e-15
# 3            Alpha       Gauteng 3.496845e-15 1.993762e-15 NA -4.108567e-16 7.404547e-15
# 4            Other       Gauteng 5.169743e-17 2.674038e-17 NA -7.127509e-19 1.041076e-16
# 5            Delta       Gauteng 3.174156e-11 1.551479e-11 NA  1.333134e-12 6.214998e-11
# 6   Omicron (BA.2)       Gauteng 1.169506e-01 3.956095e-02 NA  3.941260e-02 1.944887e-01
# 7   Omicron (BA.3)       Gauteng 1.408978e-06 1.022589e-06 NA -5.952593e-07 3.413216e-06
# 8   Omicron (BA.4)       Gauteng 7.855505e-01 8.347460e-02 NA  6.219433e-01 9.491577e-01
# 9   Omicron (BA.5)       Gauteng 9.737580e-02 6.746438e-02 NA -3.485196e-02 2.296036e-01
# 10 Omicron (BA.1*)    North West 2.244968e-03 4.522264e-04 NA  1.358620e-03 3.131315e-03
# 11            Beta    North West 3.735978e-14 1.621558e-14 NA  5.577820e-15 6.914173e-14
# 12           Alpha    North West 7.750417e-14 4.053937e-14 NA -1.951531e-15 1.569599e-13
# 13           Other    North West 4.153981e-15 1.826524e-15 NA  5.740602e-16 7.733903e-15
# 14           Delta    North West 9.909645e-10 3.992702e-10 NA  2.084092e-10 1.773520e-09
# 15  Omicron (BA.2)    North West 9.959648e-01 6.760506e-02 NA  8.634613e-01 1.128468e+00
# 16  Omicron (BA.3)    North West 8.353396e-05 5.263321e-05 NA -1.962523e-05 1.866932e-04
# 17  Omicron (BA.4)    North West 8.617630e-10 3.820170e-10 NA  1.130234e-10 1.610503e-09
# 18  Omicron (BA.5)    North West 1.706739e-03 6.776162e-02 NA -1.311036e-01 1.345171e-01
# 19 Omicron (BA.1*)    Mpumalanga 7.819253e-04 2.698220e-04 NA  2.530839e-04 1.310767e-03
# 20            Beta    Mpumalanga 4.534481e-14 2.292957e-14 NA  4.036802e-16 9.028594e-14
# 21           Alpha    Mpumalanga 7.455107e-15 5.920733e-15 NA -4.149317e-15 1.905953e-14
# 22           Other    Mpumalanga 2.132877e-15 1.095583e-15 NA -1.442646e-17 4.280181e-15
# 23           Delta    Mpumalanga 3.956304e-10 1.900531e-10 NA  2.313320e-11 7.681276e-10
# 24  Omicron (BA.2)    Mpumalanga 7.288218e-01 2.198012e-01 NA  2.980194e-01 1.159624e+00
# 25  Omicron (BA.3)    Mpumalanga 7.459760e-05 4.743832e-05 NA -1.837981e-05 1.675750e-04
# 26  Omicron (BA.4)    Mpumalanga 2.701323e-01 2.200224e-01 NA -1.611037e-01 7.013683e-01
# 27  Omicron (BA.5)    Mpumalanga 1.893983e-04 9.249534e-03 NA -1.793936e-02 1.831815e-02
# 28 Omicron (BA.1*)       Limpopo 9.005709e-05 5.932573e-05 NA -2.621921e-05 2.063334e-04
# 29            Beta       Limpopo 2.009499e-16 1.639633e-16 NA -1.204122e-16 5.223121e-16
# 30           Alpha       Limpopo 1.911195e-16 1.740481e-16 NA -1.500084e-16 5.322475e-16
# 31           Other       Limpopo 2.023369e-17 1.658640e-17 NA -1.227506e-17 5.274244e-17
# 32           Delta       Limpopo 8.083229e-12 6.452332e-12 NA -4.563110e-12 2.072957e-11
# 33  Omicron (BA.2)       Limpopo 1.766176e-01 1.120666e-01 NA -4.302881e-02 3.962640e-01
# 34  Omicron (BA.3)       Limpopo 2.412232e-06 2.160396e-06 NA -1.822066e-06 6.646531e-06
# 35  Omicron (BA.4)       Limpopo 8.039919e-01 1.427340e-01 NA  5.242383e-01 1.083745e+00
# 36  Omicron (BA.5)       Limpopo 1.929806e-02 9.239354e-02 NA -1.617899e-01 2.003861e-01
# 37 Omicron (BA.1*)  Western Cape 2.965707e-03 5.510400e-04 NA  1.885688e-03 4.045725e-03
# 38            Beta  Western Cape 7.379465e-14 2.860903e-14 NA  1.772198e-14 1.298673e-13
# 39           Alpha  Western Cape 1.003713e-13 4.763487e-14 NA  7.008672e-15 1.937339e-13
# 40           Other  Western Cape 1.518012e-15 6.039840e-16 NA  3.342251e-16 2.701799e-15
# 41           Delta  Western Cape 1.967883e-09 7.021829e-10 NA  5.916297e-10 3.344136e-09
# 42  Omicron (BA.2)  Western Cape 8.483988e-01 9.716175e-02 NA  6.579653e-01 1.038832e+00
# 43  Omicron (BA.3)  Western Cape 2.552573e-05 1.613324e-05 NA -6.094836e-06 5.714629e-05
# 44  Omicron (BA.4)  Western Cape 1.479226e-01 9.731406e-02 NA -4.280950e-02 3.386546e-01
# 45  Omicron (BA.5)  Western Cape 6.874056e-04 8.368261e-03 NA -1.571408e-02 1.708890e-02
# 46 Omicron (BA.1*)    Free State 1.909955e-03 5.039721e-04 NA  9.221882e-04 2.897722e-03
# 47            Beta    Free State 1.588340e-13 7.670602e-14 NA  8.492927e-15 3.091750e-13
# 48           Alpha    Free State 5.740563e-14 3.737156e-14 NA -1.584129e-14 1.306525e-13
# 49           Other    Free State 3.900540e-15 1.926154e-15 NA  1.253473e-16 7.675733e-15
# 50           Delta    Free State 1.434629e-09 6.538043e-10 NA  1.531964e-10 2.716062e-09
# 51  Omicron (BA.2)    Free State 9.894087e-01 1.644862e-01 NA  6.670216e-01 1.311796e+00
# 52  Omicron (BA.3)    Free State 2.308959e-24 1.564378e-24 NA -7.571652e-25 5.375083e-24
# 53  Omicron (BA.4)    Free State 9.223718e-09 4.379316e-09 NA  6.404165e-10 1.780702e-08
# 54  Omicron (BA.5)    Free State 8.681341e-03 1.648032e-01 NA -3.143269e-01 3.316896e-01
# 55 Omicron (BA.1*) KwaZulu Natal 2.958799e-04 1.312595e-04 NA  3.861608e-05 5.531437e-04
# 56            Beta KwaZulu Natal 1.073640e-14 5.993276e-15 NA -1.010201e-15 2.248301e-14
# 57           Alpha KwaZulu Natal 3.331349e-15 2.348940e-15 NA -1.272489e-15 7.935187e-15
# 58           Other KwaZulu Natal 4.463852e-16 2.518836e-16 NA -4.729754e-17 9.400679e-16
# 59           Delta KwaZulu Natal 2.490942e-10 1.327137e-10 NA -1.101991e-11 5.092084e-10
# 60  Omicron (BA.2) KwaZulu Natal 1.546545e-01 6.246007e-02 NA  3.223503e-02 2.770740e-01
# 61  Omicron (BA.3) KwaZulu Natal 1.833963e-06 1.576585e-06 NA -1.256087e-06 4.924013e-06
# 62  Omicron (BA.4) KwaZulu Natal 1.703009e-02 1.821619e-02 NA -1.867298e-02 5.273316e-02
# 63  Omicron (BA.5) KwaZulu Natal 8.280177e-01 6.959247e-02 NA  6.916189e-01 9.644164e-01
# 64 Omicron (BA.1*)  Eastern Cape 7.322681e-04 3.553831e-04 NA  3.572998e-05 1.428806e-03
# 65            Beta  Eastern Cape 1.916332e-14 1.166214e-14 NA -3.694067e-15 4.202070e-14
# 66           Alpha  Eastern Cape 7.422875e-15 5.638768e-15 NA -3.628906e-15 1.847466e-14
# 67           Other  Eastern Cape 4.286656e-16 2.640735e-16 NA -8.890902e-17 9.462402e-16
# 68           Delta  Eastern Cape 6.102616e-10 3.569503e-10 NA -8.934821e-11 1.309871e-09
# 69  Omicron (BA.2)  Eastern Cape 2.183079e-01 9.674792e-02 NA  2.868545e-02 4.079303e-01
# 70  Omicron (BA.3)  Eastern Cape 3.123864e-06 3.125139e-06 NA -3.001297e-06 9.249025e-06
# 71  Omicron (BA.4)  Eastern Cape 7.789604e-01 9.892230e-02 NA  5.850763e-01 9.728445e-01
# 72  Omicron (BA.5)  Eastern Cape 1.996316e-03 1.964476e-02 NA -3.650670e-02 4.049933e-02
# 73 Omicron (BA.1*) Northern Cape 2.121999e-03 1.074621e-03 NA  1.578083e-05 4.228217e-03
# 74            Beta Northern Cape 2.971384e-13 1.854045e-13 NA -6.624777e-14 6.605245e-13
# 75           Alpha Northern Cape 5.831014e-14 4.848215e-14 NA -3.671314e-14 1.533334e-13
# 76           Other Northern Cape 7.367474e-15 4.654110e-15 NA -1.754414e-15 1.648936e-14
# 77           Delta Northern Cape 2.607739e-09 1.567370e-09 NA -4.642509e-10 5.679728e-09
# 78  Omicron (BA.2) Northern Cape 9.434781e-01 4.275574e-01 NA  1.054810e-01 1.781475e+00
# 79  Omicron (BA.3) Northern Cape 9.413326e-06 1.170306e-05 NA -1.352425e-05 3.235090e-05
# 80  Omicron (BA.4) Northern Cape 4.572675e-05 2.919140e-05 NA -1.148734e-05 1.029408e-04
# 81  Omicron (BA.5) Northern Cape 5.434478e-02 4.285434e-01 NA -7.855848e-01 8.942743e-01





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
fit_southafrica_multi_preds_byprovince$variant = factor(fit_southafrica_multi_preds_byprovince$variant, levels=levels_VARIANTS_plot)
muller_southafrica_province_mfit = ggplot(data=fit_southafrica_multi_preds_byprovince, 
                                 aes(x=date, y=prob, group=variant)) + 
  facet_wrap(~ province) +
  scale_y_continuous(expand=c(0,0)) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  annotate("rect", xmin=max(GISAID_sel$date_num)+1, 
           xmax=as.Date(date.to+30, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.1, fill="white") + # extrapolated part
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
fit_southafrica_multi_preds_withCI$variant = factor(fit_southafrica_multi_preds_withCI$variant, 
                                                    levels=levels_VARIANTS_plot)

# multinomial model predictions with confidence intervals by province (slower)
fit_southafrica_multi_preds_byprovince_withCI = data.frame(emmeans(fit1_southafrica_province_multi,
                                                        ~ variant,
                                                        by=c("date_num","province"),
                                                        at=list(date_num=seq(date.from, date.to, by=7)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_southafrica_multi_preds_byprovince_withCI$date = as.Date(fit_southafrica_multi_preds_byprovince_withCI$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_byprovince_withCI$variant = factor(fit_southafrica_multi_preds_byprovince_withCI$variant,
                                                               levels=levels_VARIANTS_plot)

# Muller plot overall for South Africa
muller_southafrica_mfit = ggplot(data=fit_southafrica_multi_preds_withCI, 
                                   aes(x=date, y=prob, group=variant)) + 
  scale_y_continuous(expand=c(0,0)) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  annotate("rect", xmin=max(GISAID_sel$date_num)+1, 
           xmax=as.Date(date.to+30, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.1, fill="white") + # extrapolated part
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
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweek_province,
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.1, 3), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
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
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweek_province,
             aes(x=date, y=100*prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\ngenotyped", trans="sqrt",
                        range=c(0.1, 3), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
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
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
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
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
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
k=40
fit_cases = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) +
                  WEEKDAY, # + # + offset(log(testspercase)),
                  # BANKHOLIDAY,
                # + s(tests_new, bs="cs", k=8, fx=F),
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
                                                                  date.to)#,+-
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
fit_southafrica_multi_preds$variant = factor(fit_southafrica_multi_preds$variant, levels=levels_VARIANTS_plot)

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
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA","(negative binomial fit to NICD case data, with correction for weekday effects;\nvariant frequencies based on multinomial spline fit to GISAID data)") + #  &\nvariable testing intensity; marginal means calculated at constant, maximal testing effort
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
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
      
