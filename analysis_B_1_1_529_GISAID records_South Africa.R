# ANALYSIS OF GROWTH ADVANTAGE OF OMICRON (B.1.1.529/BA.XX) IN SOUTH AFRICA BASED ON GISAID SEQUENCE DATA

# T. Wenseleers
# last update 20 APRIL 2022
    
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
d6 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_20_21_subm_mar_apr_2022.tsv"), col_types = cols(.default = "c")) 

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
  (((grepl("BA.4",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
    ((grepl("L452R", GISAID$aa_substitutions)&
            grepl("486V", GISAID$aa_substitutions)&
            grepl("11F", GISAID$aa_substitutions)&
            (!grepl("D3N",GISAID$aa_substitutions)) )))) ~ "Omicron (BA.4)",
  (((grepl("BA.5",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
    ((grepl("L452R",GISAID$aa_substitutions)&
       grepl("486V",GISAID$aa_substitutions)&
     (!grepl("11F", GISAID$aa_substitutions))&
       grepl("D3N",GISAID$aa_substitutions))))) ~ "Omicron (BA.5)",
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
# contrast    estimate          SE df   asymp.LCL  asymp.UCL       p.value
# 1           Beta - (Omicron (BA.1*)) 0.040336695 0.002454835 NA  0.03552531 0.04514808  1.137484e-60
# 2          Alpha - (Omicron (BA.1*)) 0.003486441 0.009786234 NA -0.01569423 0.02266711  7.216460e-01
# 3          Other - (Omicron (BA.1*)) 0.115794720 0.001891357 NA  0.11208773 0.11950171  0.000000e+00
# 4          Delta - (Omicron (BA.1*)) 0.056230708 0.003048785 NA  0.05025520 0.06220622  5.865176e-76
# 5 Omicron (BA.2) - (Omicron (BA.1*)) 0.047628759 0.003564553 NA  0.04064236 0.05461515  1.011199e-40
# 6 Omicron (BA.3) - (Omicron (BA.1*)) 0.059494796 0.002649594 NA  0.05430169 0.06468791 1.161598e-111
# 7 Omicron (BA.4) - (Omicron (BA.1*)) 0.151708857 0.006707136 NA  0.13856311 0.16485460 2.816585e-113
# 8 Omicron (BA.5) - (Omicron (BA.1*)) 0.179976039 0.012802632 NA  0.15488334 0.20506874  6.905893e-45

# corresponding transmission advantage of Omicron BA.5 over BA.1* (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_avg[8,5]*2.2) # 1.40x
exp(delta_r_southafrica_avg[8,2]*2.2) # Omicron 1.49x transmission advantage over BA.1*
exp(delta_r_southafrica_avg[8,6]*2.2) # 1.57x


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
# contrast    estimate          SE df   asymp.LCL   asymp.UCL       p.value
# 1           Beta - (Omicron (BA.1*)) -0.13419510 0.002451514 NA -0.13899997 -0.12939022  0.000000e+00
# 2          Alpha - (Omicron (BA.1*)) -0.12446764 0.002568730 NA -0.12950225 -0.11943302  0.000000e+00
# 3          Other - (Omicron (BA.1*)) -0.14121884 0.002457287 NA -0.14603503 -0.13640264  0.000000e+00
# 4          Delta - (Omicron (BA.1*)) -0.10183039 0.002394044 NA -0.10652263 -0.09713815  0.000000e+00
# 5 Omicron (BA.2) - (Omicron (BA.1*))  0.07642427 0.001463509 NA  0.07355584  0.07929269  0.000000e+00
# 6 Omicron (BA.3) - (Omicron (BA.1*))  0.03743688 0.003178189 NA  0.03120774  0.04366601  4.989786e-32
# 7 Omicron (BA.4) - (Omicron (BA.1*))  0.19260087 0.008397668 NA  0.17614174  0.20906000 2.078064e-116
# 8 Omicron (BA.5) - (Omicron (BA.1*))  0.18422093 0.013507597 NA  0.15774653  0.21069534  2.369553e-42

# corresponding transmission advantage of Omicron BA.5 over BA.1* (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_province_avg[8,5]*2.2) # 1.42x
exp(delta_r_southafrica_province_avg[8,2]*2.2) # Omicron 1.50x transmission advantage over BA.1*
exp(delta_r_southafrica_province_avg[8,6]*2.2) # 1.59x


emtrsouthafrica_province_avg_pairw = emtrends(fit1_southafrica_province_multi, pairwise ~ variant,  
                                        var="date_num",  mode="latent",
                                        at=list(date_num=seq(today_num, today_num, by=1)),
                                        adjust="none", df=NA)
delta_r_southafrica_province_avg_pairw = data.frame(confint(emtrsouthafrica_province_avg_pairw, 
                                                      adjust="none", df=NA)$contrasts, 
                                              p.value=as.data.frame(emtrsouthafrica_province_avg_pairw$contrasts,
                                                                    adjust="none", df=NA)$p.value)
delta_r_southafrica_province_avg_pairw
#                         contrast     estimate           SE df    asymp.LCL    asymp.UCL       p.value
# 1            (Omicron (BA.1*)) - Beta  0.134195095 0.0024515140 NA  0.129390216  0.138999975  0.000000e+00
# 2           (Omicron (BA.1*)) - Alpha  0.124467635 0.0025687298 NA  0.119433017  0.129502253  0.000000e+00
# 3           (Omicron (BA.1*)) - Other  0.141218835 0.0024572874 NA  0.136402640  0.146035030  0.000000e+00
# 4           (Omicron (BA.1*)) - Delta  0.101830392 0.0023940436 NA  0.097138153  0.106522631  0.000000e+00
# 5  (Omicron (BA.1*)) - Omicron (BA.2) -0.076424267 0.0014635092 NA -0.079292693 -0.073555842  0.000000e+00
# 6  (Omicron (BA.1*)) - Omicron (BA.3) -0.037436880 0.0031781886 NA -0.043666015 -0.031207744  4.989786e-32
# 7  (Omicron (BA.1*)) - Omicron (BA.4) -0.192600870 0.0083976679 NA -0.209059996 -0.176141743 2.078064e-116
# 8  (Omicron (BA.1*)) - Omicron (BA.5) -0.184220935 0.0135075969 NA -0.210695338 -0.157746531  2.369553e-42
# 9                        Beta - Alpha -0.009727460 0.0008231339 NA -0.011340773 -0.008114148  3.166316e-32
# 10                       Beta - Other  0.007023740 0.0001789758 NA  0.006672954  0.007374526  0.000000e+00
# 11                       Beta - Delta -0.032364703 0.0005162886 NA -0.033376610 -0.031352796  0.000000e+00
# 12              Beta - Omicron (BA.2) -0.210619363 0.0028678390 NA -0.216240224 -0.204998502  0.000000e+00
# 13              Beta - Omicron (BA.3) -0.171631975 0.0040219644 NA -0.179514880 -0.163749070  0.000000e+00
# 14              Beta - Omicron (BA.4) -0.326795965 0.0087524483 NA -0.343950449 -0.309641482 4.020505e-305
# 15              Beta - Omicron (BA.5) -0.318416030 0.0137309605 NA -0.345328218 -0.291503842 5.792576e-119
# 16                      Alpha - Other  0.016751200 0.0008358665 NA  0.015112932  0.018389468  2.441935e-89
# 17                      Alpha - Delta -0.022637243 0.0009238774 NA -0.024448010 -0.020826477 1.391557e-132
# 18             Alpha - Omicron (BA.2) -0.200891903 0.0029686616 NA -0.206710372 -0.195073433  0.000000e+00
# 19             Alpha - Omicron (BA.3) -0.161904515 0.0040944658 NA -0.169929520 -0.153879509  0.000000e+00
# 20             Alpha - Omicron (BA.4) -0.317068505 0.0087860004 NA -0.334288749 -0.299848260 3.508361e-285
# 21             Alpha - Omicron (BA.5) -0.308688570 0.0137523717 NA -0.335642723 -0.281734416 1.393714e-111
# 22                      Other - Delta -0.039388443 0.0005430940 NA -0.040452888 -0.038323998  0.000000e+00
# 23             Other - Omicron (BA.2) -0.217643103 0.0028727758 NA -0.223273640 -0.212012565  0.000000e+00
# 24             Other - Omicron (BA.3) -0.178655715 0.0040254860 NA -0.186545522 -0.170765907  0.000000e+00
# 25             Other - Omicron (BA.4) -0.333819705 0.0087540672 NA -0.350977361 -0.316662049  0.000000e+00
# 26             Other - Omicron (BA.5) -0.325439770 0.0137319924 NA -0.352353980 -0.298525559 3.658756e-124
# 27             Delta - Omicron (BA.2) -0.178254659 0.0028188686 NA -0.183779540 -0.172729778  0.000000e+00
# 28             Delta - Omicron (BA.3) -0.139267272 0.0039871963 NA -0.147082033 -0.131452511 2.734601e-267
# 29             Delta - Omicron (BA.4) -0.294431262 0.0087365251 NA -0.311554536 -0.277307987 5.553925e-249
# 30             Delta - Omicron (BA.5) -0.286051327 0.0137208162 NA -0.312943632 -0.259159021  1.589881e-96
# 31    Omicron (BA.2) - Omicron (BA.3)  0.038987388 0.0032034785 NA  0.032708685  0.045266090  4.472757e-34
# 32    Omicron (BA.2) - Omicron (BA.4) -0.116176602 0.0082542473 NA -0.132354630 -0.099998575  5.428401e-45 ##
# 33    Omicron (BA.2) - Omicron (BA.5) -0.107796667 0.0134183621 NA -0.134096174 -0.081497161  9.471598e-16 ##
# 34    Omicron (BA.3) - Omicron (BA.4) -0.155163990 0.0088660686 NA -0.172541165 -0.137786815  1.410763e-68
# 35    Omicron (BA.3) - Omicron (BA.5) -0.146784055 0.0138028086 NA -0.173837063 -0.119731047  2.062367e-26
# 36    Omicron (BA.4) - Omicron (BA.5)  0.008379935 0.0151514918 NA -0.021316443  0.038076314  5.802110e-01


# fitted prop of different variantS today (on average across provinces)
# based on model with province included as a factor
multinom_preds_today_avg_province = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|1,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg_province
#           variant         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1*) 5.888078e-04 1.141869e-04 NA  3.650055e-04 8.126101e-04
# 2            Beta 2.290904e-15 1.092715e-15 NA  1.492217e-16 4.432587e-15
# 3           Alpha 1.017004e-15 5.534660e-16 NA -6.776976e-17 2.101777e-15
# 4           Other 4.773828e-17 2.292540e-17 NA  2.805325e-18 9.267124e-17
# 5           Delta 5.564542e-11 2.380463e-11 NA  8.989214e-12 1.023016e-10
# 6  Omicron (BA.2) 3.783219e-01 4.650397e-02 NA  2.871758e-01 4.694680e-01
# 7  Omicron (BA.3) 4.893913e-04 2.040835e-04 NA  8.939489e-05 8.893876e-04
# 8  Omicron (BA.4) 5.179112e-01 4.777973e-02 NA  4.242647e-01 6.115578e-01
# 9  Omicron (BA.5) 1.026887e-01 1.580568e-02 NA  7.171012e-02 1.336673e-01

# fitted prop of different variantS today (by province)
multinom_preds_today_byprovince = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|province,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_byprovince
#            variant      province         prob           SE df     asymp.LCL    asymp.UCL
# 1  Omicron (BA.1*)       Gauteng 3.446148e-05 1.056688e-05 NA  1.375078e-05 5.517219e-05
# 2             Beta       Gauteng 2.252385e-17 1.174041e-17 NA -4.869268e-19 4.553463e-17
# 3            Alpha       Gauteng 5.472549e-17 3.193373e-17 NA -7.863468e-18 1.173145e-16
# 4            Other       Gauteng 5.148416e-19 2.717363e-19 NA -1.775185e-20 1.047435e-18
# 5            Delta       Gauteng 8.914512e-13 4.412503e-13 NA  2.661653e-14 1.756286e-12
# 6   Omicron (BA.2)       Gauteng 4.518788e-02 1.175804e-02 NA  2.214255e-02 6.823322e-02
# 7   Omicron (BA.3)       Gauteng 4.786061e-05 2.207243e-05 NA  4.599440e-06 9.112177e-05
# 8   Omicron (BA.4)       Gauteng 9.006208e-01 3.122101e-02 NA  8.394288e-01 9.618129e-01
# 9   Omicron (BA.5)       Gauteng 5.410899e-02 2.602822e-02 NA  3.094619e-03 1.051234e-01
# 10 Omicron (BA.1*)    North West 3.660349e-04 1.831017e-04 NA  7.162219e-06 7.249076e-04
# 11            Beta    North West 3.138486e-16 2.117826e-16 NA -1.012376e-16 7.289348e-16
# 12           Alpha    North West 7.620646e-16 5.626215e-16 NA -3.406533e-16 1.864782e-15
# 13           Other    North West 2.697547e-17 1.831968e-17 NA -8.930441e-18 6.288139e-17
# 14           Delta    North West 1.709412e-11 1.114999e-11 NA -4.759464e-12 3.894771e-11
# 15  Omicron (BA.2)    North West 2.130873e-01 9.890272e-02 NA  1.924155e-02 4.069331e-01
# 16  Omicron (BA.3)    North West 1.102631e-03 6.602356e-04 NA -1.914066e-04 2.396669e-03
# 17  Omicron (BA.4)    North West 7.843305e-01 1.004989e-01 NA  5.873562e-01 9.813047e-01
# 18  Omicron (BA.5)    North West 1.113523e-03 1.288591e-02 NA -2.414239e-02 2.636944e-02
# 19 Omicron (BA.1*)    Mpumalanga 2.027450e-04 1.030463e-04 NA  7.779674e-07 4.047120e-04
# 20            Beta    Mpumalanga 4.426738e-16 3.009278e-16 NA -1.471338e-16 1.032481e-15
# 21           Alpha    Mpumalanga 8.848127e-17 8.130637e-17 NA -7.087628e-17 2.478388e-16
# 22           Other    Mpumalanga 1.586389e-17 1.089053e-17 NA -5.481166e-18 3.720894e-17
# 23           Delta    Mpumalanga 8.352517e-12 5.496580e-12 NA -2.420582e-12 1.912562e-11
# 24  Omicron (BA.2)    Mpumalanga 2.868534e-01 1.374921e-01 NA  1.737385e-02 5.563329e-01
# 25  Omicron (BA.3)    Mpumalanga 1.096353e-03 6.455560e-04 NA -1.689140e-04 2.361619e-03
# 26  Omicron (BA.4)    Mpumalanga 7.115877e-01 1.382943e-01 NA  4.405359e-01 9.826394e-01
# 27  Omicron (BA.5)    Mpumalanga 2.598425e-04 6.631760e-03 NA -1.273817e-02 1.325785e-02
# 28 Omicron (BA.1*)       Limpopo 4.014037e-05 2.602183e-05 NA -1.086147e-05 9.114221e-05
# 29            Beta       Limpopo 3.232851e-18 2.841943e-18 NA -2.337255e-18 8.802958e-18
# 30           Alpha       Limpopo 3.690130e-18 3.577737e-18 NA -3.322105e-18 1.070237e-17
# 31           Other       Limpopo 2.644347e-19 2.335111e-19 NA -1.932386e-19 7.221080e-19
# 32           Delta       Limpopo 2.747015e-13 2.363434e-13 NA -1.885230e-13 7.379260e-13
# 33  Omicron (BA.2)       Limpopo 1.089245e-01 6.792391e-02 NA -2.420391e-02 2.420529e-01
# 34  Omicron (BA.3)       Limpopo 8.478297e-05 6.386562e-05 NA -4.039135e-05 2.099573e-04
# 35  Omicron (BA.4)       Limpopo 8.908435e-01 6.818818e-02 NA  7.571972e-01 1.024490e+00
# 36  Omicron (BA.5)       Limpopo 1.070322e-04 4.698192e-03 NA -9.101255e-03 9.315319e-03
# 37 Omicron (BA.1*)  Western Cape 1.131368e-03 2.340224e-04 NA  6.726923e-04 1.590043e-03
# 38            Beta  Western Cape 1.459281e-15 6.689626e-16 NA  1.481383e-16 2.770424e-15
# 39           Alpha  Western Cape 2.443205e-15 1.316063e-15 NA -1.362309e-16 5.022641e-15
# 40           Other  Western Cape 2.297820e-17 1.075286e-17 NA  1.902980e-18 4.405341e-17
# 41           Delta  Western Cape 8.484834e-11 3.629544e-11 NA  1.371059e-11 1.559861e-10
# 42  Omicron (BA.2)  Western Cape 4.810898e-01 7.024648e-02 NA  3.434092e-01 6.187704e-01
# 43  Omicron (BA.3)  Western Cape 4.725328e-04 2.063747e-04 NA  6.804579e-05 8.770199e-04
# 44  Omicron (BA.4)  Western Cape 4.268828e-01 7.794219e-02 NA  2.741190e-01 5.796467e-01
# 45  Omicron (BA.5)  Western Cape 9.042347e-02 4.805182e-02 NA -3.756373e-03 1.846033e-01
# 46 Omicron (BA.1*)    Free State 1.362159e-03 2.786459e-04 NA  8.160228e-04 1.908295e-03
# 47            Beta    Free State 6.042541e-15 3.125338e-15 NA -8.300837e-17 1.216809e-14
# 48           Alpha    Free State 2.634305e-15 1.791393e-15 NA -8.767614e-16 6.145371e-15
# 49           Other    Free State 1.171730e-16 6.188968e-17 NA -4.128536e-18 2.384746e-16
# 50           Delta    Free State 1.131396e-10 5.519369e-11 NA  4.961994e-12 2.213173e-10
# 51  Omicron (BA.2)    Free State 9.980711e-01 4.798545e-04 NA  9.971306e-01 9.990116e-01
# 52  Omicron (BA.3)    Free State 2.971752e-04 3.119930e-04 NA -3.143198e-04 9.086701e-04
# 53  Omicron (BA.4)    Free State 6.266249e-11 2.056279e-11 NA  2.236017e-11 1.029648e-10
# 54  Omicron (BA.5)    Free State 2.695560e-04 1.378948e-04 NA -7.128706e-07 5.398248e-04
# 55 Omicron (BA.1*) KwaZulu Natal 1.613631e-04 5.801676e-05 NA  4.765231e-05 2.750738e-04
# 56            Beta KwaZulu Natal 3.155401e-16 1.709469e-16 NA -1.950968e-17 6.505898e-16
# 57           Alpha KwaZulu Natal 1.207482e-16 8.404136e-17 NA -4.396987e-17 2.854662e-16
# 58           Other KwaZulu Natal 1.009043e-17 5.537028e-18 NA -7.619433e-19 2.094281e-17
# 59           Delta KwaZulu Natal 1.607980e-11 8.217273e-12 NA -2.575722e-14 3.218536e-11
# 60  Omicron (BA.2) KwaZulu Natal 1.181846e-01 3.657346e-02 NA  4.650195e-02 1.898673e-01
# 61  Omicron (BA.3) KwaZulu Natal 8.241500e-05 5.093415e-05 NA -1.741410e-05 1.822441e-04
# 62  Omicron (BA.4) KwaZulu Natal 1.051197e-01 4.820095e-02 NA  1.064755e-02 1.995918e-01
# 63  Omicron (BA.5) KwaZulu Natal 7.764519e-01 6.923234e-02 NA  6.407590e-01 9.121448e-01
# 64 Omicron (BA.1*)  Eastern Cape 4.615070e-04 1.970260e-04 NA  7.534310e-05 8.476708e-04
# 65            Beta  Eastern Cape 6.950208e-16 4.220201e-16 NA -1.321234e-16 1.522165e-15
# 66           Alpha  Eastern Cape 3.289435e-16 2.507133e-16 NA -1.624455e-16 8.203325e-16
# 67           Other  Eastern Cape 1.217505e-17 7.491739e-18 NA -2.508489e-18 2.685859e-17
# 68           Delta  Eastern Cape 4.825358e-11 2.797029e-11 NA -6.567180e-12 1.030743e-10
# 69  Omicron (BA.2)  Eastern Cape 1.692338e-01 6.472252e-02 NA  4.237995e-02 2.960876e-01
# 70  Omicron (BA.3)  Eastern Cape 1.682572e-04 1.200638e-04 NA -6.706363e-05 4.035780e-04
# 71  Omicron (BA.4)  Eastern Cape 8.299807e-01 6.507405e-02 NA  7.024380e-01 9.575235e-01
# 72  Omicron (BA.5)  Eastern Cape 1.557280e-04 3.559887e-03 NA -6.821522e-03 7.132978e-03
# 73 Omicron (BA.1*) Northern Cape 1.539492e-03 5.718873e-04 NA  4.186136e-04 2.660370e-03
# 74            Beta Northern Cape 1.132348e-14 6.452324e-15 NA -1.322846e-15 2.396980e-14
# 75           Alpha Northern Cape 2.716870e-15 2.162493e-15 NA -1.521538e-15 6.955278e-15
# 76           Other Northern Cape 2.236092e-16 1.294466e-16 NA -3.010148e-17 4.773199e-16
# 77           Delta Northern Cape 2.118747e-10 1.145176e-10 NA -1.257567e-11 4.363250e-10
# 78  Omicron (BA.2) Northern Cape 9.842645e-01 2.912567e-01 NA  4.134119e-01 1.555117e+00
# 79  Omicron (BA.3) Northern Cape 1.052514e-03 7.810652e-04 NA -4.783459e-04 2.583373e-03
# 80  Omicron (BA.4) Northern Cape 1.183535e-02 2.817086e-01 NA -5.403033e-01 5.639740e-01
# 81  Omicron (BA.5) Northern Cape 1.308106e-03 7.920553e-02 NA -1.539319e-01 1.565481e-01





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
           xmax=as.Date(date.to+30, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.2, fill="white") + # extrapolated part
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

fit_cases_testing = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) +
                  WEEKDAY + # + # + offset(log(testspercase)),
                # BANKHOLIDAY,
                s(tests_new, bs="cs", k=8, fx=F),
                family=nb(), data=cases_tot,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+1*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-1*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 

# calculate average instantaneous growth rates & 95% CLs & Re values using emtrends ####
# based on the slope of the GAM fit on a log link scale
date.from = as.numeric(as.Date("2020-03-14")) # as.numeric(min(GISAID_sel$date_num))
date.to = today_num+7 # max(GISAID_sel$date_num)+extrapolate

avg_r_cases = as.data.frame(emtrends(fit_cases_testing, ~ date_num, var="date_num", 
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
  geom_line() + theme_hc() + xlab("Date of infection") + ylab("Week-on-week growth in infections (%)") + 
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
  geom_line() + theme_hc() + xlab("Date of infection") + ylab("Week-on-week growth in infections (%)") + 
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
cases_emmeans = as.data.frame(emmeans(fit_cases_testing, ~ date_num, at=list(date_num=seq(date.from, date.to, by=1),
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA","(negative binomial fit to NICD case data, with correction for weekday effects &\nvariable testing intensity; marginal means calculated at constant, maximal testing effort\nvariant frequencies based on multinomial spline fit to GISAID data)") + #  
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
      
