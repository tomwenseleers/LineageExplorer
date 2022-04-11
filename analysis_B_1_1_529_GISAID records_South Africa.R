# ANALYSIS OF GROWTH ADVANTAGE OF OMICRON (B.1.1.529/BA.XX) IN SOUTH AFRICA BASED ON GISAID SEQUENCE DATA

# T. Wenseleers
# last update 11 APRIL 2022
    
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
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.1*)")] = "red2" # "magenta"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.2)")] = "red3"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.3)")] = "green4" 
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
d3 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_01_30_19_subm_june_aug_2021.tsv"), col_types = cols(.default = "c")) 
d4 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_16_subm_sept_2021_dec_2021.tsv"), col_types = cols(.default = "c")) 
d5 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_02_04_15_subm_jan_feb_2022.tsv"), col_types = cols(.default = "c")) 
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
                                                   "Northern Cape", "Northern Cape", NA,
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
#                             contrast    estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1           Beta - (Omicron (BA.1*)) -0.05935633 0.003020570 NA -0.065276540 -0.053436123 5.701345e-86
# 2          Alpha - (Omicron (BA.1*)) -0.05977602 0.006650182 NA -0.072810136 -0.046741901 2.503338e-19
# 3          Other - (Omicron (BA.1*))  0.01468573 0.002448623 NA  0.009886521  0.019484946 2.003183e-09
# 4          Delta - (Omicron (BA.1*)) -0.05368651 0.002945835 NA -0.059460241 -0.047912780 3.295918e-74
# 5 Omicron (BA.2) - (Omicron (BA.1*))  0.03244221 0.001675215 NA  0.029158849  0.035725572 1.494430e-83
# 6 Omicron (BA.3) - (Omicron (BA.1*)) -0.01024168 0.005846031 NA -0.021699691  0.001216329 7.979044e-02
# 7 Omicron (BA.4) - (Omicron (BA.1*))  0.13203043 0.009435617 NA  0.113536959  0.150523900 1.725567e-44
# 8 Omicron (BA.5) - (Omicron (BA.1*))  0.21246507 0.020683275 NA  0.171926597  0.253003544 9.392533e-25

# corresponding transmission advantage of Omicron BA.5 over BA.1* (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_avg[8,5]*2.2) # 1.45x
exp(delta_r_southafrica_avg[8,2]*2.2) # Omicron 1.59x transmission advantage over Delta
exp(delta_r_southafrica_avg[8,6]*2.2) # 1.74x


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
#                             contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1           Beta - (Omicron (BA.1*)) -0.124401014 0.002201763 NA -0.128716391 -0.12008564 0.000000e+00
# 2          Alpha - (Omicron (BA.1*)) -0.115591118 0.002329019 NA -0.120155912 -0.11102632 0.000000e+00
# 3          Other - (Omicron (BA.1*)) -0.131693448 0.002208757 NA -0.136022533 -0.12736436 0.000000e+00
# 4          Delta - (Omicron (BA.1*)) -0.093057472 0.002139603 NA -0.097251016 -0.08886393 0.000000e+00
# 5 Omicron (BA.2) - (Omicron (BA.1*))  0.078609315 0.001704293 NA  0.075268962  0.08194967 0.000000e+00
# 6 Omicron (BA.3) - (Omicron (BA.1*))  0.004984473 0.005246546 NA -0.005298569  0.01526752 3.420877e-01
# 7 Omicron (BA.4) - (Omicron (BA.1*))  0.189641851 0.012787962 NA  0.164577905  0.21470580 9.413049e-50
# 8 Omicron (BA.5) - (Omicron (BA.1*))  0.206291742 0.021218861 NA  0.164703540  0.24787995 2.427392e-22

# corresponding transmission advantage of Omicron BA.5 over BA.1* (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_province_avg[8,5]*2.2) # 1.44x
exp(delta_r_southafrica_province_avg[8,2]*2.2) # Omicron 1.57x transmission advantage over Delta
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
# 1            (Omicron (BA.1*)) - Beta  0.124401014 0.0022017633 NA  0.120085637  0.128716391  0.000000e+00
# 2           (Omicron (BA.1*)) - Alpha  0.115591118 0.0023290194 NA  0.111026324  0.120155912  0.000000e+00
# 3           (Omicron (BA.1*)) - Other  0.131693448 0.0022087574 NA  0.127364363  0.136022533  0.000000e+00
# 4           (Omicron (BA.1*)) - Delta  0.093057472 0.0021396026 NA  0.088863928  0.097251016  0.000000e+00
# 5  (Omicron (BA.1*)) - Omicron (BA.2) -0.078609315 0.0017042932 NA -0.081949669 -0.075268962  0.000000e+00
# 6  (Omicron (BA.1*)) - Omicron (BA.3) -0.004984473 0.0052465465 NA -0.015267515  0.005298569  3.420877e-01
# 7  (Omicron (BA.1*)) - Omicron (BA.4) -0.189641851 0.0127879625 NA -0.214705797 -0.164577905  9.413049e-50
# 8  (Omicron (BA.1*)) - Omicron (BA.5) -0.206291742 0.0212188607 NA -0.247879945 -0.164703540  2.427392e-22
# 9                        Beta - Alpha -0.008809896 0.0008114085 NA -0.010400227 -0.007219565  1.836390e-27
# 10                       Beta - Other  0.007292434 0.0001854635 NA  0.006928932  0.007655936  0.000000e+00
# 11                       Beta - Delta -0.031343542 0.0005069229 NA -0.032337093 -0.030349992  0.000000e+00
# 12              Beta - Omicron (BA.2) -0.203010329 0.0027976617 NA -0.208493646 -0.197527013  0.000000e+00
# 13              Beta - Omicron (BA.3) -0.129385487 0.0056533838 NA -0.140465916 -0.118305058 6.350015e-116
# 14              Beta - Omicron (BA.4) -0.314042865 0.0129790514 NA -0.339481338 -0.288604391 2.443289e-129
# 15              Beta - Omicron (BA.5) -0.330692756 0.0213345647 NA -0.372507735 -0.288877778  3.451243e-54
# 16                      Alpha - Other  0.016102330 0.0008250100 NA  0.014485340  0.017719320  7.759733e-85
# 17                      Alpha - Delta -0.022533646 0.0009122135 NA -0.024321552 -0.020745741 1.013737e-134
# 18             Alpha - Omicron (BA.2) -0.194200434 0.0028988756 NA -0.199882125 -0.188518742  0.000000e+00
# 19             Alpha - Omicron (BA.3) -0.120575591 0.0057041438 NA -0.131755508 -0.109395675  3.539954e-99
# 20             Alpha - Omicron (BA.4) -0.305232969 0.0130012440 NA -0.330714939 -0.279750999 6.972664e-122
# 21             Alpha - Omicron (BA.5) -0.321882860 0.0213480730 NA -0.363724315 -0.280041406  2.265553e-51
# 22                      Other - Delta -0.038635976 0.0005365544 NA -0.039687603 -0.037584349  0.000000e+00
# 23             Other - Omicron (BA.2) -0.210302763 0.0028031694 NA -0.215796874 -0.204808652  0.000000e+00
# 24             Other - Omicron (BA.3) -0.136677921 0.0056561115 NA -0.147763696 -0.125592146 5.238746e-129
# 25             Other - Omicron (BA.4) -0.321335299 0.0129802397 NA -0.346776101 -0.295894496 2.689757e-135
# 26             Other - Omicron (BA.5) -0.337985190 0.0213352876 NA -0.379801586 -0.296168795  1.606681e-56
# 27             Delta - Omicron (BA.2) -0.171666787 0.0027490067 NA -0.177054742 -0.166278833  0.000000e+00
# 28             Delta - Omicron (BA.3) -0.098041945 0.0056294579 NA -0.109075480 -0.087008410  6.252382e-68
# 29             Delta - Omicron (BA.4) -0.282699323 0.0129686507 NA -0.308117411 -0.257281234 2.388630e-105
# 30             Delta - Omicron (BA.5) -0.299349214 0.0213282390 NA -0.341151795 -0.257546634  9.473294e-45
# 31    Omicron (BA.2) - Omicron (BA.3)  0.073624842 0.0054116753 NA  0.063018154  0.084231531  3.749152e-42
# 32    Omicron (BA.2) - Omicron (BA.4) -0.111032535 0.0126590615 NA -0.135843840 -0.086221231  1.770980e-18
# 33    Omicron (BA.2) - Omicron (BA.5) -0.127682427 0.0211400763 NA -0.169116215 -0.086248639  1.542785e-09
# 34    Omicron (BA.3) - Omicron (BA.4) -0.184657378 0.0137796042 NA -0.211664906 -0.157649849  5.983496e-41
# 35    Omicron (BA.3) - Omicron (BA.5) -0.201307269 0.0218303009 NA -0.244093873 -0.158520666  2.930792e-20
# 36    Omicron (BA.4) - Omicron (BA.5) -0.016649892 0.0241058628 NA -0.063896515  0.030596731  4.897548e-01


# fitted prop of different variantS today (on average across provinces)
# based on model with province included as a factor
multinom_preds_today_avg_province = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|1,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg_province
#           variant         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1*) 1.633749e-03 2.819449e-04 NA  1.081147e-03 2.186350e-03
# 2            Beta 8.360751e-14 3.437107e-14 NA  1.624146e-14 1.509736e-13
# 3           Alpha 3.210500e-14 1.526967e-14 NA  2.176995e-15 6.203300e-14
# 4           Other 1.797550e-15 7.317809e-16 NA  3.632862e-16 3.231815e-15
# 5           Delta 1.292240e-09 4.737560e-10 NA  3.636951e-10 2.220785e-09
# 6  Omicron (BA.2) 5.960412e-01 4.802139e-02 NA  5.019210e-01 6.901614e-01
# 7  Omicron (BA.3) 3.305916e-05 2.083036e-05 NA -7.767591e-06 7.388591e-05
# 8  Omicron (BA.4) 2.993131e-01 4.520635e-02 NA  2.107103e-01 3.879159e-01
# 9  Omicron (BA.5) 1.029789e-01 2.572646e-02 NA  5.255596e-02 1.534018e-01

# fitted prop of different variantS today (by province)
multinom_preds_today_byprovince = data.frame(emmeans(fit1_southafrica_province_multi, ~ variant|province,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_byprovince
#            variant      province         prob           SE df     asymp.LCL    asymp.UCL
# 1  Omicron (BA.1*)       Gauteng 1.636194e-04 6.076136e-05 NA  4.452927e-05 2.827094e-04
# 2             Beta       Gauteng 1.670399e-15 8.749728e-16 NA -4.451633e-17 3.385314e-15
# 3            Alpha       Gauteng 2.772148e-15 1.616915e-15 NA -3.969472e-16 5.941243e-15
# 4            Other       Gauteng 3.682069e-17 1.953948e-17 NA -1.475986e-18 7.511738e-17
# 5            Delta       Gauteng 3.776009e-11 1.887654e-11 NA  7.627598e-13 7.475742e-11
# 6   Omicron (BA.2)       Gauteng 1.359245e-01 4.413317e-02 NA  4.942511e-02 2.224239e-01
# 7   Omicron (BA.3)       Gauteng 1.813174e-06 1.448734e-06 NA -1.026291e-06 4.652640e-06
# 8   Omicron (BA.4)       Gauteng 7.626942e-01 8.607558e-02 NA  5.939891e-01 9.313992e-01
# 9   Omicron (BA.5)       Gauteng 1.012159e-01 6.820190e-02 NA -3.245739e-02 2.348892e-01
# 10 Omicron (BA.1*)    North West 2.751333e-03 5.796821e-04 NA  1.615177e-03 3.887489e-03
# 11            Beta    North West 4.095672e-14 1.852950e-14 NA  4.639574e-15 7.727386e-14
# 12           Alpha    North West 7.224903e-14 3.898096e-14 NA -4.152238e-15 1.486503e-13
# 13           Other    North West 3.249099e-15 1.493006e-15 NA  3.228615e-16 6.175336e-15
# 14           Delta    North West 1.235949e-09 5.194105e-10 NA  2.179228e-10 2.253974e-09
# 15  Omicron (BA.2)    North West 9.971265e-01 6.029765e-04 NA  9.959446e-01 9.983083e-01
# 16  Omicron (BA.3)    North West 1.094920e-04 7.679216e-05 NA -4.101791e-05 2.600018e-04
# 17  Omicron (BA.4)    North West 6.822369e-39 3.015233e-39 NA  9.126218e-40 1.273212e-38
# 18  Omicron (BA.5)    North West 1.271349e-05 9.276201e-06 NA -5.467528e-06 3.089451e-05
# 19 Omicron (BA.1*)    Mpumalanga 8.648010e-04 2.827582e-04 NA  3.106052e-04 1.418997e-03
# 20            Beta    Mpumalanga 4.284651e-14 2.179859e-14 NA  1.220661e-16 8.557095e-14
# 21           Alpha    Mpumalanga 6.113068e-15 4.878583e-15 NA -3.448778e-15 1.567491e-14
# 22           Other    Mpumalanga 1.472607e-15 7.632525e-16 NA -2.334075e-17 2.968554e-15
# 23           Delta    Mpumalanga 4.369207e-10 2.104836e-10 NA  2.438040e-11 8.494610e-10
# 24  Omicron (BA.2)    Mpumalanga 7.548802e-01 2.010522e-01 NA  3.608251e-01 1.148935e+00
# 25  Omicron (BA.3)    Mpumalanga 1.151975e-04 7.872351e-05 NA -3.909772e-05 2.694928e-04
# 26  Omicron (BA.4)    Mpumalanga 2.432525e-01 2.010529e-01 NA -1.508039e-01 6.373088e-01
# 27  Omicron (BA.5)    Mpumalanga 8.873298e-04 1.925981e-02 NA -3.686120e-02 3.863586e-02
# 28 Omicron (BA.1*)       Limpopo 1.068906e-04 6.854592e-05 NA -2.745695e-05 2.412381e-04
# 29            Beta       Limpopo 2.577210e-16 2.122575e-16 NA -1.582960e-16 6.737381e-16
# 30           Alpha       Limpopo 2.230984e-16 2.046117e-16 NA -1.779331e-16 6.241299e-16
# 31           Other       Limpopo 1.679852e-17 1.393644e-17 NA -1.051640e-17 4.411344e-17
# 32           Delta       Limpopo 1.259520e-11 1.013286e-11 NA -7.264839e-12 3.245523e-11
# 33  Omicron (BA.2)       Limpopo 2.106116e-01 1.287045e-01 NA -4.164461e-02 4.628678e-01
# 34  Omicron (BA.3)       Limpopo 4.055948e-06 3.807149e-06 NA -3.405926e-06 1.151782e-05
# 35  Omicron (BA.4)       Limpopo 7.878043e-01 1.310753e-01 NA  5.309014e-01 1.044707e+00
# 36  Omicron (BA.5)       Limpopo 1.473102e-03 2.632738e-02 NA -5.012762e-02 5.307382e-02
# 37 Omicron (BA.1*)  Western Cape 3.837383e-03 7.247023e-04 NA  2.416992e-03 5.257773e-03
# 38            Beta  Western Cape 7.891120e-14 3.206733e-14 NA  1.606039e-14 1.417620e-13
# 39           Alpha  Western Cape 8.538676e-14 4.213327e-14 NA  2.807073e-15 1.679665e-13
# 40           Other  Western Cape 1.187064e-15 4.957987e-16 NA  2.153169e-16 2.158812e-15
# 41           Delta  Western Cape 2.458611e-09 9.209336e-10 NA  6.536144e-10 4.263608e-09
# 42  Omicron (BA.2)  Western Cape 8.649024e-01 8.667422e-02 NA  6.950240e-01 1.034781e+00
# 43  Omicron (BA.3)  Western Cape 4.653801e-05 3.158851e-05 NA -1.537434e-05 1.084504e-04
# 44  Omicron (BA.4)  Western Cape 1.310106e-01 8.699664e-02 NA -3.949969e-02 3.015209e-01
# 45  Omicron (BA.5)  Western Cape 2.031316e-04 4.340877e-03 NA -8.304831e-03 8.711095e-03
# 46 Omicron (BA.1*)    Free State 1.995033e-03 4.826559e-04 NA  1.049045e-03 2.941022e-03
# 47            Beta    Free State 1.363611e-13 6.716210e-14 NA  4.725842e-15 2.679964e-13
# 48           Alpha    Free State 4.208269e-14 2.774649e-14 NA -1.229942e-14 9.646481e-14
# 49           Other    Free State 2.633660e-15 1.328429e-15 NA  2.998752e-17 5.237332e-15
# 50           Delta    Free State 1.486969e-09 6.896839e-10 NA  1.352133e-10 2.838725e-09
# 51  Omicron (BA.2)    Free State 9.977288e-01 5.605296e-04 NA  9.966302e-01 9.988274e-01
# 52  Omicron (BA.3)    Free State 2.093739e-15 1.572222e-15 NA -9.877591e-16 5.175236e-15
# 53  Omicron (BA.4)    Free State 4.241672e-21 1.943756e-21 NA  4.319801e-22 8.051363e-21
# 54  Omicron (BA.5)    Free State 2.761823e-04 2.038370e-04 NA -1.233308e-04 6.756954e-04
# 55 Omicron (BA.1*) KwaZulu Natal 4.459398e-04 1.959164e-04 NA  6.195079e-05 8.299288e-04
# 56            Beta KwaZulu Natal 1.293011e-14 7.320988e-15 NA -1.418763e-15 2.727898e-14
# 57           Alpha KwaZulu Natal 3.444102e-15 2.457649e-15 NA -1.372800e-15 8.261005e-15
# 58           Other KwaZulu Natal 3.986998e-16 2.285505e-16 NA -4.925088e-17 8.466505e-16
# 59           Delta KwaZulu Natal 3.803827e-10 2.050095e-10 NA -2.142867e-11 7.821940e-10
# 60  Omicron (BA.2) KwaZulu Natal 1.706218e-01 6.573051e-02 NA  4.179237e-02 2.994512e-01
# 61  Omicron (BA.3) KwaZulu Natal 2.984631e-06 2.805743e-06 NA -2.514524e-06 8.483785e-06
# 62  Omicron (BA.4) KwaZulu Natal 1.618757e-02 1.718437e-02 NA -1.749317e-02 4.986831e-02
# 63  Omicron (BA.5) KwaZulu Natal 8.127417e-01 7.208353e-02 NA  6.714606e-01 9.540228e-01
# 64 Omicron (BA.1*)  Eastern Cape 1.500409e-03 7.263464e-04 NA  7.679673e-05 2.924022e-03
# 65            Beta  Eastern Cape 4.818319e-14 3.040244e-14 NA -1.140450e-14 1.077709e-13
# 66           Alpha  Eastern Cape 1.298117e-14 1.031765e-14 NA -7.241053e-15 3.320340e-14
# 67           Other  Eastern Cape 7.713578e-16 4.930630e-16 NA -1.950280e-16 1.737744e-15
# 68           Delta  Eastern Cape 1.840725e-09 1.116310e-09 NA -3.472019e-10 4.028651e-09
# 69  Omicron (BA.2)  Eastern Cape 2.455862e-01 1.040513e-01 NA  4.164932e-02 4.495230e-01
# 70  Omicron (BA.3)  Eastern Cape 2.757704e-27 2.374026e-27 NA -1.895303e-27 7.410710e-27
# 71  Omicron (BA.4)  Eastern Cape 7.528690e-01 1.047181e-01 NA  5.476253e-01 9.581127e-01
# 72  Omicron (BA.5)  Eastern Cape 4.441976e-05 2.962563e-03 NA -5.762098e-03 5.850938e-03
# 73 Omicron (BA.1*) Northern Cape 3.038328e-03 9.397304e-04 NA  1.196490e-03 4.880165e-03
# 74            Beta Northern Cape 3.903506e-13 1.910538e-13 NA  1.589213e-14 7.648091e-13
# 75           Alpha Northern Cape 6.369291e-14 4.696351e-14 NA -2.835388e-14 1.557397e-13
# 76           Other Northern Cape 6.411847e-15 3.215752e-15 NA  1.090881e-16 1.271461e-14
# 77           Delta Northern Cape 3.740247e-09 1.714349e-09 NA  3.801853e-10 7.100309e-09
# 78  Omicron (BA.2) Northern Cape 9.869887e-01 1.851015e-01 NA  6.241965e-01 1.349781e+00
# 79  Omicron (BA.3) Northern Cape 1.745119e-05 2.087796e-05 NA -2.346886e-05 5.837124e-05
# 80  Omicron (BA.4) Northern Cape 2.354694e-17 1.159249e-17 NA  8.260734e-19 4.626780e-17
# 81  Omicron (BA.5) Northern Cape 9.955487e-03 1.856728e-01 NA -3.539565e-01 3.738675e-01





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
      
