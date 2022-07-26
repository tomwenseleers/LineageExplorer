# ANALYSIS OF GROWTH ADVANTAGE OF OMICRON SUBVARIANTS IN SOUTH AFRICA BASED ON GISAID SEQUENCE DATA

# T. Wenseleers
# last update 15 MAY 2022
    
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
  sel_reference_VOC = "Omicron (BA.1)"
  levels_VARIANTS = c(sel_reference_VOC, "Beta", "Alpha", "Other", "Delta", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.1.1, with S:R346K)", "Omicron (BA.4)", sel_target_VOC)
  levels_VARIANTS_plot = c("Other", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.1.1, with S:R346K)", "Omicron (BA.4)", "Omicron (BA.5)")
  
  n = length(levels_VARIANTS)
  lineage_cols = hcl(h = seq(0, 180, length = n), l = 60, c = 180)
  lineage_cols[which(levels_VARIANTS=="Alpha")] = "#0085FF"
  lineage_cols[which(levels_VARIANTS=="Beta")] = "green4"
  lineage_cols[which(levels_VARIANTS=="Delta")] = "mediumorchid"
  # lineage_cols[which(levels_VARIANTS=="C.1.2")] = "darkorange"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.1)")] = "red" # "magenta"
  lineage_cols[which(levels_VARIANTS=="Omicron (BA.1.1, with S:R346K)")] = "black"
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
# downloaded by submission date
d1 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_15_subm_2020.tsv"), col_types = cols(.default = "c")) 
d2 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_16_subm_jan_may_2021.tsv"), col_types = cols(.default = "c")) 
d3 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_12_19_subm_june_aug_2021.tsv"), col_types = cols(.default = "c")) 
d4 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_11_16_subm_sept_2021_dec_2021.tsv"), col_types = cols(.default = "c")) 
d5 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_04_12_19_subm_jan_feb_2022.tsv"), col_types = cols(.default = "c")) 
d6 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2022_05_15_19_subm_mar_may_2022.tsv"), col_types = cols(.default = "c")) 

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
range(GISAID$date) # "2020-03-06" "2022-05-08"
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

# unique(GISAID$pango_lineage[grepl("^BA\\.1\\.1$|BA\\.1\\.1\\.", GISAID$pango_lineage)&grepl("R346K", GISAID$aa_substitutions)])

GISAID$variant = case_when(
  # grepl("N679K", aa_substitutions) & grepl("H655Y", aa_substitutions) & grepl("P681H", aa_substitutions) ~ "B.1.1.529",
  (grepl("^BA\\.1\\.1$|BA\\.1\\.1\\.", GISAID$pango_lineage)&grepl("R346K", GISAID$aa_substitutions)) ~ "Omicron (BA.1.1, with S:R346K)",
  (grepl("^BA\\.1$|BA\\.1\\.", GISAID$pango_lineage)) ~ "Omicron (BA.1)",
  (grepl("^BA\\.3$|BA\\.3\\.", GISAID$pango_lineage)) ~ "Omicron (BA.3)",
  (((grepl("^BA\\.4",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
                                             ((grepl("L452R", GISAID$aa_substitutions)&
                                                 grepl("486V", GISAID$aa_substitutions)&
                                                 grepl("11F", GISAID$aa_substitutions)&
                                                 (!grepl("D3N",GISAID$aa_substitutions)) )))) ~ "Omicron (BA.4)",
  (((grepl("^BA\\.5",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
                                             ((grepl("L452R",GISAID$aa_substitutions)&
                                                 grepl("486V",GISAID$aa_substitutions)&
                                                 (!grepl("11F", GISAID$aa_substitutions))&
                                                 grepl("D3N",GISAID$aa_substitutions))))) ~ "Omicron (BA.5)",
  (grepl("^BA\\.2",GISAID$pango_lineage)) ~ "Omicron (BA.2)",
  grepl("B.1.617.2", GISAID$pango_lineage, fixed=T) | grepl("AY", GISAID$pango_lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", GISAID$pango_lineage) ~ "Alpha",
  grepl("B.1.351", GISAID$pango_lineage, fixed=T) ~ "Beta",
  # grepl("C.1.2", GISAID$pango_lineage, fixed=T) ~ "C.1.2",
  T ~ "Other"
)

# GISAID$variant = case_when(
#   # grepl("N679K", aa_substitutions) & grepl("H655Y", aa_substitutions) & grepl("P681H", aa_substitutions) ~ "B.1.1.529",
#   (grepl("BA.1", GISAID$pango_lineage) ~ "Omicron (BA.1)"),
#   (GISAID$pango_lineage=="BA.3") ~ "Omicron (BA.3)",
#   (((grepl("BA.4",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
#     ((grepl("L452R", GISAID$aa_substitutions)&
#             grepl("486V", GISAID$aa_substitutions)&
#             grepl("11F", GISAID$aa_substitutions)&
#             (!grepl("D3N",GISAID$aa_substitutions)) )))) ~ "Omicron (BA.4)",
#   (((grepl("BA.5",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
#     ((grepl("L452R",GISAID$aa_substitutions)&
#        grepl("486V",GISAID$aa_substitutions)&
#      (!grepl("11F", GISAID$aa_substitutions))&
#        grepl("D3N",GISAID$aa_substitutions))))) ~ "Omicron (BA.5)",
#   (grepl("BA.2",GISAID$pango_lineage)) ~ "Omicron (BA.2)",
#   grepl("B.1.617.2", GISAID$pango_lineage, fixed=T) | grepl("AY", GISAID$pango_lineage)  ~ "Delta",
#   grepl("B.1.1.7", GISAID$pango_lineage, fixed=T) ~ "Alpha",
#   grepl("B.1.351", GISAID$pango_lineage, fixed=T) ~ "Beta",
#   # grepl("C.1.2", GISAID$pango_lineage, fixed=T) ~ "C.1.2",
#   T ~ "Other"
# )

# GISAID[GISAID$variant == "Omicron" & GISAID$date<as.Date("2021-10-01"),"variant"] = "Other" # miscoded date??

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

write.csv(GISAID[,"Accession ID"],
          paste0(".\\plots\\",plotdir,
                 "\\accession IDs South Africa.csv"), 
          row.names=F)


# ANALYSIS OF VOCs IN SOUTH AFRICA ####
GISAID_sel = GISAID
nrow(GISAID_sel) # 39541
sum(GISAID_sel$variant==sel_target_VOC) # 348
table(GISAID_sel$variant)
range(GISAID_sel$date)

GISAID_sel$variant = factor(GISAID_sel$variant)
GISAID_sel$variant = relevel(GISAID_sel$variant, ref=sel_reference_VOC) # we code Omicron (BA.2) as the reference

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
fit1_southafrica_province_multi = nnet::multinom(variant ~ scale(date_num)+province, weights=count, 
                                                 data=data_agbyweek_province, maxit=1000)
fit2_southafrica_province_multi = nnet::multinom(variant ~ ns(date_num, df=2)+province, weights=count, 
                                                 data=data_agbyweek_province, maxit=1000)
BIC(fit1_southafrica_province_multi, fit2_southafrica_province_multi) 
# fit2_southafrica_province_multi  best

# and also a multinomial 2 df spline fit to the data aggregated over all provinces
fit1_southafrica_multi = nnet::multinom(variant ~ scale(date_num), weights=count, 
                                        data=data_agbyday, maxit=1000)
fit2_southafrica_multi = nnet::multinom(variant ~ ns(date_num, df=2), weights=count, 
                                        data=data_agbyday, maxit=1000)
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
# contrast    estimate          SE df    asymp.LCL   asymp.UCL       p.value
# 1                           Beta - Omicron (BA.1)  0.02537294 0.003691644 NA  0.018137451 0.032608429  6.283327e-12
# 2                          Alpha - Omicron (BA.1) -0.01863739 0.011203273 NA -0.040595401 0.003320621  9.619907e-02
# 3                          Other - Omicron (BA.1)  0.10817242 0.003344060 NA  0.101618182 0.114726658 1.498122e-229
# 4                          Delta - Omicron (BA.1)  0.04082221 0.004279675 NA  0.032434199 0.049210217  1.447411e-21
# 5                 Omicron (BA.2) - Omicron (BA.1)  0.03639528 0.003962204 NA  0.028629498 0.044161053  4.091896e-20
# 6                 Omicron (BA.3) - Omicron (BA.1)  0.01671015 0.004176869 NA  0.008523641 0.024896668  6.317113e-05
# 7 Omicron (BA.1.1, with S:R346K) - Omicron (BA.1) -0.01174025 0.007147440 NA -0.025748974 0.002268476  1.004697e-01
# 8                 Omicron (BA.4) - Omicron (BA.1)  0.13013289 0.004146742 NA  0.122005423 0.138260351 3.566500e-216
# 9                 Omicron (BA.5) - Omicron (BA.1)  0.15693168 0.005937989 NA  0.145293432 0.168569920 6.460317e-154

# corresponding transmission advantage of Omicron BA.5 over BA.1 (ie how much higher effective reproduction number is at any timepoint) 
# here I am using a short 2.2 day generation time
exp(delta_r_southafrica_avg[9,5]*2.2) # 1.38x
exp(delta_r_southafrica_avg[9,2]*2.2) # Omicron 1.41x transmission advantage over BA.1
exp(delta_r_southafrica_avg[9,6]*2.2) # 1.45x


# avg growth rate advantage of Omicron BA.2 over BA.1 (difference in growth rate per day) (today)
# based on multinomial model with province included as factor
emtrsouthafrica_province_avg = emtrends(fit2_southafrica_province_multi, trt.vs.ctrl ~ variant,  
                   var="date_num",  mode="latent",
                   at=list(date_num=seq(today_num, today_num, by=1)),
                   adjust="none", df=NA)
delta_r_southafrica_province_avg = data.frame(confint(emtrsouthafrica_province_avg, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrsouthafrica_province_avg$contrasts,
                                               adjust="none", df=NA)$p.value)
delta_r_southafrica_province_avg
# contrast     estimate          SE df   asymp.LCL   asymp.UCL       p.value
# 1                           Beta - Omicron (BA.1)  0.045262977 0.003980578 NA  0.03746119 0.053064766  5.834414e-30
# 2                          Alpha - Omicron (BA.1) -0.003820912 0.011146960 NA -0.02566855 0.018026729  7.317669e-01
# 3                          Other - Omicron (BA.1)  0.128778664 0.003615548 NA  0.12169232 0.135865007 7.370226e-278
# 4                          Delta - Omicron (BA.1)  0.060239741 0.004573000 NA  0.05127683 0.069202657  1.256515e-39
# 5                 Omicron (BA.2) - Omicron (BA.1)  0.066649659 0.004183682 NA  0.05844979 0.074849524  3.869790e-57
# 6                 Omicron (BA.3) - Omicron (BA.1)  0.054878917 0.004463642 NA  0.04613034 0.063627495  9.677593e-35
# 7 Omicron (BA.1.1, with S:R346K) - Omicron (BA.1) -0.007550789 0.007930174 NA -0.02309364 0.007992066  3.410162e-01
# 8                 Omicron (BA.4) - Omicron (BA.1)  0.179214771 0.004502076 NA  0.17039086 0.188038679  0.000000e+00
# 9                 Omicron (BA.5) - Omicron (BA.1)  0.192432873 0.006479189 NA  0.17973390 0.205131850 7.645084e-194

# corresponding transmission advantage of Omicron BA.5 over BA.1 (ie how much higher effective reproduction number is at any timepoint) 
exp(delta_r_southafrica_province_avg[9,5]*2.2) # 1.49x
exp(delta_r_southafrica_province_avg[9,2]*2.2) # Omicron 1.53x transmission advantage over BA.1
exp(delta_r_southafrica_province_avg[9,6]*2.2) # 1.57x

# pairwase growth rate advantages for multinomial spline model today
emtrsouthafrica_province_avg_pairw = emtrends(fit2_southafrica_province_multi, pairwise ~ variant,  
                                        var="date_num",  mode="latent",
                                        at=list(date_num=seq(today_num, today_num, by=1)),
                                        adjust="none", df=NA)
delta_r_southafrica_province_avg_pairw = data.frame(confint(emtrsouthafrica_province_avg_pairw, 
                                                      adjust="none", df=NA)$contrasts, 
                                              p.value=as.data.frame(emtrsouthafrica_province_avg_pairw$contrasts,
                                                                    adjust="none", df=NA)$p.value)
delta_r_southafrica_province_avg_pairw
# contrast     estimate          SE df    asymp.LCL    asymp.UCL       p.value
# 1                            Omicron (BA.1) - Beta -0.045262977 0.003980578 NA -0.053064766 -0.037461188  5.834414e-30
# 2                           Omicron (BA.1) - Alpha  0.003820912 0.011146960 NA -0.018026729  0.025668552  7.317669e-01
# 3                           Omicron (BA.1) - Other -0.128778664 0.003615548 NA -0.135865007 -0.121692320 7.370226e-278
# 4                           Omicron (BA.1) - Delta -0.060239741 0.004573000 NA -0.069202657 -0.051276826  1.256515e-39
# 5                  Omicron (BA.1) - Omicron (BA.2) -0.066649659 0.004183682 NA -0.074849524 -0.058449794  3.869790e-57
# 6                  Omicron (BA.1) - Omicron (BA.3) -0.054878917 0.004463642 NA -0.063627495 -0.046130339  9.677593e-35
# 7  Omicron (BA.1) - Omicron (BA.1.1, with S:R346K)  0.007550789 0.007930174 NA -0.007992066  0.023093643  3.410162e-01
# 8                  Omicron (BA.1) - Omicron (BA.4) -0.179214771 0.004502076 NA -0.188038679 -0.170390864  0.000000e+00
# 9                  Omicron (BA.1) - Omicron (BA.5) -0.192432873 0.006479189 NA -0.205131850 -0.179733895 7.645084e-194
# 10                                    Beta - Alpha  0.049083889 0.010623045 NA  0.028263103  0.069904674  3.827973e-06
# 11                                    Beta - Other -0.083515687 0.001775643 NA -0.086995884 -0.080035490  0.000000e+00
# 12                                    Beta - Delta -0.014976764 0.003108367 NA -0.021069051 -0.008884477  1.448519e-06
# 13                           Beta - Omicron (BA.2) -0.021386682 0.003692502 NA -0.028623854 -0.014149510  6.958590e-09
# 14                           Beta - Omicron (BA.3) -0.009615940 0.007021894 NA -0.023378600  0.004146719  1.708672e-01
# 15           Beta - Omicron (BA.1.1, with S:R346K)  0.052813766 0.008405642 NA  0.036339011  0.069288521  3.318168e-10
# 16                           Beta - Omicron (BA.4) -0.133951794 0.004467168 NA -0.142707282 -0.125196307 1.501066e-197
# 17                           Beta - Omicron (BA.5) -0.147169896 0.006317242 NA -0.159551461 -0.134788330 4.805678e-120
# 18                                   Alpha - Other -0.132599575 0.010587455 NA -0.153350605 -0.111848546  5.503196e-36
# 19                                   Alpha - Delta -0.064060653 0.010980899 NA -0.085582819 -0.042538487  5.417062e-09
# 20                          Alpha - Omicron (BA.2) -0.070470571 0.011057209 NA -0.092142302 -0.048798840  1.850404e-10
# 21                          Alpha - Omicron (BA.3) -0.058699829 0.012591888 NA -0.083379476 -0.034020181  3.135811e-06
# 22          Alpha - Omicron (BA.1.1, with S:R346K)  0.003729877 0.013387204 NA -0.022508560  0.029968314  7.805402e-01
# 23                          Alpha - Omicron (BA.4) -0.183035683 0.011339460 NA -0.205260617 -0.160810749  1.303752e-58
# 24                          Alpha - Omicron (BA.5) -0.196253784 0.012187254 NA -0.220140364 -0.172367205  2.422484e-58
# 25                                   Other - Delta  0.068538922 0.002172024 NA  0.064281833  0.072796011 1.515540e-218
# 26                          Other - Omicron (BA.2)  0.062129005 0.003290361 NA  0.055680016  0.068577993  1.600132e-79
# 27                          Other - Omicron (BA.3)  0.073899747 0.006846207 NA  0.060481428  0.087318065  3.663963e-27
# 28          Other - Omicron (BA.1.1, with S:R346K)  0.136329452 0.008241602 NA  0.120176208  0.152482696  1.840081e-61
# 29                          Other - Omicron (BA.4) -0.050436108 0.004141814 NA -0.058553914 -0.042318302  4.106540e-34
# 30                          Other - Omicron (BA.5) -0.063654209 0.006089871 NA -0.075590137 -0.051718280  1.427551e-25
# 31                          Delta - Omicron (BA.2) -0.006409918 0.004255180 NA -0.014749918  0.001930082  1.319697e-01
# 32                          Delta - Omicron (BA.3)  0.005360824 0.007495270 NA -0.009329635  0.020051284  4.744684e-01
# 33          Delta - Omicron (BA.1.1, with S:R346K)  0.067790530 0.008705429 NA  0.050728202  0.084852858  6.853462e-15
# 34                          Delta - Omicron (BA.4) -0.118975030 0.004949426 NA -0.128675726 -0.109274334 1.110528e-127
# 35                          Delta - Omicron (BA.5) -0.132193131 0.006654578 NA -0.145235864 -0.119150399  8.176584e-88
# 36                 Omicron (BA.2) - Omicron (BA.3)  0.011770742 0.006915483 NA -0.001783355  0.025324839  8.873937e-02
# 37 Omicron (BA.2) - Omicron (BA.1.1, with S:R346K)  0.074200448 0.008290588 NA  0.057951195  0.090449701  3.556031e-19
# 38                 Omicron (BA.2) - Omicron (BA.4) -0.112565112 0.003459613 NA -0.119345830 -0.105784395 3.207494e-232
# 39                 Omicron (BA.2) - Omicron (BA.5) -0.125783214 0.005665369 NA -0.136887133 -0.114679294 3.276192e-109
# 40 Omicron (BA.3) - Omicron (BA.1.1, with S:R346K)  0.062429706 0.010068490 NA  0.042695828  0.082163584  5.628291e-10
# 41                 Omicron (BA.3) - Omicron (BA.4) -0.124335854 0.006911289 NA -0.137881733 -0.110789976  2.322895e-72
# 42                 Omicron (BA.3) - Omicron (BA.5) -0.137553955 0.008401631 NA -0.154020849 -0.121087062  3.016343e-60
# 43 Omicron (BA.1.1, with S:R346K) - Omicron (BA.4) -0.186765560 0.008606990 NA -0.203634950 -0.169896170 2.083569e-104
# 44 Omicron (BA.1.1, with S:R346K) - Omicron (BA.5) -0.199983661 0.009756332 NA -0.219105720 -0.180861602  2.250956e-93
# 45                 Omicron (BA.4) - Omicron (BA.5) -0.013218101 0.005661848 NA -0.024315120 -0.002121082  1.956479e-02

write.csv(delta_r_southafrica_province_avg_pairw,
          paste0(".\\plots\\",plotdir,
                 "\\pairwise growth rate advantages per day South Africa_multinomial 2 df spline fit with province.csv"), 
          row.names=F)

# using simpler plain multinomial model with province included
emtrsouthafrica_province_avg_pairw1 = emtrends(fit1_southafrica_province_multi, pairwise ~ variant,  
                                              var="date_num",  mode="latent",
                                              at=list(date_num=seq(today_num, today_num, by=1)),
                                              adjust="none", df=NA)
delta_r_southafrica_province_avg_pairw1 = data.frame(confint(emtrsouthafrica_province_avg_pairw1, 
                                                            adjust="none", df=NA)$contrasts, 
                                                    p.value=as.data.frame(emtrsouthafrica_province_avg_pairw1$contrasts,
                                                                          adjust="none", df=NA)$p.value)
delta_r_southafrica_province_avg_pairw1
# contrast     estimate           SE df   asymp.LCL     asymp.UCL       p.value
# 1                            Omicron (BA.1) - Beta  0.093695904 0.0012306339 NA  0.09128391  0.0961079026  0.000000e+00
# 2                           Omicron (BA.1) - Alpha  0.085016570 0.0014193727 NA  0.08223465  0.0877984894  0.000000e+00
# 3                           Omicron (BA.1) - Other  0.099543818 0.0012392504 NA  0.09711493  0.1019727037  0.000000e+00
# 4                           Omicron (BA.1) - Delta  0.067140962 0.0011571745 NA  0.06487294  0.0694089826  0.000000e+00
# 5                  Omicron (BA.1) - Omicron (BA.2) -0.065731367 0.0011874824 NA -0.06805879 -0.0634039447  0.000000e+00
# 6                  Omicron (BA.1) - Omicron (BA.3) -0.007037444 0.0037491273 NA -0.01438560  0.0003107103  6.050595e-02
# 7  Omicron (BA.1) - Omicron (BA.1.1, with S:R346K) -0.013970995 0.0015301062 NA -0.01696995 -0.0109720416  6.803554e-20
# 8                  Omicron (BA.1) - Omicron (BA.4) -0.157957905 0.0037986635 NA -0.16540315 -0.1505126614  0.000000e+00
# 9                  Omicron (BA.1) - Omicron (BA.5) -0.172879651 0.0058833154 NA -0.18441074 -0.1613485644 8.606279e-190
# 10                                    Beta - Alpha -0.008679334 0.0007503599 NA -0.01015001 -0.0072086560  6.063781e-31
# 11                                    Beta - Other  0.005847913 0.0001662647 NA  0.00552204  0.0061737859 5.300519e-271
# 12                                    Beta - Delta -0.026554942 0.0004047190 NA -0.02734818 -0.0257617074  0.000000e+00
# 13                           Beta - Omicron (BA.2) -0.159427272 0.0017163560 NA -0.16279127 -0.1560632761  0.000000e+00
# 14                           Beta - Omicron (BA.3) -0.100733349 0.0039201114 NA -0.10841663 -0.0930500716 1.277333e-145
# 15           Beta - Omicron (BA.1.1, with S:R346K) -0.107666899 0.0019393346 NA -0.11146793 -0.1038658731  0.000000e+00
# 16                           Beta - Omicron (BA.4) -0.251653810 0.0039957940 NA -0.25948542 -0.2438221971  0.000000e+00
# 17                           Beta - Omicron (BA.5) -0.266575555 0.0060124755 NA -0.27835979 -0.2547913197  0.000000e+00
# 18                                   Alpha - Other  0.014527248 0.0007606786 NA  0.01303635  0.0160181502  2.636146e-81
# 19                                   Alpha - Delta -0.017875608 0.0008145230 NA -0.01947204 -0.0162791720 9.434479e-107
# 20                          Alpha - Omicron (BA.2) -0.150747937 0.0018563494 NA -0.15438632 -0.1471095596  0.000000e+00
# 21                          Alpha - Omicron (BA.3) -0.092054014 0.0039833873 NA -0.09986131 -0.0842467186 3.717613e-118
# 22          Alpha - Omicron (BA.1.1, with S:R346K) -0.098987565 0.0020642665 NA -0.10303345 -0.0949416767  0.000000e+00
# 23                          Alpha - Omicron (BA.4) -0.242974475 0.0040578968 NA -0.25092781 -0.2350211435  0.000000e+00
# 24                          Alpha - Omicron (BA.5) -0.257896221 0.0060539258 NA -0.26976170 -0.2460307442  0.000000e+00
# 25                                   Other - Delta -0.032402855 0.0004303038 NA -0.03324624 -0.0315594754  0.000000e+00
# 26                          Other - Omicron (BA.2) -0.165275185 0.0017225444 NA -0.16865131 -0.1618990600  0.000000e+00
# 27                          Other - Omicron (BA.3) -0.106581262 0.0039228254 NA -0.11426986 -0.0988926653 1.489284e-162
# 28          Other - Omicron (BA.1.1, with S:R346K) -0.113514812 0.0019448136 NA -0.11732658 -0.1097030476  0.000000e+00
# 29                          Other - Omicron (BA.4) -0.257501723 0.0039984562 NA -0.26533855 -0.2496648926  0.000000e+00
# 30                          Other - Omicron (BA.5) -0.272423468 0.0060142451 NA -0.28421117 -0.2606357647  0.000000e+00
# 31                          Delta - Omicron (BA.2) -0.132872330 0.0016644754 NA -0.13613464 -0.1296100180  0.000000e+00
# 32                          Delta - Omicron (BA.3) -0.074178407 0.0038976650 NA -0.08181769 -0.0665391235  9.353895e-81
# 33          Delta - Omicron (BA.1.1, with S:R346K) -0.081111957 0.0018935796 NA -0.08482330 -0.0774006091  0.000000e+00
# 34                          Delta - Omicron (BA.4) -0.225098867 0.0039737856 NA -0.23288734 -0.2173103908  0.000000e+00
# 35                          Delta - Omicron (BA.5) -0.240020613 0.0059978716 NA -0.25177623 -0.2282650007  0.000000e+00
# 36                 Omicron (BA.2) - Omicron (BA.3)  0.058693923 0.0038403302 NA  0.05116701  0.0662208321  9.841204e-53
# 37 Omicron (BA.2) - Omicron (BA.1.1, with S:R346K)  0.051760373 0.0016952877 NA  0.04843767  0.0550830757 9.830663e-205
# 38                 Omicron (BA.2) - Omicron (BA.4) -0.092226538 0.0035837584 NA -0.09925058 -0.0852025001 4.796326e-146
# 39                 Omicron (BA.2) - Omicron (BA.5) -0.107148283 0.0057470973 NA -0.11841239 -0.0958841795  1.415471e-77
# 40 Omicron (BA.3) - Omicron (BA.1.1, with S:R346K) -0.006933550 0.0039731296 NA -0.01472074  0.0008536406  8.096563e-02
# 41                 Omicron (BA.3) - Omicron (BA.4) -0.150920461 0.0052686994 NA -0.16124692 -0.1405939996 1.864753e-180
# 42                 Omicron (BA.3) - Omicron (BA.5) -0.165842207 0.0069240583 NA -0.17941311 -0.1522713015 8.894317e-127
# 43 Omicron (BA.1.1, with S:R346K) - Omicron (BA.4) -0.143986910 0.0039872896 NA -0.15180185 -0.1361719665 1.498166e-285
# 44 Omicron (BA.1.1, with S:R346K) - Omicron (BA.5) -0.158908656 0.0060063009 NA -0.17068079 -0.1471365228 3.031773e-154
# 45                 Omicron (BA.4) - Omicron (BA.5) -0.014921746 0.0058358389 NA -0.02635978 -0.0034837117  1.056049e-02

write.csv(delta_r_southafrica_province_avg_pairw1,
          paste0(".\\plots\\",plotdir,
                 "\\pairwise growth rate advantages per day South Africa_multinomial fit with province.csv"), 
          row.names=F)


# fitted prop of different variantS today (on average across provinces)
# based on model with province included as a factor
multinom_preds_today_avg_province = data.frame(emmeans(fit2_southafrica_province_multi, ~ variant|1,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg_province
# variant         prob           SE df     asymp.LCL    asymp.UCL
# 1 Omicron (BA.1) 6.919546e-04 1.654332e-04 NA  3.677115e-04 1.016198e-03
# 2            Beta 1.092857e-09 5.294747e-10 NA  5.510549e-11 2.130608e-09
# 3           Alpha 8.913647e-15 2.083082e-14 NA -3.191401e-14 4.974131e-14
# 4           Other 1.255755e-01 2.202183e-02 NA  8.241350e-02 1.687375e-01
# 5           Delta 6.131224e-05 2.659282e-05 NA  9.191271e-06 1.134332e-04
# 6  Omicron (BA.2) 1.957366e-01 2.535679e-02 NA  1.460382e-01 2.454350e-01
# 7  Omicron (BA.3) 5.987460e-03 2.382789e-03 NA  1.317281e-03 1.065764e-02
# 8  Omicron (BA.4) 5.684798e-01 2.698617e-02 NA  5.155879e-01 6.213717e-01
# 9  Omicron (BA.5) 1.034673e-01 1.330131e-02 NA  7.739723e-02 1.295374e-01

# fitted prop of different variantS today (by province)
multinom_preds_today_byprovince = data.frame(emmeans(fit2_southafrica_province_multi, ~ variant|province,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_byprovince
#            variant      province         prob           SE df     asymp.LCL    asymp.UCL
# 1  Omicron (BA.1)       Gauteng 3.204233e-05 1.178136e-05 NA  8.951286e-06 5.513338e-05
# 2             Beta       Gauteng 1.335817e-11 7.883256e-12 NA -2.092730e-12 2.880907e-11
# 3            Alpha       Gauteng 3.773083e-16 8.893357e-16 NA -1.365758e-15 2.120374e-15
# 4            Other       Gauteng 1.716202e-03 6.778596e-04 NA  3.876215e-04 3.044782e-03
# 5            Delta       Gauteng 7.858067e-07 4.477449e-07 NA -9.175722e-08 1.663371e-06
# 6   Omicron (BA.2)       Gauteng 1.876008e-02 5.228968e-03 NA  8.511487e-03 2.900866e-02
# 7   Omicron (BA.3)       Gauteng 4.861098e-04 2.186500e-04 NA  5.756369e-05 9.146560e-04
# 8   Omicron (BA.4)       Gauteng 9.353804e-01 2.538635e-02 NA  8.856241e-01 9.851368e-01
# 9   Omicron (BA.5)       Gauteng 4.362435e-02 2.309144e-02 NA -1.634032e-03 8.888274e-02
# 10 Omicron (BA.1)    North West 2.779876e-04 1.508299e-04 NA -1.763353e-05 5.736087e-04
# 11            Beta    North West 1.716621e-10 1.244049e-10 NA -7.216700e-11 4.154911e-10
# 12           Alpha    North West 5.988831e-15 1.422022e-14 NA -2.188229e-14 3.385995e-14
# 13           Other    North West 4.193481e-02 2.407251e-02 NA -5.246456e-03 8.911607e-02
# 14           Delta    North West 1.175528e-05 8.380362e-06 NA -4.669929e-06 2.818049e-05
# 15  Omicron (BA.2)    North West 7.285574e-02 3.540441e-02 NA  3.464370e-03 1.422471e-01
# 16  Omicron (BA.3)    North West 9.735562e-03 5.774160e-03 NA -1.581583e-03 2.105271e-02
# 17  Omicron (BA.4)    North West 8.751842e-01 6.119574e-02 NA  7.552427e-01 9.951256e-01
# 18  Omicron (BA.5)    North West 1.037954e-46 7.289320e-47 NA -3.907269e-47 2.466634e-46
# 19 Omicron (BA.1)    Mpumalanga 2.048457e-04 1.259937e-04 NA -4.209751e-05 4.517889e-04
# 20            Beta    Mpumalanga 1.627059e-10 1.291656e-10 NA -9.045396e-11 4.158658e-10
# 21           Alpha    Mpumalanga 4.599441e-16 1.147427e-15 NA -1.788972e-15 2.708860e-15
# 22           Other    Mpumalanga 2.660489e-02 1.744808e-02 NA -7.592717e-03 6.080249e-02
# 23           Delta    Mpumalanga 5.787177e-06 4.521892e-06 NA -3.075568e-06 1.464992e-05
# 24  Omicron (BA.2)    Mpumalanga 1.261739e-01 7.142156e-02 NA -1.380984e-02 2.661576e-01
# 25  Omicron (BA.3)    Mpumalanga 1.383161e-02 8.973313e-03 NA -3.755765e-03 3.141898e-02
# 26  Omicron (BA.4)    Mpumalanga 8.331790e-01 9.444254e-02 NA  6.480750e-01 1.018283e+00
# 27  Omicron (BA.5)    Mpumalanga 4.631669e-46 3.489895e-46 NA -2.208399e-46 1.147174e-45
# 28 Omicron (BA.1)       Limpopo 3.529993e-05 2.183750e-05 NA -7.500785e-06 7.810064e-05
# 29            Beta       Limpopo 1.088192e-12 1.351152e-12 NA -1.560017e-12 3.736402e-12
# 30           Alpha       Limpopo 1.327125e-17 3.447377e-17 NA -5.429609e-17 8.083859e-17
# 31           Other       Limpopo 2.153021e-04 2.500589e-04 NA -2.748043e-04 7.054084e-04
# 32           Delta       Limpopo 9.116605e-08 1.123836e-07 NA -1.291019e-07 3.114339e-07
# 33  Omicron (BA.2)       Limpopo 4.412133e-02 2.493669e-02 NA -4.753682e-03 9.299633e-02
# 34  Omicron (BA.3)       Limpopo 9.815639e-04 6.832600e-04 NA -3.576011e-04 2.320729e-03
# 35  Omicron (BA.4)       Limpopo 9.546456e-01 2.562425e-02 NA  9.044230e-01 1.004868e+00
# 36  Omicron (BA.5)       Limpopo 8.180782e-07 6.129298e-07 NA -3.832421e-07 2.019399e-06
# 37 Omicron (BA.1)  Western Cape 1.464573e-03 4.169096e-04 NA  6.474448e-04 2.281700e-03
# 38            Beta  Western Cape 1.462987e-09 7.543685e-10 NA -1.554820e-11 2.941522e-09
# 39           Alpha  Western Cape 3.380991e-14 7.899052e-14 NA -1.210087e-13 1.886285e-13
# 40           Other  Western Cape 9.786776e-02 2.629608e-02 NA  4.632839e-02 1.494071e-01
# 41           Delta  Western Cape 1.202245e-04 5.752743e-05 NA  7.472776e-06 2.329761e-04
# 42  Omicron (BA.2)  Western Cape 2.912871e-01 5.145697e-02 NA  1.904333e-01 3.921409e-01
# 43  Omicron (BA.3)  Western Cape 6.990773e-03 2.763884e-03 NA  1.573661e-03 1.240789e-02
# 44  Omicron (BA.4)  Western Cape 5.112148e-01 8.190815e-02 NA  3.506778e-01 6.717518e-01
# 45  Omicron (BA.5)  Western Cape 9.105477e-02 5.071593e-02 NA -8.346639e-03 1.904562e-01
# 46 Omicron (BA.1)    Free State 2.074239e-03 5.604180e-04 NA  9.758397e-04 3.172638e-03
# 47            Beta    Free State 5.107055e-09 2.648733e-09 NA -8.436533e-11 1.029848e-08
# 48           Alpha    Free State 2.792513e-14 6.621399e-14 NA -1.018519e-13 1.577022e-13
# 49           Other    Free State 3.501967e-01 9.818812e-02 NA  1.577515e-01 5.426419e-01
# 50           Delta    Free State 1.155960e-04 5.772583e-05 NA  2.455443e-06 2.287365e-04
# 51  Omicron (BA.2)    Free State 6.417628e-01 9.789571e-02 NA  4.498907e-01 8.336349e-01
# 52  Omicron (BA.3)    Free State 5.809549e-03 6.064555e-03 NA -6.076760e-03 1.769586e-02
# 53  Omicron (BA.4)    Free State 5.031511e-08 1.743531e-08 NA  1.614254e-08 8.448769e-08
# 54  Omicron (BA.5)    Free State 4.106209e-05 2.273707e-05 NA -3.501751e-06 8.562593e-05
# 55 Omicron (BA.1) KwaZulu Natal 1.609288e-04 6.710749e-05 NA  2.940057e-05 2.924571e-04
# 56            Beta KwaZulu Natal 1.944358e-10 1.178150e-10 NA -3.647740e-11 4.253491e-10
# 57           Alpha KwaZulu Natal 1.080217e-15 2.580378e-15 NA -3.977230e-15 6.137665e-15
# 58           Other KwaZulu Natal 2.818121e-02 1.155933e-02 NA  5.525344e-03 5.083708e-02
# 59           Delta KwaZulu Natal 2.467697e-05 1.398245e-05 NA -2.728122e-06 5.208206e-05
# 60  Omicron (BA.2) KwaZulu Natal 5.878744e-02 1.991200e-02 NA  1.976064e-02 9.781424e-02
# 61  Omicron (BA.3) KwaZulu Natal 1.075692e-03 6.265381e-04 NA -1.523000e-04 2.303684e-03
# 62  Omicron (BA.4) KwaZulu Natal 1.152851e-01 5.474961e-02 NA  7.977814e-03 2.225924e-01
# 63  Omicron (BA.5) KwaZulu Natal 7.964850e-01 7.143315e-02 NA  6.564786e-01 9.364914e-01
# 64 Omicron (BA.1)  Eastern Cape 4.063217e-04 1.972883e-04 NA  1.964375e-05 7.929996e-04
# 65            Beta  Eastern Cape 3.350881e-10 2.284475e-10 NA -1.126608e-10 7.828371e-10
# 66           Alpha  Eastern Cape 2.102455e-15 5.077285e-15 NA -7.848841e-15 1.205375e-14
# 67           Other  Eastern Cape 3.809006e-02 1.962218e-02 NA -3.687057e-04 7.654883e-02
# 68           Delta  Eastern Cape 4.594079e-05 2.961414e-05 NA -1.210185e-05 1.039834e-04
# 69  Omicron (BA.2)  Eastern Cape 6.819123e-02 2.858991e-02 NA  1.215603e-02 1.242264e-01
# 70  Omicron (BA.3)  Eastern Cape 1.837276e-03 1.300420e-03 NA -7.115000e-04 4.386053e-03
# 71  Omicron (BA.4)  Eastern Cape 8.914292e-01 4.610173e-02 NA  8.010714e-01 9.817869e-01
# 72  Omicron (BA.5)  Eastern Cape 3.309349e-48 2.183432e-48 NA -9.700978e-49 7.588797e-48
# 73 Omicron (BA.1) Northern Cape 1.571354e-03 4.100621e-04 NA  7.676468e-04 2.375061e-03
# 74            Beta Northern Cape 2.387330e-09 1.177781e-09 NA  7.892151e-11 4.695739e-09
# 75           Alpha Northern Cape 8.465762e-15 2.044815e-14 NA -3.161188e-14 4.854340e-14
# 76           Other Northern Cape 5.453726e-01 8.800016e-02 NA  3.728954e-01 7.178497e-01
# 77           Delta Northern Cape 2.269525e-04 9.643145e-05 NA  3.795034e-05 4.159547e-04
# 78  Omicron (BA.2) Northern Cape 4.396901e-01 8.713852e-02 NA  2.689018e-01 6.104785e-01
# 79  Omicron (BA.3) Northern Cape 1.313901e-02 8.684101e-03 NA -3.881515e-03 3.015954e-02
# 80  Omicron (BA.4) Northern Cape 5.715074e-13 2.013098e-13 NA  1.769474e-13 9.660674e-13
# 81  Omicron (BA.5) Northern Cape 4.408291e-47 2.453973e-47 NA -4.014083e-48 9.217990e-47





# PLOT MULTINOMIAL FIT

extrapolate = 30
date.from = as.numeric(min(GISAID_sel$date_num))
date.to = max(GISAID_sel$date_num)+extrapolate

# multinomial model predictions by province (fastest, but no confidence intervals)
predgrid = expand.grid(list(date_num=seq(date.from, date.to),
                            province=levels_PROVINCES))
fit_southafrica_multi_preds_byprovince = data.frame(predgrid, as.data.frame(predict(fit2_southafrica_province_multi, 
                                                                                    newdata=predgrid, type="prob")),check.names=F)
fit_southafrica_multi_preds_byprovince = gather(fit_southafrica_multi_preds_byprovince, variant, prob, all_of(levels_VARIANTS), factor_key=TRUE)
fit_southafrica_multi_preds_byprovince$date = as.Date(fit_southafrica_multi_preds_byprovince$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_byprovince$variant = factor(fit_southafrica_multi_preds_byprovince$variant, levels=levels_VARIANTS)
fit_southafrica_multi_preds_byprovince$province = factor(fit_southafrica_multi_preds_byprovince$province, levels=levels_PROVINCES)

# multinomial model predictions overall for South Africa (with model without province) (fastest, but no confidence intervals)
predgrid = expand.grid(list(date_num=seq(date.from, date.to)))
fit_southafrica_multi_preds = data.frame(predgrid, as.data.frame(predict(fit2_southafrica_multi, 
                                                                         newdata=predgrid, type="prob")),check.names=F)
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



# Muller plot overall for South Africa
fit_southafrica_multi_preds$variant = factor(fit_southafrica_multi_preds$variant, levels=levels_VARIANTS_plot)
muller_southafrica_mfit = ggplot(data=fit_southafrica_multi_preds, 
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



# MULTINOMIAL MODEL PREDICTIONS with confidence intervals for all of South Africa combined (slower) ####
# using fit with province included fit2_southafrica_province_multi
# I WILL USE THIS MODEL
fit_southafrica_multi_preds_withCI = data.frame(emmeans(fit2_southafrica_province_multi, # fit2_southafrica_multi,
                                                        ~ variant,
                                                        by=c("date_num"),
                                                        at=list(date_num=seq(date.from, date.to, by=7)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_southafrica_multi_preds_withCI$date = as.Date(fit_southafrica_multi_preds_withCI$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_withCI$variant = factor(fit_southafrica_multi_preds_withCI$variant, 
                                                    levels=levels_VARIANTS_plot)


# CALCULATE PAIRWISE GROWTH RATE ADVANTAGES OF EACH VARIANT AT TIMEPOINT WHERE IT FIRST EXCEEDED 1% OF THE CASES ####
# BASED ON MULTINOMIAL MULTINOMIAL 2 DF SPLINE FIT WITH PROVINCE INCLUDED AS MAIN EFFECT
# (PS should still be subsetted to LINEAGES THAT THEN CO-OCCURRED AT REASONABLE FREQUENCY)
library(dplyr)
# fit_southafrica_multi_preds_maxima = as.data.frame(fit_southafrica_multi_preds_withCI[fit_southafrica_multi_preds_withCI$date<=today,]
#                                                 %>% group_by(variant) %>% top_n(1, prob))
# fit_southafrica_multi_preds_maxima = fit_southafrica_multi_preds_maxima[fit_southafrica_multi_preds_maxima$variant!="Other",]
date_at1perc = as.data.frame(fit_southafrica_multi_preds_withCI[fit_southafrica_multi_preds_withCI$date<=today,]
                             %>% group_by(variant) %>% summarise(DATE_1PERCENT=first(date[prob>=0.01])))
date_at1perc = date_at1perc[date_at1perc$variant!="Other"&!is.na(date_at1perc$DATE_1PERCENT),]

growth_advantages_at1perc = do.call(rbind, lapply(date_at1perc$variant,
                                                  function (LINEAGE) { 
                                                    t = as.numeric(date_at1perc[date_at1perc$variant==LINEAGE,"DATE_1PERCENT"]) 
                                                    emtrSA_pw = emtrends(fit2_southafrica_province_multi, 
                                                                         pairwise ~ variant,  
                                                                              var="date_num",  mode="latent",
                                                                              at=list(date_num=t))
                                                    delta_r_SA_pw = data.frame(LINEAGE_AT1PERCENT = LINEAGE,
                                                                                    DATE_NUM_AT1PERCENT=t, 
                                                                                    DATE_AT1PERCENT=as.Date(t, origin="1970-01-01"),
                                                                                    confint(emtrSA_pw, 
                                                                                            adjust="none", df=NA)$contrasts, 
                                                                                    p.value=as.data.frame(emtrSA_pw$contrasts)$p.value) 
                                                    return(delta_r_SA_pw)
                                                  }) )

write.csv(growth_advantages_at1perc, ".//data//GISAID//south africa//pairwise daily growth rate advantage when lineages at 1 percent_south africa_multinomial 2 df spline fit with province.csv",
          row.names=F)


# PLOT OF MULTINOMIAL FIT ON LOGIT SCALE ####
fit_SA_multi_pr = fit_southafrica_multi_preds_withCI
ymin = 0.001
ymax = 0.99
fit_SA_multi_pr$asymp.LCL[fit_SA_multi_pr$asymp.LCL<ymin] = ymin
fit_SA_multi_pr$asymp.UCL[fit_SA_multi_pr$asymp.UCL<ymin] = ymin
fit_SA_multi_pr$asymp.UCL[fit_SA_multi_pr$asymp.UCL>ymax] = ymax
fit_SA_multi_pr$prob[fit_SA_multi_pr$prob<ymin] = ymin
# data_agbyweek2$prop[data_agbyweek2$collection_date<=as.Date("2021-09-01")&grepl("Omicron", data_agbyweek2$LINEAGE)] = 0 # fix above

plot_southafrica_mfit_logit = qplot(data=fit_SA_multi_pr, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN SOUTH AFRICA\n(GISAID data, multinomial spline fit VARIANT ~ ns(DATE, df=2)+PROVINCE)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today+extrapolate, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today+extrapolate, by="month")),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
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
                        range=c(0.01, 5), limits=c(1,10^4), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-03-01"),NA), ylim=c(0.001, 0.9901), expand=c(0,0))
plot_southafrica_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom 2 df spline fit with province_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\southafrica_multinom fit_logit scale.pdf"), width=8, height=6)









# multinomial model predictions with confidence intervals by province (slower)
fit_southafrica_multi_preds_byprovince_withCI = data.frame(emmeans(fit2_southafrica_province_multi,
                                                        ~ variant,
                                                        by=c("date_num","province"),
                                                        at=list(date_num=seq(date.from, date.to, by=7)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_southafrica_multi_preds_byprovince_withCI$date = as.Date(fit_southafrica_multi_preds_byprovince_withCI$date_num, origin="1970-01-01")
fit_southafrica_multi_preds_byprovince_withCI$variant = factor(fit_southafrica_multi_preds_byprovince_withCI$variant,
                                                               levels=levels_VARIANTS_plot)




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




# 2. PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

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
tail(hosp_prov)
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
data_prov_wide$daily_cases[data_prov_wide$daily_cases<0] = 0
# fill in missing dates with NAs
data_prov_wide2 = expand.grid(date=seq(min(data_prov_wide$date), max(data_prov_wide$date), by=1), province=levels(data_prov_wide$province)) 
hosp_prov_wide = merge(data_prov_wide2, data_prov_wide, all.x=TRUE)

hosp_prov_wide = as.data.frame(hosp_prov_wide %>% # test to see if diff(AdmissionstoDate) has less reporting lag than AdmissionsinPreviousDay
                      group_by(province) %>%
                      mutate(NewAdmissions = AdmissionstoDate - lag(AdmissionstoDate)))
data_prov_wide = hosp_prov_wide


# AdmissionsinPreviousDay has underreporting compared to diff(AdmissionstoDate)
# but the latter has a lot of reporting backlogs & reporting artefacts
qplot(data=data_prov_wide, x=date, y=log10(abs(NewAdmissions)+1), geom="col", fill=I("blue"), width=I(1.2)) + 
  facet_wrap(~ province, scale="free_y")  + 
  geom_col(aes(x=date, y=-log10(abs(AdmissionsinPreviousDay)+1)), fill=I("red"), width=I(1.2)) +
  xlab("reporting date") +
  ylab("diff(AdmissionstoDate) (blue) and AdmissionsinPreviousDay (red)") + 
  ggtitle("diff(AdmissionstoDate) (blue) and AdmissionsinPreviousDay (red), data NICD") + xaxis + # + scale_y_log10()
  coord_cartesian(xlim=c(as.Date("2020-03-01"), NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\hosps AdmissionsinPreviousDay diff NewAdmissions by province.png"), width=10, height=6)


# plots of cases, hospitalisations & deaths by province
qplot(data=data_prov_wide, x=date, y=daily_cases, geom="col", colour=I("blue")) + facet_wrap(~ province, scale="free_y")  + 
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

qplot(data=data_prov_wide, x=date, y=log10(daily_cases+1), geom="col", fill=I("blue"), width=I(1.2)) + 
  facet_wrap(~ province, scale="free_y")  + 
  geom_col(aes(x=date-6, y=-log10((AdmissionsinPreviousDay+1)*7)), fill=I("red"), width=I(1.2)) +
  xlab("case reporting date") +
  ylab("log10(cases per day) (top) &\nlog10(hospitalisations per week 7 days later) (bottom)") + 
  ggtitle("Confirmed Covid cases & hospitalisations in South Africa",
          "Log10 of new confirmed Covid cases per day (blue)\n& hospitalised patients testing positive per week one week later (red), data NICD") + xaxis + # + scale_y_log10()
  coord_cartesian(xlim=c(as.Date("2020-11-01"), NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day and hosps per week mirrored by province_log10 Y scale.png"), width=10, height=6)


qplot(data=data_prov_wide, x=date, y=AdmissionsinPreviousDay, geom="col", width=I(1.3), colour=I("red")) + facet_wrap(~ province, scale="free_y") + # + # + scale_y_log10()
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



qplot(data=hosp_prov[hosp_prov$variable=="CurrentlyAdmitted",], x=date, y=value, geom="col", width=I(1.3), colour=I("blue")) + facet_wrap(~ province, scale="free_y") +
  ylab("currently admitted") + ggtitle("Covid patients currently admitted in hospital\nin South Africa") + xaxis # + # + scale_y_log10()
# facet_wrap(~ variable)
ggsave(file=paste0(".\\plots\\",plotdir,"\\currently admitted by province.png"), width=10, height=6)

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

# daily cases & estimated infections per day overall for South Africa
cases_tot = cases_cum[,c("date","cases_daily", "tests_daily")] 
# TO DO: replace by sum of per province data / data_prov_wide (appears more up to date)
colnames(cases_tot) = c("date", "cases_new", "tests_new")
cases_tot = rbind(data.frame(date=seq.Date(from=min(cases_tot$date)-30,
                                           to=min(cases_tot$date)-1, by=1), 
                             cases_new=0,
                             tests_new=0), cases_tot)
cases_tot$tests_new = abs(cases_tot$tests_new)
cases_tot$tests_new[is.na(cases_tot$tests_new)] = cases_tot$cases_new[is.na(cases_tot$tests_new)]
cases_tot$date_num = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)

cases_tot$testspercase = (cases_tot$tests_new)/(cases_tot$cases_new+1)
cases_tot$posratio = (cases_tot$cases_new)/(cases_tot$tests_new+1)
qplot(data=cases_tot, x=date, y=cases_new, geom="col")
qplot(data=cases_tot, x=date, y=tests_new, geom="col")
qplot(data=cases_tot, x=date, y=testspercase, geom="col")
qplot(data=cases_tot, x=date, y=posratio, geom="col") + ylab("positivity ratio") + ylim(c(0,1)) # positivity ratio
tail(cases_tot$posratio) # at 17% now

# also add smoothed case, growth in cases, smoothed tests, growth in tests, smoothed positivity ratio & growth in positivity ratio
library(mgcv)
k=50
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

date.from = as.numeric(as.Date("2020-03-14")) # as.numeric(min(GISAID_sel$date_num))
date.to = today_num+7 # max(GISAID_sel$date_num)+extrapolate

growth_cases = as.data.frame(emtrends(fit_cases, ~ date_num, var="date_num", 
                                     at=list(date_num = seq(date.from,
                                                            date.to, by=1)),
                                     type="link"))
colnames(growth_cases)[2] = "growth_cases"
colnames(growth_cases)[5] = "growth_cases_LOWER"
colnames(growth_cases)[6] = "growth_cases_UPPER"
growth_cases$DATE = as.Date(growth_cases$date_num, origin="1970-01-01")

cases_smoothed = as.data.frame(emmeans(fit_cases, ~ date_num, 
                                        at=list(date_num = seq(date.from,
                                                               date.to, by=1)),
                                        type="response"))
colnames(cases_smoothed)[2] = "cases_smoothed"
colnames(cases_smoothed)[5] = "cases_smootheds_LOWER"
colnames(cases_smoothed)[6] = "cases_smoothed_UPPER"
cases_smoothed$DATE = as.Date(cases_smoothed$date_num, origin="1970-01-01")


fit_tests = gam(tests_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) +
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
BIC(fit_tests) # 9953.436

growth_tests = as.data.frame(emtrends(fit_tests, ~ date_num, var="date_num", 
                                      at=list(date_num = seq(date.from,
                                                             date.to, by=1)),
                                      type="link"))
colnames(growth_tests)[2] = "growth_tests"
colnames(growth_tests)[5] = "growth_tests_LOWER"
colnames(growth_tests)[6] = "growth_testss_UPPER"
growth_tests$DATE = as.Date(growth_tests$date_num, origin="1970-01-01")

tests_smoothed = as.data.frame(emmeans(fit_tests, ~ date_num, 
                                       at=list(date_num = seq(date.from,
                                                              date.to, by=1)),
                                       type="response"))
colnames(tests_smoothed)[2] = "tests_smoothed"
colnames(tests_smoothed)[5] = "tests_smoothed_LOWER"
colnames(tests_smoothed)[6] = "tests_smoothed_UPPER"
tests_smoothed$DATE = as.Date(tests_smoothed$date_num, origin="1970-01-01")

fit_posratio = gam(cbind(cases_new, tests_new-cases_new) ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) +
                  WEEKDAY, # + # + offset(log(testspercase)),
                # BANKHOLIDAY,
                # + s(tests_new, bs="cs", k=8, fx=F),
                family=binomial, data=cases_tot,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+1*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-1*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 
BIC(fit_posratio) 

growth_posratio = as.data.frame(emtrends(fit_posratio, ~ date_num, var="date_num", 
                                      at=list(date_num = seq(date.from,
                                                             date.to, by=1)),
                                      type="link"))
colnames(growth_posratio)[2] = "growth_posratio"
colnames(growth_posratio)[5] = "growth_posratio_LOWER"
colnames(growth_posratio)[6] = "growth_posratio_UPPER"
growth_posratio$DATE = as.Date(growth_posratio$date_num, origin="1970-01-01")

posratio_smoothed = as.data.frame(emmeans(fit_posratio, ~ date_num, 
                                       at=list(date_num = seq(date.from,
                                                              date.to, by=1)),
                                       type="response"))
colnames(posratio_smoothed)[2] = "posratio_smoothed"
colnames(posratio_smoothed)[5] = "posratio_smoothed_LOWER"
colnames(posratio_smoothed)[6] = "posratio_smoothed_UPPER"
posratio_smoothed$DATE = as.Date(posratio_smoothed$date_num, origin="1970-01-01")

cases_tot$growth_cases = growth_cases$growth_cases[match(cases_tot$date, growth_cases$DATE)]
cases_tot$cases_smoothed = cases_smoothed$cases_smoothed[match(cases_tot$date, cases_smoothed$DATE)]
cases_tot$growth_tests = growth_tests$growth_tests[match(cases_tot$date, growth_tests$DATE)]
cases_tot$tests_smoothed = tests_smoothed$tests_smoothed[match(cases_tot$date, tests_smoothed$DATE)]
cases_tot$growth_posratio = growth_posratio$growth_posratio[match(cases_tot$date, growth_posratio$DATE)]
cases_tot$posratio_smoothed = posratio_smoothed$posratio_smoothed[match(cases_tot$date, posratio_smoothed$DATE)]

# estimated true infections (IHME)
# (https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00484-6/fulltext?utm_campaign=lancetcovid22&utm_source=twitter&utm_medium=social#sec1, https://github.com/ihmeuw/covid-historical-model, https://github.com/ihmeuw/covid-model-infections, downloaded from https://ourworldindata.org/grapher/daily-new-estimated-covid-19-infections-ihme-model?country=~BEL)
infections_IHME = read.csv(".//data//daily-new-estimated-covid-19-infections-ihme-model.csv")
colnames(infections_IHME) = c("COUNTRY","ISO3","DATE","infections_IHME","infections_IHME_lower","infections_IHME_upper","infections_IHME_rolling_7d_mean")
infections_IHME$DATE = as.Date(infections_IHME$DATE)
infections_IHME$DATE_NUM = as.numeric(infections_IHME$DATE)
infections_IHME_SA = infections_IHME[infections_IHME$ISO3=="ZAF",]
cases_tot$infections_IHME = infections_IHME_SA$infections_IHME[match(cases_tot$date, infections_IHME$DATE)]

# estimate new infections from cases & pos ratio
library(Hmisc)
library(splines)
fit_infections = glm(Lag(infections_IHME, shift=8) ~ ns(date_num, df=3) * log(cases_smoothed), # + 
                       # log(tests_smoothed),
                       # ns(log(tests_smoothed), df=2), 
                     family=poisson, data=cases_tot[!is.na(cases_tot$infections_IHME),]) # good model

# fit_infections = glm(Lag(infections_IHME, shift=8) ~ ns(date_num, df=3) + 
#                        log(cases_smoothed) + 
#                        ns(growth_posratio, df=2),  
#                      family=poisson, data=cases_tot[!is.na(cases_tot$infections_IHME),])


plot(model.frame(fit_infections)[,1], col="black", type="l")
lines(fitted(fit_infections), col="red")
sum(fitted(fit_infections)) # 62650373, SA population is 60,652,262 as of Thursday, April 21, 2022, based on Worldometer elaboration of the latest United Nations data
sum(cases_tot$infections_IHME, na.rm=T) # 62843552

plot(predict(fit_infections, newdata=cases_tot, type="response"), type="l", col="red")
sum(predict(fit_infections, newdata=cases_tot, type="response"), na.rm=T) # 78 343 441

cases_tot$est_infections = predict(fit_infections, newdata=cases_tot, type="response")
cases_tot$est_infections[is.na(cases_tot$est_infections)] = 0
cases_tot$est_infections_cumsum = cumsum(cases_tot$est_infections)
populationSA = 60652262
cases_tot$est_infections_cumsum_perpop = 100*cases_tot$est_infections_cumsum/populationSA
plot(cases_tot$est_infections_cumsum_perpop, type="l", ylab="Cumulative infections (IHME & NICD case data) (% of population)")

# also add variable with estimated excess deaths
# TO DO : plot infection fatality rate over time based on IHME true infection estimates & excess mortality ests of the Economist
excess_mort = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/export_country_cumulative.csv")
excess_mort$date = as.Date(excess_mort$date)
excess_mort_SA = excess_mort[excess_mort$iso3c=="ZAF",]
# TO DO: interpolate, as data now only given per week
cases_tot$cumulative_estimated_daily_excess_deaths_7d_later = excess_mort_SA$cumulative_estimated_daily_excess_deaths[match(cases_tot$date+7,excess_mort_SA$date)]      
cases_tot$infection_fatality_rate = 100*cases_tot$cumulative_estimated_daily_excess_deaths_7d_later/cases_tot$est_infections_cumsum
plot(x=cases_tot$date, y=cases_tot$infection_fatality_rate, type="p", ylim=c(0,1))
sort(unique(excess_mort$iso3c))
plot(x=excess_mort[excess_mort$iso3c=="BEL",]$date,
     y=excess_mort[excess_mort$iso3c=="BEL",]$cumulative_estimated_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="BEL",]$date,
     y=excess_mort[excess_mort$iso3c=="BEL",]$cumulative_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="BEL",]$date,
     y=excess_mort[excess_mort$iso3c=="BEL",]$cumulative_daily_covid_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="ZAF",]$date,
     y=excess_mort[excess_mort$iso3c=="ZAF",]$cumulative_estimated_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="ZAF",]$date,
     y=excess_mort[excess_mort$iso3c=="ZAF",]$cumulative_daily_excess_deaths,
     type="l")
plot(x=excess_mort[excess_mort$iso3c=="ZAF",]$date,
     y=excess_mort[excess_mort$iso3c=="ZAF",]$cumulative_daily_covid_deaths,
     type="l")


# CALCULATE Re VALUES THROUGH TIME

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time, Nishiura et al. 2020: mean 4.7 d & SD 2.9 d)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

# smooth out weekday effects in case nrs using negative binomial GAM (possibility to also correct for variable testing intensity)
library(mgcv)
k=50
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

# model with correction for variable testing intensity
fit_cases_testing = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) +
                  WEEKDAY + # + # + offset(log(testspercase)),
                # BANKHOLIDAY,
                s(tests_new, bs="cs", k=8, fx=F),
                family=nb(), data=cases_tot,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+0.75*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-0.75*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 

# calculate average instantaneous growth rates & 95% CLs & Re values using emtrends ####
# based on the slope of the GAM fit on a log link scale
date.from = as.numeric(as.Date("2020-03-14")) # as.numeric(min(GISAID_sel$date_num))
date.to = today_num+7 # max(GISAID_sel$date_num)+extrapolate

avg_r_cases = as.data.frame(emtrends(fit_cases_testing, ~ date_num, var="date_num", 
                                     at=list(date_num = seq(date.from,
                                                          date.to, by=1),
                                             tests_new = max(cases_tot$tests_new,na.rm=T),
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
library(ggthemes)
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
k=50
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

qplot(data=fit_southafrica_multi_preds, x=date, y=cases, geom="line", colour=variant) +
  xlim(c(as.Date("2021-10-01"), today)) + scale_y_log10() + 
  geom_line(aes(y=smoothed_cases), lwd=I(2)) +
  ylab("new Omicron cases per day") +
  scale_colour_manual(values=lineage_cols_plot) +
  ggtitle("Inferred nr. of new Omicron cases per day in SA", "(black=observed,\nred=GAM negative binomial spline smooth with correction for\nvariable testing intensity & weekday effects & marginal means\ncalculated at uniform, maximal testing effort)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\variant cases observed and fitted.png"), width=8, height=6)

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



# STACKED AREA CHART OF NEW INFECTIONS BY VARIANT (MULTINOMIAL FIT MAPPED ONTO ESTIMATED INFECTIONS) ####

fit_southafrica_multi_preds$totinfections = cases_tot$est_infections[match(round(fit_southafrica_multi_preds$date_num),cases_tot$date_num)]
fit_southafrica_multi_preds$totinfections[is.na(fit_southafrica_multi_preds$totinfections)] = 0
fit_southafrica_multi_preds$infections = fit_southafrica_multi_preds$totinfections * fit_southafrica_multi_preds$prob
fit_southafrica_multi_preds$infections[fit_southafrica_multi_preds$infections<=0.001] = 0

fit_southafrica_multi_preds$variant = factor(fit_southafrica_multi_preds$variant, levels=levels_VARIANTS_plot)

ggplot(data=fit_southafrica_multi_preds, 
       aes(x=date, y=infections*100/populationSA, group=variant)) + 
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant, width=I(1.1)), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("Estimated infections per day (% of population)") + xlab("Date of diagnosis") +
  ggtitle("ESTIMATED SARS-CoV2 INFECTIONS PER DAY BY VARIANT IN SOUTH AFRICA","(case data NICD mapped to true number of daily infections as estimated by IHME using model\nglm(Lag(infections_IHME, shift=8) ~ ns(date, df=3)*log(cases_smoothed), family=poisson(log)),\ncombined with multinomial 2 df spline fit to GISAID lineage frequency data)") +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
# coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\estimated infections IHME per day by variant_stacked area multinomial fit.png"), width=8, height=6)

# plot of cumulative nr of infections by variant
library(dplyr)
fit_southafrica_multi_preds = fit_southafrica_multi_preds %>% 
  group_by(variant) %>% 
  arrange(date) %>% 
  mutate(cuminfections = cumsum(ifelse(is.na(infections), 0, infections)) + infections*0,
         cuminfectionsperpop = 100*(cumsum(ifelse(is.na(infections), 0, infections)) + infections*0)/populationSA)

ggplot(data=fit_southafrica_multi_preds[fit_southafrica_multi_preds$date<=today,], 
       aes(x=date, y=cuminfectionsperpop, group=variant)) + 
  geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant, width=I(1.1)), 
           position = position_stack(reverse = T)) +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("Cumulative infections (% of population)") + xlab("Date of diagnosis") +
  ggtitle("CUMULATIVE SARS-CoV2 INFECTIONS BY VARIANT IN SOUTH AFRICA","(case data NICD mapped to true number of daily infections as estimated by IHME using model\nglm(Lag(infections_IHME, shift=8) ~ ns(date, df=3)*log(cases_smoothed), family=poisson(log)),\ncombined with multinomial 2 df spline fit to GISAID lineage frequency data)") +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylim(c(0,150))
# coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cumulative estimated infections IHME by variant_stacked area multinomial fit raw case data.png"), width=8, height=6)
write.csv(fit_southafrica_multi_preds, file=paste0(".\\plots\\",plotdir,"\\cases infections per day by variant South Africa 5 Jan 2022.csv"), row.names=F)



# STACKED AREA PLOTS OF NEW CASES BY VARIANT PER PROVINCE

fit_southafrica_multi_preds_byprovince$totcases = data_prov_wide$daily_cases[match(interaction(fit_southafrica_multi_preds_byprovince$province, round(fit_southafrica_multi_preds_byprovince$date_num)),
                                                                                 interaction(data_prov_wide$province, as.numeric(data_prov_wide$date)))]
fit_southafrica_multi_preds_byprovince$cases = fit_southafrica_multi_preds_byprovince$totcases * fit_southafrica_multi_preds_byprovince$prob
fit_southafrica_multi_preds_byprovince$cases[fit_southafrica_multi_preds_byprovince$cases<=0.001] = NA
k=50
data_prov_wide$date_num = as.numeric(data_prov_wide$date)
data_prov_wide$WEEKDAY = factor(weekdays(data_prov_wide$date))
fit_cases_byprovince = gam(daily_cases ~ s(date_num, bs="cs", k=k, m=c(2), fx=F, by=province) + province +
                             WEEKDAY, 
                        family=nb(), data=data_prov_wide,
                        method = "REML",
                        knots = list(date_num = c(min(data_prov_wide$date_num)-14,
                                                  seq(min(data_prov_wide$date_num)+1*diff(range(data_prov_wide$date_num))/(k-2), 
                                                      max(data_prov_wide$date_num)-1*diff(range(data_prov_wide$date_num))/(k-2), length.out=k-2),
                                                  max(data_prov_wide$date_num)+14))
) 
cases_emmeans_byprov = as.data.frame(emmeans(fit_cases_byprovince, 
                                             ~ 1, by=c("date_num","province"), 
                                             at=list(date_num=seq(date.from, date.to, by=1)),
                                             rg.limit=100000, 
                                             type="response"))
fit_southafrica_multi_preds_byprovince$smoothed_totcases = cases_emmeans_byprov$response[match(interaction(fit_southafrica_multi_preds_byprovince$date_num,fit_southafrica_multi_preds_byprovince$province),
                                                                                               interaction(cases_emmeans_byprov$date_num, cases_emmeans_byprov$province))]
fit_southafrica_multi_preds_byprovince$smoothed_cases = fit_southafrica_multi_preds_byprovince$smoothed_totcases * fit_southafrica_multi_preds_byprovince$prob
fit_southafrica_multi_preds_byprovince$smoothed_cases[fit_southafrica_multi_preds_byprovince$smoothed_cases<=0.001] = NA
fit_southafrica_multi_preds_byprovince$variant = factor(fit_southafrica_multi_preds_byprovince$variant, levels=levels_VARIANTS_plot)

# fit_southafrica_multi_preds_byprovince[fit_southafrica_multi_preds$date==as.Date("2021-12-03"),] # this date is not plotted in plot below FIX

qplot(data=fit_southafrica_multi_preds_byprovince[fit_southafrica_multi_preds_byprovince$smoothed_cases>=10,], 
      x=date, y=smoothed_cases, geom="line", lwd=I(1.2), colour=variant) +
  facet_wrap(~province) +
  xlim(c(as.Date("2020-03-01"), today)) + scale_y_log10() + 
  ylab("new Omicron cases per day") +
  scale_colour_manual(values=lineage_cols_plot) +
  ggtitle("Inferred nr. of new Omicron cases per day in SA", "(black=observed,\nred=GAM negative binomial spline smooth with correction for\nvariable testing intensity & weekday effects & marginal means\ncalculated at uniform, maximal testing effort)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\variant cases log10 Y scale by province.png"), width=8, height=6)

ggplot(data=fit_southafrica_multi_preds_byprovince, 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~province, scale="free_y") +
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
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data_by province.png"), width=8, height=6)

ggplot(data=fit_southafrica_multi_preds_byprovince[fit_southafrica_multi_preds_byprovince$date<=today,],  # or (today-1) or (today-2) to be on the safe side
       aes(x=date, y=smoothed_cases, group=variant)) + 
  facet_wrap(~province, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA","(negative binomial fit to NICD case data, with correction for weekday effects &\nbut not for variable testing intensity;\nvariant frequencies based on multinomial spline fit to GISAID data)") + #  
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit case data_by province.png"), width=8, height=6)



# TO DO : plot growth rate in cases by variant





# DIDN'T TRY TO RUN / UPDATE THE PART BELOW YET

# EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for simplest model fit1_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from+16,
                                                  date.to-extrapolate), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit1_southafrica_multi, ~ variant , 
                                                                           at=list(date_num=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_VARIANTS)
                                                cons = lapply(seq_along(wt), function (i) { 
                                                  con = -wt; con[i] = 1 + con[i]; con })
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
      
