# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN ENGLAND BASED ON SANGER INSTITUTE BASELINE SURVEILLANCE SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & most surge testing/active surveillance)

# last update 25 JULY 2022

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) 
# today = as.Date("2021-08-26")
# today # "2022-01-03"
today_num = as.numeric(today)
plotdir = "UK_SANGER"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))



sanger = read_tsv("https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv")
sanger = as.data.frame(sanger)

# code ONS regions
LTLAs_regions = read.csv(".//data/ONS/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")
sanger$REGION = LTLAs_regions$RGN20NM[match(sanger$LTLA, LTLAs_regions$LAD20CD)]
sanger = sanger[!is.na(sanger$REGION),]
sanger$REGION = factor(sanger$REGION)
levels(sanger$REGION)
levels_REGION = c("London","North West","South West","South East","East of England","East Midlands","West Midlands",
                  "North East","Yorkshire and The Humber")
sanger$REGION = factor(sanger$REGION, levels=levels_REGION)
head(sanger)
# code NHS England regions
LTLAs_NHSregions = read.csv(".//data/ONS/LTLA_to_NHSregion_lookup.csv")
sanger$NHSREGION = LTLAs_NHSregions$NHSREGION[match(sanger$LTLA, LTLAs_NHSregions$LTLA)]
sanger = sanger[!is.na(sanger$REGION),]
sanger$NHSREGION = factor(sanger$NHSREGION)
levels(sanger$NHSREGION)
levels_NHSREGION = c("London","North West","South West","South East","East of England","Midlands",
                     "North East and Yorkshire")
sanger$NHSREGION = factor(sanger$NHSREGION, levels=levels_NHSREGION)
head(sanger)

sanger = sanger[!(sanger$Lineage=="None"|sanger$Lineage=="Lineage data suppressed"),]
range(sanger$WeekEndDate) # "2020-09-05" "2022-07-16"
# sanger$Week = lubridate::week(sanger$WeekEndDate)
sanger$DATE_NUM = as.numeric(sanger$WeekEndDate) # using start of week
colnames(sanger)

library(dplyr)
sanger$LINEAGE = case_when(
  sanger$Lineage=="BA.2.75" ~ "Omicron (BA.2.75)",
  sanger$Lineage=="BA.2.74" ~ "Omicron (BA.2.74)",
  sanger$Lineage=="BA.2.76" ~ "Omicron (BA.2.76)",
  (sanger$Lineage=="B.1.617.2")|grepl("AY", sanger$Lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", sanger$Lineage) ~ "Alpha",
  grepl("B.1.351", sanger$Lineage, fixed=T) ~ "Beta",
  grepl("^BA\\.1$|BA\\.1\\.", sanger$Lineage) ~ "Omicron (BA.1)",
  grepl("^BA\\.3$|BA\\.3\\.", sanger$Lineage) ~ "Omicron (BA.3)",
  grepl("^BA\\.4",sanger$Lineage) ~ "Omicron (BA.4)",
  grepl("^BA\\.5",sanger$Lineage)|grepl("BE|BF",sanger$Lineage) ~ "Omicron (BA.5)",
  grepl("^BA\\.2",sanger$Lineage) ~ "Omicron (BA.2)",
  sanger$Lineage!="Unassigned" ~ "Other" # remaining Unassigned will be assigned as NA
)

sanger = sanger[rep(seq_len(nrow(sanger)), sanger$Count),] # convert to long format
sanger$Count = NULL
nrow(sanger) # 1542126

table(sanger$LINEAGE)

sel_target_VOC = "Omicron (BA.2.75)"
sel_reference_VOC = "Omicron (BA.5)"
levels_LINEAGE = c(sel_reference_VOC, "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", "Omicron (BA.2.74)", "Omicron (BA.2.76)", sel_target_VOC)
levels_LINEAGE_plot = c("Other", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", "Omicron (BA.5)", "Omicron (BA.2.74)", "Omicron (BA.2.76)", "Omicron (BA.2.75)")

n = length(levels_LINEAGE)
lineage_cols = hcl(h = seq(0, 180, length = n), l = 60, c = 180)
lineage_cols[which(levels_LINEAGE=="Alpha")] = "#0085FF"
lineage_cols[which(levels_LINEAGE=="Beta")] = "green4"
lineage_cols[which(levels_LINEAGE=="Delta")] = "mediumorchid"
# lineage_cols[which(levels_LINEAGE=="C.1.2")] = "darkorange"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.1)")] = "red" # "magenta"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.2)")] = "red3"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.3)")] = "red4" 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.4)")] = "darkorange" 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.5)")] = "darkorange3" 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.2.74)")] = "black"
lineage_cols[which(levels_LINEAGE=="Omicron (BA.2.76)")] = muted(muted("magenta")) 
lineage_cols[which(levels_LINEAGE=="Omicron (BA.2.75)")] = "magenta" 
lineage_cols[which(levels_LINEAGE=="Other")] = "grey65"

lineage_cols_plot = lineage_cols[match(levels_LINEAGE_plot,levels_LINEAGE)]


# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of England
data_agbyweek1 = as.data.frame(table(sanger$WeekEndDate, sanger$LINEAGE))
colnames(data_agbyweek1) = c("WeekEndDate", "LINEAGE", "count")
data_agbyweek1_sum = aggregate(count ~ WeekEndDate, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$WeekEndDate, data_agbyweek1_sum$WeekEndDate)]
data_agbyweek1$WeekEndDate = as.Date(as.character(data_agbyweek1$WeekEndDate))
data_agbyweek1$collection_date = data_agbyweek1$WeekEndDate # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweek1$LINEAGE = factor(data_agbyweek1$LINEAGE, levels=levels_LINEAGE)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and ONS region
data_agbyweekregion1 = as.data.frame(table(sanger$WeekEndDate, sanger$REGION, sanger$LINEAGE))
colnames(data_agbyweekregion1) = c("WeekEndDate", "REGION", "LINEAGE", "count")
data_agbyweekregion1_sum = aggregate(count ~ WeekEndDate + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$WeekEndDate,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$WeekEndDate,data_agbyweekregion1_sum$REGION))]
data_agbyweekregion1$WeekEndDate = as.Date(as.character(data_agbyweekregion1$WeekEndDate))
data_agbyweekregion1$collection_date = data_agbyweekregion1$WeekEndDate # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweekregion1$LINEAGE = factor(data_agbyweekregion1$LINEAGE, levels=levels_LINEAGE)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion1[data_agbyweekregion1$collection_date==max(data_agbyweekregion1$collection_date),]

# aggregated by week and NHS region
data_agbyweeknhsregion1 = as.data.frame(table(sanger$WeekEndDate, sanger$NHSREGION, sanger$LINEAGE))
colnames(data_agbyweeknhsregion1) = c("WeekEndDate", "NHSREGION", "LINEAGE", "count")
data_agbyweeknhsregion1_sum = aggregate(count ~ WeekEndDate + NHSREGION, data=data_agbyweeknhsregion1, sum)
data_agbyweeknhsregion1$total = data_agbyweeknhsregion1_sum$count[match(interaction(data_agbyweeknhsregion1$WeekEndDate,data_agbyweeknhsregion1$NHSREGION), 
                                                                  interaction(data_agbyweeknhsregion1_sum$WeekEndDate,data_agbyweeknhsregion1_sum$NHSREGION))]
data_agbyweeknhsregion1$WeekEndDate = as.Date(as.character(data_agbyweeknhsregion1$WeekEndDate))
data_agbyweeknhsregion1$collection_date = data_agbyweeknhsregion1$WeekEndDate-3.5 # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweeknhsregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweeknhsregion1$LINEAGE = factor(data_agbyweeknhsregion1$LINEAGE, levels=levels_LINEAGE)
data_agbyweeknhsregion1$NHSREGION = factor(data_agbyweeknhsregion1$NHSREGION, levels=levels_NHSREGION)
data_agbyweeknhsregion1$collection_date_num = as.numeric(data_agbyweeknhsregion1$collection_date)
data_agbyweeknhsregion1$prop = data_agbyweeknhsregion1$count/data_agbyweeknhsregion1$total
data_agbyweeknhsregion1 = data_agbyweeknhsregion1[data_agbyweeknhsregion1$total!=0,]

data_agbyweeknhsregion1[data_agbyweeknhsregion1$collection_date==max(data_agbyweeknhsregion1$collection_date),]



# MULLER PLOT (RAW DATA)
unique(sanger$LINEAGE)

library(scales)

data_agbyweek1$LINEAGE = factor(data_agbyweek1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweekregion1$LINEAGE = factor(data_agbyweekregion1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweeknhsregion1$LINEAGE = factor(data_agbyweeknhsregion1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweeknhsregion1$NHSREGION = factor(data_agbyweeknhsregion1$NHSREGION, levels=levels_NHSREGION)


library(ggplot2)
library(ggthemes)
muller_sanger_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sanger_raw1


muller_sangerbyregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)") +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbyregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_raw data.pdf"), width=8, height=6)


muller_sangerbynhsregion_raw1 = ggplot(data=data_agbyweeknhsregion1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  facet_wrap(~ NHSREGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN ENGLAND") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)") +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_sangerbynhsregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
sanger$LINEAGE = relevel(sanger$LINEAGE, ref="Omicron (BA.2)")
# sanger$DATE_NUM_WEEK = sanger$DATE_NUM/7

data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
data_agbyweekregion1$LINEAGE = relevel(data_agbyweekregion1$LINEAGE, ref="Omicron (BA.2)")
data_agbyweeknhsregion1$DATE_NUM = data_agbyweeknhsregion1$collection_date_num
data_agbyweeknhsregion1$LINEAGE = relevel(data_agbyweeknhsregion1$LINEAGE, ref="Omicron (BA.2)")

# by ONS region
set.seed(1)
fit1_sanger_multi = nnet::multinom(LINEAGE ~ REGION + DATE_NUM, data=data_agbyweekregion1, weights=count, maxit=1000)
fit2_sanger_multi = nnet::multinom(LINEAGE ~ REGION * DATE_NUM, data=data_agbyweekregion1, weights=count, maxit=1000)
fit3_sanger_multi = nnet::multinom(LINEAGE ~ REGION + ns(DATE_NUM, df=2), data=data_agbyweekregion1, weights=count, maxit=1000)
fit4_sanger_multi = nnet::multinom(LINEAGE ~ REGION * ns(DATE_NUM, df=2), data=data_agbyweekregion1, weights=count, maxit=1000)
fit5_sanger_multi = nnet::multinom(LINEAGE ~ REGION + ns(DATE_NUM, df=3), data=data_agbyweekregion1, weights=count, maxit=1000)
fit6_sanger_multi = nnet::multinom(LINEAGE ~ REGION * ns(DATE_NUM, df=3), data=data_agbyweekregion1, weights=count, maxit=1000)
BIC(fit1_sanger_multi, fit2_sanger_multi, fit3_sanger_multi, fit4_sanger_multi, fit5_sanger_multi, fit6_sanger_multi) 
# fit6_sanger_multi fits best (lowest BIC)

# by NHS region
set.seed(1)
fit1_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION + DATE_NUM, data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit2_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION * DATE_NUM, data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit3_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION + ns(DATE_NUM, df=2), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit4_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION * ns(DATE_NUM, df=2), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit5_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION + ns(DATE_NUM, df=3), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
fit6_sanger_nhs_multi = nnet::multinom(LINEAGE ~ NHSREGION * ns(DATE_NUM, df=3), data=data_agbyweeknhsregion1, weights=count, maxit=1000)
BIC(fit1_sanger_nhs_multi, fit2_sanger_nhs_multi, fit3_sanger_nhs_multi, fit4_sanger_nhs_multi, fit5_sanger_nhs_multi, fit6_sanger_nhs_multi) 
# fit6_sanger_multi fits best (lowest BIC)

# growth rate advantage of different lineages over BA.2 based on model fit6_sanger_multi
emtrsanger6 = emtrends(fit6_sanger_multi, trt.vs.ctrl1 ~ LINEAGE,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger6 = data.frame(confint(emtrsanger6, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger6$contrasts)$p.value)
delta_r_sanger6
# contrast    estimate          SE df   asymp.LCL   asymp.UCL      p.value
# 1             Other - Omicron (BA.2)  0.10016595 0.001027839 NA  0.09815142  0.10218048 0.000000e+00
# 2             Alpha - Omicron (BA.2) -0.13060802 0.009927530 NA -0.15006562 -0.11115042 0.000000e+00
# 3              Beta - Omicron (BA.2) -0.58193609 0.022138774 NA -0.62532729 -0.53854489 0.000000e+00
# 4             Gamma - Omicron (BA.2) -0.41417375 0.045411525 NA -0.50317871 -0.32516880 0.000000e+00
# 5             Delta - Omicron (BA.2)  0.05297104 0.001266519 NA  0.05048871  0.05545338 0.000000e+00
# 6    Omicron (BA.1) - Omicron (BA.2) -0.11620039 0.001191832 NA -0.11853634 -0.11386445 0.000000e+00
# 7    Omicron (BA.3) - Omicron (BA.2)  0.02024411 0.002439222 NA  0.01546332  0.02502490 2.109424e-15
# 8  Omicron (BA.1.1) - Omicron (BA.2) -0.08844363 0.001050212 NA -0.09050201 -0.08638525 0.000000e+00
# 9    Omicron (BA.4) - Omicron (BA.2)  0.42743305 0.027397165 NA  0.37373559  0.48113051 0.000000e+00
# 10   Omicron (BA.5) - Omicron (BA.2)  0.36762081 0.035401923 NA  0.29823431  0.43700730 0.000000e+00

# pairwise growth rate advantages for all strain comparisons (i.e. pairwise differences in growth rate per day among the different lineages)
emtrsanger_pairw6 = emtrends(fit6_sanger_multi, pairwise ~ LINEAGE,  
                            var="DATE_NUM",  mode="latent",
                            at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_pairw6 = data.frame(confint(emtrsanger_pairw6, 
                                          adjust="none", df=NA)$contrasts, 
                                  p.value=as.data.frame(emtrsanger_pairw6$contrasts)$p.value)
delta_r_sanger_pairw6
# contrast    estimate          SE df   asymp.LCL   asymp.UCL      p.value
# 1         Delta - other  0.13143581 0.021564483 NA  0.08917019  0.17370142 1.038571e-07
# 2    Delta - (B.1.177+)  1.65647626 0.025240149 NA  1.60700647  1.70594604 0.000000e+00
# 3         Delta - Alpha  0.10597787 0.014520829 NA  0.07751757  0.13443817 1.152013e-10
# 4          Delta - Beta  0.80431957 0.025145670 NA  0.75503496  0.85360418 0.000000e+00
# 5         Delta - Gamma  0.69438136 0.040310542 NA  0.61537415  0.77338857 0.000000e+00
# 6       Delta - Omicron -0.23634806 0.002666635 NA -0.24157456 -0.23112155 0.000000e+00
# 7    other - (B.1.177+)  1.52504045 0.032257174 NA  1.46181755  1.58826335 0.000000e+00
# 8         other - Alpha -0.02545793 0.020340701 NA -0.06532497  0.01440911 8.728342e-01
# 9          other - Beta  0.67288377 0.031640158 NA  0.61087020  0.73489733 0.000000e+00
# 10        other - Gamma  0.56294555 0.045179793 NA  0.47439479  0.65149632 0.000000e+00
# 11      other - Omicron -0.36778386 0.021728718 NA -0.41037137 -0.32519636 0.000000e+00
# 12   (B.1.177+) - Alpha -1.55049838 0.027154349 NA -1.60371993 -1.49727684 0.000000e+00
# 13    (B.1.177+) - Beta -0.85215669 0.035585801 NA -0.92190357 -0.78240980 0.000000e+00
# 14   (B.1.177+) - Gamma -0.96209490 0.054137124 NA -1.06820171 -0.85598808 0.000000e+00
# 15 (B.1.177+) - Omicron -1.89282431 0.025380582 NA -1.94256934 -1.84307929 0.000000e+00
# 16         Alpha - Beta  0.69834170 0.028145169 NA  0.64317818  0.75350522 0.000000e+00
# 17        Alpha - Gamma  0.58840349 0.042146719 NA  0.50579744  0.67100954 0.000000e+00
# 18      Alpha - Omicron -0.34232593 0.014763631 NA -0.37126211 -0.31338974 0.000000e+00
# 19         Beta - Gamma -0.10993821 0.050338496 NA -0.20859985 -0.01127657 3.085473e-01
# 20       Beta - Omicron -1.04066763 0.025286642 NA -1.09022853 -0.99110672 0.000000e+00
# 21      Gamma - Omicron -0.93072942 0.040398636 NA -1.00990929 -0.85154954 0.000000e+00



# transm advantage over BA.2  with fixed generation time of 4.7 days (Nishiura et al. 2020)
transmadv_sanger6 =  sign(delta_r_sanger6[,c(2,5,6)])*100*(exp(abs(delta_r_sanger6[,c(2,5,6)])*4.7)-1)
transmadv_sanger6 =  data.frame(contrast=delta_r_sanger6$contrast, transmadv_sanger6)
transmadv_sanger6
# contrast      estimate     asymp.LCL     asymp.UCL
# 1      other - Delta     -85.47470    -126.23332     -52.05923
# 2 (B.1.177+) - Delta -240432.04138 -303393.83354 -190532.08717
# 3      Alpha - Delta     -64.55851     -88.11051     -43.95529
# 4       Beta - Delta   -4282.92754   -5425.39391   -3376.68494
# 5      Gamma - Delta   -2514.32856   -3689.90590   -1703.39934
# 6    Omicron - Delta     203.68958     196.32044     211.24197

emtrsanger_byregion6 = emtrends(fit6_sanger_multi, trt.vs.ctrl1 ~ LINEAGE|REGION,  
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(sanger$DATE_NUM)))
delta_r_sanger_byregion6 = data.frame(confint(emtrsanger_byregion6, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrsanger_byregion6$contrasts)$p.value)
delta_r_sanger_byregion6
#              contrast                   REGION      estimate          SE df    asymp.LCL     asymp.UCL      p.value
# 1       other - Delta                   London -0.2030419074 0.029561692 NA -0.260981759 -0.1451020555 4.085211e-10
# 2  (B.1.177+) - Delta                   London -1.6165483490 0.024751074 NA -1.665059562 -1.5680371361 0.000000e+00
# 3       Alpha - Delta                   London -0.0117891309 0.010973291 NA -0.033296385  0.0097181234 7.321369e-01
# 4        Beta - Delta                   London  0.0094558538 0.017639544 NA -0.025117018  0.0440287252 9.577026e-01
# 5       Gamma - Delta                   London -0.0024804661 0.023099385 NA -0.047754429  0.0427934970 9.996297e-01
# 6     Omicron - Delta                   London  0.2329669016 0.003280587 NA  0.226537070  0.2393967334 0.000000e+00
# 7       other - Delta               North West -0.0003647597 0.008440503 NA -0.016907841  0.0161783218 9.999758e-01
# 8  (B.1.177+) - Delta               North West  0.0496731780 0.045513507 NA -0.039531656  0.1388780120 7.220419e-01
# 9       Alpha - Delta               North West -0.0287183883 0.014829389 NA -0.057783457  0.0003466802 2.271950e-01
# 10       Beta - Delta               North West -1.1603813115 0.069541044 NA -1.296679253 -1.0240833701 0.000000e+00
# 11      Gamma - Delta               North West -1.9622415488 0.033298833 NA -2.027506062 -1.8969770356 0.000000e+00
# 12    Omicron - Delta               North West  0.2673371037 0.004946273 NA  0.257642588  0.2770316199 0.000000e+00
# 13      other - Delta               South West -0.0812938367 0.076716774 NA -0.231655950  0.0690682769 7.407321e-01
# 14 (B.1.177+) - Delta               South West -1.9910109107 0.034429014 NA -2.058490539 -1.9235312824 0.000000e+00
# 15      Alpha - Delta               South West -0.0927497376 0.045784860 NA -0.182486415 -0.0030130603 1.908832e-01
# 16       Beta - Delta               South West -0.4140407623 0.082035196 NA -0.574826792 -0.2532547322 5.661955e-06
# 17      Gamma - Delta               South West -1.0680725583 0.074372838 NA -1.213840642 -0.9223044743 0.000000e+00
# 18    Omicron - Delta               South West  0.2462318474 0.004950658 NA  0.236528737  0.2559349578 0.000000e+00
# 19      other - Delta               South East  0.0237536826 0.005969352 NA  0.012053968  0.0354533973 5.475680e-04
# 20 (B.1.177+) - Delta               South East -2.0306624179 0.043621441 NA -2.116158872 -1.9451659641 0.000000e+00
# 21      Alpha - Delta               South East  0.0174068211 0.011489407 NA -0.005112002  0.0399256447 4.535568e-01
# 22       Beta - Delta               South East  0.0412514279 0.017044649 NA  0.007844530  0.0746583259 7.927311e-02
# 23      Gamma - Delta               South East -0.9566399167 0.092262906 NA -1.137471890 -0.7758079434 0.000000e+00
# 24    Omicron - Delta               South East  0.2131750254 0.003688099 NA  0.205946484  0.2204035670 0.000000e+00
# 25      other - Delta          East of England  0.0014299630 0.007191517 NA -0.012665151  0.0155250771 9.976714e-01
# 26 (B.1.177+) - Delta          East of England -2.5776091937 0.058285476 NA -2.691846628 -2.4633717593 0.000000e+00
# 27      Alpha - Delta          East of England -0.0137096415 0.014121574 NA -0.041387418  0.0139681355 7.904439e-01
# 28       Beta - Delta          East of England -1.4935008074 0.072276314 NA -1.635159780 -1.3518418349 0.000000e+00
# 29      Gamma - Delta          East of England -0.2808708744 0.196701096 NA -0.666397939  0.1046561900 5.086685e-01
# 30    Omicron - Delta          East of England  0.1998047563 0.003947408 NA  0.192067979  0.2075415339 0.000000e+00
# 31      other - Delta            East Midlands -0.2402797584 0.061384299 NA -0.360590775 -0.1199687423 7.029407e-04
# 32 (B.1.177+) - Delta            East Midlands -0.5775604887 0.029260530 NA -0.634910074 -0.5202109037 0.000000e+00
# 33      Alpha - Delta            East Midlands -0.1714308751 0.031632366 NA -0.233429174 -0.1094325766 9.503416e-07
# 34       Beta - Delta            East Midlands -2.0561361501 0.034708234 NA -2.124163039 -1.9881092608 0.000000e+00
# 35      Gamma - Delta            East Midlands -1.6107157294 0.080338898 NA -1.768177076 -1.4532543828 0.000000e+00
# 36    Omicron - Delta            East Midlands  0.2270992781 0.005045595 NA  0.217210093  0.2369884629 0.000000e+00
# 37      other - Delta            West Midlands -0.4624217349 0.073879806 NA -0.607223494 -0.3176199758 1.231836e-08
# 38 (B.1.177+) - Delta            West Midlands -3.1933515549 0.027340415 NA -3.246937784 -3.1397653262 0.000000e+00
# 39      Alpha - Delta            West Midlands -0.0190703932 0.027747161 NA -0.073453830  0.0353130439 9.153804e-01
# 40       Beta - Delta            West Midlands -0.6340701665 0.040569430 NA -0.713584788 -0.5545555449 0.000000e+00
# 41      Gamma - Delta            West Midlands -1.5366163741 0.088235439 NA -1.709554657 -1.3636780913 0.000000e+00
# 42    Omicron - Delta            West Midlands  0.2426045180 0.005366932 NA  0.232085525  0.2531235110 0.000000e+00
# 43      other - Delta               North East -0.0828759407 0.127266919 NA -0.332314519  0.1665626373 9.270499e-01
# 44 (B.1.177+) - Delta               North East -0.0170235293 0.032091433 NA -0.079921582  0.0458745233 9.589418e-01
# 45      Alpha - Delta               North East -0.3829959004 0.080184831 NA -0.540155281 -0.2258365200 1.952772e-05
# 46       Beta - Delta               North East -0.8104033916 0.047629729 NA -0.903755944 -0.7170508390 0.000000e+00
# 47      Gamma - Delta               North East  0.5046758650 0.023099385 NA  0.459401902  0.5499498281 0.000000e+00
# 48    Omicron - Delta               North East  0.2403503493 0.019867960 NA  0.201409863  0.2792908358 0.000000e+00
# 49      other - Delta Yorkshire and The Humber -0.1378279584 0.059782109 NA -0.254998739 -0.0206571782 1.041826e-01
# 50 (B.1.177+) - Delta Yorkshire and The Humber -2.9541930353 0.032387201 NA -3.017670783 -2.8907152878 0.000000e+00
# 51      Alpha - Delta Yorkshire and The Humber -0.2507436069 0.061792242 NA -0.371854175 -0.1296330386 4.028764e-04
# 52       Beta - Delta Yorkshire and The Humber -0.7210508284 0.035502403 NA -0.790634260 -0.6514673968 0.000000e+00
# 53      Gamma - Delta Yorkshire and The Humber  0.6635293623 0.023099385 NA  0.618255399  0.7088033254 0.000000e+00
# 54    Omicron - Delta Yorkshire and The Humber  0.2575627200 0.006168628 NA  0.245472431  0.2696530091 0.000000e+00




# predicted incidences on average over all ONS regions from multinomial fit
# fitted prop of different LINEAGES in England today
multinom_preds_today_avg6 = data.frame(emmeans(fit6_sanger_multi, ~ LINEAGE|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg6
# LINEAGE         prob           SE df     asymp.LCL    asymp.UCL
# 1    Delta 2.717537e-02 4.987916e-03 NA  1.739923e-02 3.695150e-02
# 2    other 3.053303e-07 1.620503e-07 NA -1.228242e-08 6.229431e-07
# 3 B.1.177+ 6.457009e-08 2.488948e-07 NA -4.232547e-07 5.523949e-07
# 4    Alpha 3.427218e-08 4.118860e-08 NA -4.645599e-08 1.150004e-07
# 5     Beta 3.643311e-08 7.608304e-08 NA -1.126869e-07 1.855531e-07
# 6    Gamma 1.237640e-09 3.803893e-09 NA -6.217854e-09 8.693134e-09
# 7  Omicron 9.728242e-01 4.987935e-03 NA  9.630480e-01 9.826004e-01


# here given by region:
multinom_preds_today_byregion6 = data.frame(emmeans(fit6_sanger_multi, ~ LINEAGE|DATE_NUM, by=c("REGION"),
                                                   at=list(DATE_NUM=today_num), 
                                                   mode="prob", df=NA))
multinom_preds_today_byregion6
# LINEAGE DATE_NUM                   REGION          prob            SE df      asymp.LCL     asymp.UCL
# 1     Delta    18995                   London  2.676098e-03  2.025275e-04 NA   2.279151e-03  3.073044e-03
# 2     other    18995                   London  2.732980e-20  1.240113e-19 NA  -2.157279e-19  2.703875e-19
# 3  B.1.177+    18995                   London 1.256917e-131 4.780791e-131 NA -8.113261e-131 1.062710e-130
# 4     Alpha    18995                   London  2.110154e-08  3.199519e-08 NA  -4.160787e-08  8.381096e-08
# 5      Beta    18995                   London  2.603604e-08  6.542315e-08 NA  -1.021910e-07  1.542631e-07
# 6     Gamma    18995                   London  1.113876e-08  3.423504e-08 NA  -5.596069e-08  7.823820e-08
# 7   Omicron    18995                   London  9.973238e-01  2.025319e-04 NA   9.969269e-01  9.977208e-01
# 8     Delta    18995               North West  4.706646e-03  4.827694e-04 NA   3.760436e-03  5.652857e-03
# 9     other    18995               North West  1.816530e-07  1.965149e-07 NA  -2.035090e-07  5.668151e-07
# 10 B.1.177+    18995               North West  2.346483e-08  1.794132e-07 NA  -3.281785e-07  3.751082e-07
# 11    Alpha    18995               North West  5.554701e-09  1.111600e-08 NA  -1.623227e-08  2.734167e-08
# 12     Beta    18995               North West  1.783651e-90  2.100273e-89 NA  -3.938094e-89  4.294824e-89
# 13    Gamma    18995               North West 3.008837e-133 9.356385e-133 NA -1.532934e-132 2.134701e-132
# 14  Omicron    18995               North West  9.952931e-01  4.827910e-04 NA   9.943469e-01  9.962394e-01
# 15    Delta    18995               South West  1.234432e-02  1.181594e-03 NA   1.002844e-02  1.466020e-02
# 16    other    18995               South West  4.006909e-13  4.549369e-12 NA  -8.515908e-12  9.317289e-12
# 17 B.1.177+    18995               South West 2.970175e-159 1.360960e-158 NA -2.370415e-158 2.964450e-158
# 18    Alpha    18995               South West  3.072914e-13  1.948750e-12 NA  -3.512189e-12  4.126771e-12
# 19     Beta    18995               South West  7.264691e-37  1.018009e-35 NA  -1.922613e-35  2.067907e-35
# 20    Gamma    18995               South West  1.394700e-87  1.306400e-86 NA  -2.421028e-86  2.699968e-86
# 21  Omicron    18995               South West  9.876557e-01  1.181594e-03 NA   9.853398e-01  9.899716e-01
# 22    Delta    18995               South East  1.156430e-02  9.005571e-04 NA   9.799241e-03  1.332936e-02
# 23    other    18995               South East  1.564852e-06  1.157185e-06 NA  -7.031879e-07  3.832892e-06
# 24 B.1.177+    18995               South East 4.566457e-165  0.000000e+00 NA  4.566457e-165 4.566457e-165
# 25    Alpha    18995               South East  2.172712e-07  3.490339e-07 NA  -4.668227e-07  9.013652e-07
# 26     Beta    18995               South East  3.018619e-07  6.815664e-07 NA  -1.033984e-06  1.637708e-06
# 27    Gamma    18995               South East  1.118100e-73  1.618703e-72 NA  -3.060789e-72  3.284409e-72
# 28  Omicron    18995               South East  9.884336e-01  9.007186e-04 NA   9.866682e-01  9.901990e-01
# 29    Delta    18995          East of England  1.759752e-02  1.434858e-03 NA   1.478525e-02  2.040979e-02
# 30    other    18995          East of England  1.001447e-06  8.624976e-07 NA  -6.890170e-07  2.691912e-06
# 31 B.1.177+    18995          East of England 4.918001e-203  0.000000e+00 NA  4.918001e-203 4.918001e-203
# 32    Alpha    18995          East of England  5.993042e-08  1.146840e-07 NA  -1.648462e-07  2.847070e-07
# 33     Beta    18995          East of England 4.156135e-110 4.872183e-109 NA -9.133690e-109 9.964917e-109
# 34    Gamma    18995          East of England  1.378003e-28  4.393603e-27 NA  -8.473504e-27  8.749105e-27
# 35  Omicron    18995          East of England  9.824014e-01  1.434944e-03 NA   9.795890e-01  9.852139e-01
# 36    Delta    18995            East Midlands  1.916293e-02  1.894923e-03 NA   1.544894e-02  2.287691e-02
# 37    other    18995            East Midlands  5.959631e-22  5.684293e-21 NA  -1.054505e-20  1.173697e-20
# 38 B.1.177+    18995            East Midlands  6.247988e-51  2.849839e-50 NA  -4.960784e-50  6.210381e-50
# 39    Alpha    18995            East Midlands  9.329635e-17  4.220381e-16 NA  -7.338832e-16  9.204759e-16
# 40     Beta    18995            East Midlands 1.439183e-144 7.282246e-144 NA -1.283376e-143 1.571212e-143
# 41    Gamma    18995            East Midlands 3.259080e-123 1.083517e-122 NA -1.797746e-122 2.449562e-122
# 42  Omicron    18995            East Midlands  9.808371e-01  1.894923e-03 NA   9.771231e-01  9.845511e-01
# 43    Delta    18995            West Midlands  1.133280e-02  1.207883e-03 NA   8.965392e-03  1.370020e-02
# 44    other    18995            West Midlands  4.797827e-39  5.621459e-38 NA  -1.053807e-37  1.149764e-37
# 45 B.1.177+    18995            West Midlands 1.389470e-246  0.000000e+00 NA  1.389470e-246 1.389470e-246
# 46    Alpha    18995            West Midlands  4.591456e-09  1.811992e-08 NA  -3.092294e-08  4.010585e-08
# 47     Beta    18995            West Midlands  9.661073e-52  6.553906e-51 NA  -1.187931e-50  1.381153e-50
# 48    Gamma    18995            West Midlands 2.800818e-126 9.616116e-126 NA -1.604642e-125 2.164806e-125
# 49  Omicron    18995            West Midlands  9.886672e-01  1.207883e-03 NA   9.862998e-01  9.910346e-01
# 50    Delta    18995               North East  1.546132e-01  4.476671e-02 NA   6.687207e-02  2.423544e-01
# 51    other    18995               North East  2.017164e-11  4.151787e-10 NA  -7.935637e-10  8.339070e-10
# 52 B.1.177+    18995               North East  5.576659e-07  2.232254e-06 NA  -3.817472e-06  4.932804e-06
# 53    Alpha    18995               North East  6.829151e-29  7.867790e-28 NA  -1.473767e-27  1.610350e-27
# 54     Beta    18995               North East  3.266089e-57  2.374494e-56 NA  -4.327315e-56  4.980533e-56
# 55    Gamma    18995               North East  6.231460e-69  1.923153e-68 NA  -3.146164e-68  4.392456e-68
# 56  Omicron    18995               North East  8.453862e-01  4.476687e-02 NA   7.576448e-01  9.331277e-01
# 57    Delta    18995 Yorkshire and The Humber  1.058051e-02  1.253081e-03 NA   8.124511e-03  1.303650e-02
# 58    other    18995 Yorkshire and The Humber  3.149033e-16  2.775081e-15 NA  -5.124156e-15  5.753963e-15
# 59 B.1.177+    18995 Yorkshire and The Humber 1.325970e-228  0.000000e+00 NA  1.325970e-228 1.325970e-228
# 60    Alpha    18995 Yorkshire and The Humber  6.025821e-24  5.419757e-23 NA  -1.001995e-22  1.122511e-22
# 61     Beta    18995 Yorkshire and The Humber  1.510827e-59  8.001277e-59 NA  -1.417139e-58  1.719304e-58
# 62    Gamma    18995 Yorkshire and The Humber 5.570002e-104 1.712703e-103 NA -2.799836e-103 3.913836e-103
# 63  Omicron    18995 Yorkshire and The Humber  9.894195e-01  1.253081e-03 NA   9.869635e-01  9.918755e-01


# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME ####

gentime = 4.7 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantage of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit6_sanger_multi = nnet::multinom(LINEAGE ~ REGION * ns(DATE_NUM, df=3), data=sanger, maxit=1000) here
emtrsanger6 = emtrends(fit6_sanger_multi, trt.vs.ctrl1 ~ LINEAGE, by=c("DATE_NUM","REGION"), 
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=seq(as.numeric(as.Date("2021-12-01")),as.numeric(as.Date("2022-01-03")))))
delta_r_sanger6 = data.frame(confint(emtrsanger6, 
                                    adjust="none", df=NA)$contrasts)
delta_r_sanger6 = delta_r_sanger6[delta_r_sanger6$contrast=="Omicron - Delta",]
delta_r_sanger6

delta_r_sanger6$collection_date = as.Date(delta_r_sanger6$DATE_NUM, origin="1970-01-01")
delta_r_sanger6$REGION = factor(delta_r_sanger6$REGION, levels=levels_REGION)

plot_sanger_mfit6_growthadv = qplot(data=delta_r_sanger6, 
      x=collection_date, y=estimate, ymin=asymp.LCL, ymax=asymp.UCL, colour=REGION, fill=REGION, group=REGION,
      geom="blank") +
  facet_wrap(~ REGION) +
  xaxis +
  scale_y_continuous(limits=c(0,max(delta_r_sanger6$asymp.UCL))) +
  geom_ribbon(aes(colour=NULL), alpha=I(0.6)) + # fill=I("steelblue"), 
  geom_line(aes(fill=NULL)) + # colour=I("steelblue"), 
  ylab("Growth rate advantage of Omicron over Delta per day") +
  theme_hc() + xlab("") +
  ggtitle("GROWTH RATE ADVANTAGE OF OMICRON OVER DELTA\n(Sanger Institute baseline surveillance data, multinomial spline fit)") +
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
  theme(legend.position = "none") +
  xlab("") +
  scale_colour_hue(h=c(10, 310), c=100) +
  scale_fill_hue(h=c(10, 310), c=100)
# theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) # +
# coord_cartesian(xlim=c(as.Date("2021-05-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit6_growthadv

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit6_growth rate advantage Omicron.png"), width=8, height=6)


# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_sanger4 = delta_r_sanger4
transmadv_sanger4$estimate =  100*(exp(transmadv_sanger4$estimate*gentime)-1)
transmadv_sanger4$asymp.LCL =  100*(exp(transmadv_sanger4$asymp.LCL*gentime)-1)
transmadv_sanger4$asymp.UCL =  100*(exp(transmadv_sanger4$asymp.UCL*gentime)-1)
transmadv_sanger4$collection_date = as.Date(transmadv_sanger4$DATE_NUM, origin="1970-01-01")
transmadv_sanger4$REGION = factor(transmadv_sanger4$REGION, levels=levels_REGION)

plot_sanger_mfit4_transmadv = qplot(data=transmadv_sanger4, 
                               x=collection_date, y=estimate, ymin=asymp.LCL, ymax=asymp.UCL, colour=REGION, fill=REGION, group=REGION,
                               geom="blank") +
  facet_wrap(~ REGION) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous(limits=c(0,max(transmadv_sanger4$asymp.UCL))) +
  geom_ribbon(aes(fill=I("steelblue"), colour=NULL), alpha=I(0.3)) +
  geom_line(aes(colour=I("steelblue"), fill=NULL)) +
  ylab("Transmission advantage of B.1.617.2 over B.1.1.7 (%)") +
  theme_hc() + xlab("") +
  ggtitle("TRANSMISSION ADVANTAGE OF B.1.617.2 OVER B.1.1.7\n(Sanger Institute baseline surveillance data, multinomial fit)") +
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
  theme(legend.position = "none") +
  xlab("") # +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) # +
  # coord_cartesian(xlim=c(as.Date("2021-05-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_sanger_mfit4_transmadv

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit3_transm advantage B16172.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit3_transm advantage B16172.pdf"), width=8, height=6)




# PLOT MULTINOMIAL FIT ####

# extrapolate = 60
date.from = as.numeric(as.Date("2020-09-01"))
date.to = today_num # max(sanger$DATE_NUM)+extrapolate

# predictions by ONS region (without CIs) ####
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to), 
                            REGION=levels_REGION))
fit_multi_preds = data.frame(predgrid, as.data.frame(predict(fit6_sanger_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fit_multi_preds = gather(fit_multi_preds, LINEAGE, prob, all_of(levels_LINEAGE), factor_key=TRUE)
fit_multi_preds$collection_date = as.Date(fit_multi_preds$DATE_NUM, origin="1970-01-01")
fit_multi_preds$LINEAGE = factor(fit_multi_preds$LINEAGE, levels=levels_LINEAGE) 

# library(glm.predict)
# fit_multi_preds_withCIs = basepredict.multinom(model=fit6_sanger_multi, values=predgrid, type="simulation",
#                                            sim.count=10, conf.int=0.95, sigma=NULL, set.seed=NULL)

muller_sangerbyregion_mfit = ggplot(data=fit_multi_preds, 
                                      aes(x=collection_date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbyregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit.png"), width=8, height=6)


library(ggpubr)
ggarrange(muller_sangerbyregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbyregion_mfit+ggtitle("Multinomial fit (LINEAGE~REGION*ns(DATE,df=3))"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.png"), width=8, height=10)





# predictions by ONS region (with CIs) ####

fit_sanger_multi_predsbyregion = data.frame(emmeans(fit8_sanger_multi, 
                                                    ~ LINEAGE,
                                                    by=c("DATE_NUM", "REGION"),
                                                    at=list(DATE_NUM=seq(as.numeric(as.Date("2021-12-01")), date.to, by=1)), 
                                                    mode="prob", df=NA))
fit_sanger_multi_predsbyregion$collection_date = as.Date(fit_sanger_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion$LINEAGE = factor(fit_sanger_multi_predsbyregion$LINEAGE, levels=levels_LINEAGE_plot) 
fit_sanger_multi_predsbyregion$REGION = factor(fit_sanger_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of England today

muller_sangerbyregion_mfit = ggplot(data=fit_sanger_multi_predsbyregion, 
                                    aes(x=collection_date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2021-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2021-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2021-12-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbyregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit.pdf"), width=8, height=6)

library(ggpubr)
ggarrange(muller_sangerbyregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbyregion_mfit+ggtitle("Multinomial fit"), ncol=1)

# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.png"), width=8, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by ONS region_multinom fit multipanel.pdf"), width=8, height=6)




# predictions by NHS region

fit_sanger_multi_predsbynhsregion = data.frame(emmeans(fit4_sanger_nhs_multi, 
                                                  ~ LINEAGE,
                                                  by=c("DATE_NUM", "NHSREGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to, by=1)), 
                                                  mode="prob", df=NA))
fit_sanger_multi_predsbynhsregion$collection_date = as.Date(fit_sanger_multi_predsbynhsregion$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbynhsregion$LINEAGE = factor(fit_sanger_multi_predsbynhsregion$LINEAGE, levels=levels_LINEAGE_plot) 
fit_sanger_multi_predsbynhsregion$NHSREGION = factor(fit_sanger_multi_predsbynhsregion$NHSREGION, levels=levels_NHSREGION) 
# predicted incidence in different parts of England today
# fit_sanger_multi_predsbynhsregion[fit_sanger_multi_predsbynhsregion$collection_date==today,]

fit_sanger_multi_preds = data.frame(emmeans(fit4_sanger_nhs_multi, 
                                           ~ LINEAGE,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to, by=4)), 
                                           mode="prob", df=NA))
fit_sanger_multi_preds$collection_date = as.Date(fit_sanger_multi_preds$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_preds$LINEAGE = factor(fit_sanger_multi_preds$LINEAGE, levels=levels_LINEAGE_plot) 

muller_sanger_mfit = ggplot(data=fit_sanger_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data)")
muller_sanger_mfit

library(ggpubr)
ggarrange(muller_sanger_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sanger_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots_multinom fit multipanel.png"), width=8, height=8)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots_multinom fit multipanel.pdf"), width=8, height=6)


# predictions by NHS region
muller_sangerbynhsregion_mfit = ggplot(data=fit_sanger_multi_predsbynhsregion, 
                                  aes(x=collection_date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ NHSREGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(sanger$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)")
muller_sangerbynhsregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit.pdf"), width=8, height=6)

library(ggpubr)
ggarrange(muller_sangerbynhsregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_sangerbynhsregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit multipanel.png"), width=8, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_muller plots by NHS region_multinom fit multipanel.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & 95% CONFIDENCE INTERVALS

# # PHE SGTF data to superimpose on model predictions as independent estimate of prop of B.1.1.7 (these data are not up to date anymore)
# sgtf_UK = read.csv(".\\data\\SGTF_UK\\sgtf_2021-05-31_by adm region.csv") # sgtf_2021-05-11_by NHS region, TO DO : ask Nick for update
# colnames(sgtf_UK) = c("collection_date", "ereg_code", "sgtf", "non_sgtf", "REGION")
# sgtf_UK$ereg_code = NULL
# sgtf_UK$collection_date = as.Date(sgtf_UK$collection_date)
# sgtf_UK$DATE_NUM = as.numeric(sgtf_UK$collection_date)
# sgtf_UK$REGION = factor(sgtf_UK$REGION, levels=levels_REGION)
# sgtf_UK = sgtf_UK[sgtf_UK$collection_date>="2021-01-01",]
# sgtf_UK$total = sgtf_UK$sgtf + sgtf_UK$non_sgtf
# sgtf_UK = sgtf_UK[sgtf_UK$total!=0,]
# sgtf_UK$Week = lubridate::week(sgtf_UK$collection_date)
# sgtf_UK$LINEAGE = "B.1.1.7 (S dropout)"
# sgtf_UK$LINEAGE = factor(sgtf_UK$LINEAGE, levels=c(levels_LINEAGE_plot, "B.1.1.7 (S dropout)"))
# sgtf_UK2 = sgtf_UK
# sgtf_UK2$count = sgtf_UK2$sgtf 
# sgtf_UK2$collection_date_num = sgtf_UK2$DATE_NUM
# sgtf_UK2$DATE_NUM = NULL
# sgtf_UK2$sgtf = NULL
# sgtf_UK2$non_sgtf = NULL
# sgtf_UK2$prop = sgtf_UK2$count / sgtf_UK2$total
# head(sgtf_UK2)  
# range(sgtf_UK$collection_date) # "2021-01-01" "2021-05-30"
# library(tidyr)
# sgtf_UK_long = gather(sgtf_UK, Sdropout, count, sgtf:non_sgtf, factor_key=TRUE)
# sgtf_UK_long$collection_date
# 
# 
# # plot of SGTF data by region
# ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-01-01"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
#   facet_wrap(~REGION, scales="free_y") +
#   geom_col(position="fill") +
#   scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                  lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                    lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   ylab("Share") +
#   ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN ENGLAND\n(data PHE)") +
#   scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
#                      labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
#                      limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
#   xlab("Collection date")
#   
# ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by ONS region_stacked as proportion.png"), width=10, height=8)
# # ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data by ONS region_stacked as proportion.pdf"), width=10, height=8)
# 
# ggplot(data=sgtf_UK_long[sgtf_UK_long$collection_date>=as.Date("2021-01-01"),], aes(x=collection_date, y=count, fill=Sdropout, colour=Sdropout)) +
#   # facet_wrap(~REGION, scales="free_y") +
#   geom_col(position="fill") +
#   scale_fill_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                  lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   scale_colour_manual("", values=c(lineage_cols1[which(levels_LINEAGE_plot=="B.1.1.7")],
#                                    lineage_cols1[which(levels_LINEAGE_plot=="B.1.617.2")]), labels=c("S dropout (B.1.1.7 Kent variant)", "S positive (non-B.1.1.7, mostly B.1.617.2)")) +
#   ylab("Share") +
#   ggtitle("INCREASE IN S GENE POSITIVITY (SHARE OF NON-KENT VARIANTS) IN ENGLAND\n(data PHE)") +
#   scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
#                      labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
#                      limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
#   xlab("Collection date")
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.png"), width=10, height=8)
# # ggsave(file=paste0(".\\plots\\",plotdir,"\\SGTF data_stacked as proportion.pdf"), width=10, height=8)
# 
# 
# 
# # logistic spline fit to SGTF data
# 
# fit_SGTF = glm(cbind(sgtf,non_sgtf) ~ REGION * ns(DATE_NUM, df=3), family=binomial, data=sgtf_UK) # at level of ONS regions
# BIC(fit_SGTF)
# 
# SGTF_predsbyregion = data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
#                                                       at=list(DATE_NUM=seq(date.from, date.to)), 
#                                                       type="response"))
# SGTF_predsbyregion$collection_date = as.Date(SGTF_predsbyregion$DATE_NUM, origin="1970-01-01")
# SGTF_predsbyregion$LINEAGE = "B.1.1.7 (S dropout)"
# SGTF_predsbyregion$LINEAGE = factor(SGTF_predsbyregion$LINEAGE, levels=c(levels_LINEAGE_plot, "B.1.1.7 (S dropout)"))
# 
# # fitted S dropout in different parts of England today
# 1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM,
#         at=list(DATE_NUM=today_num), 
#         type="response"))[,c(2,6,5)]
# #        prob asymp.UCL asymp.LCL
# # 1 0.9997718 0.9996945 0.9998295
# 
# # here given by NHS region:
# data.frame(REGION=data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
#                    at=list(DATE_NUM=today_num), 
#                    type="response"))[,2], 
#            1-as.data.frame(emmeans(fit_SGTF, ~ DATE_NUM, by=c("REGION"),
#                         at=list(DATE_NUM=today_num), 
#                         type="response"))[,c(3,7,6)])
# # REGION      prob asymp.UCL asymp.LCL
# # 1                   London 0.9924462 0.9883269 0.9951190
# # 2               North West 0.9999878 0.9999832 0.9999911
# # 3               South West 0.9998023 0.9984491 0.9999748
# # 4               South East 0.9998247 0.9996719 0.9999063
# # 5          East of England 0.9996771 0.9994547 0.9998088
# # 6            East Midlands 0.9995856 0.9993293 0.9997439
# # 7            West Midlands 0.9997396 0.9995751 0.9998404
# # 8               North East 0.9998467 0.9996065 0.9999403
# # 9 Yorkshire and The Humber 0.9999026 0.9998289 0.9999445
# 
# 

# plot predictions multinomial fit by ONS region
# on logit scale:

fit_sanger_multi_preds2 = fit_sanger_multi_predsbyregion # fit_sanger_multi_preds
ymin = 0.001
ymax = 0.9998
fit_sanger_multi_preds2$asymp.LCL[fit_sanger_multi_preds2$asymp.LCL<ymin] = ymin
fit_sanger_multi_preds2$asymp.UCL[fit_sanger_multi_preds2$asymp.UCL<ymin] = ymin
fit_sanger_multi_preds2$asymp.UCL[fit_sanger_multi_preds2$asymp.UCL>ymax] = ymax
fit_sanger_multi_preds2$prob[fit_sanger_multi_preds2$prob<ymin] = ymin

fit_sanger_multi_predsbyregion2 = fit_sanger_multi_predsbyregion # fit_sanger_multi_predsbyregion
fit_sanger_multi_predsbyregion2$asymp.LCL[fit_sanger_multi_predsbyregion2$asymp.LCL<ymin] = ymin
fit_sanger_multi_predsbyregion2$asymp.UCL[fit_sanger_multi_predsbyregion2$asymp.UCL<ymin] = ymin
fit_sanger_multi_predsbyregion2$asymp.UCL[fit_sanger_multi_predsbyregion2$asymp.UCL>ymax] = ymax
fit_sanger_multi_predsbyregion2$prob[fit_sanger_multi_predsbyregion2$prob<ymin] = ymin

# plots of multinomial fit to Sanger Inst data by ONS region

plot_sanger_mfit_logit = qplot(data=fit_sanger_multi_predsbyregion2, 
                               x=collection_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=LINEAGE
  # ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data,\nmultinomial fit LINEAGE~REGION*ns(DATE,df=3))") +
  scale_x_continuous(breaks=seq(as.Date("2021-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2021-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2021-12-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(0.4), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.001, 3), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-12-01"),NA), ylim=c(0.001, 0.9991), expand=c(0,0))
plot_sanger_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.pdf"), width=8, height=6)
library(svglite)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_logit scale.svg"), width=8, height=6)

# on response scale:
plot_sanger_mfit = qplot(data=fit_sanger_multi_predsbyregion2, 
                         x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-09-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(1), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_sanger_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_multinom fit by ONS region_response scale.pdf"), width=8, height=6)


# plots of multinomial fit to Sanger Inst data by ONS region with PHE S dropout data overlaid

plot_sanger_mfit_logit = qplot(data=rbind(fit_sanger_multi_predsbyregion2[,1:9], SGTF_predsbyregion[,1:9]), 
                               x=collection_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(1), pch=I(16)) +
  geom_point(data=sgtf_UK2,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(0.5), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.01, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA), ylim=c(0.001, 0.9991), expand=c(0,0))
plot_sanger_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.pdf"), width=8, height=6)
library(svglite)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logit scale.svg"), width=8, height=6)
# write.csv(fit_sanger_multi_predsbyregion2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_multinomial fit sanger inst data.csv"), row.names=F)
# write.csv(SGTF_predsbyregion, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_logistic fit S dropout data.csv"), row.names=F)
# write.csv(data_agbyweekregion1, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_sanger inst data aggregated by week.csv"), row.names=F)
# write.csv(sgtf_UK2, file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit by ONS region_S dropout data.csv"), row.names=F)

# on response scale:
plot_sanger_mfit = qplot(data=rbind(fit_sanger_multi_predsbyregion2[,1:9], SGTF_predsbyregion[,1:9]), 
                         x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN ENGLAND\n(Sanger Institute baseline surveillance & PHE S dropout data, multinomial fit)") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2020-09-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(1), pch=I(16)) +
  geom_point(data=sgtf_UK2,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(0.5), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(0.1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_sanger_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger S dropout_multinom fit_response scale.pdf"), width=8, height=6)


# PROJECT MULTINOMIAL FIT ONTO CASE DATA ####

# case data from https://coronavirus.data.gov.uk/details/download

# cases_uk_ltla = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesByPublishDateRollingRate&metric=newCasesBySpecimenDate&format=csv")
# cases_uk_ltla$date = as.Date(cases_uk_ltla$date)
# cases_uk_ltla$DATE_NUM = as.numeric(cases_uk_ltla$date)
# cases_uk_ltla$LTLA = factor(cases_uk_ltla$areaCode)
# cases_uk_ltla$areaCode = NULL
# cases_uk_ltla$REGION = LTLAs_regions$RGN20NM[match(cases_uk_ltla$LTLA, LTLAs_regions$LAD20CD)]
# cases_uk_ltla = cases_uk_ltla[!is.na(cases_uk_ltla$REGION),]
# cases_uk_ltla$REGION = factor(cases_uk_ltla$REGION, levels=levels_REGION)

# daily new cases by region :
cases_uk_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDate&format=csv")
cases_uk_region$date = as.Date(cases_uk_region$date)
cases_uk_region$DATE_NUM = as.numeric(cases_uk_region$date)
cases_uk_region$REGION = factor(cases_uk_region$areaName, levels=levels_REGION)
cases_uk_region$areaName = NULL
cases_uk_region = cases_uk_region[cases_uk_region$date<=(max(cases_uk_region$date)-3),] # cut off data from last 3 days (incomplete)

ggplot(data=cases_uk_region, 
       aes(x=date, y=newCasesBySpecimenDate, group=REGION, colour=REGION, fill=REGION)) +
  facet_wrap(~REGION) +
  geom_point(cex=I(0.2)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  scale_y_log10() + ylab("New cases (per day)") + ggtitle("CONFIRMED SARS-CoV2 CASES PER DAY IN ENGLAND BY ONS REGION")

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases England by ONS region.png"), width=8, height=6)

# daily  new cases by ONS region & age group :
cases_uk_age_region = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDateAgeDemographics&format=csv")
cases_uk_age_region$date = as.Date(cases_uk_age_region$date)
cases_uk_age_region$DATE_NUM = as.numeric(cases_uk_age_region$date)
cases_uk_age_region$REGION = factor(cases_uk_age_region$areaName, levels=levels_REGION)
cases_uk_age_region$areaName = NULL

ggplot(data=cases_uk_age_region[cases_uk_age_region$age %in% c("00_59","60+"),], 
        aes(x=date, y=cases, group=age, colour=age, fill=age)) +
          geom_point(cex=I(0.2)) +
          geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
          facet_wrap(~REGION) + scale_y_log10()

ggplot(data=cases_uk_age_region[cases_uk_age_region$age %in% c("00_59","60+"),], 
       aes(x=date, y=rollingSum/7, group=age, colour=age, fill=age)) +
  geom_point(cex=I(0.2)) +
  # geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  facet_wrap(~REGION) + scale_y_log10()

# daily new cases by age group for England :
cases_ENG_age = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDateAgeDemographics&format=csv")
cases_ENG_age$date = as.Date(cases_ENG_age$date)
cases_ENG_age$DATE_NUM = as.numeric(cases_ENG_age$date)
cases_ENG_age$areaName = NULL

ggplot(data=cases_ENG_age[cases_ENG_age$age %in% c("00_59","60+"),], 
       aes(x=date, y=cases, group=age, colour=age, fill=age)) +
  geom_point(cex=I(0.2)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  scale_y_log10()

# daily hospital admissions by NHS region
hosps_uk_nhsregion = read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=newAdmissions&format=csv") # &metric=newCasesBySpecimenDate gives NAs
hosps_uk_nhsregion$date = as.Date(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$DATE_NUM = as.numeric(hosps_uk_nhsregion$date)
hosps_uk_nhsregion$REGION = factor(hosps_uk_nhsregion$areaName, levels=levels_NHSREGION)

ggplot(data=hosps_uk_nhsregion, 
       aes(x=date, y=newAdmissions, group=REGION, colour=REGION, fill=REGION)) +
  facet_wrap(~REGION) +
  geom_point(cex=I(0.2)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 30)) + 
  scale_y_log10() + 
  ylab("New hospital admissions (per day)") + 
  ggtitle("HOSPITAL ADMISSIONS PER DAY IN ENGLAND BY NHS REGION")

ggsave(file=paste0(".\\plots\\",plotdir,"\\hospital admissions England by NHS region.png"), width=8, height=6)


# daily hospital admission also see 
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/COVID-19-daily-admissions-and-beds-20210611.xlsx

# daily hospital admission by age for England
  # https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/Covid-Publication-10-06-2021-Supplementary-Data.xlsx
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/07/Covid-Publication-08-07-2021-Supplementary-Data.xlsx
# PS Next publication: Thursday 8 July 2021
library(rio)
# dat = rio::import(file = "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/07/Covid-Publication-08-07-2021-Supplementary-Data.xlsx",
#            which = 1)
dat = rio::import(file = "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/08/Covid-Publication-12-08-2021-Supplementary-Data.xlsx",
                  which = 1)
dat[,1] = NULL
dates = as.Date(as.vector(unlist(dat[11,-1])), origin="1899-12-30")
ages = gsub("Total reported admissions and diagnoses  ", "", dat[14:(nrow(dat)-3),1])
dat = dat[14:(nrow(dat)-3),-1]
colnames(dat) = dates
dat = data.frame(age=ages, dat, check.names=F)
library(tidyr)
hosps_ENG_age = gather(dat, date, newAdmissions, all_of(as.character(dates)), factor_key=TRUE)
hosps_ENG_age = hosps_ENG_age[hosps_ENG_age$age!="Unknown age",]
hosps_ENG_age$age = factor(hosps_ENG_age$age, levels=c("0-5","6-17","18-24","25-34", "35-44", "45-54","55-64","65-74","75-84","85+"))
hosps_ENG_age$date = as.Date(as.character(hosps_ENG_age$date))

ggplot(data=hosps_ENG_age, 
       aes(x=date, y=newAdmissions, group=age, colour=age, fill=age)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  geom_point(cex=I(0.8)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 35)) + 
  # facet_wrap(~ type) + 
  scale_y_log10() + ylab("Daily new confirmed hospital admissions") + 
  ggtitle("Daily new SARS-CoV2 hospitalisations in England\n(data NHS)") +
  scale_colour_discrete(guide = guide_legend(reverse = TRUE) ) +
  scale_fill_discrete(guide = guide_legend(reverse = TRUE) )

# combined with daily new cases by age (matched to closest age category)
cases_hosps_ENG_age = rbind(data.frame(type="hospitalisations",hosps_ENG_age), 
                            data.frame(type="cases",hosps_ENG_age))
colnames(cases_hosps_ENG_age)[4] = "count"
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="00_04"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                    cases_ENG_age[cases_ENG_age$age=="00_04","date"])]  
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="6-17"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="05_09"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                    cases_ENG_age[cases_ENG_age$age=="05_09","date"])] +
                                                                                               cases_ENG_age$cases[cases_ENG_age$age=="10_14"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                     cases_ENG_age[cases_ENG_age$age=="10_14","date"])] +
                                                                                               cases_ENG_age$cases[cases_ENG_age$age=="15_19"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                     cases_ENG_age[cases_ENG_age$age=="15_19","date"])]
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="18-54"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="20_24"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="20_24","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="25_29"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="25_29","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="30_34"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="30_34","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="35_39"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="35_39","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="40_44"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="40_44","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="45_49"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="45_49","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="50_54"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="50_54","date"])] 
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="55-64"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="55_59"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="55_59","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="60_64"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="60_64","date"])]
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="65-74"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="65_69"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="65_69","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="70_74"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="70_74","date"])] 
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="75-84"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="75_79"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="75_79","date"])] +
                                                                                                cases_ENG_age$cases[cases_ENG_age$age=="80_84"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="80_84","date"])] 
cases_hosps_ENG_age$count[cases_hosps_ENG_age$age=="85+"&cases_hosps_ENG_age$type=="cases"] = cases_ENG_age$cases[cases_ENG_age$age=="85_89"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                      cases_ENG_age[cases_ENG_age$age=="85_89","date"])] +
                                                                                              cases_ENG_age$cases[cases_ENG_age$age=="90+"][match(cases_hosps_ENG_age[cases_hosps_ENG_age$age=="0-5"&cases_hosps_ENG_age$type=="cases","date"],
                                                                                                                                                    cases_ENG_age[cases_ENG_age$age=="90+","date"])]
cases_hosps_ENG_age$age = factor(cases_hosps_ENG_age$age, levels=c("0-5","6-17","18-54","55-64","65-74","75-84","85+"))
ggplot(data=cases_hosps_ENG_age, 
       aes(x=date, y=count, group=age, colour=age, fill=age)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  geom_point(cex=I(0.8)) +
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 35)) + 
  facet_wrap(~ type) + scale_y_log10() + ylab("Daily new confirmed cases & hospital admissions") + 
  ggtitle("Daily new SARS-CoV2 cases & hospitalisations in England\n(data gov.uk & NHS)")

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.pdf"), width=8, height=6)
# write.csv(cases_hosps_ENG_age, file=paste0(".\\plots\\",plotdir,"\\confirmed cases hospitalisations by age England.csv"), row.names=F)

cases_hosps_ENG_age2 = cases_hosps_ENG_age
lag = 7 # lag = 7 from cases to hosps gives best BIC in Poisson GLM
cases_hosps_ENG_age2$date[cases_hosps_ENG_age2$type=="cases"] = cases_hosps_ENG_age2$date[cases_hosps_ENG_age2$type=="cases"]-lag
cases_hosps_ENG_age_wide = spread(cases_hosps_ENG_age2, type, count)
cases_hosps_ENG_age_wide$DATE_NUM = as.numeric(cases_hosps_ENG_age_wide$date)
cases_hosps_ENG_age_wide = cases_hosps_ENG_age_wide[complete.cases(cases_hosps_ENG_age_wide),]
cases_hosps_ENG_age_wide$CHR = cases_hosps_ENG_age_wide$hospitalisations / cases_hosps_ENG_age_wide$cases
fit = glm(hospitalisations~cases*age*DATE_NUM, data=cases_hosps_ENG_age_wide, family=poisson)
BIC(fit)
qplot(data=cases_hosps_ENG_age_wide, x=date, y=CHR*100, group=age, colour=age, fill=age, geom=c("blank")) + 
  geom_smooth(method="gam", se=TRUE, formula = y ~ s(x, k = 10)) + 
  # facet_wrap(~age) +
  # scale_x_log10() + 
  scale_y_log10() +
  ylab("Case (-7 days) to hospitalisation ratio (%)") +
  ggtitle("CASE TO HOSPITALISATION RATIO IN ENGLAND\n(data gov.uk & NHS)")
ggsave(file=paste0(".\\plots\\",plotdir,"\\case to hospitalisation ratio England.png"), width=8, height=6)
cases_hosps_ENG_age_wide[cases_hosps_ENG_age_wide$date==max(cases_hosps_ENG_age_wide$date),]


# MAP MULTINOMIAL FIT ONTO CASE NUMBERS ####

newdat = expand.grid(REGION=levels_REGION,
                     DATE_NUM=seq(date.from, date.to))
fit_sanger_multi_predsbyregion_day = data.frame(newdat,
                                          predict(fit4_sanger_multi, 
                                                  newdata = newdat,
                                                  type = "prob"), check.names=F)  
fit_sanger_multi_predsbyregion_day = gather(fit_sanger_multi_predsbyregion_day, LINEAGE, prob, all_of(levels_LINEAGE_plot))
fit_sanger_multi_predsbyregion_day$collection_date = as.Date(fit_sanger_multi_predsbyregion_day$DATE_NUM, origin="1970-01-01")
fit_sanger_multi_predsbyregion_day$LINEAGE = factor(fit_sanger_multi_predsbyregion_day$LINEAGE, levels=levels_LINEAGE_plot)

fit_sanger_multi_predsbyregion_day$totcases_rollingmean = cases_uk_region$newCasesBySpecimenDateRollingSum[match(interaction(fit_sanger_multi_predsbyregion_day$collection_date,
                                                                                                                         fit_sanger_multi_predsbyregion_day$REGION),
                                                                                        interaction(cases_uk_region$date,
                                                                                                    cases_uk_region$REGION))]/7
fit_sanger_multi_predsbyregion_day$totcases = cases_uk_region$newCasesBySpecimenDate[match(interaction(fit_sanger_multi_predsbyregion_day$collection_date,
                                                                                                       fit_sanger_multi_predsbyregion_day$REGION),
                                                                                                 interaction(cases_uk_region$date,
                                                                                                             cases_uk_region$REGION))]
fit_sanger_multi_predsbyregion_day$WEEKDAY = factor(weekdays(fit_sanger_multi_predsbyregion_day$collection_date))
library(mgcv)
gamfittotcases = gam(totcases ~ s(DATE_NUM, bs="cs", k=25, by=REGION) + REGION + WEEKDAY, family=poisson(log), 
                     data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$collection_date<=max(sanger$WeekEndDate),])
BIC(gamfittotcases) 
gamfitemmeans = as.data.frame(emmeans(gamfittotcases, ~ DATE_NUM, by=c("REGION","DATE_NUM"),
                                                                         at=list(REGION=levels(fit_sanger_multi_predsbyregion_day$REGION),
                                                                           DATE_NUM=seq(min(fit_sanger_multi_predsbyregion_day$DATE_NUM),
                                                         max(fit_sanger_multi_predsbyregion_day$DATE_NUM)+14)), type="response"))
fit_sanger_multi_predsbyregion_day$totcases_smoothed = gamfitemmeans$rate[match(interaction(fit_sanger_multi_predsbyregion_day$REGION,
                                                                                            fit_sanger_multi_predsbyregion_day$DATE_NUM),
                                                                            interaction(gamfitemmeans$REGION,gamfitemmeans$DATE_NUM))]
fit_sanger_multi_predsbyregion_day$cases = fit_sanger_multi_predsbyregion_day$totcases * fit_sanger_multi_predsbyregion_day$prob
fit_sanger_multi_predsbyregion_day$cases_rollingmean = fit_sanger_multi_predsbyregion_day$totcases_rollingmean * fit_sanger_multi_predsbyregion_day$prob
fit_sanger_multi_predsbyregion_day$cases_smoothed = fit_sanger_multi_predsbyregion_day$totcases_smoothed * fit_sanger_multi_predsbyregion_day$prob
fit_sanger_multi_predsbyregion_day$REGION = factor(fit_sanger_multi_predsbyregion_day$REGION, levels=levels_REGION)
fit_sanger_multi_predsbyregion_day$cases[fit_sanger_multi_predsbyregion_day$cases<=0.001] = NA
fit_sanger_multi_predsbyregion_day$cases_smoothed[fit_sanger_multi_predsbyregion_day$cases_smoothed<=0.001] = NA

cases_uk_region$totcases_smoothed = gamfitemmeans$rate[match(interaction(cases_uk_region$REGION,cases_uk_region$DATE_NUM),
                                                             interaction(gamfitemmeans$REGION,gamfitemmeans$DATE_NUM))]

fit_sanger_multi_predsbyregion_day$cases[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2"&
                                           fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA
fit_sanger_multi_predsbyregion_day$cases_smoothed[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2"&
                                           fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA
fit_sanger_multi_predsbyregion_day$cases_rollingmean[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2"&
                                                    fit_sanger_multi_predsbyregion_day$collection_date<as.Date("2021-03-14")] = NA

# plot new cases per day by region
ggplot(data=fit_sanger_multi_predsbyregion_day,
       aes(x=collection_date, y=cases_smoothed, 
           group=LINEAGE)) +
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_line(aes(lwd=I(1.3), colour=LINEAGE)) +
#  geom_smooth(aes(lwd=I(1.3), colour=LINEAGE), method = "gam", family = "poisson", formula = y ~ s(x, k=4), se=FALSE) +
# geom_smooth(aes(lwd=I(1), colour=LINEAGE), method="loess", span=0.8, se=FALSE) +
  geom_line(data=data.frame(collection_date=cases_uk_region$date,
                            cases_smoothed=cases_uk_region$totcases_smoothed, 
                              REGION=cases_uk_region$REGION,
                              LINEAGE="total"), aes(lwd=I(1.9), alpha=I(0.5)), 
              colour=alpha(I("black"),0.7)) +
#  geom_smooth(data=data.frame(collection_date=cases_uk_region$date,
#                              cases=cases_uk_region$newCasesBySpecimenDate, 
#                              REGION=cases_uk_region$REGION,
#                              LINEAGE="total"), aes(lwd=I(1.9), alpha=I(0.5)), 
#              method = "gam", family = "poisson", formula = y ~ s(x, k=4), se=FALSE, colour=alpha(I("black"),0.7)) +
  #  geom_smooth(data=data.frame(collection_date=cases_uk_region$date,
#                              cases=cases_uk_region$newCasesBySpecimenDate, 
#                              REGION=cases_uk_region$REGION,
#                              LINEAGE="total"), aes(lwd=I(1.5)), method="loess", span=0.8, se=FALSE, colour=I("black")) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_y_log10() +
  scale_fill_manual("variant", values=c(lineage_cols1,I("black"))) +
  scale_colour_manual("variant", values=c(lineage_cols1,I("black"))) 

# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day by lineage multinomial fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day by lineage multinomial fit.pdf"), width=8, height=6)


ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) 

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases_smoothed+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit smoothed case data_log10 Y scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases_rollingmean, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data.pdf"), width=8, height=6)

ggplot(data=fit_sanger_multi_predsbyregion_day, 
       aes(x=collection_date, y=cases_rollingmean+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND\n(blue=B.1.1.7 (alpha), magenta=B.1.617.2 (delta))") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_uk_region$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit rolling mean case data_log10 Y scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day stacked area multinomial fit raw case data.pdf"), width=8, height=6)

# ggplot(data=fit_sanger_multi_predsbyregion_day, 
#        aes(x=collection_date, y=cases_smoothed, group=LINEAGE)) + 
#   facet_wrap(~ REGION, scale="free", ncol=3) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                      limits=as.Date(c("2021-03-01","2021-06-30")), expand=c(0,0)) +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("New confirmed cases per day") +
#   ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN ENGLAND") +
#   scale_fill_manual("variant", values=lineage_cols1) +
#   scale_colour_manual("variant", values=lineage_cols1) 
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_cases per day smoothed stacked area multinomial fit.pdf"), width=8, height=6)



# CALCULATE Re VALUES PER VARIANT ####

# Function to calculate Re values from intrinsic growth rate
# cf. https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (assuming gamma distributed gen time)
Re.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) { # Nishiura et al. 2020, or use values from Ferretti et al. 2020 (gamma_mean=5.5, gamma_sd=1.8)
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}
library(mgcv)

qplot(data=fit_sanger_multi_predsbyregion_day, x=collection_date, y=totcases, geom="line") +
  facet_wrap(~ REGION) + scale_y_log10()

fit_sanger_multi_predsbyregion_day$weekday = factor(weekdays(fit_sanger_multi_predsbyregion_day$collection_date))
k=15
fit_cases = gam(totcases ~ s(DATE_NUM, bs="cs", k=k, m=c(2), by=REGION) + REGION + weekday,
                method = "REML",
                knots = list(DATE_NUM = c(min(fit_sanger_multi_predsbyregion_day$DATE_NUM)-14,
                                          seq(min(fit_sanger_multi_predsbyregion_day$DATE_NUM)+1*diff(range(fit_sanger_multi_predsbyregion_day$DATE_NUM))/(k-2), 
                                              max(fit_sanger_multi_predsbyregion_day$DATE_NUM)-1*diff(range(fit_sanger_multi_predsbyregion_day$DATE_NUM))/(k-2), length.out=k-2),
                                          max(fit_sanger_multi_predsbyregion_day$DATE_NUM)+14)),
                family=poisson(log), data=fit_sanger_multi_predsbyregion_day[fit_sanger_multi_predsbyregion_day$LINEAGE==levels_LINEAGE[[1]],],
) # TO DO: add testing
BIC(fit_cases)

# calculate average instantaneous growth rates & 95% CLs using emtrends ####
# based on the slope of the GAM fit on a log link scale
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, var="DATE_NUM", by="REGION",
                              at=list(DATE_NUM=seq(date.from,
                                                   date.to),
                                      REGION=levels_REGION
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
avg_r_cases$REGION = factor(avg_r_cases$REGION, levels=levels_REGION)
qplot(data=avg_r_cases, x=DATE, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue"), group=REGION) +
  facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re VALUES IN ENGLAND BY REGION BASED ON NEW CONFIRMED CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
  # coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants4 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
       function (dat) { 
  wt = as.data.frame(emmeans(fit4_sanger_multi, ~ LINEAGE, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
  # wt = rep(1/length(levels_LINEAGE), length(levels_LINEAGE)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
  cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
  names(cons) = seq_along(cons)
  EMT = emtrends(fit3_sanger_multi,  ~ LINEAGE, by=c("DATE_NUM","REGION"),
               var="DATE_NUM", mode="latent",
               at=list(DATE_NUM=dat,
                       REGION=levels(SGTF_predsbyregion$REGION)))
  out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
  # sum(out$estimate*wt) # should sum to zero
  return(out) } ))

# # for simpler model fit1_sanger_multi
# above_avg_r_variants1 = do.call(rbind, lapply(seq(date.from,
#                                                   date.to), 
#                                               function (dat) { 
#                                                 wt = as.data.frame(emmeans(fit3_sanger_multi, ~ LINEAGE, at=list(DATE_NUM=dat), type="response"))$prob   # important: these should sum to 1
#                                                 # wt = rep(1/length(levels_LINEAGE), length(levels_LINEAGE)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
#                                                 cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
#                                                 names(cons) = seq_along(cons)
#                                                 EMT = emtrends(fit1_sanger_multi,  ~ LINEAGE, by=c("DATE_NUM","REGION"),
#                                                                var="DATE_NUM", mode="latent",
#                                                                at=list(DATE_NUM=dat,
#                                                                        REGION=levels(SGTF_predsbyregion$REGION)))
#                                                 out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
#                                                 # sum(out$estimate*wt) # should sum to zero
#                                                 return(out) } ))
above_avg_r_variants = above_avg_r_variants4
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_LINEAGE), 
                                       labels=levels(data_agbyweeknhsregion1$LINEAGE))
above_avg_r_variants$LINEAGE = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$DATE_NUM = above_avg_r_variants$DATE_NUM-0 # the -7 to calculate back to time of infection
above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
range(above_avg_r_variants$collection_date) # "2020-09-01" "2021-07-31"
above_avg_r_variants$avg_r = avg_r_cases$r[match(interaction(above_avg_r_variants$REGION,above_avg_r_variants$collection_date),
                                                 interaction(avg_r_cases$REGION,avg_r_cases$DATE))]  # average growth rate of all lineages calculated from case nrs
above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
df = data.frame(contrast=NA,
                DATE_NUM=avg_r_cases$DATE_NUM-0, # -7 to calculate back to time of infection
                REGION=avg_r_cases$REGION,
                estimate=NA,
                SE=NA,
                df=NA,
                asymp.LCL=NA,
                asymp.UCL=NA,
                # p.value=NA,
                collection_date=avg_r_cases$DATE-0,
                LINEAGE="avg",
                avg_r=avg_r_cases$r,
                r=avg_r_cases$r,
                r_LOWER=avg_r_cases$r_LOWER,
                r_UPPER=avg_r_cases$r_UPPER,
                Re=avg_r_cases$Re,
                Re_LOWER=avg_r_cases$Re_LOWER,
                Re_UPPER=avg_r_cases$Re_UPPER)
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$LINEAGE = factor(above_avg_r_variants$LINEAGE, levels=c(levels_LINEAGE_plot,"avg"))
above_avg_r_variants$prob = fit_sanger_multi_predsbyregion_day$prob[match(interaction(above_avg_r_variants$DATE_NUM,
                                                                                      above_avg_r_variants$LINEAGE,
                                                                                      above_avg_r_variants$REGION),
                                                                          interaction(fit_sanger_multi_predsbyregion_day$DATE_NUM,
                                                                                      fit_sanger_multi_predsbyregion_day$LINEAGE,
                                                                                      fit_sanger_multi_predsbyregion_day$REGION))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 2
ymin = 1/ymax
above_avg_r_variants2$Re[above_avg_r_variants2$Re>=ymax] = NA
above_avg_r_variants2$Re[above_avg_r_variants2$Re<=ymin] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER>=ymax] = ymax
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$Re_LOWER<=ymin] = ymin
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER>=ymax] = ymax
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$Re_UPPER<=ymin] = ymin
above_avg_r_variants2$Re[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_LOWER[above_avg_r_variants2$prob<0.01] = NA
above_avg_r_variants2$Re_UPPER[above_avg_r_variants2$prob<0.01] = NA
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$LINEAGE %in% c("other"))|(above_avg_r_variants2$collection_date>max(cases_uk_region$date))),], 
      x=collection_date-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE, fill=LINEAGE, alpha=I(0.5),
      group=LINEAGE, linetype=I(0)) +
  facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=LINEAGE, colour=LINEAGE), alpha=I(0.5))
  geom_line(aes(colour=LINEAGE), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=seq(as.Date("2019-12-01"), today, by="month"),
                     labels=substring(months(seq(as.Date("2019-12-01"), today, by="month")),1,1),
                     limits=as.Date(c("2020-09-01",NA)), expand=c(0,0)) +  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN ENGLAND BY REGION\n(case data gov.uk & multinomial fit to Sanger Institute sequence data)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),max(cases_uk_region$date))) +
  scale_fill_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  scale_colour_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
  theme(legend.position="right") # ,  
# axis.title.x=element_blank()

ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant_with clipping.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\sanger_Re values per variant_with clipping.pdf"), width=8, height=6)


