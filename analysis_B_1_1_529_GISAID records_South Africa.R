# ANALYSIS OF GROWTH ADVANTAGE OF B.1.1.529 IN SOUTH AFRICA 
# (GISAID RECORDS & 
# SGTF data (now proxy for B.1.1.529), read off graph by John Burn Murdoch from briefing by Tulio de Oliveira & Richard Lessells, 
# https://www.youtube.com/watch?v=Vh4XMueP1zQ, https://twitter.com/chrischirp/status/1463885565221384202/photo/1)

# T. Wenseleers
# last update 25 NOVEMBER 2021
    
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
    
  today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
  # today = as.Date("2021-09-06")
  today_num = as.numeric(today)
  plotdir = "South Africa_GISAID"
  suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import GISAID records for South Africa
d1 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_11_25_20_subm_2020.tsv"), col_types = cols(.default = "c")) 
d2 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_11_25_20_subm_jan_may_2021.tsv"), col_types = cols(.default = "c")) 
d3 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_11_25_20_subm_june_aug_2021.tsv"), col_types = cols(.default = "c")) 
d4 = read_tsv(file(".//data//GISAID//South Africa//gisaid_hcov-19_2021_11_25_20_sept_nov_2021.tsv"), col_types = cols(.default = "c")) 

# import SGTF data (now proxy for B.1.1.529)
sgtf = read.csv(".//data//GISAID//South Africa//South Africa SGTF.csv") # SGTF data, from graph shown at press conference, now proxy for B.1.1.529
sgtf$date = as.Date(sgtf$date)
plot(sgtf$date, sgtf$SGTF, type="l")
sgtf_nov = sgtf[sgtf$date>=as.Date("2021-11-01"),] # SGTF data for november only (there dominated by B.1.1.529)
B_1_1_529_nov = data.frame(variant="B.1.1.529", date=sgtf_nov$date, "count"=sgtf_nov$SGTF-2) # we subtract a 2% baseline, "count" is actually a percentage here

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
range(GISAID$date) # "2020-03-06" "2021-11-22"
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
GISAID$floor_date = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$date_num = as.numeric(GISAID$date)
# GISAID = GISAID[GISAID$pango_lineage!="None",]
attach(GISAID)

GISAID$variant = case_when(
  grepl("N679K", aa_substitutions) & grepl("H655Y", aa_substitutions) & grepl("P681H", aa_substitutions) ~ "B.1.1.529",
  grepl("B.1.617.2", pango_lineage, fixed=T) | grepl("AY", pango_lineage)  ~ "Delta",
  grepl("B.1.1.7", pango_lineage, fixed=T) ~ "Alpha",
  grepl("B.1.351", pango_lineage, fixed=T) ~ "Beta",
  grepl("C.1.2", pango_lineage, fixed=T) ~ "C.1.2",
  T ~ "Other"
)

table(GISAID[GISAID$variant=="Other"&GISAID$date>as.Date("2021-10-01"),]$pango_lineage)

GISAID = GISAID[!(GISAID$variant=="Other"&GISAID$pango_lineage=="None"),]


# ANALYSIS OF VOCs IN SOUTH AFRICA ####

GISAID_sel = GISAID
nrow(GISAID_sel) # 22333
sum(GISAID_sel$variant==sel_target_VOC) # 66
table(GISAID_sel$variant)
# Alpha B.1.1.529      Beta     C.1.2     Delta     Other 
# 219        66      6763       256     10690      4339  
range(GISAID_sel$date) # "2020-03-06" "2021-11-22"

sel_target_VOC = "B.1.1.529"
sel_reference_VOC = "Delta"

GISAID_sel$variant = factor(GISAID_sel$variant)
GISAID_sel$variant = relevel(GISAID_sel$variant, ref=sel_reference_VOC) # we code Delta as the reference
levels_variant = c(sel_reference_VOC, "Beta", "Alpha", "C.1.2", "Other", sel_target_VOC)
GISAID_sel$variant = factor(GISAID_sel$variant, levels=levels_variant)
table(GISAID_sel$variant, GISAID_sel$`Sampling strategy`)

# for november I extrapolate a multinomial fit to the GISAID & SGTF data and use that to fill in the missing values for the frequencies of the other non-B.1.1.529 lineages
# this would be similar to using an EM algorithm to fill in missing data
set.seed(1)
fit1_southafrica_multi0 = nnet::multinom(variant ~ scale(date_num), data=GISAID_sel[GISAID_sel$variant!="B.1.1.529",], maxit=1000)
fit2_southafrica_multi0 = nnet::multinom(variant ~ ns(date_num, df=2), data=GISAID_sel[GISAID_sel$variant!="B.1.1.529",], maxit=1000)
BIC(fit1_southafrica_multi0, fit2_southafrica_multi0) 
#                         df      BIC
# fit1_southafrica_multi0  8 21196.06
# fit2_southafrica_multi0 12 17802.09 # best

# multinomial model predictions for november, merged with SGTF data for november as proxy for B.1.1.529
predgrid = expand.grid(list(date_num=seq(as.Date("2021-11-01"), max(sgtf$date), by=1)))
data_nov = data.frame(predgrid, as.data.frame(predict(fit2_southafrica_multi0, newdata=predgrid, type="prob")),check.names=F)
data_nov = gather(data_nov, variant, prob, all_of(levels_variant[levels_variant!="B.1.1.529"]), factor_key=TRUE)
data_nov$date = as.Date(data_nov$date_num, origin="1970-01-01")
data_nov$date_num = NULL
data_nov$count = 100*data_nov$prob*(1-B_1_1_529_nov$count/100)
data_nov$prob = NULL
head(data_nov)
data_nov = rbind(data_nov, B_1_1_529_nov)
data_nov$variant = factor(data_nov$variant, levels=levels_variant) 
data_nov$total = 100

# aggregated data to make Muller plots of raw data
# aggregate by day to identify days on which INSA performed (days with a lot of sequences)
# we subset the data to just those days to avoid sampling biases (delta infection clusters etc)
GISAID_sel2 = GISAID_sel[GISAID_sel$date<as.Date("2021-11-01"), c("date","variant")]
data_agbyday = as.data.frame(table(GISAID_sel2$date, GISAID_sel2$variant))
colnames(data_agbyday) = c("date", "variant", "count")
data_agbyday_sum = aggregate(count ~ date, data=data_agbyday, sum)
data_agbyday$total = data_agbyday_sum$count[match(data_agbyday$date, data_agbyday_sum$date)]
data_agbyday$date = as.Date(as.character(data_agbyday$date))
data_agbyday = rbind(data_agbyday, data_nov) # merge with SGTF+extrapolated GISAID for november
data_agbyday$variant = factor(data_agbyday$variant, levels=levels_variant)
data_agbyday$date_num = as.numeric(data_agbyday$date)
data_agbyday$prop = data_agbyday$count/data_agbyday$total
data_agbyday$floor_date = NULL
# GISAID_sel$total_sequenced_on_that_day = data_agbyday$total[match(GISAID_sel$date, data_agbyday$date)]

# # aggregated by week for selected variant lineages
# data_agbyweek = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$variant))
# colnames(data_agbyweek) = c("floor_date", "variant", "count")
# data_agbyweek_sum = aggregate(count ~ floor_date, data=data_agbyweek, sum)
# data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$floor_date, data_agbyweek_sum$floor_date)]
# sum(data_agbyweek[data_agbyweek$variant=="Beta","total"]) == nrow(GISAID_sel) # TRUE
# data_agbyweek$date = as.Date(as.character(data_agbyweek$floor_date))
# data_agbyweek$variant = factor(data_agbyweek$variant, levels=levels_variant)
# data_agbyweek$date_num = as.numeric(data_agbyweek$date)
# data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
# data_agbyweek$floor_date = NULL
# data_agbyweek$date_num = as.numeric(data_agbyweek$date)


# MULLER PLOT (RAW DATA)
  n2 = length(levels(GISAID_sel$variant))
  lineage_cols2 = hcl(h = seq(0, 180, length = n2), l = 60, c = 180)
  lineage_cols2[which(levels(GISAID_sel$variant)=="Alpha")] = "#0085FF"
  lineage_cols2[which(levels(GISAID_sel$variant)=="Beta")] = "green4"
  lineage_cols2[which(levels(GISAID_sel$variant)=="Delta")] = "mediumorchid"
  lineage_cols2[which(levels(GISAID_sel$variant)=="C.1.2")] = "darkorange"
  lineage_cols2[which(levels(GISAID_sel$variant)==sel_target_VOC)] = "red2" # "magenta"
  lineage_cols2[which(levels(GISAID_sel$variant)=="Other")] = "grey65"

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                             labels=substring(months(firststofmonth),1,1),
                             expand=c(0,0))
  
muller_southafrica_raw2 = ggplot(data=data_agbyday, aes(x=date, y=count, group=variant)) + 
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
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID & SGTF data)") 
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
# fit1_southafrica_multi 10 22689.37
# fit2_southafrica_multi 15 19244.09

# growth rate advantage compared to Delta (difference in growth rate per day) 
emtrsouthafrica = emtrends(fit2_southafrica_multi, trt.vs.ctrl ~ variant,  
                   var="date_num",  mode="latent",
                   at=list(date_num=today_num),
                   adjust="none", df=NA)
delta_r_southafrica = data.frame(confint(emtrsouthafrica, 
                                 adjust="none", df=NA)$contrasts, 
                         p.value=as.data.frame(emtrsouthafrica$contrasts,
                                               adjust="none", df=NA)$p.value)
delta_r_southafrica
#            contrast     estimate          SE df    asymp.LCL    asymp.UCL      p.value
# 1      Beta - Delta -0.092397180 0.004389190 NA -0.100999834 -0.083794527 2.235305e-98
# 2     Alpha - Delta -0.092365861 0.007764417 NA -0.107583839 -0.077147883 1.240889e-32
# 3     C.1.2 - Delta -0.002049636 0.002703825 NA -0.007349037  0.003249764 4.484207e-01
# 4     Other - Delta -0.036490645 0.004066734 NA -0.044461298 -0.028519993 2.886498e-19
# 5 B.1.1.529 - Delta  0.380062489 0.024951131 NA  0.331159170  0.428965807 2.159452e-52

exp(0.38*4.7) # B.1.1.529 6x higher R value than Delta, R0 would be 39 then (5.97 x 6.5!!

# fitted prop of different variantS today
multinom_preds_today_avg = data.frame(emmeans(fit2_southafrica_multi, ~ variant|1,
                                              at=list(date_num=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
#     variant         prob           SE df     asymp.LCL    asymp.UCL
# 1     Delta 7.090931e-02 1.650922e-02 NA  3.855182e-02 1.032668e-01
# 2      Beta 1.085318e-08 6.206581e-09 NA -1.311494e-09 2.301786e-08
# 3     Alpha 3.375923e-09 3.517119e-09 NA -3.517503e-09 1.026935e-08
# 4     C.1.2 1.744099e-03 5.058077e-04 NA  7.527338e-04 2.735463e-03
# 5     Other 1.347741e-05 6.970240e-06 NA -1.840066e-07 2.713883e-05
# 6 B.1.1.529 9.273331e-01 1.691559e-02 NA  8.941792e-01 9.604871e-01



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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID & SGTF data, multinomial fit)")
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID & SGTF data, multinomial fit)") +
  xaxis +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  # geom_point(data=data_agbyday,
  #            aes(x=date, y=prop, size=total,
  #                colour=variant
  #            ),
  #            alpha=I(1)) +
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
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN SOUTH AFRICA\n(GISAID & SGTF data, multinomial fit)") +
  xaxis +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",NA)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  geom_point(data=data_agbyday,
             aes(x=date, y=100*prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\ngenotyped", trans="sqrt",
                        range=c(0.5, 5), limits=c(1,max(data_agbyweek$total)), breaks=c(10,100,1000,10000)) +
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

cases_tot = as.data.frame(get_national_data(countries = "South Africa"))
# cases_tot = cases_tot[cases_tot$date>=as.Date("2020-01-01"),]
cases_tot$date_num = as.numeric(cases_tot$date)
# cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
range(cases_tot$date) # "2020-01-03" "2021-11-22"

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
library(mgcv)
k=55
fit_cases = gam(cases_new ~ s(date_num, bs="cs", k=k, m=c(2), fx=F) + 
                  WEEKDAY, # + 
                # BANKHOLIDAY,
                # s(TESTS_ALL, bs="cs", k=8, fx=F),
                family=poisson(log), data=cases_tot,
                method = "REML",
                knots = list(date_num = c(min(cases_tot$date_num)-14,
                                          seq(min(cases_tot$date_num)+0.7*diff(range(cases_tot$date_num))/(k-2), 
                                              max(cases_tot$date_num)-0.7*diff(range(cases_tot$date_num))/(k-2), length.out=k-2),
                                          max(cases_tot$date_num)+14))
) 
BIC(fit_cases)

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_southafrica_multi_preds_withCI$totcases = cases_tot$cases_new[match(round(fit_southafrica_multi_preds_withCI$date_num),cases_tot$date_num)]
fit_southafrica_multi_preds_withCI$cases = fit_southafrica_multi_preds_withCI$totcases * fit_southafrica_multi_preds_withCI$prob
fit_southafrica_multi_preds_withCI$cases[fit_southafrica_multi_preds_withCI$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ date_num, at=list(date_num=seq(date.from, date.to, by=0.5), BANHOLIDAY="no"), type="response"))
fit_southafrica_multi_preds_withCI$smoothed_totcases = cases_emmeans$rate[match(fit_southafrica_multi_preds_withCI$date_num,cases_emmeans$date_num)]
fit_southafrica_multi_preds_withCI$smoothed_cases = fit_southafrica_multi_preds_withCI$smoothed_totcases * fit_southafrica_multi_preds_withCI$prob
fit_southafrica_multi_preds_withCI$smoothed_cases[fit_southafrica_multi_preds_withCI$smoothed_cases<=0.001] = NA
fit_southafrica_multi_preds_withCI$variant = factor(fit_southafrica_multi_preds_withCI$variant, levels=levels_variant)

ggplot(data=fit_southafrica_multi_preds_withCI, 
       aes(x=date, y=cases, group=variant)) + 
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA\n(case data & multinomial fit to GISAID & SGTF data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) # +
  # coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_southafrica_multi_preds_withCI[fit_southafrica_multi_preds_withCI$date<=today,], 
       aes(x=date-7, y=smoothed_cases, group=variant)) + 
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  xaxis +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN SOUTH AFRICA\n(case data & multinomial fit to GISAID & SGTF data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit case data.png"), width=8, height=6)



# DIDN'T TRY TO RUN / UPDATE THE PART BELOW

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
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ date_num, var="date_num", 
                                     at=list(date_num=seq(date.from+16,
                                                          date.to-extrapolate)#,
                                             # BANKHOLIDAY="no"
                                     ), # weekday="Wednesday",
                                     type="link"))
colnames(avg_r_cases)[2] = "r"
colnames(avg_r_cases)[5] = "r_LOWER"
colnames(avg_r_cases)[6] = "r_UPPER"
avg_r_cases$DATE = as.Date(avg_r_cases$date_num, origin="1970-01-01") 
avg_r_cases$Re = Re.from.r(avg_r_cases$r)
avg_r_cases$Re_LOWER = Re.from.r(avg_r_cases$r_LOWER)
avg_r_cases$Re_UPPER = Re.from.r(avg_r_cases$r_UPPER)
avg_r_cases = avg_r_cases[complete.cases(avg_r_cases),]
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) + # -7 TO CALCULATE BACK TO INFECTION DATE
  # facet_wrap(~ REGION) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  xaxis +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN SOUTH AFRICA AT MOMENT OF INFECTION BASED ON NEW CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants0 = do.call(rbind, lapply(seq(date.from+16,
                                                  date.to-extrapolate), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit3_southafrica_multi, ~ variant , at=list(date_num=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_variant)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(fit3_southafrica_multi,  ~ variant , by=c("date_num"),
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
      
