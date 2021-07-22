# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN DENMARK (ECDC TESSy GENOMIC EPIDEMIOLOGY METADATA)
# T. Wenseleers
# last update 22 July 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)
library(lubridate)
library(tsibble)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-07-22")
today_num = as.numeric(today)
plotdir = "Greece_TESSy"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

# import ECDC TESSy variant data from https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea
# https://opendata.ecdc.europa.eu/covid19/virusvariant/csv/data.csv

TESSY = read.csv("https://opendata.ecdc.europa.eu/covid19/virusvariant/csv/data.csv")
colnames(TESSY)[1] = "country"
colnames(TESSY)[6] = "total"
colnames(TESSY)[10] = "count"
TESSY = TESSY[TESSY$source=="TESSy",] # just keep ECDC TESSy data & discard ECDC TESSy data
TESSY = TESSY[TESSY$variant!="UNK",] # remove unknown
TESSY = TESSY[!is.na(TESSY$count),]
TESSY = TESSY[TESSY$valid_denominator=="Yes",]
TESSY$date = as.Date(yearweek(gsub("-","-W",TESSY$year_week)))+3.5 # week midpoint
TESSY$DATE_NUM = as.numeric(TESSY$date)
unique(TESSY$variant)
# "B.1.1.7"         "B.1.1.7+E484K"   "Other"           "B.1.351"         "B.1.525"         "P.1"             "B.1.617.2"       "B.1.617.1"      
# "B.1.621"         "B.1.617"         "B.1.616"         "B.1.620"         "B.1.617.3"       "B.1.427/B.1.429" "P.3"
levels_VARIANTS = c("B.1.1.7","B.1.1.7+E484K","B.1.427/B.1.429","B.1.525",
                    "B.1.351","P.1","P.3","B.1.616","B.1.620","B.1.621",
                    "B.1.617","B.1.617.1","B.1.617.2","B.1.617.3","Other")
TESSY$variant = factor(TESSY$variant, levels=levels_VARIANTS)
TESSY = TESSY[rep(seq_len(nrow(TESSY)), TESSY$count),] # convert to long format
TESSY$count = NULL
head(TESSY)

# ANALYSIS OF VOCs IN GREECE ####

selected_countries = c("Greece")
TESSY_sel = TESSY[TESSY$country %in% selected_countries,]
# use data from Jan  1 onwards
TESSY_sel = TESSY_sel[TESSY_sel$date>=as.Date("2021-01-01"),]
range(TESSY_sel$date) # "2021-01-07" "2021-07-15"
TESSY_sel$variant = droplevels(TESSY_sel$variant)
levels_VARIANTS = levels(TESSY_sel$variant)


# AGGREGATE DATA
# aggregated by week
data_agbyweek = as.data.frame(table(TESSY_sel$date, TESSY_sel$variant))
colnames(data_agbyweek) = c("date", "variant", "count")
data_agbyweek_sum = aggregate(count ~ date, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$date, data_agbyweek_sum$date)]
sum(data_agbyweek[data_agbyweek$variant=="B.1.617.1","total"]) == nrow(TESSY_sel) # correct
data_agbyweek$date = as.Date(as.character(data_agbyweek$date))
data_agbyweek$variant = factor(data_agbyweek$variant, levels=levels_VARIANTS)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
data_agbyweek$DATE_NUM = as.numeric(data_agbyweek$date)

# MULLER PLOT (RAW DATA)
library(scales)
n2 = length(levels(TESSY_sel$variant))
lineage_cols2 = hcl(h = seq(15, 320, length = n2), l = 65, c = 200)
lineage_cols2[which(levels(TESSY_sel$variant)=="B.1.1.7")] = "#0085FF"
lineage_cols2[which(levels(TESSY_sel$variant)=="B.1.351")] = "#9A9D00"
lineage_cols2[which(levels(TESSY_sel$variant)=="P.1")] = "cyan3"
lineage_cols2[which(levels(TESSY_sel$variant)=="B.1.617.1")] = muted("magenta")
lineage_cols2[which(levels(TESSY_sel$variant)=="B.1.617.2")] = "magenta"
lineage_cols2[which(levels(TESSY_sel$variant)=="Other")] = "grey75"

muller_greece_raw2 = ggplot(data=data_agbyweek, aes(x=date, y=count, group=variant)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols2) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 variant LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN GREECE\n(ECDC TESSy data)") 
# +
# coord_cartesian(xlim=c(1,max(TESSY_sel$Week)))
muller_greece_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_muller plots_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_muller plots_raw data.pdf"), width=8, height=6)


# multinomial fits

library(nnet)
library(splines)
set.seed(1)
fit1_greece_multi = nnet::multinom(variant ~ scale(DATE_NUM), weights=count, data=data_agbyweek, maxit=1000)
fit2_greece_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2), weights=count, data=data_agbyweek, maxit=1000)
fit3_greece_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3), weights=count, data=data_agbyweek, maxit=1000)
BIC(fit1_greece_multi, fit2_greece_multi, fit3_greece_multi) 
# df      BIC
# fit1_greece_multi 14 28876.60
# fit2_greece_multi 21 27951.83
# fit3_greece_multi 28 27531.69

# growth rate advantage compared to UK type B.1.1.7 (difference in growth rate per day) 
emtrgreece = emtrends(fit3_greece_multi, trt.vs.ctrl ~ variant,  
                   var="DATE_NUM",  mode="latent",
                   at=list(DATE_NUM=max(TESSY_sel$DATE_NUM)))
delta_r_greece = data.frame(confint(emtrgreece, 
                             adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_greece) = delta_r_greece[,"contrast"]
delta_r_greece = delta_r_greece[,-1]
delta_r_greece
# estimate    asymp.LCL    asymp.UCL
# B.1.525 - B.1.1.7   -0.065643366 -0.172726491  0.041439758
# B.1.351 - B.1.1.7    0.142672451  0.129660800  0.155684102
# P.1 - B.1.1.7        0.053456150 -0.020159755  0.127072054
# B.1.621 - B.1.1.7   -0.403752139 -0.677527926 -0.129976352
# B.1.617.1 - B.1.1.7 -3.430288449 -3.797915818 -3.062661079
# B.1.617.2 - B.1.1.7  0.177823147  0.164699827  0.190946466
# Other - B.1.1.7     -0.003699635 -0.007978155  0.000578884

# implied increase in infectiousness (due to combination of increased transmissibility and/or immune escape)
# assuming generation time of 4.7 days (Nishiura et al. 2020)
exp(delta_r_greece*4.7) 
# estimate    asymp.LCL    asymp.UCL
# B.1.525 - B.1.1.7   7.345305e-01 4.440516e-01 1.215028e+00
# B.1.351 - B.1.1.7   1.955333e+00 1.839338e+00 2.078643e+00
# P.1 - B.1.1.7       1.285624e+00 9.095995e-01 1.817094e+00
# B.1.621 - B.1.1.7   1.499228e-01 4.140386e-02 5.428681e-01
# B.1.617.1 - B.1.1.7 9.957490e-08 1.769110e-08 5.604605e-07
# B.1.617.2 - B.1.1.7 2.306587e+00 2.168616e+00 2.453335e+00
# Other - B.1.1.7     9.827620e-01 9.631970e-01 1.002724e+00

# fitted prop of different LINEAGES in the Greece today
multinom_preds_today_avg = data.frame(emmeans(fit3_greece_multi, ~ variant|1,
                                              at=list(DATE_NUM=today_num), 
                                              mode="prob", df=NA))
multinom_preds_today_avg
# variant          prob            SE df      asymp.LCL     asymp.UCL
# 1   B.1.1.7  1.961844e-02  2.743277e-03 NA   1.424171e-02  2.499516e-02
# 2   B.1.525  1.550954e-06  4.030417e-06 NA  -6.348518e-06  9.450426e-06
# 3   B.1.351  1.066047e-01  1.807852e-02 NA   7.117146e-02  1.420380e-01
# 4       P.1  2.787717e-04  3.364344e-04 NA  -3.806276e-04  9.381710e-04
# 5   B.1.621  2.628011e-09  1.238298e-08 NA  -2.164219e-08  2.689821e-08
# 6 B.1.617.1 8.483729e-106 1.448931e-105 NA -1.991481e-105 3.688226e-105
# 7 B.1.617.2  8.648352e-01  1.986121e-02 NA   8.259079e-01  9.037624e-01
# 8     Other  8.661361e-03  1.321740e-03 NA   6.070797e-03  1.125192e-02

# % non-B.1.1.7
colSums(multinom_preds_today_avg[-1, c("prob","asymp.LCL","asymp.UCL")])
#      prob asymp.LCL asymp.UCL 
# 0.9803816 0.9027632 1.0580000 


# PLOT MULTINOMIAL FIT

# extrapolate = 30
date.from = as.numeric(as.Date("2021-01-01"))
date.to = as.numeric(as.Date("2021-07-31")) # max(TESSY_sel$DATE_NUM)+extrapolate

# multinomial model predictions (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to)))
fit_greece_multi_preds = data.frame(predgrid, as.data.frame(predict(fit3_greece_multi, newdata=predgrid, type="prob")),check.names=F)
library(tidyr)
library(tidyselect)
fit_greece_multi_preds = gather(fit_greece_multi_preds, variant, prob, all_of(levels_VARIANTS), factor_key=TRUE)
fit_greece_multi_preds$date = as.Date(fit_greece_multi_preds$DATE_NUM, origin="1970-01-01")
fit_greece_multi_preds$variant = factor(fit_greece_multi_preds$variant, levels=levels_VARIANTS) 

muller_greece_mfit = ggplot(data=fit_greece_multi_preds, 
                                   aes(x=date, y=prob, group=variant)) + 
  # facet_wrap(~ STATE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols2) +
  annotate("rect", xmin=max(TESSY_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN GREECE\n(ECDC TESSy data, multinomial fit)")
muller_greece_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_muller plots_multinom fit.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_muller plots_multinom fit.pdf"), width=10, height=6)


library(ggpubr)
ggarrange(muller_greece_raw2 + coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_greece_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_muller plots multipanel_multinom fit.png"), width=10, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_muller plots multipanel_multinom fit.pdf"), width=10, height=10)





# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# multinomial model predictions by state with confidence intervals (but slower)
fit_greece_multi_preds_withCI = data.frame(emmeans(fit3_greece_multi,
                                                        ~ variant,
                                                        by=c("DATE_NUM"),
                                                        at=list(DATE_NUM=seq(date.from, date.to, by=1)),  # by=XX to speed up things a bit
                                                        mode="prob", df=NA))
fit_greece_multi_preds_withCI$date = as.Date(fit_greece_multi_preds_withCI$DATE_NUM, origin="1970-01-01")
fit_greece_multi_preds_withCI$variant = factor(fit_greece_multi_preds_withCI$variant, levels=levels_variant)
fit_greece_multi_preds2 = fit_greece_multi_preds_withCI


# on logit scale:

ymin = 0.001
ymax = 0.999
fit_greece_multi_preds2$asymp.LCL[fit_greece_multi_preds2$asymp.LCL<ymin] = ymin
fit_greece_multi_preds2$asymp.UCL[fit_greece_multi_preds2$asymp.UCL<ymin] = ymin
fit_greece_multi_preds2$asymp.UCL[fit_greece_multi_preds2$asymp.UCL>ymax] = ymax
fit_greece_multi_preds2$prob[fit_greece_multi_preds2$prob<ymin] = ymin

plot_greece_mfit_logit = qplot(data=fit_greece_multi_preds2, x=date, y=prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN GREECE\n(ECDC TESSy data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")), ylim=c(0.001, 0.991), expand=c(0,0))
plot_greece_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_multinom fit_logit scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_multinom fit_logit scale.pdf"), width=10, height=6)


# on response scale:
plot_greece_mfit = qplot(data=fit_greece_multi_preds2, x=date, y=100*prob, geom="blank") +
  # facet_wrap(~ STATE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN GREECE\n(ECDC TESSy data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
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
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 5), limits=c(1,max(data_agbyweek$total)), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_greece_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_multinom fit_response scale.png"), width=10, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\greece_multinom fit_response scale.pdf"), width=10, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot = as.data.frame(get_national_data(countries = "Greece"))
cases_tot = cases_tot[cases_tot$date>=as.Date("2020-08-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
# cases_tot$BANKHOLIDAY = bankholiday(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)
range(cases_tot$date)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
library(mgcv)
k=20
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=k, m=c(2), fx=F) + 
                  WEEKDAY, # + 
                # BANKHOLIDAY,
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

fit_greece_multi_preds_withCI$totcases = cases_tot$cases_new[match(round(fit_greece_multi_preds_withCI$DATE_NUM),cases_tot$DATE_NUM)]
fit_greece_multi_preds_withCI$cases = fit_greece_multi_preds_withCI$totcases * fit_greece_multi_preds_withCI$prob
fit_greece_multi_preds_withCI$cases[fit_greece_multi_preds_withCI$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM, at=list(DATE_NUM=seq(date.from, date.to, by=0.5), BANHOLIDAY="no"), type="response"))
fit_greece_multi_preds_withCI$smoothed_totcases = cases_emmeans$rate[match(fit_greece_multi_preds_withCI$DATE_NUM,cases_emmeans$DATE_NUM)]
fit_greece_multi_preds_withCI$smoothed_cases = fit_greece_multi_preds_withCI$smoothed_totcases * fit_greece_multi_preds_withCI$prob
fit_greece_multi_preds_withCI$smoothed_cases[fit_greece_multi_preds_withCI$smoothed_cases<=0.001] = NA

ggplot(data=fit_greece_multi_preds_withCI[fit_greece_multi_preds_withCI$date>=as.Date("2021-02-18"),], 
       aes(x=date, y=cases, group=variant)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),max(cases_tot$date)), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN GREECE\n(case data & multinomial fit to ECDC TESSy data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-02-18"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_greece_multi_preds_withCI[fit_greece_multi_preds_withCI$date>=as.Date("2021-02-18")&
                                              fit_greece_multi_preds_withCI$date<=max(cases_tot$date),], 
       aes(x=date-7, y=smoothed_cases, group=variant)) + 
  # facet_wrap(~ REGION, scale="free", ncol=3) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     # limits=c(as.Date("2021-03-01"),today), 
                     expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN GREECE\n(case data & multinomial fit to ECDC TESSy data)") +
  scale_fill_manual("variant", values=lineage_cols2) +
  scale_colour_manual("variant", values=lineage_cols2) +
  coord_cartesian(xlim=c(as.Date("2021-02-18"),max(cases_tot$date)))

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
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN GREECE AT MOMENT OF INFECTION BASED ON NEW CASES") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants3 = do.call(rbind, lapply(seq(date.from,
                                                  date.to), 
                                              function (d) { 
                                                wt = as.data.frame(emmeans(fit3_greece_multi, ~ variant , at=list(DATE_NUM=d), type="response"))$prob   # important: these should sum to 1
                                                # wt = rep(1/length(levels_variantS), length(levels_variantS)) # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_variant)
                                                cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                                names(cons) = seq_along(cons)
                                                EMT = emtrends(fit3_greece_multi,  ~ variant , by=c("DATE_NUM"),
                                                               var="DATE_NUM", mode="latent",
                                                               at=list(DATE_NUM=d))
                                                out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                                # sum(out$estimate*wt) # should sum to zero
                                                return(out) } ))
above_avg_r_variants = above_avg_r_variants3
above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                       levels=1:length(levels_VARIANTS), 
                                       labels=levels_VARIANTS)
above_avg_r_variants$variant = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
above_avg_r_variants$date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
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
                DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
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
# df = df[df$DATE_NUM<=max(above_avg_r_variants$DATE_NUM)&df$DATE_NUM>=(min(above_avg_r_variants$DATE_NUM)+7),]
above_avg_r_variants = rbind(above_avg_r_variants, df)
above_avg_r_variants$variant = factor(above_avg_r_variants$variant, levels=c(levels_VARIANTS,"avg"))
above_avg_r_variants$prob = fit_greece_multi_preds_withCI$prob[match(interaction(above_avg_r_variants$DATE_NUM,
                                                                      above_avg_r_variants$variant),
                                                          interaction(fit_greece_multi_preds_withCI$DATE_NUM,
                                                                      fit_greece_multi_preds_withCI$variant))]
above_avg_r_variants2 = above_avg_r_variants
ymax = 4
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
qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$variant %in% c("other"))|above_avg_r_variants2$date>max(cases_tot$DATE)),], 
      x=date-7, # -7 to calculate back to date of infection
      y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=variant, fill=variant, alpha=I(0.5),
      group=variant, linetype=I(0)) +
  # facet_wrap(~ REGION) +
  # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
  geom_line(aes(colour=variant), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN GREECE\nAT MOMENT OF INFECTION\n(based on case data & multinomial fit to ECDC TESSy data)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date))) +
  scale_fill_manual("variant", values=c(lineage_cols2,"black")) +
  scale_colour_manual("variant", values=c(lineage_cols2,"black")) +
  theme(legend.position="right") 

ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)
        
