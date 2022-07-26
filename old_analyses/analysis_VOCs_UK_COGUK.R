# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE UK BASED ON COG-UK SEQUENCING DATA ####

# last update 25 July 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
# today = as.Date("2021-08-19")
today_num = as.numeric(today)
plotdir = "UK_COGUK"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

library(data.table)
library(R.utils)
coguk = fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz")
table(coguk$lineage)
sum(coguk$lineage=="BA.2.75") # 26 # 19 in GISAID
levels_REGION = c("England","Scotland","Wales","Northern Ireland")
coguk$REGION = factor(coguk$adm1, levels=c("UK-ENG","UK-SCT","UK-WLS","UK-NIR"), labels=levels_REGION)
coguk = coguk[!is.na(coguk$REGION),]
levels(coguk$REGION)
coguk$date = as.Date(coguk$sample_date)
head(coguk)

# use data from sept 2020 onwards & Pillar 2 (community as opposed to hospitalised) only
#cogukp1 = coguk[coguk$date>=as.Date("2020-09-01")&coguk$is_pillar_2=="N",]
#cogukp2 = coguk[coguk$date>=as.Date("2020-09-01")&coguk$is_pillar_2=="Y",]
cogukp2 = coguk[coguk$date>=as.Date("2020-09-01"),]

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
# cogukp2 = cogukp2[-which(grepl("B.1.617", cogukp2$lineage, fixed=TRUE)&cogukp2$date<=as.Date("2021-04-14")),]  

nrow(cogukp2) # 2745422
range(cogukp2$date) # "2020-09-01" "2022-07-21"
library(lubridate)
cogukp2$WeekStartDate = floor_date(cogukp2$date,unit="week") # lubridate::week(cogukp2$date)
cogukp2$DATE_NUM = as.numeric(cogukp2$date) 
colnames(cogukp2)

sort(unique(cogukp2$lineage))
sel_target_VOC = "BA.2.75"

cogukp2$LINEAGE = case_when(
  cogukp2$lineage=="BA.2.75" ~ "Omicron (BA.2.75)",
  cogukp2$lineage=="BA.2.74" ~ "Omicron (BA.2.74)",
  cogukp2$lineage=="BA.2.76" ~ "Omicron (BA.2.76)",
  (cogukp2$lineage=="B.1.617.2")|grepl("AY", cogukp2$lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", cogukp2$lineage) ~ "Alpha",
  grepl("B.1.351", cogukp2$lineage, fixed=T) ~ "Beta",
  grepl("^BA\\.1$|BA\\.1\\.", cogukp2$lineage) ~ "Omicron (BA.1)",
  grepl("^BA\\.3$|BA\\.3\\.", cogukp2$lineage) ~ "Omicron (BA.3)",
  grepl("^BA\\.4",cogukp2$lineage) ~ "Omicron (BA.4)",
  grepl("^BA\\.5",cogukp2$lineage)|grepl("BE|BF",cogukp2$lineage) ~ "Omicron (BA.5)",
  grepl("^BA\\.2",cogukp2$lineage) ~ "Omicron (BA.2)",
  cogukp2$lineage!="Unassigned" ~ "Other" # remaining Unassigned will be assigned as NA
)

cogukp2 = cogukp2[!is.na(cogukp2$LINEAGE),]
table(cogukp2$LINEAGE)
nrow(cogukp2) # 2729444
sum(cogukp2$LINEAGE == "Omicron (BA.2.75)", na.rm=T) 
# 26 BA.2.75 with valid dates

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

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))




# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of UK
data_agbyweek1 = as.data.frame(table(cogukp2$WeekStartDate, cogukp2$LINEAGE))
colnames(data_agbyweek1) = c("WeekStartDate", "LINEAGE", "count")
data_agbyweek1_sum = aggregate(count ~ WeekStartDate, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$WeekStartDate, data_agbyweek1_sum$WeekStartDate)]
data_agbyweek1$WeekStartDate = as.Date(data_agbyweek1$WeekStartDate)
data_agbyweek1$collection_date = data_agbyweek1$WeekStartDate # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweek1$LINEAGE = factor(data_agbyweek1$LINEAGE, levels=levels_LINEAGE)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and region
data_agbyweekregion1 = as.data.frame(table(cogukp2$WeekStartDate, cogukp2$REGION, cogukp2$LINEAGE))
colnames(data_agbyweekregion1) = c("WeekStartDate", "REGION", "LINEAGE", "count")
data_agbyweekregion1_sum = aggregate(count ~ WeekStartDate + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$WeekStartDate,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$WeekStartDate,data_agbyweekregion1_sum$REGION))]
data_agbyweekregion1$WeekStartDate = as.Date(data_agbyweekregion1$WeekStartDate) # as.numeric(as.character(data_agbyweekregion1$WeekStartDate))
data_agbyweekregion1$collection_date = data_agbyweekregion1$WeekStartDate # lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$WeekStartDate - 1 ) - 3.5 # we use the week midpoint
data_agbyweekregion1$LINEAGE = factor(data_agbyweekregion1$LINEAGE, levels=levels_LINEAGE)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion1[data_agbyweekregion1$collection_date==max(data_agbyweekregion1$collection_date),]

# MULLER PLOT (RAW DATA)
unique(cogukp2$LINEAGE)

library(scales)

data_agbyweek1$LINEAGE = factor(data_agbyweek1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweekregion1$LINEAGE = factor(data_agbyweekregion1$LINEAGE, levels=levels_LINEAGE_plot)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)

library(ggplot2)
library(ggthemes)
muller_cogukp2_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VOCs IN THE UK\n(COG-UK data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_cogukp2_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plot avg_raw data.png"), width=8, height=6)


muller_cogukp2byregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VOCs IN THE UK\n(COG-UK data)") +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_cogukp2byregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
data_agbyweekregion1$LINEAGE = relevel(data_agbyweekregion1$LINEAGE, ref=sel_reference_VOC) # we take BA.5 as baseline / reference level
data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
set.seed(1)
fit1_cogukp2_multi = nnet::multinom(LINEAGE ~ REGION + DATE_NUM, weights=count, data=data_agbyweekregion1, maxit=1000)
fit2_cogukp2_multi = nnet::multinom(LINEAGE ~ REGION * DATE_NUM, weights=count, data=data_agbyweekregion1, maxit=1000)
fit3_cogukp2_multi = nnet::multinom(LINEAGE ~ REGION + ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion1, maxit=1000)
fit4_cogukp2_multi = nnet::multinom(LINEAGE ~ REGION * ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion1, maxit=1000)
BIC(fit1_cogukp2_multi, fit2_cogukp2_multi, fit3_cogukp2_multi, fit4_cogukp2_multi) 
# fit4_cogukp2_multi fits best (lowest BIC), but fit3_cogukp2 almost as good

# growth rate advantages of different VOCs compared to BA.5 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2 = emtrends(fit3_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2 = data.frame(confint(emtrcogukp2, 
                                   adjust="none")$contrasts, 
                            p.value=as.data.frame(emtrcogukp2$contrasts, adjust="sidak")$p.value)
delta_r_cogukp2

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_cogukp2 =  sign(delta_r_cogukp2[,c(2,5,6)])*100*(exp(abs(delta_r_cogukp2[,c(2,5,6)])*4.7)-1)
transmadv_cogukp2 =  data.frame(contrast=delta_r_cogukp2$contrast, transmadv_cogukp2)
transmadv_cogukp2


# growth rate & transmission advantages of different VOCs compared to UK type B.1.1.7 by REGION (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2_region = emtrends(fit4_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE,  
                     var="DATE_NUM",  by=c("REGION"), mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2_region = data.frame(confint(emtrcogukp2_region, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrcogukp2_region$contrasts)$p.value)
delta_r_cogukp2_region


# PLOT MULTINOMIAL FIT ####

extrapolate = 30
date.from = as.numeric(as.Date("2020-09-01"))
date.to = today_num+extrapolate # max(cogukp2$DATE_NUM)+extrapolate

# without CIs
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to),
                            REGION=levels_REGION))
fit_preds = data.frame(predgrid, as.data.frame(predict(fit3_cogukp2_multi, 
                                                       newdata=predgrid, type="prob")),check.names=F)
fit_preds = gather(fit_preds, LINEAGE, prob, all_of(levels_LINEAGE), factor_key=TRUE)
fit_preds$date = as.Date(fit_preds$DATE_NUM, origin="1970-01-01")
fit_preds$LINEAGE = factor(fit_preds$LINEAGE, levels=levels_LINEAGE_plot)
fit_preds$REGION = factor(fit_preds$REGION, levels=levels_REGION)

# # with CIs
# fit_cogukp2_multi_predsbyregion = data.frame(emmeans(fit3_cogukp2_multi, 
#                                                   ~ LINEAGE,
#                                                   by=c("DATE_NUM", "REGION"),
#                                                   at=list(DATE_NUM=seq(date.from, date.to, by=1)), 
#                                                   mode="prob", df=NA))

muller_cogukp2byregion_mfit = ggplot(data=fit_preds, 
                           aes(x=date, y=prob, group=LINEAGE)) + 
  facet_wrap(~ REGION, ncol=2) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  annotate("rect", xmin=max(cogukp2$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VOCs IN THE UK\n(COG-UK data)")
muller_cogukp2byregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plot avg_multinom fit.png"), width=8, height=6)


library(ggpubr)
ggarrange(muller_cogukp2byregion_raw1+coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_cogukp2byregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.png"), width=8, height=10)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.pdf"), width=8, height=6)


# PLOT MODEL FIT WITH DATA & 95% CONFIDENCE INTERVALS

# on logit scale:
fit_cogukp2_multi_preds2 = fit_preds
ymin = 0.00001
ymax = 0.998
fit_cogukp2_multi_preds2$asymp.LCL[fit_cogukp2_multi_preds2$asymp.LCL<ymin] = ymin
fit_cogukp2_multi_preds2$asymp.UCL[fit_cogukp2_multi_preds2$asymp.UCL<ymin] = ymin
fit_cogukp2_multi_preds2$asymp.UCL[fit_cogukp2_multi_preds2$asymp.UCL>ymax] = ymax
fit_cogukp2_multi_preds2$prob[fit_cogukp2_multi_preds2$prob<ymin] = ymin

# plots of multinomial fit to cogukp2 data

plot_cogukp2_mfit_logit = qplot(data=fit_cogukp2_multi_preds2, 
                               x=date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  #geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
  #                fill=LINEAGE
  #), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit variant~ns(date,df=2)+country)") +
  xaxis +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(1), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(1, 6), limits=c(1,10^6), breaks=c(1,10,100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA), 
                  ylim=c(0.00001, 0.998), expand=0)
plot_cogukp2_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.pdf"), width=8, height=6)
# library(svglite)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.svg"), width=8, height=6)

# on response scale:
plot_cogukp2_mfit = qplot(data=fit_cogukp2_multi_preds2, 
                         x=date, y=100*prob, geom="blank") +
  facet_wrap(~REGION) +
  # geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
  #                 fill=LINEAGE
  # ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit variant~ns(date,df=2)+country)") +
  xaxis +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",date.to-1)),
                  ylim=c(0,100), expand=0) +
  scale_fill_manual("variant", values=lineage_cols_plot) +
  scale_colour_manual("variant", values=lineage_cols_plot) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE, 
             ),
             alpha=I(1), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_cogukp2_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_response scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_response scale.pdf"), width=8, height=6)




# PLOTS OF NEW CASES PER DAY BY VARIANT & EFFECTIVE REPRODUCTION NUMBER BY VARIANT THROUGH TIME ####

# load case data
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)  

cases_tot = as.data.frame(get_regional_data(country = "United Kingdom", regions=c("England","Scotland","Wales","Northern Ireland")))
cases_tot = cases_tot[cases_tot$date>=as.Date("2020-08-01"),]
cases_tot$DATE_NUM = as.numeric(cases_tot$date)
cases_tot$WEEKDAY = weekdays(cases_tot$date)
cases_tot$region = factor(cases_tot$region, levels=c("England","Scotland","Wales","Northern Ireland"))
cases_tot = cases_tot[cases_tot$date<=(max(cases_tot$date)-3),] # cut off data from last 3 days (incomplete)

# smooth out weekday effects in case nrs using GAM (if testing data is available one could correct for testing intensity as well)
fit_cases = gam(cases_new ~ s(DATE_NUM, bs="cs", k=20, fx=F, by=region) + region +
                  WEEKDAY # +
                  # s(tested_new, bs="cs", k=8, fx=F, by=region)
                ,
                family=poisson(log), data=cases_tot,
) 
BIC(fit_cases)

# STACKED AREA CHART OF NEW CASES BY VARIANT (MULTINOMIAL FIT MAPPED ONTO CASE DATA) ####

fit_cogukp2_multi_predsbyregion$totcases = cases_tot$cases_new[match(interaction(fit_cogukp2_multi_predsbyregion$DATE_NUM,fit_cogukp2_multi_predsbyregion$REGION),
                                                                     interaction(cases_tot$DATE_NUM,cases_tot$region))]
fit_cogukp2_multi_predsbyregion$cases = fit_cogukp2_multi_predsbyregion$totcases * fit_cogukp2_multi_predsbyregion$prob
fit_cogukp2_multi_predsbyregion$cases[fit_cogukp2_multi_predsbyregion$cases<=0.001] = NA
cases_emmeans = as.data.frame(emmeans(fit_cases, ~ DATE_NUM|region, at=list(DATE_NUM=seq(date.from, date.to, by=0.5)
                                                                     ), type="response"))
fit_cogukp2_multi_predsbyregion$smoothed_totcases = cases_emmeans$rate[match(interaction(fit_cogukp2_multi_predsbyregion$DATE_NUM,fit_cogukp2_multi_predsbyregion$REGION),
                                                                             interaction(cases_emmeans$DATE_NUM,cases_emmeans$region))]
fit_cogukp2_multi_predsbyregion$smoothed_cases = fit_cogukp2_multi_predsbyregion$smoothed_totcases * fit_cogukp2_multi_predsbyregion$prob
fit_cogukp2_multi_predsbyregion$smoothed_cases[fit_cogukp2_multi_predsbyregion$smoothed_cases<=0.001] = NA

ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=cases+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_cogukp2_multi_predsbyregion,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE)) +
  facet_wrap(~ REGION, scale="free") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE, group=LINEAGE), position="stack") +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=smoothed_cases+1, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE=="B.1.1.7",], aes(lwd=I(1.2), colour=NULL), fill=I("#0085FF"), position="identity", alpha=I(0.7)) +
  geom_area(data=fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$LINEAGE=="B.1.617.2",], aes(lwd=I(1.2), colour=NULL), fill=I("magenta"), position="identity", alpha=I(0.7)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  scale_y_log10() +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),max(cases_tot$date)))

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_smoothed_stacked area multinomial fit raw case data_log10 Y scale.png"), width=8, height=6)


ggplot(data=fit_cogukp2_multi_predsbyregion, 
       aes(x=collection_date, y=cases, group=LINEAGE)) + 
  facet_wrap(~ REGION, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE, group=LINEAGE)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),max(cases_tot$date)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") + 
  ylab("New confirmed cases per day") + xlab("Date of diagnosis") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit raw case data.png"), width=8, height=6)

ggplot(data=fit_cogukp2_multi_predsbyregion,
       aes(x=collection_date-7, y=smoothed_cases, group=LINEAGE)) +
  facet_wrap(~ REGION, scale="free") +
  geom_line(aes(lwd=I(1.2), colour=LINEAGE, group=LINEAGE)) +
  scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01"))),1,1),
                     limits=c(as.Date("2021-01-01"),today), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right") +
  ylab("New confirmed cases per day (smoothed)") + xlab("Date of infection") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT\nIN THE UK\n(case data gov.uk & multinomial fit to COG-UK data)") +
  scale_fill_manual("variant", values=lineage_cols1) +
  scale_colour_manual("variant", values=lineage_cols1) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today)) +
  scale_y_log10()

ggsave(file=paste0(".\\plots\\",plotdir,"\\cases per day_log10 y scale_multinomial fit smoothed case data.png"), width=8, height=6)






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
avg_r_cases = as.data.frame(emtrends(fit_cases, ~ DATE_NUM, by="region", var="DATE_NUM", 
                                     at=list(DATE_NUM=seq(date.from,
                                                          date.to, by=3)
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
qplot(data=avg_r_cases, x=DATE-7, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("steelblue")) +
  facet_wrap(~ region) +
  geom_line() + theme_hc() + xlab("Date of infection") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
  # scale_y_continuous(limits=c(1/2, 2), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re IN THE UK AT MOMENT OF INFECTION BASED ON NEW CASES\n(data gov.uk)") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# coord_cartesian(xlim=c(as.Date("2020-01-01"),NA))

# calculate above-average intrinsic growth rates per day of each variant over time based on multinomial fit using emtrends weighted effect contrasts ####
# for best model fit3_sanger_multi
above_avg_r_variants4 = do.call(rbind, lapply(levels_REGION, function(region) { do.call(rbind, 
                                                                                      lapply(seq(date.from,
                                                 date.to, by=3), 
                                             function (d) { 
                                               wt = as.data.frame(emmeans(fit4_cogukp2_multi, ~ LINEAGE , by="REGION", 
                                                                          at=list(DATE_NUM=d, REGION=region), type="response"))$prob   # important: these should sum to 1
                                               # wt = rep(1/length(levels_VARIANTS), length(levels_VARIANTS)) 
                                               # this would give equal weights, equivalent to emmeans:::eff.emmc(levs=levels_LINEAGE)
                                               cons = lapply(seq_along(wt), function (i) { con = -wt; con[i] = 1 + con[i]; con })
                                               names(cons) = seq_along(cons)
                                               EMT = emtrends(fit4_cogukp2_multi,  ~ LINEAGE , by=c("DATE_NUM", "REGION"),
                                                              var="DATE_NUM", mode="latent",
                                                              at=list(DATE_NUM=d, REGION=region))
                                               out = as.data.frame(confint(contrast(EMT, cons), adjust="none", df=NA))
                                               # sum(out$estimate*wt) # should sum to zero
                                               return(out) } )) } ))
  above_avg_r_variants = above_avg_r_variants4
  above_avg_r_variants$contrast = factor(above_avg_r_variants$contrast, 
                                         levels=1:length(levels_LINEAGE), 
                                         labels=levels(data_agbyweekregion1$LINEAGE))
  above_avg_r_variants$LINEAGE = above_avg_r_variants$contrast # gsub(" effect|\\(|\\)","",above_avg_r_variants$contrast)
  above_avg_r_variants$collection_date = as.Date(above_avg_r_variants$DATE_NUM, origin="1970-01-01")
  range(above_avg_r_variants$collection_date) # "2020-09-01" "2021-07-31"
  # average growth rate of all lineages calculated from case nrs
  above_avg_r_variants$avg_r = avg_r_cases$r[match(interaction(above_avg_r_variants$collection_date,above_avg_r_variants$REGION),
                                                   interaction(avg_r_cases$DATE,avg_r_cases$region))]  
  above_avg_r_variants$r = above_avg_r_variants$avg_r+above_avg_r_variants$estimate
  above_avg_r_variants$r_LOWER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.LCL
  above_avg_r_variants$r_UPPER = above_avg_r_variants$avg_r+above_avg_r_variants$asymp.UCL
  above_avg_r_variants$Re = Re.from.r(above_avg_r_variants$r)
  above_avg_r_variants$Re_LOWER = Re.from.r(above_avg_r_variants$r_LOWER)
  above_avg_r_variants$Re_UPPER = Re.from.r(above_avg_r_variants$r_UPPER)
  df = data.frame(contrast=NA,
                  DATE_NUM=avg_r_cases$DATE_NUM, # -7 to calculate back to time of infection
                  REGION=avg_r_cases$region,
                  estimate=NA,
                  SE=NA,
                  df=NA,
                  asymp.LCL=NA,
                  asymp.UCL=NA,
                  # p.value=NA,
                  collection_date=avg_r_cases$DATE,
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
  above_avg_r_variants$prob = fit_cogukp2_multi_predsbyregion$prob[match(interaction(round(above_avg_r_variants$DATE_NUM),
                                                                        as.character(above_avg_r_variants$LINEAGE),
                                                                        as.character(above_avg_r_variants$REGION)),
                                                            interaction(round(fit_cogukp2_multi_predsbyregion$DATE_NUM),
                                                                        as.character(fit_cogukp2_multi_predsbyregion$LINEAGE),
                                                                        as.character(fit_cogukp2_multi_predsbyregion$REGION)))]
  above_avg_r_variants2 = above_avg_r_variants
  ymax = 2
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
  qplot(data=above_avg_r_variants2[!((above_avg_r_variants2$LINEAGE %in% c("other"))|(above_avg_r_variants2$collection_date>=max(cases_tot$date))),], # |above_avg_r_variants2$collection_date>max(cases_tot$DATE)
        x=collection_date-7, # -7 to calculate back to date of infection
        y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", colour=LINEAGE, fill=LINEAGE, alpha=I(0.5),
        group=LINEAGE, linetype=I(0)) +
    facet_wrap(~ REGION) +
    # geom_ribbon(aes(fill=variant, colour=variant), alpha=I(0.5))
    geom_line(aes(colour=LINEAGE), lwd=I(0.72)) + theme_hc() + xlab("Date of infection") +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01","2021-07-01")),
                       labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J")) +
    # scale_y_continuous(limits=c(1/ymax,ymax), trans="log2") +
    geom_hline(yintercept=1, colour=I("red")) +
    ggtitle("Re VALUES OF SARS-CoV2 VARIANTS IN THE UK\nAT MOMENT OF INFECTION\n(based on gov.uk case data & multinomial fit to\nCOG-UK lineage frequencies)") +
    # labs(tag = tag) +
    # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
    theme(plot.tag.position = "bottomright",
          plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
    coord_cartesian(xlim=c(as.Date("2020-11-01"),max(cases_tot$date))) +
    scale_fill_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
    scale_colour_manual("variant", values=c(tail(lineage_cols1,-1),"black")) +
    theme(legend.position="right") 
  
ggsave(file=paste0(".\\plots\\",plotdir,"\\Re values per variant_avgRe_from_cases_with clipping.png"), width=8, height=6)


