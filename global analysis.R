# GLOBAL ANALYSIS OF SARS-Cov2 VARIANTS OF CONCERN & INTEREST USING MULTINOMIAL SPLINE FITS ####
# DATA: GISAID OR NCBI DATA, ACCESSED VIA COVSPECTRUM API ####

# T. Wenseleers
# last update 8 JANUARY 2023

# for similar analysis see https://nbviewer.org/github/gerstung-lab/SARS-CoV-2-International/blob/main/genomicsurveillance-int.ipynb#Check-some-fast-growing-lineages

# note: script below is fairly memory hungry - best to run this on workstation 
# with 64 Gb RAM - it runs quite fast though - just ca 30 mins including
# downloading all the lastest data from GISAID & COG-UK


# rm(list = ls()) # clear workspace
gc()

# set data source: open NCBI data or GISAID data ####

source = "GISAID" # or "NCBI", using GISAID data here


# set COVSPECTRUM credentials if you would like to use private GISAID data ####

# do this by setting REACT_APP_LAPIS_ACCESS_KEY as an environment variable 
# using Sys.setenv(REACT_APP_LAPIS_ACCESS_KEY = "XXXX") 
# cf https://cov-spectrum.org/static/js/main.b4ab11d7.js
# (using open NCBI data is also possible by setting source="NCBI" and then no key is needed)
if (file.exists("..//set_COVSPECTRUM_credentials.R")) source("..//set_COVSPECTRUM_credentials.R") 


# load (& if needed install) required packages & load some utility function ###
# install.packages("pacman")
library(pacman)
pacman::p_load(devtools, nnet, splines, mclogit, readr, ggplot2, ggthemes, scales,
               archive, dplyr, stringr, lubridate, tidyr, countrycode,
               memoise, readxl, covidregionaldata, tidyquant, data.table, R.utils,
               locatexec, pals, inspectdf, zoo, RSelenium, jsonlite)
Sys.setenv(GITHUB_PAT = "")
if (!require(marginaleffects)) devtools::install_github("tomwenseleers/marginaleffects")
require(marginaleffects)
pacman::p_load_gh("melff/mclogit/pkg", "rvlenth/emmeans",
                  "epiforecasts/covidregionaldata")

# load some utility functions to get aggregated NextcladePangolin lineage frequencies from covSPECTRUM
source(".//download_covSpectrum.R") # utility function download_covSpectrum to download data from covSpectrum


# select target variant (CovSpectrum nextcladePangoLineage notation) ####
target_variant = "XBB.1.5*"

# select countries to use in analyses, here those with >=5 XBB.1.5* sequences ####
# over the past 30 days & 100 sequenced genomes in total over that period
minseqs = 5
mintotalseqs = 100
lastdays = 30
toplist = countrieswithvariant(target_variant="XBB.1.5*", 
                                 minseqs=minseqs,
                                 mintotalseqs=mintotalseqs,
                                 lastdays=lastdays)
toplist
countries = unique(toplist$country)

# define countries for which we would like to divide by division (state/province) ####
countries_bydivision = c("USA", "India") 

# download variant data ####
date_from = "2019-12-01"
data_percountry = bind_rows(lapply(countries[!countries %in% countries_bydivision],
                         function (country) download_covSpectrum(source=source, 
                            date_from=date_from,
                            country=country, # vector of country names, NA for all
                            nextcladePangoLineages=NA, # get all lineages
                            bydate=TRUE, 
                            bypangolineage=TRUE,
                            bydivision=FALSE)
                         ))
data_percountry$division = data_percountry$country
data_percountry_perdivision = bind_rows(lapply(countries[countries %in% countries_bydivision],
                                   function (country) download_covSpectrum(source=source, 
                                                                           date_from=date_from,
                                                                           country=country, # vector of country names, NA for all
                                                                           nextcladePangoLineages=NA, # get all lineages
                                                                           bydate=TRUE, 
                                                                           bypangolineage=TRUE,
                                                                           bydivision=TRUE)
))
sort(unique(data_percountry_perdivision$division))
data_percountry_perdivision$division = gsub("Stockton", "California", data_percountry_perdivision$division, fixed=T)
data_percountry_perdivision$division = gsub("Washington DC", "Washington", data_percountry_perdivision$division, fixed=T)
data_percountry_perdivision = data_percountry_perdivision[!data_percountry_perdivision$division %in% 
                                                            c(NA, "Un", "USA", "Stockton", "American Samoa",
                                                              "Guam", "Northern Mariana Islands",
                                                              "Puerto Rico", "Virgin Islands"),]
sort(unique(data_percountry_perdivision$division))

data = bind_rows(data_percountry, data_percountry_perdivision)


# convert date to date format ####
data$date = as.Date(fast_strptime(data$date, "%Y-%m-%d")) # faster than as.Date(data$date)

# define output directories & info to put on plots if desired ####
today = as.Date(Sys.time())
today_num = as.numeric(today)
target_dir = "./data/GISAID"
suppressWarnings(dir.create(target_dir))
plotdir = "./plots/GISAID"
suppressWarnings(dir.create(plotdir))
tag = paste("@TWenseleers\n",today)


# define main variant lineages ####
# I need to code this in a generic way as in 
# https://nbviewer.org/github/gerstung-lab/SARS-CoV-2-International/blob/main/genomicsurveillance-int.ipynb#Check-some-fast-growing-lineages
# lineage_aliases = fromJSON("https://cov-spectrum.org/api/v2/resource/pango-lineage-alias")

# some regular expression helper functions
# lineage
lin = function (lineage, x=data$nextcladePangoLineage) { pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$")
grepl(pat, x, fixed=F, perl=T) }

# lineage plus sublineages
linplus = function (lineage, x=data$nextcladePangoLineage) { pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$","|",paste0("^", gsub(".","\\.",lineage, fixed=T)))
grepl(pat, x, fixed=F, perl=T) }

# one of X lineages
lin_oneof = function (lineages, x=data$nextcladePangoLineage) { 
  rowSums(sapply(lineages, function (lineage) { 
    pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$")
    grepl(pat, x, fixed=F, perl=T) }))>=1
}
# e.g. lin_oneof(c("BA.5","BF.7"))

# one of X lineages plus sublineages
linplus_oneof = function (lineages, x=data$nextcladePangoLineage) { 
  rowSums(sapply(lineages, function (lineage) { 
    pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$","|",paste0("^", gsub(".","\\.",lineage, fixed=T)))
    grepl(pat, x, fixed=F, perl=T) }))>=1
}
# linplus_oneof(c("BA.5","BF.7"))

# date from XX  
datefrom = function (date, x=data$date) { x >= as.Date(date) }
# datefrom(c("2021-01-01"))

# define lineages
levels_VARIANTS = c("Other", "B.1.177 (EU1)", "B.1.160 (EU2)", 
                    "B.1.221 (20A/S:98F)", "Beta", "Alpha", "Delta", 
                    "Omicron (BA.1*)", "Omicron (BA.2*)", 
                    "Omicron (BA.2.12.1*)", 
                    "Omicron (BA.2.38*)", 
                    "Omicron (BA.4*)", "Omicron (BA.4.6*)", 
                    "Omicron (BA.5*)", "Omicron (BA.5.2*)", "Omicron (BF.7*)", 
                    "Omicron (BA.2.76*)", "Omicron (BA.2.75*)", 
                    "Omicron (BQ.1*)", "Omicron (XAY.2*)",
                    "Omicron (CH.1.1*)", "Omicron (XBB*)", 
                    "Omicron (XBB.1.5*)"
)
target_variant = "Omicron (XBB.1.5*)"
baseline = "Omicron (BQ.1*)"
# I am using this order in plots, baseline is in fits coded as reference level

data$variant <- case_when(
  linplus("XBB.1.5") ~ "Omicron (XBB.1.5*)",
  linplus("XBB") ~ "Omicron (XBB*)",
  linplus("CH.1.1") ~ "Omicron (CH.1.1*)",
  linplus("XAY.2") ~ "Omicron (XAY.2*)",
  linplus("BQ.1")&
    datefrom("2022-02-01") ~ "Omicron (BQ.1*)",
  linplus("BA.2.75")&
    (!(linplus("BN")|linplus("BM")|linplus("BR")|linplus("BL")|linplus("CH")|linplus("BQ")))&
    datefrom("2022-04-01") ~ "Omicron (BA.2.75*)",
  linplus("BA.2.76") ~ "Omicron (BA.2.76*)",
  linplus("BF.7")&
    datefrom("2022-02-01") ~ "Omicron (BF.7*)",
  linplus("BA.5.2")&
    datefrom("2022-02-01") ~ "Omicron (BA.5.2*)",
  linplus("BA.5")&
    datefrom("2021-12-01") ~ "Omicron (BA.5*)",
  linplus("BA.4.6")&
    datefrom("2021-12-01") ~ "Omicron (BA.4.6*)",
  linplus("BA.4")&
    datefrom("2021-12-01") ~ "Omicron (BA.4*)",
  linplus("BA.2.12.1")&
    datefrom("2021-12-01") ~ "Omicron (BA.2.12.1*)",
  linplus("BA.2.38")&
    datefrom("2021-09-01") ~ "Omicron (BA.2.38*)",         
  linplus("BA.2")&
    datefrom("2021-09-01") ~ "Omicron (BA.2*)",
  linplus("BA.1")&
    datefrom("2021-09-01") ~ "Omicron (BA.1*)",
  linplus("B.1.617.2")|linplus("AY")&
    datefrom("2020-10-30") ~ "Delta",
  linplus("B.1.1.7")&
    datefrom("2020-09-20") ~ "Alpha",
  linplus("B.1.351")&
    datefrom("2020-08-10") ~ "Beta",
  linplus("B.1.221")&
    datefrom("2020-03-01") ~ "B.1.221 (20A/S:98F)",
  linplus("B.1.160")&
    datefrom("2020-02-15") ~ "B.1.160 (EU2)",
  linplus("B.1.177")&
    datefrom("2020-05-27") ~ "B.1.177 (EU1)",
  TRUE ~ "Other"
  )

length(unique(data$variant)) == length(levels_VARIANTS) # correct

# note: in India BA.2.38 & BA.2.38.1 caused an infection wave in some states - hence separated out above
# B.1.177+B.1.160+B.1.221 were behind the 2020 wave in fall in Europe & each had one spike mutations & a small growth rate advantage relative to predominant B.1.1

# TO DO: maybe also add XBF, XBK and BR.2.1 - the main contenders in Europe besides CH.1.1 and XBB.1.5*

# earliest realistic dates were taken from
# https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/clade_emergence_dates.tsv

# define lineage colours
n = length(levels_VARIANTS)
lineage_cols = case_when(
    levels_VARIANTS=="Other" ~ "grey65",
    levels_VARIANTS=="B.1.177 (EU1)" ~ "darkorange4",
    levels_VARIANTS=="B.1.160 (EU2)" ~ "darkorange3",
    levels_VARIANTS=="B.1.221 (20A/S:98F)" ~ "darkorange2",
    levels_VARIANTS=="Beta" ~ "green4",
    levels_VARIANTS=="Alpha" ~ "#0085FF",
    levels_VARIANTS=="Delta" ~ "mediumorchid",
    levels_VARIANTS=="Omicron (BA.1*)" ~ "red",
    levels_VARIANTS=="Omicron (BA.2*)" ~ "red3",
    levels_VARIANTS=="Omicron (BA.2.12.1*)" ~ "black",
    levels_VARIANTS=="Omicron (BA.2.38*)" ~ "red4",
    levels_VARIANTS=="Omicron (BA.4*)" ~ "green3",
    levels_VARIANTS=="Omicron (BA.4.6*)" ~ "green2",
    levels_VARIANTS=="Omicron (BA.5*)" ~ "blue4",
    levels_VARIANTS=="Omicron (BA.5.2*)" ~ "blue3",
    levels_VARIANTS=="Omicron (BF.7*)" ~ "dodgerblue",
    levels_VARIANTS=="Omicron (BA.2.76*)" ~ "magenta4",
    levels_VARIANTS=="Omicron (BA.2.75*)" ~ "magenta",
    levels_VARIANTS=="Omicron (BQ.1*)" ~ "orange",
    levels_VARIANTS=="Omicron (XAY.2*)" ~ "yellow2",
    levels_VARIANTS=="Omicron (CH.1.1*)" ~ "yellow3",  
    levels_VARIANTS=="Omicron (XBB*)" ~ "cyan3",
    levels_VARIANTS=="Omicron (XBB.1.5*)" ~ "cyan"
  )
names(lineage_cols) = levels_VARIANTS
length(levels_VARIANTS) == length(lineage_cols) # correct
pal.bands(lineage_cols)

data$variant = factor(data$variant, levels=levels_VARIANTS)

table(data$variant)
table(data$region, data$variant)


# add week, year & start of week ####
# TO DO: add these columns in download functions?

data$Week = lubridate::week(data$date)
data$Year = lubridate::year(data$date)
data$Year_Week = interaction(data$Year,data$Week)
data$floor_date = fast_strptime(as.character(cut(data$date, "week")), "%Y-%m-%d") # start of week
data$date_num = as.numeric(data$date)

# sort countries by incidence of target variant
data$country = factor(data$country, levels=toplist$country)

# sort divisions/states & countries by incidence of target variant
total_count_target_variant = data %>%
  filter(date>=(as.Date(Sys.time())-lastdays)) %>%
  group_by(variant, division) %>%
  summarise(total_count_target_variant = sum(count)) %>%
  filter(variant==target_variant)
total_n = data %>%
  filter(date>=(as.Date(Sys.time())-lastdays)) %>%
  group_by(division) %>%
  summarise(total_n = sum(count))
total_count_target_variant = left_join(total_count_target_variant, total_n) %>%
  mutate(prop_target_variant = total_count_target_variant / total_n) %>%
  arrange(desc(prop_target_variant))

data$division = factor(data$division, levels=c(total_count_target_variant$division, 
                                               sort(unique(data$division[!data$division %in% total_count_target_variant$division]))))
table(data$division, data$variant)

levels_regions = c("Asia","North America","Europe","Africa","South America","Oceania")
data$region = factor(data$region, levels=levels_regions)
data$region = droplevels(data$region)



# 2. COMBINED GLOBAL ANALYSIS USING MULTINOMIAL SPLINE FITS ####

# AGGREGATED DATA BY WEEK ####
data_agbyweek1 = data %>%
  group_by(variant, floor_date) %>%
  summarise(count = sum(count))
data_agbyweek1_total = data %>%
  group_by(floor_date) %>%
  summarise(total = sum(count))
data_agbyweek1 = left_join(data_agbyweek1, data_agbyweek1_total) %>%
  rename(date = floor_date) %>%
  mutate(prop = count/total)
data_agbyweek1 = as.data.frame(data_agbyweek1)
data_agbyweek1$date = as.Date(data_agbyweek1$date)
data_agbyweek1$date_num = as.numeric(data_agbyweek1$date)
# write.csv(data_agbyweek1, file=".//data//COVSPECTRUM//COVSPECTRUM aggregated counts by week_all.csv", row.names=F)


# AGGREGATED DATA BY WEEK & COUNTRY/DIVISION ####
data_agbyweekcountry1 = data %>%
  group_by(division, variant, floor_date) %>%
  summarise(count = sum(count))
data_agbyweekcountry1_total = data %>%
  group_by(division, floor_date) %>%
  summarise(total = sum(count))
data_agbyweekcountry1 = left_join(data_agbyweekcountry1, data_agbyweekcountry1_total) %>%
  rename(date = floor_date) %>%
  mutate(prop = count/total)
data_agbyweekcountry1 = as.data.frame(data_agbyweekcountry1)
data_agbyweekcountry1$date = as.Date(data_agbyweekcountry1$date)
data_agbyweekcountry1$date_num = as.numeric(data_agbyweekcountry1$date)
data_agbyweekcountry1$region = data$region[match(data_agbyweekcountry1$division, data$division)]
# write.csv(data_agbyweekcountry1, file="./data/COVSPECTRUM/COVSPECTRUM aggregated counts by start of week and lineage_all.csv", row.names=F)

gc()



# MULLER PLOT (RAW DATA, selected countries pooled, but with big sampling biases across countries/divisions) ####
data_agbyweek1$variant = factor(data_agbyweek1$variant, levels=levels_VARIANTS)
muller_raw_all = ggplot(data=data_agbyweek1, aes(x=date, y=count, group=variant)) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES",
       subtitle=paste0("Raw GISAID data up to ",today," (using NextClade Pangolin lineages)")) +
  coord_cartesian(xlim=c(min(data$date),max(data$date)), expand=c(0)) +
  labs(tag = tag) + theme(plot.tag.position = "bottomright", plot.tag = element_text(vjust = 1, hjust = 1, size=8))
muller_raw_all

ggsave(file=file.path(plotdir,"muller plot_raw data.png"), width=12, height=5)




# FIT NNET::MULTINOM MULTINOMIAL SPLINE MODEL ####

# code uses baseline lineage as reference level
data_agbyweekcountry1$variant = relevel(data_agbyweekcountry1$variant, ref=baseline)

set.seed(1)
gc()
system.time(fit <- nnet::multinom(variant ~ ns(date_num, df=2)+ns(date_num, df=2):region+division, 
                                       weights=count, 
                                       data=data_agbyweekcountry1, 
                                       maxit=10000, MaxNWts=100000)) # iter 2680, 649s
# syntax of model to put on plot legend
model = "variant ~ ns(date, df=2)+ns(date, df=2):region+division" 

# TO DO: change to mclogit::mblogit fit (can take into account overdispersion &
# latest github version should run - previously it was giving fitting errors) or
# the MGLM package - but that one also gave fitting errors

# model to use below - I just fitted 1 possible model now - one could vary df etc
fit_best = fit

# we calculate the Hessian using my own faster Rcpp Kronecker-product based function
source(".//fastmultinomHess.R") # faster way to calculation Hessian of multinomial fits
gc()
system.time(fit_best$Hessian <- fastmultinomHess(fit_best, model.matrix(fit_best))) # 81s
# we add variance-covariance matrix as extra slot to be re-used later
system.time(fit_best$vcov <- vcov(fit_best))  # 6s
gc()

save.image("~/Github/LineageExplorer/environment_2023_9_1.RData")x
# load("~/Github/LineageExplorer/environment_2023_9_1.RData")

# save multinom fit
saveRDS(fit_best, file="~/Github/LineageExplorer/LineageExplorer/fits/multinom_fit.rds")


# CALCULATE GROWTH RATE ADVANTAGE OVER BASELINE REFERENCE LEVEL BQ.1* TODAY ####

## TO DO still a bug here with groth rate advantages coming out negative - need to fix this

# with new faster marginaleffects code
library(marginaleffects)
system.time(meffects <- marginaleffects(fit_best, 
                                               type = "link", # = additive log-ratio = growth rate advantage relative to BQ.1*
                                               variables = c("date_num"),
                                               by = c("group"),
                                               vcov = fit_best$vcov,
                                               newdata = datagrid(date_num = today_num 
                                               ))) # 30s

# growth rate advantage compared to reference level BA.5.2 by region
system.time(meffects_byregion <- marginaleffects(fit_best, 
                                               type = "link", # = additive log-ratio = growth rate advantage relative to BQ.1*
                                               variables = c("date_num"),
                                               by = c("group", "region"),
                                               vcov = fit_best$vcov,
                                               newdata = datagrid(date_num = today_num,
                                                                  region = unique(data_agbyweekcountry1$region)
                                               ))) # 30s

# for all pairwise growth rate differences:
# growth_differences = comparisons(
#   fit_best,
#   newdata = datagrid(date_num = today_num),
#   variables = "date_num",
#   by = "region",
#   type = "clr", # here we could either use "clr" (centered logratio) or "link" (additive logratio) - this gives same result
#   hypothesis = "pairwise")


# old emtrends code to calculate pairwise growth rate differences
# system.time(emtr_pairw <- emtrends(fit_best, revpairwise ~ variant,
#                                    by="region",
#                                    var="date_num",  mode="latent",
#                                    at=list(date_num=today_num))) # 
# delta_r_pairw = data.frame(confint(emtr_pairw,
#                                    adjust="none", df=NA)$contrasts,
#                            p.value=as.data.frame(emtr_pairw$contrasts)$p.value)
# delta_r_pairw
# write.csv(delta_r_pairw, file.path(plotdir, "growth rate advantage all variants vs BA_5_2.csv"), row.names=F)


# plot of growth rate advantage of last n newest variants
# TO DO: order by selective advantage and then take top lastn
lastn = 11 # last n variants to show - change in top n ?
sel_variants = tail(levels_VARIANTS,lastn)
sel_variants = sel_variants[!sel_variants %in% c(baseline, "Omicron (BA.4.6*)", "Omicron (BA.5*)")]
meffects_sel1 = meffects[meffects$group %in% sel_variants,]
meffects_sel1$group = factor(meffects_sel1$group, levels=meffects_sel1$group[order(meffects_sel1$dydx, decreasing=T)])
cols = colorRampPalette(c("red3", "blue3"))(length(levels(meffects_sel1$group)))
qplot(data=meffects_sel1, 
      x=group, y=dydx*100, ymin=conf.low*100, ymax=conf.high*100, fill=group, geom="col", 
      width=I(0.7)) +
  geom_linerange(aes(lwd=I(0.4))) + ylab("Growth rate advantage\nrelative to BQ.1* (% per day)") +
  scale_fill_manual(values=cols) +
  theme(legend.position="none") + xlab("") +
  ggtitle("GROWTH RATE ADVANTAGE OF SARS-CoV2 VARIANTS",
          subtitle=paste0("based on multinomial fit ", model, "\nGISAID data with NextcladePangolin lineage definition,\nusing data from countries with >=", minseqs, " XBB.1.5* sequences") ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
  
ggsave(file=file.path(plotdir,"growth rate advantage VOCs_overall.png"), width=7, height=5)

# plot of growth rate advantage of last X newest variants by continent
# TO DO: order by selective advantage and then take top n
sel_variants = tail(levels_VARIANTS,lastn)
sel_variants = sel_variants[!sel_variants %in% c(baseline, "Omicron (BA.4.6*)", "Omicron (BA.5*)")]
sel_regions = unique(data_agbyweekcountry1$region)
# sel_regions = sel_regions[!sel_regions %in% c("Africa")] # too little data
meffects_sel2 = meffects_byregion[meffects_byregion$region %in% sel_regions,]
meffects_sel2 = meffects_sel2[meffects_sel2$group %in% sel_variants,]
meffects_sel2$group = factor(meffects_sel2$group, levels=levels(meffects_sel1$group))
# outlier = (abs(meffects_sel2$dydx)>=0.35)|(meffects_sel2$dydx<0) # typically due to there being too little data
# meffects_sel2 = meffects_sel2[!outlier,]

tbl = as.data.frame(table(data[data$variant %in% sel_variants,"region"], 
                          data[data$variant %in% sel_variants,"variant"]))
colnames(tbl) = c("region", "variant", "count")
meffects_sel2$count = tbl$count[match(interaction(meffects_sel2$region, meffects_sel2$group),
                                      interaction(tbl$region, tbl$variant))]
# retain only estimates with total count > minvariantseqs per region/continet
minvariantseqs = 50
meffects_sel2 = meffects_sel2[meffects_sel2$count>=minvariantseqs, ]

qplot(data=meffects_sel2, 
      x=group, y=dydx*100, ymin=conf.low*100, ymax=conf.high*100, fill=group, geom="col", 
      width=I(0.7)) +
  facet_wrap(~ region) +
  geom_linerange(aes(lwd=I(0.4))) + ylab("Growth rate advantage\nrelative to BQ.1* (% per day)") +
  scale_fill_manual(values=cols) +
  theme(legend.position="none") + xlab("") +
  ggtitle("GROWTH RATE ADVANTAGE OF SARS-CoV2 VARIANTS",
          subtitle=paste0("based on multinomial fit ", model, "\nGISAID data with NextcladePangolin lineage definition,\nusing data from countries with >=", minseqs, " XBB.1.5* sequences", 
                          "\nEstimates shown for regions with >=",minvariantseqs," seqs of each variant") ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(file=file.path(plotdir,"growth rate advantage VOCs_by continent.png"), width=7, height=5)

reg = "North America"
qplot(data=meffects_sel2[meffects_sel2$region==reg,], 
      x=group, y=dydx*100, ymin=conf.low*100, ymax=conf.high*100, fill=group, geom="col", 
      width=I(0.7)) +
  geom_linerange(aes(lwd=I(0.4))) + ylab("Growth rate advantage\nrelative to BQ.1* (% per day)") +
  scale_fill_manual(values=cols) +
  theme(legend.position="none") + xlab("") +
  ggtitle(paste0("GROWTH RATE ADVANTAGE OF SARS-CoV2 VARIANTS\nIN ", toupper(reg)),
          subtitle=paste0("based on multinomial fit ", model, "\nGISAID data") 
          ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(file=file.path(plotdir,"growth rate advantage VOCs_NAmerica.png"), width=7, height=5)



# PLOT MULTINOMIAL FIT ####

extrapolate = 14
date.from = as.numeric(as.Date("2020-01-01"))
date.to = today_num+extrapolate

# multinomial model predictions by country/divisions with CIs calculated using margineffects::predictions

step=3
predgrid = expand.grid(list(date_num=as.numeric(seq(date.from, date.to, by=step)),
                            division=unique(data_agbyweekcountry1$division)))
predgrid$region = data_agbyweekcountry1$region[match(predgrid$division,
                                                     data_agbyweekcountry1$division)]


# note: now using Delta method on response scale, better to
# calculate CIs as in Effects package on link scale (type="link") or
# on isometric logratio scale (type="ilr") & then backtransform
# but still having some problems with over/underflows with
# type="link" and type="ilr"=isometric logratio is a bit more hassle to backtransform

# rm(fit_preds)
gc()
system.time(fit_preds <- data.frame(predictions(fit_best, 
                       newdata = predgrid,
                       type = "probs",
                       vcov = fit_best$vcov))) # %>% # 810s
             # transform(conf.low = predicted - 1.96 * std.error,
             #          conf.high = predicted + 1.96 * std.error) %>%
             # group_by(rowid) |>
             #mutate_at(c("predicted", "conf.low", "conf.high"), function (x) plogis(x)))
range(fit_preds$predicted)
hist(fit_preds$predicted) # comes out as either 0 or 1 - fix bug
gc()
fit_preds$conf.high[fit_preds$conf.high>0.99999] = 0.99999 # slight artefact of Delta method on response scale
# fit_preds$conf.high[fit_preds$conf.high<1E-10] = 1E-10
fit_preds$conf.low[fit_preds$conf.low<1E-10] = 1E-10
# fit_preds$conf.low[fit_preds$conf.low>0.99999-10] = 0.99999
# fit_preds$predicted[fit_preds$predicted>0.99999] = 0.99999
# fit_preds$predicted[fit_preds$predicted<1E-10] = 1E-10

# replace NAs by 0
# fit_preds <- fit_preds %>% mutate(predicted = ifelse(is.na(predicted), 0, predicted),
#                                   conf.low = ifelse(is.na(conf.low), 0, conf.low),
#                                   conf.high = ifelse(is.na(conf.high), 0, conf.high))
fit_preds$date = as.Date(fit_preds$date_num, origin="1970-01-01")
fit_preds$variant = NULL
colnames(fit_preds)[which(colnames(fit_preds)=="group")] = "variant"
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS)

# TO DO: fix bug with type="link" where some predictions come out as NA,
# and/or switch to type="ilr" - check in my marginaleffects fork

write_csv(fit_preds, file=file.path(plotdir, "COVSPECTRUM fitted lineage frequencies global multinomial fit.csv"))
# saveRDS(fit_preds, file.path(plotdir, "COVSPECTRUM fitted lineage frequencies global multinomial fit.rds"))



# PLOT MULTINOMIAL FIT ON LOGIT SCALE ####

ncls = round(sqrt(length(unique(data$division))))
pl = qplot(data=fit_preds[fit_preds$variant!="Other",], 
                         x=date, y=predicted, geom="blank") +
  facet_wrap(~ division, ncol=ncls) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, ", multinomial fit ", model, ",\nall countries with >=", minseqs, " XBB.1.5* sequences included, using NextCladePangolin lineages")) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +   
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other",],
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.1, 3), limits=c(1,max(data_agbyweekcountry1$total)), 
                        breaks=c(10,100,1000, 10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)") +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  labs(tag = tag) + theme(plot.tag.position = "bottomright", plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
  # theme(strip.text.x = element_text(size = 20))
pl

ggsave(file=file.path(plotdir, "predicted lineage freqs_logit scale.png"), 
       width=2*ncls, 
       height=(2/1.5)*ncls)


# zoomed in on last 6 months
pl = qplot(data=fit_preds[fit_preds$variant!="Other",], 
                         x=date, y=predicted, geom="blank") +
  facet_wrap(~ division) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, ", multinomial fit ", model, ",\nall countries with >=", minseqs, " XBB.1.5* sequences included, using NextCladePangolin lineages")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand=FALSE) +   
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other",],
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.1, 3), limits=c(1,max(data_agbyweekcountry1$total)), 
                        breaks=c(10,100,1000, 10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)") +
  coord_cartesian(xlim=c(today-30*6,NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
pl

ggsave(file=file.path(plotdir, "predicted lineage freqs_last6months_logit scale.png"), 
       width=2*ncls, 
       height=(2/1.5)*ncls)


# plot for selected divisions
divisions = c("New York","Connecticut","New Hampshire","Massachusetts","Rhode Island","New Jersey","Pennsylvania",
              "Israel","Netherlands","United Kingdom","Canada","Denmark","Germany","France")

ncls = round(length(divisions)/2) # round(sqrt(length(divisions)))
pl = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$division %in% divisions,], 
           x=date, y=predicted, geom="blank") +
  facet_wrap(~ division, ncol=ncls) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, "\nusing NextCladePangolin lineages\nmultinomial fit ", model)) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +   
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$division %in% divisions,],
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.1, 3), limits=c(1,max(data_agbyweekcountry1$total)), 
                        breaks=c(10,100,1000, 10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)") +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  labs(tag = tag) + theme(plot.tag.position = "bottomright", plot.tag = element_text(vjust = 1, hjust = 1, size=8))
pl

ggsave(file=file.path(plotdir, "predicted lineage freqs_logit scale_selected states countries.png"), 
       width=2*ncls, 
       height=0.75*ncls)

# zoomed in on last 6 months
pl = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$division %in% divisions,], 
           x=date, y=predicted, geom="blank") +
  facet_wrap(~ division, ncol=ncls) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, "\nusing NextCladePangolin lineages\nmultinomial fit ", model)) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand=FALSE) +   
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$division %in% divisions,],
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.1, 3), limits=c(1,max(data_agbyweekcountry1$total)), 
                        breaks=c(10,100,1000, 10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)") +
  coord_cartesian(xlim=c(today-30*6,NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
pl

ggsave(file=file.path(plotdir, "predicted lineage freqs_last6months_logit scale_selected states countries.png"), 
       width=3*ncls, 
       height=1*ncls)



for (division in divisions) {

# plot just for given division
pl = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$division==division,], 
           x=date, y=predicted, geom="blank") +
  # facet_wrap(~ country) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed\ninfections (%)") +
  theme_hc() + xlab("") +
  ggtitle(paste0("SARS-CoV2 LINEAGE FREQUENCIES IN ", toupper(division)),
          subtitle=paste0("GISAID data up to ",today, ", multinomial fit ", model, ", using NextCladePangolin lineages")) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +   
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&
                                          data_agbyweekcountry1$division==division,],
             aes(x=date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.5, 4), limits=c(1,max(data_agbyweekcountry1$total)), 
                        breaks=c(10,100,1000, 10000), guide=F) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "bottom") +
  xlab("Collection date (start of week)") +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# theme(plot.title=element_text(size=25)) +
# theme(plot.subtitle=element_text(size=20))
pl

ggsave(file=file.path(plotdir, paste0("predicted lineage freqs_",division,"_logit scale.png")), width=12, height=5)

}




# TO DO: DIDN'T UPDATE PART BELOW YET
# (map share of variants on cases & hospitalisations etc - should get US data by state)


# plot predicted values as Muller plot
pl = ggplot(data=fit_preds[fit_preds$date>=as.Date("2021-01-01"),], 
                    aes(x=date, y=predicted, group=variant)) +
  facet_wrap(~ division, ncol=ncls) +
  geom_area(aes(width=I(10), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +   
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES",
       subtitle=paste0("GISAID data up to ",today, ", multinomial fit ", model, ",\nall countries with >=", minseqs, " XBB.1.5* sequences shown")) +
  annotate("rect", 
           xmin=max(data$date)+1, 
           xmax=as.Date(date.to, origin="1970-01-01")+5, 
           ymin=0, ymax=1, alpha=0.5, fill="white") + # extrapolated part
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
pl

ggsave(file=file.path(plotdir, "muller plot.png"), 
       width=4*ncls, 
       height=(4/1.2)*ncls)


# MAP VARIANT SHARE ONTO CASE NUMBERS ####

country_data = get_national_data(countries=sel_countries,
                                 source="who",
                                 level=1)
# sel_countries_fig = c("Austria", "Belgium", "Denmark", "France", "Germany", 
#                       "Italy", "Netherlands", "New Zealand", "Singapore", "United Kingdom")
# country_data = get_national_data(countries=sel_countries_fig,
#                                  source="who",
#                                  level=1)
country_data$country[country_data$country=="United States"] = "USA"
country_data$country = factor(country_data$country, levels=levels(fit_preds$country))
country_data = country_data[country_data$date<(max(country_data$date)-3),] # drop data last 3 days
qplot(data=country_data, x=date, y=cases_new, group=country, geom="blank", colour=country) +
  facet_wrap(~country, scale="free_y") +
  geom_ma(ma_fun = SMA, n = 7, lty=1) +
  # scale_colour_hue("") +
  ylab("New cases per day") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  coord_trans(y="sqrt") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  theme(legend.position="none") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY",
          subtitle="7 day simple moving average, using\nWHO case data accessed via covidregionaldata package")
# ggsave(file=file.path(plotdir, "new cases countries with more than 10 sequenced BA_2_75 cases.png"), width=7, height=5)

fit_preds$totnewcases = 
  country_data$cases_new[match(interaction(fit_preds$country,
                                          fit_preds$date),
                                      interaction(country_data$country, 
                                                  country_data$date))]
fit_preds = fit_preds %>% 
  group_by(country) %>% 
  mutate(totnewcases_smoothed = rollmean(totnewcases, 7, na.pad = T))
fit_preds$cases = fit_preds$totnewcases_smoothed*fit_preds$predicted
fit_preds$cases[fit_preds$cases==0] = NA
fit_preds$cases[fit_preds$predicted<0.001] = NA

fit_preds_sel = fit_preds
# fit_preds_sel = fit_preds_sel[fit_preds_sel$country %in% sel_countries_fig,]
fit_preds_sel = fit_preds_sel[fit_preds_sel$country %in% sel_countries,]
fit_preds_sel = fit_preds_sel[!fit_preds_sel$country %in%
                                c("Hong Kong","Czech Republic","Reunion"),]
fit_preds_sel$country = factor(fit_preds_sel$country)

fit_preds2 = fit_preds_sel
fit_preds2$cases[fit_preds2$cases==0] = NA
fit_preds2$cases[fit_preds2$cases<=1] = NA
fit_preds2$country = factor(fit_preds2$country)
fit_preds2$variant = factor(fit_preds2$variant, levels=levels_VARIANTS)

ggplot(data=fit_preds2, 
       aes(x=date, y=cases)) + 
  facet_wrap(~ country, scale="free", ncol=ncls) +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_preds2, aes(x=date, y=totnewcases_smoothed, lwd=I(1.5)), 
            colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial fit ", model, ",\nselected countries with >=", min_n, " level5 or level6+ variant sequences shown")) +
  scale_colour_manual("lineage", values=lineage_cols) +
  scale_y_log10() +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) 
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$date)-20))

ggsave(file=file.path(plotdir,"confirmed cases multinomial fit_by country_log10 scale.png"), 
       width=4*ncls, height=4*ncls/1.2)

ggplot(data=fit_preds2, 
       aes(x=date, y=cases)) + 
  facet_wrap(~ country, scale="free", ncol=ncls) +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_preds2, aes(x=date, y=totnewcases_smoothed, lwd=I(1.5)), 
            colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial fit ", model, ",\nselected countries with >=", min_n, " level5 or level6+ variant sequences shown")) +
  scale_colour_manual("lineage", values=lineage_cols) +
  scale_y_log10() +
  # coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2021-01-01"),NA)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(xlim=c(as.Date("2022-01-01"),NA))

ggsave(file=file.path(plotdir,"confirmed cases multinomial fit_by country_log10 scale_zoomed.png"), 
       width=4*ncls, height=4*ncls/1.2)

# stacked area chart
ggplot(data=fit_preds2, 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y", ncol=ncls) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial fit ", model, ",\nselected countries with >=", min_n, " level5 or level6+ variant sequences shown")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot by country.png"), 
       width=4*ncls, height=4*ncls/1.2)


# stacked area chart, some selected countries
ggplot(data=fit_preds2[fit_preds2$country %in% 
                         c("France", "Belgium", "Denmark", "Germany"),], 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y", ncol=2) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial fit ", model, ",\na few selected countries shown")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot by country_sel countries.png"), 
       width=6*2, height=6*2/2)


# stacked area chart, just for Belgium
ggplot(data=fit_preds2[fit_preds2$country=="Belgium",], 
       aes(x=date, y=cases, group=variant)) + 
  # facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN BELGIUM",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, "\nand multinomial fit ", model, ",\nHere only showing data for Belgium")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot_BELGIUM.png"), width=12, height=6)

# stacked area chart, just for France
ggplot(data=fit_preds2[fit_preds2$country=="France",], 
       aes(x=date, y=cases, group=variant)) + 
  # facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN FRANCE",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, "\nand multinomial fit ", model, ",\nHere only showing data for France")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot_FRANCE.png"), width=12, height=6)


# stacked area chart, just for Switzerland
ggplot(data=fit_preds2[fit_preds2$country=="Switzerland",], 
       aes(x=date, y=cases, group=variant)) + 
  # facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN SWITZERLAND",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, "\nand multinomial fit ", model, ",\nHere only showing data for Switzerland")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot_SWITZERLAND.png"), width=12, height=6)


# stacked area chart, just for United Kingdom
ggplot(data=fit_preds2[fit_preds2$country=="United Kingdom",], 
       aes(x=date, y=cases, group=variant)) + 
  # facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN THE UK",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " & COG-UK data\nand multinomial fit ", model, ",\nHere only showing data for the UK")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot_UK.png"), width=12, height=6)


# stacked area chart, just for Singapore
ggplot(data=fit_preds2[fit_preds2$country=="Singapore",], 
       aes(x=date, y=cases, group=variant)) + 
  # facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN SINGAPORE",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, "\nand multinomial fit ", model, ",\nHere only showing data for Singapore")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot_SINGAPORE.png"), width=12, height=6)

# stacked area chart, just for Israel
ggplot(data=fit_preds2[fit_preds2$country=="Israel",],
       aes(x=date, y=cases, group=variant)) +
  # facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant),
            position="stack") +
  scale_fill_manual("", values=lineage_cols) +
  # annotate("rect", xmin=max(GISAID_india$date_num)+1,
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN ISRAEL",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, "\nand multinomial fit ", model, ",\nHere only showing data for Israel")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2020-01-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot_ISRAEL.png"), width=12, height=6)



# zoomed in since june 2022
ggplot(data=fit_preds2[fit_preds2$date>=as.Date("2022-06-01"),], 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y", ncol=ncls) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols[-c(2,3,4,5,6,7,8)]) + # TO DO: determine dropped levels automatically
  # annotate("rect", xmin=max(GISAID_india$date_num)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial fit ", model, ",\nselected countries with >=", min_n, " level5 or level6+ variant sequences shown")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2022-06-01"),NA))

ggsave(file=file.path(plotdir,"\\confirmed cases multinomial fit_stacked area plot by country_zoomed.png"), 
       width=4*ncls, height=4*ncls/1.2)

# TO DO : get IHME infection estimates from links below,
# convolve these to cases & map them onto variant frequencies
# https://www.healthdata.org/covid/data-downloads
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2020.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2021.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2022.csv

# TO DO : get mortality data from Eurostat with Eurostat package & https://www.mortality.org/Data/STMF
# & for the US CDC Wonder data or usmortality.com
# convolve case data by variant to mortality data & calculate
# death toll of each variant

# TO DO : finish spatial multinomial tensor spline fit in function of latitude & longitude & time
# (especially good for countries with limited or no sequencing data)


# FUNCTION TO CALCULATE EFFECTIVE REPRODUCTION NUMBER R FROM
# THE INSTANTANEOUS GROWTH RATE r, ASSUMING A GAMMA DISTRIBUTED GENERATION TIME
# from epiforecasts package growth_to_R
# https://github.com/epiforecasts/EpiNow2/blob/5015e75f7048c2580b2ebe83e46d63124d014861/R/utilities.R#L109
# https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
# (better if distribution is known to be gamma, based on the full integral)
Re_from_r <- function(r, gamma_mean=4.7, gamma_sd=2.9) {
  k <- (gamma_sd / gamma_mean)^2
  R <- (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

Re_from_r(0.18, gamma_mean=4.7, gamma_sd=2.9) # from fit Moritz Gerstung for Germany
# 2.1/6=35% of the population susceptible to BQ.1.1
gc()


