# GLOBAL ANALYSIS OF SARS-Cov2 VARIANTS OF CONCERN & INTEREST
# T. Wenseleers
# last update 6 OCTOBER 2022

# set GISAID credentials ####
# set them first using 
# Sys.setenv(GISAIDR_USERNAME = "XXXXX")
# Sys.setenv(GISAIDR_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_USERNAME = "XXXXX") # not needed for this script
# Sys.setenv(GISAIDJSON_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_STREAM = "XXXXX")
if (file.exists("..//set_GISAID_credentials.R")) source("..//set_GISAID_credentials.R") # using my GISAID credentials


# load some utility functions & install required packages 
library(devtools)
remotes::install_github("Wytamma/GISAIDR")
library(GISAIDR)
source(".//download_GISAID.R") # load function to download GISAID metadata download package (lacking records from last few days)
source(".//download_GISAID_records.R") # load functions to download most recent GISAID records
source(".//download_COGUK.R") # load function to download COG-UK metadata
use_coguk = TRUE # use COG-UK data instead of GISAID data for UK?

# load required packages
library(nnet)
library(splines)
library(devtools)
remotes::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
remotes::install_github("rvlenth/emmeans", dependencies = FALSE)
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)
library(archive)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(countrycode)
library(memoise)
library(readxl)
# install.packages("covidregionaldata",
#                   repos = "https://epiforecasts.r-universe.dev")
remotes::install_github("epiforecasts/covidregionaldata")
library(covidregionaldata)
library(tidyquant)
library(data.table)
library(R.utils)
library(locatexec)
library(tidyr)
library(pals)
library(devtools)
install_github("tomwenseleers/marginaleffects")
library(marginaleffects)
library(inspectdf)
library(zoo)
library(RSelenium)


# 1. LOAD DATA ####

today = as.Date(Sys.time())
today_num = as.numeric(today)
today # "2021-09-25"
plotdir = "GISAID_COGUK_GLOBAL_ANALYSIS"
target_dir = "C:/Users/bherr/OneDrive - KU Leuven/Documents/Github/LineageExplorer/LineageExplorer/data/GISAID" # target download directory GISAID data
suppressWarnings(dir.create(file.path("plots",plotdir)))
tag = paste("@TWenseleers\n",today)

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2023-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

# import GISAID metadata ####

# download latest GISAID metadata ####

# note: make sure to have a working installation of RSelenium & chrome
# browser installed
system.time(GISAID <- download_GISAD_meta(target_dir = target_dir,
                                          headless = FALSE,
                                          usr = Sys.getenv("GISAIDR_USERNAME"),
                                          psw = Sys.getenv("GISAIDR_PASSWORD"))) # 194s
download = tail(list.files(target_dir, pattern=".tar.xz"), 1)
download # "metadata_tsv_2022_10_03.tar.xz"  

# records go up to submission date
GISAID_max_submdate = as.Date(max(GISAID$submission_date, na.rm=T))
GISAID_max_submdate # "2022-10-01"
nrow(GISAID) # 13341604


# add some extra recently submitted records not available in download package ####

# I'll use https://github.com/Wytamma/GISAIDR
# but in current implementation field AA substitutions still missing, see 
# https://github.com/Wytamma/GISAIDR/issues/27
# https://github.com/Wytamma/GISAIDR/issues/26
# so added extra function in download_GSIAD_records.R that fixes this
credentials = login(username = Sys.getenv("GISAIDR_USERNAME"),
                    password = Sys.getenv("GISAIDR_PASSWORD"),
                    database = "EpiCoV")
# dataframe with recently uploaded records
# (not included in metadata download package)
recent_records = as.vector(query(
  credentials = credentials, 
  from_subm = as.character(GISAID_max_submdate), 
  to_subm = as.character(today),
  fast = TRUE
))$accession_id
length(recent_records) # 49648   records
recent_records = recent_records[!recent_records %in% GISAID$accession_id]
length(recent_records) # 48110 record not available in GISAID download package

# dataframe with recently submitted records that are not yet in GISAID metadata package download
d_extra = download_GISAID_records(accession_ids = recent_records,
                                  get_sequence = FALSE, 
                                  clean_up = FALSE,
                                  target_dir = target_dir, # TO DO: check if directory exists, and if not make it
                                  max_batch_size = 10000, # maximum batch size, usually either 10000 or 5000
                                  headless = FALSE,
                                  usr = Sys.getenv("GISAIDR_USERNAME"),
                                  psw = Sys.getenv("GISAIDR_PASSWORD"))
dim(d_extra) # 48110                                               19

# merge GISAID download package & recently submitted records
GISAID = dplyr::bind_rows(GISAID, d_extra)
nrow(GISAID) # 13389714

# LOAD COG-UK DATA FOR THE UNITED KINGDOM ####
if (use_coguk) { coguk = download_COGUK_meta()
  # MERGE GISAID (MINUS UK GISAID DATA) & COG-UK DATA FOR UK
  GISAID = dplyr::bind_rows(GISAID[GISAID$country!="United Kingdom",], 
                           coguk)
  GISAID$country = factor(GISAID$country)
  GISAID$country = factor(GISAID$country)
  }
nrow(GISAID) # 13312997

levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
GISAID$continent= factor(GISAID$continent, levels=levels_continents)
GISAID$country = factor(GISAID$country)
levels_countries = levels(GISAID$country)
length(levels_countries) # 217
GISAID$location = factor(GISAID$location)
levels_locations = levels(GISAID$location)
length(levels_locations) # 587


# PARSE GISAID DATA ####

# parse date & check dates are valid
# records with valid date
GISAID$date_isvalid = (str_count(GISAID$collection_date,
                                 pattern = "-")==2)
GISAID$date = as.Date(NA)
GISAID$date[which(GISAID$date_isvalid)] = as.Date(fast_strptime(GISAID$collection_date[which(GISAID$date_isvalid)], "%Y-%m-%d")) # faster than as.Date(GISAID$collection_date)


# CODE VARIANT LINEAGES ####

sum(is.na(GISAID$aa_substitutions)) # 11323
GISAID$aa_substitutions[is.na(GISAID$aa_substitutions)] = ""

# # convert AA substitions to nested column "muts"
# system.time(GISAID <- GISAID %>% 
#   # convert mutations to nested list column, you can unnest this again using unnest(aa_substitutions)
#   mutate(muts = strsplit(aa_substitutions, ","))) # 344 s
# splitting unique mutations into separate columns first did not speed up things
# https://stackoverflow.com/questions/73758344/fast-way-to-split-comma-separated-strings-into-sparse-boolean-matrix-in-r

# some helper functions to simplify regular expressions below ####
# is a mutation present?
mut = function (mutation, x=GISAID$aa_substitutions) { str_detect(x, fixed(mutation)) }
# system.time(output <- mut("NSP3_S403L")) # 2.6s

# is one of a number of mutations present?
mut_oneof = function (mutations, x=GISAID$aa_substitutions) { rowSums(sapply(mutations, function (mutation) str_detect(x, fixed(mutation))))>=1 }
# system.time(output <- mut_oneof(c("NSP3_S403L","NSP8_N118S"))) # 5s

# are all of a given nr of mutations present?
mut_allof = function (mutations, x=GISAID$aa_substitutions) { rowSums(sapply(mutations, function (mutation) str_detect(x, fixed(mutation))))==length(mutations) }
# system.time(output <- mut_allof(c("NSP3_S403L","NSP8_N118S"))) # 5s

# lineage
lin = function (lineage, x=GISAID$pango_lineage) { pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$")
                                                   grepl(pat, x, fixed=F, perl=T) }

# lineage plus sublineages
linplus = function (lineage, x=GISAID$pango_lineage) { pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$","|",paste0("^", gsub(".","\\.",lineage, fixed=T)))
                                                       grepl(pat, x, fixed=F, perl=T) }

# one of X lineages
lin_oneof = function (lineages, x=GISAID$pango_lineage) { 
  rowSums(sapply(lineages, function (lineage) { 
    pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$")
    grepl(pat, x, fixed=F, perl=T) }))>=1
   }
# system.time(output <- lin_oneof(c("BA.5","BF.7"))) # 1.2s

# one of X lineages plus sublineages
linplus_oneof = function (lineages, x=GISAID$pango_lineage) { 
  rowSums(sapply(lineages, function (lineage) { 
    pat = paste0("^", gsub(".","\\.",lineage, fixed=T), "$","|",paste0("^", gsub(".","\\.",lineage, fixed=T)))
    grepl(pat, x, fixed=F, perl=T) }))>=1
  }
# system.time(output <- linplus_oneof(c("BA.5","BF.7"))) # 1.4s

# date from XX  
datefrom = function (date, x=GISAID$date) { x >= as.Date(date) }
# system.time(output <- datefrom(c("2021-01-01"))) # 0.14s

# save.image("~/Github/LineageExplorer/LineageExplorer/environment_5 oct 2022.RData")
# load("~/Github/LineageExplorer/LineageExplorer/environment_5 oct 2022.RData")

# DEFINE VARIANT LINEAGES OF INTEREST

# PS this currently crashes my R session - not sure why...
system.time(GISAID$variant <- case_when(
  (linplus_oneof(c("BA.2.75.2","BL.1"))|
    mut_allof(c("NSP3_S403L","NSP8_N118S","Spike_R346T","Spike_F486S"))|
    mut_allof(c("NSP3_S403L","E_T11A","Spike_R346T","Spike_F486S")))&
    datefrom("2022-04-01") ~ "Omicron (BA.2.75.2)", # also add CA.1 ( BA.2.75.2 + Spike L452R)?
  (linplus("BA.2.75")|
    mut_allof(c("NSP3_S403L","NSP8_N118S"))|
    mut_allof(c("NSP3_S403L","E_T11A")))&
    datefrom("2022-04-01") ~ "Omicron (BA.2.75)", # =22D, all other BA.2.75
  linplus("BA.2.76")|
    mut_allof(c("Spike_Y248N","Spike_R346T")) ~ "Omicron (BA.2.76)",
  linplus("BA.2.3")&
    mut_allof(c("Spike_K444R","Spike_L452M","Spike_N460K")) ~ "Omicron (BA.2.3.20)",
  linplus("BA.4.6")&
    datefrom("2021-12-01") ~ "Omicron (BA.4.6)",
  lin_oneof(c("BE.1.1"))&
    mut_allof(c("Spike_K444T","Spike_N460K","Spike_R346T"))&
    datefrom("2022-02-01") ~ "Omicron (BQ.1.1)",
  lin_oneof(c("BE.1.1"))&
    mut_allof(c("Spike_K444T","Spike_N460K"))&
    datefrom("2022-02-01") ~ "Omicron (BQ.1)",
  lin_oneof(c("BJ.1","BA.2.10"))&mut_allof(c("M_D3Y","N_T282I")) ~ "Omicron (BJ.1)", # PS some BJ.1 get assigned as BA.2.10
  lin_oneof(c("BA.2","BA.2.10"))&mut_allof(c("Spike_H146Q","Spike_G252V","E_T11A")) ~ "Omicron (XBB)",
  linplus_oneof(c("BF.7"))&
    datefrom("2022-02-01") ~ "Omicron (BF.7)", # BF.7=BA.5.2.1.7, also potentially interesting: BF.5=BA.5.2.1.5, BF.10=BA.5.2.1.10, BF.14=BA.5.2.1.14
  linplus_oneof("BA.5.2")&
    datefrom("2022-02-01") ~ "Omicron (BA.5.2)",
  linplus_oneof(c("B.1.617.2","AY"))&
    datefrom("2020-10-30") ~ "Delta",
  linplus("B.1.1.7")&
    datefrom("2020-09-20") ~ "Alpha",
  linplus("B.1.351")&
    datefrom("2020-08-10") ~ "Beta",
  linplus("BA.1")&
    datefrom("2021-09-01") ~ "Omicron (BA.1)", # = 21K
  (linplus("BA.4")|
    (linplus("BA.2")&
    mut_allof(c("Spike_L452R","Spike_F486V","NS7b_L11F"))&
    (!mut("M_D3N"))))&
    datefrom("2021-12-01") ~ "Omicron (BA.4)", # = 22A
  (linplus_oneof(c("BA.5","BE","BF"))|
     (!lin("Unassigned")&mut("M_D3N")))&
    datefrom("2021-12-01") ~ "Omicron (BA.5)", # =22B, cf pattern used by Alex Selby
  lin("BA.2.12.1")&
    datefrom("2021-12-01") ~ "Omicron (BA.2.12.1)", # =22C
  linplus("BA.2.38")&
    datefrom("2021-09-01") ~ "Omicron (BA.2.38)",
  linplus("BA.2")&
    datefrom("2021-09-01") ~ "Omicron (BA.2)", # =21L
  linplus("B.1.177")&
    datefrom("2020-05-27") ~ "B.1.177 (EU1)",
  linplus("B.1.160")&
    datefrom("2020-02-15") ~ "B.1.160 (EU2)",
  linplus("B.1.221") ~ "B.1.221 (20A/S:98F)",
  !lin("Unassigned") ~ "Other" # assigns NA to remaining Unassigned & remove them later on
  # TRUE ~ "Other" # alternative: to assign Unassigned lineages to category Other
)) # 140s - note: could be sped up by using multidplyr & parallelization

# for all Spike_N460K lineages together:
# cf https://twitter.com/rquiroga777/status/1575564008501182464
# S:R346T | S:K356T | S:K444T | S:F486V | S:F486S
# mut(Spike_N460K) & mut_oneof(Spike_R346T, Spike_K356T, Spike_K444T,Spike_F486V,Spike_F486S)
# "& datefrom("2022-07-01")

# earliest realistic dates were taken from
# https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/clade_emergence_dates.tsv

# note: in India BA.2.38 & BA.2.38.1 caused an infection wave in some states - hence separated out above

# B.1.177+B.1.160+B.1.221 were behind the 2020 wave in fall in Europe & each had one spike mutations & a small growth rate advantage relative to predominant B.1.1

# variants to watch & maybe add in due time:
# https://cov-spectrum.org/collections/1
# https://cov-spectrum.org/collections/32

table(GISAID$variant)

# GISAID = GISAID[!is.na(GISAID$variant),]
nrow(GISAID) # 13164288

# define variant lineages and colours ####
sel_reference_VOC = "Omicron (BA.5.2)"
levels_VARIANTS = c(sel_reference_VOC, "B.1.177 (EU1)", "B.1.160 (EU2)", "B.1.221 (20A/S:98F)", "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.38)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.2.76)", "Omicron (BA.2.75)", "Omicron (BA.2.75.2)", "Omicron (BA.5)", "Omicron (BF.7)", "Omicron (BQ.1)", "Omicron (BQ.1.1)", "Omicron (BJ.1)", "Omicron (XBB)", "Omicron (BA.2.3.20)")
levels_VARIANTS_plot = c("Other", "B.1.177 (EU1)", "B.1.160 (EU2)", "B.1.221 (20A/S:98F)", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.38)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.5)", "Omicron (BA.5.2)", "Omicron (BF.7)", "Omicron (BA.2.76)", "Omicron (BA.2.75)", "Omicron (BA.2.75.2)", "Omicron (BJ.1)", "Omicron (XBB)", "Omicron (BA.2.3.20)", "Omicron (BQ.1)", "Omicron (BQ.1.1)")

n = length(levels_VARIANTS)
lineage_cols_plot = case_when(
  levels_VARIANTS_plot=="Other" ~ "grey65",
  levels_VARIANTS_plot=="B.1.177 (EU1)" ~ "darkorange4",
  levels_VARIANTS_plot=="B.1.160 (EU2)" ~ "darkorange3",
  levels_VARIANTS_plot=="B.1.221 (20A/S:98F)" ~ "darkorange2",
  levels_VARIANTS_plot=="Beta" ~ "green4",
  levels_VARIANTS_plot=="Alpha" ~ "#0085FF",
  levels_VARIANTS_plot=="Delta" ~ "mediumorchid",
  levels_VARIANTS_plot=="Omicron (BA.1)" ~ "red",
  levels_VARIANTS_plot=="Omicron (BA.2)" ~ "red3",
  levels_VARIANTS_plot=="Omicron (BA.2.12.1)" ~ "black",
  levels_VARIANTS_plot=="Omicron (BA.2.38)" ~ "red4",
  levels_VARIANTS_plot=="Omicron (BA.4)" ~ "green3",
  levels_VARIANTS_plot=="Omicron (BA.4.6)" ~ "green2",
  levels_VARIANTS_plot=="Omicron (BA.5)" ~ "blue4",
  levels_VARIANTS_plot=="Omicron (BA.5.2)" ~ "blue3",
  levels_VARIANTS_plot=="Omicron (BF.7)" ~ "dodgerblue",
  levels_VARIANTS_plot=="Omicron (BA.2.76)" ~ "magenta4",
  levels_VARIANTS_plot=="Omicron (BA.2.75)" ~ "magenta3",
  levels_VARIANTS_plot=="Omicron (BA.2.75.2)" ~ "magenta",
  levels_VARIANTS_plot=="Omicron (BJ.1)" ~ "yellow3",
  levels_VARIANTS_plot=="Omicron (XBB)" ~ "yellow2",
  levels_VARIANTS_plot=="Omicron (BA.2.3.20)" ~ "orange",  
  levels_VARIANTS_plot=="Omicron (BQ.1)" ~ "cyan3",
  levels_VARIANTS_plot=="Omicron (BQ.1.1)" ~ "cyan"
)
pal.bands(lineage_cols_plot)
# pal.volcano(lineage_cols_plot)
# pal.zcurve(lineage_cols_plot)
lineage_cols = lineage_cols_plot[match(levels_VARIANTS,levels_VARIANTS_plot)]



# GISAID SELECTION ####

GISAID_sel = GISAID

# remove records with invalid/incomplete dates ####
GISAID_sel = GISAID_sel[GISAID_sel$date_isvalid,]
GISAID_sel = GISAID_sel[which(GISAID_sel$host=="Human"),]
nrow(GISAID_sel) # 13078490

# filter to desired date range ####
# start_date = "2020-06-01"
# end_date = today
# GISAID_sel = GISAID_sel[GISAID_sel$date>=as.Date(start_date)&
#                          GISAID_sel$date<=as.Date(end_date),]

# range(GISAID_sel$date, na.rm=T) # "2020-06-01" "2022-08-09"
# nrow(GISAID_sel) # 11958660

# add week, year & start of week
# fix a wrong date
GISAID_sel[which(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$date==as.Date("2022-01-07")),"date"] = as.Date("2022-07-01") # had day & month flipped around
# min(GISAID$date[GISAID$variant=="Omicron (BA.2.75)"],na.rm=T) # "2022-05-26"

# TO DO: add these columns in download functions??
GISAID_sel$Week = lubridate::week(GISAID_sel$date)
GISAID_sel$Year = lubridate::year(GISAID_sel$date)
GISAID_sel$Year_Week = interaction(GISAID_sel$Year,GISAID_sel$Week)
GISAID_sel$floor_date = fast_strptime(as.character(cut(GISAID_sel$date, "week")), "%Y-%m-%d") # start of week
GISAID_sel$DATE_NUM = as.numeric(GISAID_sel$date)

table(GISAID_sel$variant)
table(GISAID_sel$continent, GISAID_sel$variant)
# table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","country"])

# nBA_2_75 = sum(GISAID_sel$variant=="Omicron (BA.2.75)", na.rm=T)
# nBA_2_75 # 8337 BA.2.75 records so far with valid date
# 
# nBA_2_75_india = sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="India", na.rm=T)
# nBA_2_75_india # 5489 BA.2.75 records so far for India with valid date

# sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="United Kingdom", na.rm=T)
# 206 BA.2.75 in UK so far

# sum(GISAID_sel$variant=="Omicron (BQ.1)"&GISAID_sel$country=="United Kingdom", na.rm=T)
# 69 BQ.1 in UK so far

maxsubmdate = today


# selected countries to include, here those with min 5 BA.2.75.2 or 5 BQ.1* or 5 BA.2.3.20 sequences
tab = as.data.frame(table(GISAID_sel$country, GISAID_sel$variant))
sel_countries = sort(unique(c(tab[tab$Var2=="Omicron (BA.2.75.2)"&tab$Freq>=5,"Var1"],
                              tab[tab$Var2=="Omicron (BQ.1)"&tab$Freq>=5,"Var1"],
                              tab[tab$Var2=="Omicron (BQ.1.1)"&tab$Freq>=5,"Var1"],
                              tab[tab$Var2=="Omicron (BJ.1)"&tab$Freq>=5,"Var1"],
                              tab[tab$Var2=="Omicron (XBB)"&tab$Freq>=5,"Var1"],
                              tab[tab$Var2=="Omicron (BA.2.3.20)"&tab$Freq>=5,"Var1"])))
sel_countries
# [1] Australia      Austria        Bangladesh     Belgium        Brunei         Canada         Denmark       
# [8] France         Germany        India          Ireland        Israel         Italy          Japan         
# [15] Netherlands    New Zealand    Nigeria        Portugal       Singapore      Slovakia       South Korea   
# [22] Spain          Sweden         Switzerland    Thailand       United Kingdom USA 

GISAID_sel = as.data.frame(GISAID_sel)
GISAID_sel = GISAID_sel[GISAID_sel$country %in% sel_countries,]
GISAID_sel$country = droplevels(GISAID_sel$country)
GISAID_sel$continent = factor(GISAID_sel$continent, levels=unique(GISAID_sel$continent))


# 2. GLOBAL ANALYSIS OF THE GISAID+COGUK DATA USING MULTINOMIAL FITS ####

# TO DO: replace this with more elegant/tidy dplyr code? 

# AGGREGATED DATA BY DATE & COUNTRY ####
data_agbydatecountry1 = as.data.frame(table(GISAID_sel$date, GISAID_sel$country, GISAID_sel$variant))
colnames(data_agbydatecountry1) = c("date", "country", "variant", "count")
data_agbydatecountry1_sum = aggregate(count ~ date + country, data=data_agbydatecountry1, sum)
data_agbydatecountry1$total = data_agbydatecountry1_sum$count[match(interaction(data_agbydatecountry1$date,data_agbydatecountry1$country),
                                                                    interaction(data_agbydatecountry1_sum$date,data_agbydatecountry1_sum$country))]
data_agbydatecountry1$date = as.Date(as.character(data_agbydatecountry1$date))
data_agbydatecountry1$variant = factor(data_agbydatecountry1$variant, levels=levels_VARIANTS)
data_agbydatecountry1$date_num = as.numeric(data_agbydatecountry1$date)
data_agbydatecountry1$prop = data_agbydatecountry1$count/data_agbydatecountry1$total
data_agbydatecountry1 = data_agbydatecountry1[data_agbydatecountry1$total!=0,]
data_agbydatecountry1$country = factor(data_agbydatecountry1$country)
data_agbydatecountry1$continent = GISAID_sel$continent[match(data_agbydatecountry1$country, GISAID_sel$country)]
data_agbydatecountry1$continent = factor(data_agbydatecountry1$continent)
# write.csv(data_agbydatecountry1, file=".//data//GISAID//GISAID aggregated counts by date and lineage_all.csv", row.names=F)

# AGGREGATED DATA BY WEEK ####
data_agbyweek1 = as.data.frame(table(GISAID_sel$date, GISAID_sel$variant))
colnames(data_agbyweek1) = c("date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ date, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(interaction(data_agbyweek1$date),
                                                                    interaction(data_agbyweek1_sum$date))]
data_agbyweek1$date = as.Date(as.character(data_agbyweek1$date))
data_agbyweek1$variant = factor(data_agbyweek1$variant, levels=levels_VARIANTS_plot)
data_agbyweek1$date_num = as.numeric(data_agbyweek1$date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1 = data_agbyweek1[data_agbyweek1$total!=0,]
# write.csv(data_agbyweek1, file=".//data//GISAID//GISAID aggregated counts by week_all.csv", row.names=F)

# AGGREGATED DATA BY WEEK & COUNTRY ####
data_agbyweekcountry1 = as.data.frame(table(GISAID_sel$floor_date, GISAID_sel$country, GISAID_sel$variant))
colnames(data_agbyweekcountry1) = c("floor_date", "country", "variant", "count")
data_agbyweekcountry1_sum = aggregate(count ~ floor_date + country, data=data_agbyweekcountry1, sum)
data_agbyweekcountry1$total = data_agbyweekcountry1_sum$count[match(interaction(data_agbyweekcountry1$floor_date,data_agbyweekcountry1$country),
                                                                  interaction(data_agbyweekcountry1_sum$floor_date,data_agbyweekcountry1_sum$country))]
data_agbyweekcountry1$collection_date = as.Date(as.character(data_agbyweekcountry1$floor_date))
data_agbyweekcountry1$variant = factor(data_agbyweekcountry1$variant, levels=levels_VARIANTS)
data_agbyweekcountry1$collection_date_num = as.numeric(data_agbyweekcountry1$collection_date)
data_agbyweekcountry1$DATE_NUM = data_agbyweekcountry1$collection_date_num
data_agbyweekcountry1$prop = data_agbyweekcountry1$count/data_agbyweekcountry1$total
data_agbyweekcountry1 = data_agbyweekcountry1[data_agbyweekcountry1$total!=0,]
data_agbyweekcountry1$floor_date = NULL
# unmark to start here & read in aggregated counts
# data_agbyweekcountry1 = read.csv(file="./data/GISAID/GISAID aggregated counts by start of week and lineage.csv") 
data_agbyweekcountry1$country = factor(data_agbyweekcountry1$country)
if (is.null(data_agbyweekcountry1$continent)) data_agbyweekcountry1$continent = GISAID_sel$continent[match(data_agbyweekcountry1$country, GISAID_sel$country)]
data_agbyweekcountry1$continent = factor(data_agbyweekcountry1$continent)
data_agbyweekcountry1$collection_date = as.Date(data_agbyweekcountry1$collection_date)
data_agbyweekcountry1$variant = factor(data_agbyweekcountry1$variant, levels=levels_VARIANTS)
# write.csv(data_agbyweekcountry1, file="./data/GISAID/GISAID aggregated counts by start of week and lineage_all.csv", row.names=F)


# MULLER PLOT (RAW DATA, selected countries pooled, but with big sampling biases across countries)
muller_raw_all = ggplot(data=data_agbyweek1, aes(x=date, y=count, group=variant)) +
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES",
       subtitle=paste0("Raw GISAID data up to ",today," plus COG-UK data, pooled over all countries")) +
 coord_cartesian(xlim=c(min(GISAID_sel$date),max(GISAID_sel$date)), expand=c(0))
muller_raw_all

ggsave(file=file.path("plots", plotdir,"global analysis_muller plot_raw data all countries pooled.png"), width=7, height=5)


# FIT NNET::MULTINOM MULTINOMIAL SPLINE MODEL ####

set.seed(1)
# best model:
system.time(fit_global_multi <- nnet::multinom(variant ~ 
                                    ns(DATE_NUM, df=2)+
                                    ns(DATE_NUM, df=2):continent+
                                    country, 
                                  weights=count, 
                                  data=data_agbyweekcountry1, 
                                  maxit=10000, MaxNWts=100000)) # 935s; 5h if we use all GISAID data from all countries
# we calculate the Hessian using my own faster Rcpp Kronecker-product based function
source(".//fastmultinomHess.R") # faster way to calculation Hessian of multinomial fits
system.time(fit_global_multi$Hessian <- fastmultinomHess(fit_global_multi, model.matrix(fit_global_multi))) # 54s
# we add variance-covariance matrix as extra slot to be re-used later
system.time(fit_global_multi$vcov <- vcov(fit_global_multi)) # 1.75s

# with saved environment you can start from here ####

# clean up some memory
rm(d_extra, recent_records, coguk) 
# mem_usage = inspect_mem(GISAID_sel)
# print(mem_usage, n=100)
GISAID_sel$aa_substitutions = NULL # takes up 3 Gb
GISAID_sel$virus_name = NULL # takes up 1 Gb
gc()

# save.image("~/Github/LineageExplorer/environment_2022_10_06_small.RData")
# load("~/Github/LineageExplorer/environment_2022_10_06_small.RData")




# CALCULATE GROWTH RATE ADVANTAGE OVER REFERENCE LEVEL BA.5.2 ####

# with new faster marginaleffects code
system.time(meffects <- marginaleffects(fit_global_multi, 
                                               type = "link", # = additive log-ratio = growth rate advantage relative to BA.5.2
                                               variables = c("DATE_NUM"),
                                               by = c("group"),
                                               vcov = fit_global_multi$vcov,
                                               newdata = datagrid(DATE_NUM = today_num 
                                               ))) # 21s

# growth rate advantage compared to reference level BA.5.2 by continent
system.time(meffects_bycontinent <- marginaleffects(fit_global_multi, 
                                               type = "link", # = additive log-ratio = growth rate advantage relative to BA.5.2
                                               variables = c("DATE_NUM"),
                                               by = c("group", "continent"),
                                               vcov = fit_global_multi$vcov,
                                               newdata = datagrid(DATE_NUM = today_num,
                                                                  continent = unique(data_agbyweekcountry1$continent)
                                               ))) # 23s

# for all pairwise growth rate differences:
# growth_differences = comparisons(
#   fit_global_multi,
#   newdata = datagrid(DATE_NUM = today_num),
#   variables = "DATE_NUM",
#   by = "continent",
#   type = "clr", # here we could either use "clr" (centered logratio) or "link" (additive logratio) - this gives same result
#   hypothesis = "pairwise")


# old emtrends code to calculate pairwise growth rate differences
# system.time(emtr_pairw <- emtrends(fit_global_multi, revpairwise ~ variant,
#                                    by="continent",
#                                    var="DATE_NUM",  mode="latent",
#                                    at=list(DATE_NUM=today_num))) # 
# delta_r_pairw = data.frame(confint(emtr_pairw,
#                                    adjust="none", df=NA)$contrasts,
#                            p.value=as.data.frame(emtr_pairw$contrasts)$p.value)
# delta_r_pairw
# write.csv(delta_r_pairw, file.path("plots", plotdir, "growth rate advantage all variants vs BA_5_2.csv"), row.names=F)


# plot of growth rate advantage of last n newest variants
# TO DO: order by selective advantage and then take top n
lastn = 11
sel_variants = tail(levels_VARIANTS,lastn)
sel_variants = sel_variants[!sel_variants %in% c(sel_reference_VOC, "Omicron (BA.5)", "Omicron (BA.2.76)")]
meffects_sel1 = meffects[meffects$group %in% sel_variants,]
meffects_sel1$group = factor(meffects_sel1$group, levels=meffects_sel1$group[order(meffects_sel1$dydx, decreasing=T)])
cols = colorRampPalette(c("red3", "blue3"))(length(levels(meffects_sel1$group)))
qplot(data=meffects_sel1, 
      x=group, y=dydx*100, ymin=conf.low*100, ymax=conf.high*100, fill=group, geom="col", 
      width=I(0.7)) +
  geom_linerange(aes(lwd=I(0.4))) + ylab("Growth rate advantage\nrelative to BA.5.2 (% per day)") +
  scale_fill_manual(values=cols) +
  theme(legend.position="none") + xlab("") +
  ggtitle("GROWTH RATE ADVANTAGE OF SARS-CoV2 VARIANTS",
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nGISAID & COG-UK data, using data from countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2\n(",
                   paste0(sel_countries[1:10], collapse=", "), "\n",
                   paste0(sel_countries[11:18], collapse=", "), "\n",
                   paste0(sel_countries[19:length(sel_countries)], collapse=", "), ")") ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
  
ggsave(file=file.path("plots", plotdir,"growth rate advantage VOCs_global fit.png"), width=9, height=7)

# plot of growth rate advantage of last X newest variants by continent
# TO DO: order by selective advantage and then take top n
lastn = 11
sel_variants = tail(levels_VARIANTS,lastn)
sel_variants = sel_variants[!sel_variants %in% c(sel_reference_VOC, "Omicron (BA.5)", "Omicron (BA.2.76)")]
sel_continents = unique(data_agbyweekcountry1$continent)
sel_continents = sel_continents[!sel_continents %in% c("Africa")] # too little data
meffects_sel2 = meffects_bycontinent[meffects_bycontinent$continent %in% sel_continents,]
meffects_sel2 = meffects_sel2[meffects_sel2$group %in% sel_variants,]
meffects_sel2$group = factor(meffects_sel2$group, levels=levels(meffects_sel1$group))
outlier = abs(meffects_sel2$dydx)>=0.40 # typically due to there being too little data
meffects_sel2 = meffects_sel2[!outlier,]
qplot(data=meffects_sel2, 
      x=group, y=dydx*100, ymin=conf.low*100, ymax=conf.high*100, fill=group, geom="col", 
      width=I(0.7)) +
  facet_wrap(~ continent) +
  geom_linerange(aes(lwd=I(0.4))) + ylab("Growth rate advantage\nrelative to BA.5.2 (% per day)") +
  scale_fill_manual(values=cols) +
  theme(legend.position="none") + xlab("") +
  ggtitle("GROWTH RATE ADVANTAGE OF SARS-CoV2 VARIANTS",
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nGISAID & COG-UK data, using data from countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2\n(",
                          paste0(sel_countries[1:10], collapse=", "), "\n",
                          paste0(sel_countries[11:18], collapse=", "), "\n",
                          paste0(sel_countries[19:length(sel_countries)], collapse=", "), ")") ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(file=file.path("plots", plotdir,"growth rate advantage VOCs_global fit_by continent.png"), width=12, height=8)

qplot(data=meffects_sel2[meffects_sel2$continent=="Europe",], 
      x=group, y=dydx*100, ymin=conf.low*100, ymax=conf.high*100, fill=group, geom="col", 
      width=I(0.7)) +
  # facet_wrap(~ continent) +
  geom_linerange(aes(lwd=I(0.4))) + ylab("Growth rate advantage\nrelative to BA.5.2 (% per day)") +
  scale_fill_manual(values=cols) +
  theme(legend.position="none") + xlab("") +
  ggtitle("GROWTH RATE ADVANTAGE OF SARS-CoV2 VARIANTS IN EUROPE",
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nGISAID & COG-UK data, using data from countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2\n(",
                          paste0(sel_countries[1:10], collapse=", "), "\n",
                          paste0(sel_countries[11:18], collapse=", "), "\n",
                          paste0(sel_countries[19:length(sel_countries)], collapse=", "), ")") ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(file=file.path("plots", plotdir,"growth rate advantage VOCs_global fit_europe.png"), width=12, height=8)


# PLOT MULTINOMIAL FIT ####

extrapolate = 60
date.from = as.numeric(as.Date("2021-01-01"))
date.to = today_num+extrapolate

# multinomial model predictions by country with CIs calculated using margineffects::predictions

predgrid = expand.grid(list(DATE_NUM=as.numeric(seq(date.from, date.to)),
                            country=unique(data_agbyweekcountry1$country))) # unique(GISAID_sel$country)
predgrid$continent = data_agbyweekcountry1$continent[match(predgrid$country,
                                                                data_agbyweekcountry1$country)]
# note: now using Delta method on response scale, better to
# calculate CIs as in Effects package on link scale (type="link") or
# on isometric logratio scale (type="ilr") & then backtransform
# but still having some problems with over/underflows with
# type="link" and type="ilr" is a bit more hassle to backtransform

system.time(fit_preds <- data.frame(predictions(fit_global_multi, 
                       newdata = predgrid,
                       type = "probs",
                       vcov = fit_global_multi$vcov))) # %>% # 
             # transform(conf.low = predicted - 1.96 * std.error,
             #          conf.high = predicted + 1.96 * std.error) %>%
             # group_by(rowid) |>
             #mutate_at(c("predicted", "conf.low", "conf.high"), function (x) plogis(x)))
             # 442s
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
fit_preds$date = as.Date(fit_preds$DATE_NUM, origin="1970-01-01")
fit_preds$variant = NULL
colnames(fit_preds)[which(colnames(fit_preds)=="group")] = "variant"
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)

# TO DO: fix bug with type="link" where some predictions come out as NA,
# and/or switch to type="ilr" - check in my marginaleffects fork

fit_preds$country = factor(fit_preds$country)
levels(fit_preds$country)
fit_preds$continent = factor(fit_preds$continent)

write_csv(fit_preds, file=file.path("plots", plotdir, "GISAID fitted lineage frequencies global multinomial spline fit.csv"))


# PLOT MULTINOMIAL FIT ON LOGIT SCALE ####

plot_preds_logit = qplot(data=fit_preds[fit_preds$variant!="Other",], 
                         x=date, y=predicted, geom="blank") +
  facet_wrap(~ country) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&
                                          data_agbyweekcountry1$country %in% sel_countries,],
             aes(x=collection_date, y=prop, size=total,
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
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
plot_preds_logit

ggsave(file=file.path("plots", plotdir, "global multinom fit_all data_predictions_logit scale.png"), width=20, height=12)

# zoomed in on last 6 months
sel_variants = tail(levels_VARIANTS_plot,-11)
plot_preds_logit = qplot(data=fit_preds[fit_preds$variant %in% sel_variants,], 
                         x=date, y=predicted, geom="blank") +
  facet_wrap(~ country) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-11)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-11)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant %in% sel_variants&
                                          data_agbyweekcountry1$country %in% sel_countries,],
             aes(x=collection_date, y=prop, size=total,
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
plot_preds_logit

ggsave(file=file.path("plots", plotdir, "global multinom fit_all data_predictions_last6months_logit scale.png"), width=20, height=12)


# plot just for Belgium
plot_preds_logit_be = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$country=="Belgium",], 
                         x=date, y=predicted, geom="blank") +
  # facet_wrap(~ country) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN BELGIUM",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nonly Belgium shown here")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$country=="Belgium",],
             aes(x=collection_date, y=prop, size=total,
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
plot_preds_logit_be

ggsave(file=file.path("plots", plotdir, "global multinom fit_belgium_logit scale.png"), width=12, height=5)

# plot just for UK
plot_preds_logit_uk = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$country=="United Kingdom",], 
                            x=date, y=predicted, geom="blank") +
  # facet_wrap(~ country) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant), alpha=I(0.3)) +
  geom_line(aes(colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN UNITED KINGDOM",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nonly UK shown here")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$country=="United Kingdom",],
             aes(x=collection_date, y=prop, size=total,
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
plot_preds_logit_uk

ggsave(file=file.path("plots", plotdir, "global multinom fit_UK_logit scale.png"), width=12, height=5)


# plot just for DK
plot_preds_logit_dk = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$country=="Denmark",], 
                            x=date, y=predicted, geom="blank") +
  # facet_wrap(~ country) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                   fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=predicted,
                colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN DENMARK",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nonly Denmark shown here")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$country=="Denmark",],
             aes(x=collection_date, y=prop, size=total,
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
plot_preds_logit_dk

ggsave(file=file.path("plots", plotdir, "global multinom fit_Denmark_logit scale.png"), width=12, height=5)


# plot just for Austria
plot_preds_logit_au = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$country=="Austria",], 
                            x=date, y=predicted, geom="blank") +
  # facet_wrap(~ country) +
  geom_ribbon(aes(y=predicted, ymin=conf.low, ymax=conf.high, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=predicted,
                colour=variant), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES IN AUSTRIA",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nonly Denmark shown here")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$country=="Denmark",],
             aes(x=collection_date, y=prop, size=total,
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
plot_preds_logit_au

ggsave(file=file.path("plots", plotdir, "global multinom fit_Austria_logit scale.png"), width=12, height=5)




# plot predicted values as Muller plot
muller_fit = ggplot(data=fit_preds[fit_preds$date>=as.Date("2021-01-01"),], 
                    aes(x=date, y=predicted, group=variant)) +
  facet_wrap(~ country, ncol=5) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES",
       subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01")+5, 
           ymin=0, ymax=1, alpha=0.5, fill="white") + # extrapolated part
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
muller_fit

ggsave(file=file.path("plots", plotdir, "global multinom fit_all data_predictions_muller plot.png"), width=20, height=12)


# save.image("~/Github/LineageExplorer/environment_2022_09_29_small.RData")
# load("~/Github/LineageExplorer/environment_2022_09_29_small.RData")

# MAP VARIANT SHARE ONTO CASE NUMBERS ####

# sel_countries = c("Austria","Belgium","Denmark","France","Germany","Israel","Netherlands","Singapore","United Kingdom") # "Bangladesh",
sel_countries_fig = c("Austria", "Bangladesh", "Belgium", "Brunei", "Denmark", "France", "Germany", "Ireland", "Israel", "Italy", "Netherlands", "New Zealand", "Singapore", "Slovakia", "Switzerland", "United Kingdom")
country_data = get_national_data(countries=sel_countries_fig, # sel_countries_top3,
                                 source="who",
                                 level=1)
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY IN\nCOUNTRIES WITH >=5 SEQUENCED BA.2.75.2 or BQ.1* CASES",
          subtitle="7 day simple moving average, using\nWHO case data accessed via covidregionaldata package")
# ggsave(file=file.path("plots", plotdir, "new cases countries with more than 10 sequenced BA_2_75 cases.png"), width=7, height=5)

  
fit_preds[(fit_preds$date==today)&(fit_preds$variant=="Omicron (BQ.1)"),]
fit_preds[(fit_preds$date==today)&(fit_preds$variant=="Omicron (BQ.1)"),][fit_preds[(fit_preds$date==today)&(fit_preds$variant=="Omicron (BQ.1)"),]$predicted>0.2,]

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
fit_preds_sel = fit_preds_sel[fit_preds_sel$country %in% sel_countries_fig,]
fit_preds_sel$country = factor(fit_preds_sel$country)

fit_preds2 = fit_preds_sel
fit_preds2$cases[fit_preds2$cases==0] = NA
fit_preds2$cases[fit_preds2$cases<=1] = NA
fit_preds2$country = factor(fit_preds2$country)
fit_preds2$variant = factor(fit_preds2$variant, levels=levels_VARIANTS_plot)

ggplot(data=fit_preds2, 
       aes(x=date, y=cases)) + 
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_preds2, aes(x=date, y=totnewcases_smoothed, lwd=I(1.5)), 
            colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nselected countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  scale_colour_manual("lineage", values=lineage_cols_plot) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2021-01-01"),NA)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) 
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=file.path("plots", plotdir,"global multinom fit_all data_predictions_confirmed cases multinomial fit by country.png"), width=20, height=12)

ggplot(data=fit_preds2, 
       aes(x=date, y=cases)) + 
  facet_wrap(~ country, scale="free") +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_preds2, aes(x=date, y=totnewcases_smoothed, lwd=I(1.5)), 
            colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nselected countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  scale_colour_manual("lineage", values=lineage_cols_plot) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2021-01-01"),NA)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(xlim=c(as.Date("2022-01-01"),NA))

ggsave(file=file.path("plots", plotdir,"global multinom fit_all data_predictions_confirmed cases multinomial fit by country_zoomed.png"), width=20, height=12)

# stacked area chart
ggplot(data=fit_preds2, 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nselected countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=file.path("plots", plotdir,"\\global multinom fit_all data_predictions_confirmed cases stacked area multinomial fit by country.png"), width=20, height=12)

# zoomed in since jan 2022
ggplot(data=fit_preds2, 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nselected countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2022-01-01"),NA))

ggsave(file=file.path("plots", plotdir,"\\global multinom fit_all data_predictions_confirmed cases stacked area multinomial fit by country_zoomed.png"), width=20, height=12)

# zoomed in since june 2022
ggplot(data=fit_preds2[fit_preds2$date>=as.Date("2022-06-01"),], 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
            position="stack") +
  scale_fill_manual("", values=lineage_cols_plot[-c(2,3,4,5,6,7,8)]) + # TO DO: determine dropped levels automatically
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT",
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nselected countries with >=5 BQ.1, BQ.1.1, BJ.1, XBB, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2022-06-01"),NA))

ggsave(file=file.path("plots", plotdir,"\\global multinom fit_all data_predictions_confirmed cases stacked area multinomial fit by country_zoomed2.png"), width=20, height=12)

# save.image("~/Github/LineageExplorer/environment_2022_09_29_small.RData")
# load("~/Github/LineageExplorer/environment_2022_09_29_small.RData")

# TO DO : get IHME infection estimates from links below,
# convolve these to cases & map them onto variant frequencies
# https://www.healthdata.org/covid/data-downloads
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2020.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2021.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2022.csv

# TO DO : get mortality data from Eurostat with Eurostat package & https://www.mortality.org/Data/STMF
# convolve case data by variant to mortality data & calculate
# death toll of each variant

# TO DO : finish spatial multinomial tensor spline fit in function of latitude & longitude & time
# (especially good for countries with limited or no sequencing data)