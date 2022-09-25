# GLOBAL ANALYSIS OF SARS-Cov2 VARIANTS OF CONCERN & INTEREST
# T. Wenseleers
# last update 25 SEPTEMBER 2022

# set GISAID credentials ####
# set them first using 
# Sys.setenv(GISAIDR_USERNAME = "XXXXX")
# Sys.setenv(GISAIDR_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_USERNAME = "XXXXX") # not needed for this script
# Sys.setenv(GISAIDJSON_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_STREAM = "XXXXX")

# load some utility functions & install required packages 
library(devtools)
remotes::install_github("Wytamma/GISAIDR")
library(GISAIDR)
source(".//set_GISAID_credentials.R") # set GISAID credentials
source(".//download_GISAID.R") # load function to download GISAID metadata download package (lacking records from last few days)
source(".//download_GISAID_records.R") # load functions to download most recent GISAID records
source(".//download_COGUK.R") # load function to download COG-UK metadata
use_coguk = TRUE # use COG-UK data instead of GISAID data for UK?
source(".//fastmultinomHess.R") # faster way to calculation Hessian of multinomial fits

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

system.time(GISAID <- download_GISAD_meta(target_dir = target_dir,
                                          headless = FALSE,
                                          usr = Sys.getenv("GISAIDR_USERNAME"),
                                          psw = Sys.getenv("GISAIDR_PASSWORD"))) # 194s
download = tail(list.files(target_dir, pattern=".tar.xz"), 1)
download # "metadata_tsv_2022_09_24.tar.xz"  

# records go up to submission date
GISAID_max_submdate = as.Date(max(GISAID$submission_date, na.rm=T))
GISAID_max_submdate # "2022-09-22"
nrow(GISAID) # 13217765


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
length(recent_records) # 34680  records
recent_records = recent_records[!recent_records %in% GISAID$accession_id]
length(recent_records) # 25540 record not available in GISAID download package

# dataframe with recently submitted records that are not yet in GISAID metadata package download
d_extra = download_GISAID_records(accession_ids = recent_records,
                                  get_sequence = FALSE, 
                                  clean_up = FALSE,
                                  target_dir = target_dir, # TO DO: check if directory exists, and if not make it
                                  max_batch_size = 10000, # maximum batch size, usually either 10000 or 5000
                                  headless = FALSE,
                                  usr = Sys.getenv("GISAIDR_USERNAME"),
                                  psw = Sys.getenv("GISAIDR_PASSWORD"))
dim(d_extra) # 25540                                       19

# merge GISAID download package & recently submitted records
GISAID = dplyr::bind_rows(GISAID, d_extra)
nrow(GISAID) # 13243305

# LOAD COG-UK DATA FOR THE UNITED KINGDOM ####
if (use_coguk) { coguk = download_COGUK_meta()
  # MERGE GISAID (MINUS UK GISAID DATA) & COG-UK DATA FOR UK
  GISAID = dplyr::bind_rows(GISAID[GISAID$country!="United Kingdom",], 
                           coguk)
  GISAID$country = factor(GISAID$country)
  GISAID$country = factor(GISAID$country)
  }
nrow(GISAID) # 13164288

levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
GISAID$continent= factor(GISAID$continent, levels=levels_continents)
GISAID$country = factor(GISAID$country)
levels_countries = levels(GISAID$country)
length(levels_countries) # 217
GISAID$location = factor(GISAID$location)
levels_locations = levels(GISAID$location)
length(levels_locations) # 289


# PARSE GISAID DATA ####
# GISAID = as.data.frame(GISAID) 
# TO DO: keep as more efficient data.table or as tibble
# & make sure everything below still works

# parse date & check dates are valid
# records with valid date
GISAID$date_isvalid = (str_count(GISAID$collection_date,
                                 pattern = "-")==2)
GISAID$date = as.Date(NA)
GISAID$date[which(GISAID$date_isvalid)] = as.Date(fast_strptime(GISAID$collection_date[which(GISAID$date_isvalid)], "%Y-%m-%d")) # faster than as.Date(GISAID$collection_date)


# CODE VARIANT LINEAGES ####

sum(is.na(GISAID$aa_substitutions)) # 11134
GISAID$aa_substitutions[is.na(GISAID$aa_substitutions)] = ""

# # convert AA substitions to nested column "muts"
# system.time(GISAID <- GISAID %>% 
#   # convert mutations to nested list column, you can unnest this again using unnest(aa_substitutions)
#   mutate(muts = strsplit(aa_substitutions, ","))) # 344 s
# splitting unique mutations into separate columns first did not speed up things
# https://stackoverflow.com/questions/73758344/fast-way-to-split-comma-separated-strings-into-sparse-boolean-matrix-in-r

# some helper functions to simplify regular expressions below
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

# DEFINE VARIANT LINEAGES OF INTEREST
system.time(GISAID$variant <- case_when(
  (linplus("BA.2.75.2")|
    mut_allof(c("NSP3_S403L","NSP8_N118S","Spike_R346T","Spike_F486S"))|
    mut_allof(c("NSP3_S403L","E_T11A","Spike_R346T","Spike_F486S")))&
    datefrom("2022-04-01") ~ "Omicron (BA.2.75.2)",
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
  linplus_oneof(c("BE.1.1"))&
    mut_allof(c("Spike_K444T","Spike_N460K"))&
    datefrom("2022-02-01") ~ "Omicron (BQ.1)", # BQ.1 (BE.1.1=BA.5.3.1.1.1 + K444T + N460K), for BQ.1.1 +R346T
  linplus_oneof(c("BA.5.2.1","BF"))&
    datefrom("2022-02-01") ~ "Omicron (BA.5.2.1)", # includes BF.7=BA.5.2.1.7, BF.5=BA.5.2.1.5, BF.10=BA.5.2.1.10, BF.14=BA.5.2.1.14
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
)) # 86s - note: could be sped up by using multidplyr & parallelization

# earliest realistic dates were taken from
# https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/clade_emergence_dates.tsv

# note: in India BA.2.38 & BA.2.38.1 caused an infection wave in some states - hence separated out above

# B.1.177+B.1.160+B.1.221 were behind the 2020 wave in fall in Europe & each had one spike mutations & a small growth rate advantage relative to predominant B.1.1

# variants to watch & maybe add in due time:
# https://cov-spectrum.org/collections/1
# https://cov-spectrum.org/collections/32
# BJ.1, https://github.com/cov-lineages/pango-designation/issues/915
# (just 11 seqs for now, https://cov-spectrum.org/explore/World/AllSamples/Past6M/variants?variantQuery=%5B4-of%3A+ORF1a%3A47R%2C+S%3A83A%2C+S%3A146Q%2C+S%3A213E%2C+S%3A339H%2C+S%3A445P%2C+S%3A483A%2C+S%3A1003I%2C+M%3A3Y%2C+ORF7a%3A110T%2C+N%3A282I%2C+15738T%2C+15939C%5D&)

# sum(is.na(GISAID$variant)) # 224476 unassigned
# sum(GISAID$pango_lineage=="Unassigned", na.rm=T) # 224635 originally unassigned in GISAID
# sum(GISAID$variant=="Omicron (BA.2.75)", na.rm=T) # 8337 BA.2.75
# sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$country=="India", na.rm=T) # 5489 BA.2.75 for India
# sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$date_isvalid, na.rm=T) # 8337 BA.2.75 with valid date
# sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$country=="India"&GISAID$date_isvalid, na.rm=T) # 5489 BA.2.75 for India with valid date
# sum(GISAID$pango_lineage=="BA.2.75", na.rm=T) # 8786
# sum(GISAID$variant=="Omicron (BA.5)", na.rm=T) # 474004
# sum(GISAID$variant=="Omicron (BA.5.2)", na.rm=T) # 188559
# sum(GISAID$variant=="Omicron (BA.5.2.1)", na.rm=T) # 281604
# sum(GISAID$variant=="Omicron (BQ.1)", na.rm=T) # 381
table(GISAID$variant)

# GISAID = GISAID[!is.na(GISAID$variant),]
nrow(GISAID) # 13164288

# define variant lineages and colours ####
sel_reference_VOC = "Omicron (BA.5.2)"
levels_VARIANTS = c(sel_reference_VOC, "B.1.177 (EU1)", "B.1.160 (EU2)", "B.1.221 (20A/S:98F)", "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.38)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.2.76)", "Omicron (BA.5)", "Omicron (BA.5.2.1)", "Omicron (BQ.1)", "Omicron (BA.2.75)", "Omicron (BA.2.75.2)","Omicron (BA.2.3.20)")
levels_VARIANTS_plot = c("Other", "B.1.177 (EU1)", "B.1.160 (EU2)", "B.1.221 (20A/S:98F)", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.38)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.5)", "Omicron (BA.5.2)", "Omicron (BA.5.2.1)", "Omicron (BQ.1)", "Omicron (BA.2.76)", "Omicron (BA.2.75)", "Omicron (BA.2.75.2)","Omicron (BA.2.3.20)")

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
  levels_VARIANTS_plot=="Omicron (BA.5.2)" ~ "blue2",
  levels_VARIANTS_plot=="Omicron (BA.5.2.1)" ~ "cyan4",
  levels_VARIANTS_plot=="Omicron (BQ.1)" ~ "cyan3",
  levels_VARIANTS_plot=="Omicron (BA.2.76)" ~ "magenta4",
  levels_VARIANTS_plot=="Omicron (BA.2.75)" ~ "magenta3",
  levels_VARIANTS_plot=="Omicron (BA.2.75.2)" ~ "magenta",
  levels_VARIANTS_plot=="Omicron (BA.2.3.20)" ~ "orange"  
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
nrow(GISAID_sel) # 12933396

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

GISAID_sel$Week = lubridate::week(GISAID_sel$date)
GISAID_sel$Year = lubridate::year(GISAID_sel$date)
GISAID_sel$Year_Week = interaction(GISAID_sel$Year,GISAID_sel$Week)
GISAID_sel$floor_date = fast_strptime(as.character(cut(GISAID_sel$date, "week")), "%Y-%m-%d") # start of week
GISAID_sel$DATE_NUM = as.numeric(GISAID_sel$date)

table(GISAID_sel$variant)
table(GISAID_sel$continent, GISAID_sel$variant)
table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","country"])

# nBA_2_75 = sum(GISAID_sel$variant=="Omicron (BA.2.75)", na.rm=T)
# nBA_2_75 # 8337 BA.2.75 records so far with valid date
# 
# nBA_2_75_india = sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="India", na.rm=T)
# nBA_2_75_india # 5489 BA.2.75 records so far for India with valid date

sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="United Kingdom", na.rm=T)
# 206 BA.2.75 in UK so far

sum(GISAID_sel$variant=="Omicron (BQ.1)"&GISAID_sel$country=="United Kingdom", na.rm=T)
# 69 BQ.1 in UK so far

maxsubmdate = today


# these countries are used as a selection for plotting later on
tab = as.data.frame(table(GISAID_sel$country, GISAID_sel$variant))

# countries with min 5 BA.2.75.2 or 5 BQ.1* or 5 BA.2.3.20 sequences
sel_countries_min5 = sort(unique(c(as.character(tab[tab$Var2=="Omicron (BA.2.75.2)"&tab$Freq>=5,"Var1"]),
                                   as.character(tab[tab$Var2=="Omicron (BQ.1)"&tab$Freq>=5,"Var1"]),
                                   as.character(tab[tab$Var2=="Omicron (BA.2.3.20)"&tab$Freq>=5,"Var1"]))))
sel_countries_min5
# [1] "Australia"      "Austria"        "Bangladesh"     "Belgium"        "Canada"        
# [6] "Denmark"        "France"         "Germany"        "India"          "Ireland"       
# [11] "Israel"         "Italy"          "Japan"          "Netherlands"    "New Zealand"   
# [16] "Singapore"      "South Korea"    "Switzerland"    "United Kingdom" "USA"  

sel_countries = sel_countries_min5


# 2. GLOBAL ANALYSIS OF ALL OF THE GISAID+COGUK DATA ####

# TO DO: replace this with more elegant/tidy dplyr coce? 

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

data_agbyweek1$variant = factor(data_agbyweek1$variant, levels=levels_VARIANTS_plot)


# MULLER PLOT (RAW DATA, all countries pooled, but with big sampling biases across countries)
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


# fit multinomial spline model ####
data_agbyweekcountry1_subs = data_agbyweekcountry1[data_agbyweekcountry1$country %in% 
                                                     sel_countries,]
data_agbyweekcountry1_subs$country = droplevels(data_agbyweekcountry1_subs$country)
data_agbyweekcountry1_subs$variant = factor(data_agbyweekcountry1_subs$variant, levels=levels_VARIANTS)

data_agbyweekcountry1_subs2 = data_agbyweekcountry1_subs[data_agbyweekcountry1_subs$continent=="Europe",]
data_agbyweekcountry1_subs2$country = droplevels(data_agbyweekcountry1_subs2$country)


set.seed(1)
# best model:
system.time(fit_global_multi <- nnet::multinom(variant ~ 
                                    ns(DATE_NUM, df=2)+
                                    ns(DATE_NUM, df=2):continent+
                                    country, 
                                  weights=count, 
                                  data=data_agbyweekcountry1_subs, 
                                  maxit=10000, MaxNWts=100000)) # 597s; 5h if we use all GISAID data from all countries
system.time(fit_global_multi$Hessian <- fastmultinomHess(fit_global_multi, model.matrix(fit_global_multi))) # 17s
system.time(fit_global_multi$vcov <- vcov(fit_global_multi)) # 0.6s



# calculate current pairwise growth rate differences compared to ref level BA.5.2 ####

# for all pairwise growth rate differences:
# growth_differences = comparisons(
#   fit_global_multi,
#   newdata = datagrid(DATE_NUM = today_num),
#   variables = "DATE_NUM",
#   by = "continent",
#   type = "clr",
#   hypothesis = "pairwise")


# WITH SAVED ENVIRONMENT YOU CAN START FROM HERE ####

# save.image("~/Github/LineageExplorer/environment_2022_09_25.RData")
# load("~/Github/LineageExplorer/environment_2022_09_25.RData")
 
# clean up some memory
rm(GISAID, d_extra, recent_records, coguk) 
# mem_usage = inspect_mem(GISAID_sel)
# print(mem_usage, n=100)
GISAID_sel$aa_substitutions = NULL # takes up 3 Gb
GISAID_sel$virus_name = NULL # takes up 1 Gb
gc()

# save.image("~/Github/LineageExplorer/environment_2022_09_25_small.RData")
# load("~/Github/LineageExplorer/environment_2022_09_25_small.RData")


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

# average growth rate advantage compared to reference level BA.5.2
# with new faster marginaleffects code
system.time(meffects <- marginaleffects(fit_global_multi, 
                                               type = "link", # = additive log-ratio = growth rate advantage relative to BA.5.2
                                               variables = c("DATE_NUM"),
                                               by = c("group"),
                                               vcov = fit_global_multi$vcov,
                                               newdata = datagrid(DATE_NUM = today_num 
                                               ))) # 13s


# growth rate advantage compared to reference level BA.5.2 by continent
system.time(meffects_bycontinent <- marginaleffects(fit_global_multi, 
                                               type = "link", # = additive log-ratio = growth rate advantage relative to BA.5.2
                                               variables = c("DATE_NUM"),
                                               by = c("group", "continent"),
                                               vcov = fit_global_multi$vcov,
                                               newdata = datagrid(DATE_NUM = today_num,
                                                                  continent = unique(data_agbyweekcountry1$continent)
                                               ))) # 16s

# plot of growth rate advantage of last X newest variants
lastn = 4
sel_variants = tail(levels_VARIANTS,lastn)
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
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nGISAID & COG-UK data, using data from countries with >=5 BQ.1*, BA.2.3.20 or BA.2.75.2\n(",
                   paste0(sel_countries_min5[1:10], collapse=", "), "\n",
                   paste0(sel_countries_min5[11:18], collapse=", "), "\n",
                   paste0(sel_countries_min5[19:length(sel_countries_min5)], collapse=", "), ")") ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
  # theme(axis.text.x = element_text(angle = 45, hjust=1))
  
ggsave(file=file.path("plots", plotdir,"growth rate advantage VOCs_global fit.png"), width=7, height=5)

# plot of growth rate advantage of last X newest variants by continent
lastn = 4
sel_variants = tail(levels_VARIANTS,lastn)
sel_continents = unique(data_agbyweekcountry1_subs$continent)
meffects_sel2 = meffects_bycontinent[meffects_bycontinent$continent %in% sel_continents,]
meffects_sel2 = meffects_sel2[meffects_sel2$group %in% sel_variants,]
meffects_sel2$group = factor(meffects_sel2$group, levels=levels(meffects_sel1$group))
qplot(data=meffects_sel2, 
      x=group, y=dydx*100, ymin=conf.low*100, ymax=conf.high*100, fill=group, geom="col", 
      width=I(0.7)) +
  facet_wrap(~ continent) +
  geom_linerange(aes(lwd=I(0.4))) + ylab("Growth rate advantage\nrelative to BA.5.2 (% per day)") +
  scale_fill_manual(values=cols) +
  theme(legend.position="none") + xlab("") +
  ggtitle("GROWTH RATE ADVANTAGE OF SARS-CoV2 VARIANTS",
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nGISAID & COG-UK data, using data from countries with >=5 BQ.1*, BA.2.3.20 or BA.2.75.2\n(",
                          paste0(sel_countries_min5[1:10], collapse=", "), "\n",
                          paste0(sel_countries_min5[11:18], collapse=", "), "\n",
                          paste0(sel_countries_min5[19:length(sel_countries_min5)], collapse=", "), ")") ) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) # +
# theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(file=file.path("plots", plotdir,"growth rate advantage VOCs_global fit_by continent.png"), width=12, height=8)



# PLOT MULTINOMIAL FIT

extrapolate = 60
date.from = as.numeric(as.Date("2020-01-01"))
date.to = today_num+extrapolate # max(GISAID_sel$DATE_NUM, na.rm=T)

# plot predicted values


# multinomial model predictions by province (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=as.numeric(seq(date.from, date.to)),
                            country=unique(data_agbyweekcountry1_subs$country))) # unique(GISAID_sel$country)
predgrid$continent = data_agbyweekcountry1_subs$continent[match(predgrid$country,
                                                                data_agbyweekcountry1_subs$country)]
# note: CIs are backtransformed from a logit scale as in the Effects package
# also possible (and slightly better) to use isometric logratio scale type="ilr", 
# but more hassle to backtransform
system.time(fit_preds <- data.frame(predictions(fit_global_multi, 
                       newdata = predgrid,
                       type = "logit",
                       vcov = fit_global_multi$vcov) %>%
             transform(conf.low = predicted - 1.96 * std.error,
                       conf.high = predicted + 1.96 * std.error) %>%
             group_by(rowid) |>
             mutate_at(c("predicted", "conf.low", "conf.high"), function (x) plogis(x)))) # 238s
# replace NAs by 0
fit_preds <- fit_preds %>% mutate(predicted = ifelse(is.na(predicted), 0, predicted),
                                  conf.low = ifelse(is.na(conf.low), 0, conf.low),
                                  conf.high = ifelse(is.na(conf.high), 0, conf.high))
fit_preds$date = as.Date(fit_preds$DATE_NUM, origin="1970-01-01")
fit_preds$variant = fit_preds$group
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)

# TO DO: check why some predictions before jan 2021 come out as NA,
# presumably due to some overflow??

fit_preds$country = factor(fit_preds$country) # , levels=levels_country
levels(fit_preds$country)
fit_preds$continent = factor(fit_preds$continent)

write_csv(fit_preds, file=file.path("plots", plotdir, "GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage.csv"))

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
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1*, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&
                                          data_agbyweekcountry1$country %in% sel_countries_min5,],
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
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1*, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-11)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-11)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant %in% sel_variants&
                                          data_agbyweekcountry1$country %in% sel_countries_min5,],
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
       subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1*, BA.2.30.20 or BA.2.75.2 sequences shown")) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01")+5, 
           ymin=0, ymax=1, alpha=0.5, fill="white") + # extrapolated part
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
muller_fit

ggsave(file=file.path("plots", plotdir, "global multinom fit_all data_predictions_muller plot.png"), width=20, height=12)




# map variant share onto case numbers

# sel_countries_top3 = c("India", "Nepal", "Singapore")
sel_countries = c("Austria","Bangladesh","Belgium","Denmark","France","Germany","Israel","Singapore","United Kingdom")
country_data = get_national_data(countries=sel_countries, # sel_countries_top3,
                                 source="who",
                                 level=1)
country_data$country[country_data$country=="United States"] = "USA"
country_data$country = factor(country_data$country, levels=levels(fit_preds$country))
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

  
fit_preds[(fit_preds$date==today)&(fit_preds$variant=="Omicron (BA.2.75)"),]

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
fit_preds_sel = fit_preds_sel[fit_preds_sel$country %in% sel_countries,]
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
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1*, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  scale_colour_manual("lineage", values=lineage_cols_plot) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2021-01-01"),NA)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) 
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=file.path("plots", plotdir,"global multinom fit_all data_predictions_confirmed cases multinomial fit by country.png"), width=20, height=12)


fit_preds3 = fit_preds2
# fit_preds3$cases[fit_preds3$date<=as.Date("2022-03-01")] = 0

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
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >=5 BQ.1*, BA.2.3.20 or BA.2.75.2 sequences shown")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20)) +
  coord_cartesian(ylim=c(1,NA), xlim=c(as.Date("2021-01-01"),NA))

ggsave(file=file.path("plots", plotdir,"\\global multinom fit_all data_predictions_confirmed cases stacked area multinomial fit by country.png"), width=20, height=12)

# TO DO : get IHME infection estimates from links below,
# convolve these to cases & map them onto variant frequencies
# https://www.healthdata.org/covid/data-downloads
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2020.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2021.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2022.csv
# get mortality data from Eurostat with Eurostat package, convolve
# case data by variant to mortality data & figure out
# death toll of each variant