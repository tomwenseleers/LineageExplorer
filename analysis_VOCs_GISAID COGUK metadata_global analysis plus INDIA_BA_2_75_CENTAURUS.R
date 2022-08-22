# ANALYSIS OF GROWTH ADVANTAGE OF BA.2.75, AKA CENTAURUS
# T. Wenseleers
# last update 21 AUGUST 2022

# set GISAID credentials ####
# set them first using 
# Sys.setenv(GISAIDR_USERNAME = "XXXXX")
# Sys.setenv(GISAIDR_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_USERNAME = "XXXXX") # not needed for this script
# Sys.setenv(GISAIDJSON_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_STREAM = "XXXXX")

# devtools::install_github("Wytamma/GISAIDR")
library(GISAIDR)
# setwd("~/Github/newcovid_belgium")
source(".//set_GISAID_credentials.R") # set GISAID credentials
source(".//download_GISAID.R") # load function to download GISAID metadata download package (lacking records from last few days)
source(".//download_GISAID_records.R") # load functions to download most recent GISAID records
source(".//download_COGUK.R") # load function to download COG-UK metadata
use_coguk = TRUE # use COG-UK data instead of GISAID data for UK?

# load required packages
library(nnet)
library(splines)
library(devtools)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
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
#                  repos = "https://epiforecasts.r-universe.dev"
# )
library(covidregionaldata)
library(tidyquant)
library(data.table)
library(R.utils)
library(locatexec)


# 1. LOAD DATA ####

today = as.Date(Sys.time())
today_num = as.numeric(today)
today # "2021-08-13"
plotdir = "GISAID_BA_2_75"
suppressWarnings(dir.create(file.path("plots",plotdir)))
tag = paste("@TWenseleers\n",today)

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2023-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

# import GISAID metadata ####

# download latest GISAID metadata ####
target_dir = "C:/Users/bherr/OneDrive - KU Leuven/Documents/Github/newcovid_belgium/data/GISAID" # target download directory
system.time(GISAID <- download_GISAD_meta(target_dir = target_dir,
                                          headless = FALSE,
                                          usr = Sys.getenv("GISAIDR_USERNAME"),
                                          psw = Sys.getenv("GISAIDR_PASSWORD"))) # 187s
download = tail(list.files(target_dir, pattern=".tar.xz"), 1)
download # "metadata_tsv_2022_08_13.tar.xz"  

# records go up to submission date
GISAID_max_submdate = as.Date(max(GISAID$submission_date, na.rm=T))
GISAID_max_submdate # "2022-08-11"
# GISAID_max_submdate = as.Date("2022-08-11")
nrow(GISAID) # 12513323




# add some extra manually downloaded data from the last few days ####

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
length(recent_records) # 155853  records
recent_records = recent_records[!recent_records %in% GISAID$accession_id]
length(recent_records) # 119080 record not available in GISAID download package

# dataframe with recently submitted records that are not yet in GISAID metadata package download
d_extra = download_GISAID_records(accession_ids = recent_records,
                                  get_sequence=FALSE, 
                                  clean_up=FALSE,
                                  target_dir="C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID/BA_2_75/extra",
                                  max_batch_size=10000, # maximum batch size, usually either 10000 or 5000
                                  headless = FALSE,
                                  usr=Sys.getenv("GISAIDR_USERNAME"),
                                  psw=Sys.getenv("GISAIDR_PASSWORD"))
dim(d_extra) # 127666                       17

# merge GISAID download package & recently submitted records
GISAID = dplyr::bind_rows(GISAID, d_extra)
nrow(GISAID) # 12640989

# with GISAIDR functions this should normally have worked, but since AA substitutions field is missing
# I am using my own function above (this is being fixed in GISAIDR)
# # function to split vector in chunks of max size chunk_length
# chunk = function(x, chunk_length=4000) split(x, ceiling(seq_along(x)/chunk_length))
# chunks = chunk(df_recent$accession_id)
# downloads = do.call(rbind, lapply(1:length(chunks),
#                    function (chunk) {
#                      message(paste0("Downloading batch ", chunk, " out of ", length(chunks)))
#                      Sys.sleep(3)
#                      return(download(credentials = credentials, 
#                               list_of_accession_ids = chunks[[chunk]],
#                               clean_up=FALSE)) } ))
# dim(downloads)
# names(downloads)

  
# # search in which countries at least 1 BA.2.75 was picked up in the last few days
# df_countries_BA_2_75 = query(
#   credentials = credentials, 
#   lineage = "BA.2.75", 
#   from_subm = as.character(GISAID_max_submdate), 
#   to_subm = as.character(today),
#   fast = TRUE
# )
# full_df_countries_BA_2_75 = download(credentials = credentials, 
#                            list_of_accession_ids = df_countries_BA_2_75$accession_id)
# full_df_countries_BA_2_75$sequence = NULL
# countries_BA_2_75 = unique(full_df_countries_BA_2_75$country)
# countries_BA_2_75



# parse location field ####
# parse continent / country / location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID$location), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
# unique(loc[,1]) # continent
# unique(loc[,2]) # country
# unique(loc[,3]) # city or province/state

levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
GISAID$continent = factor(loc[,1], levels=levels_continents)
GISAID$country = factor(loc[,2])
levels_countries = levels(GISAID$country)
GISAID$location = factor(loc[,3])
levels_locations = levels(GISAID$location)


# LOAD COG-UK DATA FOR THE UNITED KINGDOM ####
if (use_coguk) { coguk = download_COGUK_meta()
  # MERGE GISAID (MINUS UK GISAID DATA) & COG-UK DATA FOR UK
  GISAID = dplyr::bind_rows(GISAID[GISAID$country!="United Kingdom",], 
                           coguk)
  }
nrow(GISAID) # 12562374

# PARSE GISAID DATA ####
GISAID = as.data.frame(GISAID) 
# TO DO: keep as more efficient data.table or as tibble
# & make sure everything below still works

# parse date & check dates are valid
# records with valid date
GISAID$date_isvalid = (str_count(GISAID$collection_date,
                                 pattern = "-")==2)
GISAID$date = as.Date(NA)
GISAID$date[which(GISAID$date_isvalid)] = as.Date(fast_strptime(GISAID$collection_date[which(GISAID$date_isvalid)], "%Y-%m-%d")) # faster than as.Date(GISAID$collection_date)

# define variant lineages and colours ####
sel_target_VOC = "Omicron (BA.2.75)"
sel_reference_VOC = "Omicron (BA.5)"
levels_VARIANTS = c(sel_reference_VOC, "B.1.177 (EU1)", "B.1.160 (EU2)", "B.1.221 (20A/S:98F)", "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.38)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.2.76)", "Omicron (BA.5.2)", "Omicron (BA.5.2.1.7)", sel_target_VOC)
levels_VARIANTS_plot = c("Other", "B.1.177 (EU1)", "B.1.160 (EU2)", "B.1.221 (20A/S:98F)", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.38)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.5)", "Omicron (BA.5.2)", "Omicron (BA.5.2.1.7)", "Omicron (BA.2.76)", "Omicron (BA.2.75)")

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
  levels_VARIANTS_plot=="Omicron (BA.5.2.1.7)" ~ "cyan3",
  levels_VARIANTS_plot=="Omicron (BA.2.76)" ~ "magenta4",
  levels_VARIANTS_plot=="Omicron (BA.2.75)" ~ "magenta",
)
# library(pals)
# pal.bands(lineage_cols_plot)
# pal.volcano(lineage_cols_plot)
# pal.zcurve(lineage_cols_plot)
lineage_cols = lineage_cols_plot[match(levels_VARIANTS,levels_VARIANTS_plot)]


# CODE MAIN LINEAGES ####
GISAID$aa_substitutions[is.na(GISAID$aa_substitutions)] = ""

# PS this is slow - optimize this, maybe use multiple cores / multidplyr?
# or avoid having to use all these grepl patterns

GISAID$variant = case_when(
  (grepl("BA.2.75", GISAID$pango_lineage, fixed=T)|(grepl("NSP3_S403L",GISAID$aa_substitutions)& 
                                                      grepl("NSP8_N118S",GISAID$aa_substitutions))|
     (grepl("NSP3_S403L",GISAID$aa_substitutions)& 
        grepl("E_T11A",GISAID$aa_substitutions))) ~ "Omicron (BA.2.75)",
   (grepl("BA.2.76", GISAID$pango_lineage, fixed=T)|(grepl("Spike_Y248N",GISAID$aa_substitutions)&
      grepl("Spike_R346T",GISAID$aa_substitutions))) ~ "Omicron (BA.2.76)", 
  grepl("^BA\\.4\\.6$", GISAID$pango_lineage) ~ "Omicron (BA.4.6)",
  grepl("^BA\\.5\\.2\\.1$|^BA\\.5\\.2\\.1|^BF\\.7", GISAID$pango_lineage)&
    grepl("Spike_R346T",GISAID$aa_substitutions)&
    (GISAID$date>as.Date("2022-02-01")) ~ "Omicron (BA.5.2.1.7)",
  (grepl("^BA\\.5\\.2$|^BA\\.5\\.2", GISAID$pango_lineage)|
     (grepl("BF", GISAID$pango_lineage)) )&
    (GISAID$date>as.Date("2022-02-01")) ~ "Omicron (BA.5.2)",
  grepl("B.1.617.2", GISAID$pango_lineage, fixed=T) | grepl("AY", GISAID$pango_lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", GISAID$pango_lineage) ~ "Alpha",
  grepl("B.1.351", GISAID$pango_lineage, fixed=T) ~ "Beta",
  (grepl("^BA\\.1$|BA\\.1\\.", GISAID$pango_lineage)&(GISAID$date>as.Date("2021-11-01"))) ~ "Omicron (BA.1)",
  (((grepl("^BA\\.4",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
                                                  ((grepl("L452R", GISAID$aa_substitutions)&
                                                    grepl("486V", GISAID$aa_substitutions)&
                                                    grepl("11F", GISAID$aa_substitutions)&
                                                    (!grepl("D3N",GISAID$aa_substitutions)) )))) ~ "Omicron (BA.4)",
  (grepl("^BA\\.5",GISAID$pango_lineage)|grepl("BE|BF",GISAID$pango_lineage)|((GISAID$pango_lineage!="Unassigned")&grepl("M_D3N",GISAID$aa_substitutions)))&(GISAID$date>as.Date("2022-02-01")) ~ "Omicron (BA.5)", # cf pattern used by Alex Selby
  (grepl("^BA\\.2\\.12\\.1$", GISAID$pango_lineage))&(GISAID$date>as.Date("2022-02-01")) ~ "Omicron (BA.2.12.1)",
  (grepl("^BA\\.2\\.38$|^BA\\.2\\.38\\.", GISAID$pango_lineage)) ~ "Omicron (BA.2.38)",
  (grepl("^BA\\.2",GISAID$pango_lineage)) ~ "Omicron (BA.2)",
  (grepl("^B\\.1\\.177$|^B\\.1\\.177\\.", GISAID$pango_lineage)) ~ "B.1.177 (EU1)",
  (grepl("^B\\.1\\.160$|^B\\.1\\.160\\.", GISAID$pango_lineage)) ~ "B.1.160 (EU2)",
  (grepl("^B\\.1\\.221$|^B\\.1\\.221\\.", GISAID$pango_lineage)) ~ "B.1.221 (20A/S:98F)",
  GISAID$pango_lineage!="Unassigned" ~ "Other" # assigns NA to remaining Unassigned & remove them later on
  # TRUE ~ "Other" # to assign Unassigned lineages to category Other
)

# note: in India BA.2.38 & BA.2.38.1 caused an infection wave in some states - hence separated out above
# B.1.177+B.1.160+B.1.221 were behind the 2020 wave in fall in Europe & each had one spike mutations & a small growth rate advantage relative to predominant B.1.1
# BA.5.2.1.7 = BF.7

# variants to watch & maybe add in due time:
# https://cov-spectrum.org/collections/1
# BJ.1, https://github.com/cov-lineages/pango-designation/issues/915
# (just 11 seqs for now, https://cov-spectrum.org/explore/World/AllSamples/Past6M/variants?variantQuery=%5B4-of%3A+ORF1a%3A47R%2C+S%3A83A%2C+S%3A146Q%2C+S%3A213E%2C+S%3A339H%2C+S%3A445P%2C+S%3A483A%2C+S%3A1003I%2C+M%3A3Y%2C+ORF7a%3A110T%2C+N%3A282I%2C+15738T%2C+15939C%5D&)

sum(is.na(GISAID$variant)) # 197922 unassigned
sum(GISAID$pango_lineage=="Unassigned", na.rm=T) # 197996 originally unassigned in GISAID
sum(GISAID$variant=="Omicron (BA.2.75)", na.rm=T) # 3852 BA.2.75
sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$country=="India", na.rm=T) # 3036 BA.2.75 for India
sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$date_isvalid, na.rm=T) # 3307 BA.2.75 with valid date
sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$country=="India"&GISAID$date_isvalid, na.rm=T) # 2493 BA.2.75 for India with valid date
sum(GISAID$pango_lineage=="BA.2.75", na.rm=T) # 3144
sum(GISAID$variant=="Omicron (BA.5)", na.rm=T) # 327070
sum(GISAID$variant=="Omicron (BA.5.2.1.7)", na.rm=T) # 2089
table(GISAID$variant)

# GISAID = GISAID[!is.na(GISAID$variant),]
nrow(GISAID) # 12562374



# GISAID SELECTION ####

GISAID_sel = GISAID

# remove records with invalid/incomplete dates ####
GISAID_sel = GISAID_sel[GISAID_sel$date_isvalid,]
GISAID_sel = GISAID_sel[which(GISAID_sel$host=="Human"),]
nrow(GISAID_sel) # 12343167

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

nBA_2_75 = sum(GISAID_sel$variant=="Omicron (BA.2.75)", na.rm=T)
nBA_2_75 # 3293 BA.2.75 records so far with valid date

nBA_2_75_india = sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="India", na.rm=T)
nBA_2_75_india # 2493 BA.2.75 records so far for India with valid date

sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="United Kingdom", na.rm=T)
# 59 BA.2.75 in UK so far
maxsubmdate = today


# these countries are used as a selection for plotting later on
tab = as.data.frame(table(GISAID_sel$country, GISAID_sel$variant))

sel_countries_min10 = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&tab$Freq>=10,"Var1"])
sel_countries_min10
# [1] "Australia"      "Austria"        "Canada"         "Denmark"       
# [5] "Germany"        "India"          "Israel"         "Japan"         
# [9] "Nepal"          "Russia"         "Singapore"      "South Korea"   
# [13] "United Kingdom" "USA"   

sel_countries = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&
                                   tab$Freq>1,"Var1"])
sel_countries



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
set.seed(1)
# best model:
fit_global_multi = nnet::multinom(variant ~ 
                                    ns(DATE_NUM, df=2)+
                                    ns(DATE_NUM, df=2):continent+
                                    country, 
                                  weights=count, 
                                  data=data_agbyweekcountry1, 
                                  maxit=10000, MaxNWts=100000)
# BIC(fit_global_multi) 
# saveRDS(fit_global_multi, file="./fits/GISAID_global_fit.rds")
# fit_global_multi = readRDS(file="./fits/GISAID_global_fit.rds")

# # calculate current pairwise growth rate differences
# emtr_pairw = emtrends(bestfit_multi, revpairwise ~ variant, by="continent", 
#                 var="DATE_NUM",  mode="latent",
#                 at=list(DATE_NUM=max(GISAID_sel$DATE_NUM, na.rm=T)))
# delta_r_pairw = data.frame(confint(emtr_pairw, 
#                              adjust="none", df=NA)$contrasts, 
#                      p.value=as.data.frame(emtr_pairw$contrasts)$p.value)
# delta_r_pairw
# write.csv(delta_r_pairw, file.path("plots", plotdir, "growth advantage BA.2.75 vs all other strains by continent.csv"), row.names=F)
# 
# # plot growth rate advantage of BA.2.75
# delta_r_pairw2 = delta_r_pairw[delta_r_pairw$contrast %in%
#                                  c("Omicron (BA.2.75) - Omicron (BA.2)",
#                                    "Omicron (BA.2.75) - Omicron (BA.2.38)",
#                                    "Omicron (BA.2.75) - Omicron (BA.2.12.1)",
#                                    "Omicron (BA.2.75) - Omicron (BA.5)",
#                                    "Omicron (BA.2.75) - Omicron (BA.5.2)",
#                                    "Omicron (BA.2.75) - Omicron (BA.4)",
#                                    "Omicron (BA.2.75) - Omicron (BA.4.6)"),]
# delta_r_pairw2$contrast = droplevels(delta_r_pairw2$contrast)
# delta_r_pairw2$contrast = factor(delta_r_pairw2$contrast,
#                                  levels=c("Omicron (BA.2.75) - Omicron (BA.2)",
#                                           "Omicron (BA.2.75) - Omicron (BA.2.12.1)",
#                                           "Omicron (BA.2.75) - Omicron (BA.2.38)",
#                                           "Omicron (BA.2.75) - Omicron (BA.4)",
#                                           "Omicron (BA.2.75) - Omicron (BA.5)",
#                                           "Omicron (BA.2.75) - Omicron (BA.5.2)",
#                                           "Omicron (BA.2.75) - Omicron (BA.4.6)"))
# delta_r_pairw2 = delta_r_pairw2[!delta_r_pairw2$continent %in%
#                                   c("South America", "Africa"),]
# 
# library(RColorBrewer)
# cols = brewer.pal(n = 8, "Blues")[-1]
# qplot(data=delta_r_pairw2,
#       x=continent, group=contrast,
#       y=estimate,
#       ymin=asymp.LCL,
#       ymax=asymp.UCL,
#       geom="blank") +
#   geom_col(position=position_dodge(), aes(fill=contrast)) +
#   geom_linerange(position=position_dodge(.9)) +
#   labs(tag = tag) +
#   theme(plot.tag.position = "bottomright",
#         plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
#   scale_fill_manual(name="", values=cols) +
#   # scale_fill_manual("", values=c(lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.2.12.1)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.4.6)")])) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   xlab("") +
#   ylab("Difference in growth rate (/day)") +
#   ggtitle("GROWTH RATE ADVANTAGE OF BA.2.75 SARS-CoV2 VARIANT",
#           subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)*continent+country,\nusing GISAID data with submission data up to ", maxsubmdate, " &\ndata from countries with at least 1 BA.2.75 sequence plus COG-UK data\n(n=", nBA_2_75, " BA.2.75 sequences in total)"))
# 
# ggsave(file=file.path("plots", plotdir,"growth rate advantage BA_2_75_by continent.png"), width=7, height=5)
# 
# # plot growth rate advantage of BA.4.6
# emtr_pairw_BA_4_6 = emtrends(bestfit_multi, trt.vs.ctrl ~ variant, 
#                        ref = which(bestfit_multi$lab=="Omicron (BA.4.6)"),
#                        reverse = TRUE,
#                        by="continent", 
#                       var="DATE_NUM",  mode="latent",
#                       at=list(DATE_NUM=max(GISAID_sel$DATE_NUM, na.rm=T)))
# delta_r_BA_4_6 = data.frame(confint(emtr_pairw_BA_4_6, 
#                                    adjust="none", df=NA)$contrasts, 
#                            p.value=as.data.frame(emtr_pairw_BA_4_6$contrasts)$p.value)
# 
# delta_r_BA_4_6_2 = delta_r_BA_4_6[delta_r_BA_4_6$contrast %in% 
#                                     c("Omicron (BA.4.6) - Omicron (BA.5)",
#                                       "Omicron (BA.4.6) - Omicron (BA.5.2)",
#                                       "Omicron (BA.4.6) - Omicron (BA.4)",
#                                       "Omicron (BA.4.6) - Omicron (BA.2.75)",
#                                       "Omicron (BA.4.6) - Omicron (BA.2.12.1)",
#                                       "Omicron (BA.4.6) - Omicron (BA.2.38)",
#                                       "Omicron (BA.4.6) - Omicron (BA.2)"
#                                       ),]
# delta_r_BA_4_6_2$contrast = droplevels(delta_r_BA_4_6_2$contrast)
# delta_r_BA_4_6_2$contrast = factor(delta_r_BA_4_6_2$contrast,
#                                  levels=c("Omicron (BA.4.6) - Omicron (BA.2)",
#                                           "Omicron (BA.4.6) - Omicron (BA.2.12.1)",
#                                           "Omicron (BA.4.6) - Omicron (BA.2.38)",
#                                           "Omicron (BA.4.6) - Omicron (BA.4)",
#                                           "Omicron (BA.4.6) - Omicron (BA.5)",
#                                           "Omicron (BA.4.6) - Omicron (BA.5.2)",
#                                           "Omicron (BA.4.6) - Omicron (BA.2.75)"))
# delta_r_BA_4_6_2 = delta_r_BA_4_6_2[!delta_r_BA_4_6_2$continent %in%
#                                   c("South America", "Africa"),]
# 
# library(RColorBrewer)
# cols = brewer.pal(n = 8, "Blues")[-1]
# qplot(data=delta_r_BA_4_6_2,
#       x=continent, group=contrast,
#       y=estimate,
#       ymin=asymp.LCL,
#       ymax=asymp.UCL,
#       geom="blank") +
#   geom_col(position=position_dodge(), aes(fill=contrast)) +
#   geom_linerange(position=position_dodge(.9)) +
#   labs(tag = tag) +
#   theme(plot.tag.position = "bottomright",
#         plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
#   scale_fill_manual(name="", values=cols) +
#   # scale_fill_manual("", values=c(lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.2.12.1)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
#   #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.4.6)")])) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   xlab("") +
#   ylab("Difference in growth rate (/day)") +
#   ggtitle("GROWTH RATE ADVANTAGE OF BA.4.6 SARS-CoV2 VARIANT",
#           subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)*continent+country,\nusing GISAID data with submission data up to ", maxsubmdate, " &\ndata from countries with at least 1 BA.2.75 sequence plus COG-UK data\n(n=", nBA_2_75, " BA.2.75 sequences in total)"))
# 
# ggsave(file=file.path("plots", plotdir,
#                       "growth rate advantage BA_4_6_by continent.png"), width=7, height=5)




# PLOT MULTINOMIAL FIT

extrapolate = 150
date.from = as.numeric(as.Date("2020-01-01"))
date.to = today_num+extrapolate # max(GISAID_sel$DATE_NUM, na.rm=T)

# plot predicted values
# multinomial model predictions by province (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to),
                            country=unique(GISAID_sel$country)))
predgrid$continent = GISAID_sel$continent[match(predgrid$country,
                                                        GISAID_sel$country)]

fit_preds = data.frame(predgrid, as.data.frame(predict(fit_global_multi, 
                                                       newdata=predgrid, type="prob")),check.names=F)
fit_preds = gather(fit_preds, variant, prob, all_of(levels_VARIANTS), factor_key=TRUE)
fit_preds$date = as.Date(fit_preds$DATE_NUM, origin="1970-01-01")
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)
levels_country = rev(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","country"][order(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","prob"])])
as.character(levels_country)
levels_country = c("India", "Nepal", "China", 
                   as.character(levels_country)[!as.character(levels_country) %in% c("India","Nepal","China")])
fit_preds$country = factor(fit_preds$country, levels=levels_country)
levels(fit_preds$country)
fit_preds$continent = factor(fit_preds$continent)
write_csv(fit_preds, file=file.path("plots", plotdir, "GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage.csv"))




# PLOT MULTINOMIAL FIT ON LOGIT SCALE ####
ymin = 0.0001
ymax = 0.999
if (!is.null(fit_preds$asymp.LCL)) fit_preds$asymp.LCL[fit_preds$asymp.LCL<ymin] = ymin
if (!is.null(fit_preds$asymp.UCL)) fit_preds$asymp.UCL[fit_preds$asymp.UCL<ymin] = ymin
if (!is.null(fit_preds$asymp.UCL)) fit_preds$asymp.UCL[fit_preds$asymp.UCL>ymax] = ymax
fit_preds$prob[fit_preds$prob<ymin] = ymin
fit_preds$prob[fit_preds$prob>ymax] = ymax

data_agbyweekcountry1$variant = factor(data_agbyweekcountry1$variant, levels=levels_VARIANTS_plot)
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)

plot_preds_logit = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$country %in% sel_countries,], 
                         x=date, y=prob, geom="blank") +
  facet_wrap(~ country) +
  # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
  #                 fill=variant
  # ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >1 BA.2.75 sequences shown (n=", nBA_2_75," BA.2.75 sequences in total)")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$country %in% sel_countries,],
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
  coord_cartesian(xlim=c(as.Date("2020-05-01"),NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
plot_preds_logit

ggsave(file=file.path("plots", plotdir, "global multinom fit_all data_predictions_logit scale.png"), width=20, height=12)

# plot predicted values as Muller plot
muller_fit = ggplot(data=fit_preds[fit_preds$country %in% sel_countries&
                                     fit_preds$date>=as.Date("2020-05-01"),], 
                    aes(x=date, y=prob, group=variant)) +
  facet_wrap(~ country, ncol=4) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES",
       subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >1 BA.2.75 sequences shown (n=", nBA_2_75," BA.2.75 sequences in total)")) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01")+5, 
           ymin=0, ymax=1, alpha=0.5, fill="white") + # extrapolated part
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))
muller_fit

ggsave(file=file.path("plots", plotdir, "global multinom fit_all data_predictions_muller plot.png"), width=20, height=12)




# map variant share onto case numbers

# sel_countries_top3 = c("India", "Nepal", "Singapore")
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY IN\nCOUNTRIES WITH >1 SEQUENCED BA.2.75 CASES",
          subtitle="7 day simple moving average, using\nWHO case data accessed via covidregionaldata package")
# ggsave(file=file.path("plots", plotdir, "new cases countries with more than 10 sequenced BA_2_75 cases.png"), width=7, height=5)

  
fit_preds[(fit_preds$date==today)&(fit_preds$variant=="Omicron (BA.2.75)"),]

# fit_multi_predsbycountry = data.frame(emmeans(fit8_multi,
#                                                   ~ variant,
#                                                   by=c("DATE_NUM",
#                                                        "continent",
#                                                        "country"),
#                                                   at=list(DATE_NUM=today_num), # by=7 just to speed up things a bit
#                                                   mode="prob", df=NA,
#                                                   rg.limit=100000))
# 
# fit_multi_predsbycountry[fit_multi_predsbycountry$variant=="Omicron (BA.2.75)",]


fit_preds$totnewcases = 
  country_data$cases_new[match(interaction(fit_preds$country,
                                          fit_preds$date),
                                      interaction(country_data$country, 
                                                  country_data$date))]
library(zoo)
fit_preds = fit_preds %>% 
  group_by(country) %>% 
  mutate(totnewcases_smoothed = rollmean(totnewcases, 7, na.pad = T))
fit_preds$cases = fit_preds$totnewcases_smoothed*fit_preds$prob
fit_preds$cases[fit_preds$cases==0] = NA
fit_preds$cases[fit_preds$prob<0.001] = NA

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
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >1 BA.2.75 sequences shown (n=", nBA_2_75," BA.2.75 sequences in total)")) +
  scale_colour_manual("lineage", values=lineage_cols_plot) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA)) +
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
          subtitle=paste0("case data accessed via the covidregionaldata package\nlineage frequencies based on GISAID data up to ",today, " plus COG-UK data\nand multinomial spline fit variant ~ ns(date, df=2)+ns(date, df=2):continent+country,\nall countries with >1 BA.2.75 sequences shown (n=", nBA_2_75," BA.2.75 sequences in total)")) +
  # coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))
  theme(plot.title=element_text(size=25)) +
  theme(plot.subtitle=element_text(size=20))

ggsave(file=file.path("plots", plotdir,"\\global multinom fit_all data_predictions_confirmed cases stacked area multinomial fit by country.png"), width=20, height=12)

# TO DO : get IHME infection estimates from links below,
# convolute these to cases & map them onto variant frequencies
# https://www.healthdata.org/covid/data-downloads
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2020.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2021.csv
# https://ihmecovid19storage.blob.core.windows.net/archive/2022-07-19/data_download_file_reference_2022.csv




# 3. ANALYSIS OF BA.2.75 GROWTH RATE ADVANTAGE (CODE NOT FINISHED/WRITTEN YET) ####
# (here just using countries with at least 1 sequenced BA.2.75 case)

# TO DO: replace this with more elegant/tidy dplyr coce? 

# we first subset to countries where at least 1 BA.2.75 was detected
tab = as.data.frame(table(GISAID_sel$country, GISAID_sel$variant))

sel_countries_min10 = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&tab$Freq>=10,"Var1"])
sel_countries_min10
# [1] "Australia"      "Austria"        "Canada"         "Denmark"        "Germany"        "India"         
# [7] "Israel"         "Japan"          "Nepal"          "Russia"         "Singapore"      "United Kingdom"
# [13] "USA" 

sel_countries = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&tab$Freq>=1,"Var1"])
sel_countries

# GISAID selection : subsetted to countries with >=1 BA.2.75 sequence ####
GISAID_sel_BA_2_75 = GISAID_sel[as.character(GISAID_sel$country) %in% sel_countries,]
GISAID_sel_BA_2_75 = GISAID_sel_BA_2_75[!is.na(GISAID_sel_BA_2_75$variant),] # we remove records with unassigned lineage
GISAID_sel_BA_2_75$country = factor(GISAID_sel_BA_2_75$country)
table(GISAID_sel_BA_2_75$country, GISAID_sel_BA_2_75$variant)
table(GISAID_sel_BA_2_75$continent, GISAID_sel_BA_2_75$variant)

df_cont=as.data.frame(table(GISAID_sel_BA_2_75[GISAID_sel_BA_2_75$variant=="Omicron (BA.2.75)","continent"]))
colnames(df_cont)=c("continent","BA.2.75")
df_cont[order(df_cont$BA.2.75,decreasing=T),]

df=as.data.frame(table(GISAID_sel_BA_2_75[GISAID_sel_BA_2_75$variant=="Omicron (BA.2.75)","location"]))
df=df[df$Freq!=0,]
df=df[order(df$Freq, decreasing=T),]
df

df=as.data.frame(table(GISAID_sel_BA_2_75[GISAID_sel_BA_2_75$variant=="Omicron (BA.2.75)","country"]))
df=df[df$Freq!=0,]
df=df[order(df$Freq, decreasing=T),]
colnames(df)=c("country","BA.2.75")
df


# AGGREGATE DATA BY DATE & COUNTRY ####
data_agbydatecountry1 = as.data.frame(table(GISAID_sel_BA_2_75$date, GISAID_sel_BA_2_75$country, GISAID_sel_BA_2_75$variant))
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
data_agbydatecountry1$continent = GISAID_sel_BA_2_75$continent[match(data_agbydatecountry1$country, GISAID_sel_BA_2_75$country)]
data_agbydatecountry1$continent = factor(data_agbydatecountry1$continent)
# write.csv(data_agbydatecountry1, file=".//data//GISAID//GISAID aggregated counts by date and lineage.csv", row.names=F)

# AGGREGATE DATA BY WEEK & COUNTRY ####
data_agbyweekcountry1 = as.data.frame(table(GISAID_sel_BA_2_75$floor_date, GISAID_sel_BA_2_75$country, GISAID_sel_BA_2_75$variant))
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
if (is.null(data_agbyweekcountry1$continent)) data_agbyweekcountry1$continent = GISAID_sel_BA_2_75$continent[match(data_agbyweekcountry1$country, GISAID_sel_BA_2_75$country)]
data_agbyweekcountry1$continent = factor(data_agbyweekcountry1$continent)
data_agbyweekcountry1$collection_date = as.Date(data_agbyweekcountry1$collection_date)
data_agbyweekcountry1$variant = factor(data_agbyweekcountry1$variant, levels=levels_VARIANTS)
# write.csv(data_agbyweekcountry1, file="./data/GISAID/GISAID aggregated counts by start of week and lineage.csv", row.names=F)

# fit multinomial spline model ####
set.seed(1)
#fit0_multi = nnet::multinom(variant ~ scale(DATE_NUM), weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit1_multi = nnet::multinom(variant ~ scale(DATE_NUM)+continent, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit2_multi = nnet::multinom(variant ~ scale(DATE_NUM)+country, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit3_multi = nnet::multinom(variant ~ scale(DATE_NUM)*continent, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit4_multi = nnet::multinom(variant ~ scale(DATE_NUM)*continent+country, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit5_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+continent, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit6_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+country, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit7_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*continent, weights=count, data=data_agbyweekcountry1, maxit=1000)

# best model:
fit8_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*continent+country, weights=count, data=data_agbyweekcountry1, maxit=10000)
# fit8b_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+ns(DATE_NUM, df=2):continent+country, weights=count, data=data_agbyweekcountry1, maxit=10000)

#fit10_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*country, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit11_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2), weights=count, data=data_agbyweekcountry1, maxit=1000)
# BIC(fit8_multi, fit8b_multi) # fit8_multi has best BIC

bestfit_multi = fit8_multi # I will use this model

# PS model.matrix(fit8_multi) gives the model matrix

# tests with some other multinomial fitting functions
# library(mclogit)
# fit8b_mblogit = mblogit(variant ~ ns(DATE_NUM, df=2)+ns(DATE_NUM, df=2):continent+country, weights=count, estimator="REML", dispersion=FALSE, from.table=TRUE, data=data_agbyweekcountry1, maxit=10000)
# fit9b_mblogit = mblogit(variant ~ ns(DATE_NUM, df=3)+ns(DATE_NUM, df=3):continent+country, weights=count, estimator="REML", dispersion=FALSE, from.table=TRUE, data=data_agbyweekcountry1, maxit=10000)
# library(mlogit)
# mldata = mlogit.data(data_agbyweekcountry1, choice="variant", shape="wide")
# fit8_mlogit = mlogit(variant ~ 0 | ns(DATE_NUM, df=2)*continent+country, weights=count,
#                        data=mldata)
# library(mnlogit)
# fit8_mnlogit = mnlogit(variant ~ 0 | ns(DATE_NUM, df=2)*continent+country, weights=count,
#                        data=mldata, ncores=15)
# library(gmnl)
# check

# growth rate advantage compared to BA.5
# emtr = emtrends(bestfit_multi, trt.vs.ctrl ~ variant, by="continent", 
#                      var="DATE_NUM",  mode="latent",
#                      at=list(DATE_NUM=max(GISAID_sel$DATE_NUM, na.rm=T)))
# # TO DO: use faster https://github.com/vincentarelbundock/marginaleffects package?
# # but first need to make sure that nnet:multinom predict method
# # both supports type="link" and type="prob"/"response"
# delta_r = data.frame(confint(emtr, 
#                                    adjust="none", df=NA)$contrasts, 
#                            p.value=as.data.frame(emtr$contrasts)$p.value)
# delta_r
# write.csv(delta_r, file.path(plots,"growth advantage BA.2.75 vs BA.5 by continent.csv"), row.names=F)

emtr_pairw = emtrends(bestfit_multi, revpairwise ~ variant, by="continent", 
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=max(GISAID_sel_BA_2_75$DATE_NUM, na.rm=T)))
delta_r_pairw = data.frame(confint(emtr_pairw, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtr_pairw$contrasts)$p.value)
delta_r_pairw
write.csv(delta_r_pairw, file.path("plots", plotdir, "growth advantage BA.2.75 vs all other strains by continent.csv"), row.names=F)

# plot growth rate advantage of BA.2.75
delta_r_pairw2 = delta_r_pairw[delta_r_pairw$contrast %in%
                                 c("Omicron (BA.2.75) - Omicron (BA.2)",
                                   "Omicron (BA.2.75) - Omicron (BA.2.38)",
                                   "Omicron (BA.2.75) - Omicron (BA.2.12.1)",
                                   "Omicron (BA.2.75) - Omicron (BA.5)",
                                   "Omicron (BA.2.75) - Omicron (BA.5.2)",
                                   "Omicron (BA.2.75) - Omicron (BA.4)",
                                   "Omicron (BA.2.75) - Omicron (BA.4.6)"),]
delta_r_pairw2$contrast = droplevels(delta_r_pairw2$contrast)
delta_r_pairw2$contrast = factor(delta_r_pairw2$contrast,
                                 levels=c("Omicron (BA.2.75) - Omicron (BA.2)",
                                          "Omicron (BA.2.75) - Omicron (BA.2.12.1)",
                                          "Omicron (BA.2.75) - Omicron (BA.2.38)",
                                          "Omicron (BA.2.75) - Omicron (BA.4)",
                                          "Omicron (BA.2.75) - Omicron (BA.5)",
                                          "Omicron (BA.2.75) - Omicron (BA.5.2)",
                                          "Omicron (BA.2.75) - Omicron (BA.4.6)"))
delta_r_pairw2 = delta_r_pairw2[!delta_r_pairw2$continent %in%
                                  c("South America", "Africa"),]

library(RColorBrewer)
cols = brewer.pal(n = 8, "Blues")[-1]
qplot(data=delta_r_pairw2,
      x=continent, group=contrast,
      y=estimate,
      ymin=asymp.LCL,
      ymax=asymp.UCL,
      geom="blank") +
  geom_col(position=position_dodge(), aes(fill=contrast)) +
  geom_linerange(position=position_dodge(.9)) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  scale_fill_manual(name="", values=cols) +
  # scale_fill_manual("", values=c(lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.2.12.1)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.4.6)")])) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("Difference in growth rate (/day)") +
  ggtitle("GROWTH RATE ADVANTAGE OF BA.2.75 SARS-CoV2 VARIANT",
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)*continent+country,\nusing GISAID data with submission data up to ", maxsubmdate, " &\ndata from countries with at least 1 BA.2.75 sequence plus COG-UK data\n(n=", nBA_2_75, " BA.2.75 sequences in total)"))

ggsave(file=file.path("plots", plotdir,"growth rate advantage BA_2_75_by continent.png"), width=7, height=5)

# plot growth rate advantage of BA.4.6
emtr_pairw_BA_4_6 = emtrends(bestfit_multi, trt.vs.ctrl ~ variant, 
                             ref = which(bestfit_multi$lab=="Omicron (BA.4.6)"),
                             reverse = TRUE,
                             by="continent", 
                             var="DATE_NUM",  mode="latent",
                             at=list(DATE_NUM=max(GISAID_sel_BA_2_75$DATE_NUM, na.rm=T)))
delta_r_BA_4_6 = data.frame(confint(emtr_pairw_BA_4_6, 
                                    adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtr_pairw_BA_4_6$contrasts)$p.value)

delta_r_BA_4_6_2 = delta_r_BA_4_6[delta_r_BA_4_6$contrast %in% 
                                    c("Omicron (BA.4.6) - Omicron (BA.5)",
                                      "Omicron (BA.4.6) - Omicron (BA.5.2)",
                                      "Omicron (BA.4.6) - Omicron (BA.4)",
                                      "Omicron (BA.4.6) - Omicron (BA.2.75)",
                                      "Omicron (BA.4.6) - Omicron (BA.2.12.1)",
                                      "Omicron (BA.4.6) - Omicron (BA.2.38)",
                                      "Omicron (BA.4.6) - Omicron (BA.2)"
                                    ),]
delta_r_BA_4_6_2$contrast = droplevels(delta_r_BA_4_6_2$contrast)
delta_r_BA_4_6_2$contrast = factor(delta_r_BA_4_6_2$contrast,
                                   levels=c("Omicron (BA.4.6) - Omicron (BA.2)",
                                            "Omicron (BA.4.6) - Omicron (BA.2.12.1)",
                                            "Omicron (BA.4.6) - Omicron (BA.2.38)",
                                            "Omicron (BA.4.6) - Omicron (BA.4)",
                                            "Omicron (BA.4.6) - Omicron (BA.5)",
                                            "Omicron (BA.4.6) - Omicron (BA.5.2)",
                                            "Omicron (BA.4.6) - Omicron (BA.2.75)"))
delta_r_BA_4_6_2 = delta_r_BA_4_6_2[!delta_r_BA_4_6_2$continent %in%
                                      c("South America", "Africa"),]

library(RColorBrewer)
cols = brewer.pal(n = 8, "Blues")[-1]
qplot(data=delta_r_BA_4_6_2,
      x=continent, group=contrast,
      y=estimate,
      ymin=asymp.LCL,
      ymax=asymp.UCL,
      geom="blank") +
  geom_col(position=position_dodge(), aes(fill=contrast)) +
  geom_linerange(position=position_dodge(.9)) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  scale_fill_manual(name="", values=cols) +
  # scale_fill_manual("", values=c(lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.2.12.1)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.5.2)")],
  #                                lineage_cols_plot[which(levels_VARIANTS_plot=="Omicron (BA.4.6)")])) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("Difference in growth rate (/day)") +
  ggtitle("GROWTH RATE ADVANTAGE OF BA.4.6 SARS-CoV2 VARIANT",
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=2)*continent+country,\nusing GISAID data with submission data up to ", maxsubmdate, " &\ndata from countries with at least 1 BA.2.75 sequence plus COG-UK data\n(n=", nBA_2_75, " BA.2.75 sequences in total)"))

ggsave(file=file.path("plots", plotdir,
                      "growth rate advantage BA_4_6_by continent.png"), width=7, height=5)




# PLOT MULTINOMIAL FIT

extrapolate = 150
date.from = as.numeric(as.Date("2020-01-01"))
date.to = today_num+extrapolate # max(GISAID_sel_BA_2_75$DATE_NUM, na.rm=T)

# plot predicted values
# multinomial model predictions by province (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to),
                            country=levels(GISAID_sel$country)))
predgrid$continent = GISAID_sel_BA_2_75$continent[match(predgrid$country,
                                                        GISAID_sel_BA_2_75$country)]

fit_preds = data.frame(predgrid, as.data.frame(predict(bestfit_multi, 
                                                       newdata=predgrid, type="prob")),check.names=F)
fit_preds = gather(fit_preds, variant, prob, all_of(levels_VARIANTS), factor_key=TRUE)
fit_preds$date = as.Date(fit_preds$DATE_NUM, origin="1970-01-01")
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)
levels_country = rev(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","country"][order(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","prob"])])
as.character(levels_country)
levels_country = c("India", "Nepal", "China", 
                   as.character(levels_country)[!as.character(levels_country) %in% c("India","Nepal","China")])
fit_preds$country = factor(fit_preds$country, levels=levels_country)
fit_preds$continent = factor(fit_preds$continent)
# write.csv(fit_preds, file=".//data//GISAID//GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage.csv", row.names=F)


# PLOT MULTINOMIAL FIT ON LOGIT SCALE ####
ymin = 0.0001
ymax = 0.999
if (!is.null(fit_preds$asymp.LCL)) fit_preds$asymp.LCL[fit_preds$asymp.LCL<ymin] = ymin
if (!is.null(fit_preds$asymp.UCL)) fit_preds$asymp.UCL[fit_preds$asymp.UCL<ymin] = ymin
if (!is.null(fit_preds$asymp.UCL)) fit_preds$asymp.UCL[fit_preds$asymp.UCL>ymax] = ymax
fit_preds$prob[fit_preds$prob<ymin] = ymin
fit_preds$prob[fit_preds$prob>ymax] = ymax

data_agbyweekcountry1$variant = factor(data_agbyweekcountry1$variant, levels=levels_VARIANTS_plot)
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)

plot_preds_logit = qplot(data=fit_preds[fit_preds$variant!="Other"&fit_preds$country %in% sel_countries_min10,], 
                         x=date, y=prob, geom="blank") +
  facet_wrap(~ country) +
  # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
  #                 fill=variant
  # ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SARS-CoV2 LINEAGE FREQUENCIES",
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=2)*continent+country,\nall countries with >=10 BA.2.75 sequences shown (n=", nBA_2_75," BA.2.75 sequences in total)")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  scale_colour_manual("variant", values=tail(as.vector(lineage_cols_plot),-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$country %in% sel_countries_min10,],
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
  coord_cartesian(xlim=c(as.Date("2020-05-01"),NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_preds_logit

ggsave(file=file.path("plots", plotdir, "predictions global multinom fit_logit scale.png"), width=16, height=8)

# plot predicted values as Muller plot
muller_fit = ggplot(data=fit_preds[fit_preds$country %in% sel_countries_min10,], 
                    aes(x=date, y=prob, group=variant)) +
  facet_wrap(~ country, ncol=4) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SARS-CoV2 LINEAGE FREQUENCIES",
       subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=2)*continent+country,\nall countries with >=10 BA.2.75 sequences shown (n=", nBA_2_75," BA.2.75 sequences in total)")) +
  annotate("rect", xmin=max(GISAID_sel$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01")+5, 
           ymin=0, ymax=1, alpha=0.5, fill="white") # extrapolated part
muller_fit

ggsave(file=file.path("plots", plotdir,"predictions global multinom fit_muller plot.png"), width=14, height=7)



# map variant share onto case numbers

# sel_countries_top3 = c("India", "Nepal", "Singapore")
country_data = get_national_data(countries=sel_countries_min10, # sel_countries_top3,
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY IN\nTOP 3 COUNTRIES WITH MOST BA.2.75 CASES",
          subtitle="7 day simple moving average, using\nWHO case data accessed via covidregionaldata package")
# ggsave(file=file.path("plots", plotdir, "new cases countries with more than 10 sequenced BA_2_75 cases.png"), width=7, height=5)

fit_preds[(fit_preds$date==today)&(fit_preds$variant=="Omicron (BA.2.75)"),]

# fit_multi_predsbycountry = data.frame(emmeans(fit8_multi,
#                                                   ~ variant,
#                                                   by=c("DATE_NUM",
#                                                        "continent",
#                                                        "country"),
#                                                   at=list(DATE_NUM=today_num), # by=7 just to speed up things a bit
#                                                   mode="prob", df=NA,
#                                                   rg.limit=100000))
# 
# fit_multi_predsbycountry[fit_multi_predsbycountry$variant=="Omicron (BA.2.75)",]


fit_preds$totnewcases = 
  country_data$cases_new[match(interaction(fit_preds$country,
                                           fit_preds$date),
                               interaction(country_data$country, 
                                           country_data$date))]
library(zoo)
fit_preds = fit_preds %>% 
  group_by(country) %>% 
  mutate(totnewcases_smoothed = rollmean(totnewcases, 7, na.pad = T))
fit_preds$cases = fit_preds$totnewcases_smoothed*fit_preds$prob
fit_preds$cases[fit_preds$cases==0] = NA
fit_preds$cases[fit_preds$prob<0.001] = NA
fit_preds = fit_preds[fit_preds$country %in% sel_countries_min10,]
fit_preds$country = factor(fit_preds$country)

fit_preds2 = fit_preds
fit_preds2$cases[fit_preds2$cases==0] = NA
fit_preds2$cases[fit_preds2$cases<=1] = NA
fit_preds2$country = factor(fit_preds2$country)
fit_preds2$variant = factor(fit_preds2$variant, levels=levels_VARIANTS_plot)

ggplot(data=fit_preds2, 
       aes(x=date, y=cases)) + 
  facet_wrap(~ country, scale="free", ncol=4) +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_preds2, aes(x=date, y=totnewcases_smoothed, lwd=I(1.5)), 
            colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN\nCOUNTRIES WITH >=10 SEQUENCED BA.2.75 CASES") +
  scale_colour_manual("lineage", values=lineage_cols_plot) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=file.path("plots", plotdir,"confirmed cases multinomial fit by country.png"), width=18, height=8)


fit_preds3 = fit_preds2
# fit_preds3$cases[fit_preds3$date<=as.Date("2022-03-01")] = 0

ggplot(data=fit_preds2, 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y", ncol=4) +
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN\nCOUNTRIES WITH >=10 SEQUENCED BA.2.75 CASES") +
  coord_cartesian(xlim=c(as.Date("2021-06-01"),NA))

ggsave(file=file.path("plots", plotdir,"\\confirmed cases stacked area multinomial fit by country.png"), width=18, height=8)




# 4. ANALYSIS OF VOC LINEAGE FREQUENCIES IN INDIA ####

# load case data India per state
cases_india_bystate = read.csv("https://data.covid19bharat.org/csv/latest/states.csv") # cumulative cases
cases_india_bystate$Date = as.Date(cases_india_bystate$Date)
cases_india_bystate = cases_india_bystate[cases_india_bystate$Date >= as.Date("2020-06-01"),]
tail(cases_india_bystate)

cases_india_bystate = do.call(rbind,lapply(unique(cases_india_bystate$State), function (state) { df =  cases_india_bystate[cases_india_bystate$State==state,]
df$newcases = c(NA, diff(df$Confirmed))
df$newtests = c(NA, diff(df$Tested))
return(df)
} ))
cases_india_bystate$posratio = cases_india_bystate$newcases*100/cases_india_bystate$newtests
cases_india_bystate$posratio[cases_india_bystate$posratio>100] = NA
cases_india_bystate = cases_india_bystate[cases_india_bystate$State!="State Unassigned",]
cases_india_bystate = cases_india_bystate[cases_india_bystate$State!="India",]

levels_STATES = sort(unique(cases_india_bystate$State))
levels_STATES

# plot new cases per day by state
ggplot(data=cases_india_bystate[cases_india_bystate$Date>=as.Date("2022-01-01"),],
       aes(x=Date, y=newcases, 
           group=State)) +
  facet_wrap(~ State, scale="free", ncol=5) +
  geom_smooth(aes(lwd=I(1), colour=State), method="loess", span=0.1, se=FALSE) +
  # geom_line(aes(lwd=I(1), colour=State)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY STATE IN INDIA") +
  scale_y_log10() +
  theme(legend.position = "none") # +
#  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=file.path("plots", plotdir,"\\india_cases per day by state.png"), width=12, height=12)

# plot positivity ratios per day by state
ggplot(data=cases_india_bystate[cases_india_bystate$Date>=as.Date("2022-01-01"),],
       aes(x=Date, y=posratio, 
           group=State)) +
  facet_wrap(~ State, scale="free", ncol=5) +
  geom_smooth(aes(lwd=I(1), colour=State), method="loess", span=0.2, se=FALSE) +
  # geom_line(aes(lwd=I(1), colour=State)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("Positivity ratios") +
  ggtitle("POSITIVITY RATIOS BY STATE IN INDIA") +
  scale_y_log10() +
  theme(legend.position = "none") # +
#  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=file.path("plots", plotdir,"india_pos ratios by state.png"), width=12, height=12)


# GISAID data India

GISAID_india = GISAID_sel[as.character(GISAID_sel$country)=="India",]
GISAID_india$location = droplevels(GISAID_india$location)
GISAID_india = GISAID_india[!is.na(GISAID_india$location),]
GISAID_india = GISAID_india[!toupper(as.character(GISAID_india$location))=="NAIROBI",]
nrow(GISAID_india) # 202556
table(GISAID_india$variant, GISAID_india$location)
sort(unique(toupper(as.character(GISAID_india$location))))

GISAID_india$location = factor(toupper(as.character(GISAID_india$location)),
                               levels = c("ANDAMAN AND NICOBAR ISLANDS", "ANDHRA PRADESH", "ARUNACHAL PRADESH", "ASSAM", "BANGLADESH", "BIHAR", "CHANDIGARH", "CHHATISGARH",  "CHHATTISGARH", "DADRA AND NAGAR HAVELI",                   "DADRA AND NAGAR HAVELI AND DAMAN AND DIU", "DELHI", "GOA", "GUJARAT", "GUJART",  "GUJRAT",  "HARAYANA", "HARYANA", "HIMACHAL PRADESH", "JAMMU",             "JAMMU & KASHMIR",   "JAMMU AND KASHMIR", "JAMMU AND KASHMĪR", "JARKHAND",  "JHARKHAND", "KARNATAKA", "KERALA", "LADAKH", "LAKSHADWEEP", "MADHYA PRADESH", "MAHARASHTRA", "MAHARASTHRA", "MAHARSHTRA",  "MANIPUR", "MEGHALAYA", "MIZORAM", "MUMBAI",      "NAGALAND", "NEW DELHI", "ODISHA", "PONDICHERRY", "PUDUCHERRY", "PUNJAB", "RAJASTHAN", "SIKKIM", "TAMIL NADU", "TAMILNADU",  "TELANGANA", "TRIPURA", "UTTAR PRADESH", "UTTARAKHAND", "WEST BENGAL" ),
                               labels = c("Andaman and Nicobar Islands", "Andhra Pradesh", "Arunachal Pradesh", "Assam", "",           "Bihar", "Chandigarh", "Chhattisgarh", "Chhattisgarh", "Dadra and Nagar Haveli and Daman and Diu", "Dadra and Nagar Haveli and Daman and Diu", "Delhi", "Goa", "Gujarat", "Gujarat", "Gujarat", "Haryana",  "Haryana", "Himachal Pradesh", "Jammu and Kashmir", "Jammu and Kashmir", "Jammu and Kashmir", "Jammu and Kashmir", "Jharkhand", "Jharkhand", "Karnataka", "Kerala", "Ladakh", "Lakshadweep", "Madhya Pradesh", "Maharashtra", "Maharashtra", "Maharashtra", "Manipur", "Meghalaya", "Mizoram", "Maharashtra", "Nagaland", "Delhi",     "Odisha", "Puducherry",  "Puducherry", "Punjab", "Rajasthan", "Sikkim", "Tamil Nadu", "Tamil Nadu", "Telangana", "Tripura", "Uttar Pradesh", "Uttarakhand", "West Bengal" ))
GISAID_india = GISAID_india[as.character(GISAID_india$location)!="",]
GISAID_india$location = droplevels(GISAID_india$location)

# Muller plot of aggregated data using all data from India using all states
nBA_2_75_india = sum(GISAID_india$variant=="Omicron (BA.2.75)", na.rm=T)
nBA_2_75_india # 2492

# PS this was for Ewen Callaway's Nature piece on BA.2.75
# aggregated by week for selected variant lineages for whole of India using data from all states
data_agbyweek1 = as.data.frame(table(GISAID_india$floor_date,
                                     GISAID_india$variant))
colnames(data_agbyweek1) = c("floor_date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date,
                               data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$variant=="Omicron (BA.1)","total"]) == nrow(GISAID_india) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$variant = factor(data_agbyweek1$variant,
                                levels=levels_VARIANTS)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

# MULLER PLOT (RAW DATA)
levels_VARIANTS_plot_india = levels_VARIANTS_plot[levels_VARIANTS_plot %in% unique(data_agbyweek1$variant)]
data_agbyweek1$variant2 = factor(data_agbyweek1$variant, 
                                 levels=levels_VARIANTS_plot_india)
lineage_cols_plot_india = lineage_cols_plot[levels_VARIANTS_plot %in% unique(data_agbyweek1$variant)]
muller_india_raw1_all = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=variant2)) +
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant2, group=variant2), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
       subtitle=paste0("Raw GISAID data up to ",today,"\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total with valid date)"))
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_india_raw1_all

ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data all states.png"), width=7, height=5)
# ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data all states.pdf"), width=7, height=5)
# ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data all states.svg"), width=7, height=5)
# saveRDS(muller_india_raw1_all, file.path("plots", plotdir,"india_muller plots_raw data all states.rds"))
# data_agbyweek2=data_agbyweek1
# data_agbyweek2$variant2=NULL
# write.csv(data_agbyweek2, file=file.path("plots", plotdir,"india_muller plots_raw data all states.csv"), row.names=F)


# per state analysis, using data from states with at least 5 BA.2.75 sequences
tab = as.data.frame(table(GISAID_india$variant, GISAID_india$location))
tab_BA_2_75 = tab[tab$Var1=="Omicron (BA.2.75)",]
tab_BA_2_75 = tab_BA_2_75[order(tab_BA_2_75$Freq, decreasing=T),]
tab_BA_2_75[tab_BA_2_75$Freq!=0,]
sel_states = as.character(tab_BA_2_75$Var2[tab_BA_2_75$Freq>=5])
sel_states

GISAID_india = GISAID_india[as.character(GISAID_india$location) %in%
                              sel_states, ]
GISAID_india$location = droplevels(GISAID_india$location)

# GISAID_india = GISAID_india[GISAID_india$location=="Maharashtra",]
table(GISAID_india$variant)
table(GISAID_india[GISAID_india$date>=as.Date("2022-06-01"),"variant"],
      GISAID_india[GISAID_india$date>=as.Date("2022-06-01"),"pango_lineage"])

nBA_2_75_india = sum(GISAID_india$variant=="Omicron (BA.2.75)", na.rm=T)
nBA_2_75_india # 2479

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of India
data_agbyweek1 = as.data.frame(table(GISAID_india$floor_date, 
                                     GISAID_india$variant))
colnames(data_agbyweek1) = c("floor_date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date, 
                               data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$variant=="Omicron (BA.1)","total"]) == nrow(GISAID_india) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$variant = factor(data_agbyweek1$variant, 
                                levels=levels_VARIANTS)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

# aggregated by week and state
data_agbyweekregion1 = as.data.frame(table(GISAID_india$floor_date, 
                                           GISAID_india$location, 
                                           GISAID_india$variant))
colnames(data_agbyweekregion1) = c("floor_date", "division", "variant", "count")
data_agbyweekregion1_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$floor_date,data_agbyweekregion1$division),
                                                                  interaction(data_agbyweekregion1_sum$floor_date,data_agbyweekregion1_sum$division))]
data_agbyweekregion1$collection_date = as.Date(as.character(data_agbyweekregion1$floor_date))
data_agbyweekregion1$variant = factor(data_agbyweekregion1$variant, levels=levels_VARIANTS)
data_agbyweekregion1$division = factor(data_agbyweekregion1$division, levels=levels_STATES)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]
data_agbyweekregion1$floor_date = NULL
data_agbyweekregion1$location = data_agbyweekregion1$division
data_agbyweekregion1$location = droplevels(data_agbyweekregion1$location)
data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
write.csv(data_agbyweekregion1, file=".//data//GISAID//GISAID aggregated counts by start of week and lineage_india_bystate_allstateswithatleast_5_BA_2_75.csv", row.names=F)


# MULLER PLOT (RAW DATA)
data_agbyweek1$variant2 = factor(data_agbyweek1$variant, 
                                 levels=levels_VARIANTS_plot_india)
muller_india_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=variant2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant2, group=variant2), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VOCs IN INDIA",
       subtitle=paste0("Raw GISAID data up to ",today, ", subsetted to states with at least 5 BA.2.75 sequences\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_india_raw1

ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data.png"), width=7, height=5)

data_agbyweekregion1$variant2 = factor(data_agbyweekregion1$variant, levels=levels_VARIANTS_plot_india)
muller_indiabystate_raw2 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=variant2)) +
  facet_wrap(~ division, ncol=4) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant2, group=variant2), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
       subtitle=paste0("Raw GISAID data up to ",today, ", subsetted to states with at least 5 BA.2.75 sequences\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) 

# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_indiabystate_raw2

ggsave(file=file.path("plots", plotdir,"india_muller plots by state_raw data.png"), width=14, height=7)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
GISAID_india$variant = factor(GISAID_india$variant, levels=levels_VARIANTS)
GISAID_india$variant = droplevels(GISAID_india$variant)
# fit1_india_multi = nnet::multinom(variant ~ scale(DATE_NUM)+location, data=GISAID_india, maxit=1000)
# fit2_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+location, data=GISAID_india, maxit=1000)
# fit3_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)+location, data=GISAID_india, maxit=1000)
# fit1_india_multi = nnet::multinom(variant ~ scale(DATE_NUM)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
fit2_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
# fit3_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
# fit4_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*location, weights=count, data=data_agbyweekregion1, maxit=10000)
# fit5_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)*location, weights=count, data=data_agbyweekregion1, maxit=10000)
# BIC(fit1_india_multi, fit2_india_multi, fit3_india_multi, fit4_india_multi, fit5_india_multi) 
# fit3_india_multi has best BIC, but fit2_india_multi almost as good & since data is sparse I will use
# that simpler model

bestfit_india_multi = fit2_india_multi
  
# growth rate advantage compared to BA.5
# max(GISAID_india$date) # 2022-07-03
# emtrindia = emtrends(fit3_india_multi, trt.vs.ctrl ~ variant, 
#                      by="DATE_NUM",
#                      var="DATE_NUM",  
#                      mode="latent",
#                      at=list(location="Delhi",
#                              DATE_NUM=
#                                as.numeric(max(GISAID_india$date,na.rm=T))))
# delta_r_india = data.frame(confint(emtrindia, 
#                                    adjust="none", df=NA)$contrasts, 
#                            p.value=as.data.frame(emtrindia$contrasts)$p.value)
# delta_r_india
# write.csv(delta_r_india, ".//plots//GISAID_BA_2_75//growth advantage BA.2.75 vs BA.5 India.csv", row.names=F)


# pairwise growth rate difference (differences in growth rate per day) 
emtrindia_pairw = emtrends(bestfit_india_multi, revpairwise ~ variant,  
                           var="DATE_NUM",  by="DATE_NUM",
                           mode="latent",
                           at=list(location="Delhi",
                                   DATE_NUM=max(GISAID_india$DATE_NUM, na.rm=T)))
delta_r_india_pairw = data.frame(confint(emtrindia_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindia_pairw$contrasts)$p.value)
delta_r_india_pairw
write.csv(delta_r_india_pairw, ".//plots//GISAID_BA_2_75//growth advantage BA.2.75 vs all other strains India.csv", row.names=F)


# PLOT MULTINOMIAL FIT
extrapolate = 60
date.from = as.numeric(as.Date("2020-06-01"))
date.to = today_num+extrapolate


# PLOT MODEL FIT ####

# multinomial model predictions by province (fastest, but no confidence intervals, use emmeans to get conf intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to),
                            location=levels(GISAID_india$location)))
fit_preds = data.frame(predgrid, as.data.frame(predict(bestfit_india_multi, 
                                                       newdata=predgrid, type="prob")),check.names=F)
fit_preds = gather(fit_preds, variant, prob, all_of(levels_VARIANTS_plot_india), factor_key=TRUE)
fit_preds$date = as.Date(fit_preds$DATE_NUM, origin="1970-01-01")
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)
levels_state = rev(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","location"][order(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","prob"])])
as.character(levels_state)
fit_preds$location = factor(fit_preds$location, levels=levels_state)

# write.csv(fit_preds, file=".//data//GISAID//GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage_india_bystate_allstateswithatleast_5_BA_2_75.csv", row.names=F)

# to get predictions with confidence intervals (but slow) :
# fit_india_multi_predsbystate = data.frame(emmeans(bestfit_india_multi,
#                                                   ~ variant,
#                                                   by=c("DATE_NUM", 
#                                                        "location"),
#                                                   at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
#                                                   mode="prob", df=NA,
#                                                   rg.limit=100000))
# fit_india_multi_predsbystate$date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
# fit_india_multi_predsbystate$variant = factor(fit_india_multi_predsbystate$variant, levels=levels_VARIANTS_plot)
# fit_india_multi_predsbystate$location = factor(fit_india_multi_predsbystate$location, levels=levels_STATES)


# plot model prediction as line plot on logit scale:

fit_india_multi_preds2 = fit_preds
ymin = 0.0001
ymax = 0.999
if (!is.null(fit_india_multi_preds2$asymp.LCL)) fit_india_multi_preds2$asymp.LCL[fit_india_multi_preds2$asymp.LCL<ymin] = ymin
if (!is.null(fit_india_multi_preds2$asymp.UCL)) fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL<ymin] = ymin
if (!is.null(fit_india_multi_preds2$asymp.UCL)) fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL>ymax] = ymax
fit_india_multi_preds2$prob[fit_india_multi_preds2$prob<ymin] = ymin

plot_india_mfit_logit = qplot(data=fit_india_multi_preds2[fit_india_multi_preds2$variant!="Other",], 
                              x=date, y=prob, geom="blank") +
  facet_wrap(~location, ncol=4) +
  # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
  #                 fill=variant
  # ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
          subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=2)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
  xaxis +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(lineage_cols_plot_india,-1)) +
  scale_colour_manual("variant", values=tail(lineage_cols_plot_india,-1)) +
  geom_point(data=data_agbyweekregion1[data_agbyweekregion1$variant!="Other",],
             aes(x=collection_date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.001, 2), limits=c(1,max(data_agbyweekregion1$total)), breaks=c(10,100, 1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today_num+extrapolate+1), 
                  ylim=c(0.001, 0.99), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_india_mfit_logit

ggsave(file=file.path("plots", plotdir,"india_multinom fit_logit scale.png"), width=14, height=8)


# plot model prediction as line plot on response scale:
plot_india_mfit = qplot(data=fit_india_multi_preds2[fit_india_multi_preds2$variant!="Other",], x=date, 
                        y=100*prob, geom="blank") +
  facet_wrap(~location, ncol=4) +
  # geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, 
  #                 colour=NULL,
  #                 fill=variant
  # ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among new infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
          subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=2)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
  xaxis +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today_num+extrapolate+1), 
                  ylim=c(0, 100)) +
  scale_fill_manual("variant", values=tail(lineage_cols_plot_india,-1)) +
  scale_colour_manual("variant", values=tail(lineage_cols_plot_india,-1)) +
  geom_point(data=data_agbyweekregion1[data_agbyweekregion1$variant!="Other",],
             aes(x=collection_date, y=100*prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.001, 2), limits=c(1,max(data_agbyweekregion1$total)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_india_mfit

ggsave(file=file.path("plots", plotdir,"india_multinom fit_response scale.png"), width=14, height=8)


# plot predicted values as Muller plot
muller_indiabystate_fit = ggplot(data=fit_preds, 
                                 aes(x=date, y=prob, group=variant)) +
  facet_wrap(~ location, ncol=4) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
       subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=2)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, # extrapolated part
                   xmax=as.Date(date.to, origin="1970-01-01")+5, 
           ymin=0, ymax=1, alpha=0.5, fill="white") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) 
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_indiabystate_fit

ggsave(file=file.path("plots", plotdir,"india_muller plots by state_multinom fit.png"), width=14, height=7)



# project multinomial fit onto case data ####

# changed to new source
# https://github.com/DataKind-BLR/covid19bharat_data/
# https://data.covid19bharat.org/csv/latest/states.csv

cases_india_bystate2 = cases_india_bystate[cases_india_bystate$State %in% 
                                             sel_states,]
colnames(cases_india_bystate2)[2]="STATE"

newdat = expand.grid(DATE_NUM=seq(as.numeric(min(cases_india_bystate2$Date)),
                                  as.numeric(max(cases_india_bystate2$Date))),
                     location=unique(as.character(cases_india_bystate2$STATE)))
fit_india_multi_predsbystate = data.frame(newdat,
                        predict(bestfit_india_multi, 
                                  newdata = newdat,
                                  type = "prob"), check.names=F) 
library(tidyr)
fit_india_multi_predsbystate = gather(fit_india_multi_predsbystate, variant, prob, all_of(levels_VARIANTS_plot_india))
fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
fit_india_multi_predsbystate$variant = factor(fit_india_multi_predsbystate$variant, levels=levels_VARIANTS)
# colnames(fit_india_multi_predsbystate)[2] = "location"
# fit_india_multi_predsbystate$location = factor(fit_india_multi_predsbystate$location)
fit_india_multi_predsbystate$totnewcases = 
  cases_india_bystate2$newcases[match(interaction(fit_india_multi_predsbystate$location,
                                                  fit_india_multi_predsbystate$collection_date),
                                      interaction(cases_india_bystate2$STATE, 
                                                  cases_india_bystate2$Date))]                                                    
fit_india_multi_predsbystate$totnewcases[which(fit_india_multi_predsbystate$totnewcases>5&fit_india_multi_predsbystate$collection_date>as.Date("2021-08-01")&fit_india_multi_predsbystate$collection_date<as.Date("2021-09-01")&fit_india_multi_predsbystate$location=="Chandigarh")] = 5 # outliers

library(zoo)
fit_india_multi_predsbystate = fit_india_multi_predsbystate %>% 
  group_by(location) %>% 
  mutate(totnewcases_smoothed = rollmean(totnewcases, 7, na.pad = T))
fit_india_multi_predsbystate$totnewtests = 
  cases_india_bystate2$newtests[match(interaction(fit_india_multi_predsbystate$location,
                                                  fit_india_multi_predsbystate$collection_date),
                                      interaction(cases_india_bystate2$STATE, 
                                                  cases_india_bystate2$Date))]                                                                                                                                                                                                                                       
fit_india_multi_predsbystate$totnewtests[fit_india_multi_predsbystate$location=="Haryana"&fit_india_multi_predsbystate$collection_date=="2022-06-23"] = NA # likely error in data
fit_india_multi_predsbystate$totnewtests[fit_india_multi_predsbystate$location=="Odisha"&fit_india_multi_predsbystate$collection_date=="2022-07-04"] = NA # likely error in data
fit_india_multi_predsbystate = fit_india_multi_predsbystate %>% 
  group_by(location) %>% 
  mutate(totnewtests_smoothed = rollmean(totnewtests, 7, na.pad = T))
fit_india_multi_predsbystate$cases = fit_india_multi_predsbystate$totnewcases_smoothed*fit_india_multi_predsbystate$prob
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$cases==0] = NA
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$location=="Chandigarh"&fit_india_multi_predsbystate$cases>2000] = NA
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$prob<0.001] = NA
fit_india_multi_predsbystate$posratio = 100*fit_india_multi_predsbystate$cases/fit_india_multi_predsbystate$totnewtests_smoothed
fit_india_multi_predsbystate$posratio[fit_india_multi_predsbystate$posratio<0] = NA
fit_india_multi_predsbystate$posratio[fit_india_multi_predsbystate$posratio>100] = NA
fit_india_multi_predsbystate$posratio[fit_india_multi_predsbystate$variant=="Delta"&fit_india_multi_predsbystate$posratio>10&fit_india_multi_predsbystate$collection_date>as.Date("2021-08-01")&fit_india_multi_predsbystate$collection_date<as.Date("2021-09-01")&fit_india_multi_predsbystate$location=="Chandigarh"] = 0 # outliers

# fit_india_multi_predsbystate$collection_date[which((fit_india_multi_predsbystate$posratio>10)&(fit_india_multi_predsbystate$variant=="Omicron (BA.2.76)")&(fit_india_multi_predsbystate$location=="Odisha"))]
# "2022-07-04" # outlier
fit_india_multi_predsbystate$location = factor(fit_india_multi_predsbystate$location, levels=levels_state)

fit_india_multi_predsbystate2 = fit_india_multi_predsbystate
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases==0] = NA
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases<=1] = NA
fit_india_multi_predsbystate2$location = factor(fit_india_multi_predsbystate2$location, levels=levels_state)
cases_india_bystate2$location = factor(cases_india_bystate2$STATE, levels=levels_state)
fit_india_multi_predsbystate2$variant = factor(fit_india_multi_predsbystate2$variant, levels=levels_VARIANTS_plot)

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_india_multi_predsbystate2, aes(x=collection_date, y=totnewcases_smoothed, lwd=I(1.5)), colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  scale_colour_manual("lineage", values=lineage_cols_plot_india) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA)) +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_confirmed cases multinomial fit.png"), width=16, height=8)

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))


ggsave(file=file.path("plots", plotdir,"india_confirmed cases stacked area multinomial fit.png"), width=16, height=8)

fit_india_multi_predsbystate3 = fit_india_multi_predsbystate2
fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date<as.Date("2022-05-01"),"cases"] = 0
fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date<as.Date("2022-05-01"),"posratio"] = 0

ggplot(data=fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date>=as.Date("2022-04-30"),], 
       aes(x=collection_date, y=cases, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_confirmed cases stacked area multinomial fit_ZOOMED.png"), width=12, height=6)

# graph with positivity ratios by variant
ggplot(data=fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date>=as.Date("2022-04-30"),], 
       aes(x=collection_date, y=posratio, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("% of tests that were positive") +
  ggtitle("SARS-CoV2 POSITIVITY RATIOS BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
          " plus records from last days\ncase & testing data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_positivity ratios stacked area multinomial fit_ZOOMED.png"), width=12, height=6)
# write.csv(fit_india_multi_predsbystate, file="./data/GISAID/GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage_india_bystate__plus_positivity_ratios_by_state_allstateswithatleast_5_BA_2_75.csv", row.names=F)


ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=posratio, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot_india) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("% of tests that were positive") +
  ggtitle("SARS-CoV2 POSITIVITY RATIOS BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase & testing data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_positivity ratios stacked area multinomial fit.png"), width=16, height=8)


# 5. ANALYSIS OF LINEAGE FREQUENCIES IN BELGIUM (CODE NOT FINISHED/WRITTEN YET) ####

# load case data India per state
cases_india_bystate = read.csv("https://data.covid19bharat.org/csv/latest/states.csv") # cumulative cases
cases_india_bystate$Date = as.Date(cases_india_bystate$Date)
cases_india_bystate = cases_india_bystate[cases_india_bystate$Date >= as.Date("2020-06-01"),]
tail(cases_india_bystate)

cases_india_bystate = do.call(rbind,lapply(unique(cases_india_bystate$State), function (state) { df =  cases_india_bystate[cases_india_bystate$State==state,]
df$newcases = c(NA, diff(df$Confirmed))
df$newtests = c(NA, diff(df$Tested))
return(df)
} ))
cases_india_bystate$posratio = cases_india_bystate$newcases*100/cases_india_bystate$newtests
cases_india_bystate$posratio[cases_india_bystate$posratio>100] = NA
cases_india_bystate = cases_india_bystate[cases_india_bystate$State!="State Unassigned",]
cases_india_bystate = cases_india_bystate[cases_india_bystate$State!="India",]

levels_STATES = sort(unique(cases_india_bystate$State))
levels_STATES

# plot new cases per day by state
ggplot(data=cases_india_bystate[cases_india_bystate$Date>=as.Date("2022-01-01"),],
       aes(x=Date, y=newcases, 
           group=State)) +
  facet_wrap(~ State, scale="free", ncol=5) +
  geom_smooth(aes(lwd=I(1), colour=State), method="loess", span=0.1, se=FALSE) +
  # geom_line(aes(lwd=I(1), colour=State)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY STATE IN INDIA") +
  scale_y_log10() +
  theme(legend.position = "none") # +
#  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=file.path("plots", plotdir,"\\india_cases per day by state.png"), width=12, height=12)

# plot positivity ratios per day by state
ggplot(data=cases_india_bystate[cases_india_bystate$Date>=as.Date("2022-01-01"),],
       aes(x=Date, y=posratio, 
           group=State)) +
  facet_wrap(~ State, scale="free", ncol=5) +
  geom_smooth(aes(lwd=I(1), colour=State), method="loess", span=0.2, se=FALSE) +
  # geom_line(aes(lwd=I(1), colour=State)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right",
                     axis.title.x=element_blank()) +
  ylab("Positivity ratios") +
  ggtitle("POSITIVITY RATIOS BY STATE IN INDIA") +
  scale_y_log10() +
  theme(legend.position = "none") # +
#  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=file.path("plots", plotdir,"india_pos ratios by state.png"), width=12, height=12)


# GISAID data India

GISAID_india = GISAID_sel[as.character(GISAID_sel$country)=="India",]
GISAID_india$location = droplevels(GISAID_india$location)
GISAID_india = GISAID_india[!is.na(GISAID_india$location),]
GISAID_india = GISAID_india[!toupper(as.character(GISAID_india$location))=="NAIROBI",]
nrow(GISAID_india) # 191323
table(GISAID_india$variant, GISAID_india$location)
sort(unique(toupper(as.character(GISAID_india$location))))

GISAID_india$location = factor(toupper(as.character(GISAID_india$location)),
                               levels = c("ANDAMAN AND NICOBAR ISLANDS", "ANDHRA PRADESH", "ARUNACHAL PRADESH", "ASSAM", "BANGLADESH", "BIHAR", "CHANDIGARH", "CHHATISGARH",  "CHHATTISGARH", "DADRA AND NAGAR HAVELI",                   "DADRA AND NAGAR HAVELI AND DAMAN AND DIU", "DELHI", "GOA", "GUJARAT", "GUJART",  "GUJRAT",  "HARAYANA", "HARYANA", "HIMACHAL PRADESH", "JAMMU",             "JAMMU & KASHMIR",   "JAMMU AND KASHMIR", "JAMMU AND KASHMĪR", "JARKHAND",  "JHARKHAND", "KARNATAKA", "KERALA", "LADAKH", "LAKSHADWEEP", "MADHYA PRADESH", "MAHARASHTRA", "MAHARASTHRA", "MAHARSHTRA",  "MANIPUR", "MEGHALAYA", "MIZORAM", "MUMBAI",      "NAGALAND", "NEW DELHI", "ODISHA", "PONDICHERRY", "PUDUCHERRY", "PUNJAB", "RAJASTHAN", "SIKKIM", "TAMIL NADU", "TAMILNADU",  "TELANGANA", "TRIPURA", "UTTAR PRADESH", "UTTARAKHAND", "WEST BENGAL" ),
                               labels = c("Andaman and Nicobar Islands", "Andhra Pradesh", "Arunachal Pradesh", "Assam", "",           "Bihar", "Chandigarh", "Chhattisgarh", "Chhattisgarh", "Dadra and Nagar Haveli and Daman and Diu", "Dadra and Nagar Haveli and Daman and Diu", "Delhi", "Goa", "Gujarat", "Gujarat", "Gujarat", "Haryana",  "Haryana", "Himachal Pradesh", "Jammu and Kashmir", "Jammu and Kashmir", "Jammu and Kashmir", "Jammu and Kashmir", "Jharkhand", "Jharkhand", "Karnataka", "Kerala", "Ladakh", "Lakshadweep", "Madhya Pradesh", "Maharashtra", "Maharashtra", "Maharashtra", "Manipur", "Meghalaya", "Mizoram", "Maharashtra", "Nagaland", "Delhi",     "Odisha", "Puducherry",  "Puducherry", "Punjab", "Rajasthan", "Sikkim", "Tamil Nadu", "Tamil Nadu", "Telangana", "Tripura", "Uttar Pradesh", "Uttarakhand", "West Bengal" ))
GISAID_india = GISAID_india[as.character(GISAID_india$location)!="",]
GISAID_india$location = droplevels(GISAID_india$location)

# Muller plot of aggregated data using all data from India using all states
nBA_2_75_india = sum(GISAID_india$variant=="Omicron (BA.2.75)", na.rm=T)
nBA_2_75_india # 1521

# PS this was for Ewen Callaway's Nature piece on BA.2.75
# aggregated by week for selected variant lineages for whole of India using data from all states
data_agbyweek1 = as.data.frame(table(GISAID_india$floor_date,
                                     GISAID_india$variant))
colnames(data_agbyweek1) = c("floor_date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date,
                               data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$variant=="Omicron (BA.1)","total"]) == nrow(GISAID_india) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$variant = factor(data_agbyweek1$variant,
                                levels=levels_VARIANTS)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

# MULLER PLOT (RAW DATA)
data_agbyweek1$variant2 = factor(data_agbyweek1$variant, levels=levels_VARIANTS_plot)
muller_india_raw1_all = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=variant2)) +
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant2, group=variant2), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
       subtitle=paste0("Raw GISAID data up to ",today,"\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)"))
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_india_raw1_all

ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data all states.png"), width=7, height=5)
ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data all states.pdf"), width=7, height=5)
ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data all states.svg"), width=7, height=5)
saveRDS(muller_india_raw1_all, file.path("plots", plotdir,"india_muller plots_raw data all states.rds"))
data_agbyweek2=data_agbyweek1
data_agbyweek2$variant2=NULL
write.csv(data_agbyweek2, file=file.path("plots", plotdir,"india_muller plots_raw data all states.csv"), row.names=F)


# per state analysis, using data from states with at least 5 BA.2.75 sequences
tab = as.data.frame(table(GISAID_india$variant, GISAID_india$location))
tab_BA_2_75 = tab[tab$Var1=="Omicron (BA.2.75)",]
tab_BA_2_75 = tab_BA_2_75[order(tab_BA_2_75$Freq, decreasing=T),]
tab_BA_2_75[tab_BA_2_75$Freq!=0,]
sel_states = as.character(tab_BA_2_75$Var2[tab_BA_2_75$Freq>=5])
sel_states
# [1] "West Bengal"      "Maharashtra"      "Odisha"           "Delhi"           
# [5] "Telangana"        "Rajasthan"        "Gujarat"          "Himachal Pradesh"
# [9] "Haryana"          "Tamil Nadu"       "Sikkim"           "Assam"           
# [13] "Manipur"          "Chhattisgarh"     "Madhya Pradesh"   "Chandigarh"   

GISAID_india = GISAID_india[as.character(GISAID_india$location) %in%
                              sel_states, ]
GISAID_india$location = droplevels(GISAID_india$location)

# GISAID_india = GISAID_india[GISAID_india$location=="Maharashtra",]
table(GISAID_india$variant)
table(GISAID_india[GISAID_india$date>=as.Date("2022-06-01"),"variant"],
      GISAID_india[GISAID_india$date>=as.Date("2022-06-01"),"pango_lineage"])

nBA_2_75_india = sum(GISAID_india$variant=="Omicron (BA.2.75)", na.rm=T)
nBA_2_75_india # 1505

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of India
data_agbyweek1 = as.data.frame(table(GISAID_india$floor_date, 
                                     GISAID_india$variant))
colnames(data_agbyweek1) = c("floor_date", "variant", "count")
data_agbyweek1_sum = aggregate(count ~ floor_date, 
                               data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
sum(data_agbyweek1[data_agbyweek1$variant=="Omicron (BA.1)","total"]) == nrow(GISAID_india) # correct
data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
data_agbyweek1$variant = factor(data_agbyweek1$variant, 
                                levels=levels_VARIANTS)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
data_agbyweek1$floor_date = NULL

# aggregated by week and state
data_agbyweekregion1 = as.data.frame(table(GISAID_india$floor_date, 
                                           GISAID_india$location, 
                                           GISAID_india$variant))
colnames(data_agbyweekregion1) = c("floor_date", "division", "variant", "count")
data_agbyweekregion1_sum = aggregate(count ~ floor_date + division, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$floor_date,data_agbyweekregion1$division),
                                                                  interaction(data_agbyweekregion1_sum$floor_date,data_agbyweekregion1_sum$division))]
data_agbyweekregion1$collection_date = as.Date(as.character(data_agbyweekregion1$floor_date))
data_agbyweekregion1$variant = factor(data_agbyweekregion1$variant, levels=levels_VARIANTS)
data_agbyweekregion1$division = factor(data_agbyweekregion1$division, levels=levels_STATES)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]
data_agbyweekregion1$floor_date = NULL
data_agbyweekregion1$location = data_agbyweekregion1$division
data_agbyweekregion1$location = droplevels(data_agbyweekregion1$location)
data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
write.csv(data_agbyweekregion1, file=".//data//GISAID//GISAID aggregated counts by start of week and lineage_india_bystate_allstateswithatleast_5_BA_2_75.csv", row.names=F)


# MULLER PLOT (RAW DATA)
data_agbyweek1$variant2 = factor(data_agbyweek1$variant, levels=levels_VARIANTS_plot)
muller_india_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=variant2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant2, group=variant2), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VOCs IN INDIA",
       subtitle=paste0("Raw GISAID data up to ",today, ", subsetted to states with at least 5 BA.2.75 sequences\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_india_raw1

ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data.png"), width=7, height=5)

data_agbyweekregion1$variant2 = factor(data_agbyweekregion1$variant, levels=levels_VARIANTS_plot)
muller_indiabystate_raw2 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=variant2)) +
  facet_wrap(~ division, ncol=4) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant2, group=variant2), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
       subtitle=paste0("Raw GISAID data up to ",today, ", subsetted to states with at least 5 BA.2.75 sequences\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) 

# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_indiabystate_raw2

ggsave(file=file.path("plots", plotdir,"india_muller plots by state_raw data.png"), width=14, height=7)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
GISAID_india$variant = factor(GISAID_india$variant, levels=levels_VARIANTS)
# fit1_india_multi = nnet::multinom(variant ~ scale(DATE_NUM)+location, data=GISAID_india, maxit=1000)
# fit2_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+location, data=GISAID_india, maxit=1000)
# fit3_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)+location, data=GISAID_india, maxit=1000)
# fit1_india_multi = nnet::multinom(variant ~ scale(DATE_NUM)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
fit2_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
# fit3_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
# fit4_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*location, weights=count, data=data_agbyweekregion1, maxit=10000)
# fit5_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)*location, weights=count, data=data_agbyweekregion1, maxit=10000)
# BIC(fit1_india_multi, fit2_india_multi, fit3_india_multi, fit4_india_multi, fit5_india_multi) 
# fit3_india_multi has best BIC, but fit2_india_multi almost as good & since data is sparse I will use
# that simpler model

bestfit_india_multi = fit2_india_multi

# growth rate advantage compared to BA.5
# max(GISAID_india$date) # 2022-07-03
# emtrindia = emtrends(fit3_india_multi, trt.vs.ctrl ~ variant, 
#                      by="DATE_NUM",
#                      var="DATE_NUM",  
#                      mode="latent",
#                      at=list(location="Delhi",
#                              DATE_NUM=
#                                as.numeric(max(GISAID_india$date,na.rm=T))))
# delta_r_india = data.frame(confint(emtrindia, 
#                                    adjust="none", df=NA)$contrasts, 
#                            p.value=as.data.frame(emtrindia$contrasts)$p.value)
# delta_r_india
# write.csv(delta_r_india, ".//plots//GISAID_BA_2_75//growth advantage BA.2.75 vs BA.5 India.csv", row.names=F)


# pairwise growth rate difference (differences in growth rate per day) 
emtrindia_pairw = emtrends(bestfit_india_multi, revpairwise ~ variant,  
                           var="DATE_NUM",  by="DATE_NUM",
                           mode="latent",
                           at=list(location="Delhi",
                                   DATE_NUM=max(GISAID_india$DATE_NUM, na.rm=T)))
delta_r_india_pairw = data.frame(confint(emtrindia_pairw, 
                                         adjust="none", df=NA)$contrasts, 
                                 p.value=as.data.frame(emtrindia_pairw$contrasts)$p.value)
delta_r_india_pairw
write.csv(delta_r_india_pairw, ".//plots//GISAID_BA_2_75//growth advantage BA.2.75 vs all other strains India.csv", row.names=F)


# PLOT MULTINOMIAL FIT
extrapolate = 60
date.from = as.numeric(as.Date("2020-06-01"))
date.to = today_num+extrapolate


# PLOT MODEL FIT ####

# multinomial model predictions by province (fastest, but no confidence intervals, use emmeans to get conf intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to),
                            location=levels(GISAID_india$location)))
fit_preds = data.frame(predgrid, as.data.frame(predict(bestfit_india_multi, 
                                                       newdata=predgrid, type="prob")),check.names=F)
fit_preds = gather(fit_preds, variant, prob, all_of(levels_VARIANTS), factor_key=TRUE)
fit_preds$date = as.Date(fit_preds$DATE_NUM, origin="1970-01-01")
fit_preds$variant = factor(fit_preds$variant, levels=levels_VARIANTS_plot)
levels_state = rev(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","location"][order(fit_preds[fit_preds$date==today&fit_preds$variant=="Omicron (BA.2.75)","prob"])])
as.character(levels_state)
fit_preds$location = factor(fit_preds$location, levels=levels_state)

write.csv(fit_preds, file=".//data//GISAID//GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage_india_bystate_allstateswithatleast_5_BA_2_75.csv", row.names=F)

# to get predictions with confidence intervals (but slow) :
# fit_india_multi_predsbystate = data.frame(emmeans(bestfit_india_multi,
#                                                   ~ variant,
#                                                   by=c("DATE_NUM", 
#                                                        "location"),
#                                                   at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
#                                                   mode="prob", df=NA,
#                                                   rg.limit=100000))
# fit_india_multi_predsbystate$date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
# fit_india_multi_predsbystate$variant = factor(fit_india_multi_predsbystate$variant, levels=levels_VARIANTS_plot)
# fit_india_multi_predsbystate$location = factor(fit_india_multi_predsbystate$location, levels=levels_STATES)


# plot model prediction as line plot on logit scale:

fit_india_multi_preds2 = fit_preds
ymin = 0.0001
ymax = 0.999
if (!is.null(fit_india_multi_preds2$asymp.LCL)) fit_india_multi_preds2$asymp.LCL[fit_india_multi_preds2$asymp.LCL<ymin] = ymin
if (!is.null(fit_india_multi_preds2$asymp.UCL)) fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL<ymin] = ymin
if (!is.null(fit_india_multi_preds2$asymp.UCL)) fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL>ymax] = ymax
fit_india_multi_preds2$prob[fit_india_multi_preds2$prob<ymin] = ymin

plot_india_mfit_logit = qplot(data=fit_india_multi_preds2[fit_india_multi_preds2$variant!="Other",], 
                              x=date, y=prob, geom="blank") +
  facet_wrap(~location, ncol=4) +
  # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
  #                 fill=variant
  # ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
          subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=2)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
  xaxis +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(lineage_cols_plot,-1)) +
  scale_colour_manual("variant", values=tail(lineage_cols_plot,-1)) +
  geom_point(data=data_agbyweekregion1[data_agbyweekregion1$variant!="Other",],
             aes(x=collection_date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.001, 2), limits=c(1,max(data_agbyweekregion1$total)), breaks=c(10,100, 1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today_num+extrapolate+1), 
                  ylim=c(0.001, 0.99), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_india_mfit_logit

ggsave(file=file.path("plots", plotdir,"india_multinom fit_logit scale.png"), width=14, height=8)


# plot model prediction as line plot on response scale:
plot_india_mfit = qplot(data=fit_india_multi_preds2[fit_india_multi_preds2$variant!="Other",], x=date, 
                        y=100*prob, geom="blank") +
  facet_wrap(~location, ncol=4) +
  # geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, 
  #                 colour=NULL,
  #                 fill=variant
  # ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among new infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
          subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=2)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
  xaxis +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),today_num+extrapolate+1), 
                  ylim=c(0, 100)) +
  scale_fill_manual("variant", values=tail(lineage_cols_plot,-1)) +
  scale_colour_manual("variant", values=tail(lineage_cols_plot,-1)) +
  geom_point(data=data_agbyweekregion1[data_agbyweekregion1$variant!="Other",],
             aes(x=collection_date, y=100*prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(0.001, 2), limits=c(1,max(data_agbyweekregion1$total)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_india_mfit

ggsave(file=file.path("plots", plotdir,"india_multinom fit_response scale.png"), width=14, height=8)


# plot predicted values as Muller plot
muller_indiabystate_fit = ggplot(data=fit_preds, 
                                 aes(x=date, y=prob, group=variant)) +
  facet_wrap(~ location, ncol=4) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant, group=variant), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
       subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=2)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
  annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, # extrapolated part
           xmax=as.Date(date.to, origin="1970-01-01")+5, 
           ymin=0, ymax=1, alpha=0.5, fill="white") +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) 
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_indiabystate_fit

ggsave(file=file.path("plots", plotdir,"india_muller plots by state_multinom fit.png"), width=14, height=7)



# project multinomial fit onto case data ####

# changed to new source
# https://github.com/DataKind-BLR/covid19bharat_data/
# https://data.covid19bharat.org/csv/latest/states.csv

cases_india_bystate2 = cases_india_bystate[cases_india_bystate$State %in% 
                                             sel_states,]
colnames(cases_india_bystate2)[2]="STATE"

newdat = expand.grid(DATE_NUM=seq(as.numeric(min(cases_india_bystate2$Date)),
                                  as.numeric(max(cases_india_bystate2$Date))),
                     location=unique(as.character(cases_india_bystate2$STATE)))
fit_india_multi_predsbystate = data.frame(newdat,
                                          predict(bestfit_india_multi, 
                                                  newdata = newdat,
                                                  type = "prob"), check.names=F) 
library(tidyr)
fit_india_multi_predsbystate = gather(fit_india_multi_predsbystate, variant, prob, all_of(levels_VARIANTS))
fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
fit_india_multi_predsbystate$variant = factor(fit_india_multi_predsbystate$variant, levels=levels_VARIANTS)
# colnames(fit_india_multi_predsbystate)[2] = "location"
# fit_india_multi_predsbystate$location = factor(fit_india_multi_predsbystate$location)
fit_india_multi_predsbystate$totnewcases = 
  cases_india_bystate2$newcases[match(interaction(fit_india_multi_predsbystate$location,
                                                  fit_india_multi_predsbystate$collection_date),
                                      interaction(cases_india_bystate2$STATE, 
                                                  cases_india_bystate2$Date))]                                                    
fit_india_multi_predsbystate$totnewcases[which(fit_india_multi_predsbystate$totnewcases>5&fit_india_multi_predsbystate$collection_date>as.Date("2021-08-01")&fit_india_multi_predsbystate$collection_date<as.Date("2021-09-01")&fit_india_multi_predsbystate$location=="Chandigarh")] = 5 # outliers

library(zoo)
fit_india_multi_predsbystate = fit_india_multi_predsbystate %>% 
  group_by(location) %>% 
  mutate(totnewcases_smoothed = rollmean(totnewcases, 7, na.pad = T))
fit_india_multi_predsbystate$totnewtests = 
  cases_india_bystate2$newtests[match(interaction(fit_india_multi_predsbystate$location,
                                                  fit_india_multi_predsbystate$collection_date),
                                      interaction(cases_india_bystate2$STATE, 
                                                  cases_india_bystate2$Date))]                                                                                                                                                                                                                                       
fit_india_multi_predsbystate$totnewtests[fit_india_multi_predsbystate$location=="Haryana"&fit_india_multi_predsbystate$collection_date=="2022-06-23"] = NA # likely error in data
fit_india_multi_predsbystate$totnewtests[fit_india_multi_predsbystate$location=="Odisha"&fit_india_multi_predsbystate$collection_date=="2022-07-04"] = NA # likely error in data
fit_india_multi_predsbystate = fit_india_multi_predsbystate %>% 
  group_by(location) %>% 
  mutate(totnewtests_smoothed = rollmean(totnewtests, 7, na.pad = T))
fit_india_multi_predsbystate$cases = fit_india_multi_predsbystate$totnewcases_smoothed*fit_india_multi_predsbystate$prob
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$cases==0] = NA
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$location=="Chandigarh"&fit_india_multi_predsbystate$cases>2000] = NA
fit_india_multi_predsbystate$cases[fit_india_multi_predsbystate$prob<0.001] = NA
fit_india_multi_predsbystate$posratio = 100*fit_india_multi_predsbystate$cases/fit_india_multi_predsbystate$totnewtests_smoothed
fit_india_multi_predsbystate$posratio[fit_india_multi_predsbystate$posratio<0] = NA
fit_india_multi_predsbystate$posratio[fit_india_multi_predsbystate$posratio>100] = NA
fit_india_multi_predsbystate$posratio[fit_india_multi_predsbystate$variant=="Delta"&fit_india_multi_predsbystate$posratio>10&fit_india_multi_predsbystate$collection_date>as.Date("2021-08-01")&fit_india_multi_predsbystate$collection_date<as.Date("2021-09-01")&fit_india_multi_predsbystate$location=="Chandigarh"] = 0 # outliers

# fit_india_multi_predsbystate$collection_date[which((fit_india_multi_predsbystate$posratio>10)&(fit_india_multi_predsbystate$variant=="Omicron (BA.2.76)")&(fit_india_multi_predsbystate$location=="Odisha"))]
# "2022-07-04" # outlier
fit_india_multi_predsbystate$location = factor(fit_india_multi_predsbystate$location, levels=levels_state)

fit_india_multi_predsbystate2 = fit_india_multi_predsbystate
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases==0] = NA
fit_india_multi_predsbystate2$cases[fit_india_multi_predsbystate2$cases<=1] = NA
fit_india_multi_predsbystate2$location = factor(fit_india_multi_predsbystate2$location, levels=levels_state)
cases_india_bystate2$location = factor(cases_india_bystate2$STATE, levels=levels_state)
fit_india_multi_predsbystate2$variant = factor(fit_india_multi_predsbystate2$variant, levels=levels_VARIANTS_plot)

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_india_multi_predsbystate2, aes(x=collection_date, y=totnewcases_smoothed, lwd=I(1.5)), colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  scale_colour_manual("lineage", values=lineage_cols_plot) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA)) +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_confirmed cases multinomial fit.png"), width=16, height=8)

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))


ggsave(file=file.path("plots", plotdir,"india_confirmed cases stacked area multinomial fit.png"), width=16, height=8)

fit_india_multi_predsbystate3 = fit_india_multi_predsbystate2
fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date<as.Date("2022-05-01"),"cases"] = 0
fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date<as.Date("2022-05-01"),"posratio"] = 0

ggplot(data=fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date>=as.Date("2022-04-30"),], 
       aes(x=collection_date, y=cases, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_confirmed cases stacked area multinomial fit_ZOOMED.png"), width=12, height=6)

# graph with positivity ratios by variant
ggplot(data=fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date>=as.Date("2022-04-30"),], 
       aes(x=collection_date, y=posratio, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("% of tests that were positive") +
  ggtitle("SARS-CoV2 POSITIVITY RATIOS BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase & testing data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_positivity ratios stacked area multinomial fit_ZOOMED.png"), width=12, height=6)
write.csv(fit_india_multi_predsbystate, file="./data/GISAID/GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage_india_bystate__plus_positivity_ratios_by_state_allstateswithatleast_5_BA_2_75.csv", row.names=F)


ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=posratio, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=4) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), position="stack") +
  scale_fill_manual("", values=lineage_cols_plot) +
  # annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
  #          xmax=as.Date("2021-05-31"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("% of tests that were positive") +
  ggtitle("SARS-CoV2 POSITIVITY RATIOS BY VARIANT IN INDIA",
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=2) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase & testing data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=file.path("plots", plotdir,"india_positivity ratios stacked area multinomial fit.png"), width=16, height=8)
