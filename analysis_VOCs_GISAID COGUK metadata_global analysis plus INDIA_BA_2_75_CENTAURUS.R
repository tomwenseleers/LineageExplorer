# ANALYSIS OF GROWTH ADVANTAGE OF BA.2.75, AKA CENTAURUS
# T. Wenseleers
# last update 9 AUGUST 2022

# set GISAID credentials ####
# set them first using 
# Sys.setenv(GISAIDR_USERNAME = "XXXXX")
# Sys.setenv(GISAIDR_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_USERNAME = "XXXXX") # not needed for this script
# Sys.setenv(GISAIDJSON_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_STREAM = "XXXXX")
# for me:

# devtools::install_github("Wytamma/GISAIDR")
library(GISAIDR)
source(".//set_GISAID_credentials.R") # set GISAID credentials
source(".//download_GISAID.R") # load function to download GISAID metadata download package (lacking records from last few days)
source(".//download_GISAID_records.R") # load functions to download most recent GISAID records

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
today # "2021-08-09"
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
target_dir = "C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID" # target download directory
system.time(GISAID <- download_GISAD_meta(target_dir = target_dir,
                                          headless = FALSE,
                                          usr = Sys.getenv("GISAIDR_USERNAME"),
                                          psw = Sys.getenv("GISAIDR_PASSWORD"))) # 214s
# TO DO: check version available for download & if already downloaded don't download it again
# allow downloading FASTA or other available download packages

# note: if you get a message This version of ChromeDriver only supports Chrome version xxx
# then make sure that you have specified the right chromedriver version 
# (check if installed chrome version matches some version available under binman/binman_chromedriver/XXX which ones get installed)
# it may be necessary to downgrade your Chrome browser to version xxx using instructions at
# https://browserhow.com/how-to-downgrade-and-install-older-version-of-chrome/#download-the-older-chrome-version
# download browser install from https://filehippo.com/download_google-chrome/history/
# and disable chrome updates

colnames(GISAID)
# [1] "virus_name"                      "type"                            "accession_id"                   
# [4] "collection_date"                 "location"                        "additional_location_information"
# [7] "sequence_length"                 "host"                            "patient_age"                    
# [10] "gender"                          "clade"                           "pango_lineage"                  
# [13] "pangolin_version"                "variant"                         "aa_substitutions"               
# [16] "submission_date"                 "is_reference"                    "is_complete"                    
# [19] "is_high_coverage"                "is_low_coverage"                 "n_content"                      
# [22] "gc_content" 

download = tail(list.files(target_dir, pattern=".tar.xz"), 1)
download # "metadata_tsv_2022_08_08.tar.xz"  

# records go up to submission date
GISAID_max_submdate = as.Date(max(GISAID$submission_date, na.rm=T))
GISAID_max_submdate # "2022-08-06"
nrow(GISAID) # 12372529

# PS fread is slightly faster (multicore), but requires file to be unzipped first
# system.time(archive_extract(archive=file.path(target_dir,download),
#                 dir=target_dir)) # 33s
# system.time(GISAID <- fread(file.path(target_dir,"metadata.tsv"))) # 51s


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
length(recent_records) # 28369
recent_records = recent_records[!recent_records %in% GISAID$accession_id]
length(recent_records) # 25972

# dataframe with recently submitted records that are not yet in GISAID metadata package download
d_extra = download_GISAID_records(accession_ids = recent_records,
                                  get_sequence=FALSE, 
                                  clean_up=FALSE,
                                  target_dir="C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID/BA_2_75/extra",
                                  max_batch_size=10000, # maximum batch size
                                  headless = FALSE,
                                  usr=Sys.getenv("GISAIDR_USERNAME"),
                                  psw=Sys.getenv("GISAIDR_PASSWORD"))
dim(d_extra) # 25972        17
colnames(d_extra) 
# [1] "virus_name"                      "accession_id"                    "collection_date"                
# [4] "location"                        "host"                            "additional_location_information"
# [7] "sampling_strategy"               "gender"                          "patient_age"                    
# [10] "patient_status"                  "last_vaccinated"                 "passage"                        
# [13] "specimen"                        "additional_host_information"     "pango_lineage"                        
# [16] "clade"                           "aa_substitutions"  

GISAID = dplyr::bind_rows(GISAID, d_extra)
nrow(GISAID) # 12398501

sum(grepl("BA.2.75", GISAID$pango_lineage, fixed=T)) # 1702
sum(grepl("BA.2.75", GISAID$pango_lineage, fixed=T)&grepl("India", GISAID$location)) # 1318


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
# # [1] "strain"                "virus"                 "accession_id"         
# # [4] "genbank_accession"     "date"                  "region"               
# # [7] "country"               "division"              "location"             
# # [10] "region_exposure"       "country_exposure"      "division_exposure"    
# # [13] "segment"               "length"                "host"                 
# # [16] "age"                   "sex"                   "Nextstrain_clade"     
# # [19] "pangolin_lineage"      "GISAID_clade"          "originating_lab"      
# # [22] "submitting_lab"        "authors"               "url"                  
# # [25] "title"                 "paper_url"             "date_submitted"       
# # [28] "purpose_of_sequencing" "sequence"  
  
  
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
coguk = fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz")
colnames(coguk)
colnames(coguk)[which(colnames(coguk) %in% c("sample_date", "lineage", "mutations"))] = c("collection_date", "pango_lineage", "aa_substitutions") # code as in GISAID
coguk$aa_substitutions = gsub(":","_",gsub("S:","Spike_",coguk$aa_substitutions, fixed=T), fixed=T)
coguk$aa_substitutions = gsub("orf1ab_S1221L", "NSP3_S403L",coguk$aa_substitutions, fixed=T)
coguk$aa_substitutions = gsub("orf1ab_N4060S", "NSP8_N118S",coguk$aa_substitutions, fixed=T)

# TO DO: convert syntax to that used in GISAID genomic epidemiology metadata download
# now orf1ab_ still needs to be converted to NSP_ notation, rest should be OK
# see https://github.com/theosanderson/Codon2Nucleotide/blob/main/src/App.js
# https://codon2nucleotide.theo.io/
# (currently only still problems for NSPs <> orf1ab in coguk)
# NSP8_N118S.
# Look up the start of NSP8 in the nsps array: 3943.
# Add 3943+118-1 = 4060 (subtract 1 because they count from 1).
# So NSP8_N118S = orf1ab:N4060S

# NSP3_S403L start of NSP3=819, 819+403-1=1221 -> = orf1ab_S1221L   
# NSP8_N118S start of NSP8=3943, 3943+118-1=4060 -> = orf1ab_N4060S
# I now just converted these 2 mutations to GISAID notation, as I am using these later on
# to help identify BA.2.75

coguk$host = "Human"
coguk$country = "United Kingdom"
coguk$continent = "Europe"
coguk$collection_date = as.character(coguk$collection_date)
nrow(coguk) # 2818167

# MERGE GISAID (MINUS UK GISAID DATA) & COG-UK DATA FOR UK ####
GISAID = dplyr::bind_rows(GISAID[GISAID$country!="United Kingdom",], 
                           coguk)
nrow(GISAID) # 12320979

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
levels_VARIANTS = c(sel_reference_VOC, "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.2.76)", "Omicron (BA.5.2)", sel_target_VOC)
levels_VARIANTS_plot = c("Other", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.2.12.1)", "Omicron (BA.4)", "Omicron (BA.4.6)", "Omicron (BA.5)", "Omicron (BA.5.2)", "Omicron (BA.2.76)", "Omicron (BA.2.75)")

n = length(levels_VARIANTS)
lineage_cols_plot = case_when(
  levels_VARIANTS_plot=="Other" ~ "grey65",
  levels_VARIANTS_plot=="Beta" ~ "green4",
  levels_VARIANTS_plot=="Alpha" ~ "#0085FF",
  levels_VARIANTS_plot=="Delta" ~ "mediumorchid",
  levels_VARIANTS_plot=="Omicron (BA.1)" ~ "red",
  levels_VARIANTS_plot=="Omicron (BA.2)" ~ "red3",
  levels_VARIANTS_plot=="Omicron (BA.2.12.1)" ~ "red4",
  levels_VARIANTS_plot=="Omicron (BA.4)" ~ "green3",
  levels_VARIANTS_plot=="Omicron (BA.4.6)" ~ "green2",
  levels_VARIANTS_plot=="Omicron (BA.5)" ~ "blue4",
  levels_VARIANTS_plot=="Omicron (BA.5.2)" ~ "blue2",
  levels_VARIANTS_plot=="Omicron (BA.2.76)" ~ "magenta4",
  levels_VARIANTS_plot=="Omicron (BA.2.75)" ~ "magenta",
)
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
  # (grepl("BA.2.74", GISAID$pango_lineage, fixed=T)|(grepl("Spike_L452M",GISAID$aa_substitutions)&
  #    grepl("Spike_R346T",GISAID$aa_substitutions))) ~ "Omicron (BA.2.74)",
   (grepl("BA.2.76", GISAID$pango_lineage, fixed=T)|(grepl("Spike_Y248N",GISAID$aa_substitutions)&
      grepl("Spike_R346T",GISAID$aa_substitutions))) ~ "Omicron (BA.2.76)", 
  grepl("^BA\\.4\\.6$", GISAID$pango_lineage) ~ "Omicron (BA.4.6)",
  grepl("^BA\\.5\\.2$|^BA\\.5\\.2\\.|BF", GISAID$pango_lineage)&(GISAID$date>as.Date("2022-02-01")) ~ "Omicron (BA.5.2)",
  # grepl("^BA\\.5\\.3$|^BA\\.5\\.3\\.|BE", GISAID$pango_lineage) ~ "Omicron (BA.5.3)",
  grepl("B.1.617.2", GISAID$pango_lineage, fixed=T) | grepl("AY", GISAID$pango_lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", GISAID$pango_lineage) ~ "Alpha",
  grepl("B.1.351", GISAID$pango_lineage, fixed=T) ~ "Beta",
  (grepl("^BA\\.1$|BA\\.1\\.", GISAID$pango_lineage)&(GISAID$date>as.Date("2021-11-01"))) ~ "Omicron (BA.1)",
  # (grepl("^BA\\.3$|BA\\.3\\.", GISAID$pango_lineage)) ~ "Omicron (BA.3)",
  (((grepl("^BA\\.4",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
                                                ((grepl("L452R", GISAID$aa_substitutions)&
                                                    grepl("486V", GISAID$aa_substitutions)&
                                                    grepl("11F", GISAID$aa_substitutions)&
                                                    (!grepl("D3N",GISAID$aa_substitutions)) )))) ~ "Omicron (BA.4)",
  (grepl("^BA\\.5",GISAID$pango_lineage)|grepl("BE|BF",GISAID$pango_lineage)|((GISAID$pango_lineage!="Unassigned")&grepl("M_D3N",GISAID$aa_substitutions)))&(GISAID$date>as.Date("2022-02-01")) ~ "Omicron (BA.5)", # cf pattern used by Alex Selby
  (grepl("^BA\\.2\\.12\\.1$", GISAID$pango_lineage))&(GISAID$date>as.Date("2022-02-01")) ~ "Omicron (BA.2.12.1)",
  (grepl("^BA\\.2",GISAID$pango_lineage)) ~ "Omicron (BA.2)",
  GISAID$pango_lineage!="Unassigned" ~ "Other" # assigns NA to remaining Unassigned & remove them later on
  # TRUE ~ "Other"
)


sum(is.na(GISAID$variant)) # 189574 unassigned
sum(GISAID$pango_lineage=="Unassigned", na.rm=T) # 189615 originally unassigned in GISAID
sum(GISAID$variant=="Omicron (BA.2.75)", na.rm=T) # 2111 BA.2.75
sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$country=="India", na.rm=T) # 1720 BA.2.75 for India
sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$date_isvalid, na.rm=T) # 1912 BA.2.75 with valid date
sum(GISAID$variant=="Omicron (BA.2.75)"&GISAID$country=="India"&GISAID$date_isvalid, na.rm=T) # 1522 BA.2.75 for India with valid date
sum(GISAID$pango_lineage=="BA.2.75", na.rm=T) # 1701
sum(GISAID$variant=="Omicron (BA.5)", na.rm=T) # 242018
table(GISAID$variant)

# GISAID = GISAID[!is.na(GISAID$variant),]
nrow(GISAID) # 12320979



# GISAID SELECTION ####
GISAID_sel = GISAID

# remove records with invalid/incomplete dates ####
GISAID_sel = GISAID_sel[GISAID_sel$date_isvalid,]
GISAID_sel = GISAID_sel[which(GISAID_sel$host=="Human"),]
nrow(GISAID_sel) # 12104208

# filter to desired date range ####
start_date = "2020-06-01"
end_date = today
GISAID_sel = GISAID_sel[GISAID_sel$date>=as.Date(start_date)&
                          GISAID_sel$date<=as.Date(end_date),]

range(GISAID_sel$date, na.rm=T) # "2020-06-01" "2022-08-07"
nrow(GISAID_sel) # 11958660

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
nBA_2_75 # 1912 BA.2.75 records so far with valid date

nBA_2_75_india = sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="India", na.rm=T)
nBA_2_75_india # 1522 BA.2.75 records so far for India with valid date

sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="United Kingdom", na.rm=T)
# 33 BA.2.75 in UK so far
maxsubmdate = today


# 2. ANALYSIS OF GLOBAL BA.2.75 GROWTH RATE ADVANTAGE ####

# TO DO: replace this with more elegant/tidy dplyr coce? 

# we first subset to countries where at least 1 BA.2.75 was detected
tab = as.data.frame(table(GISAID_sel$country, GISAID_sel$variant))

sel_countries_min10 = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&tab$Freq>=10,"Var1"])
sel_countries_min10
# [1] "Australia"      "Canada"         "Denmark"        "India"          "Israel"        
# [6] "Japan"          "Nepal"          "Singapore"      "United Kingdom" "USA"  

sel_countries = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&tab$Freq>=1,"Var1"])
sel_countries

# GISAID selection : subsetted to countries with >=1 BA.2.75 sequence ####
GISAID_sel = GISAID_sel[as.character(GISAID_sel$country) %in% sel_countries,]
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$variant),] # we remove records with unassigned lineage
GISAID_sel$country = factor(GISAID_sel$country)
table(GISAID_sel$country, GISAID_sel$variant)
table(GISAID_sel$continent, GISAID_sel$variant)

df_cont=as.data.frame(table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","continent"]))
colnames(df_cont)=c("continent","BA.2.75")
df_cont[order(df_cont$BA.2.75,decreasing=T),]
#       continent BA.2.75
# 1          Asia    1664
# 3 North America     122
# 2        Europe      81
# 4       Oceania      42
# 5 South America       3

df=as.data.frame(table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","location"]))
df=df[df$Freq!=0,]
df=df[order(df$Freq, decreasing=T),]
df

df=as.data.frame(table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","country"]))
df=df[df$Freq!=0,]
df=df[order(df$Freq, decreasing=T),]
colnames(df)=c("country","BA.2.75")
df
#           country BA.2.75
# 10          India    1522
# 28            USA      87
# 14          Japan      53
# 21      Singapore      43
# 1       Australia      36
# 4          Canada      34
# 27 United Kingdom      33
# 17          Nepal      17
# 7         Denmark      16
# 12         Israel      11
# 2         Austria       9
# 11      Indonesia       9
# 8          France       7
# 9         Germany       6
# 19    New Zealand       6
# 24    South Korea       5
# 15     Luxembourg       4
# 13          Italy       2
# 20           Peru       2
# 25       Thailand       2
# 3        Cambodia       1
# 5           Chile       1
# 6           China       1
# 16     Martinique       1
# 18    Netherlands       1
# 22       Slovakia       1
# 23       Slovenia       1
# 26         Turkey       1


# AGGREGATE DATA BY DATE & COUNTRY ####
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
write.csv(data_agbydatecountry1, file=".//data//GISAID//GISAID aggregated counts by date and lineage.csv", row.names=F)

# AGGREGATE DATA BY WEEK & COUNTRY ####
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
write.csv(data_agbyweekcountry1, file="./data/GISAID/GISAID aggregated counts by start of week and lineage.csv", row.names=F)

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
                at=list(DATE_NUM=max(GISAID_sel$DATE_NUM, na.rm=T)))
delta_r_pairw = data.frame(confint(emtr_pairw, 
                             adjust="none", df=NA)$contrasts, 
                     p.value=as.data.frame(emtr_pairw$contrasts)$p.value)
delta_r_pairw
write.csv(delta_r_pairw, file.path("plots", plotdir, "growth advantage BA.2.75 vs all other strains by continent.csv"), row.names=F)

# plot growth rate advantage of BA.2.75
delta_r_pairw2 = delta_r_pairw[delta_r_pairw$contrast %in%
                                 c("Omicron (BA.2.75) - Omicron (BA.2)",
                                   "Omicron (BA.2.75) - Omicron (BA.2.12.1)",
                                   "Omicron (BA.2.75) - Omicron (BA.5)",
                                   "Omicron (BA.2.75) - Omicron (BA.5.2)",
                                   "Omicron (BA.2.75) - Omicron (BA.4)",
                                   "Omicron (BA.2.75) - Omicron (BA.4.6)"),]
delta_r_pairw2$contrast = droplevels(delta_r_pairw2$contrast)
delta_r_pairw2$contrast = factor(delta_r_pairw2$contrast,
                                 levels=c("Omicron (BA.2.75) - Omicron (BA.2)",
                                          "Omicron (BA.2.75) - Omicron (BA.2.12.1)",
                                          "Omicron (BA.2.75) - Omicron (BA.4)",
                                          "Omicron (BA.2.75) - Omicron (BA.5)",
                                          "Omicron (BA.2.75) - Omicron (BA.5.2)",
                                          "Omicron (BA.2.75) - Omicron (BA.4.6)"))
delta_r_pairw2 = delta_r_pairw2[!delta_r_pairw2$continent %in%
                                  c("South America", "Africa"),]

library(RColorBrewer)
cols = brewer.pal(n = 7, "Blues")[-1]
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
                      at=list(DATE_NUM=max(GISAID_sel$DATE_NUM, na.rm=T)))
delta_r_BA_4_6 = data.frame(confint(emtr_pairw_BA_4_6, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtr_pairw_BA_4_6$contrasts)$p.value)

delta_r_BA_4_6_2 = delta_r_BA_4_6[delta_r_BA_4_6$contrast %in% 
                                    c("Omicron (BA.4.6) - Omicron (BA.5)",
                                      "Omicron (BA.4.6) - Omicron (BA.5.2)",
                                      "Omicron (BA.4.6) - Omicron (BA.4)",
                                      "Omicron (BA.4.6) - Omicron (BA.2.75)",
                                      "Omicron (BA.4.6) - Omicron (BA.2.12.1)",
                                      "Omicron (BA.4.6) - Omicron (BA.2)"
                                      ),]
delta_r_BA_4_6_2$contrast = droplevels(delta_r_BA_4_6_2$contrast)
delta_r_BA_4_6_2$contrast = factor(delta_r_BA_4_6_2$contrast,
                                 levels=c("Omicron (BA.4.6) - Omicron (BA.2)",
                                          "Omicron (BA.4.6) - Omicron (BA.2.12.1)",
                                          "Omicron (BA.4.6) - Omicron (BA.4)",
                                          "Omicron (BA.4.6) - Omicron (BA.5)",
                                          "Omicron (BA.4.6) - Omicron (BA.5.2)",
                                          "Omicron (BA.4.6) - Omicron (BA.2.75)"))
delta_r_BA_4_6_2 = delta_r_BA_4_6_2[!delta_r_BA_4_6_2$continent %in%
                                  c("South America", "Africa"),]

library(RColorBrewer)
cols = brewer.pal(n = 7, "Blues")[-1]
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
date.from = as.numeric(as.Date("2020-06-01"))
date.to = today_num+extrapolate # max(GISAID_sel$DATE_NUM, na.rm=T)

# plot predicted values
# multinomial model predictions by province (fastest, but no confidence intervals)
predgrid = expand.grid(list(DATE_NUM=seq(date.from, date.to),
                            country=levels(GISAID_sel$country)))
predgrid$continent = GISAID_sel$continent[match(predgrid$country,
                                                GISAID_sel$country)]

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
write.csv(fit_preds, file=".//data//GISAID//GISAID fitted lineage frequencies multinomial spline fit by start of week and lineage.csv", row.names=F)


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
  coord_cartesian(xlim=c(as.Date("2020-11-01"),NA), 
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




# 3. ANALYSIS OF VOC LINEAGE FREQUENCIES IN INDIA ####

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
nBA_2_75_india = sum(GISAID_india$variant=="Omicron (BA.2.75)")
nBA_2_75_india # 1521

# # PS this was for Ewen Callaway's Nature piece on BA.2.75
# # aggregated by week for selected variant lineages for whole of India using data from all states
# data_agbyweek1 = as.data.frame(table(GISAID_india$floor_date, 
#                                      GISAID_india$variant))
# colnames(data_agbyweek1) = c("floor_date", "variant", "count")
# data_agbyweek1_sum = aggregate(count ~ floor_date, 
#                                data=data_agbyweek1, sum)
# data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$floor_date, data_agbyweek1_sum$floor_date)]
# sum(data_agbyweek1[data_agbyweek1$variant=="Omicron (BA.1)","total"]) == nrow(GISAID_india) # correct
# data_agbyweek1$collection_date = as.Date(as.character(data_agbyweek1$floor_date))
# data_agbyweek1$variant = factor(data_agbyweek1$variant, 
#                                 levels=levels_VARIANTS)
# data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
# data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total
# data_agbyweek1$floor_date = NULL
# 
# # MULLER PLOT (RAW DATA)
# data_agbyweek1$variant2 = factor(data_agbyweek1$variant, levels=levels_VARIANTS_plot)
# muller_india_raw1_all = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=variant2)) + 
#   # facet_wrap(~ STATE, ncol=1) +
#   # geom_col(aes(lwd=I(1.2), colour=NULL, fill=variant), width=1, position="fill") +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant2, group=variant2), position="fill") +
#   scale_fill_manual("", values=lineage_cols_plot) +
#   xaxis +
#   # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
#   theme_hc() +
#   ylab("Share") + 
#   theme(legend.position="right",  
#         axis.title.x=element_blank()) +
#   labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
#        subtitle=paste0("Raw GISAID data up to ",today,"\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)")) 
# # +
# # coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
# muller_india_raw1_all
# 
# ggsave(file=file.path("plots", plotdir,"india_muller plots_raw data all states.png"), width=7, height=5)
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

nBA_2_75_india = sum(GISAID_india$variant=="Omicron (BA.2.75)")
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
BIC(fit1_india_multi, fit2_india_multi, fit3_india_multi, fit4_india_multi, fit5_india_multi) 
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
