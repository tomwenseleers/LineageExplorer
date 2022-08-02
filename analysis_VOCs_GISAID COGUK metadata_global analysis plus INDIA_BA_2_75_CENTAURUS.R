# ANALYSIS OF GROWTH ADVANTAGE OF BA.2.75, AKA CENTAURUS
# T. Wenseleers
# last update 1 AUGUST 2022

# set GISAID credentials ####
# set them first using 
# Sys.setenv(GISAIDR_USERNAME = "XXXXX")
# Sys.setenv(GISAIDR_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_USERNAME = "XXXXX") # not needed for this script
# Sys.setenv(GISAIDJSON_PASSWORD = "XXXXX")
# Sys.setenv(GISAIDJSON_STREAM = "XXXXX")
# for me:
source(".//set_GISAID_credentials.R")
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
# devtools::install_github("Wytamma/GISAIDR")
library(GISAIDR)

# 1. LOAD DATA ####

today = as.Date(Sys.time())
today_num = as.numeric(today)
today # "2021-07-29"
plotdir = "GISAID_BA_2_75"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))
tag = paste("@TWenseleers\n",today)

# X axis for plots
firststofmonth = seq(as.Date("2020-01-01"), as.Date("2022-12-01"), by="month")
xaxis = scale_x_continuous(breaks=firststofmonth,
                           labels=substring(months(firststofmonth),1,1),
                           expand=c(0,0))

# import GISAID metadata ####

# download latest GISAID metadata ####
target_dir = "C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID" # target download directory
system.time(GISAID <- download_GISAD_meta(target_dir = target_dir)) # 150s
download = tail(list.files(target_dir, pattern=".tar.xz"), 1)
download # "metadata_tsv_2022_08_01.tar.xz"  

# records go up to submission date
GISAID_max_submdate = as.Date(max(GISAID$`Submission date`, na.rm=T))
GISAID_max_submdate # "2022-07-30"
nrow(GISAID) # 12242947

# PS fread is slightly faster (multicore), but requires file to be unzipped first
# system.time(archive_extract(archive=paste0(target_dir, "//", download),
#                 dir=target_dir)) # 33s
# system.time(GISAID <- fread(paste0(target_dir, "//", "metadata.tsv"))) # 51s


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
length(recent_records) # 23031
recent_records = recent_records[!recent_records %in% GISAID$`Accession ID`]
length(recent_records) # 20378

write.csv(recent_records[1:10000],file=".//data//GISAID//BA_2_75//extra//ids1.csv", row.names=F)
write.csv(recent_records[10001:20000],file=".//data//GISAID//BA_2_75//extra//ids2.csv", row.names=F)
write.csv(recent_records[20001:length(recent_records)],file=".//data//GISAID//BA_2_75//extra//ids3.csv", row.names=F)
d1 = read_tsv(".//data//GISAID//BA_2_75//extra//gisaid_hcov-19_2022_08_01_14.tsv", col_types = cols(.default = "c"))
d2 = read_tsv(".//data//GISAID//BA_2_75//extra//gisaid_hcov-19_2022_08_01_14 (1).tsv", col_types = cols(.default = "c"))
d3 = read_tsv(".//data//GISAID//BA_2_75//extra//gisaid_hcov-19_2022_08_01_14 (2).tsv", col_types = cols(.default = "c"))

# normally this should have worked, but AA field missing, so still did it manually as above
# extra_records = download_GISAID_records(
#     accession_ids=recent_records,
#     get_sequence=FALSE, 
#     clean_up=FALSE,
#     target_dir="C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID/BA_2_75/extra",
#     max_batch_size=10000, # maximum batch size
#     usr=Sys.getenv("GISAIDR_USERNAME"),
#     psw=Sys.getenv("GISAIDR_PASSWORD"))

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
# dim(downloads) # 103356     29
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


d_extra = rbind(d1, d2, d3)
d_extra = d_extra[!d_extra$`Accession ID` %in% GISAID$`Accession ID`,] # remove duplicates if present
d_extra$duplicated = duplicated(d_extra$`Accession ID`)
d_extra = d_extra[!d_extra$duplicated,]
d_extra$duplicated = NULL
nrow(d_extra) # 20378

GISAID = dplyr::bind_rows(GISAID, d_extra)
nrow(GISAID) # 12263325

# parse location field ####
# parse continent / country / location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID$Location), "/", TRUE)) # parse location info
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
coguk$`Collection date` = as.character(coguk$sample_date)
coguk$`Pango lineage` = coguk$lineage
coguk$aa_substitutions = "" 
# TO DO: convert syntax to that used in GISAID genomic epidemiology metadata download
# see https://github.com/theosanderson/Codon2Nucleotide/blob/main/src/App.js
# https://codon2nucleotide.theo.io/
coguk$Host = "Human"
coguk$country = "United Kingdom"
coguk$continent = "Europe"
nrow(coguk) # 2802865

# MERGE GISAID (MINUS UK GISAID DATA) & COG-UK DATA FOR UK ####
GISAID = dplyr::bind_rows(GISAID[GISAID$country!="United Kingdom",], 
                           coguk)
nrow(GISAID) # 12186865

# PARSE GISAID DATA ####
GISAID = as.data.frame(GISAID) 
# TO DO: keep as more efficient data.table or as tibble
# & make sure everything below still works

# define variant lineages and colours ####
sel_target_VOC = "Omicron (BA.2.75)"
sel_reference_VOC = "Omicron (BA.5)"
levels_VARIANTS = c(sel_reference_VOC, "Beta", "Alpha", "Other", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", "Omicron (BA.2.74)", "Omicron (BA.2.76)", sel_target_VOC)
levels_VARIANTS_plot = c("Other", "Beta", "Alpha", "Delta", "Omicron (BA.1)", "Omicron (BA.2)", "Omicron (BA.3)", "Omicron (BA.4)", "Omicron (BA.5)", "Omicron (BA.2.74)", "Omicron (BA.2.76)", "Omicron (BA.2.75)")

n = length(levels_VARIANTS)
lineage_cols = hcl(h = seq(0, 180, length = n), l = 60, c = 180)
lineage_cols[which(levels_VARIANTS=="Alpha")] = "#0085FF"
lineage_cols[which(levels_VARIANTS=="Beta")] = "green4"
lineage_cols[which(levels_VARIANTS=="Delta")] = "mediumorchid"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.1)")] = "red" # "magenta"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2)")] = "red3"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.3)")] = "red4" 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.4)")] = "darkorange" 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.5)")] = "darkorange3" 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2.74)")] = "black"
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2.76)")] = muted(muted("magenta")) 
lineage_cols[which(levels_VARIANTS=="Omicron (BA.2.75)")] = "magenta" 
lineage_cols[which(levels_VARIANTS=="Other")] = "grey65"

lineage_cols_plot = lineage_cols[match(levels_VARIANTS_plot,levels_VARIANTS)]

# CODE MAIN LINEAGES ####
GISAID$pango_lineage = GISAID$`Pango lineage`
GISAID$aa_substitutions = GISAID$`AA Substitutions`
GISAID$aa_substitutions[is.na(GISAID$aa_substitutions)] = ""

# PS this is slow - optimize this, maybe use multiple cores / multidplyr?
# or avoid having to use all these grepl patterns
GISAID$variant = case_when(
  (grepl("BA.2.75", GISAID$pango_lineage, fixed=T)|(grepl("NSP3_S403L",GISAID$aa_substitutions)& 
                                                      grepl("NSP8_N118S",GISAID$aa_substitutions))|
     (grepl("NSP3_S403L",GISAID$aa_substitutions)& 
        grepl("E_T11A",GISAID$aa_substitutions))) ~ "Omicron (BA.2.75)",
  (grepl("BA.2.74", GISAID$pango_lineage, fixed=T)|(grepl("Spike_L452M",GISAID$aa_substitutions)&
      grepl("Spike_R346T",GISAID$aa_substitutions))) ~ "Omicron (BA.2.74)",
  (grepl("BA.2.76", GISAID$pango_lineage, fixed=T)|(grepl("Spike_Y248N",GISAID$aa_substitutions)&
      grepl("Spike_R346T",GISAID$aa_substitutions))) ~ "Omicron (BA.2.76)",  
  grepl("B.1.617.2", GISAID$pango_lineage, fixed=T) | grepl("AY", GISAID$pango_lineage)  ~ "Delta",
  grepl("^B\\.1\\.1\\.7$", GISAID$pango_lineage) ~ "Alpha",
  grepl("B.1.351", GISAID$pango_lineage, fixed=T) ~ "Beta",
  (grepl("^BA\\.1$|BA\\.1\\.", GISAID$pango_lineage)) ~ "Omicron (BA.1)",
  (grepl("^BA\\.3$|BA\\.3\\.", GISAID$pango_lineage)) ~ "Omicron (BA.3)",
  (((grepl("^BA\\.4",GISAID$pango_lineage)))|((grepl("BA.2",GISAID$pango_lineage))&
                                                ((grepl("L452R", GISAID$aa_substitutions)&
                                                    grepl("486V", GISAID$aa_substitutions)&
                                                    grepl("11F", GISAID$aa_substitutions)&
                                                    (!grepl("D3N",GISAID$aa_substitutions)) )))) ~ "Omicron (BA.4)",
  # (((grepl("^BA\\.5",GISAID$pango_lineage)|(grepl("BE|BF",GISAID$pango_lineage))))|((grepl("BA.2",GISAID$pango_lineage))&
  #                                               ((grepl("L452R",GISAID$aa_substitutions)&
  #                                                   grepl("486V",GISAID$aa_substitutions)&
  #                                                   (!grepl("11F", GISAID$aa_substitutions))&
  #                                                   grepl("D3N",GISAID$aa_substitutions))))) ~ "Omicron (BA.5)", # old pattern I was using
  (grepl("^BA\\.5",GISAID$pango_lineage)|grepl("BE|BF",GISAID$pango_lineage)|((GISAID$pango_lineage!="Unassigned")&grepl("M_D3N",GISAID$aa_substitutions))) ~ "Omicron (BA.5)", # cf pattern used by Alex Selby
  (grepl("^BA\\.2",GISAID$pango_lineage)) ~ "Omicron (BA.2)",
  GISAID$pango_lineage!="Unassigned" ~ "Other" # assign NA to remaining Unassigned & remove them
  # TRUE ~ "Other"
)

sum(is.na(GISAID$variant)) # 196644 unassigned, in GISAID somewhat less - check where discrepancy comes from
sum(GISAID$pango_lineage=="Unassigned", na.rm=T) # 182498
sum(GISAID$variant=="Omicron (BA.2.75)", na.rm=T) # 1312 BA.2.75 
sum(GISAID$pango_lineage=="BA.2.75", na.rm=T) # 808
sum(GISAID$variant=="Omicron (BA.5)", na.rm=T) # 318579
table(GISAID$variant)

# GISAID = GISAID[!is.na(GISAID$variant),]
nrow(GISAID) # 12186865


# FILTER OUT SOME INVALID GISAID RECORDS ####
GISAID_sel = GISAID

# remove records with invalid/incomplete dates ####
GISAID_sel$date = GISAID_sel$`Collection date` 
date_isvalid = (str_count(GISAID_sel$date,
                          pattern = "-")==2)
GISAID_sel = GISAID_sel[which(date_isvalid),]
GISAID_sel$date = as.Date(fast_strptime(GISAID_sel$date, "%Y-%m-%d")) # faster than as.Date(GISAID$date)

# keep only Human records ####
# unique(GISAID$Host)
# GISAID[GISAID$Host!="Human","Lineage"]
GISAID_sel = GISAID_sel[which(GISAID_sel$Host=="Human"),]
nrow(GISAID_sel) # 11971603

# filter to desired date range ####
start_date = "2020-06-01"
end_date = today
GISAID_sel = GISAID_sel[GISAID_sel$date>=as.Date(start_date)&
                          GISAID_sel$date<=as.Date(end_date),]

range(GISAID_sel$date, na.rm=T) # "2020-06-01" "2022-07-28"
nrow(GISAID_sel) # 11826153

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
nBA_2_75 # 1118 BA.2.75 so far
sum(GISAID_sel$variant=="Omicron (BA.2.75)"&GISAID_sel$country=="United Kingdom", na.rm=T) # 26
# 30 BA.2.75 in UK so far
maxsubmdate = today


# 2. ANALYSIS OF GLOBAL BA.2.75 GROWTH RATE ADVANTAGE ####

# TO DO: replace this with more elegant/tidy dplyr coce? 

# we first subset to countries where at least 1 BA.2.75 was detected
tab = as.data.frame(table(GISAID_sel$country, GISAID_sel$variant))

sel_countries_min10 = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&tab$Freq>=10,"Var1"])
sel_countries_min10
# [1] "Australia"      "Canada"         "Denmark"        "India"          "Japan"          "Nepal"         
# [7] "Singapore"      "United Kingdom" "USA"

sel_countries = as.character(tab[tab$Var2=="Omicron (BA.2.75)"&tab$Freq>=1,"Var1"])
sel_countries

# GISAID selection : subsetted to countries with >=BA.2.75 sequence ####
GISAID_sel = GISAID_sel[as.character(GISAID_sel$country) %in% sel_countries,]
GISAID_sel = GISAID_sel[!is.na(GISAID_sel$variant),]
GISAID_sel$country = factor(GISAID_sel$country)
GISAID_sel$country = droplevels(GISAID_sel$country)
table(GISAID_sel$country, GISAID_sel$variant)
table(GISAID_sel$continent, GISAID_sel$variant)

df_cont=as.data.frame(table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","continent"]))
colnames(df_cont)=c("continent","BA.2.75")
df_cont[order(df_cont$BA.2.75,decreasing=T),]

df=as.data.frame(table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","location"]))
df=df[df$Freq!=0,]
df=df[order(df$Freq, decreasing=T),]
df

df=as.data.frame(table(GISAID_sel[GISAID_sel$variant=="Omicron (BA.2.75)","country"]))
df=df[df$Freq!=0,]
df=df[order(df$Freq, decreasing=T),]
colnames(df)=c("country","BA.2.75")
df

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
# data_agbyweekcountry1 = read.csv(file=".//data//GISAID//GISAID aggregated counts by week and lineage.csv") 
data_agbyweekcountry1$country = factor(data_agbyweekcountry1$country)
if (is.null(data_agbyweekcountry1$continent)) data_agbyweekcountry1$continent = GISAID_sel$continent[match(data_agbyweekcountry1$country, GISAID_sel$country)]
data_agbyweekcountry1$continent = factor(data_agbyweekcountry1$continent)
data_agbyweekcountry1$collection_date = as.Date(data_agbyweekcountry1$collection_date)
data_agbyweekcountry1$variant = factor(data_agbyweekcountry1$variant, levels=levels_VARIANTS)
write.csv(data_agbyweekcountry1, file=".//data//GISAID//GISAID aggregated counts by week and lineage.csv", row.names=F)

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

# best models:
# fit8_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*continent+country, weights=count, data=data_agbyweekcountry1, maxit=10000)
#fit8b_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+ns(DATE_NUM, df=2):continent+country, weights=count, data=data_agbyweekcountry1, maxit=10000)
fit9_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)*continent+country, weights=count, data=data_agbyweekcountry1, maxit=10000)
fit9b_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)+ns(DATE_NUM, df=3):continent+country, weights=count, data=data_agbyweekcountry1, maxit=10000)

#fit10_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*country, weights=count, data=data_agbyweekcountry1, maxit=1000)
#fit11_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2), weights=count, data=data_agbyweekcountry1, maxit=1000)

BIC(fit8_multi, fit8b_multi,fit9_multi, fit9b_multi) 
# fit9_multi has best BIC, fit9b_multi & fit8b_multi close

bestfit_multi = fit9b_multi # I will use this model

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
# write.csv(delta_r, paste0(".\\plots\\",plotdir,"\\growth advantage BA.2.75 vs BA.5 by continent.csv"), row.names=F)

emtr_pairw = emtrends(bestfit_multi, revpairwise ~ variant, by="continent", 
                var="DATE_NUM",  mode="latent",
                at=list(DATE_NUM=max(GISAID_sel$DATE_NUM, na.rm=T)))
delta_r_pairw = data.frame(confint(emtr_pairw, 
                             adjust="none", df=NA)$contrasts, 
                     p.value=as.data.frame(emtr_pairw$contrasts)$p.value)
delta_r_pairw
write.csv(delta_r_pairw, paste0(".\\plots\\",plotdir,"\\growth advantage BA.2.75 vs all other strains by continent.csv"), row.names=F)

delta_r_pairw2 = delta_r_pairw[delta_r_pairw$contrast %in%
                                 c("Omicron (BA.2.75) - Omicron (BA.2)",
                                   "Omicron (BA.2.75) - Omicron (BA.5)"),]
delta_r_pairw2$contrast = droplevels(delta_r_pairw2$contrast)
delta_r_pairw2$contrast = factor(delta_r_pairw2$contrast,
                                 rev(levels(delta_r_pairw$contrast)))
delta_r_pairw2 = delta_r_pairw2[!delta_r_pairw2$continent %in%
                                  c("South America"),]

qplot(data=delta_r_pairw2,
      x=continent, group=contrast,
      y=estimate,
      ymin=asymp.LCL,
      ymax=asymp.UCL,
      geom="blank") +
  geom_col(position=position_dodge(), aes(fill=contrast)) +
  geom_linerange(position=position_dodge(.9)) +
  scale_fill_manual("", values=c("red","blue")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("Difference in growth rate (/day)") +
  ggtitle("GROWTH RATE ADVANTAGE OF BA.2.75 SARS-CoV2 VARIANT",
          subtitle=paste0("based on multinomial spline fit variant ~ ns(date, df=3)+ns(date, df=3):continent+country,\nusing GISAID data with submission data up to ", maxsubmdate, " &\ndata from countries with at least 1 BA.2.75 sequence plus COG-UK data\n(n=", nBA_2_75, " BA.2.75 sequences in total)")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
ggsave(file=paste0(".\\plots\\",plotdir,"\\growth rate advantage BA_2_75_by continent.png"), width=7, height=5)



# PLOT MULTINOMIAL FIT

extrapolate = 30
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

# PLOT OF MULTINOMIAL FIT ON LOGIT SCALE ####
ymin = 0.0001
ymax = 0.999
if (!is.null(fit_preds$asymp.LCL)) fit_preds$asymp.LCL[fit_preds$asymp.LCL<ymin] = ymin
if (!is.null(fit_preds$asymp.UCL)) fit_preds$asymp.UCL[fit_preds$asymp.UCL<ymin] = ymin
if (!is.null(fit_preds$asymp.UCL)) fit_preds$asymp.UCL[fit_preds$asymp.UCL>ymax] = ymax
fit_preds$prob[fit_preds$prob<ymin] = ymin
fit_preds$prob[fit_preds$prob>ymax] = ymax

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
          subtitle=paste0("GISAID data up to ",today, " plus COG-UK data, multinomial spline fit variant ~ ns(date, df=3)+ns(date, df=3):continent+country,\nall countries with >=10 BA.2.75 sequences shown (n=", nBA_2_75," BA.2.75 sequences in total)")) +
  xaxis +  
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=tail(lineage_cols_plot,-1)) +
  scale_colour_manual("variant", values=tail(lineage_cols_plot,-1)) +
  geom_point(data=data_agbyweekcountry1[data_agbyweekcountry1$variant!="Other"&data_agbyweekcountry1$country %in% sel_countries_min10,],
             aes(x=collection_date, y=prop, size=total,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(1, 5), limits=c(1,max(data_agbyweekcountry1$total)), 
                        breaks=c(10,100,1000, 10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date (start of week)") +
  coord_cartesian(xlim=c(as.Date("2022-01-01"),NA), 
                  ylim=c(0.0001, 0.99901), expand=0) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_preds_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\predictions global multinom fit_logit scale.png"), width=11, height=8)

# map variant share onto case numbers

sel_countries_top3 = c("India", "Nepal", "Japan")
country_data = get_national_data(countries=sel_countries_top3,
                                 source="who",
                                 level=1)
country_data$country[country_data$country=="United States"] = "USA"
country_data$country = factor(country_data$country, levels=levels(fit_preds$country))
qplot(data=country_data, x=date, y=cases_new, group=country, geom="blank", colour=country) +
  facet_wrap(~country, scale="free_y") +
  geom_ma(ma_fun = SMA, n = 7, lty=1) +
  scale_colour_hue("") +
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
ggsave(file=paste0(".\\plots\\",plotdir,"\\new cases in top 3 countries with most BA_2_75 cases.png"), width=7, height=5)

fit_preds[(fit_preds$date==today)&(fit_preds$variant=="Omicron (BA.2.75)"),]
# 81% [] in Nepal, 65% [] in India, 6.5% [] in Japan 

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
fit_preds = fit_preds[fit_preds$country %in% sel_countries_top3,]
fit_preds$country = factor(fit_preds$country)

fit_preds2 = fit_preds
fit_preds2$cases[fit_preds2$cases==0] = NA
fit_preds2$cases[fit_preds2$cases<=1] = NA
fit_preds2$country = factor(fit_preds2$country)
fit_preds2$variant = factor(fit_preds2$variant, levels=levels_VARIANTS_plot)

ggplot(data=fit_preds2, 
       aes(x=date, y=cases)) + 
  facet_wrap(~ country, scale="free", ncol=1) +
  geom_line(aes(lwd=I(1), colour=variant, group=variant)) +
  geom_line(data=fit_preds2, aes(x=date, y=totnewcases_smoothed, lwd=I(1.5)), 
            colour=alpha("black",0.3)) +
  # geom_line(aes(lwd=I(1), colour=LINEAGE2, group=LINEAGE2)) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("New confirmed cases per day") +
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES PER DAY BY VARIANT IN\nTOP 3 COUNTRIES WITH MOST BA.2.75 CASES") +
  scale_colour_manual("lineage", values=lineage_cols_plot) +
  scale_y_log10() +
  coord_cartesian(ylim=c(1,NA)) # +
# coord_cartesian(xlim=c(as.Date("2021-01-01"),max(fit_india_multi_predsbystate2$collection_date)-20))

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases multinomial fit by country.png"), width=8, height=10)


fit_preds3 = fit_preds2
fit_preds3$cases[fit_preds3$date<=as.Date("2022-03-01")] = 0

ggplot(data=fit_preds3, 
       aes(x=date, y=cases, group=variant)) + 
  facet_wrap(~ country, scale="free_y", ncol=1) +
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
  ggtitle("NEW CONFIRMED SARS-CoV2 CASES BY VARIANT IN\nTOP 3 COUNTRIES WITH MOST SEQUENCED BA.2.75 CASES") +
  coord_cartesian(xlim=c(as.Date("2022-03-01"),NA))

ggsave(file=paste0(".\\plots\\",plotdir,"\\confirmed cases stacked area multinomial fit by country.png"), width=8, height=10)




# predictions with 95% confidence intervals (slower)

# fit_multi_predsbycountry = data.frame(emmeans(fit8_multi,
#                                                   ~ variant,
#                                                   by=c("DATE_NUM", 
#                                                        "continent",
#                                                        "country"),
#                                                   at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
#                                                   mode="prob", df=NA,
#                                                   rg.limit=100000))
# fit_multi_predsbycountry$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
# fit_multi_predsbycountry$LINEAGE2 = factor(fit_india_multi_predsbystate$LINEAGE2, levels=levels_LINEAGE2)
# fit_multi_predsbycountry$STATE = factor(fit_india_multi_predsbystate$division, levels=levels_STATES)

# # test with multinomial mixed model
# 
# data_continent_country_location = GISAID_sel  %>% count(date, variant, country, location)
# # add zero counts for missing dates
# library(tidyr)
# data_continent_country_location = data_continent_country_location %>% 
#   group_by(country, variant) %>%
#   complete(date = seq(min(GISAID_sel$date), 
#                       max(GISAID_sel$date),1), 
#            fill = list(n = 0))
# data_continent_country_location = as.data.frame(data_continent_country_location)
# data_continent_country_location$DATE_NUM = as.numeric(data_continent_country_location$date)
# data_continent_country_location$continent = GISAID_sel$continent[match(data_continent_country_location$country, GISAID_sel$country)]
# data_continent_country_location$continent = droplevels(data_continent_country_location$continent)
# data_continent_country_location$variant = factor(data_continent_country_location$variant, levels=levels_VARIANTS)
# 
# fit = mblogit(variant ~ ns(DATE_NUM, df=2)*continent, 
#               # random=~1|country/location,
#               weight=n, 
#               data=data_continent_country_location, 
#               from.table=TRUE # dispersion=TRUE
#               )
# emtr_pairw_mblogit = emtrends(fit, pairwise ~ variant, by="continent", 
#                       var="DATE_NUM",  mode="latent",
#                       at=list(DATE_NUM=today_num))
# delta_r_pairw_mblogit = data.frame(confint(emtr_pairw_mblogit, 
#                                    adjust="none", df=NA)$contrasts, 
#                            p.value=as.data.frame(emtr_pairw_mblogit$contrasts)$p.value)
# delta_r_pairw_mblogit
# write.csv(delta_r_pairw, ".//data//GISAID//growth advantage BA.2.75 vs all other strains by continent.csv", row.names=F)






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

ggplot(data=cases_india_bystate[cases_india_bystate$State=="Chandigarh",],
       aes(x=Date, y=newcases))+geom_line()+coord_cartesian(xlim=c(as.Date("2022-04-01"),NA),
                                                            ylim=c(0,200))

ggplot(data=cases_india_bystate[cases_india_bystate$State=="Chandigarh",],
       aes(x=Date, y=posratio))+geom_line()+coord_cartesian(xlim=c(as.Date("2022-04-01"),NA),
                                                            ylim=c(0,30))

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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_cases per day by state.png"), width=12, height=12)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_cases per day by state.pdf"), width=12, height=12)

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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_pos ratios by state.png"), width=12, height=12)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_pos ratios by state.pdf"), width=12, height=12)


# GISAID data India

GISAID_india = GISAID_sel[as.character(GISAID_sel$country)=="India",]
GISAID_india$location = droplevels(GISAID_india$location)
GISAID_india = GISAID_india[!is.na(GISAID_india$location),]
GISAID_india = GISAID_india[!toupper(as.character(GISAID_india$location))=="NAIROBI",]
nrow(GISAID_india) # 187814
table(GISAID_india$variant, GISAID_india$location)
sort(unique(toupper(as.character(GISAID_india$location))))

GISAID_india$location = factor(toupper(as.character(GISAID_india$location)),
                               levels = c("ANDAMAN AND NICOBAR ISLANDS", "ANDHRA PRADESH", "ARUNACHAL PRADESH", "ASSAM", "BANGLADESH", "BIHAR", "CHANDIGARH", "CHHATISGARH",  "CHHATTISGARH", "DADRA AND NAGAR HAVELI",                   "DADRA AND NAGAR HAVELI AND DAMAN AND DIU", "DELHI", "GOA", "GUJARAT", "GUJART",  "GUJRAT",  "HARAYANA", "HARYANA", "HIMACHAL PRADESH", "JAMMU",             "JAMMU & KASHMIR",   "JAMMU AND KASHMIR", "JAMMU AND KASHMĪR", "JARKHAND",  "JHARKHAND", "KARNATAKA", "KERALA", "LADAKH", "LAKSHADWEEP", "MADHYA PRADESH", "MAHARASHTRA", "MAHARASTHRA", "MAHARSHTRA",  "MANIPUR", "MEGHALAYA", "MIZORAM", "MUMBAI",      "NAGALAND", "NEW DELHI", "ODISHA", "PONDICHERRY", "PUDUCHERRY", "PUNJAB", "RAJASTHAN", "SIKKIM", "TAMIL NADU", "TAMILNADU",  "TELANGANA", "TRIPURA", "UTTAR PRADESH", "UTTARAKHAND", "WEST BENGAL" ),
                               labels = c("Andaman and Nicobar Islands", "Andhra Pradesh", "Arunachal Pradesh", "Assam", "",           "Bihar", "Chandigarh", "Chhattisgarh", "Chhattisgarh", "Dadra and Nagar Haveli and Daman and Diu", "Dadra and Nagar Haveli and Daman and Diu", "Delhi", "Goa", "Gujarat", "Gujarat", "Gujarat", "Haryana",  "Haryana", "Himachal Pradesh", "Jammu and Kashmir", "Jammu and Kashmir", "Jammu and Kashmir", "Jammu and Kashmir", "Jharkhand", "Jharkhand", "Karnataka", "Kerala", "Ladakh", "Lakshadweep", "Madhya Pradesh", "Maharashtra", "Maharashtra", "Maharashtra", "Manipur", "Meghalaya", "Mizoram", "Maharashtra", "Nagaland", "Delhi",     "Odisha", "Puducherry",  "Puducherry", "Punjab", "Rajasthan", "Sikkim", "Tamil Nadu", "Tamil Nadu", "Telangana", "Tripura", "Uttar Pradesh", "Uttarakhand", "West Bengal" ))
GISAID_india = GISAID_india[as.character(GISAID_india$location)!="",]
GISAID_india$location = droplevels(GISAID_india$location)

# select only states with at least 5 BA.2.75 sequences
tab = as.data.frame(table(GISAID_india$variant, GISAID_india$location))
tab_BA_2_75 = tab[tab$Var1=="Omicron (BA.2.75)",]
tab_BA_2_75 = tab_BA_2_75[order(tab_BA_2_75$Freq, decreasing=T),]
tab_BA_2_75[tab_BA_2_75$Freq!=0,]
sel_states = as.character(tab_BA_2_75$Var2[tab_BA_2_75$Freq>=5])
sel_states
# [1] "West Bengal"      "Maharashtra"      "Rajasthan"        "Delhi"            "Gujarat"          "Odisha"          
# [7] "Himachal Pradesh" "Haryana"          "Tamil Nadu"       "Telangana"        "Chhattisgarh"     "Madhya Pradesh"  
# [13] "Chandigarh"

GISAID_india = GISAID_india[as.character(GISAID_india$location) %in%
                              sel_states, ]
GISAID_india$location = droplevels(GISAID_india$location)

# GISAID_india = GISAID_india[GISAID_india$location=="Maharashtra",]
table(GISAID_india$variant)
table(GISAID_india[GISAID_india$date>=as.Date("2022-06-01"),"variant"],
      GISAID_india[GISAID_india$date>=as.Date("2022-06-01"),"pango_lineage"])

nBA_2_75_india = sum(GISAID_india$variant=="Omicron (BA.2.75)")
nBA_2_75_india # 855

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
       subtitle=paste0("Raw GISAID data up to ",today, ", subsetted to states with at least 5 BA.2.75 sequences\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)")) 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_india_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_raw data.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_raw data.pdf"), width=8, height=6)

data_agbyweekregion1$variant2 = factor(data_agbyweekregion1$variant, levels=levels_VARIANTS_plot)
muller_indiabystate_raw2 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=variant2)) +
  facet_wrap(~ division, ncol=3) +
  geom_col(aes(width=I(10), colour=NULL, fill=variant2, group=variant2), position="fill") +
  scale_fill_manual("", values=lineage_cols_plot) +
  xaxis +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  ylab("Share") +
  theme(legend.position="right",
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
       subtitle=paste0("Raw GISAID data up to ",today, ", subsetted to states with at least 5 BA.2.75 sequences\n(n=", nBA_2_75_india," BA.2.75 sequences in India in total)")) 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india$Week)))
muller_indiabystate_raw2

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.png"), width=10, height=8)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_raw data.pdf"), width=10, height=8)



# multinomial fits
library(nnet)
library(splines)
set.seed(1)
GISAID_india$variant = factor(GISAID_india$variant, levels=levels_VARIANTS)
# fit1_india_multi = nnet::multinom(variant ~ scale(DATE_NUM)+location, data=GISAID_india, maxit=1000)
# fit2_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+location, data=GISAID_india, maxit=1000)
# fit3_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)+location, data=GISAID_india, maxit=1000)
fit1_india_multi = nnet::multinom(variant ~ scale(DATE_NUM)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
fit2_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
fit3_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)+location, weights=count, data=data_agbyweekregion1, maxit=10000)
fit4_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=2)*location, weights=count, data=data_agbyweekregion1, maxit=10000)
fit5_india_multi = nnet::multinom(variant ~ ns(DATE_NUM, df=3)*location, weights=count, data=data_agbyweekregion1, maxit=10000)
BIC(fit1_india_multi, fit2_india_multi, fit3_india_multi, fit4_india_multi, fit5_india_multi) 
# fit3_india_multi has best BIC

bestfit_india_multi = fit3_india_multi
  
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
emtrindia_pairw = emtrends(fit3_india_multi, revpairwise ~ variant,  
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
extrapolate = 30
date.from = as.numeric(as.Date("2020-06-01"))
date.to = today_num+extrapolate

# fit_india_multi_predsbystate = data.frame(emmeans(fit3_india_multi,
#                                                   ~ variant,
#                                                   by=c("DATE_NUM", 
#                                                        "location"),
#                                                   at=list(DATE_NUM=seq(date.from, date.to, by=7)), # by=7 just to speed up things a bit
#                                                   mode="prob", df=NA,
#                                                   rg.limit=100000))
# fit_india_multi_predsbystate$collection_date = as.Date(fit_india_multi_predsbystate$DATE_NUM, origin="1970-01-01")
# fit_india_multi_predsbystate$LINEAGE2 = factor(fit_india_multi_predsbystate$LINEAGE2, levels=levels_LINEAGE2)
# fit_india_multi_predsbystate$STATE = factor(fit_india_multi_predsbystate$division, levels=levels_STATES)
# 
# fit_india_multi_preds = data.frame(emmeans(fit3_india_multi, 
#                                            ~ variant,
#                                            by=c("DATE_NUM"),
#                                            at=list(DATE_NUM=seq(date.from, date.to, by=1)), # by=7 just to speed up things a bit
#                                            mode="prob", df=NA,
#                                            rg.limit=100000))
# fit_india_multi_preds$collection_date = as.Date(fit_india_multi_preds$DATE_NUM, origin="1970-01-01")
# fit_india_multi_preds$date = as.Date(fit_india_multi_preds$DATE_NUM, origin="1970-01-01")
# fit_india_multi_preds$variant = factor(fit_india_multi_preds$variant, levels=levels_VARIANTS_plot) 
# 
# fit_india_multi_preds[fit_india_multi_preds$date==today,]
# # BA.2.75 today prop 9.196145e-01
# fit_india_multi_preds[fit_india_multi_preds$date==(today-7),]
# # 70% B1.2.75 1 week ago
# 
# muller_india_mfit = ggplot(data=fit_india_multi_preds, 
#                            aes(x=collection_date, y=prob, group=variant)) + 
#   # facet_wrap(~ STATE, ncol=1) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant, group=variant), 
#             position="stack") +
#   scale_fill_manual("", values=lineage_cols_plot) +
#   annotate("rect", xmin=max(GISAID_india$DATE_NUM,na.rm=T)+1, 
#            xmax=max(GISAID_india$DATE_NUM,na.rm=T)+1+extrapolate, 
#            ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
#   xaxis +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("Share") +
#   ggtitle("SPREAD OF SARS-CoV2 VOCs IN MAHARASHTRA, INDIA\n(multinomial 2 df spline fit)")
# muller_india_mfit
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_multinom fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots_multinom fit.pdf"), width=8, height=6)


# library(ggpubr)
# ggarrange(muller_india_raw1+coord_cartesian(xlim=c(as.Date("2020-06-01"),max(GISAID_india$date, na.rm=T)+extrapolate+1))+
#             theme(legend.background = element_rect(fill = alpha("white", 0)),
#                   legend.key = element_rect(fill = alpha("white", 0)),
#                   legend.text=element_text(color = "white")) +
#             guides(colour = guide_legend(override.aes = list(alpha = 0)),
#                    fill = guide_legend(override.aes = list(alpha = 0))) +
#             ggtitle("SPREAD OF SARS-CoV2 VOCs in MAHARASHTRA, INDIA", "raw GISAID data"), 
#           muller_india_mfit+ggtitle("", "multinomial 2 df spline fit"), ncol=1)
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots multipanel_multinom fit.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots multipanel_multinom fit.pdf"), width=8, height=6)


# muller_indiabystate_mfit = ggplot(data=fit_india_multi_predsbystate, 
#                                   aes(x=collection_date, y=prob, group=LINEAGE2)) + 
#   facet_wrap(~ STATE, ncol=2) +
#   geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
#   scale_fill_manual("", values=lineage_cols2) +
#   annotate("rect", xmin=max(GISAID_india$DATE_NUM)+1, 
#            xmax=as.Date("2021-06-30"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
#   scale_x_continuous(breaks=as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                      labels=substring(months(as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                      limits=as.Date(c("2020-06-01",NA)), expand=c(0,0)) +
#   # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
#   theme_hc() + theme(legend.position="right", 
#                      axis.title.x=element_blank()) + 
#   ylab("Share") +
#   ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\n(multinomial fit)")
# muller_indiabystate_mfit
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.png"), width=6, height=8)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state_multinom fit.pdf"), width=6, height=8)
# 
# ggarrange(muller_indiabystate_raw2+
#             coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-30")))+
#             theme(legend.background = element_rect(fill = alpha("white", 0)),
#                   legend.key = element_rect(fill = alpha("white", 0)),
#                   legend.text=element_text(color = "white")) +
#             guides(colour = guide_legend(override.aes = list(alpha = 0)),
#                    fill = guide_legend(override.aes = list(alpha = 0)))+
#             ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN INDIA\nRaw GISAID data"), 
#           muller_indiabystate_mfit+ggtitle("\nMultinomial fit")+
#             coord_cartesian(xlim=c(as.Date("2020-06-01"),as.Date("2021-06-30"))), nrow=2)
# 
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state multipanel_multinom fit.png"), width=7, height=11)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_muller plots by state multipanel_multinom fit.pdf"), width=7, height=11)

# PLOT MODEL FIT WITH DATA

# plot predicted values
# multinomial model predictions by province (fastest, but no confidence intervals)
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

# on logit scale:

fit_india_multi_preds2 = fit_preds
ymin = 0.0001
ymax = 0.999
if (!is.null(fit_india_multi_preds2$asymp.LCL)) fit_india_multi_preds2$asymp.LCL[fit_india_multi_preds2$asymp.LCL<ymin] = ymin
if (!is.null(fit_india_multi_preds2$asymp.UCL)) fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL<ymin] = ymin
if (!is.null(fit_india_multi_preds2$asymp.UCL)) fit_india_multi_preds2$asymp.UCL[fit_india_multi_preds2$asymp.UCL>ymax] = ymax
fit_india_multi_preds2$prob[fit_india_multi_preds2$prob<ymin] = ymin

plot_india_mfit_logit = qplot(data=fit_india_multi_preds2[fit_india_multi_preds2$variant!="Other",], 
                              x=date, y=prob, geom="blank") +
  facet_wrap(~location, ncol=3) +
  # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
  #                 fill=variant
  # ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN INDIA",
          subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=3)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_multinom fit_logit scale.png"), width=10, height=6)


# on response scale:
plot_india_mfit = qplot(data=fit_india_multi_preds2[fit_india_multi_preds2$variant!="Other",], x=date, 
                        y=100*prob, geom="blank") +
  facet_wrap(~location, ncol=3) +
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
          subtitle=paste0("GISAID data up to ",today, ", multinomial spline fit variant ~ ns(date, df=3)+state,\nall states with >=5 BA.2.75 sequences shown (n=", nBA_2_75_india," BA.2.75 sequences in total)")) +
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

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_multinom fit_response scale.png"), width=10, height=6)


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
  facet_wrap(~ location, scale="free", ncol=3) +
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
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=3) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases multinomial fit.png"), width=11, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases multinomial fit.pdf"), width=8, height=10)

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=cases, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=3) +
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
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=3) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))


ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases stacked area multinomial fit.png"), width=11, height=6)
# ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases stacked area multinomial fit.pdf"), width=8, height=10)

fit_india_multi_predsbystate3 = fit_india_multi_predsbystate2
fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date<as.Date("2022-05-01"),"cases"] = 0
fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date<as.Date("2022-05-01"),"posratio"] = 0

ggplot(data=fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date>=as.Date("2022-04-30"),], 
       aes(x=collection_date, y=cases, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=3) +
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
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=3) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_confirmed cases stacked area multinomial fit_ZOOMED.png"), width=10, height=6)

# graph with positivity ratios by variant
ggplot(data=fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date>=as.Date("2022-04-30"),], 
       aes(x=collection_date, y=posratio, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=3) +
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
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=3) + state) to GISAID data\nfile version ",download,
          " plus records from last days\ncase & testing data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_positivity ratios stacked area multinomial fit_ZOOMED.png"), width=10, height=6)
#write.csv(fit_india_multi_predsbystate3[fit_india_multi_predsbystate3$collection_date>=as.Date("2022-04-30"),],
#                                        "multinom fit India.csv", row.names=F)

ggplot(data=fit_india_multi_predsbystate2, 
       aes(x=collection_date, y=posratio, group=variant)) + 
  facet_wrap(~ location, scale="free", ncol=3) +
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
          subtitle=paste0("variant frequencies based on multinomial spline fit variant ~ ns(date, df=3) + state) to GISAID data\nfile version ",download,
                          " plus records from last days\ncase & testing data from data.covid19bharat.org/csv/latest/states.csv")) +
  labs(tag = tag) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8))

ggsave(file=paste0(".\\plots\\",plotdir,"\\india_positivity ratios stacked area multinomial fit.png"), width=11, height=6)
