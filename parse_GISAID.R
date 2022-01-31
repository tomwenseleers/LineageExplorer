library(readr)

# import GISAID metadata (file version metadata_tsv_2022_01_28.tar.xz)
gc()
message("Parsing GISAID metadata...")
GISAID = read_tsv(".//data//GISAID//metadata.tsv", col_types = cols(.default = "c"), quote = "\"")
GISAID = as.data.frame(GISAID)
# dim(GISAID)

# only keep human samples
GISAID = GISAID[GISAID$Host=="Human",] 
GISAID = GISAID[!is.na(GISAID$Type),] 
# nrow(GISAID) 
# names(GISAID)
# head(GISAID$Location)
# head(GISAID$"Collection date")
# head(GISAID$"Pango lineage")
# unique(GISAID$"Virus name")
# unique(GISAID$"Type")
# unique(GISAID$"Pango lineage")
# unique(GISAID$"Pangolin version") # "2022-01-05"
# unique(GISAID$"Variant")
# unique(GISAID$"AA Substitutions")

# parse dates & remove records with invalid dates
library(stringr)
date_isvalid = (str_count(GISAID$"Collection date", pattern = "-")==2)
GISAID = GISAID[date_isvalid,]
GISAID$date = as.Date(GISAID$"Collection date") 
GISAID = GISAID[!is.na(GISAID$date),]
# nrow(GISAID)

# only keep complete genomes
GISAID = GISAID[!is.na(GISAID$"Is complete?"),]
# nrow(GISAID) 

# add numeric version of date, week midpoint & week & year
# range(GISAID$date) # "2019-12-24" "2022-01-16"
# GISAID[GISAID$date==min(GISAID$date),] # first available Wuhan genome from Dec 24 (submitted 11 Jan 2020)
GISAID$Week = lubridate::week(GISAID$date)
GISAID$Year = lubridate::year(GISAID$date)
GISAID$Year_Week = interaction(GISAID$Year,GISAID$Week)
library(lubridate)
GISAID$week_midpoint = as.Date(as.character(cut(GISAID$date, "week")))+3.5 # week midpoint date
GISAID$DATE_NUM = as.numeric(GISAID$date)
# names(GISAID)

# parse continent / country / location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID$Location), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
# unique(loc[,1]) # continents : "South America" "North America" "Europe"        "Oceania"       "Asia"          "Africa" 

library(countrycode)
# canonicalize country names
loc[,2] = countrycode( countrycode(loc[,2], "country.name", "iso3c",
                                            custom_match = c('Bonaire' = 'BES',
                                                             'Canary Islands' = 'ESP',
                                                             'Crimea' = 'UKR',
                                                             'Kosovo' = 'SRB',
                                                             'Saint Martin' = 'SXM', 
                                                             'Sint Eustatius' = 'BES') ), 
                                   "iso3c", "country.name" ) 

# sort(unique(loc[,2])) # country

loc2 = loc # recoded version of country with Scotland, England, Wales & Northern Ireland recoded separately for UK
loc2[,2][loc2[,3]=="Scotland"] = "Scotland"
loc2[,2][loc2[,3]=="England"] = "England"
loc2[,2][loc2[,3]=="Wales"] = "Wales"
loc2[,2][loc2[,3]=="Northern Ireland"] = "Northern Ireland"
loc2[,2][loc2[,2]=="United Kingdom"] = "England"

# sort(unique(loc2[,2])) # country, with United Kingdom split up in England, Scotland, Wales & Northern Ireland

# sort(unique(loc[,3])) # city or province
# unique(loc[,4][loc[,2]=="England"])
sort(unique(loc[,3][loc[,2]=="United States"])) # US states - needs some tidying up
sort(unique(loc[,3][loc[,2]=="South Africa"])) # SA states - needs some tidying up
sort(unique(loc[,3][loc[,2]=="India"])) # provinces of India - needs some tidying up
sort(unique(loc[,3][loc[,2]=="Canada"])) # for Canada - needs some tidying up
sort(unique(loc[,3][loc[,2]=="Russia"])) # for Russia - needs some tidying up # vast majority from Saint-Petersburg              Saint Petersburg & Moscow                 Moscow region                 Moscow Region
sort(unique(loc[,3][loc[,2]=="Brazil"])) # for Brazil - needs some tidying up
sort(unique(loc[,3][loc[,2]=="Australia"])) # for Australia - needs some tidying up

# unique(GISAID$'Additional location information'[loc[,2]=="England"]) # sometimes has ZIP code, can maybe be converted to ONS region?
# unique(loc[,3][loc[,2]=="Belgium"]) # usually city
# unique(loc[,3][loc[,2]=="Netherlands"]) # province, incl Amsterdam
# unique(loc[,4][loc[,2]=="Netherlands"]) # city, incl Amsterdam
# unique(GISAID$'Additional location information'[loc[,2]=="Netherlands"]) # ZIP code
# sum(grepl("1012|1011|1015|1016|1017|1052|1018|1013|1051|1031|1021|1053|1072|1074|1073|1071|1054|1091|1093|1056", GISAID$'Additional location information'[loc[,2]=="Netherlands"]), na.rm=T)
# sum((loc[,4][loc[,2]=="Netherlands"]=="Amsterdam")|(loc[,3][loc[,2]=="Netherlands"]=="Amsterdam"),na.rm=T) # 167
# unique(loc[,3][loc[,2]=="Switzerland"]) # city, incl Geneva / Genève
# sum((loc[,4][loc[,2]=="Switzerland"]=="Geneva")|(loc[,3][loc[,2]=="Switzerland"]=="Geneva"),na.rm=T) # 10459
# sum((loc[,4][loc[,2]=="Switzerland"]=="Genève")|(loc[,3][loc[,2]=="Switzerland"]=="Genève"),na.rm=T) # 8

levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
GISAID$continent = factor(loc[,1], levels=levels_continents)
GISAID$country = factor(loc[,2])
levels_countries = levels(GISAID$country)
GISAID$iso3c = countrycode(loc[,2], "country.name", "iso3c")
GISAID$country2 = factor(loc2[,2])
levels_countries2 = levels(GISAID$country2)
GISAID$location = factor(loc[,3])
levels_locations = levels(GISAID$location)

# code VOCs based on pango lineages & Variant field
GISAID = GISAID[!(GISAID$"Pango lineage"=="None"|GISAID$"Pango lineage"==""|is.na(GISAID$"Pango lineage")),]
#nrow(GISAID) 
unique(GISAID$"Pango lineage")
length(unique(GISAID$"Pango lineage")) # 1576 unique pango lineages

# sel=GISAID[GISAID$`Pango lineage`=="B.1.1.529",]
# sel[order(sel$date),]

library(dplyr)
GISAID$VARIANT = case_when(
  grepl("Alpha", GISAID$Variant) ~ "Alpha", 
  grepl("Beta", GISAID$Variant) ~ "Beta", 
  grepl("Gamma", GISAID$Variant) ~ "Gamma",
  grepl("Delta", GISAID$Variant) ~ "Delta",
  grepl("Omicron", GISAID$Variant) ~ "Omicron",
  T ~ "Other"
)
# table(GISAID$VARIANT)
# pie(table(GISAID$VARIANT))
VARIANT_levels = c("Other", "Alpha", "Beta", "Gamma", "Delta", "Omicron")
GISAID$VARIANT = factor(GISAID$VARIANT, levels=VARIANT_levels)
VARIANT_colours = c("grey70","#0085FF","#9A9D00","cyan3","magenta","blue")

# with Omicron subtypes coded separately
GISAID$VARIANT2 = case_when(
  grepl("Alpha", GISAID$Variant) ~ "Alpha", 
  grepl("Beta", GISAID$Variant) ~ "Beta", 
  grepl("Gamma", GISAID$Variant) ~ "Gamma",
  grepl("Delta", GISAID$Variant) ~ "Delta",
  (GISAID$"Pango lineage"=="BA.1") ~ "Omicron (BA.1)",
  (GISAID$"Pango lineage"=="BA.1.1") ~ "Omicron (BA.1.1)",
  grepl("BA.2", GISAID$"Pango lineage") ~ "Omicron (BA.2)",
  grepl("BA.3", GISAID$"Pango lineage") ~ "Omicron (BA.3)",
  T ~ "Other"
)
# table(GISAID$VARIANT2)
# pie(table(GISAID$VARIANT2))
VARIANT2_levels = c("Other", "Alpha", "Beta", "Gamma", "Delta", "Omicron (BA.1)", "Omicron (BA.1.1)", "Omicron (BA.2)", "Omicron (BA.3)")
GISAID$VARIANT2 = factor(GISAID$VARIANT2, levels=VARIANT2_levels)
VARIANT2_colours = c("grey70","#0085FF","#9A9D00","cyan3","magenta",
                     colorRampPalette(c("red", "orange", "blue"))(4)
                     )
# code sgtf (S gene target failure) based on presence of Spike_H69del / Spike_V70del
GISAID$sgtf = case_when(
  grepl("Spike_H69del",GISAID$"AA Substitutions", fixed=T)|grepl("Spike_V70del",GISAID$"AA Substitutions", fixed=T) ~ "sgtf", 
  T ~ "non_sgtf"
)
table(GISAID$VARIANT2,GISAID$sgtf) # Alpha, B.1.1.529, BA.1, BA.1.1 and BA.3 generally have S dropout, while Beta, Delta, Gamma, BA.2 & Other strains generally do not

saveRDS(GISAID, ".//data//GISAID//GISAID.rds")


# TO DO ALSO ADD LAT & LONGITUDE OF EACH COUNTRY?


# convert to frequency table by country, date (week midpoint) & variant ####
# for VARIANT
GISAID_ftable = as.data.frame(ftable(formula = VARIANT ~ week_midpoint+country, data=GISAID))
head(GISAID_ftable)
# add total count over all variants per date
GISAID_ftable_total = as.data.frame(GISAID_ftable %>% 
                                      group_by(week_midpoint, country) %>% 
                                      summarise(across(c(Freq), list(Total = sum))))
GISAID_ftable_total$week_midpoint = as.Date(as.character(GISAID_ftable_total$week_midpoint))
# library(ggplot2)
# library(ggthemes)
# qplot(data=GISAID_ftable_total, x=week_midpoint, y=Freq_Total, geom="line")+facet_wrap(~country, scale="free_y")
GISAID_ftable$Total = GISAID_ftable_total$Freq_Total[match(interaction(GISAID_ftable$week_midpoint, GISAID_ftable$country),
                                                           interaction(GISAID_ftable_total$week_midpoint, GISAID_ftable_total$country))]
GISAID_ftable = GISAID_ftable[GISAID_ftable$Total!=0,]
GISAID_ftable$Prop = GISAID_ftable$Freq/GISAID_ftable$Total
# add latitude & longitude of country centroids
# source(".//get_longitudes_latitudes_country_centroids.R") # returns dataframe countries_long_lat, manually curated version below
countries_long_lat = read.csv(".//data//GISAID//countries_long_lat_centroids_curated.csv")
GISAID_ftable$longitude = countries_long_lat$longitude[match(GISAID_ftable$country, countries_long_lat$country)]
GISAID_ftable$latitude = countries_long_lat$latitude[match(GISAID_ftable$country, countries_long_lat$country)]
# add continent
GISAID_ftable$continent = GISAID$continent[match(GISAID_ftable$country, GISAID$country)]
# add numeric version of date
GISAID_ftable$week_midpoint = as.Date(as.character(GISAID_ftable$week_midpoint))
GISAID_ftable$date_num = as.numeric(GISAID_ftable$week_midpoint)
GISAID_ftable$colour = VARIANT_colours[match(GISAID_ftable$VARIANT, VARIANT_levels)]
write.csv(GISAID_ftable, ".//data//GISAID//GISAID_freq_VARIANT_byweek_country.csv", row.names=F)

# for VARIANT2
GISAID_ftable2 = as.data.frame(ftable(formula = VARIANT2 ~ week_midpoint+country, data=GISAID))
head(GISAID_ftable2)
# add total count over all variants per date
GISAID_ftable2_total = as.data.frame(GISAID_ftable2 %>% 
                                      group_by(week_midpoint, country) %>% 
                                      summarise(across(c(Freq), list(Total = sum))))
GISAID_ftable2_total$week_midpoint = as.Date(as.character(GISAID_ftable2_total$week_midpoint))
GISAID_ftable2$Total = GISAID_ftable2_total$Freq_Total[match(interaction(GISAID_ftable2$week_midpoint, GISAID_ftable2$country),
                                                           interaction(GISAID_ftable2_total$week_midpoint, GISAID_ftable2_total$country))]
GISAID_ftable2 = GISAID_ftable2[GISAID_ftable2$Total!=0,]
GISAID_ftable2$Prop = GISAID_ftable2$Freq/GISAID_ftable2$Total
GISAID_ftable2$longitude = countries_long_lat$longitude[match(GISAID_ftable2$country, countries_long_lat$country)]
GISAID_ftable2$latitude = countries_long_lat$latitude[match(GISAID_ftable2$country, countries_long_lat$country)]
# add continent
GISAID_ftable2$continent = GISAID$continent[match(GISAID_ftable2$country, GISAID$country)]
# add numeric version of date
GISAID_ftable2$week_midpoint = as.Date(as.character(GISAID_ftable2$week_midpoint))
GISAID_ftable2$date_num = as.numeric(GISAID_ftable2$week_midpoint)
GISAID_ftable2$colour = VARIANT2_colours[match(GISAID_ftable2$VARIANT2, VARIANT2_levels)]
write.csv(GISAID_ftable2, ".//data//GISAID//GISAID_freq_VARIANT2_byweek_country.csv", row.names=F)



# broken down by sgtf status & variant
# using country2 (with UK split up in England, Scotland, Wales & Northern Ireland)
GISAID_ftable_sgtf = as.data.frame(ftable(formula = sgtf ~ week_midpoint+country2, data=GISAID))
library(tidyr)
GISAID_ftable_sgtf = spread(GISAID_ftable_sgtf, sgtf, Freq)
GISAID_ftable_VARIANT2 = as.data.frame(ftable(formula = VARIANT2 ~ week_midpoint+country2, data=GISAID))
GISAID_ftable_VARIANT2 = spread(GISAID_ftable_VARIANT2, VARIANT2, Freq)
GISAID_ftable_sgtf_VARIANTS = cbind(GISAID_ftable_sgtf, GISAID_ftable_VARIANT2[,-c(1:2)])
GISAID_ftable_sgtf_VARIANTS$Omicron = GISAID_ftable_sgtf_VARIANTS$`Omicron (B.1.1.529)` +
                                      GISAID_ftable_sgtf_VARIANTS$`Omicron (BA.1)` +
                                      GISAID_ftable_sgtf_VARIANTS$`Omicron (BA.2)` +
                                      GISAID_ftable_sgtf_VARIANTS$`Omicron (BA.3)`
GISAID_ftable_sgtf_VARIANTS$Total = GISAID_ftable_sgtf_VARIANTS$sgtf + GISAID_ftable_sgtf_VARIANTS$non_sgtf
GISAID_ftable_sgtf_VARIANTS$longitude = countries_long_lat$longitude[match(GISAID_ftable_sgtf_VARIANTS$country2, countries_long_lat$country)] 
GISAID_ftable_sgtf_VARIANTS$latitude = countries_long_lat$latitude[match(GISAID_ftable_sgtf_VARIANTS$country2, countries_long_lat$country)]
# add continent
GISAID_ftable_sgtf_VARIANTS$continent = GISAID$continent[match(GISAID_ftable_sgtf_VARIANTS$country2, GISAID$country2)]
# add numeric version of date
GISAID_ftable_sgtf_VARIANTS$week_midpoint = as.Date(as.character(GISAID_ftable_sgtf_VARIANTS$week_midpoint))
GISAID_ftable_sgtf_VARIANTS$date_num = as.numeric(GISAID_ftable_sgtf_VARIANTS$week_midpoint)
head(GISAID_ftable_sgtf_VARIANTS)
write.csv(GISAID_ftable_sgtf_VARIANTS, ".//data//GISAID//GISAID_freq_sgtf_VARIANTS_byweek_country2.csv", row.names=F)

# using country (with UK not split up into England, Scotland, Wales & Northern Ireland)
GISAID_ftable_sgtf = as.data.frame(ftable(formula = sgtf ~ week_midpoint+country, data=GISAID))
library(tidyr)
GISAID_ftable_sgtf = spread(GISAID_ftable_sgtf, sgtf, Freq)
GISAID_ftable_VARIANT2 = as.data.frame(ftable(formula = VARIANT2 ~ week_midpoint+country, data=GISAID))
GISAID_ftable_VARIANT2 = spread(GISAID_ftable_VARIANT2, VARIANT2, Freq)
GISAID_ftable_sgtf_VARIANTS = cbind(GISAID_ftable_sgtf, GISAID_ftable_VARIANT2[,-c(1:2)])
GISAID_ftable_sgtf_VARIANTS$Omicron = GISAID_ftable_sgtf_VARIANTS$`Omicron (B.1.1.529)` +
  GISAID_ftable_sgtf_VARIANTS$`Omicron (BA.1)` +
  GISAID_ftable_sgtf_VARIANTS$`Omicron (BA.2)` +
  GISAID_ftable_sgtf_VARIANTS$`Omicron (BA.3)`
GISAID_ftable_sgtf_VARIANTS$Total = GISAID_ftable_sgtf_VARIANTS$sgtf + GISAID_ftable_sgtf_VARIANTS$non_sgtf
GISAID_ftable_sgtf_VARIANTS$longitude = countries_long_lat$longitude[match(GISAID_ftable_sgtf_VARIANTS$country, countries_long_lat$country)]
GISAID_ftable_sgtf_VARIANTS$latitude = countries_long_lat$latitude[match(GISAID_ftable_sgtf_VARIANTS$country, countries_long_lat$country)]
# add continent
GISAID_ftable_sgtf_VARIANTS$continent = GISAID$continent[match(GISAID_ftable_sgtf_VARIANTS$country, GISAID$country)]
# add numeric version of date
GISAID_ftable_sgtf_VARIANTS$week_midpoint = as.Date(as.character(GISAID_ftable_sgtf_VARIANTS$week_midpoint))
GISAID_ftable_sgtf_VARIANTS$date_num = as.numeric(GISAID_ftable_sgtf_VARIANTS$week_midpoint)
head(GISAID_ftable_sgtf_VARIANTS)
write.csv(GISAID_ftable_sgtf_VARIANTS, ".//data//GISAID//GISAID_freq_sgtf_VARIANTS_byweek_country.csv", row.names=F)


# old code to delete

# library(tidyr)
# # proportion of sgtf samples that are omicron in function of time by country & date to run logistic regression on
# GISAID_sgtf_isomicron = spread(GISAID_sgtf_isomicron, Omicron, Freq) 
# GISAID_sgtf_isomicron$date = as.Date(GISAID_sgtf_isomicron$date)
# GISAID_sgtf_isomicron$DATE_NUM = as.numeric(GISAID_sgtf_isomicron$date)
# saveRDS(GISAID_sgtf_isomicron, file=".//data//GISAID//GISAID_sgtf_isomicron.rds")
# 
# # code Omicron lineage
# GISAID$Omicron = "Other"
# GISAID$Omicron[GISAID$VARIANT=="Omicron"] = "Omicron"
# 
# # code Delta lineage
# GISAID$Delta = "Other"
# GISAID$Delta[grepl("Delta",GISAID$Variant)] = "Delta"
# 
# GISAID_sgtf_isomicron = as.data.frame(ftable(formula=Omicron~country+date, data=GISAID[GISAID$sgtf=="sgtf",]))
# library(tidyr)
# # proportion of sgtf samples that are omicron in function of time by country & date to run logistic regression on
# GISAID_sgtf_isomicron = spread(GISAID_sgtf_isomicron, Omicron, Freq) 
# GISAID_sgtf_isomicron$date = as.Date(GISAID_sgtf_isomicron$date)
# GISAID_sgtf_isomicron$DATE_NUM = as.numeric(GISAID_sgtf_isomicron$date)
# saveRDS(GISAID_sgtf_isomicron, file=".//data//GISAID//GISAID_sgtf_isomicron.rds")







