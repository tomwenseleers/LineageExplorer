# PARSE GISAID METADATA
# T. Wenseleers, last update 9 February 2022

library(readr)
library(stringr)
library(lubridate)
library(stringi)
library(countrycode)
# devtools::install_github("SymbolixAU/googleway")
library(ggmap)
library(googleway)
# PS first register your Google Maps API key using
# register_google(key="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX", account_type="premium", write=TRUE)
library(futile.logger)
library(utils)
library(plyr)
library(curl)
library(sp)
library(rgdal)
library(rgeos)
library(rworldmap)
library(rworldxtra)


# import GISAID metadata (file version metadata_tsv_2022_02_04.tar.xz) (or go directly to line 65)
# gc()
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

date_isvalid = (str_count(GISAID$"Collection date", pattern = "-")==2)
GISAID = GISAID[date_isvalid,]

GISAID$date = as.Date(fast_strptime(GISAID$"Collection date", "%Y-%m-%d")) 
GISAID = GISAID[!is.na(GISAID$date),]
GISAID$date_num = as.numeric(GISAID$date) # numeric version of date
# nrow(GISAID)

# only keep complete genomes
GISAID = GISAID[!is.na(GISAID$"Is complete?"),]
# nrow(GISAID) 

# add numeric version of date, week midpoint & week & year
# range(GISAID$date) # "2019-12-24" "2022-02-02"
# GISAID[GISAID$date==min(GISAID$date),] # first available Wuhan genome from Dec 24 (submitted 11 Jan 2020)
GISAID$week = lubridate::week(GISAID$date)
GISAID$year = lubridate::year(GISAID$date)
GISAID$week_startdate = floor_date(GISAID$date, "weeks", week_start = 1) 
# names(GISAID)


# read GISAID, here from rds

GISAID = readRDS(".//data//GISAID//GISAID.rds")

# parse continent / country / location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(GISAID$Location), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
# unique(loc[,1]) # continents : "South America" "North America" "Europe"        "Oceania"       "Asia"          "Africa" 

loc = apply(loc, 2, stri_enc_toutf8) # convert to UTF-8
loc[,3][loc[,3] %in% c("Unknown", "NA", "unknown", "?")] = NA



# canonicalize country names
# countries = data.frame(country_GISAID = sort(unique(loc[,2]))) # original country names in GISAID
# countries$iso3c = countrycode(countries$country_GISAID, "country.name", "iso3c",
#                                            custom_match = c('Bonaire' = 'BES',
#                                                             'Canary Islands' = 'ESP',
#                                                             'Crimea' = 'UKR',
#                                                             'Kosovo' = 'KSV',
#                                                             'Saint Martin' = 'SXM', 
#                                                             'Sint Eustatius' = 'BES'))
# countries$country = countrycode(countries$iso3c, "iso3c", "country.name",
#                                 custom_match = c('SOL' = 'Somaliland',
#                                                  'KSV' = 'Kosovo'))  # canonicalized country name
# countries$continent = loc[,1][match(countries$country_GISAID, loc[,2])]
# countries$un.regionsub.name = countrycode(countries$iso3c, "iso3c", "un.regionsub.name", 
#                                   custom_match = c('KSV' = 'Southern Europe',
#                                                   'TWN' = 'Eastern Asia'))
# length(sort(unique(countries$country))) # 209

# sort(unique(loc[,2])) # country

# file with latitudes & longitudes of countries & canonicalized country names & continents
# (country names are as given by the countrycode package)
countries_long_lat = read.csv(".//data//GISAID//countries_long_lat_centroids_curated_final.csv") 

GISAID$country_GISAID = loc[,2]
GISAID$iso3c = countrycode(GISAID$country_GISAID, "country.name", "iso3c",
                                            custom_match = c('Bonaire' = 'BES',
                                                             'Canary Islands' = 'ESP',
                                                             'Crimea' = 'UKR',
                                                             'Kosovo' = 'SRB',
                                                             'Saint Martin' = 'SXM', 
                                                             'Sint Eustatius' = 'BES') )
# GISAID$country = countrycode(GISAID$iso3c, "iso3c", "country.name",
#                              custom_match = c('SOL' = 'Somaliland',
#                                               'KSV' = 'Kosovo'))  # canonicalized country name
GISAID$country = countries_long_lat$country[match(GISAID$iso3c, countries_long_lat$iso3c)]  # canonicalized country name
GISAID$continent = countries_long_lat$continent[match(GISAID$iso3c, countries_long_lat$iso3c)]
GISAID$un.regionsub.name = countrycode(GISAID$iso3c, "iso3c", "un.regionsub.name", 
                                       custom_match = c('KSV' = 'Southern Europe',
                                                        'TWN' = 'Eastern Asia'))
GISAID$region = loc[,3]
GISAID$continent.country.region = factor(tidyr::unite_(GISAID, col="x", from=c("continent","country","region"), sep=" / ")$x)

GISAID_ftable_countries_regions = as.data.frame(xtabs(formula = ~ continent.country.region, data=GISAID))
GISAID_ftable_countries_regions$continent = GISAID$continent[match(GISAID_ftable_countries_regions$continent.country.region, GISAID$continent.country.region)]
GISAID_ftable_countries_regions$country = GISAID$country[match(GISAID_ftable_countries_regions$continent.country.region, GISAID$continent.country.region)]
GISAID_ftable_countries_regions$region = GISAID$region[match(GISAID_ftable_countries_regions$continent.country.region, GISAID$continent.country.region)]
GISAID_ftable_countries_regions$iso2c = countrycode(GISAID_ftable_countries_regions$country, "country.name", "iso2c")
GISAID_ftable_countries_regions$iso3c = countrycode(GISAID_ftable_countries_regions$country, "country.name", "iso3c")

# write.csv(GISAID_ftable_countries_regions, ".//data//GISAID//GISAID_countries_regions.csv", row.names=F)
sum(GISAID_ftable_countries_regions$Freq>100) # 1239 regions with >100 sequences

GISAID_ftable_countries = as.data.frame(xtabs(formula = ~ country, data=GISAID))
sum(GISAID_ftable_countries$Freq>1000) # 93 countries with >1000 sequences
sort(GISAID_ftable_countries$country[GISAID_ftable_countries$Freq>1000])

apiKey = google_key()

# retry utility function

retry = function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=100, sleep=0.02) {
  attempts = 0
  retval = try(eval(expr))
  while (isError(retval)) {
    attempts = attempts + 1
    if (attempts >= maxErrors) {
      msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
      flog.fatal(msg)
      stop(msg)
    } else {
      msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts,  maxErrors, 
                    capture.output(str(retval)))
      flog.error(msg)
      warning(msg)
      
    }
    if (sleep > 0) Sys.sleep(sleep)
    retval = try(eval(expr))
  }
  return(retval)
}



# look up latitude & longitude of regions using googleway::google_geocode
message("Looking up region longitudes & latitudes")
get_geocode_google = function (i) { 
  Sys.sleep(0.02)  
  query = stri_enc_toutf8(GISAID_ftable_countries_regions$country[i])
  if (!is.na(GISAID_ftable_countries_regions$region[i])) { query = paste0(stri_enc_toutf8(GISAID_ftable_countries_regions$region[i]), ", ", query) } 
  # while(!has_internet()) {}
  res = google_geocode(address = query, 
                       region = substring(countrycode(GISAID_ftable_countries_regions$country[i], "country.name", "cctld"),2), 
                       key = apiKey)$results$geometry$location[1,]
  if (is.null(res)) res = data.frame(lat=NA, lng=NA)
  res = data.frame(row=i, res);
  res
}

# GISAID_ftable_countries_regions_latlongs =  read.csv(".//data//GISAID//GISAID_ftable_countries_regions_latlongs.csv")

GISAID_ftable_countries_regions_latlongs = do.call(rbind, llply(1:nrow(GISAID_ftable_countries_regions), 
                                                  function (i) retry(get_geocode_google(i)), .progress = "text"))

# write.csv(GISAID_ftable_countries_regions_latlongs, ".//data//GISAID//GISAID_ftable_countries_regions_latlongs.csv", row.names=F)

GISAID_ftable_countries_regions$latitude_country = countries_long_lat$latitude[match(GISAID_ftable_countries_regions$iso3c, countries_long_lat$iso3c)]
GISAID_ftable_countries_regions$longitude_country = countries_long_lat$longitude[match(GISAID_ftable_countries_regions$iso3c, countries_long_lat$iso3c)]

GISAID_ftable_countries_regions$latitude_region = GISAID_ftable_countries_regions_latlongs$lat
GISAID_ftable_countries_regions$longitude_region = GISAID_ftable_countries_regions_latlongs$lng

# use latitude & longitude of country centroid for those regions that returned an NA for latitude & longitude or that had an NA for region
# (alternative would be to remove those records)
# GISAID_ftable_countries_regions[is.na(GISAID_ftable_countries_regions$latitude_region),c("continent.country.region","Freq")]
# sum(GISAID_ftable_countries_regions[is.na(GISAID_ftable_countries_regions$latitude_region),c("Freq")]) # 219

GISAID_ftable_countries_regions$latitude_region[is.na(GISAID_ftable_countries_regions$latitude_region)|is.na(GISAID_ftable_countries_regions$region)] = GISAID_ftable_countries_regions$latitude_country[is.na(GISAID_ftable_countries_regions$latitude_region)|is.na(GISAID_ftable_countries_regions$region)]
GISAID_ftable_countries_regions$longitude_region[is.na(GISAID_ftable_countries_regions$longitude_region)|is.na(GISAID_ftable_countries_regions$region)] = GISAID_ftable_countries_regions$longitude_country[is.na(GISAID_ftable_countries_regions$longitude_region)|is.na(GISAID_ftable_countries_regions$region)]

write.csv(GISAID_ftable_countries_regions, ".//data//GISAID//GISAID_countries_regions_withlatlong.csv", row.names=F)

# save.image("~/GitHub/newcovid_belgium/environment_7_2_2022.RData")

# reverse geocode to locality, administrative_area_level_1, administrative_area_level_2 & country using googleway::google_reverse_geocode
message("Reverse geocoding to locality, administrative_area_level_1, administrative_area_level_2 & country using Google API")

get_revgeocode_google = function (i) { Sys.sleep(0.02); 
  # print(i); 
  # while(!has_internet()) {}
  res = googleway::google_reverse_geocode(
    location = c(GISAID_ftable_countries_regions$latitude_region[i], 
                 GISAID_ftable_countries_regions$longitude_region[i]),
    result_type = "political|country|administrative_area_level_1|administrative_area_level_2|locality",
    language = "en", # return result in English
    #location_type = "GEOMETRIC_CENTER|APPROXIMATE",
    key = apiKey
  )
  addr_comps = geocode_address_components(res)
  locality = addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="locality",]$long_name 
  if (length(locality)==0) locality=NA 
  administrative_area_level_1 = addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="administrative_area_level_1",]$long_name 
  if (length(administrative_area_level_1)==0) administrative_area_level_1=NA 
  administrative_area_level_2 = addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="administrative_area_level_2",]$long_name 
  if (length(administrative_area_level_2)==0) administrative_area_level_2=NA 
  country = addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="country",]$long_name  
  if (length(country)==0) country=NA 
  res = data.frame(row=i, 
                   locality_revgeo=locality, 
                   administrative_area_level_1_revgeo=administrative_area_level_1,
                   administrative_area_level_2_revgeo=administrative_area_level_2,
                   country_revgeo=country);
  res
}

# GISAID_ftable_revgeo = read.csv(".//data//GISAID//GISAID_ftable_revgeo.csv")

GISAID_ftable_revgeo = do.call(rbind, llply(1:nrow(GISAID_ftable_countries_regions), 
                                                                function (i) retry(get_revgeocode_google(i)), .progress = "text"))

# write.csv(GISAID_ftable_revgeo, ".//data//GISAID//GISAID_ftable_revgeo.csv", row.names=F)

sum(is.na(GISAID_ftable_revgeo$country_revgeo)) # 90 record where reverse geocoding failed - we use country centroid for those

GISAID_ftable_countries_regions$locality_revgeo = GISAID_ftable_revgeo$locality_revgeo
GISAID_ftable_countries_regions$administrative_area_level_1_revgeo = GISAID_ftable_revgeo$administrative_area_level_1_revgeo
GISAID_ftable_countries_regions$administrative_area_level_2_revgeo = GISAID_ftable_revgeo$administrative_area_level_2_revgeo
GISAID_ftable_countries_regions$country_revgeo = countrycode(countrycode(GISAID_ftable_revgeo$country_revgeo, "country.name", "iso3c", custom_match = c('Bonaire' = 'BES',
                                                                                                                                            'Canary Islands' = 'ESP',
                                                                                                                                            'Crimea' = 'UKR',
                                                                                                                                            'Kosovo' = 'KSV',
                                                                                                                                            'Saint Martin' = 'SXM', 
                                                                                                                                            'Sint Eustatius' = 'BES')),
                                                             "iso3c", "country.name")

sort(unique(GISAID_ftable_countries_regions[GISAID_ftable_countries_regions$country_revgeo!=GISAID_ftable_countries_regions$country,"country"]))
GISAID_ftable_countries_regions[GISAID_ftable_countries_regions$country_revgeo!=GISAID_ftable_countries_regions$country,c("country","country_revgeo","Freq")]
sum(GISAID_ftable_countries_regions[GISAID_ftable_countries_regions$country_revgeo!=GISAID_ftable_countries_regions$country,c("country","country_revgeo","Freq")]$Freq, na.rm=T)
# 891 = 0.01% of all records return the wrong country
sum(GISAID_ftable_countries_regions[is.na(GISAID_ftable_countries_regions$country_revgeo),"Freq"]) # 41215 = 0.5% of all records where reverse geocoding failed
sum(GISAID_ftable_countries_regions$Freq, na.rm=T) # 7589371 total records


# for records where reverse geocoded country did not match original country or where reverse geocoded country returned NA we 
# replace latitude_region and longitude_region by country centroid & country_revgeo by country

write.csv(GISAID_ftable_countries_regions, ".//data//GISAID//GISAID_countries_regions_withlatlong.csv", row.names=F)

# TO DO: geocode administrative_area_level_1 & administrative_area_level_2 to latitude & longitude


# saveRDS(GISAID_ftable_countries_regions_latlongs, ".//data//GISAID//GISAID_ftable_countries_regions_latlongs.rds")
# saveRDS(GISAID_ftable_revgeo, ".//data//GISAID//GISAID_ftable_revgeo.rds")
# saveRDS(GISAID_ftable_countries_regions, ".//data//GISAID//GISAID_ftable_countries_regions.rds")

GISAID_ftable_countries_regions_latlongs = readRDS(".//data//GISAID//GISAID_ftable_countries_regions_latlongs.rds")
GISAID_ftable_revgeo = readRDS(".//data//GISAID//GISAID_ftable_revgeo.rds")
GISAID_ftable_countries_regions = readRDS(".//data//GISAID//GISAID_ftable_countries_regions.rds")


# reverse geocode to administrative_area_level_1 / state/province and country using sp / rgdal & Natural Earth state map & rworldxtra world map


# countries_map = readOGR(dsn=".//data//GISAID//countries_map", layer="ne_10m_admin_0_countries", use_iconv = TRUE, encoding = "UTF-8") # from https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip
countries_map = getMap(resolution='high') # from rworldxtra 
countries_map@data$ISO3 = as.character(countries_map@data$ISO3)
countries_map@data$ISO3[countries_map@data$ADMIN=="South Sudan"] = "SSD" # we use SSD instead of SDS for South Sudan as in dataset The Economist
countries_map@data$ISO3 = factor(countries_map@data$ISO3)
countries_map@data$REGION = as.character(countries_map@data$REGION)
countries_map@data$REGION[countries_map@data$ADMIN=="Clipperton Island"] = "Australia"
countries_map@data$REGION[countries_map@data$ADMIN=="Heard Island and McDonald Islands"] = "Antarctica"
countries_map@data$REGION[countries_map@data$ADMIN=="French Southern and Antarctic Lands"] = "Antarctica"
countries_map@data$REGION[countries_map@data$REGION=="Australia"] = "Oceania"
countries_map@data$REGION = factor(countries_map@data$REGION)

states_map = readOGR(dsn=".//data//GISAID//states_map", layer="ne_10m_admin_1_states_provinces", use_iconv = TRUE, encoding = "UTF-8") # from https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip
#  states_map@data$adm0_a3 = isoc3 code
# unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3]
# states_map@data[states_map@data$adm0_a3 %in% unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3],]
# states_map@data[states_map@data$adm0_a3 %in% unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3],"name"]
# states_map@data[states_map@data$adm0_a3 %in% unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3],"geonunit"]
# states_map@data[states_map@data$adm0_a3 %in% unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3],"adm0_a3"]
# unique(states_map@data[states_map@data$adm0_a3 %in% unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3],"adm0_a3"])
# # "PSX" "SDS" "ALD" "PGA" "ATC"

# "Argentina" "Uruguay"   "Indonesia" "Malaysia"  "Chile"

# unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3]
# countrycode(unique(states_map@data$adm0_a3)[!unique(states_map@data$adm0_a3) %in% countries_map@data$ISO3], "iso3c", "country.name")
# states_map@data$type_en # Province/Department/State/Region/...
# states_map@data$latitude
# states_map@data$longitude

# country name & iso & continent dictionary to use, with latitude, longitude & population (from rworldxtra world map)
custdict = data.frame(country.name=as.character(countries_map@data$ADMIN), 
                      iso3c=as.character(countries_map@data$ISO3), 
                      iso2c=as.character(countries_map@data$ISO_A2),
                      continent=as.character(countries_map@data$REGION),
                      latitude=countries_map@data$LAT,
                      longitude=countries_map@data$LON,
                      population=countries_map@data$POP_EST)
custdict = rbind(custdict,
                 data.frame(country.name=c("Caribbean Netherlands", "Tokelau"),
                            iso3c=c("BES","TKL"),
                            iso2c=c("BQ", "TK"),
                            continent=c("South America and the Caribbean", "Oceania"),
                            latitude=c(12.178361, -8.967363),
                            longitude=c(-68.238534, -171.8559),
                            population=c(25157, 1411)))

# ISO codes & country names as used in rworldxtra countries_map = getMap(resolution='high')
country_to_iso3c = function (country.names) countrycode(country.names, "country.name", "iso3c", custom_dict=custdict)
iso3c_to_country = function (iso3c) countrycode(iso3c, "iso3c", "country.name", custom_dict=custdict)

                                                     custom_match = c('Ashm' = 'Ashmore and Cartier Islands',
                                                                      'CLP' = 'Clipperton Island',
                                                                      'CNM' = 'Cyprus No Mans Area',
                                                                      'CSI' = 'Coral Sea Islands',
                                                                      'CYN' = 'Northern Cyprus',
                                                                      'ESB' = 'Dhekelia Sovereign Base Area',
                                                                      'Gaza' = 'Gaza',
                                                                      'IOA' = 'Indian Ocean Territories',
                                                                      'KAB' = 'Baykonur Cosmodrome',
                                                                      'KAS' = 'Kashmir',
                                                                      'KNM' = 'Korea No Mans Area',
                                                                      'KOS' = 'Kosovo',
                                                                      'MAF' = 'Saint Martin',
                                                                      'SAH' = 'Western Sahara',
                                                                      'SDS' = 'South Sudan',
                                                                      'SOL' = 'Somaliland',
                                                                      'USG' = 'US Naval Base Guantanamo Bay',
                                                                      'WSB' = 'Akrotiri Sovereign Base Area'))

iso3c_economist[!iso3c_economist %in% countries_map@data$ISO3] # "BES", "SSD" and "TKL" missing

countries_map_centroids = data.frame(as.data.frame(gCentroid(countries_map, byid=TRUE)), 
                                                   iso3c = countries_map@data$ISO3,
                                                   country = iso3c_to_country(countries_map@data$ISO3),
                                                   country_rworldxtra = countries_map@data$ADMIN) 
countries_map_centroids$country[is.na(countries_map_centroids$country)] = countries_map_centroids$country_natearth[is.na(countries_map_centroids$country)]
countries_map_centroids$iso3c[countries_map_centroids$iso3c=="-99"] = NA
countries_map_centroids$iso3c[countries_map_centroids$country=="Somaliland"] = "SOL"
countries_map_centroids$iso3c[countries_map_centroids$country=="Kosovo"] = "KSV"
colnames(countries_map_centroids)[1:2] = c("longitude","latitude")
states_map_centroids = data.frame(as.data.frame(gCentroid(states_map, byid=TRUE)), 
                                  state=states_map@data$name_en) # or name, also had latitude longitude
colnames(states_map_centroids)[1:2] = c("longitude","latitude")

message("Reverse geocoding to administrative_area_level_1 & country using Natureal Earth maps")

library(rworldmap)
library(rworldxtra)

countriesSP <- getMap(resolution='low')


coords2country = function(points)
{  

  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}


# function to reverse geocode to country based on coordinates, see https://rstudio-pubs-static.s3.amazonaws.com/619723_d259330340dc4fc9820c227f5f6291b4.html
get_revgeo_natearth = function(locations){
  lats = locations[,1]
  longs = locations[,2]
  
  # First the coordinates are transformed to spatialpoints
  points = SpatialPoints(matrix(c(longs,
                                 lats),ncol=2))
  # Creating a projection of the coordinates on the map of countries
  proj4string(points) = suppressWarnings(proj4string(countries_map))
  # To see where the name of the country is stored in the map object, you need to explore it in R and see the “data” element. In this case, “NAME” has the information that we want. The function over returns the name of the country given the coordinates projected in the countries_map
  res = over(points, countries_map)
  iso3c = as.character(res$ISO_A3) # NAME_EN
  iso3c[iso3c=="-99"] = NA
  country = countrycode(iso3c, "iso3c", "country.name")
  country_natearth = res$NAME_EN
  country[is.na(country)] = country_natearth[is.na(country)]
  iso3c[country=="Somaliland"] = "SOL"
  iso3c[country=="Kosovo"] = "KSV"
  
  # The same for state
  proj4string(points) = suppressWarnings(proj4string(states_map))
  state = as.character(over(points, states_map)$name_en)
  if (length(state)==0) state=NA
  
  return(data.frame(row=1:nrow(locations), state=state, iso3c=iso3c, country=country, country_natearth=country_natearth)) 
}
# revgeo_natearth(matrix(c(50.78102, 5.464813), nrow=1)) # Limburg
# revgeo_natearth(matrix(c(33.25888, -86.82953), nrow=1)) # Alabama # for US admin1 is at state level
# revgeo_natearth(matrix(c(53.383331, -1.466667), nrow=1)) # Sheffield # natural earth does admin1 for UK at county level
# revgeo_natearth(matrix(c(0.5998746, 37.79594437), nrow=1))

GISAID_ftable_revgeo_natearth = get_revgeo_natearth(locations=data.frame(GISAID_ftable_countries_regions$latitude_region,
                                                                     GISAID_ftable_countries_regions$longitude_region))

GISAID_ftable_countries_regions$state_revgeo_natearth = GISAID_ftable_revgeo_natearth$state
GISAID_ftable_countries_regions$country_revgeo_natearth = countries_map_centroids$country[match(GISAID_ftable_revgeo_natearth$iso3c, countries_map_centroids$iso3c)]

write.csv(GISAID_ftable_countries_regions, ".//data//GISAID//GISAID_countries_regions_withlatlong.csv", row.names=F)

# countries for which to aggregate by state/province/district
# (rest will be aggregated per country using country centroid as latitude/longitude)
countries_per_state = c("South Africa", "India", "Brazil", "United Kingdom", "Canada", "United States", "Australia", "Argentina", "China", "Russia")
# GISAID_ftable_revgeo_natearth looks good for these
# United Kingdom?



# then aggregate by 
# VARIANT, week & country
# VARIANT, week & administrative_area_level_1
# VARIANT, week & administrative_area_level_2
# VARIANT, week & locality

# SUBVARIANT, week & country
# SUBVARIANT, week & administrative_area_level_1
# SUBVARIANT, week & administrative_area_level_2
# SUBVARIANT, week & locality





# TO DO: try to geocode this using
# devtools::install_github("hrbrmstr/nominatim")
# library(nominatim)
# OR
# devtools::install_github(repo = 'rCarto/photon')  
# https://github.com/rCarto/photon
# OR
# https://jessecambon.github.io/tidygeocoder/index.html
# install.packages('tidygeocoder')
library(tidygeocoder)
# ?api_parameter_reference
Sys.setenv(GOOGLEGEOCODE_API_KEY = "AIzaSyADpSYvVxz0FOTF7yVPO6dwYHP6Qrgx2aI") # my Google API key from https://mapsplatform.google.com/, twenseleers@gmail.com 
tidygeocoder::geo(address = "Tongeren, Belgium", method = "osm")

tidygeocoder::geo(address = "New York, USA", method = "osm")
tidygeocoder::get_api_query(address = "New York, USA", method = "google")
tidygeocoder::geocode(address = "New York, USA", method = "google")

tidygeocoder::geo(address = "New York, USA", method = "osm")

res = ggmap::revgeocode(location = c(-97.358112, 37.683829), output="all")
geocode_address(res)

# devtools::install_github("SymbolixAU/googleway")
library(googleway)
apiKey = "AIzaSyADpSYvVxz0FOTF7yVPO6dwYHP6Qrgx2aI"
google_geocode(address = "Leuven, Belgium", key = apiKey)$results$geometry$location
google_geocode(address = "Adrar, Algeria", key = apiKey)$results$geometry$location


library(googleway)
res = googleway::google_reverse_geocode(
  location = c(50.88229, 4.713764), # c(40.7, -74.0), 
  key = apiKey
)
addr_comps = geocode_address_components(res)
addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="locality",] # Leuven
addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="administrative_area_level_1",] # Vlaams Gewest
addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="administrative_area_level_2",] # Vlaams-Brabant
addr_comps[sapply(addr_comps$types, function (x) x[[1]])=="country",] # Belgium




library(dplyr)
df3 <- data.frame(longitude=-74.0, latitude=40.7) %>%
  rowwise() %>%
  do(revgeocode(c(.$longitude[1], .$latitude[1]), output = "all")) %>%
  ungroup()
df3




# and then reverse geocode to provinces/states (administrative_area_level_1) using
# sp / rgdal
# as in https://rstudio-pubs-static.s3.amazonaws.com/619723_d259330340dc4fc9820c227f5f6291b4.html
# using https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip
# or using
# https://stackoverflow.com/questions/46515102/extract-address-components-from-coordiantes
# https://stackoverflow.com/questions/8751497/latitude-longitude-coordinates-to-state-code-in-r

library(ggmap)
geo_information = revgeocode(c(50.88229, 4.713764), output = "all")
geo_information$administrative_area_level_1
# Google API key from https://mapsplatform.google.com/
# twenseleers@gmail.com 
# API KEY AIzaSyADpSYvVxz0FOTF7yVPO6dwYHP6Qrgx2aI
register_google(key="AIzaSyADpSYvVxz0FOTF7yVPO6dwYHP6Qrgx2aI", account_type="premium", write=TRUE)


# using MapQuest Nominatim / OpenStreetMap
# devtools::install_github("hrbrmstr/nominatim")
# https://github.com/hrbrmstr/nominatim
library(nominatim)
osm_search("Alabama", country_codes="us", key="JrURMuiggQT2rfgKETl0BKJrK5KFECTg", limit=1) # my username TWenseleers, password Tompie99, tom.wenseleers@kuleuven.be, https://developer.mapquest.com/user/me/apps
reverse_geocode_coords(lat=33.25888, lon=-86.82953, zoom=6, email="tom.wenseleers@kuleuven.be", key="JrURMuiggQT2rfgKETl0BKJrK5KFECTg")  

osm_search("Tongeren", country_codes="be", key="JrURMuiggQT2rfgKETl0BKJrK5KFECTg", limit=1) # my username TWenseleers, password Tompie99, tom.wenseleers@kuleuven.be, https://developer.mapquest.com/user/me/apps
reverse_geocode_coords(lat=50.78102, lon=5.464813, zoom=6, email="tom.wenseleers@kuleuven.be", key="JrURMuiggQT2rfgKETl0BKJrK5KFECTg")  

reverse_geocode_coords(lat=50.88229, lon=4.713764, zoom=6, email="tom.wenseleers@kuleuven.be", accept_language="en", key="JrURMuiggQT2rfgKETl0BKJrK5KFECTg")


# library(revgeo) # no good
# revgeo(longitude=5.464813, 
#        latitude=50.78102, 
#        provider = 'photon', output="frame")

library(sp)
library(rgdal)
states_map <- readOGR(dsn=".//data//GISAID//states_map", layer="ne_10m_admin_1_states_provinces") # from https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip

# This is a function to reverse geocode based on coordinates, see https://rstudio-pubs-static.s3.amazonaws.com/619723_d259330340dc4fc9820c227f5f6291b4.html
rev_geo<-function(lat,long){
  #First the coordinates are transformed to spatialpoints
  points<-SpatialPoints(matrix(c(long,
                                 lat),ncol=2,nrow=1))
  # #Creating a projection of the coordinates on the map of countries
  # proj4string(points) <- proj4string(countries_map)
  # #To see where the name of the country is stored in the map object, you need to explore it in R and see the “data” element. In this case, “NAME” has the information that we want. The function over returns the name of the country given the coordinates projected in the countries_map
  # country<-as.character(over(points, countries_map)$NAME)
  
  #The same for state
  proj4string(points) <- suppressWarnings(proj4string(states_map))
  state <- as.character(over(points, states_map)$name)
  
  # #The same for LGA (I have only the map for NSW LGAs)
  # proj4string(points) <- proj4string(lgas_map)
  # LGA<-as.character(over(points, lgas_map)$NSW_LGA__3)
  
  return(as.vector(c(state))) # c(country,state,LGA)
}
rev_geo(50.78102, 5.464813) # Limburg
rev_geo(33.25888, -86.82953) # Alabama
rev_geo(53.383331, -1.466667) # Sheffield






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







