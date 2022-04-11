library(rgeos)
library(rworldmap)

# get world map
wmap <- getMap(resolution="high")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with longitudes & latitudes of country centroids
countries_long_lat <- as.data.frame(centroids)
colnames(countries_long_lat) = c("longitude","latitude")
countries_long_lat$country_rworldmap = rownames(countries_long_lat)
rownames(countries_long_lat) = NULL
countries_long_lat$iso3c = countrycode(countries_long_lat$country_rworldmap, "country.name", "iso3c",
                                       custom_match = c('Ashmore and Cartier Islands' = NA,
                                                          'Gaza' = NA, # just counting West Bank under Palestinian Territories / PSE
                                                          'Indian Ocean Territories' = NA,
                                                          'Kosovo' = 'KSV',
                                                          'Northern Cyprus' = NA, # just counting the rest of Cyprus as Cyprus
                                                          'Saint Martin' = NA, # just counting Sint Maarten as SXM
                                                          'Siachen Glacier' = NA,
                                                          'Somaliland' = 'SOL'))
countries_long_lat = countries_long_lat[!is.na(countries_long_lat$iso3c), ]
countries_long_lat$country = countrycode(countries_long_lat$iso3c, "iso3c", "country.name", custom_match = c('SOL' = 'Somaliland',
                                                                                                             'KSV' = 'Kosovo'))

# rownames(table(countries_long_lat$iso3c))[table(countries_long_lat$iso3c)>1] 
# countries_long_lat[countries_long_lat$iso3c=="SXM",]



# make sure to have longitudes & latitudes for all countries in the analysis of The Economist
# fit_economist = read.csv("https://raw.githubusercontent.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/output-for-interactive/second_map.csv")
# iso3c_economist = sort(unique(fit_economist$iso3c))
# iso3c_economist[!iso3c_economist %in% countries_long_lat$iso3c] # "BES" "GIB" "TKL" missing

countries_long_lat = rbind(countries_long_lat,
                           data.frame(longitude=c(-68.238534, -5.35257, -171.855881),
                                      latitude=c(12.178361, 36.14474, -8.967363),
                                      country_rworldmap=c(NA, NA, NA),
                                      iso3c=c("BES","GIB","TKL"),
                                      country=countrycode(c("BES","GIB","TKL"), "iso3c", "country.name")))

# make sure to have longitudes & latitudes for all countries included in GISAID
# countries$iso3c[!countries$iso3c %in% countries_long_lat$iso3c] # "GLP" "MTQ" "MYT" "REU" missing

countries_long_lat = rbind(countries_long_lat,
                           data.frame(longitude=c(-61.580002, -61.024174, 45.130741, 55.448101),
                                      latitude=c(16.270000, 14.641528, -12.809645, -20.878901),
                                      country_rworldmap=c(NA, NA, NA, NA),
                                      iso3c=c("GLP","MTQ","MYT", "REU"),
                                      country=countrycode(c("GLP","MTQ","MYT", "REU"), "iso3c", "country.name")))

countries_long_lat$country_GISAID = countries$country_GISAID[match(countries_long_lat$iso3c, countries$iso3c)]
countries_long_lat$continent = countries$continent[match(countries_long_lat$iso3c, countries$iso3c)]

write.csv(countries_long_lat, ".//data//GISAID//countries_long_lat_centroids.csv", row.names=F)


