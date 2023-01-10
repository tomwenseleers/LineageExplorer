# DOWNLOAD AGGREGATED nextcladePangoLineage COUNTS FROM COVSPECTRUM WITH CHOICE TO USE OPEN NCBI DATA OR GISAID DATA ####
# T. Wenseleers, 8 JANUARY 2023

# 1. function to download aggregated nextcladePangoLineage counts from covSpectrum ####
download_covSpectrum = function(source="GISAID",  # or "NCBI"
                                date_from="2019-12-01",
                                country=NA, # country, NA for all
                                nextcladePangoLineages=NA, # e.g. "XBB.1.5*", NA for all lineages
                                bydate=TRUE, # by date - TRUE or FALSE
                                bypangolineage=TRUE, # by nextcladePangoLineage?
                                bydivision=FALSE # by division TRUE or FALSE (e.g. by state for USA)
                                ) {
  require(jsonlite)
  
  # set REACT_APP_LAPIS_ACCESS_KEY as environment variable using Sys.setenv(REACT_APP_LAPIS_ACCESS_KEY = "XXXX") 
  # to be able to access COVSPECTRUM private GISAID part, cf. https://cov-spectrum.org/static/js/main.b4ab11d7.js
  if (file.exists("..//set_COVSPECTRUM_credentials.R")) source("..//set_COVSPECTRUM_credentials.R") 
  
  # query the COVSPECTRUM API
  if (source=="NCBI") { base_url = "https://lapis.cov-spectrum.org/open/v1/sample/aggregated?"
                        key = NA
  } else { base_url = "https://lapis.cov-spectrum.org/gisaid/v1/sample/aggregated?"
           key = Sys.getenv("REACT_APP_LAPIS_ACCESS_KEY") # REACT_APP_LAPIS_ACCESS_KEY key to access private GISAID part of COVSPECTRUM
    }
  
  query = paste0(
    base_url,
    "dateFrom=", date_from,
    "&host=Human",
    "&fields=region,country"
  )
  if (bypangolineage) query = paste0(query, ",nextcladePangoLineage")
  if (bydate) query = paste0(query, ",date")
  if (bydivision) query = paste0(query, ",division")
  if (!is.na(nextcladePangoLineages[[1]])) query = paste0(query, "&nextcladePangoLineage=", paste0(nextcladePangoLineages,collapse=",") ) 
  if (!is.na(country[[1]])) query = paste0(query, "&country=", country)
  if (!is.na(key)) query = paste0(query, "&accessKey=", key)
  
  response = fromJSON(URLencode(query))

  # Check for errors
  errors = response$errors
  if (length(errors) > 0) {
    stop("Query returned errors")
  }

  # Check for deprecation
  deprecationDate = response$info$deprecationDate
  if (!is.null(deprecationDate)) {
    warning(paste0("This version of the API will be deprecated on ", deprecationDate,
                 ". Message: ", response$info$deprecationInfo))
  }

  # The data is good to be used!
  data = response$data
  
  return(data) 
}

# e.g. to get country breakdown over past month
# download_covSpectrum(date_from="2022-12-01", country="Sweden", bydate=F) %>% arrange(count)


# 2. function to check countries where at least minseqs of a variant have been detected over the past lastdays days ####
countrieswithvariant = function (target_variant="XBB.1.5*", 
                             minseqs=10,
                             mintotalseqs=100,
                             lastdays=30) { 
  require(dplyr)
  date_from = format(Sys.Date() - as.difftime(lastdays, unit = "days"), "%Y-%m-%d")
  data = download_covSpectrum(source=source, 
                              date_from=date_from,
                              country=NA, # vector of country names, NA for all
                              nextcladePangoLineages=target_variant,
                              bydate=FALSE,
                              bypangolineage=TRUE,
                              bydivision=FALSE # by division YES or NO (e.g. by state for USA)
  )
  data_totals = download_covSpectrum(source=source, 
                                     date_from=date_from,
                                     country=NA, # vector of country names, NA for all
                                     nextcladePangoLineages=NA, # NA for all lineages
                                     bydate=FALSE,
                                     bypangolineage=FALSE,
                                     bydivision=FALSE # by division YES or NO (e.g. by state for USA)
  ) %>% rename(total = count)
  data = left_join(data_totals, data) %>% 
             filter((!is.na(count))&(count>=minseqs)&(total>=mintotalseqs)) %>% 
             mutate(prop_variant = count/total) %>%
             arrange(desc(prop_variant))
  return(data)
}






# PS for info on API see https://lapis.cov-spectrum.org/open/docs/#aggregation 
# https://lapis-docs.readthedocs.io/en/latest/reference/sars_cov_2.html 
# https://lapis-docs.readthedocs.io/en/latest/concepts/variant_query.html
# for all available fields & queries that are possible

# see also https://github.com/GenSpectrum/cov-spectrum-server/blob/develop/docs/api.md
# for pango lineage aliases: fromJSON("https://cov-spectrum.org/api/v2/resource/pango-lineage-alias")
# for countries: fromJSON("https://cov-spectrum.org/api/v2/resource/country")
# for cases & deaths: fromJSON("https://cov-spectrum.org/api/v2/resource/case?country=Belgium&fields=date")
# for collections: fromJSON("https://cov-spectrum.org/api/v2/resource/collection")
# for collection: fromJSON("https://cov-spectrum.org/api/v2/resource/collection/{1}")
