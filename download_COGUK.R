# DOWNLOAD COG-UK METADATA ####
# T. Wenseleers, 13 AUGUST 2022

# function to download & read COG-UK metadata

download_COGUK_meta = function() {
  require(data.table)
  require(dplyr)
  
  message("Downloading & reading COG-UK metadata...")
  
  coguk = tibble(fread("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz"))
  # rename some variables as in GISAID metadata
  coguk = coguk %>% rename(collection_date = sample_date,
                           pango_lineage = lineage,
                           aa_substitutions = mutations)
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

return(coguk)
}
