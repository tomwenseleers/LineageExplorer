# DOWNLOAD GISAID METADATA DOWNLOAD PACKAGE ####
# T. Wenseleers, 2 AUGUST 2022

# note: set GISAID access credentials first using ####
# Sys.setenv(GISAIDR_USERNAME = "XXX") # GISAID username
# Sys.setenv(GISAIDR_PASSWORD = "XXX") # GISAID password
# source("..//set_GISAID_credentials.R")

download_GISAD_meta = function(target_dir = getwd(), # target download directory
                               clean_up = FALSE,
                               headless = FALSE,
                               genom_epidem = FALSE, # if TRUE use Genomic Epidemiology Metadata package download, else use regular Metadata package download  
                               chromedriver_version = as.character(unlist(binman::list_versions("chromedriver")))[grepl(as.character(locatexec::exec_version("chrome")[[1, 1]]), as.character(unlist(binman::list_versions("chromedriver"))))][[1]], # "105.0.5195.127" or before "104.0.5112.79"
                               usr = Sys.getenv("GISAIDR_USERNAME"),  
                               psw = Sys.getenv("GISAIDR_PASSWORD")) {
  # TO DO: also implement arguments clean_up=TRUE to delete downloaded file (default best set to FALSE though)
  # and get_sequence=TRUE to also download FASTA with sequences & add those to outputted dataframe
  
  require(RSelenium)
  require(readr)
  require(archive)
  require(data.table)
  require(dplyr)
  require(stringr)
  
  if (!dir.exists(target_dir)) dir.create(target_dir)
                              
arg = c('--no-sandbox', 
        '--disable-dev-shm-usage', 
        '--disable-blink-features=AutomationControlled',
        'user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/61.0.3163.100 Safari/537.36'
)
if (headless) { arg = c('--headless', arg) } # for headless operation
eCaps = list(chromeOptions = list(
  args = arg,
  prefs = list(
    "profile.default_content_settings.popups" = 0L,
    "profile.default_content_setting_values.automatic_downloads" = 1L,
    "download.prompt_for_download" = FALSE,
    "download.directory_upgrade" = TRUE,
    "safebrowsing.enabled" = TRUE,
    "safebrowsing.disable_download_protection" = TRUE,
    "useAutomationExtension" = FALSE,
    "default_directory" = normalizePath(target_dir)
  )
))
browser = wdman::chrome(port = 4570L, version = chromedriver_version, check = TRUE) 
remDr = remoteDriver(port = 4570L, 
                     version=chromedriver_version,
                     browserName = "chrome", 
                     extraCapabilities = eCaps)
remDr$open()
# note: if you get a message This version of ChromeDriver only supports Chrome version xxx
# then make sure that you have specified the right chromedriver version 
# (check if installed chrome version matches some version available under binman/binman_chromedriver/XXX which ones get installed)
# it may be necessary to downgrade your Chrome browser to version xxx using instructions at
# https://browserhow.com/how-to-downgrade-and-install-older-version-of-chrome/#download-the-older-chrome-version
# download browser install from https://filehippo.com/download_google-chrome/history/
# and disable chrome updates



# code below needed in headless mode
# https://stackoverflow.com/questions/35504731/specify-download-folder-in-rselenium
remDr$queryRD(
  ipAddr = paste0(remDr$serverURL, "/session/", remDr$sessionInfo[["id"]], "/chromium/send_command"),
  method = "POST",
  qdata = list(
    cmd = "Page.setDownloadBehavior",
    params = list(
      behavior = "allow",
      downloadPath = normalizePath(target_dir)
    )
  )
)

if (genom_epidem) { message("Checking GISAID genomic epidemiology metadata...") } else { message("Checking GISAID metadata...") }

remDr$navigate("https://www.epicov.org/epi3/start") 
remDr$setImplicitWaitTimeout(milliseconds = 100)

# enter credentials
username = NULL
while (length(username)==0) username = remDr$findElement(using = "xpath", "//input[@id='elogin']")
username$sendKeysToElement(list(usr))
password = remDr$findElement(using = "xpath", "//input[@id='epassword']")
password$sendKeysToElement(list(psw))

# click Login buttom
login_button = remDr$findElement(using = "xpath", "//input[@value='Login']")
login_button$clickElement() 

epicov_tab = NULL
while (length(epicov_tab)==0) { 
  suppressMessages(tryCatch({
    epicov_tab = remDr$findElement("xpath", "//a[contains(text(),'EpiCoV™')]") 
  }, error = function( err ) { epicov_tab = NULL })) }
epicov_tab$click()

downloads_tab = NULL
while (length(downloads_tab)==0) downloads_tab = remDr$findElements("class", "sys-actionbar-action-ni")[[3]]
downloads_tab$clickElement()

# switch to right frame
frames = NULL
while (length(frames)==0) frames = remDr$findElements("tag name", "iframe")
remDr$switchToFrame(frames[[1]])
remDr$setImplicitWaitTimeout(milliseconds = 10)

# available download buttons
download_buttons = NULL
while (length(download_buttons)==0) download_buttons = remDr$findElements("class", "kachel75")
# length(download_buttons) # 27 downloads available in total for me

# sapply(download_buttons, function(d) d$findChildElement("class", "downicon"))

downicons = remDr$findElements("class", "downicon")
downicons_titles = unlist(sapply(downicons, function (d) d$getElementAttribute('title')))

linkicons = remDr$findElements("class", "linkicon")
linkicons_titles = unlist(sapply(linkicons, function (d) d$getElementAttribute('title')))

# TO DO: figure out all titles, in the order in which they appear on the Downloads page
# all_titles = XXX
  
download_nr_gisaid_meta = which(grepl("TSV-File", downicons_titles))+
  sum(linkicons_titles %in% c("Audacity files archive","LANL"))
download_nr_gisaid_genom_epidem_meta = which(grepl("metadata_", downicons_titles)) +
  sum(linkicons_titles %in% c("Audacity files archive", 
                              "LANL",
                              "FASTA and Metadata per clade",
                              "FASTA and Metadata per lineage",      
                              "Global Phylogeny",
                              "Select input for the Augur pipeline"))

# # TO DO : get titles of <div class="downicon" onclick="sys.call('c_rg1kzy_13e','DownloadFile',new Object({'id':'gisaid:metadata_tsv.tar.xz'}));" title="TSV-File (2022-08-01)">
# <img src="/epi3/app_entities/entities/corona2020/download_other2.png"><div>metadata</div>
#   </div>

# check version available for download & if already downloaded don't download it again

# TO DO: allow downloading FASTA or other available download packages

downl_title = downicons_titles[which(grepl("TSV-File", downicons_titles))]
patt = "(?i)((?:(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z ]*\\.?)|(?:\\d{1,2}))[\\/ ,-](\\d{1,2})(?:[\\/ ,-]\\s*(\\d{4}|\\d{2}))?"
available_metadata_todownload = paste0("metadata_tsv_20", gsub("-", "_", str_match_all(downl_title, patt)[[1]][[1]], fixed=T), ".tar.xz")
available_genomepidem_todownload = downicons_titles[which(grepl("metadata_", downicons_titles))]
available_todownload = ifelse(genom_epidem, available_genomepidem_todownload, available_metadata_todownload)

message(paste0("Metadata file version available for download: ", available_todownload))

metadata_already_downloaded = tail(list.files(target_dir, pattern=".tar.xz"),1)
genomepidem_already_downloaded = tail(list.files(target_dir, pattern=".tsv.gz"),1)
already_downloaded = ifelse(genom_epidem, genomepidem_already_downloaded, metadata_already_downloaded)
if (is.na(already_downloaded)) already_downloaded="NA"

message(paste0("Metadata file version already downloaded available in target directory: ", already_downloaded))
  
if (available_todownload != already_downloaded) {
# DOWNLOAD PATIENT METADATA (EITHER REGULAR GISAID METADATA OR GENOMIC EPIDEMIOLOGY METADATA)
download_nr = ifelse(genom_epidem, download_nr_gisaid_genom_epidem_meta, download_nr_gisaid_meta)
metadata_button = download_buttons[[download_nr]] # patient metadata

# TO DO: check available version & if that file is already present in target_dir don't bother downloading it again
metadata_button$clickElement()

# CLICK CHECKBOX NOTICE AND REMINDER OF TERMS OF USE & PRESS DOWNLOAD (THIS ONE DOES NOT ALWAYS SHOW UP)
frames = NULL
while (length(frames)==0) frames = remDr$findElements("tag name", "iframe")
remDr$switchToFrame(frames[[1]])
remDr$setImplicitWaitTimeout(milliseconds = 10)

checkbox_iagree=NULL
while (length(checkbox_iagree)==0) checkbox_iagree = remDr$findElements("class", "sys-event-hook")[[1]]
checkbox_iagree$clickElement() 

Sys.sleep(3)

download_button = NULL
while (length(download_button)==0) download_button = remDr$findElements("class", "sys-form-button")[[2]]
suppressMessages(tryCatch({
  download_button$clickElement() 
}, error = function( err ) { message("") }))

Sys.sleep(3)

# wait until download finishes
while (length(list.files(target_dir, pattern="crdownload", full.names=T))>=1) {
  Sys.sleep(1)
}

if (genom_epidem) pat = ".tsv.gz" else pat = ".tar.xz"
# df = file.info(list.files(target_dir, pattern=pat, full.names = T))
# gsub(paste0(target_dir,"/"), "", rownames(df)[which.max(df$mtime)])
download = tail(list.files(target_dir, pattern=pat, full.names = F), 1) 

message(paste0("Downloaded GISAID metadata file version ", download))

} else { message(paste0("No new metadata package available to download, using previously downloaded metadata available in target directory"))
         download = already_downloaded
}

remDr$close()
browser$stop()
remDr$quit()

message(paste0("Reading GISAID metadata file version ", download))

if (!genom_epidem) { output = read_tsv( # we directly read from archive
  archive_read(file.path(target_dir,download), file=2), 
  col_types = cols(.default = "c"))
# PS fread is slightly faster (multicore), but requires file to be unzipped first
# system.time(archive_extract(archive=file.path(target_dir,download),
#                 dir=target_dir)) # 33s
# system.time(GISAID <- fread(file.path(target_dir,"metadata.tsv")))
colnames(output) = gsub("-", "_", gsub("?", "", gsub(" ", "_", tolower(colnames(output))), fixed=T), fixed=T)
output$pango_lineage = gsub(" (marker override based on Emerging Variants AA substitutions)", "",  output$pango_lineage, fixed=T)
# for me the regular GISAID metadata package download then has the following column names:
# colnames(output)
# [1] "virus_name"                      "type"                           
# [3] "accession_id"                    "collection_date"                
# [5] "location"                        "additional_location_information"
# [7] "sequence_length"                 "host"                           
# [9] "patient_age"                     "gender"                         
# [11] "clade"                           "pango_lineage"                  
# [13] "pangolin_version"                "variant"                        
# [15] "aa_substitutions"                "submission_date"                
# [17] "is_reference"                    "is_complete"                    
# [19] "is_high_coverage"                "is_low_coverage"                
# [21] "n_content"                       "gc_content" 

# remove leading and trailing round brackets from aa_substitutions
output$aa_substitutions = stringr::str_sub(output$aa_substitutions, 2, -2)

# parse location field ####
# parse continent / country / location (PS: location is sometimes city & sometimes province)
loc = do.call(cbind, data.table::tstrsplit(unlist(output$location), "/", TRUE)) # parse location info
loc = trimws(loc, whitespace = " ")
levels_continents = c("Asia","North America","Europe","Africa","South America","Oceania")
output$continent = factor(loc[,1], levels=levels_continents)
output$country = factor(loc[,2])
output$location = factor(loc[,3]) # city or province/state

} else {
  # system.time(output <- read_tsv( # we directly read from archive
  #   gzfile(file.path(target_dir,download)), 
  #   col_types = cols(.default = "c"))) # 61s
  system.time(output <- tibble(fread( # we directly read from archive, fread 2x faster than read_tsv
    file.path(target_dir,download), 
    colClasses = c("character")))) # 28s
  colnames(output) = tolower(colnames(output))
  # the GISAID genomic epidemiology metadata package download has the following column names:
  # TO DO: make some names consistent across metadata files?
  # colnames(output)
  # [1] "strain"                "virus"                 "gisaid_epi_isl"        "genbank_accession"    
  # [5] "date"                  "region"                "country"               "division"             
  # [9] "location"              "region_exposure"       "country_exposure"      "division_exposure"    
  # [13] "segment"               "length"                "host"                  "age"                  
  # [17] "sex"                   "nextstrain_clade"      "pango_lineage"         "gisaid_clade"         
  # [21] "originating_lab"       "submitting_lab"        "authors"               "url"                  
  # [25] "title"                 "paper_url"             "date_submitted"        "purpose_of_sequencing"
} 

colnames(output) = gsub("region", "continent", colnames(output), fixed=TRUE)
colnames(output) = gsub("location", "city", colnames(output), fixed=TRUE)
colnames(output) = gsub("division", "location", colnames(output), fixed=TRUE)

if (clean_up) unlink(file.path(target_dir,download))

return(output)

}

# example
# GISAID = download_GISAD_meta(target_dir = "C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID")
