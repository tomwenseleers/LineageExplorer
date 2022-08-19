# DOWNLOAD GISAID RECORDS ####
# T. Wenseleers, 31 JULY 2022

# note: set GISAID access credentials first using ####
# Sys.setenv(GISAIDR_USERNAME = "XXX") # GISAID username
# Sys.setenv(GISAIDR_PASSWORD = "XXX") # GISAID password
# source(".//set_GISAID_credentials.R")

library(locatexec)

# function to download given records, after having logged in first
# (max batch of 10 000 records at a time)
downl_records = function (accession_ids,
                          get_sequence=FALSE, 
                          clean_up=FALSE,
                          target_dir,
                          remDr=remDr) {
  require(readr)
  require(dplyr)
  require(stringr)
  require(dplyr)
  
  remDr$refresh()
  
  # click Select tab at the bottom
  select_tab = NULL
  while(length(select_tab)==0) select_tab = remDr$findElements("class", "sys-form-button")[[2]]
  select_tab$clickElement()
  
  # enter desired GISAID access nrs in input field 
  frames=NULL
  while (length(frames)==0) frames = remDr$findElements("tag name", "iframe")
  remDr$switchToFrame(frames[[1]])
  remDr$setImplicitWaitTimeout(milliseconds = 10)
  
  # write accession IDs to file & select them via Choose file...
  IDfile = file.path(target_dir,"IDs.csv")
  write.csv(accession_ids, IDfile, row.names=F) # make temporary file with IDs
  # clear input field
  input_field = remDr$findElements("class", "sys-form-fi-multiline")[[2]]
  input_field$clickElement()
  input_field$clearElement()
  
  remDr$refresh()
  
  frames=NULL
  while (length(frames)==0) frames = remDr$findElements("tag name", "iframe")
  remDr$switchToFrame(frames[[1]])
  remDr$setImplicitWaitTimeout(milliseconds = 10)

  frames=NULL
  while (length(frames)==0) frames = remDr$findElements("tag name", "iframe")
  remDr$switchToFrame(frames[[1]])
  remDr$setImplicitWaitTimeout(milliseconds = 10)
  
  file_field = remDr$findElements("name", "data")[[1]]
  file_field$sendKeysToElement(list(IDfile))
  
  unlink(IDfile) # remove temporary file with IDs again
  
  # # put GISAID_ids in entry form
  # # this approach only worked for small queries, so doing this via file as above 
  # input_field = remDr$findElements("class", "sys-form-fi-multiline")[[2]]
  # input_field$clickElement()
  # input_field$clearElement()
  # input_field$sendKeysToElement(list(paste0(accession_ids, collapse=", "))) 

  Sys.sleep(1)
  
  # click OK
  remDr$refresh()
  
  frames=NULL
  while (length(frames)==0) frames = remDr$findElements("tag name", "iframe")
  remDr$switchToFrame(frames[[1]])
  remDr$setImplicitWaitTimeout(milliseconds = 10)
  
  OK_button = remDr$findElements("class", "sys-form-button")[[2]]
  OK_button$clickElement()
  Sys.sleep(1) # TO DO change to while loop
  
  # click OK again to Message XXX entries selected
  buttonid = ""
  while (buttonid=="") {
    html = remDr$getPageSource()[[1]]
    buttonid = str_extract(html, "(?<=button id=\")[0-9]*") # button id appears dynamic
    # print(buttonid)
  }
  OK_button = remDr$findElement("id", buttonid)
  OK_button$clickElement()
  
  # click Download at the bottom
  remDr$refresh()
  download_tab = remDr$findElements("class", "sys-form-button")[[4]]
  download_tab$clickElement()
  # Sys.sleep(20) # TO DO change to while loop
  
  # select Patient status metadata checkbox
  frames=NULL
  while (length(frames)==0) {
    suppressMessages(tryCatch({
      frames = remDr$findElements("tag name", "iframe")
    }, error = function( err ) { frames = NULL }))
    }
  remDr$switchToFrame(frames[[1]])
  remDr$setImplicitWaitTimeout(milliseconds = 10)
  
  Sys.sleep(1)
  checkbox_metadata = NULL
  while (length(checkbox_metadata)==0) { 
    checkbox_metadata = remDr$findElements("class", "sys-event-hook")[[4]] }
  checkbox_metadata$clickElement() 
  
  # click Download
  mostrecenttsv = suppressWarnings(max(file.info(list.files(target_dir, pattern=".tsv", full.names=T))$mtime))
  Sys.sleep(1) # TO DO change to while loop
  mostrecenttsv_new = mostrecenttsv
  download_button = NULL
  while (length(download_button)==0) {
    suppressMessages(tryCatch({
      download_button = remDr$findElements("class", "sys-form-button")[[2]]
      download_button$clickElement() 
    }, error = function( err ) { download_button = NULL }))
  } 
  
  # CLICK CHECKBOX NOTICE AND REMINDER OF TERMS OF USE & PRESS DOWNLOAD (THIS ONE DOES NOT ALWAYS SHOW UP)
  frames=NULL
  while (length(frames)==0) {
    suppressMessages(tryCatch({
      frames = remDr$findElements("tag name", "iframe")
    }, error = function( err ) { frames = NULL }))
  }
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
  Sys.sleep(1)
  
  # wait until download finishes
  while (mostrecenttsv_new==mostrecenttsv) {
    mostrecenttsv_new = suppressWarnings(max(file.info(list.files(target_dir, pattern=".tsv", full.names=T))$mtime))
    Sys.sleep(1)
  }
  df = file.info(list.files(target_dir, pattern=".tsv", full.names = T))
  download = rownames(df)[which.max(df$mtime)]
  
  message(paste0("Downloaded GISAID metadata file ", download))
  
  # # TO DO finish part below if get_sequence=TRUE
  # # click Download again at the bottom to select FASTA sequences & download those
  # # if desired
  # remDr$refresh()
  # download_tab = remDr$findElements("class", "sys-form-button")[[4]]
  # download_tab$clickElement()
  # 
  # # select Nucleotide sequences (FASTA) checkbox
  # frames = remDr$findElements("tag name", "iframe")
  # remDr$switchToFrame(frames[[1]])
  # remDr$setImplicitWaitTimeout(milliseconds = 10)
  # 
  # checkbox_FASTA = remDr$findElements("class", "sys-event-hook")[[3]]
  # checkbox_FASTA$clickElement() 
  # 
  # # click Download
  # mostrecentfasta = max(file.info(list.files(target_dir, pattern=".fasta", full.names=T))$mtime)
  # mostrecentfasta_new = mostrecenttsv
  # download_button = remDr$findElements("class", "sys-form-button")[[2]]
  # download_button$clickElement() 
  # 
  # Sys.sleep(1)
  # while (mostrecentfasta_new==mostrecentfasta) {
  #   mostrecentfasta_new = max(file.info(list.files(target_dir, pattern=".fasta", full.names=T))$mtime)
  #   Sys.sleep(1)
  # }
  # df = file.info(list.files(target_dir, pattern=".fasta", full.names = T))
  # download_fasta = gsub(paste0(target_dir,"/"), "", rownames(df)[which.max(df$mtime)])
  # 
  # message(paste0("Downloaded GISAID FASTA file ", download_fasta))
  
  # read in .tsv file download
  output = read_tsv(download, 
                    col_types = cols(.default = "c"))
  
  colnames(output) = gsub(" ", "_", tolower(colnames(output)))
  colnames(output)[which(colnames(output) %in% c("lineage"))] = "pango_lineage" # code lineage as pango_lineage as in GISAID metadata download package
  output$pango_lineage = gsub(" (marker override based on Emerging Variants AA substitutions)", "",  output$pango_lineage, fixed=T)
  
  if (clean_up) unlink(download)
  
  return(output)
}


# function to download given records
# (will be split in batches of max 10 000 records each)
download_GISAID_records = function(
                            accession_ids,
                            get_sequence=FALSE, 
                            clean_up=FALSE,
                            target_dir="C:/Users/bherr/OneDrive - KU Leuven/Documents/Github/newcovid_belgium/data/GISAID",
                            max_batch_size=10000, # maximum batch size
                            headless = FALSE,
                            chromedriver_version = as.character(unlist(binman::list_versions("chromedriver")))[grepl(as.character(locatexec::exec_version("chrome")[[1, 1]]), as.character(unlist(binman::list_versions("chromedriver"))))], # "104.0.5112.79"
                            usr=Sys.getenv("GISAIDR_USERNAME"),
                            psw=Sys.getenv("GISAIDR_PASSWORD")) {
require(RSelenium)  

if (!dir.exists(target_dir)) dir.create(target_dir)  

arg = c('--no-sandbox', 
        '--disable-dev-shm-usage', 
        '--disable-blink-features=AutomationControlled',
        'user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/61.0.3163.100 Safari/537.36'
)
if (headless) { arg = c('--headless', arg) } # for headless operation

eCaps = list(chromeOptions = list(
  args = c(#'--headless', # for headless operation
           '--no-sandbox', 
           '--disable-dev-shm-usage', 
           '--disable-blink-features=AutomationControlled',
           'user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/61.0.3163.100 Safari/537.36'
  ),
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
                     version = chromedriver_version, 
                     browserName = "chrome", 
                     extraCapabilities = eCaps)
remDr$open()
# code below needed for headless operation
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

message("Downloading GISAID records...")

remDr$navigate("https://www.epicov.org/epi3/start") 
remDr$setImplicitWaitTimeout(milliseconds = 100)

# enter credentials
username = NULL
while (length(username)==0) { username = remDr$findElement(using = "xpath", "//input[@id='elogin']") }
username$sendKeysToElement(list(usr))
password = remDr$findElement(using = "xpath", "//input[@id='epassword']")
password$sendKeysToElement(list(psw))

# click Login buttom
login_button = remDr$findElement(using = "xpath", "//input[@value='Login']")
login_button$clickElement() 

# click EpiCov tab
epicov_tab = NULL
while (length(epicov_tab)==0) { 
  suppressMessages(tryCatch({
    epicov_tab = remDr$findElement("xpath", "//a[contains(text(),'EpiCoV™')]") 
  }, error = function( err ) { epicov_tab = NULL })) }
epicov_tab$click()

# click Search tab
search_tab = NULL
while (length(search_tab)==0) search_tab = remDr$findElements("class", "sys-actionbar-action-ni")[[2]]
search_tab$clickElement()
# Sys.sleep(10) # TO DO add while loop in downl_records

# download records in batches of maximum size max_batch_size

# function to split vector in chunks of max size chunk_length
chunk = function(x, chunk_length=max_batch_size) split(x, ceiling(seq_along(x)/chunk_length))

batches = chunk(accession_ids)

downloads = do.call(bind_rows, lapply(1:length(batches),
                                  function (batchnr) {
                                    message(paste0("Downloading batch ", batchnr, " out of ", length(batches)))
                                    Sys.sleep(1)
                                    output = downl_records(accession_ids=batches[[batchnr]],
                                                                       get_sequence, 
                                                                       clean_up,
                                                                       target_dir,
                                                                       remDr)
                                    return(output) } ))

remDr$close()
browser$stop()
remDr$quit()

return(downloads)

}

# example
# d=download_GISAID_records(accession_ids=c("EPI_ISL_14041928", "EPI_ISL_14042096", "EPI_ISL_14041961"),
#                           max_batch_size=2)
# dim(d)
