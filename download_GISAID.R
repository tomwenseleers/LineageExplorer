# DOWNLOAD GISAID METADATA DOWNLOAD PACKAGE ####
# T. Wenseleers, 2 AUGUST 2022

# note: set GISAID access credentials first using ####
# Sys.setenv(GISAIDR_USERNAME = "XXX") # GISAID username
# Sys.setenv(GISAIDR_PASSWORD = "XXX") # GISAID password
# source(".//set_GISAID_credentials.R")

download_GISAD_meta = function(target_dir = getwd(), # target download directory
                               usr=Sys.getenv("GISAIDR_USERNAME"),  
                               psw=Sys.getenv("GISAIDR_PASSWORD")) {
  # TO DO: also implement arguments clean_up=TRUE to delete downloaded file (default best set to FALSE though)
  # and get_sequence=TRUE to also download FASTA with sequences & add those to outputted dataframe
  
  require(RSelenium)
  require(readr)
  require(archive)
  
  if (!dir.exists(target_dir)) dir.create(target_dir)
                              
browser = wdman::chrome(port = 4570L, version="102.0.5005.61", check=FALSE)
eCaps = list(chromeOptions = list(
  args = c('--headless', # for headless operation
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
    "default_directory" = gsub("/", "\\", target_dir, fixed=T)
  )
))
remDr = remoteDriver(port = 4570L, 
                     version="102.0.5005.61", 
                     browserName = "chrome", 
                     extraCapabilities = eCaps)
remDr$open()
# code below needed in headless mode
# https://stackoverflow.com/questions/35504731/specify-download-folder-in-rselenium
remDr$queryRD(
  ipAddr = paste0(remDr$serverURL, "/session/", remDr$sessionInfo[["id"]], "/chromium/send_command"),
  method = "POST",
  qdata = list(
    cmd = "Page.setDownloadBehavior",
    params = list(
      behavior = "allow",
      downloadPath = gsub("/", "\\", target_dir, fixed=T)
    )
  )
)

message("Downloading GISAID metadata...")

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
Sys.sleep(6)

epicov_tab = remDr$findElement("xpath", "//a[contains(text(),'EpiCoV™')]")
epicov_tab$click()
Sys.sleep(1)

downloads_tab = remDr$findElements("class", "sys-actionbar-action-ni")[[3]]
downloads_tab$clickElement()
Sys.sleep(3)

# switch to right frame
frames = NULL
while (length(frames)==0) { frames = remDr$findElements("tag name", "iframe") }
remDr$switchToFrame(frames[[1]])
remDr$setImplicitWaitTimeout(milliseconds = 10)

# available download buttons
download_buttons = NULL
while (length(download_buttons)==0) { download_buttons = remDr$findElements("class", "kachel75") }
# length(download_buttons) # 26 downloads available in total

# DOWNLOAD PATIENT METADATA
metadata_button = download_buttons[[12]] # patient metadata
# TO DO: check available version & if that file is already present in target_dir don't bother downloading it again
metadata_button$clickElement()
Sys.sleep(1)

frames = NULL
while (length(frames)==0) { frames = remDr$findElements("tag name", "iframe") }
remDr$switchToFrame(frames[[1]])
remDr$setImplicitWaitTimeout(milliseconds = 10)
checkbox = remDr$findElements("class", "sys-event-hook")[[1]]
checkbox$clickElement() 
Sys.sleep(1)

download_button = remDr$findElements("class", "sys-form-button")[[2]]
download_button$clickElement() 
Sys.sleep(15)

# press OK to NOTICE AND REMINDER OF TERMS OF USE (THIS ONE DOES NOT ALWAYS SHOW UP)
download_button = remDr$findElements("class", "sys-form-button")[[2]]
suppressMessages(tryCatch({
  download_button$clickElement() 
}, error = function( err ) { message("") }))
Sys.sleep(1)
# wait until download finishes
while (length(list.files(target_dir, pattern="crdownload", full.names=T))>=1) {
  Sys.sleep(1)
}
df = file.info(list.files(target_dir, pattern=".tar.xz", full.names = T))
download = gsub(paste0(target_dir,"/"), "", rownames(df)[which.max(df$mtime)])

message(paste0("Downloaded GISAID metadata file version ", download))

remDr$close()
browser$stop()
remDr$quit()

message(paste0("Reading GISAID metadata file version ", download))

return(read_tsv( # to directly read from archive
  archive_read(paste0(target_dir, "//", download), file=2), 
  col_types = cols(.default = "c")))

}

# example
# GISAID = download_GISAD_meta(target_dir = "C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID")
