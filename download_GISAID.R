# DOWNLOAD GISAID METADATA 
# T. Wenseleers, 20 JULY 2022

library(RSelenium)
usr="u0008355" # put GISAID username here
psw="Tompie99" # put GISAID password here
target_dir = "C:/Users/bherr/Documents/Github/newcovid_belgium/data/GISAID" # target download directory

message("Downloading GISAID metadata...")
browser = wdman::chrome(port = 4570L, version="102.0.5005.61", check=FALSE)
eCaps = list(chromeOptions = list(
  args = c('--headless', # for headless operation, but didn't work for me
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
# supposedly needed in headless mode, but headless operation didn't work for me
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
remDr$navigate("https://www.epicov.org/epi3/start") 
remDr$setImplicitWaitTimeout(milliseconds = 100)

# enter credentials
username = remDr$findElement(using = "xpath", "//input[@id='elogin']")
username$sendKeysToElement(list(usr))
password = remDr$findElement(using = "xpath", "//input[@id='epassword']")
password$sendKeysToElement(list(psw))

# click Login buttom
login_button = remDr$findElement(using = "xpath", "//input[@value='Login']")
login_button$clickElement() 
Sys.sleep(5)

epicov_tab = remDr$findElement("xpath", "//a[contains(text(),'EpiCoV™')]")
epicov_tab$click()
Sys.sleep(1)

downloads_tab = remDr$findElements("class", "sys-actionbar-action-ni")[[3]]
downloads_tab$clickElement()
Sys.sleep(2)

# switch to right frame
frames = remDr$findElements("tag name", "iframe")
remDr$switchToFrame(frames[[1]])
remDr$setImplicitWaitTimeout(milliseconds = 10)

# available download buttons
download_buttons = remDr$findElements("class", "kachel75")
length(download_buttons) # 26 downloads available in total

# DOWNLOAD PATIENT METADATA
metadata_button = download_buttons[[12]] # patient metadata
metadata_button$clickElement()
frames = remDr$findElements("tag name", "iframe")
remDr$switchToFrame(frames[[1]])
remDr$setImplicitWaitTimeout(milliseconds = 10)
checkbox = remDr$findElements("class", "sys-event-hook")[[1]]
checkbox$clickElement() 
Sys.sleep(1)
download = remDr$findElements("class", "sys-form-button-icon")[[2]]
download$clickElement() # download file
Sys.sleep(1)
# wait until download finishes
l = 1
while (l!=0) {
  f = list.files(target_dir, pattern=".crdownload")
  l = length(f)
  if (l==1) download = gsub(".crdownload", "", f, fixed=T) # name of downloaded file
  Sys.sleep(1)
} 
message(paste0("Downloaded GISAID metadata file version ", download))

remDr$close()
browser$stop()
remDr$quit()

