# adapted from https://github.com/JoFAM by Joris Meys

# This script loads the data
source("scripts/checkPackages.R")
thefuns <- dir("functions", pattern = "\\.R")
for(i in thefuns){
  source(file.path("functions",i), encoding="UTF-8")
}


if(!dir.exists("Processed")) dir.create("Processed")
#--------------------------------------

# Read in the raw datasets
message("Start download data. This might take a minute.")
message("Processing cases data.")
## For cases by province, region, gender and agegroup
rawcases<- refresh("cases",
                   "COVID19BE_CASES_AGESEX.csv",
                   process = function(x){
                     add_totals(x,
                                values = "CASES",
                                groups = c("PROVINCE","REGION",
                                           "SEX","AGEGROUP"),
                                along = c("DATE"),
                                name = c("All","Belgium","All","All")) %>%
                       filter(keep_combs(REGION,PROVINCE))
                   })

## For tests by province (and region)
message("Processing test information.")
rawtest <- refresh("tests",
                   "COVID19BE_tests.csv",
                   process = function(x){
                     mutate(x,
                            PROVINCE = add_unknown(PROVINCE),
                            REGION = add_region(PROVINCE)) %>%
                     add_totals(values = c("TESTS_ALL","TESTS_ALL_POS"),
                                groups = c("PROVINCE","REGION"),
                                along = "DATE",
                                name = c("All","Belgium")) %>%
                       filter(keep_combs(REGION,PROVINCE))
                   })

## For hospitalisations by province and region
message("Processing hospitalisations.")
rawhospit <- refresh("hospitalisations",
                     "COVID19BE_HOSP.csv",
                     process = function(x){
                       add_totals(x,
                                  values = c("TOTAL_IN",
                                             "TOTAL_IN_ICU",
                                             "TOTAL_IN_RESP",
                                             "TOTAL_IN_ECMO",
                                             "NEW_IN",
                                             "NEW_OUT"),
                                  groups = c("PROVINCE",
                                             "REGION"),
                                  along = "DATE",
                                  name = c("All","Belgium")) %>%
                         filter(keep_combs(REGION,PROVINCE))
                     })

## For deaths by region, sex and agegroup
message("Processing deaths.")
rawdeaths <- refresh("deaths",
                     "COVID19BE_MORT.csv",
                     process = function(x){
                       add_totals(x,
                                  values = "DEATHS",
                                  groups = c("REGION",
                                             "SEX","AGEGROUP"),
                                  along = c("DATE"),
                                  name = c("Belgium","All","All")) 
                     })

## Numbers for each municipality
message("Processing municipalities data")
rawmunicipalities <- refresh("municipalities",
                             "COVID19BE_CASES_MUNI.csv") %>%
  mutate(DISTRICT = translate_districts(TX_ADM_DSTR_DESCR_NL,
                                        TX_ADM_DSTR_DESCR_NL,
                                        REGION),
         MUNTY = translate_munty(TX_DESCR_NL,
                                 TX_DESCR_FR,
                                 REGION),
         cd_munty_refnis = NIS5,
         binned = bin_cases(CASES)) %>%
  dplyr::select(REGION, PROVINCE, DISTRICT,MUNTY,
         cd_munty_refnis, cases = CASES,
         binned,DATE) %>%
  na.omit()
  
replacefile(rawmunicipalities,
            "binnedmunty",
            "rds",
            "Processed")

#--------------------------------------
# Combine the raw datasets into a number of separate 
# datasets that can be used for analysis.
# Smooth the data using a 7 day window where possible.

## Combine cases, tests, hospitalisations and deaths by date and region
message("Combining data sets and saving.")
# rawcases = rawcases[complete.cases(rawcases),]
casetemp <- filter(rawcases, PROVINCE == "All" & SEX == "All" &
                     AGEGROUP == "All") %>%
  dplyr::select(-c(SEX, PROVINCE, AGEGROUP))

testtemp <- filter(rawtest, PROVINCE == "All") %>%
  dplyr::select(-PROVINCE)

hosptemp <- filter(rawhospit, PROVINCE == "All") %>%
  dplyr::select(-PROVINCE)

deathtemp <- filter(rawdeaths, SEX == "All" & AGEGROUP == "All") %>%
  dplyr::select(-c(SEX, AGEGROUP))

regionalweekavg <- full_join(casetemp,
                             testtemp,
                             by = c("REGION","DATE")) %>%
  full_join(hosptemp, by = c("REGION","DATE")) %>%
  full_join(deathtemp, by = c("REGION","DATE")) %>%
  filter(as.Date(DATE) >= as.Date("2020-03-15")) %>%
  mutate(across(where(is.numeric),replaceby0)) %>%
  group_by(REGION) %>%
  mutate(across(where(is.numeric),
                ~ zoo::rollmean(., 7, align = "right", fill = NA)))  %>%
  mutate(POSRATE = TESTS_ALL/TESTS_ALL_POS,
         CHANGE_CASES = changeabsolute(CASES),
         CHANGE_TESTS = changeabsolute(TESTS_ALL),
         CHANGE_POS = changeabsolute(TESTS_ALL_POS),
         CHANGE_NEW_IN = changeabsolute(NEW_IN),
         CHANGE_DEATHS = changeabsolute(DEATHS),
         CHANGE_TOTAL_IN = changeabsolute(TOTAL_IN),
         CHANGE_TOTAL_IN_ICU = changeabsolute(TOTAL_IN_ICU),
         CHANGE_TOTAL_IN_RESP = changeabsolute(TOTAL_IN_RESP),
         CHANGE_TOTAL_IN_ECMO = changeabsolute(TOTAL_IN_ECMO)) %>%
  as.data.frame()

# Add relative changes
# NEED AGEDIST FOR PROVINCES AS WELL !!!!
# n <- read.csv("Data/AgedistPopBe.csv")
# 
# agecases <- rawcases %>%
#   group_by(DATE, REGION, PROVINCE, AGEGROUP, SEX) %>%
#   mutate(across(where(is.numeric),
#                 ~ zoo::rollmean(., 7, align = "right", fill = NA)))


replacefile(regionalweekavg,
            "regionalweekavg",
            "rds",
            "Processed")
message("Succes!")



# a few extras
# completeness percentages, from Bart Mesuere
completeness = data.frame(WEEKDAY_DAYSBEFOREPRESENT= 
                            c(paste0("Monday.",as.character(1:4)),
                              paste0("Tuesday.",as.character(1:4)),
                              paste0("Wednesday.",as.character(1:4)),
                              paste0("Thursday.",as.character(1:4)),
                              paste0("Friday.",as.character(1:4)),
                              paste0("Saturday.",as.character(1:4)),
                              paste0("Sunday.",as.character(1:4))),
                          COMPLETENESS=c(0.30, 0.86, 0.96, 0.98,
                                         0.40, 0.82, 0.97, 0.98,
                                         0.42, 0.88, 0.97, 0.99,
                                         0.37, 0.89, 0.97, 0.99,
                                         0.34, 0.91, 0.96, 0.98,
                                         0.46, 0.81, 0.98, 0.98,
                                         0.49, 0.88, 0.97, 0.99))
# read as "data from Monday are 30% complete after 2 days, 86% complete after 3 days, 96% complete after 4 days"

# total tests performed
tests_tot = rawtest[rawtest$PROVINCE=="All"&rawtest$REGION=="Belgium",]
tests_tot$DATE = as.Date(tests_tot$DATE)

# tests performed per province
tests_prov = rawtest[rawtest$PROVINCE!="All"&rawtest$REGION!="Belgium",]
tests_prov$DATE = as.Date(tests_prov$DATE)
tests_prov$PROVDATE = interaction(tests_prov$PROVINCE, tests_prov$DATE)


# age distribution of Belgium
agedist_tot = read.csv("./Data/AgedistPopBE.csv")
agedist_tot = agedist_tot[agedist_tot$REGION=="Belgium"&agedist_tot$GENDER=="All",]
agedist_tot$AGEGROUP = as.factor(agedist_tot$AGEGROUP)
agedist_tot$AGENUM = as.numeric(as.character(factor(agedist_tot$AGEGROUP, levels=c("0-9","10-19","20-29",
                                                                                   "30-39","40-49","50-59",
                                                                                   "60-69","70-79","80-89","90+"),
                                                    labels=as.character(4.5+seq(0,90,by=10)))))
agedist_tot$PROP = agedist_tot$POP/sum(agedist_tot$POP)

# infection fatality rate & infection hospitalisation rate as previously estimated for Belgium for period 15 March - 1 May 2020
IFR = function (ages) {  # IFRs as estimated by Molenberghs et al. 2020, https://www.medrxiv.org/content/10.1101/2020.06.20.20136234v1
  splfit = smooth.spline(
    x=c((0+24/2), (25+44)/2, (45+64)/2, (65+74)/2, (75+84)/2, 90),
    y=log(c(0.000005, 0.00017, 0.0021, 0.022, 0.043, 0.12)),
    all.knots=TRUE)
  exp(predict(splfit, ages)$y)  
}

IHR = function (ages) {  # IHRs as estimated by Sciensano
  splfit = smooth.spline(
    x=c((0+24/2), (25+44)/2, (45+64)/2, (65+74)/2, (75+84)/2, 90),
    y=log(c(0.0022, 0.0081, 0.022, 0.061, 0.069, 0.069)),
    all.knots=TRUE)
  exp(predict(splfit, ages)$y)  
}

agedist_tot$IFR = IFR(agedist_tot$AGENUM)
agedist_tot$IHR = IHR(agedist_tot$AGENUM)
totpop = sum(agedist_tot[,3])

d = as.Date(max(rawcases$DATE[rawcases$DATE!="unknown"])) 
# d = format(as.Date(max(rawcases$DATE[rawcases$DATE!="unknown"])), "%b %d %y")
tag = paste("@TWenseleers\n",d)

# case data (total)
cases_tot = rawcases[rawcases$PROVINCE=="All"&
                       rawcases$REGION=="Belgium"&
                       rawcases$SEX=="All"&
                       rawcases$AGEGROUP=="All",]
cases_tot$DATE = as.Date(cases_tot$DATE)
cases_tot = cases_tot[!is.na(cases_tot$DATE),]
range(cases_tot$DATE) # "2020-03-01" "2022-04-11"
cases_tot$DATE_NUM = as.numeric(cases_tot$DATE)
cases_tot$TESTS_ALL = tests_tot$TESTS_ALL[match(cases_tot$DATE,tests_tot$DATE)]
cases_tot$BANKHOLIDAY = bankholiday(cases_tot$DATE)
cases_tot$WEEKDAY = as.factor(weekdays(as.Date(cases_tot$DATE)))
cases_tot$OBS = factor(1:nrow(cases_tot))
cases_tot$DAYSBEFOREPRESENT = as.Date(Sys.time())-cases_tot$DATE
cases_tot$WEEKDAY_DAYSBEFOREPRESENT = interaction(cases_tot$WEEKDAY,cases_tot$DAYSBEFOREPRESENT)
tail(cases_tot,14)
cases_tot$completeness = completeness$COMPLETENESS[match(cases_tot$WEEKDAY_DAYSBEFOREPRESENT,completeness$WEEKDAY_DAYSBEFOREPRESENT)]
cases_tot$completeness[is.na(cases_tot$completeness)] = 1
cases_tot$CASES = round(cases_tot$CASES/cases_tot$completeness,0) # apply correction to most recent data
cases_tot$TESTS_ALL = round(cases_tot$TESTS_ALL/cases_tot$completeness,0) # apply correction to most recent data
cases_tot$CASESPERTEST = cases_tot$CASES/cases_tot$TESTS_ALL
cases_tot$TESTSPERCASES = cases_tot$TESTS_ALL/(cases_tot$CASES+1)
tail(cases_tot)
clip = 2
if (weekdays(today)=="Sunday") clip = clip+1
if (weekdays(today)=="Monday") clip = clip+2
cases_tot = cases_tot[cases_tot$DATE<(as.Date(Sys.time())-clip),] # we leave out last 3-4 days since they are incomplete
cases_tot = cases_tot[complete.cases(cases_tot),]
# tail(cases_tot,14)
# cases_tot = cases_tot[cases_tot$DATE>=as.Date("2020-07-01"),]
# range(cases_tot$DATE) # "2020-03-01" "2022-04-10"

# cases per province
cases_prov = rawcases[rawcases$PROVINCE!="All"&rawcases$REGION!="Belgium"&rawcases$SEX=="All"&rawcases$AGEGROUP!="unknown"&rawcases$AGEGROUP=="All"&rawcases$PROVINCE!="unknown",]
cases_prov$PROVINCE = as.factor(cases_prov$PROVINCE)
cases_prov$AGEGROUP = as.factor(cases_prov$AGEGROUP)
cases_prov$DATE = as.Date(cases_prov$DATE)
cases_prov$DATE_NUM = as.numeric(cases_prov$DATE)
cases_prov$PROVDATE = interaction(cases_prov$PROVINCE,cases_prov$DATE)
cases_prov$TESTS_ALL = tests_prov$TESTS_ALL[match(cases_prov$PROVDATE,tests_prov$PROVDATE)]
cases_prov$WEEKDAY = as.factor(weekdays(as.Date(cases_prov$DATE)))
cases_prov$BANKHOLIDAY = as.factor(bankholiday(as.Date(cases_prov$DATE)))
cases_prov = cases_prov[complete.cases(cases_prov),]
cases_prov$PROVINCEAGE = interaction(cases_prov$PROVINCE, cases_prov$AGEGROUP)
cases_prov = cases_prov[cases_prov$DATE!=max(cases_prov$DATE),]
cases_prov = cases_prov[cases_prov$DATE!=max(cases_prov$DATE),]
# cases_prov = cases_prov[cases_prov$DATE>=as.Date("2020-08-01"),]
cases_prov$TESTSPERCASES = (cases_prov$TESTS_ALL+1)/(cases_prov$CASES+1)

# hospitalisation data (total)
hosp_tot = rawhospit[rawhospit$PROVINCE=="All"&rawhospit$REGION=="Belgium"&rawhospit$PROVINCE!="unknown",]
hosp_tot$DATE = as.Date(hosp_tot$DATE)
hosp_tot$DATE_NUM = as.numeric(hosp_tot$DATE)
hosp_tot$BANKHOLIDAY = bankholiday(hosp_tot$DATE)
hosp_tot$WEEKDAY = as.factor(weekdays(as.Date(hosp_tot$DATE)))
hosp_tot$OBS = factor(1:nrow(hosp_tot))
hosp_tot$DAYSBEFOREPRESENT = as.Date(Sys.time())-hosp_tot$DATE
hosp_tot$WEEKDAY_DAYSBEFOREPRESENT = interaction(hosp_tot$WEEKDAY,hosp_tot$DAYSBEFOREPRESENT)
# TO DO: calculate completeness based on archived data & correct for it
#hosp_tot$completeness = completeness$COMPLETENESS[match(hosp_tot$WEEKDAY_DAYSBEFOREPRESENT,completeness$WEEKDAY_DAYSBEFOREPRESENT)] 
#hosp_tot$completeness[is.na(hosp_tot$completeness)] = 1
#hosp_tot$TOTAL_IN = round(hosp_tot$TOTAL_IN/hosp_tot$completeness,0) # apply correction to most recent data
hosp_tot = hosp_tot[complete.cases(hosp_tot),]
# tail(hosp_tot,14)

# hospitalisation data per province
hosp_prov = rawhospit[rawhospit$PROVINCE!="All"&rawhospit$REGION!="Belgium"&rawhospit$PROVINCE!="unknown",]
hosp_prov$PROVINCE = as.factor(hosp_prov$PROVINCE)
# hosp_prov$AGEGROUP = as.factor(hosp_prov$AGEGROUP) # MAYBE GET THIS FROM SCIENSANO DASHBOARD?
hosp_prov$DATE = as.Date(hosp_prov$DATE)
hosp_prov$DATE_NUM = as.numeric(hosp_prov$DATE)
hosp_prov$BANKHOLIDAY = bankholiday(hosp_prov$DATE)
hosp_prov$PROVDATE = interaction(hosp_prov$PROVINCE,hosp_prov$DATE)
hosp_prov$TESTS_ALL = tests_prov$TESTS_ALL[match(hosp_prov$PROVDATE,tests_prov$PROVDATE)]
hosp_prov$WEEKDAY = as.factor(weekdays(as.Date(hosp_prov$DATE)))
hosp_prov = hosp_prov[complete.cases(hosp_prov),]
# tail(hosp_prov)
hosp_prov = hosp_prov[hosp_prov$DATE<=(max(hosp_prov$DATE)-3),]
# range(hosp_prov$DATE)


