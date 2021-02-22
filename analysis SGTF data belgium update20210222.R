# ANALYSIS OF S-GENE TARGET FAILURE (S DROPOUT) DATA FROM BELGIUM TO INFER CONTAGIOUSNESS OF NEW VARIANT OF CONCERN B.1.1.7 / 501Y.V1 ####
# PLUS INTERNATIONAL COMPARISON (USING DATA FROM THE UK, DENMARK, SWITZERLAND & THE US)
# AND PRELIMINARY INFERENCE OF GROWTH ADVANTAGE OF THE SOUTH AFRICAN VOC 501Y.V2 BASED ON SEQUENCING DATA
# T. Wenseleers & N. Hens

# Data provided by Emmanuel André (BE), COG-UK, PHE & N. Davies (UK), 
# Statens Serum Institut & Danish Covid-19 Genome Consortium (DK, https://www.covid19genomics.dk/statistics), 
# Christian Althaus, Swiss Viollier Sequencing Consortium, Institute of Medical Virology, University of Zurich, 
# Swiss National Covid-19 Science Task Force (Switzerland, https://ispmbern.github.io/covid-19/variants/, 
# https://ispmbern.github.io/covid-19/variants/data & https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv)
# and Helix, San Mateo, CA, Karthik Gangavarapu & Kristian G. Andersen (US, https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data, https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1)

# last update 22 FEBR. 2021

library(lme4)
library(splines)
library(purrr)
library(readxl)
# unloadNamespace("emmeans")
require(devtools)
# remotes::install_github("rvlenth/emmeans") # we use the development version 1.5.4-09001 as the 1.5.4 version on CRAN had a bug resulting in incorrect results for multinom models
library(emmeans)
library(effects)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(quantreg)
library(gamm4)
# install from https://github.com/tomwenseleers/export
# library(devtools)
# devtools::install_github("tomwenseleers/export")
library(export) 
library(afex)
library(dfoptim)
library(optimx)
library(mclogit)
library(lubridate)
library(zoo)
library(gridExtra)
library(sf)
library(broom)
library(nnet)


dat="2021_02_22" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today_num = as.numeric(today)
today # "2021-02-22"

# 1. ESTIMATE PROPORTION OF S DROPOUT SAMPLES THAT ARE B.1.1.7 / 501Y.V1 IN BELGIUM IN FUNCTION OF TIME BASED ON SEQUENCING DATA ####
# SEQUENCING DATA PROVIDED BY EMMANUEL ANDRÉ

datBE_b117 = read.csv(paste0(".//data//", dat, "//sequencing_Sdropouts.csv"), check.names=F) # n_b117/n_sgtf_seq = prop of S dropout samples that were B.1.1.7
datBE_b117$country = "Belgium"
datBE_b117$collection_date = as.Date(datBE_b117$collection_date)
datBE_b117$collection_date_num = as.numeric(datBE_b117$collection_date)
datBE_b117$propB117 = datBE_b117$n_b117/datBE_b117$n_sgtf_seq
datBE_b117$obs = factor(1:nrow(datBE_b117))
datBE_b117

write.csv(datBE_b117, file=".\\data\\be_latest\\sequencing_Sdropouts.csv", row.names=FALSE)

fit_seq = glmer(cbind(n_b117,n_sgtf_seq-n_b117) ~ (1|obs)+scale(collection_date_num), family=binomial(logit), data=datBE_b117)
summary(fit_seq)

# implied growth rate advantage of B.1.1.7 over other earlier strains showing S dropout:
as.data.frame(emtrends(fit_seq, ~ 1, var="collection_date_num"))[,c(2,5,6)]
#   collection_date_num.trend  asymp.LCL asymp.UCL
# 1             0.1107948 0.07844251  0.143147

# with a generation time of 4.7 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of
exp(4.7*as.data.frame(emtrends(fit_seq, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1              1.683265  1.445825  1.959699

plot(fit_seq)

# PLOT MODEL FIT
extrapolate = 90 # nr of days to extrapolate fit into the future
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                   at=list(collection_date_num=seq(as.numeric(min(datBE_b117$collection_date)),
                                                               as.numeric(max(datBE_b117$collection_date))+extrapolate)), 
                                   type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")

# prop of S dropout samples among newly diagnosed infections that are now estimated to be B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$collection_date==today,]
#    collection_date_num      prob         SE  df asymp.LCL asymp.UCL collection_date
# 62           18680 0.997214 0.002161191 Inf 0.9873354 0.9993923      2021-02-22

# prop of S dropout samples among new infections that are now estimated to be B.1.1.7 / 501Y.V1 (using 7 days for time from infection to diagnosis)
fitseq_preds[fitseq_preds$collection_date==(today+7),]
#    collection_date_num     prob          SE  df asymp.LCL asymp.UCL collection_date
# 69           18687 0.9987151 0.001143097 Inf  0.992683 0.9997756      2021-03-01

# from 13th of Jan 2021 >80% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.80,"collection_date"][1]

# from 20th of Jan 2021 >90% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.90,"collection_date"][1]

# from 11th of Feb 2021 >99% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.99,"collection_date"][1]

# on logit scale:
plot_fitseq = qplot(data=fitseq_preds, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~laboratory) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), fill=I("#b0c4de"), alpha=I(1)) +
  geom_line(aes(y=prob), colour=I("steelblue"), alpha=I(1)) +
  ylab("S dropout samples that are B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("REPRESENTATION OF 501Y.V1 AMONG S DROPOUT SAMPLES") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-12-01"),today), 
    ylim=c(0.01,0.999002), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=datBE_b117, 
             aes(x=collection_date, y=propB117, size=n_sgtf_seq), colour=I("steelblue"), alpha=I(1)) +
  scale_size_continuous("number of S dropout\nsamples sequenced", trans="sqrt", 
                        range=c(1, 6), limits=c(1,
                                                   10^(round(log10(max(datBE_b117$n_sgtf_seq)),0)+1) ), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right") + xlab("Collection date")
plot_fitseq

# saveRDS(plot_fitseq, file = paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_link scale.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_link scale.pptx"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_link scale.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_link scale.pdf"), width=8, height=6)


# same on response scale:
plot_fitseq_response = qplot(data=fitseq_preds, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~laboratory) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL), fill=I("#b0c4de"), alpha=I(1)) +
  geom_line(aes(y=100*prob), colour=I("steelblue"), alpha=I(1)) +
  ylab("S dropout samples that are B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("REPRESENTATION OF 501Y.V1 AMONG S DROPOUT SAMPLES") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-12-01"),today), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=datBE_b117, 
             aes(x=collection_date, y=100*propB117, size=n_sgtf_seq), colour=I("steelblue"), alpha=I(1)) +
  scale_size_continuous("number of S dropout\nsamples sequenced", trans="sqrt", 
                        range=c(1, 6), limits=c(1,10^(round(log10(max(datBE_b117$n_sgtf_seq)),0)+1) ), 
                        breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right") + xlab("Collection date")
plot_fitseq_response

saveRDS(plot_fitseq, file = paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_response scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig1_dataBE_propSdropoutB117_binomial GLMM_response scale.pdf"), width=8, height=6)




# 2. ANALYSIS OF Ct VALUES OF S-DROPOUT & NON-S-DROPOUT SAMPLES IN BELGIUM ####

# Read in Ct data of all valid PCRs
file_jan = paste0(".//data//", dat, "//PCR January 2021 complete.xlsx")
file_feb = paste0(".//data//", dat, "//PCR February 2020 until 9 Feb.xlsx")
file_feb2 = paste0(".//data//", dat, "//PCR February 2020 10 to 21 Feb.xlsx")
sheets = excel_sheets(file_jan)
ctdata_jan = map_df(sheets, ~ read_excel(file_jan, sheet = .x, skip = 1, 
                                       col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                       col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
ctdata_feb = map_df(sheets, ~ read_excel(file_feb, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
ctdata_feb2 = map_df(sheets, ~ read_excel(file_feb2, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
ctdata = bind_rows(ctdata_jan, ctdata_feb, ctdata_feb2)
unique(ctdata$Laboratory) 
# "UMons - Jolimont" "UZA"              "ULB"              "Saint LUC - UCL"  "UZ leuven"        "UZ Gent"          "Namur"           
# "ULG - FF 3.x" 
ctdata$Outcome[ctdata$Outcome=="Detected"] = "Positive"
ctdata$Outcome[ctdata$Outcome=="Not detected"] = "Negative"
unique(ctdata$Outcome)
ctdata$Analysis_date = as.Date(as.numeric(ctdata$Analysis_date), origin="1899-12-30")
sum(is.na(ctdata$Analysis_date)) # 0
range(ctdata$Analysis_date) # "2021-01-01" - "2021-02-21"
ctdata$collection_date = ctdata$Analysis_date-1 # collection date = analysis date-1 
sum(is.na(ctdata$collection_date)) # 0
ctdata$collection_date_num = as.numeric(ctdata$collection_date)
range(ctdata$collection_date) # "2020-12-31" "2021-02-20"
ctdata$group = interaction(ctdata$Outcome, ctdata$S_dropout)
ctdata$group = droplevels(ctdata$group)
unique(ctdata$group) # Positive.0 Positive.1 Negative.0
ctdata$Outcome = factor(ctdata$Outcome)
unique(ctdata$Outcome) # Positive Negative
ctdata$Laboratory = factor(ctdata$Laboratory)
ctdata$S_dropout = factor(ctdata$S_dropout)
head(ctdata)
str(ctdata)
nrow(ctdata) # 537655

ctdata_onlypos = ctdata[ctdata$Outcome=="Positive",] # subset with only the positive samples
ctdata_onlypos = bind_rows(ctdata_onlypos[ctdata_onlypos$S_dropout=="0",], ctdata_onlypos[ctdata_onlypos$S_dropout=="1",])

# ANALYSIS OF Ct VALUES OF S DROPOUT & NON-S DROPOUT SAMPLES

# plot & analysis of Ct values of all labs for dates from 13th of Jan onward when >80% of all S dropouts were B.1.1.7 / 501Y.V1

(fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1]) # "2021-01-13", from 13th of Jan >80% of all S dropouts are B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections
# this is the same criterion that was used for the SGTF analysis in the UK (N. Davies, pers. comm.)
subs = (ctdata_onlypos$collection_date > (fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1])) & 
       (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30) 
ctdata_onlypos_subs = ctdata_onlypos[subs,]
ctdata_onlypos_subs = ctdata_onlypos_subs[!(is.na(ctdata_onlypos_subs$S_dropout)|
                                              is.na(ctdata_onlypos_subs$N_cq)|
                                              is.na(ctdata_onlypos_subs$ORF1_cq)|
                                              (ctdata_onlypos_subs$ORF1_cq==0)),]

cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.76

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                        Laboratory correlation_Ct_N_ORF1ab
# Namur                       Namur             0.947866198
# Saint LUC - UCL   Saint LUC - UCL             0.981080423
# ULB                           ULB             0.981429368
# ULG - FF 3.x         ULG - FF 3.x             0.955847687
# UMons - Jolimont UMons - Jolimont            -0.282023995
# UZ Gent                   UZ Gent            -0.003990932
# UZ leuven               UZ leuven             0.976912170
# UZA                           UZA             0.881584309

ctcorplot_all_labs = qplot(data=ctdata_onlypos_subs, x=ORF1_cq, y=N_cq, group=Laboratory, fill=S_dropout, colour=S_dropout, size=I(3), shape=I(16)) +
  facet_wrap(~Laboratory) + 
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Ct value ORF1ab gene") + ylab("Ct value N gene") +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctcorplot_all_labs
# PS weird results for UMons & UZ Gent, and to some extent UZA, not sure of the cause
saveRDS(ctcorplot_all_labs, file = paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.pptx"), width=7, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.png"), width=7, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_correlation Ct values N ORF1 by lab_all labs.pdf"), width=7, height=6)


ctplot_rawdataN_all_labs = qplot(data=ctdata_onlypos_subs, x=collection_date, y=N_cq, group=S_dropout, 
                        colour=S_dropout, fill=S_dropout, geom="point", size=I(1), shape=I(16)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "N gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataN_all_labs
# PS UZA has suspect Ct values before the 21st of Jan & UMons & UGhent also have very different ranges for the Ct values
ctdata_onlypos_subs[ctdata_onlypos_subs$Laboratory=="UZA"&ctdata_onlypos_subs$N_cq<=15,"collection_date"][1,] # "2021-01-21"
ctplot_rawdataORF1_all_labs = qplot(data=ctdata_onlypos_subs, x=collection_date, y=ORF1_cq, group=S_dropout, 
                                 colour=S_dropout, fill=S_dropout, geom="point", size=I(1), shape=I(16)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.2), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "ORF1ab gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataORF1_all_labs

ctplots_rawdata_all_labs = ggarrange(ctplot_rawdataN_all_labs+xlab("")+theme(axis.text.x = element_blank()), 
                                     ctplot_rawdataORF1_all_labs,
                                     ncol=1, common.legend=TRUE, legend="right")
ctplots_rawdata_all_labs

saveRDS(ctplots_rawdata_all_labs, file = paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.pptx"), width=7, height=9)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.png"), width=7, height=9)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_all labs.pdf"), width=7, height=9)



# plot & analysis of Ct values of all labs except Mons (low S dropout counts & weird Ct value distribution), 
# ULG (low sample size), UZ Gent (heavily involved in active surveillance, so not representative, plus weird Ct value distribution) &
# UZA (heavily involved in active surveillance, so not representative, plus weird Ct value distribution)
# for dates from 13th of Jan onward when >80% of all S dropouts were B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections

# PS we could still consider including the data from UZA since the 21st of Jan, as Ct values 
# from then onwards are in the same range as the other labs

labs_to_remove = c("UMons - Jolimont", "ULG - FF 3.x", "UZ Gent", "UZA")
sel_labs = setdiff(unique(ctdata_onlypos_subs$Laboratory), labs_to_remove) 
sel_labs # "UZ leuven"       "Saint LUC - UCL" "ULB"             "Namur" 
# we use data from these 4 labs as the data distribution was comparable for these
# they also had large sample size & were not heavily involved in active surveillance
# sel_labs = unique(ctdata_onlypos$Laboratory) # to select data from all the labs, but distribution not comparable for all
# we use the subset of timepoints (from 13th Jan 2021 onwards) where >80% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
(fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1]) # "2021-01-13", from 13th of Jan >80% of all S dropouts are B.1.1.7 / 501Y.V1
# we also just use the positive samples with relatively strong signal, (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30)
# to not include pos samples with very low viral titers (indicative of old infections etc)
# this is the same criterion that was used for the SGTF analysis in the UK (N. Davies, pers. comm.)
subs = (ctdata_onlypos$collection_date > (fitseq_preds[fitseq_preds$prob>0.8,"collection_date"][1])) &
  (ctdata_onlypos$Laboratory %in% sel_labs) & (ctdata_onlypos$N_cq<30) & (ctdata_onlypos$ORF1_cq<30)
ctdata_onlypos_subs = ctdata_onlypos[subs,]
ctdata_onlypos_subs = ctdata_onlypos_subs[!(is.na(ctdata_onlypos_subs$S_dropout)|
                                              is.na(ctdata_onlypos_subs$N_cq)|
                                              is.na(ctdata_onlypos_subs$ORF1_cq)|
                                              (ctdata_onlypos_subs$ORF1_cq==0)),]
ctdata_onlypos_subs$Laboratory = droplevels(ctdata_onlypos_subs$Laboratory)
  
# make joint dataset for integrated analysis of both genes to estimate average effect across both sets of genes
ctdata_onlypos_subs_bothgenes = rbind(data.frame(ctdata_onlypos_subs, Gene="N gene", Ct=ctdata_onlypos_subs$N_cq), 
                                      data.frame(ctdata_onlypos_subs, Gene="ORF1ab gene", Ct=ctdata_onlypos_subs$ORF1_cq))
# we define a high viral load as one where the Ct value was 1.25x lower than in the non-S dropout sample group
# which was a Ct value < 15.099 for the N gene and < 16.001 for the ORF1ab gene
thresh_N = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","N_cq"]))/1.25
thresh_N # 15.099
thresh_ORF1 = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","ORF1_cq"]))/1.25
thresh_ORF1 # 16.001
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="N gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="N gene"]<thresh_N
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"]<thresh_ORF1


# check correlation between Ct values for N & ORF1ab gene
cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.97, t=549.68, p<2E-16

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                      Laboratory correlation_Ct_N_ORF1ab
# Namur                     Namur               0.9478662
# Saint LUC - UCL Saint LUC - UCL               0.9810804
# ULB                         ULB               0.9814294
# UZ leuven             UZ leuven               0.9769122


ctcorplot_sellabs = qplot(data=ctdata_onlypos_subs, x=ORF1_cq, y=N_cq, group=Laboratory, fill=S_dropout, colour=S_dropout, size=I(3), shape=I(16)) +
  facet_wrap(~Laboratory) + 
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Ct value ORF1ab gene") + ylab("Ct value N gene") +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept=thresh_N, colour=alpha("black", 1), lwd=I(0.3), lty=I(2)) +
  geom_vline(xintercept=thresh_ORF1, colour=alpha("black", 1), lwd=I(0.3), lty=I(2))
ctcorplot_sellabs
saveRDS(ctcorplot_sellabs, file = paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.pptx"), width=7, height=4.5)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.png"), width=7, height=4.5)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig2_dataBE_correlation Ct values N ORF1 by lab_4 main labs.pdf"), width=7, height=4.5)


ctplot_rawdataN = qplot(data=ctdata_onlypos_subs, x=collection_date, y=N_cq, group=S_dropout, 
                        colour=S_dropout, fill=S_dropout, geom="point", shape=I(16), size=I(2)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "N gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataN # there is no obvious temporal patterns

ctplot_rawdataORF1 = qplot(data=ctdata_onlypos_subs, x=collection_date, y=ORF1_cq, group=S_dropout, 
                           colour=S_dropout, fill=S_dropout, geom="point", shape=I(16), size=I(2)) +
  # geom_smooth(lwd=2, method="lm", alpha=I(0.4), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  # stat_smooth(geom="line", lwd=1.2, method="lm", alpha=I(1), fullrange=TRUE, expand=c(0,0)) + # formula='y ~ s(x, bs = "cs", k=3)') +
  facet_wrap(~Laboratory) +
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Collection date") + ylab("Ct value") + labs(title = "ORF1ab gene") +
  theme(axis.text.x = element_text(angle = 0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctplot_rawdataORF1 # there is no obvious temporal patterns

ctplots_rawdata = ggarrange(ctplot_rawdataN+xlab("")+theme(axis.text.x = element_blank()), 
                            ctplot_rawdataORF1,
                            ncol=1, common.legend=TRUE, legend="right")
ctplots_rawdata # there is no obvious temporal patterns

saveRDS(ctplots_rawdata, file = paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.pptx"), width=6, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.png"), width=6, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\dataBE_Ct values_raw data_temporal.pdf"), width=6, height=8)




# quantile/median regression to compare median Ct values of both genes across S dropout & non-S dropout samples in the different labs
qr_bothgenes0 = rq(Ct ~ Gene + S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes1 = rq(Ct ~ Gene * S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes2 = rq(Ct ~ Gene + S_dropout * Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes3 = rq(Ct ~ Gene * Laboratory + S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes4 = rq(Ct ~ (Gene + Laboratory + S_dropout)^2, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes5 = rq(Ct ~ Gene * Laboratory * S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
AIC(qr_bothgenes0, k=-1) # 214755.8
AIC(qr_bothgenes1, k=-1) # 214734
AIC(qr_bothgenes2, k=-1) # 214572.4 
AIC(qr_bothgenes3, k=-1) # 214767.4
AIC(qr_bothgenes4, k=-1) # 214556.7 # fits data best based on BIC criterion (lowest value, PS: here AIC with k<0 returns BIC)
AIC(qr_bothgenes5, k=-1) # 214585.4

summary(qr_bothgenes4)
# tau: [1] 0.5
# 
# Coefficients:
#   Value     Std. Error t value   Pr(>|t|) 
# (Intercept)             18.30352   0.05320  344.02131   0.00000
# Gene1                   -0.78890   0.05221  -15.11005   0.00000
# Laboratory1             -0.50165   0.09395   -5.33940   0.00000
# Laboratory2             -0.61455   0.09444   -6.50718   0.00000
# Laboratory3             -0.64705   0.08850   -7.31106   0.00000
# S_dropout1               1.22808   0.05277   23.27324   0.00000
# Gene1:Laboratory1        0.16578   0.09228    1.79652   0.07242
# Gene1:Laboratory2        0.07128   0.09594    0.74288   0.45756
# Gene1:Laboratory3       -0.18283   0.08814   -2.07426   0.03806
# Gene1:S_dropout1         0.23072   0.05343    4.31848   0.00002
# Laboratory1:S_dropout1  -0.03735   0.09378   -0.39828   0.69043
# Laboratory2:S_dropout1  -0.15855   0.09451   -1.67757   0.09344
# Laboratory3:S_dropout1  -0.72295   0.08858   -8.16187   0.00000

qr_emmeans_bylab99 = data.frame(emmeans(qr_bothgenes4, ~ Laboratory + Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans_bylab99
qr_emmeans99 = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans99
qr_emmeans = data.frame(emmeans(qr_bothgenes4, ~ Gene + S_dropout, level=0.95)) # median Ct values + 95% CLs
qr_emmeans
#          Gene S_dropout   emmean        SE    df lower.CL upper.CL
# 1      N gene         0 18.97342 0.0896920 32327 18.79763 19.14922
# 2 ORF1ab gene         0 20.08977 0.1005441 32327 19.89270 20.28685
# 3      N gene         1 16.05582 0.1206785 32327 15.81929 16.29236
# 4 ORF1ab gene         1 18.09507 0.1098416 32327 17.87978 18.31037

# mean difference in median Ct value of 2.46, which is highly significant across both genes: p<0.0001
contrast(emmeans(qr_bothgenes4, ~ S_dropout, level=0.95), method="pairwise") 
# contrast estimate    SE    df t.ratio p.value
# 0 - 1        2.46 0.106 32327 23.273  <.0001
# mean difference in median Ct value of 2.92 of rN gene and 1.99 for ORF1ab gene
confint(contrast(emmeans(qr_bothgenes4, ~ S_dropout|Gene, level=0.95), method="pairwise"))
# Gene = N gene:
#   contrast estimate   SE    df lower.CL upper.CL
# 0 - 1        2.92 0.15 32327     2.62     3.21
# 
# Gene = ORF1ab gene:
#   contrast estimate   SE    df lower.CL upper.CL
# 0 - 1        1.99 0.15 32327     1.70     2.29
# 
# Results are averaged over the levels of: Laboratory 
# Confidence level used: 0.95 


# violin plots by gene & lab & S dropout with expected marginal means+99% CLs of best fitting median regression model
ctviolinplots_bylab = ggplot(data=ctdata_onlypos_subs_bothgenes, aes(x=factor(S_dropout), y=Ct, fill=factor(S_dropout))) +
  geom_violin(alpha=1, colour=NA, trim=TRUE, draw_quantiles=TRUE, adjust=2, scale="width") +
  geom_crossbar(data=qr_emmeans_bylab99, aes(x=factor(S_dropout), y=emmean, ymin=lower.CL, ymax=upper.CL, group=Gene)) +
  # stat_summary(fun.data=data_summary,  
  #             geom="pointrange", aes(color=factor(S_dropout))) +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  # geom_point(aes(colour=factor(S_dropout))) +
  facet_wrap(~ Gene + Laboratory,ncol=4) +
  scale_colour_manual("", values=alpha(c("steelblue","lightcoral"), 1), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("steelblue","lightcoral"), 1), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("") + ylab("Ct value") + 
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("0","1"), labels=c("S pos","S dropout")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctviolinplots_bylab

saveRDS(ctviolinplots_bylab, file = paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.pptx"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.png"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots_by lab.pdf"), width=6, height=6)

# violin plots by gene & S dropout with expected marginal means+99% CLs of best fitting median regression model
ctviolinplots = ggplot(data=ctdata_onlypos_subs_bothgenes, aes(x=factor(S_dropout), y=Ct, fill=factor(S_dropout))) +
  geom_violin(alpha=1, colour=NA, trim=TRUE, draw_quantiles=TRUE, adjust=2, scale="width") +
  geom_crossbar(data=qr_emmeans99, aes(x=factor(S_dropout), y=emmean, ymin=lower.CL, ymax=upper.CL, group=Gene, lwd=I(0.1))) +
  # stat_summary(fun.data=data_summary,  
  #             geom="pointrange", aes(color=factor(S_dropout))) +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  # geom_point(aes(colour=factor(S_dropout))) +
  facet_wrap(~Gene,ncol=4) +
  scale_colour_manual("", values=muted(c("steelblue","lightcoral"), l=55), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=muted(c("steelblue","lightcoral"), l=55), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("") + ylab("Ct value") + 
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("0","1"), labels=c("S pos","S dropout")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctviolinplots

saveRDS(ctviolinplots, file = paste0(".\\plots\\",dat,"\\Ct values_violin plots.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.pptx"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.png"), width=6, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Ct values_violin plots.pdf"), width=6, height=6)




# binomial GLMMs to test for differences in prop with high viral load (Ct values a factor of 1.25 lower than median Ct in non-S dropout samples) :

fitct_highvirload_0A = glmer(high_viral_load ~ (1|Laboratory) + Gene + S_dropout, family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_1A = glmer(high_viral_load ~ (1|Laboratory) + Gene + S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_2A = glmer(high_viral_load ~ (1|Laboratory) + Gene + S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_3A = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene + S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_4A = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene + S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_0B = glmer(high_viral_load ~ (1|Laboratory) + Gene * S_dropout, family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_1B = glmer(high_viral_load ~ (1|Laboratory) + Gene * S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_2B = glmer(high_viral_load ~ (1|Laboratory) + Gene * S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_3B = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene * S_dropout + scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)
fitct_highvirload_4B = glmer(high_viral_load ~ (collection_date_num||Laboratory) + Gene * S_dropout * scale(collection_date_num), family=binomial(logit), data=ctdata_onlypos_subs_bothgenes)

BIC(fitct_highvirload_0A, fitct_highvirload_1A, fitct_highvirload_2A, fitct_highvirload_3A, fitct_highvirload_4A,
    fitct_highvirload_0B, fitct_highvirload_1B, fitct_highvirload_2B, fitct_highvirload_3B, fitct_highvirload_4B)
# df      BIC
# fitct_highvirload_0A  4 41190.92
# fitct_highvirload_1A  5 41192.97
# fitct_highvirload_2A  6 41202.62
# fitct_highvirload_3A  6 41203.28
# fitct_highvirload_4A  7 41212.94
# fitct_highvirload_0B  5 41167.42
# fitct_highvirload_1B  6 41169.44
# fitct_highvirload_2B  9 41199.54
# fitct_highvirload_3B  7 41179.76
# fitct_highvirload_4B 10 41209.86


# fitct_highvirload_0B best model
summary(fitct_highvirload_0B) # S dropout samples more frequently have high viral load based on N gene Ct values
# Random effects:
#   Groups     Name        Variance Std.Dev.
# Laboratory (Intercept) 0.03284  0.1812  
# Number of obs: 32340, groups:  Laboratory, 4
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.59899    0.09152  -6.545 5.95e-11 ***
#   Gene1             0.12286    0.01242   9.892  < 2e-16 ***
#   S_dropout1       -0.16234    0.01249 -12.993  < 2e-16 ***
#   Gene1:S_dropout1 -0.07224    0.01242  -5.816 6.02e-09 ***

plot(allEffects(fitct_highvirload_0B))


# odds to encounter high viral load samples based on Ct values of both genes (high vir load = Ct values >1.25x lower than median in non-S dropout samples)
# 1.38x [1.32-1.45x] 95% CLs increased among S dropout samples
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout, type="response"), method="revpairwise", type="response"))
# contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.38 0.0346 Inf      1.32      1.45

# odds ratio = 1.6 [1.49-1.71] for N gene & 1.2 [1.12-1.28] for ORF1ab gene 
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout|Gene, type="response"), method="revpairwise", type="response"))
# Gene = N gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.6 0.0554 Inf      1.49      1.71
# 
# Gene = ORF1ab gene:
#   contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.2 0.0429 Inf      1.12      1.28

fitct_highvirload_emmeans = as.data.frame(emmeans(fitct_highvirload_0B, ~ S_dropout+Gene, type="response"))
fitct_highvirload_emmeans$S_dropout = factor(fitct_highvirload_emmeans$S_dropout)
fitct_highvirload_emmeans$Gene = factor(fitct_highvirload_emmeans$Gene)
fitct_highvirload_emmeans
# S_dropout        Gene      prob         SE  df asymp.LCL asymp.UCL
# 1         0      N gene 0.3294412 0.02052853 Inf 0.2905225 0.3708482
# 2         1      N gene 0.4399057 0.02338181 Inf 0.3947135 0.4861166
# 3         0 ORF1ab gene 0.3074735 0.01980499 Inf 0.2700726 0.3475874
# 4         1 ORF1ab gene 0.3471130 0.02158259 Inf 0.3061012 0.3905270

ct_highvirload_binGLMM = ggplot(data=fitct_highvirload_emmeans, 
                                               aes(x=S_dropout, y=prob*100, fill=S_dropout, group=S_dropout)) +
  facet_wrap(~ Gene) +
  geom_col(colour=NA, position=position_dodge2(width=0.8, padding=0.5)) +
  geom_linerange(aes(ymin=asymp.LCL*100, ymax=asymp.UCL*100), position=position_dodge2(width=0.8, padding=0.5)) +
  scale_fill_manual("", values=muted(c("steelblue","lightcoral"), l=55), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_y_continuous(breaks=seq(0,100,by=10), expand=c(0,0)) +
  ylab("High viral load samples (%)") + xlab("Gene") + coord_cartesian(ylim=c(0,50)) +
  theme(legend.position = "none") +
  scale_x_discrete(breaks=c("0","1"), labels=c("S pos","S dropout"), expand=c(0.3,0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ct_highvirload_binGLMM

ctplots_multipanel = ggarrange(ctviolinplots+xlab(""), 
                               ct_highvirload_binGLMM,
                                ncol=1, legend=NULL, common.legend=FALSE)
ctplots_multipanel

saveRDS(ctplots_all_rawdata, file = paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig3_dataBE_Ct values_multipanel violin plot plus high viral load.pdf"), width=8, height=6)




# 3. ESTIMATE GROWTH RATE AND TRANSMISSION ADVANTAGE OF B.1.1.7 / 501Y.V1 IN BELGIUM BASED ON S-GENE TARGET FAILURE DATA ####

# we remove ULG - FF 3.x due to low sample size & also just use positive samples with Ct values for N & ORF1ab < 30 to focus on active, recent infections
sel_labs = setdiff(unique(ctdata$Laboratory), c("ULG - FF 3.x"))  # unique(ctdata$Laboratory) 
# setdiff(unique(ctdata$Laboratory), c("UMons - Jolimont", "ULG - FF 3.x", "UZ Gent", "UZA")) 
# setdiff(unique(ctdata$Laboratory), c("UMons - Jolimont", "UZ Gent","UZA","ULG - FF 3.x")) 
sel_labs  
pos_ctbelow30 = (ctdata$Laboratory %in% sel_labs) & (((ctdata$Outcome=="Positive")&(ctdata$N_cq<30)&(ctdata$ORF1_cq<30)))
pos_ctbelow30 = pos_ctbelow30[!is.na(pos_ctbelow30)]
pos_ctabove30 = (ctdata$Laboratory %in% sel_labs) & (((ctdata$Outcome=="Positive")&(ctdata$N_cq>=30)|(ctdata$ORF1_cq>=30)))
pos_ctabove30 = pos_ctabove30[!is.na(pos_ctabove30)]
100*sum(pos_ctabove30) / (sum(pos_ctabove30)+sum(pos_ctbelow30)) # 25% of positives have Ct > 30

subs = which((ctdata$Laboratory %in% sel_labs) & (((ctdata$Outcome=="Positive")&(ctdata$N_cq<30)&(ctdata$ORF1_cq<30))|
                                            (ctdata$Outcome=="Negative")))
ctdata_subs = ctdata
ctdata_subs = ctdata_subs[subs,]
nrow(ctdata) # 537655
nrow(ctdata_subs) # 511739

ctdata_subs$group = factor(ctdata_subs$group, 
                           levels=c("Negative.0","Positive.0","Positive.1"), 
                           labels=c("negative","S_pos","S_dropout"))
ctdata_subs$Laboratory = droplevels(ctdata_subs$Laboratory)

# aggregated counts by date (sample date) and Laboratory
data_ag = as.data.frame(table(ctdata_subs$collection_date, ctdata_subs$Laboratory, ctdata_subs$group), check.names=F)
colnames(data_ag) = c("collection_date", "LABORATORY", "GROUP", "COUNT")
data_ag_wide = spread(data_ag, GROUP, COUNT)
colnames(data_ag_wide)[colnames(data_ag_wide) %in% c("negative","S_pos","S_dropout")] = c("n_neg","n_spos","n_sgtf")
data_ag_wide$n_pos = data_ag_wide$n_spos+data_ag_wide$n_sgtf
data_ag_wide$total = data_ag_wide$n_neg + data_ag_wide$n_pos
data_ag_wide$collection_date = as.Date(data_ag_wide$collection_date)
data_ag_wide$collection_date_num = as.numeric(data_ag_wide$collection_date)
# calculate prop of S dropout that is actually B.1.1.7 / 501Y.V1 estimated from binomial GLMM:
# (using expected marginal mean calculated using emmeans, taking into account random effects)
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                     at=list(collection_date_num=seq(min(data_ag_wide$collection_date_num),
                                                                     max(data_ag_wide$collection_date_num))),
                                     type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")
# prob that S dropout was B.1.1.7 / 501Y.V1
data_ag_wide$TRUEPOS = fitseq_preds$prob[match(data_ag_wide$collection_date, fitseq_preds$collection_date)] 
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate
data_ag_wide$est_n_B117 = data_ag_wide$n_sgtf * data_ag_wide$TRUEPOS 
data_ag_wide$propB117 = data_ag_wide$est_n_B117 / data_ag_wide$n_pos
data_ag_wide$obs = factor(1:nrow(data_ag_wide))
data_ag_wide = data_ag_wide[data_ag_wide$total != 0, ]
head(data_ag_wide)

write.csv(data_ag_wide, file=".\\data\\be_latest\\be_B117_by lab.csv", row.names=FALSE)


# aggregated counts by date over all Laboratories
data_ag_byday = as.data.frame(table(ctdata_subs$collection_date, ctdata_subs$group), check.names=F)
colnames(data_ag_byday) = c("collection_date", "GROUP", "COUNT")
data_ag_byday_wide = spread(data_ag_byday, GROUP, COUNT)
colnames(data_ag_byday_wide)[colnames(data_ag_byday_wide) %in% c("negative","S_pos","S_dropout")] = c("n_neg","n_spos","n_sgtf")
data_ag_byday_wide$n_pos = data_ag_byday_wide$n_spos+data_ag_byday_wide$n_sgtf
data_ag_byday_wide$total = data_ag_byday_wide$n_neg + data_ag_byday_wide$n_pos
data_ag_byday_wide$collection_date = as.Date(data_ag_byday_wide$collection_date)
data_ag_byday_wide$collection_date_num = as.numeric(data_ag_byday_wide$collection_date)
# calculate prop of S dropout that is actually B.1.1.7 / 501Y.V1 estimated from binomial GLMM:
# (using expected marginal mean calculated using emmeans, taking into account random effects)
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                     at=list(collection_date_num=seq(min(data_ag_byday_wide$collection_date_num),
                                                                     max(data_ag_byday_wide$collection_date_num))),
                                     type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")
# prob that S dropout was B.1.1.7 / 501Y.V1
data_ag_byday_wide$TRUEPOS = fitseq_preds$prob[match(data_ag_byday_wide$collection_date, fitseq_preds$collection_date)] 
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate
data_ag_byday_wide$est_n_B117 = data_ag_byday_wide$n_sgtf * data_ag_byday_wide$TRUEPOS
data_ag_byday_wide$propB117 = data_ag_byday_wide$est_n_B117 / data_ag_byday_wide$n_pos
data_ag_byday_wide$obs = factor(1:nrow(data_ag_byday_wide))
data_ag_byday_wide = data_ag_byday_wide[data_ag_byday_wide$total != 0, ]
head(data_ag_byday_wide)

write.csv(data_ag_byday_wide, file=".\\data\\be_latest\\be_B117_total.csv", row.names=FALSE)



# 3.1 ESTIMATE GROWTH RATE & TRANSMISSION ADVANTAGE OF B.1.1.7 USING BINOMIAL GLMM (LOGISTIC FIT) ####

# fit common-slope and separate-slopes binomial GLM
set_sum_contrasts()
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit1)
fit1 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)+LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings)  # common slope model, with lab coded as fixed factor
fit2 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)*LABORATORY, family=binomial(logit), 
             data=data_ag_wide, subset=data_ag_wide$n_pos>0, control=glmersettings) # separate slopes model, with lab coded as fixed factor
BIC(fit1,fit2) 
#  df      BIC
# fit1  9 2389.264
# fit2 15 2407.336


# common-slope model fit1 fits best, i.e. rate at which B.1.1.7 is displacing other strains constant across regions/labs

summary(fit1)
# Random effects:
#   Groups Name        Variance Std.Dev.
# obs    (Intercept) 0.3395   0.5827  
# Number of obs: 363, groups:  obs, 363
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                -1.45051    0.03915 -37.048  < 2e-16 ***
#   scale(collection_date_num)  1.12805    0.04206  26.818  < 2e-16 ***
#   LABORATORY1                -0.22667    0.09433  -2.403 0.016264 *  
#   LABORATORY2                 0.29534    0.08659   3.411 0.000648 ***
#   LABORATORY3                 0.44740    0.08661   5.166 2.39e-07 ***
#   LABORATORY4                -0.94862    0.09829  -9.651  < 2e-16 ***
#   LABORATORY5                 0.32641    0.09237   3.534 0.000410 ***
#   LABORATORY6                -0.20639    0.08904  -2.318 0.020456 *   

# growth rate advantage (differences in growth rate between B.1.1.7 and old strains):
# results common-slope model
fit1_emtrends = as.data.frame(emtrends(fit1, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=today_num), mode="link", adjust="Tukey")$emtrends)
fit1_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL  asymp.UCL
# 1                0.07525432 0.06975435 0.08075429

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*4.7) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.424321  1.387974   1.46162

# with a generation time of 5.5 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*5.5) 
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.512704   1.46763  1.559162


# tests for differences in date of introduction
# UCL, ULB, Ghent & UZA earlier than avg, Namur, Mons & UZ Leuven later than avg
emmeans(fit1,eff~LABORATORY)$contrasts 
# contrast                  estimate     SE  df z.ratio p.value
# Namur effect                -0.227 0.0943 Inf -2.403  0.0190 
# (Saint LUC - UCL) effect     0.295 0.0866 Inf  3.411  0.0009 
# ULB effect                   0.447 0.0866 Inf  5.166  <.0001 
# (UMons - Jolimont) effect   -0.949 0.0983 Inf -9.651  <.0001 
# UZ Gent effect               0.326 0.0924 Inf  3.534  0.0007 
# UZ leuven effect            -0.206 0.0890 Inf -2.318  0.0205 
# UZA effect                   0.313 0.0862 Inf  3.625  0.0007 
# # 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: fdr method for 8 tests

# results of growth rate advantage of separate-slopes model fit2 by lab/region:                         
fit2_emtrends = emtrends(fit2, revpairwise ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$emtrends
fit2_emtrends
# LABORATORY       collection_date_num.trend      SE  df asymp.LCL asymp.UCL
# Namur                               0.0834 0.00787 Inf    0.0679    0.0988
# Saint LUC - UCL                     0.0636 0.00647 Inf    0.0509    0.0763
# ULB                                 0.0742 0.00654 Inf    0.0614    0.0870
# UMons - Jolimont                    0.1058 0.00925 Inf    0.0877    0.1239
# UZ Gent                             0.0673 0.00762 Inf    0.0523    0.0822
# UZ leuven                           0.0755 0.00696 Inf    0.0619    0.0892
# UZA                                 0.0693 0.00655 Inf    0.0564    0.0821

# for Mons we observe a sign above-average rate of spread
fit2_contrasts = emtrends(fit2, eff ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$contrasts
fit2_contrasts
# contrast                  estimate      SE  df z.ratio p.value
# Namur effect               0.00637 0.00720 Inf  0.884  0.9634 
# (Saint LUC - UCL) effect  -0.01340 0.00614 Inf -2.181  0.1873 
# ULB effect                -0.00283 0.00619 Inf -0.457  0.9993 
# (UMons - Jolimont) effect  0.02881 0.00829 Inf  3.476  0.0036 
# UZ Gent effect            -0.00972 0.00697 Inf -1.393  0.7134 
# UZ leuven effect          -0.00149 0.00648 Inf -0.230  1.0000 
# UZA effect                -0.00775 0.00618 Inf -1.253  0.8082 
# 
# P value adjustment: sidak method for 7 tests 



# PLOT MODEL FIT

# for best fitting common slope model fit1
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-05-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit1_preds = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                                         # by="LABORATORY", 
                                         at=list(collection_date_num=seq(date.from,
                                                                         date.to)), 
                                         type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_preds$collection_date = as.Date(fit1_preds$collection_date_num, origin="1970-01-01")


total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1))$sdcor, function (x) x^2))) 
fit1_preds_bylab = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                                   by="LABORATORY", 
                                   at=list(collection_date_num=seq(date.from,
                                                               date.to)), 
                                    type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_preds_bylab$collection_date = as.Date(fit1_preds_bylab$collection_date_num, origin="1970-01-01")
# order labs by estimated date of introduction (intercepts)
dfemmeanslabs = as.data.frame(emmeans(fit1,~LABORATORY))
levels_BE = as.character(dfemmeanslabs$LABORATORY[order(dfemmeanslabs$emmean,decreasing=T)])
fit1_preds_bylab$LABORATORY = factor(fit1_preds_bylab$LABORATORY, 
                                     levels=levels_BE)


# estimated share of B.1.1.7 among currently diagnosed infections based on fit1
fit1_preds[fit1_preds$collection_date==today,]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 27              18680 0.6372716 0.01608896 Inf 0.6052543 0.6682487      2021-02-22

# estimated share of B.1.1.7 among new infections (assuming time between infection & diagnosis of 7 days)
fit1_preds[fit1_preds$collection_date==(today+7),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 34               18687 0.7417011 0.01675286 Inf 0.7076555 0.7732591      2021-03-01


sum(tail(data_ag_byday_wide$est_n_B117, 14))/sum(tail(data_ag_byday_wide$n_pos,14)) 
# 47% of the samples of last 2 weeks in the dataset were estimated to be by British variant
# PS: with data 31/1 this was 15.4%  
# note: this is not the same as the estimated prop of the new infections or new diagnoses today that are of the British
# variant, which are much higher, see above)


# implied Re of wild type and B.1.1.7 given this predicted share of B.1.1.7 among all infections ####
# under a particular fitted transmission advantage
# based on the fact that the overall Re is a weighted average of the Re of the individual variants
# functions to calculate Re of wild type and of B.1.1.7 based on overall Re value and prop of positives that is B.1.1.7 propB117
# and transmission advantage of B.1.1.7 M
Re_wild_type = function (Re, propB117, M=1.5) {
  Re / (1-propB117+M*propB117)
}
Re_B117 = function (Re, propB117, M=1.5) {
  M*Re / (1-propB117+M*propB117)
}

Re_cases = read.csv(".//Re_fits//Re_cases.csv") 
# Re values calculated from instantaneous growth rate in nr of new cases 
# with instant growth rate = derivative/emtrends in function of sample date of GAM fit on new cases 
# gam(cbind(NEWCASES, totpop-NEWCASES) ~ s(DATE_NUM, bs="cs", k=23, fx=F) + WEEKDAY + s(log(TESTS_ALL), bs="cs", k=5, fx=F), family=binomial(cloglog), data=cases_tot) 
# and with Re calculated using R.from.r with gamma_mean=4.7, gamma_sd=2.9
# Re here is at time of diagnosis, we use a shift of 7 days to recalculate to Re at day of infection
Re_cases$DATE = as.Date(Re_cases$DATE_NUM, origin="1970-01-01")
Re_cases$collection_date_num = Re_cases$DATE_NUM
Re_cases$propB117 = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                      # by="LABORATORY", 
                      at=list(collection_date_num=seq(min(Re_cases$collection_date_num),
                                                      max(Re_cases$collection_date_num))), 
                      type="response"), bias.adjust = TRUE, sigma = total.SD)$prob

head(Re_cases)
Re_cases$Re_WT = Re_wild_type(Re=Re_cases$Re, propB117=Re_cases$propB117, M=1.5)
Re_cases$Re_WT_LOWER = Re_wild_type(Re=Re_cases$Re_LOWER, propB117=Re_cases$propB117, M=1.5)
Re_cases$Re_WT_UPPER = Re_wild_type(Re=Re_cases$Re_UPPER, propB117=Re_cases$propB117, M=1.5)
Re_cases$Re_B117 = Re_B117(Re=Re_cases$Re, propB117=Re_cases$propB117, M=1.5)
Re_cases$Re_B117_LOWER = Re_B117(Re=Re_cases$Re_LOWER, propB117=Re_cases$propB117, M=1.5)
Re_cases$Re_B117_UPPER = Re_B117(Re=Re_cases$Re_UPPER, propB117=Re_cases$propB117, M=1.5)

# d = as.Date(max(Re_cases$DATE))
# tag = paste("@TWenseleers\n",dat)

qplot(data=Re_cases, x=DATE, y=Re, ymin=Re_LOWER, ymax=Re_UPPER, geom="ribbon", alpha=I(0.5), fill=I("grey")) +
  geom_line() + theme_hc() + xlab("") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous(limits=c(1/3,3), trans="log2") +
  geom_hline(yintercept=1, colour=I("red")) +
  ggtitle("Re OF B.1.1.7 (red) AND WILD TYPE (blue) IN BELGIUM\nBASED ON NEW CONFIRMED CASES AND\nB.1.1.7 TRANSMISSION ADVANTAGE OF 50%") +
  # labs(tag = tag) +
  # theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) +
  theme(plot.tag.position = "bottomright",
        plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  geom_ribbon(aes(y=Re_WT, ymin=Re_WT_LOWER, ymax=Re_WT_UPPER), fill=I("blue"), alpha=I(0.5)) +
  geom_line(aes(y=Re_WT), fill=I("blue"), alpha=I(0.8)) +
  geom_ribbon(data=Re_cases[Re_cases$DATE>=as.Date("2021-01-01"),], 
              aes(y=Re_B117, ymin=Re_B117_LOWER, ymax=Re_B117_UPPER), fill=I("red"), alpha=I(0.5)) +
  geom_line(data=Re_cases[Re_cases$DATE>=as.Date("2021-01-01"),],
            aes(y=Re_B117), fill=I("blue"), alpha=I(0.8)) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"),max(Re_cases$DATE)),
                  ylim=c(0.5,1.5))
ggsave(file=paste0(".//plots//Re",dat,"_cases_Re_B117_Re_wildtype.png"), width=7, height=5)
Re_cases[Re_cases$DATE==max(Re_cases$DATE),]




# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50%, 75% 90% by B.1.1.7 :
fit1_preds$collection_date[fit1_preds[,"prob"]>=0.5][1]-7 # >50% by 7th of February [6th Feb - 9 Feb] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.5][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.5][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.75][1]-7 # >75% by 23d of February [21 Feb - 26 Feb] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.75][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.75][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.9][1]-7 # >90% by 10th of March [7 March - 14 March] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.9][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.9][1]-7



# PLOT MODEL FIT common-slope model fit1

# plot for the whole of Belgium on response scale:
plot_fit1_response = qplot(data=fit1_preds, x=collection_date, y=100*prob, geom="blank") +
  # facet_wrap(~LABORATORY) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, 
                  # fill=LABORATORY
  ), 
  fill=I("#b0c4de"), 
  alpha=I(1)) +
  geom_line(aes(y=100*prob, 
                # colour=LABORATORY
  ), 
  colour=I("steelblue"), 
  alpha=I(1)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC 202012/01 BY NHS REGION") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(data_ag_byday_wide$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_byday_wide, 
             aes(x=collection_date, y=100*propB117, size=n_pos,
                 # colour=LABORATORY
             ), 
             colour=I("steelblue"), 
             alpha=I(1)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(10,10^round(log10(max(data_ag_byday_wide$n_pos)+1),0)), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date")
plot_fit1_response


saveRDS(plot_fit1_response, file = paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_B117_Belgium_response scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_B117_Belgium_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_B117_Belgium_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig4_fit1_binomGLMM_B117_Belgium_response scale.pdf"), width=8, height=6)


# plot per lab on logit scale:
plot_fit1 = qplot(data=fit1_preds_bylab, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~LABORATORY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=LABORATORY
                  ), 
              # fill=I("steelblue"), 
              alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=LABORATORY
                ), 
            # colour=I("steelblue"), 
            alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(data_ag_byday_wide$collection_date), as.Date("2021-04-01")-1), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.01,0.99), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_wide,  
             aes(x=collection_date, y=propB117, size=n_pos,
                 colour=LABORATORY
                 ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="log10", 
                        range=c(1, 4), limits=c(10,10^round(log10(max(data_ag_wide$n_pos)+1),0)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90))
plot_fit1

saveRDS(plot_fit1, file = paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_B117_Belgium by lab.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_B117_Belgium by lab.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_B117_Belgium by lab.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig5_fit1_binomGLMM_B117_Belgium by lab.pdf"), width=8, height=6)




# 3.2 ESTIMATE GROWTH RATE & R VALUE OF B.1.1.7 & WILD TYPE STRAINS SEPARATELY USING MULTINOMIAL MODEL ####

# PS this analysis I couldn't update since last batch of data from 10/2-21/2 do not include PCR negative tests

# Function to recalculate Malthusian growth rate to effective reproduction number Rt
# (assuming generation time is gamma distributed)
# Ref: Park et al. 2020, https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0144
R.from.r <- function(r, gamma_mean=4.7, gamma_sd=2.9) {
  k = (gamma_sd / gamma_mean)^2
  R = (1 + k * r * gamma_mean)^(1 / k)
  return(R)
}

data_ag_wide2 = data_ag_wide 
data_ag_wide2$n_b117 = round(data_ag_wide2$est_n_B117, 0)
data_ag_wide2$n_spos = data_ag_wide2$n_pos-data_ag_wide2$n_b117
(data_ag_wide2$n_neg+data_ag_wide2$n_spos+data_ag_wide2$n_b117)==data_ag_wide2$total # check
data_ag_wide2 = data_ag_wide2[,c("collection_date","LABORATORY","n_neg","n_spos","n_b117")]
head(data_ag_wide2)
data_ag_long = gather(data_ag_wide2, outcome, count, n_neg:n_b117, factor_key=TRUE)
data_ag_long$outcome = factor(data_ag_long$outcome, levels=c("n_neg","n_spos","n_b117"))
data_ag_long$collection_date_num = as.numeric(data_ag_long$collection_date)

test_outcomes = ggplot(data=data_ag_long, 
      aes(x=collection_date, 
      y=count, fill=outcome, group=outcome)) +
  facet_wrap(~LABORATORY) +
  geom_area(aes(fill=outcome), position = position_fill(reverse = TRUE)) +
  theme_hc() +
  scale_fill_manual("test outcome", values=c("darkgrey","steelblue","lightcoral"), labels=c("negative","S positive","S dropout")) +
  ylab("Share") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  # ggtitle("Test outcomes") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right")
test_outcomes

saveRDS(test_outcomes, file = paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_for_multinomial spline fit.pdf"), width=7, height=5)

test_outcomes_pos = ggplot(data=data_ag_long[data_ag_long$outcome!="n_neg",], 
                       aes(x=collection_date, 
                           y=count, fill=outcome, group=outcome)) +
  facet_wrap(~LABORATORY) +
  geom_area(aes(fill=outcome), position = position_fill(reverse = TRUE)) +
  theme_hc() +
  scale_fill_manual("test outcome", values=c("steelblue","lightcoral"), labels=c("S positive","S dropout")) +
  ylab("Share") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  # ggtitle("Test outcomes") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right")
test_outcomes_pos

saveRDS(test_outcomes, file = paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\test_outcomes_share_Sdropout_positives only.pdf"), width=7, height=5)


# multinomial spline fit on test outcome data (negative / positive wild type / positive B.1.1.7
# to be able to estimate growth rate and Rt of B.1.1.7 and wild type separately

set.seed(1)
# we use data from the 14th of Jan onwards, as data has been approx randomly sampled from then on
sel_labs = unique(data_ag_long$LABORATORY)
date.from = as.Date("2021-01-14")
data_ag_long_subs = data_ag_long[(data_ag_long$LABORATORY %in% sel_labs)&(data_ag_long$collection_date>=date.from),]
data_ag_long_subs$LABORATORY = droplevels(data_ag_long_subs$LABORATORY)

mfit0 = nnet::multinom(outcome ~ scale(collection_date_num) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mfit1 = nnet::multinom(outcome ~ scale(collection_date_num) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000)
mfit2 = nnet::multinom(outcome ~ ns(collection_date_num, df=2) + LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
mfit3 = nnet::multinom(outcome ~ ns(collection_date_num, df=3) * LABORATORY, weights=count, data=data_ag_long_subs, maxit=1000) 
BIC(mfit0, mfit1, mfit2, mfit3) # mfit3 fits best, df splines tuned based on BIC
#       df      BIC
# mfit0 16 150913.8
# mfit1 28 150894.7
# mfit2 18 149513.9
# mfit3 56 149434.8
summary(mfit3)

plot(Effect("collection_date_num",mfit3), style="stacked")
plot(Effect("collection_date_num",mfit3), confint=list(style="bands"), rug=FALSE)


# mblogit mulinomial fit with possible orrections for overdispersion
# should be better, but I still have to update my emmeans analyses below to use these models instead
# and right now emmeans only has experimental support for mblogit models
# PS lab can also be coded as a random factor here
mblogitfit0 = mblogit(outcome ~ scale(collection_date_num) + LABORATORY, weights=count, data=data_ag_long_subs, 
                      dispersion=FALSE, control=mclogit.control(maxit = 100, trace=TRUE))
mblogitfit1 = mblogit(outcome ~ scale(collection_date_num) * LABORATORY, weights=count, data=data_ag_long_subs, 
                      dispersion=FALSE)
mblogitfit2 = mblogit(outcome ~ ns(collection_date_num, df=2) * LABORATORY, weights=count, data=data_ag_long_subs, 
                      dispersion=FALSE)
mblogitfit3 = mblogit(outcome ~ ns(collection_date_num, df=3) * LABORATORY, weights=count, data=data_ag_long_subs, 
                      dispersion=FALSE)
BIC(mblogitfit0, mblogitfit1, mblogitfit2, mblogitfit3)
#             df      BIC
# mblogitfit0 16 150913.8
# mblogitfit1 28 150870.9
# mblogitfit2 42 149479.6
# mblogitfit3 56 149434.8
summary(mblogitfit3)
dispersion(mblogitfit0, method="Afroz") # 36.85, i.e. >>1, there is overdispersion in the data
dispersion(mblogitfit3, method="Afroz") # 31.68, i.e. >>1, there is overdispersion in the data


# coefs & conf intervals of mfit0 & mblogitfit0 are almost identical as should be the case
mfit0_coefs = data.frame(tidy(mfit0, conf.int=TRUE))[,c("term","estimate","conf.low","conf.high")]
rownames(mfit0_coefs) = paste0(rep(mfit0$lab[-1],each=nrow(mfit0_coefs)/length(mfit0$lab[-1])),"~",mfit0_coefs$term) 
mfit0_coefs[,c("term")] = NULL
mfit0_coefs = mfit0_coefs[grepl("date",rownames(mfit0_coefs)),]
mfit0_coefs = mfit0_coefs / attr(scale(data_ag_long_subs$collection_date_num),"scaled:scale")
mfit0_coefs
#                                      estimate    conf.low   conf.high
# n_spos~scale(collection_date_num) -0.03388852 -0.03642875 -0.03134829
# n_b117~scale(collection_date_num)  0.01453788  0.01036944  0.01870632


mblogitfit_coefs = data.frame(tidy(mblogitfit0, conf.int=TRUE))[,c("term","estimate","conf.low","conf.high")]
rownames(mblogitfit_coefs) = mblogitfit_coefs$term
mblogitfit_coefs[,c("term")] = NULL
mblogitfit_coefs = mblogitfit_coefs[grepl("date",rownames(mblogitfit_coefs)),]
mblogitfit_coefs = mblogitfit_coefs / attr(scale(data_ag_long_subs$collection_date_num),"scaled:scale")
mblogitfit_coefs
#                                    estimate    conf.low   conf.high
# n_spos~scale(collection_date_num) -0.03388029 -0.03642333 -0.03133724
# n_b117~scale(collection_date_num)  0.01454377  0.01037025  0.01871729

# this is another way to calculate these multinomial slopes:
confint(emtrends(mfit0, trt.vs.ctrl~outcome|1, var="collection_date_num",  mode="latent", adjust="none"))$contrasts
# contrast       estimate      SE df lower.CL upper.CL
# n_spos - n_neg  -0.0339 0.00130 16  -0.0366  -0.0311
# n_b117 - n_neg   0.0145 0.00213 16   0.0100   0.0190


# results for best model mfit3:
# average growth rates of S-positive/wild type & S dropout/B.1.1.7 cases over period from Jan 14 till now
emtrends(mfit3, ~outcome|1, var="collection_date_num",   mode="latent")
# PS emtrends(mblogitfit3, ~outcome|1, var="collection_date_num",  mode="latent") gives Error in qr.default(t(const)) : NA/NaN/Inf in foreign function call (arg 1)
# emmeans_1.5.4 output:
# outcome collection_date_num.trend      SE df lower.CL upper.CL
# n_neg                     0.00392 0.00322 56 -0.00253   0.0104
# n_spos                   -0.02795 0.00396 56 -0.03589  -0.0200
# n_b117                    0.02403 0.00579 56  0.01244   0.0356

R.from.r(-0.02795) # Rt of S-positive/wild type = 0.87
R.from.r(0.02403) # Rt of S dropout/B.1.1.7 variant = 1.12
R.from.r(0.02403)/R.from.r(-0.02795) # Rt of B.1.1.7 = 1.28x times higher than of wild type


# implied growth rate & transmission advantage of B.1.1.7 vs wild type + 95% CLs
delta_r = data.frame(confint(contrast(emtrends(mfit3, ~outcome|1, var="collection_date_num",  
                                               at=list(outcome=c("n_spos","n_b117")), mode="latent"), method="revpairwise")))[,c(2,5,6)]
delta_r # growth advantage
#    estimate lower.CL  upper.CL
# 1 0.05197395 0.03317961 0.07076829
exp(delta_r*4.7) # transmission advantage
#    estimate lower.CL upper.CL
# 1 1.276699 1.168761 1.394605


# average growth rates of S-positive/wild type & S dropout/B.1.1.7 cases evaluated today
emtrends(mfit3, ~outcome|1, var="collection_date_num",  at=list(collection_date_num=today_num), mode="latent")
# emmeans_1.5.3 output:
# outcome collection_date_num.trend      SE df lower.CL upper.CL
# n_neg                      0.1027 0.00562 56   0.0914   0.1140
# n_spos                    -0.0725 0.00741 56  -0.0873  -0.0576
# n_b117                    -0.0302 0.00982 56  -0.0499  -0.0106
R.from.r(-0.0725) # Rt S pos / wild type = 0.69
R.from.r(-0.0873) # Rt S pos / wild type LCL = 0.64
R.from.r(-0.0576) # Rt S pos / wild type UCL = 0.75

R.from.r(-0.0302) # Rt of S dropout = 0.86
R.from.r(-0.0499) # Rt of S dropout = 0.78
R.from.r(-0.0106) # Rt of S dropout = 0.95


R.from.r(-0.0302)/R.from.r(-0.0725) # Rt of B.1.1.7 = 1.24x times higher than of wild type
# PS recent shift from active surveillance for B.1.1.7 in beginning of January to random sampling
# right now might introduce a downward bias here though

# implied growth rate & transmission advantage
delta_r = data.frame(confint(contrast(emtrends(mfit3, ~outcome|1, var="collection_date_num",  
                                               at=list(outcome=c("n_spos","n_b117"),
                                                       collection_date_num=today_num), mode="latent"), method="revpairwise")))[,c(2,5,6)]
delta_r # growth advantage
#       estimate lower.CL  upper.CL
# 1 0.04222412 0.009238266 0.07520997
exp(delta_r*4.7) # transmission advantage
#    estimate lower.CL upper.CL
# 1 1.219515 1.044376 1.424024


# growth rates and Re values of the B.1.1.7 variant and the wild type calculated over time
extrapolate = 0
r_and_Re_B117_wildtype = data.frame(emtrends(mfit3, ~outcome|1, var="collection_date_num", by=c("collection_date_num"), # by=c("collection_date_num","LABORATORY"),  
                                  at=list(collection_date_num=seq(min(data_ag_long_subs$collection_date_num),
                                                                  today_num+extrapolate),
                                          outcome=c("n_spos","n_b117")), mode="latent"))
r_and_Re_B117_wildtype$collection_date = as.Date(r_and_Re_B117_wildtype$collection_date_num, origin="1970-01-01")
colnames(r_and_Re_B117_wildtype)[colnames(r_and_Re_B117_wildtype) %in% c("collection_date_num.trend","lower.CL","upper.CL")] = c("r","r.LCL","r.UCL")
r_and_Re_B117_wildtype$Re = R.from.r(r_and_Re_B117_wildtype$r)
r_and_Re_B117_wildtype$Re.LCL = R.from.r(r_and_Re_B117_wildtype$r.LCL)
r_and_Re_B117_wildtype$Re.UCL = R.from.r(r_and_Re_B117_wildtype$r.UCL)

plot_Re_B117_WT = qplot(data=r_and_Re_B117_wildtype, x=collection_date, y=Re, ymin=Re.LCL, ymax=Re.UCL, geom="ribbon", alpha=I(0.5), 
      fill=outcome , colour=NULL, group=outcome ) +
  # facet_wrap(~LABORATORY) +
  geom_line(aes(colour=outcome )) + theme_hc() + xlab("") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous(trans="log2", breaks=c(1/seq(3,1),seq(1,4)),
                     labels=round(c(1/seq(3,1),seq(1,4)),2)) +
  coord_cartesian(xlim=c(min(data_ag_long_subs$collection_date),
                         max(r_and_Re_B117_wildtype$collection_date)), 
                  ylim=c(1/2,2), 
                  expand=c(0,0)) +
  geom_hline(yintercept=1, colour=alpha(I("black"),0.2)) +
  # theme(legend.position = "none") +
  ggtitle("Effective reproduction nr. Re of B.1.1.7 and wild type") +
  guides(colour=FALSE) +
  scale_colour_manual("", values=c("steelblue","lightcoral")) +
  scale_fill_manual("", values=c("steelblue","lightcoral"), labels=c("S positive (wild type)","S dropout (B.1.1.7)")) +
  theme(legend.position = "bottom")
  # labs(tag = tag) +
  # theme(plot.tag.position = "bottomright",
  #      plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_Re_B117_WT

# implied growth rate advantage + 95% CLs over time
extrapolate = 0
growthadvantage = data.frame(confint(contrast(emtrends(mfit3, ~outcome|1, var="collection_date_num", by="collection_date_num", 
                                               at=list(outcome=c("n_spos","n_b117"),
                                                       collection_date_num=seq(min(data_ag_long$collection_date_num),
                                                                               today_num+extrapolate)), 
                                               mode="latent"), method="revpairwise")))
colnames(growthadvantage)[colnames(growthadvantage) %in% c("estimate","lower.CL","upper.CL")] = c("delta_r","delta_r.LCL","delta_r.UCL")
growthadvantage$transmadv = (exp(4.7*growthadvantage$delta_r)-1)*100
growthadvantage$transmadv.LCL = (exp(4.7*growthadvantage$delta_r.LCL)-1)*100
growthadvantage$transmadv.UCL = (exp(4.7*growthadvantage$delta_r.UCL)-1)*100
growthadvantage$collection_date = as.Date(growthadvantage$collection_date_num, origin="1970-01-01")

plot_growthadvB117 = qplot(data=growthadvantage, x=collection_date, y=transmadv, 
                            ymin=transmadv.LCL, ymax=transmadv.UCL, geom="ribbon", alpha=I(0.5), 
      fill=I("lightcoral") , colour=NULL ) +
  geom_line(aes(colour=I("lightcoral"))) + theme_hc() + xlab("") + ylab("Transmission advantage (%)") +
  geom_hline(yintercept=0, colour=alpha(I("black"),0.2)) +
  ggtitle("Transmission advantage of B.1.1.7 over wild type") +
  coord_cartesian(xlim=c(min(data_ag_long_subs$collection_date),
                         max(r_and_Re_B117_wildtype$collection_date)), 
                  expand=c(0,0))
# labs(tag = tag) +
# theme(plot.tag.position = "bottomright",
#      plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_growthadvB117

multipanel_Re_growthadv_multinom = ggarrange(plot_Re_B117_WT, plot_growthadvB117, ncol=1)
multipanel_Re_growthadv_multinom

saveRDS(multipanel_Re_growthadv_multinom, file = paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.pptx"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.png"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\multinomial_Re_growthadv_B117_WT.pdf"), width=7, height=8)



# 3.3 PRELIMINARY ASSESSMENT OF GROWTH RATE ADVANTAGES OF VOC 501Y.V1,VOC 501Y.V2&VOC 501Y.V3 IN BELGIUM BASED ON SEQUENCING DATA ####

# data taken from latest Sciensano weekly report "COVID 19 WEKELIJKS EPIDEMIOLOGISCH BULLETIN (19 FEBRUARI 2021), p. 24"
# https://covid-19.sciensano.be/nl/covid-19-epidemiologische-situatie

# unloadNamespace("emmeans") # install latest development version of emmeans, since CRAN version 1.5.4 had bug in multinom models
library(devtools)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(splines)
library(broom)
library(nnet)
library(emmeans)
library(ggpubr)
library(export)
library(mclogit)

dat="2021_02_22" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
filedate = as.Date(gsub("_","-",dat)) # file date
filedate_num = as.numeric(filedate)
today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today_num = as.numeric(today)

be_seqdata = read.csv(".\\data\\be_latest\\sequencing_501YV1_501YV2_501YV3.csv")
# data is split up in baseline surveillance (randomly sampled) and active surveillance (from travellers, known outbreaks &
# S dropout sequencing), below I will use the randomly sampled baseline surveillance part
be_seqdata$collection_date = as.Date(be_seqdata$collection_date)
be_seqdata$baselinesurv_n_wild_type = be_seqdata$baselinesurv_total_sequenced-be_seqdata$baselinesurv_n_501Y.V1-be_seqdata$baselinesurv_n_501Y.V2-be_seqdata$baselinesurv_n_501Y.V3
be_seqdata$baselinesurv_prop501YV1 = be_seqdata$baselinesurv_n_501Y.V1 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_prop501YV2 = be_seqdata$baselinesurv_n_501Y.V2 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_prop501YV3 = be_seqdata$baselinesurv_n_501Y.V3 / be_seqdata$baselinesurv_total_sequenced
be_seqdata$baselinesurv_n_501Y.V1plusV2plusV3 = be_seqdata$baselinesurv_n_501Y.V1+be_seqdata$baselinesurv_n_501Y.V2+be_seqdata$baselinesurv_n_501Y.V3
be_seqdata$baselinesurv_n_501Y.propV1V2V3 = be_seqdata$baselinesurv_n_501Y.V1plusV2plusV3 / be_seqdata$baselinesurv_total_sequenced

be_seqdata$activesurv_n_wild_type = be_seqdata$activesurv_total_sequenced-be_seqdata$activesurv_n_501Y.V1-be_seqdata$activesurv_n_501Y.V2-be_seqdata$activesurv_n_501Y.V3
be_seqdata$activesurv_prop501YV1 = be_seqdata$activesurv_n_501Y.V1 / be_seqdata$activesurv_total_sequenced
be_seqdata$activesurv_prop501YV2 = be_seqdata$activesurv_n_501Y.V2 / be_seqdata$activesurv_total_sequenced
be_seqdata$activesurv_prop501YV3 = be_seqdata$activesurv_n_501Y.V3 / be_seqdata$activesurv_total_sequenced
be_seqdata$activesurv_n_501Y.V1plusV2plusV3 = be_seqdata$activesurv_n_501Y.V1+be_seqdata$activesurv_n_501Y.V2+be_seqdata$activesurv_n_501Y.V3
be_seqdata$activesurv_n_501Y.propV1V2V3 = be_seqdata$activesurv_n_501Y.V1plusV2plusV3 / be_seqdata$activesurv_total_sequenced


be_seqdata$basplusactivesurv_n_501Y.V1 = be_seqdata$baselinesurv_n_501Y.V1+be_seqdata$activesurv_n_501Y.V1
be_seqdata$basplusactivesurv_n_501Y.V2 = be_seqdata$baselinesurv_n_501Y.V2+be_seqdata$activesurv_n_501Y.V2
be_seqdata$basplusactivesurv_n_501Y.V3 = be_seqdata$baselinesurv_n_501Y.V3+be_seqdata$activesurv_n_501Y.V3
be_seqdata$basplusactivesurv_total_sequenced = be_seqdata$baselinesurv_total_sequenced+be_seqdata$activesurv_total_sequenced 
be_seqdata$basplusactivesurv_n_wild_type = be_seqdata$basplusactivesurv_total_sequenced-be_seqdata$basplusactivesurv_n_501Y.V1-be_seqdata$basplusactivesurv_n_501Y.V2-be_seqdata$basplusactivesurv_n_501Y.V3
be_seqdata$basplusactivesurv_prop501YV1 = be_seqdata$basplusactivesurv_n_501Y.V1 / be_seqdata$basplusactivesurv_total_sequenced
be_seqdata$basplusactivesurv_prop501YV2 = be_seqdata$basplusactivesurv_n_501Y.V2 / be_seqdata$basplusactivesurv_total_sequenced
be_seqdata$basplusactivesurv_prop501YV3 = be_seqdata$basplusactivesurv_n_501Y.V3 / be_seqdata$basplusactivesurv_total_sequenced
be_seqdata$basplusactivesurv_n_501Y.V1plusV2plusV3 = be_seqdata$basplusactivesurv_n_501Y.V1+be_seqdata$basplusactivesurv_n_501Y.V2+be_seqdata$basplusactivesurv_n_501Y.V3
be_seqdata$basplusactivesurv_n_501Y.propV1V2V3 = be_seqdata$basplusactivesurv_n_501Y.V1plusV2plusV3 / be_seqdata$basplusactivesurv_total_sequenced


head(be_seqdata)


# BASELINE SURVEILLANCE DATA ####

be_basseqdata_long = gather(be_seqdata[,c("collection_date",
                                       "baselinesurv_n_wild_type",
                                       "baselinesurv_n_501Y.V1",
                                       "baselinesurv_n_501Y.V2",
                                       "baselinesurv_n_501Y.V3",
                                       "baselinesurv_n_501Y.V1plusV2plusV3",
                                       "baselinesurv_total_sequenced")], 
                         variant, count, c("baselinesurv_n_wild_type",
                                           "baselinesurv_n_501Y.V1",
                                           "baselinesurv_n_501Y.V2",
                                           "baselinesurv_n_501Y.V3",
                                           "baselinesurv_n_501Y.V1plusV2plusV3"), factor_key=TRUE)
be_basseqdata_long$variant = factor(be_basseqdata_long$variant, 
                                 levels=c("baselinesurv_n_wild_type","baselinesurv_n_501Y.V1","baselinesurv_n_501Y.V2","baselinesurv_n_501Y.V3","baselinesurv_n_501Y.V1plusV2plusV3"), 
                                 labels=c("wild type", "501Y.V1", "501Y.V2", "501Y.V3", "501Y.V1+V2+V3"))
be_basseqdata_long$collection_date_num = as.numeric(be_basseqdata_long$collection_date)
be_basseqdata_long$prop = be_basseqdata_long$count / be_basseqdata_long$baselinesurv_total_sequenced

baseline_sequencing = ggplot(data=be_basseqdata_long[be_basseqdata_long$variant!="501Y.V1+V2+V3",], 
                       aes(x=collection_date, 
                           y=count, fill=variant, group=variant)) +
  # facet_wrap(~LABORATORY) +
  geom_area(aes(fill=variant), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("variant", values=c("darkgrey","blue","red","green3"), 
                               labels=c("wild type","501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-03-01")), expand=c(0,0)) +
  ylab("Share among newly diagnosed infections") +
  xlab("Collection date") +
  # ggtitle("Test outcomes") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right")
baseline_sequencing

saveRDS(baseline_sequencing, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3.pdf"), width=7, height=5)


# multinomial spline fit on share of each type (wild type / British / SA / Brazilian
# to be able to estimate growth rate advantage of each type compared to wild type

set.seed(1)
be_seq_mfit0 = nnet::multinom(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                        subset=be_basseqdata_long$variant!="501Y.V1+V2+V3", maxit=1000)
be_seq_mfit1 = nnet::multinom(variant ~ scale(ns(collection_date_num, df=2), center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                        subset=be_basseqdata_long$variant!="501Y.V1+V2+V3", maxit=1000) 
be_seq_mfit2 = nnet::multinom(variant ~ scale(ns(collection_date_num, df=3), center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                        subset=be_basseqdata_long$variant!="501Y.V1+V2+V3", maxit=1000) 
be_seq_mfitallVOC = nnet::multinom(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE), weights=count, data=be_basseqdata_long, 
                                              subset=be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="501Y.V1+V2+V3", maxit=1000) 

BIC(be_seq_mfit0, be_seq_mfit1, be_seq_mfit2) # mfit0 fits best
#              df      BIC
# be_seq_mfit0  6 2405.575
# be_seq_mfit1  9 2420.704
# be_seq_mfit2 12 2443.147
summary(be_seq_mfit0)

# growth rate advantage compared to wild type
delta_r_501V1_501YV2_501YV3 = data.frame(confint(emtrends(be_seq_mfit0, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_501V1_501YV2_501YV3) = delta_r_501V1_501YV2_501YV3[,"contrast"]
delta_r_501V1_501YV2_501YV3 = delta_r_501V1_501YV2_501YV3[,-1]
delta_r_501V1_501YV2_501YV3
#                       estimate   asymp.LCL  asymp.UCL
# 501Y.V1 - wild type 0.08294642 0.06919210 0.09670073
# 501Y.V2 - wild type 0.07684561 0.04674708 0.10694414
# 501Y.V3 - wild type 0.18033281 0.03023209 0.33043353

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_501V1_501YV2_501YV3*4.7) 
#                     estimate asymp.LCL asymp.UCL
# 501Y.V1 - wild type 1.476757  1.384311  1.575375
# 501Y.V2 - wild type 1.435014  1.245717  1.653075
# 501Y.V3 - wild type 2.333955  1.152681  4.725803

# for all 3 variants together
# growth rate advantage compared to wild type
delta_r_allVOCs = data.frame(confint(emtrends(be_seq_mfitallVOC, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent"), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_allVOCs) = delta_r_allVOCs[,"contrast"]
delta_r_allVOCs = delta_r_allVOCs[,-1]
delta_r_allVOCs
#                               estimate  asymp.LCL  asymp.UCL
# (501Y.V1+V2+V3) - wild type 0.08273544 0.06978487 0.09568601

# implied transmission advantage (assuming no immune evasion advantage of 501Y.V2, if there is such an advantage, transm advantage would be less)
exp(delta_r_allVOCs*4.7) 
#                             estimate asymp.LCL asymp.UCL
# (501Y.V1+V2+V3) - wild type   1.475293  1.388174   1.56788




# PS: mblogit fit would also be possible & would take into account overdispersion
be_basseqdata_long$obs = factor(1:nrow(be_basseqdata_long))
be_seq_mblogitfit = mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
                            # random = ~ 1|obs,
                            weights = count, data=be_basseqdata_long, 
                            dispersion = TRUE)
be_seq_mblogitfit1 = mblogit(variant ~ scale(ns(collection_date_num, df=2), center=TRUE, scale=FALSE),
                             # random = ~ 1|obs,
                             weights = count, data=be_basseqdata_long, 
                             dispersion = TRUE)
BIC(be_seq_mblogitfit, be_seq_mblogitfit1)
#                    df      BIC
# be_seq_mblogitfit   6 2405.575
# be_seq_mblogitfit1  9 2420.665
summary(be_seq_mblogitfit)

dispersion(mblogit(variant ~ scale(collection_date_num, center=TRUE, scale=FALSE),
                   # random = ~ 1|obs,
                   weights = count, data=be_basseqdata_long,
                   subset = be_basseqdata_long$variant=="wild type"|be_basseqdata_long$variant=="501Y.V1+V2+V3",
                   dispersion = TRUE), method="Afroz") # dispersion coefficient = 3.2

# growth rate advantage compared to wild type
delta_r_501V1_501YV2_501YV3 = data.frame(confint(emtrends(be_seq_mblogitfit, trt.vs.ctrl ~ variant|1, var="collection_date_num",  mode="latent", df=NA), adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_501V1_501YV2_501YV3) = delta_r_501V1_501YV2_501YV3[,"contrast"]
delta_r_501V1_501YV2_501YV3 = delta_r_501V1_501YV2_501YV3[,-1]
delta_r_501V1_501YV2_501YV3
#                       estimate  asymp.LCL  asymp.UCL
# 501Y.V1 - wild type 0.08294652 0.06919221 0.09670084
# 501Y.V2 - wild type 0.07684610 0.04674750 0.10694470
# 501Y.V3 - wild type 0.18032114 0.03022225 0.33042003




# plot multinomial model fit

# library(effects)
# plot(Effect("collection_date_num",be_seq_mfit0), style="stacked")

extrapolate = 90
date.from = as.numeric(as.Date("2020-11-01")) # min(be_basseqdata_long$collection_date_num)
date.to = max(be_basseqdata_long$collection_date_num)+extrapolate

be_seq_mfit0_preds = data.frame(emmeans(be_seq_mfit0, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfit0_preds$collection_date = as.Date(be_seq_mfit0_preds$collection_date_num, origin="1970-01-01")
be_seq_mfit0_preds$variant = factor(be_seq_mfit0_preds$variant, levels=c("wild type","501Y.V1","501Y.V2","501Y.V3"),
                                    labels=c("wild type","501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)"))

be_seq_mfitallVOCs_preds = data.frame(emmeans(be_seq_mfitallVOC, ~ variant+collection_date_num, at=list(collection_date_num=seq(date.from, date.to)), mode="prob", df=NA))
be_seq_mfitallVOCs_preds$collection_date = as.Date(be_seq_mfitallVOCs_preds$collection_date_num, origin="1970-01-01")
be_seq_mfitallVOCs_preds$variant = factor(be_seq_mfitallVOCs_preds$variant, levels=c("wild type","501Y.V1+V2+V3"),
                                    labels=c("wild type","501Y.V1+V2+V3"))

be_basseqdata_long2 = be_basseqdata_long[be_basseqdata_long$variant!="wild type",]
be_basseqdata_long2$variant = droplevels(be_basseqdata_long2$variant) 
be_basseqdata_long2$variant = factor(be_basseqdata_long2$variant, levels=c("501Y.V1","501Y.V2","501Y.V3","501Y.V1+V2+V3"),
                                  labels=c("501Y.V1 (British)","501Y.V2 (South African)","501Y.V3 (Brazilian)","501Y.V1+V2+V3"))

muller_be_seq_mfit0 = ggplot(data=be_seq_mfit0_preds, 
                     aes(x=collection_date, y=prob, group=variant)) + 
  # facet_wrap(~LABORATORY) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
  annotate("rect", xmin=max(be_basseqdata_long$collection_date)+1, 
                   xmax=as.Date("2021-03-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_fill_manual("variant", values=c("darkgrey","blue","red","green3")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                     limits=as.Date(c("2020-12-01","2021-03-01")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Relative abundance") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")
muller_be_seq_mfit0

saveRDS(muller_be_seq_mfit0, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.pptx"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.png"), width=7, height=5)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit.pdf"), width=7, height=5)

multinom_501YV1_501YV2_501YV3 = ggarrange(baseline_sequencing+ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)")+xlab("")+ylab("Share"), 
          muller_be_seq_mfit0+ggtitle("Multinomial fit plus extrapolation")+ylab("Share"), ncol=1, 
          common.legend = TRUE, legend="bottom")
multinom_501YV1_501YV2_501YV3
saveRDS(multinom_501YV1_501YV2_501YV3, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.pptx"), width=7, height=10)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.png"), width=7, height=10)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_multipanel.pdf"), width=7, height=10)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS
be_seq_mfit0_preds2 = rbind(be_seq_mfit0_preds[be_seq_mfit0_preds$variant!="wild type",],
                            be_seq_mfitallVOCs_preds[be_seq_mfitallVOCs_preds$variant!="wild type",])
be_seq_mfit0_preds2$variant = droplevels(be_seq_mfit0_preds2$variant)


# on response scale:
plot_multinom_501YV1_501YV2_501YV3_response = qplot(data=be_seq_mfit0_preds2, x=collection_date, y=100*prob, geom="blank") +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share of diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(be_basseqdata_long2$collection_date), as.Date("2021-03-01")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("blue","red","green3","black")) +
  scale_colour_manual("variant", values=c("blue","red","green3","black")) +
  geom_point(data=be_basseqdata_long2,
             aes(x=collection_date, y=100*prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("number\nsequenced", trans="sqrt",
                        range=c(1, 4), limits=c(1,10^3), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_multinom_501YV1_501YV2_501YV3_response

saveRDS(plot_multinom_501YV1_501YV2_501YV3_response, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2_501YV3_multinomial fit_model preds_response.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pdf"), width=8, height=6)


# on logit scale:

be_seq_mfit0_preds3 = be_seq_mfit0_preds2
ymin = 0.001
ymax = 0.900001
be_seq_mfit0_preds3$asymp.LCL[be_seq_mfit0_preds3$asymp.LCL<ymin] = ymin
be_seq_mfit0_preds3$asymp.UCL[be_seq_mfit0_preds3$asymp.UCL<ymin] = ymin
be_seq_mfit0_preds3$asymp.UCL[be_seq_mfit0_preds3$asymp.UCL>ymax] = ymax
be_seq_mfit0_preds3$prob[be_seq_mfit0_preds3$prob<ymin] = ymin

plot_multinom_501YV1_501YV2_501YV3 = qplot(data=be_seq_mfit0_preds3, x=collection_date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=variant
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=variant
  ), alpha=I(1)) +
  ylab("Share of diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("blue","red","green3","black")) +
  scale_colour_manual("variant", values=c("blue","red","green3","black")) +
  geom_point(data=be_basseqdata_long2,
             aes(x=collection_date, y=prop, size=baselinesurv_total_sequenced,
                 colour=variant
             ),
             alpha=I(1)) +
  scale_size_continuous("total number\nsequenced", trans="sqrt",
                        range=c(1, 4), limits=c(10,10^3), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-03-01")), ylim=c(0.001, 0.9000001), expand=c(0,0))
plot_multinom_501YV1_501YV2_501YV3


saveRDS(plot_multinom_501YV1_501YV2_501YV3, file = paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_sequencing_501YV1 501YV2 501YV3_multinomial fit_model preds.pdf"), width=8, height=6)


# estimated proportion of 501Y.V1 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18678 0.5937361 0.05109728 NA 0.4935873 0.6938849      2021-02-20

# estimated proportion of 501Y.V1 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1 (British)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V1 (British)               18685 0.6500041 0.1175849 NA 0.4195419 0.8804663      2021-02-27

# estimated proportion of 501Y.V2 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18678 0.08207829 0.02675573 NA 0.02963803 0.1345186      2021-02-20

# estimated proportion of 501Y.V2 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V2 (South African)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V2 (South African)               18685 0.0861002 0.03854629 NA 0.01055087 0.1616495      2021-02-27

# estimated proportion of 501Y.V3 among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18678 0.04207796 0.05265073 NA -0.06111558 0.1452715      2021-02-20

# estimated proportion of 501Y.V3 among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V3 (Brazilian)",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 501Y.V3 (Brazilian)               18685 0.09108313 0.1506244 NA -0.2041353 0.3863016      2021-02-27



# estimated proportion of one of the three VOCs among new lab diagnoses today
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==today&
                                                 be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18678 0.7100834 0.031975 NA 0.6474136 0.7727533      2021-02-20

# estimated proportion of one of the three VOCs among new infections today (counted one week before lab diagnosis)
be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date==(today+7)&
                      be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
#            variant collection_date_num      prob       SE df asymp.LCL asymp.UCL collection_date
# 2241 501Y.V1+V2+V3               18678 0.8138067 0.03018072 NA 0.7546536 0.8729598      2021-02-27




# the time at which new lab diagnoses would be by more than 50%, 75% 90% by one of the three VOCs :
be_seq_mfit0_preds2_subs = be_seq_mfit0_preds2[be_seq_mfit0_preds2$collection_date>as.Date("2021-01-01")&
                                               be_seq_mfit0_preds2$variant=="501Y.V1+V2+V3",]
# >50% by 10th of February [8th Feb - 12th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]

# >75% by 23d of February [19th Feb - 27th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]

# >90% by 8th of March [2nd March - 14th of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]


# the time at which new infections would be by more than 50%, 75% 90% by one of the three VOCs
# (counting 7 days between infection & diagnosis) :
# >50% by 3d of Feb [1st Feb - 5th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.5)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.5)[1]]-7

# >75% by 16th of February [12th Feb - 20th Feb] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.75)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.75)[1]]-7

# >90% by 1st of March [23 Febr - 7th of March] 95% CLs
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$prob>=0.90)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.UCL>=0.90)[1]]-7
be_seq_mfit0_preds2_subs$collection_date[which(be_seq_mfit0_preds2_subs$asymp.LCL>=0.90)[1]]-7





# 6. SOME INTERNATIONAL COMPARISONS ####

# 6.1. DATA UK: PILLAR 2 SGTF DATA ####

levels_UKregions = c("South East","London","East of England",
                     "South West","Midlands","North East and Yorkshire",
                     "Scotland","North West","Wales")

sgtfdata_uk = read.csv(".//data//uk//sgtf_pillar2_UK-2021-01-25.csv") # Pillar 2 S gene targeted failure data (SGTF) (S dropout)
sgtfdata_uk$other = sgtfdata_uk$other+sgtfdata_uk$sgtf
colnames(sgtfdata_uk) = c("collection_date","REGION","SGTF","TOTAL")
sgtfdata_uk_truepos = read.csv(".//data//uk//sgtf_pillar2_UK-2021-01-25_nick davies_modelled true pos rate sgtfv.csv") # modelled proportion of S dropout that was actually the VOC
# PS this could also be estimated from the COG-UK data based on the presence of deletion 69/70, which is S dropout
sgtfdata_uk$TRUEPOS = sgtfdata_uk_truepos$sgtfv[match(interaction(sgtfdata_uk$REGION, sgtfdata_uk$collection_date),
                                                      interaction(sgtfdata_uk_truepos$group, sgtfdata_uk_truepos$date))] # modelled proportion of S dropout samples that were actually the VOC
sgtfdata_uk$est_n_B117 = sgtfdata_uk$SGTF * sgtfdata_uk$TRUEPOS
sgtfdata_uk$COUNTRY = "UK"
sgtfdata_uk = sgtfdata_uk[,c("collection_date","COUNTRY","REGION","est_n_B117","TOTAL")]
colnames(sgtfdata_uk)[which(colnames(sgtfdata_uk)=="TOTAL")] = "n_pos"
range(sgtfdata_uk$collection_date) # "2020-10-01" "2021-01-24"
sgtfdata_uk$collection_date = as.Date(sgtfdata_uk$collection_date)
sgtfdata_uk$collection_date_num = as.numeric(sgtfdata_uk$collection_date)
sgtfdata_uk$REGION = factor(sgtfdata_uk$REGION, levels=levels_UKregions)
sgtfdata_uk$obs = factor(1:nrow(sgtfdata_uk))
sgtfdata_uk$propB117 = sgtfdata_uk$est_n_B117 / sgtfdata_uk$n_pos
head(sgtfdata_uk)

set_sum_contrasts()
glmersettings = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e5)) # bobyqa, PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B"))
glmersettings3 = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb"))
glmersettings4 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
fit_ukSGTF_1 = glmer(cbind(est_n_B117, n_pos-est_n_B117 ) ~ (1|obs)+scale(collection_date_num)+REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings)  # common slope model for country
fit_ukSGTF_2 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)*REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings) # separate slopes model for country
fit_ukSGTF_3 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+ns(collection_date_num,df=3)+REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings) # with additive spline term
fit_ukSGTF_4 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+ns(collection_date_num,df=3)*REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings3) # with spline term in interaction with region
BIC(fit_ukSGTF_1, fit_ukSGTF_2, fit_ukSGTF_3, fit_ukSGTF_4) 
# separate-slopes 3 df spline model fit_be_uk2_4 best
# df      BIC
# fit_ukSGTF_1  9 7048.188
# fit_ukSGTF_2 15 6574.709
# fit_ukSGTF_3 11 6929.594
# fit_ukSGTF_4 29 5868.880


# model fit_ukSGTF_4 best

summary(fit_ukSGTF_4)

# avg growth rate advantage for UK (difference in growth rate between B.1.1.7 and old strains):
fit_ukSGTF_4_emtrends = as.data.frame(emtrends(fit_ukSGTF_4, revpairwise ~ 1, 
                                               var="collection_date_num",
                                               mode="link", adjust="Tukey")$emtrends)
fit_ukSGTF_4_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                 0.11106 0.1096202 0.1124998

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_ukSGTF_4_emtrends[,c(2,5,6)]*4.7) 
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.685365  1.673998  1.696808


# PLOT MODEL FIT

# spline model fit_ukSGTF_4
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_ukSGTF_4))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_ukSGTF_4_preds = as.data.frame(emmeans(fit_ukSGTF_4, ~ collection_date_num, 
                                           by=c("REGION"), 
                                           at=list(collection_date_num=seq(date.from,
                                                                       date.to)), 
                                           type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_ukSGTF_4_preds$collection_date = as.Date(fit_ukSGTF_4_preds$collection_date_num, origin="1970-01-01")

n = length(levels(fit_ukSGTF_4_preds$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
# reg_cols[2:n] = rev(reg_cols[2:n])


# PLOT MODEL FIT (logit scale):
plot_UK_SGTF = qplot(data=fit_ukSGTF_4_preds, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~COUNTRY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("UK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_UK_SGTF


# PLOT MODEL FIT (response scale):
plot_UK_SGTF_response = qplot(data=fit_ukSGTF_4_preds, x=collection_date, y=prob*100, geom="blank") +
  # facet_wrap(~COUNTRY) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117*100, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("UK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_UK_SGTF_response



# 6.2. DATA UK: COG-UK SEQUENCING DATA ####
# (TO DO - THIS HAS NOT BEEN UPDATED YET)

data_uk = read.csv(".//data//uk//COGUKdata_agbydayregion.csv") 
data_uk = data_uk[data_uk$variant=="VOC 202012/01",]
# COG-UK sequencing data, aggregated by NHS region, from https://github.com/nicholasdavies/newcovid/tree/master/multinomial_logistic_fits/data
head(data_uk)
data_uk$COUNTRY = "UK"
data_uk = data_uk[,c("collection_date","COUNTRY","nhs_name","count","total")]
colnames(data_uk) = c("collection_date","COUNTRY","REGION","VOC","TOTAL")
# XXX




# 6.3. DATA SWITZERLAND ####

# data from https://ispmbern.github.io/covid-19/variants/ & https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv

data_geneva = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_GE.csv")
data_geneva$date = as.Date(data_geneva$date)
data_geneva$lab = "Geneva"
colnames(data_geneva)[colnames(data_geneva) %in% c("N501Y")] = c("n_B117")
head(data_geneva)
data_zurich = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_ZH.csv")
data_zurich$date = as.Date(data_zurich$date)
data_zurich$lab = "Zürich"
colnames(data_zurich)[colnames(data_zurich) %in% c("N501Y")] = c("n_B117")
head(data_zurich)
data_bern = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_BE.csv")
data_bern$date = as.Date(data_bern$date)
data_bern$lab = "Bern"
colnames(data_bern)[colnames(data_bern) %in% c("N501Y")] = c("n_B117")
head(data_bern)

data_viollier_risch = read.csv("https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv")
data_viollier_risch[is.na(data_viollier_risch)] = 0
data_viollier_risch$date = as.Date(NA)
data_viollier_risch$date[data_viollier_risch$week>=51] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_viollier_risch$week[data_viollier_risch$week>=51] - 1 ) + 1
data_viollier_risch$date[data_viollier_risch$week<51] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_viollier_risch$week[data_viollier_risch$week<51] - 1 ) + 6 # PS dates were made to match the ones given in https://ispmbern.github.io/covid-19/variants/data/variants_CH.csv
colnames(data_viollier_risch)[colnames(data_viollier_risch) %in% c("n","b117")] = c("total","n_B117")
data_viollier_risch = data_viollier_risch[,c("date","total","n_B117","lab")]

data_switzerland = rbind(data_geneva, data_zurich, data_bern, data_viollier_risch)[,c("date","lab","n_B117","total")]
data_switzerland$lab = factor(data_switzerland$lab, levels=c("Geneva","Zürich","Bern","Viollier","Risch"),
                                 labels=c("Geneva University Hospitals","University Hospital Zürich","University of Bern","Viollier lab","Risch lab"))
data_switzerland$date_num = as.numeric(data_switzerland$date)
data_switzerland$obs = factor(1:nrow(data_switzerland))
data_switzerland$propB117 = data_switzerland$n_B117 / data_switzerland$total


fit_switerland = glmer(cbind(n_B117,total-n_B117) ~ (1|obs)+lab+scale(date_num), family=binomial(logit), data=data_switzerland)
summary(fit_switerland)

# avg growth rate advantage for SWITZERLAND (difference in growth rate between B.1.1.7 and old strains):
fit_switzerland_emtrends = as.data.frame(emtrends(fit_switerland, revpairwise ~ 1, 
                                               var="date_num",
                                               mode="link", adjust="Tukey")$emtrends)
fit_switzerland_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                 0.1050824 0.09580577  0.114359

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_switzerland_emtrends[,c(2,5,6)]*4.7) 
#   date_num.trend asymp.LCL asymp.UCL
# 1       1.638674  1.568763  1.711701


# PS: note that the effect on Rt is in https://ispmbern.github.io/covid-19/variants/
# (1) assumed to be additive as opposed to multiplicative (not quite correct -
# with gamma distributed gen time and small k Rt=exp(r*GT) and multiplicative
# is the logical choice and mult effect on Rt = exp(logistic regression slope*GT);
# (for derivation see https://cmmid.github.io/topics/covid19/uk-novel-variant.html)
# with exponentially distributed GT Rt=r*GT and multiplicative effect on Rt = 1+logistic slope*GT, 
# (2) that unlike in our analysis overdispersion is ignored and that (3) a generation time 
# of 5.2 days instead of 4.7 days is used.


# PLOT MODEL FIT
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_switerland))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_switzerland_preds = as.data.frame(emmeans(fit_switerland, ~ date_num, 
                                              by=c("lab"), 
                                              at=list(date_num=seq(date.from,
                                                                   date.to)), 
                                              type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_switzerland_preds$date = as.Date(fit_switzerland_preds$date_num, origin="1970-01-01")

n = length(levels(fit_switzerland_preds$lab))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
# reg_cols[2:n] = rev(reg_cols[2:n])

# PLOT MODEL FIT (logit scale):
plot_switzerland = qplot(data=fit_switzerland_preds, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=lab
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=lab
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-11-01"),as.Date("2021-03-01")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_switzerland, 
             aes(x=date, y=propB117, size=total,
                 colour=lab
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("SWITZERLAND") +
  theme(plot.title = element_text(hjust = 0.5))
plot_switzerland


# PLOT MODEL FIT (response scale):
plot_switzerland_response = qplot(data=fit_switzerland_preds, x=date, y=prob*100, geom="blank") +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=lab
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=lab
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-11-01"),as.Date("2021-03-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_switzerland, 
             aes(x=date, y=propB117*100, size=total,
                 colour=lab
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(500,1000,2000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("SWITZERLAND") +
  theme(plot.title = element_text(hjust = 0.5))
plot_switzerland_response





# 6.4. DATA DENMARK: SEQUENCING DATA ####

# analysis of sequencing data from Denmark, split by region
# from https://www.covid19genomics.dk/statistics
data_denmark = read.csv(".//data/dk//B117_denmark_20210211.csv", sep=";", dec=",")
data_denmark = data_denmark[data_denmark$Region!="Whole Denmark",]
data_denmark$percent = NULL
data_denmark$Region = gsub("SjÃ¦lland","Sjælland",data_denmark$Region)
data_denmark = data_denmark[data_denmark$Region!="Other",]
data_denmark$WEEK = sapply(data_denmark$Week, function(s) as.numeric(strsplit(s, "W")[[1]][[2]]))
data_denmark$date = as.Date(NA)
data_denmark$date[data_denmark$WEEK>=42] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK>=42] - 1 ) + 1
data_denmark$date[data_denmark$WEEK<42] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK<42] - 1 ) + 6 
data_denmark$date_num = as.numeric(data_denmark$date)
data_denmark$obs = factor(1:nrow(data_denmark))
levels_DK = c("Syddanmark","Sjælland","Nordjylland","Hovedstaden","Midtjylland")
data_denmark$Region = factor(data_denmark$Region, levels=levels_DK)
colnames(data_denmark)[colnames(data_denmark) %in% c("yes")] = "n_B117"
data_denmark$propB117 = data_denmark$n_B117 / data_denmark$total

fit_denmark = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + Region + scale(date_num), family=binomial(logit), data=data_denmark)
summary(fit_denmark)

# avg growth rate advantage for DENMARK (difference in growth rate between B.1.1.7 and old strains):
fit_denmark_emtrends = as.data.frame(emtrends(fit_denmark, revpairwise ~ 1, 
                                                  var="date_num",
                                                  mode="link", adjust="Tukey")$emtrends)
fit_denmark_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                 0.07699966 0.0670402 0.08695911

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_denmark_emtrends[,c(2,5,6)]*4.7) 
#   date_num.trend asymp.LCL asymp.UCL
# 1       1.436053  1.370381  1.504872

# predicted prop of diagnosed samples that are B.1.1.7 today
emmeans(fit_denmark, revpairwise ~ Region, at=list(date_num=today_num), type="response")$emmeans
# Region       prob     SE  df asymp.LCL asymp.UCL
# Syddanmark  0.389 0.0676 Inf     0.267     0.526
# Sjælland    0.379 0.0674 Inf     0.258     0.517
# Nordjylland 0.376 0.0715 Inf     0.249     0.523
# Hovedstaden 0.327 0.0641 Inf     0.216     0.463
# Midtjylland 0.247 0.0542 Inf     0.156     0.367

emmeans(fit_denmark, revpairwise ~ 1, at=list(date_num=today_num), type="response")$emmeans
# 1        prob     SE  df asymp.LCL asymp.UCL
# overall 0.341 0.0449 Inf      0.26     0.434

# predicted prop of infections that are B.1.1.7 today
emmeans(fit_denmark, revpairwise ~ Region, at=list(date_num=today_num+7), type="response")$emmeans
# Region       prob     SE  df asymp.LCL asymp.UCL
# Syddanmark  0.522 0.0764 Inf     0.374     0.665
# Sjælland    0.511 0.0766 Inf     0.365     0.656
# Nordjylland 0.508 0.0819 Inf     0.352     0.662
# Hovedstaden 0.455 0.0783 Inf     0.310     0.608
# Midtjylland 0.360 0.0716 Inf     0.234     0.508

emmeans(fit_denmark, revpairwise ~ 1, at=list(date_num=today_num+7), type="response")$emmeans
# 1        prob     SE  df asymp.LCL asymp.UCL
# overall 0.47 0.0574 Inf     0.361     0.582


# PLOT MODEL FIT
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_denmark))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_denmark_preds = as.data.frame(emmeans(fit_denmark, ~ date_num, 
                                              by=c("Region"), 
                                              at=list(date_num=seq(date.from,
                                                                   date.to)), 
                                              type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_denmark_preds$date = as.Date(fit_denmark_preds$date_num, origin="1970-01-01")

n = length(levels(fit_denmark_preds$Region))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
# reg_cols[2:n] = rev(reg_cols[2:n])

# PLOT MODEL FIT (response scale)
plot_denmark = qplot(data=fit_denmark_preds, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=Region
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=Region
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-10-01"),as.Date("2021-03-01")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_denmark, 
             aes(x=date, y=propB117, size=total,
                 colour=Region
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("DENMARK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_denmark


# PLOT MODEL FIT (response scale)
plot_denmark_response = qplot(data=fit_denmark_preds, x=date, y=prob*100, geom="blank") +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=Region
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=Region
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-10-01"),as.Date("2021-03-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_denmark, 
             aes(x=date, y=propB117*100, size=total,
                 colour=Region
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("DENMARK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_denmark_response






# 6.5. DATA US ####

# US data from https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data
# https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1

helix_b117 = read_tsv("https://github.com/andersen-lab/paper_2021_early-b117-usa/raw/master/b117_frequency/data/covid_baseline_for_b117_paper.20210127_update.txt") %>%
  select(state, collection_date, n_b117, n_sgtf_seq) # n_b117/n_sgtf_seq = prop of S dropout samples that are B117

helix_sgtf = read_tsv("https://github.com/andersen-lab/paper_2021_early-b117-usa/raw/master/b117_frequency/data/covid_baseline_for_b117_paper.20210201_klados20211029_phyloseq.txt") %>%
  select(state, collection_date, n, n_sgtf) # n_sgtf/n = prop of pos tests that have S dropout
helix_sgtf = helix_sgtf[helix_sgtf$state %in% unique(helix_b117$state),]

helix_metadata = left_join(helix_sgtf, helix_b117, by=c("state", "collection_date"))

tmp = helix_metadata %>%
  group_by(collection_date) %>%
  summarise(n_sgtf = sum(n_sgtf), n = sum(n)) %>%
  mutate(state = "USA")

tmp = bind_rows(tmp, helix_metadata)
states_gt_500 = tmp %>% group_by(state) %>% summarise(n = sum(n), n_sgtf = sum(n_sgtf)) %>% filter(n > 500 & n_sgtf > 0) %>% select(state) %>% as_vector()

helix_b117$collection_date_num = as.numeric(helix_b117$collection_date)
helix_b117$obs = factor(1:nrow(helix_b117))
helix_b117$state = factor(helix_b117$state)
helix_sgtf$collection_date_num = as.numeric(helix_sgtf$collection_date)
helix_sgtf$state = factor(helix_sgtf$state)
helix_sgtf$obs = factor(1:nrow(helix_sgtf))


fit_us_propB117amongSGTF = glmer(cbind(n_b117, n_sgtf_seq-n_b117) ~ (1|state)+scale(collection_date_num), 
                                family=binomial(logit), data=helix_b117)

# implied growth rate advantage of B.1.1.7 over other earlier strains showing S dropout:
as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)]
#   collection_date_num.trend  asymp.LCL asymp.UCL
# 1                0.08999387 0.06042387 0.1195639

# with a generation time of 4.7 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of B.1.1.7 over other strains showing S dropout of
exp(4.7*as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.52649  1.328423   1.75409


# FIT FOR WHOLE US + PLOT ####

fitted_truepos = predict(fit_us_propB117amongSGTF, newdat=helix_sgtf, type="response") 
# fitted true positive rate, ie prop of S dropout samples that are B.1.1.7 for dates & states in helix_sgtf

helix_sgtf$est_n_B117 = helix_sgtf$n_sgtf*fitted_truepos # estimated nr of B.1.1.7 samples
helix_sgtf$propB117 = helix_sgtf$est_n_B117/helix_sgtf$n 
fit_us = glmer(cbind(est_n_B117, n-est_n_B117) ~ (1|state/obs)+scale(collection_date_num), 
               family=binomial(logit), data=helix_sgtf) # random intercepts by state
fit_us2 = glmer(cbind(est_n_B117, n-est_n_B117) ~ (collection_date_num||state/obs)+scale(collection_date_num), 
               family=binomial(logit), data=helix_sgtf) # random intercepts+slopes by state, with uncorrelated intercepts & slopes
BIC(fit_us, fit_us2) # random intercept model is best
# df      BIC
# fit_us   4 1171.852
# fit_us2  6 1188.542
summary(fit_us)

# growth advantage of 8.6% per day [8.0-9.2%] 95% CLs
as.data.frame(emtrends(fit_us, ~ 1, "collection_date_num"))[,c(2,5,6)]
#   date_num.trend  asymp.LCL  asymp.UCL
# 1     0.08578087 0.07995376 0.09160798

# increased infectiousness for GT=4.7 days: 49% more infectious [46-54%] 95% CLs
# PS: most logical value to use for GT is the one that has always been used to calculate Rt values in the US, what value is that?
exp(4.7*as.data.frame(emtrends(fit_us, ~ 1, "collection_date_num"))[,c(2,5,6)])
#   date_num.trend asymp.LCL asymp.UCL
# 1        1.496561  1.456131  1.538115
# increased infectiousness for GT=5 days: 54% more infectious [49-58%] 95% CLs
exp(5*as.data.frame(emtrends(fit_us, ~ 1, "collection_date_num"))[,c(2,5,6)])
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.535574   1.49148  1.580972
# increased infectiousness for GT=6.5 days: 75% more infectious [68-81%] 95% CLs
exp(6.5*as.data.frame(emtrends(fit_us, ~ 1, "collection_date_num"))[,c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.746433  1.681522   1.81385


# plot model fit fit_us

date.to = as.numeric(as.Date("2021-06-01"))
sel_states = intersect(rownames(ranef(fit_us)$state)[order(ranef(fit_us)$state[,1], decreasing=T)],states_gt_500)[1:16] # unique(helix_sgtf$state[helix_sgtf$propB117>0.03])
rem_states = c("NY","NJ","MN","IL","AL","OH","MI") # states with too few data points we don't want to show on plot
sel_states = setdiff(sel_states,rem_states)

total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_us))$sdcor, function (x) x^2))) 
fit_us_preds = as.data.frame(emmeans(fit_us, ~ collection_date_num, 
                                         # by="state", 
                                         at=list(collection_date_num=seq(min(helix_sgtf$collection_date_num),
                                                                         date.to)), 
                                         type="link"), bias.adjust = TRUE, sigma = total.SD)
fit_us_preds$collection_date = as.Date(fit_us_preds$collection_date_num, origin="1970-01-01")
fit_us_preds2 = do.call(rbind,lapply(unique(helix_sgtf$state), function(st) { ranintercs = ranef(fit_us)$state
                                raninterc = ranintercs[rownames(ranintercs)==st,]
                                data.frame(state=st, fit_us_preds, raninterc=raninterc)}))
fit_us_preds2$prob = plogis(fit_us_preds2$emmean+fit_us_preds2$raninterc)
fit_us_preds2$prob.asymp.LCL = plogis(fit_us_preds2$asymp.LCL+fit_us_preds2$raninterc)
fit_us_preds2$prob.asymp.UCL = plogis(fit_us_preds2$asymp.UCL+fit_us_preds2$raninterc)
fit_us_preds2 = fit_us_preds2[as.character(fit_us_preds2$state) %in% sel_states,]
fit_us_preds2$state = droplevels(fit_us_preds2$state)
fit_us_preds2$state = factor(fit_us_preds2$state, # we order states by random intercept, ie date of introduction
                             levels=intersect(rownames(ranef(fit_us)$state)[order(ranef(fit_us)$state[,1], decreasing=T)],
                                              sel_states))

# PLOT MODEL FIT (logit scale)
plot_us = qplot(data=fit_us_preds2, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~state) +
  geom_ribbon(aes(y=prob, ymin=prob.asymp.LCL, ymax=prob.asymp.UCL, colour=NULL, 
                  fill=state
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=state
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_us_preds$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.001,0.9990001), expand=c(0,0)) +
  scale_color_discrete("state", h=c(0, 240), c=120, l=50) +
  scale_fill_discrete("state", h=c(0, 240), c=120, l=50) +
  geom_point(data=helix_sgtf[helix_sgtf$state %in% sel_states,],  
             aes(x=collection_date, y=propB117, size=n,
                 colour=state
             ), pch=I(16),
             # colour=I("steelblue"), 
             alpha=I(0.3)) +
  scale_size_continuous("number of\npositive tests", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10^(round(log10(max(helix_sgtf$n)),0)+1)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  ggtitle("US") +
  theme(plot.title = element_text(hjust = 0.5))
plot_us

saveRDS(plot_us, file = paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.pdf"), width=8, height=6)


# PLOT MODEL FIT (response scale)
plot_us_response = qplot(data=fit_us_preds2, x=collection_date, y=prob*100, geom="blank") +
  facet_wrap(~state) +
  geom_ribbon(aes(y=prob*100, ymin=prob.asymp.LCL*100, ymax=prob.asymp.UCL*100, colour=NULL, 
                  fill=state
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=state
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_us_preds2$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("state", h=c(0, 240), c=120, l=50) +
  scale_fill_discrete("state", h=c(0, 240), c=120, l=50) +
  geom_point(data=helix_sgtf[helix_sgtf$state %in% sel_states,],  
             aes(x=collection_date, y=propB117*100, size=n,
                 colour=state
             ), pch=I(16),
             # colour=I("steelblue"), 
             alpha=I(0.3)) +
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 4), limits=c(1,10^(round(log10(max(helix_sgtf$n)),0)+1)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  ggtitle("US") +
  theme(plot.title = element_text(hjust = 0.5))

plot_us_response






# 6.6. MULTIPANEL PLOT INTERNATIONAL COMPARISONS ####

fit_uk_preds2 = fit_ukSGTF_4_preds
fit_uk_preds2$country = "UK"
colnames(fit_uk_preds2)[2] = "REGION"
colnames(fit_uk_preds2)[1] = "date_num"
colnames(fit_uk_preds2)[8] = "date"
fit_belgium_preds2 = fit1_preds
fit_belgium_preds2$country = "Belgium"
fit_belgium_preds2$REGION = "Belgium"
colnames(fit_belgium_preds2)[1] = "date_num"
colnames(fit_belgium_preds2)[7] = "date"
fit_switzerland_preds2 = fit_switzerland_preds
fit_switzerland_preds2$country = "Switzerland"
colnames(fit_switzerland_preds2)[2] = "REGION"
colnames(fit_switzerland_preds2)[1] = "date_num"
colnames(fit_switzerland_preds2)[8] = "date"
fit_switzerland_preds2$REGION = factor(fit_switzerland_preds2$REGION, 
                                  levels=c("Geneva University Hospitals","University Hospital Zürich","University of Bern","Viollier lab","Risch lab"),
                                  labels=c("Geneva","Zürich","Bern","Swiss Viollier lab","Swiss Risch lab"))
fit_denmark_preds2 = fit_denmark_preds
fit_denmark_preds2$country = "Denmark"
colnames(fit_denmark_preds2)[2] = "REGION"
colnames(fit_denmark_preds2)[1] = "date_num"
colnames(fit_denmark_preds2)[8] = "date"
fit_us_preds3 = fit_us_preds2
fit_us_preds3$country = "USA"
fit_us_preds3 = fit_us_preds3[,-which(colnames(fit_us_preds3) %in% c("asymp.LCL","asymp.UCL"))]
colnames(fit_us_preds3)[1] = "REGION"
colnames(fit_us_preds3)[2] = "date_num"
colnames(fit_us_preds3)[6] = "date"
colnames(fit_us_preds3)[9] = "asymp.LCL"
colnames(fit_us_preds3)[10] = "asymp.UCL"
fit_us_preds3 = fit_us_preds3[fit_us_preds3$REGION %in% c("FL","CA"),]
fit_us_preds3$REGION = factor(fit_us_preds3$REGION, levels=c("FL","CA"), labels=c("Florida","California"))
fit_us_preds3 = fit_us_preds3[,c("date_num","REGION","prob","SE","df","asymp.LCL","asymp.UCL","date","country")]

fits_international = rbind(fit_uk_preds2,fit_denmark_preds2,fit_belgium_preds2,fit_switzerland_preds2,fit_us_preds3)
fits_international$country = factor(fits_international$country, levels=c("UK","Denmark","Belgium","Switzerland","USA"))

sgtfdata_uk2 = sgtfdata_uk
sgtfdata_uk2$country = "UK"
colnames(sgtfdata_uk2)[colnames(sgtfdata_uk2) %in% c("collection_date","n_pos")] = c("date","total")
sgtfdata_uk2 = sgtfdata_uk2[,c("date","country","REGION","propB117","total")]

data_belgium2 = data_ag_byday_wide
data_belgium2$country = "Belgium"
data_belgium2$REGION = "Belgium"
colnames(data_belgium2)[colnames(data_belgium2) %in% c("collection_date","n_pos")] = c("date","total")
data_belgium2 = data_belgium2[,c("date","country","REGION","propB117","total")]

data_switzerland2 = data_switzerland
data_switzerland2$country = "Switzerland"
colnames(data_switzerland2)[colnames(data_switzerland2) %in% c("lab")] = c("REGION")
data_switzerland2 = data_switzerland2[,c("date","country","REGION","propB117","total")]
data_switzerland2$REGION = factor(data_switzerland2$REGION, 
                                  levels=c("Geneva University Hospitals","University Hospital Zürich","University of Bern","Viollier lab","Risch lab"),
                                  labels=c("Geneva","Zürich","Bern","Swiss Viollier lab","Swiss Risch lab"))

data_denmark2 = data_denmark
data_denmark2$country = "Denmark"
colnames(data_denmark2)[colnames(data_denmark2) %in% c("Region")] = c("REGION")
data_denmark2 = data_denmark2[,c("date","country","REGION","propB117","total")]

data_us2 = data.frame(helix_sgtf)
data_us2$country = "USA"
colnames(data_us2)[1] = "REGION"
colnames(data_us2)[2] = "date"
colnames(data_us2)[3] = "total"
data_us2 = data_us2[,c("date","country","REGION","propB117","total")]
data_us2 = data_us2[data_us2$REGION %in% c("FL","CA"),]
data_us2$REGION = factor(data_us2$REGION, levels=c("FL","CA"), labels=c("Florida","California"))

data_international = rbind(sgtfdata_uk2, data_denmark2, data_belgium2, data_switzerland2, data_us2)
data_international$country = factor(data_international$country, levels=c("UK","Denmark","Belgium","Switzerland","USA"))

# n1 = length(levels(fit_uk_preds2$REGION))
# n2 = length(levels(fit_switzerland_preds2$REGION))
# n3 = length(levels(fit_denmark_preds2$REGION))
# reg_cols = c(hcl(h = seq(290, 0, length = n1), l = 50, c = 255),
#              muted(hcl(h = seq(290, 0, length = n2+n3), l = 50, c = 255), c=200, l=40))

ymin = 0.001
ymax = 0.999
data_international$propB117[data_international$propB117>ymax] = ymax
fits_international$prob[fits_international$prob>ymax] = ymax
fits_international$asymp.LCL[fits_international$asymp.LCL>ymax] = ymax
fits_international$asymp.UCL[fits_international$asymp.UCL>ymax] = ymax

# PLOT MODEL FITS (response scale)
plot_international = qplot(data=fits_international, x=date, y=prob, geom="blank") +
  facet_wrap(~country, ncol=1, scales="fixed") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9"),
                      limits = c(ymin,ymax+1E-7)) +
  # scale_color_manual("", values=reg_cols) +
  # scale_fill_manual("", values=reg_cols) +
  scale_color_discrete("region", h=c(0, 240), c=250, l=50) +
  scale_fill_discrete("region", h=c(0, 240), c=250, l=50) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_international, 
             aes(x=date, y=propB117, size=total, shape=country,
                 colour=REGION, fill=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(10,10000), breaks=c(100,1000,10000)) +
  scale_shape_manual(values=21:25) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2),
    fill = guide_legend(order = 2),
    size = guide_legend(order = 3)
  ) + 
  coord_cartesian( 
    xlim=c(as.Date("2020-09-01"),as.Date("2021-03-01")),
    ylim=c(ymin,ymax+1E-7), 
    expand=FALSE)
  # ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
  # theme(plot.title = element_text(hjust = 0.5))
plot_international

saveRDS(plot_international, file = paste0(".\\plots\\",dat,"\\Fig6_binomGLMM_B117_fits_UK_DK_BE_CH_USA.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig6_binomGLMM_B117_fits_UK_DK_BE_CH_USA.pptx"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig6_binomGLMM_B117_fits_UK_DK_BE_CH_USA.png"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig6_binomGLMM_B117_fits_UK_DK_BE_CH_USA.pdf"), width=7, height=8)




# PLOT MODEL FITS (response scale)
plot_international_response = qplot(data=fits_international, x=date, y=prob*100, geom="blank") +
  facet_wrap(~country, ncol=1) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_color_discrete("region", h=c(0, 240), c=120, l=50) +
  scale_fill_discrete("region", h=c(0, 240), c=120, l=50) +
#   scale_color_manual("", values=reg_cols) +
#  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_international, 
             aes(x=date, y=propB117*100, size=total, shape=country,
                 colour=REGION, fill=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(10,10000), breaks=c(100,1000,10000)) +
  scale_shape_manual(values=21:25) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  guides(
  shape = guide_legend(order = 1),
  color = guide_legend(order = 2),
  fill = guide_legend(order = 2),
  size = guide_legend(order = 3)
) +
  coord_cartesian( 
    xlim=c(as.Date("2020-09-01"),as.Date("2021-03-01")),
    ylim=c(0,100), expand=c(0,0))
# +
  # ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
  # theme(plot.title = element_text(hjust = 0.5))
plot_international_response

saveRDS(plot_international_response, file = paste0(".\\plots\\",dat,"\\binomGLMM_B117_fits_UK_DK_BE_CH_USA_response.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\binomGLMM_B117_fits_UK_DK_BE_CH_USA_response.pptx"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\binomGLMM_B117_fits_UK_DK_BE_CH_USA_response.png"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\binomGLMM_B117_fits_UK_DK_BE_CH_USA_response.pdf"), width=7, height=8)


plot_us2 = plot_us + coord_cartesian(xlim=c(as.Date("2020-11-01"), as.Date("2021-03-31")),
                          ylim=c(0.001,99.9), expand=c(0,0)) + ggtitle("SPREAD OF VARIANT B.1.1.7 IN THE US")
plot_us2
saveRDS(plot_us2, file = paste0(".\\plots\\",dat,"\\Fig7_binomGLMM_B117_fit_US.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\Fig7_binomGLMM_B117_fit_US.pptx"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig7_binomGLMM_B117_fit_US.png"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\Fig7_binomGLMM_B117_fit_US.pdf"), width=7, height=8)

