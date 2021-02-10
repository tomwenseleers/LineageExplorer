# ANALYSIS OF S-GENE TARGET FAILURE (S DROPOUT) DATA FROM BELGIUM TO INFER CONTAGIOUSNESS OF NEW VARIANT OF CONCERN B.1.1.7 / 501Y.V1 ####
# PLUS INTERNATIONAL COMPARISON (USING DATA FROM THE UK, DENMARK, SWITZERLAND & THE US)
# T. Wenseleers & N. Hens

# Data provided by Emmanuel André (BE), COG-UK, PHE & N. Davies (UK), 
# Statens Serum Institut & Danish Covid-19 Genome Consortium (DK, https://www.covid19genomics.dk/statistics), 
# Christian Althaus, Swiss Viollier Sequencing Consortium, Institute of Medical Virology, University of Zurich, 
# Swiss National Covid-19 Science Task Force (Switzerland, https://ispmbern.github.io/covid-19/variants/, 
# https://ispmbern.github.io/covid-19/variants/data & https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv)
# and Helix, San Mateo, CA, Karthik Gangavarapu & Kristian G. Andersen (US, https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data, https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1)

# last update 9 FEBR. 2021

library(lme4)
library(splines)
library(purrr)
library(readxl)
library(emmeans)
library(effects)
library(ggplot2)
library(ggthemes)
library(ggpubr)
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
library(tidyverse)
library(lubridate)
library(zoo)
library(gridExtra)
library(sf)

# define emm_basis method to have emmeans support mblogit multinomial mixed models
# cf https://cran.r-project.org/web/packages/emmeans/vignettes/xtending.html
emm_basis.mblogit = function(object, ...) {
  object$coefficients = object$coefmat
  object$lev = levels(object$model[[1]])
  object$edf = Inf
  emmeans:::emm_basis.multinom(object, ...)
}

dat="2021_02_09" # desired file version for Belgian data (date/path in //data)
suppressWarnings(dir.create(paste0(".//plots//",dat)))
today = as.Date(gsub("_","-",dat))
today_num = as.numeric(today)


# 1. ESTIMATE PROPORTION OF S DROPOUT SAMPLES THAT ARE B.1.1.7 / 501Y.V1 IN BELGIUM IN FUNCTION OF TIME BASED ON SEQUENCING DATA ####
# SEQUENCING DATA PROVIDED BY EMMANUEL ANDRÉ

datBE_b117 = read.csv(paste0(".//data//", dat, "//sequencing_Sdropouts.csv"), check.names=F) # n_b117/n_sgtf_seq = prop of S dropout samples that were B.1.1.7
datBE_b117$country = "Belgium"
datBE_b117$collection_date = as.Date(datBE_b117$collection_date)
datBE_b117$collection_date_num = as.numeric(datBE_b117$collection_date)
datBE_b117$propB117 = datBE_b117$n_b117/datBE_b117$n_sgtf_seq
datBE_b117$obs = factor(1:nrow(datBE_b117))
datBE_b117

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
extrapolate = 30 # nr of days to extrapolate fit into the future
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
# 62           18667 0.9883506 0.006573705 Inf 0.9651977  0.996167      2021-02-09

# prop of S dropout samples among new infections that are now estimated to be B.1.1.7 / 501Y.V1 (using 7 days for time from infection to diagnosis)
fitseq_preds[fitseq_preds$collection_date==(today+7),]
#    collection_date_num     prob          SE  df asymp.LCL asymp.UCL collection_date
# 69           18674 0.9945993 0.003662389 Inf 0.9797699 0.9985754      2021-02-16

# from 13th of Jan 2021 >80% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.80,"collection_date"][1]

# from 20th of Jan 2021 >90% of all S dropout samples were indeed B.1.1.7 / 501Y.V1
fitseq_preds[fitseq_preds$prob>0.90,"collection_date"][1]

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
    xlim=c(as.Date("2020-12-01"),as.Date("2021-02-08")), 
    ylim=c(0.01,0.99002), expand=c(0,0)) +
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
    xlim=c(as.Date("2020-12-01"),as.Date("2021-02-08")), 
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
sheets = excel_sheets(file_jan)
ctdata_jan = map_df(sheets, ~ read_excel(file_jan, sheet = .x, skip = 1, 
                                       col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                       col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
ctdata_feb = map_df(sheets, ~ read_excel(file_feb, sheet = .x, skip = 1, 
                                         col_names=c("Analysis_date","Laboratory","Outcome","ORF1_cq","S_cq","N_cq","S_dropout"), 
                                         col_types=c("text","text","text","numeric","numeric","numeric","numeric"))) 
ctdata = bind_rows(ctdata_jan, ctdata_feb)
unique(ctdata$Laboratory) 
# "UMons - Jolimont" "UZA"              "ULB"              "Saint LUC - UCL"  "UZ leuven"        "UZ Gent"          "Namur"           
# "ULG - FF 3.x" 
ctdata$Outcome[ctdata$Outcome=="Detected"] = "Positive"
ctdata$Outcome[ctdata$Outcome=="Not detected"] = "Negative"
unique(ctdata$Outcome)
ctdata$Analysis_date = as.Date(as.numeric(ctdata$Analysis_date), origin="1899-12-30")
range(ctdata$Analysis_date) # "2021-01-01" - "2021-02-08"
ctdata$collection_date = ctdata$Analysis_date-1 # collection date = analysis date-1 
ctdata$collection_date_num = as.numeric(ctdata$collection_date)
range(ctdata$collection_date) # "2020-12-31" "2021-02-07"
ctdata$group = interaction(ctdata$Outcome, ctdata$S_dropout)
ctdata$high_viral_load_N = ifelse(ctdata$N_cq<20, 1, 0)
ctdata$high_viral_load_ORF1 = ifelse(ctdata$ORF1_cq<20, 1, 0)
ctdata$Outcome = factor(ctdata$Outcome)
ctdata$Laboratory = factor(ctdata$Laboratory)
ctdata$S_dropout = factor(ctdata$S_dropout)
head(ctdata)
str(ctdata)
nrow(ctdata) # 387653

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

cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.74

# below we continue with the data from all labs, except 
# UMons & ULG (which both had low sample size) and 
# UZ Ghent & UZA (which both were very heavily involved in active surveillance, 
# which may results in biased, nonrepresentative samples)

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                        Laboratory correlation_Ct_N_ORF1ab
# Namur                       Namur             0.945827328
# Saint LUC - UCL   Saint LUC - UCL             0.981151728
# ULB                           ULB             0.981289669
# ULG - FF 3.x         ULG - FF 3.x             0.955847687
# UMons - Jolimont UMons - Jolimont            -0.264402762
# UZ Gent                   UZ Gent            -0.008149941
# UZ leuven               UZ leuven             0.976820844
# UZA                           UZA             0.836841503

ctcorplot_all_labs = qplot(data=ctdata_onlypos_subs, x=ORF1_cq, y=N_cq, group=Laboratory, fill=S_dropout, colour=S_dropout, size=I(3), shape=I(16)) +
  facet_wrap(~Laboratory) + 
  scale_colour_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  scale_fill_manual("", values=alpha(c("blue","red"), 0.05), breaks=c("0","1"), labels=c("S positive","S dropout")) +
  xlab("Ct value ORF1ab gene") + ylab("Ct value N gene") +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5,fill=NA))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ctcorplot_all_labs
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



# plot & analysis of Ct values of all labs except Mons (low S dropout counts) & UZ Gent (heavily involved in active surveillance, so not representative)
# for dates from 13th of Jan onward when >80% of all S dropouts were B.1.1.7 / 501Y.V1
# we also just use the pos samples with Ct values < 30 to be able to focus only on new, active infections

sel_labs = setdiff(unique(ctdata_onlypos_subs$Laboratory), c("UMons - Jolimont", "ULG - FF 3.x", "UZ Gent", "UZA")) 
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
# which was a Ct value < 15.03 for the N gene and < 15.92 for the ORF1ab gene
thresh_N = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","N_cq"]))/1.25
thresh_N # 15.03
thresh_ORF1 = median(unlist(ctdata_onlypos_subs[ctdata_onlypos_subs$S_dropout=="0","ORF1_cq"]))/1.25
thresh_ORF1 # 15.92264
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="N gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="N gene"]<thresh_N
ctdata_onlypos_subs_bothgenes$high_viral_load[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"] = ctdata_onlypos_subs_bothgenes$Ct[ctdata_onlypos_subs_bothgenes$Gene=="ORF1ab gene"]<thresh_ORF1


# check correlation between Ct values for N & ORF1ab gene
cor.test(ctdata_onlypos_subs$N_cq, ctdata_onlypos_subs$ORF1_cq, method="pearson") # Pearson R=0.97, t=436.44, p<2E-16

do.call( rbind, lapply( split(ctdata_onlypos_subs, ctdata_onlypos_subs$Laboratory),
                        function(x) data.frame(Laboratory=x$Laboratory[1], correlation_Ct_N_ORF1ab=cor(x$N_cq, x$ORF1_cq)) ) )
#                      Laboratory correlation_Ct_N_ORF1ab
# Namur                     Namur               0.9458273
# Saint LUC - UCL Saint LUC - UCL               0.9811517
# ULB                         ULB               0.9812897
# UZ leuven             UZ leuven               0.9768208


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
qr_bothgenes0 = rq(ORF1_cq ~ Ggene + S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes1 = rq(ORF1_cq ~ Gene * S_dropout + Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes2 = rq(ORF1_cq ~ Gene + S_dropout * Laboratory, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes3 = rq(ORF1_cq ~ Gene * Laboratory + S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes4 = rq(ORF1_cq ~ (Gene + Laboratory + S_dropout)^2, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
qr_bothgenes5 = rq(ORF1_cq ~ Gene * Laboratory * S_dropout, data=ctdata_onlypos_subs_bothgenes, tau=0.5)
AIC(qr_bothgenes0, k=-1) # 138099.4
AIC(qr_bothgenes1, k=-1) # 138109.3
AIC(qr_bothgenes2, k=-1) # 138020.9 # fits data best based on BIC criterion (PS: here AIC with k<0 returns BIC)
AIC(qr_bothgenes3, k=-1) # 138129.2
AIC(qr_bothgenes4, k=-1) # 138060.7
AIC(qr_bothgenes5, k=-1) # 138090.5

summary(qr_bothgenes)
qr_emmeans_bylab99 = data.frame(emmeans(qr_bothgenes2, ~ Laboratory + Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans_bylab99
qr_emmeans99 = data.frame(emmeans(qr_bothgenes2, ~ Gene + S_dropout, level=0.99)) # median Ct values + 99% CLs
qr_emmeans99
qr_emmeans = data.frame(emmeans(qr_bothgenes2, ~ Gene + S_dropout, level=0.95)) # median Ct values + 95% CLs
qr_emmeans
# Gene S_dropout   emmean        SE    df lower.CL upper.CL
# 1      N gene         0 20.00580 0.1213226 20825 19.76800 20.24360
# 2 ORF1ab gene         0 20.00580 0.1213226 20825 19.76800 20.24360
# 3      N gene         1 18.25605 0.1406013 20825 17.98046 18.53164
# 4 ORF1ab gene         1 18.25605 0.1406013 20825 17.98046 18.53164
contrast(emmeans(qr_bothgenes2, ~ S_dropout, level=0.95), method="pairwise") # difference in median Ct value highly significant across both genes: p<0.0001
# contrast estimate    SE    df t.ratio p.value
# 0 - 1        1.75 0.153 20825 11.456  <.0001 
confint(contrast(emmeans(qr_bothgenes2, ~ S_dropout|Gene, level=0.95), method="pairwise"))
# gene = N gene:
#   contrast estimate    SE    df lower.CL upper.CL
# 0 - 1        1.75 0.153 20825     1.45     2.05
# 
# gene = ORF1ab gene:
#   contrast estimate    SE    df lower.CL upper.CL
# 0 - 1        1.75 0.153 20825     1.45     2.05
# 
# Results are averaged over the levels of: Laboratory 
# Confidence level used: 0.95 


# # Function to produce summary statistics (median + 25 & 75% quantiles)
# data_summary = function(x) {
#   m <- median(x) # mean(x)
#   ymin <- quantile(x,0.25) # m-sd(x)
#   ymax <- quantile(x,0.75) # m+sd(x)
#   return(data.frame(y=m,ymin=ymin,ymax=ymax))
# }

# violin plots by gene & lab & S dropout with expected marginal means+99% CLs of best fitting median regression model
ctviolinplots_bylab = ggplot(data=ctdata_onlypos_subs_bothgenes, aes(x=factor(S_dropout), y=Ct, fill=factor(S_dropout))) +
  geom_violin(alpha=1, colour=NA, trim=FALSE, draw_quantiles=TRUE, adjust=2, scale="width") +
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
  geom_violin(alpha=1, colour=NA, trim=FALSE, draw_quantiles=TRUE, adjust=2, scale="width") +
  geom_crossbar(data=qr_emmeans99, aes(x=factor(S_dropout), y=emmean, ymin=lower.CL, ymax=upper.CL, group=Gene)) +
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
#                      df      BIC
# fitct_highvirload_0A  4 26262.19
# fitct_highvirload_1A  5 26269.82
# fitct_highvirload_2A  6 26279.71
# fitct_highvirload_3A  6 26279.73
# fitct_highvirload_4A  7 26289.62
# fitct_highvirload_0B  5 26250.33
# fitct_highvirload_1B  6 26257.96
# fitct_highvirload_2B  9 26287.55
# fitct_highvirload_3B  7 26267.88
# fitct_highvirload_4B 10 26297.47


# fitct_highvirload_0B best model
summary(fitct_highvirload_0B) # S dropout samples more frequently have high viral load based on N gene Ct values
# Random effects:
#   Groups     Name        Variance Std.Dev.
# Laboratory (Intercept) 0.02795  0.1672  
# Number of obs: 20834, groups:  Laboratory, 4
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                -0.72877    0.08720  -8.357  < 2e-16 ***
#   GeneORF1ab gene            -0.10455    0.03499  -2.988  0.00281 ** 
#   S_dropout1                  0.44578    0.04583   9.726  < 2e-16 ***
#   GeneORF1ab gene:S_dropout1 -0.30741    0.06596  -4.661 3.15e-06 ***

plot(allEffects(fitct_highvirload_0B))


# odds to encounter high viral load samples based on Ct values of both genes (high vir load = Ct values >1.25x lower than median in non-S dropout samples)
# 1.34x [1.25-1.43x] 95% CLs increased among S dropout samples
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout, type="response"), method="revpairwise", type="response"))
# contrast odds.ratio     SE  df asymp.LCL asymp.UCL
# 1 / 0          1.34 0.0444 Inf      1.25      1.43

# odds ratio = 1.56 [1.43-1.71] for N gene & 1.15 [1.05-1.26] for ORF1ab gene 
confint(contrast(emmeans(fitct_highvirload_0B, ~ S_dropout|Gene, type="response"), method="revpairwise", type="response"))
# contrast odds.ratio    SE  df asymp.LCL asymp.UCL
# 1 / 0          1.15 0.055 Inf      1.04      1.26

fitct_highvirload_emmeans = as.data.frame(emmeans(fitct_highvirload_0B, ~ S_dropout+Gene, type="response"))
fitct_highvirload_emmeans$S_dropout = factor(fitct_highvirload_emmeans$S_dropout)
fitct_highvirload_emmeans$Gene = factor(fitct_highvirload_emmeans$Gene)
fitct_highvirload_emmeans
# S_dropout        Gene      prob         SE  df asymp.LCL asymp.UCL
# 1         0      N gene 0.3254656 0.01914425 Inf 0.2891163 0.3640446
# 2         1      N gene 0.4297220 0.02260911 Inf 0.3860823 0.4744819
# 3         0 ORF1ab gene 0.3029439 0.01844108 Inf 0.2680616 0.3402551
# 4         1 ORF1ab gene 0.3329331 0.02067025 Inf 0.2937238 0.3746005

ct_highvirload_binGLMM = ggplot(data=fitct_highvirload_emmeans, 
                                               aes(x=S_dropout, y=prob*100, fill=S_dropout, group=S_dropout)) +
  facet_wrap(~ Gene) +
  geom_col(colour=NA, position=position_dodge2(width=0.8, padding=0.2)) +
  geom_linerange(aes(ymin=asymp.LCL*100, ymax=asymp.UCL*100), position=position_dodge2(width=0.8, padding=0.2)) +
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

# TO DO: I DIDN'T COMPLETE / UPDATE THE REST BELOW YET ####





# 3. ESTIMATE GROWTH RATE AND TRANSMISSION ADVANTAGE OF B.1.1.7 / 501Y.V1 IN BELGIUM BASED ON S-GENE TARGET FAILURE DATA ####

# We remove ULG - FF 3.x because of low sample size
# Note: we look at last 14 days to minimise impact
#testdata_onlypos = testdata_onlypos[!testdata_onlypos$Laboratory %in% c("UZ Gent","UZA","ULG - FF 3.x"),]
excluded_labs = c("ULG - FF 3.x","ULG")
testdata = testdata[!testdata$Laboratory %in% excluded_labs,]
nrow(testdata) # 281518
testdata_onlypos = testdata_onlypos[!testdata_onlypos$Laboratory %in% excluded_labs,]


# aggregated counts by date (sample date) and Laboratory
data_ag = as.data.frame(table(testdata_onlypos$date, testdata_onlypos$Laboratory, testdata_onlypos$Sdropout), check.names=F)
colnames(data_ag) = c("collection_date", "LABORATORY", "S_DROPOUT", "COUNT")
data_ag_sum = aggregate(COUNT ~ collection_date + LABORATORY, data=data_ag, sum)
data_ag$TOTAL = data_ag_sum$COUNT[match(interaction(data_ag$collection_date,data_ag$LABORATORY),
                                        interaction(data_ag_sum$collection_date,data_ag_sum$LABORATORY))]
data_ag$collection_date = as.Date(data_ag$collection_date)
data_ag$S_DROPOUT = factor(data_ag$S_DROPOUT, levels=c(FALSE,TRUE))
data_ag = data_ag[data_ag$S_DROPOUT==TRUE,]
data_ag$S_DROPOUT = NULL
colnames(data_ag)[which(colnames(data_ag)=="COUNT")] = "S_DROPOUT"
data_ag$LABORATORY = factor(data_ag$LABORATORY)
data_ag$collection_date_num = as.numeric(data_ag$collection_date)
# calculate prop of S dropout that is actually B.1.1.7 / 501Y.V1 estimated from binomial GLMM:
# (using expected marginal mean calculated using emmeans, taking into account random effects)
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                     at=list(collection_date_num=seq(min(data_ag$collection_date_num),
                                                                 max(data_ag$collection_date_num))),
                                     type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")
data_ag$TRUEPOS = fitseq_preds$prob[match(data_ag$collection_date, fitseq_preds$collection_date)] # prob that S dropout was B.1.1.7 / 501Y.V1
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate
data_ag$n_b117 = data_ag$S_DROPOUT*data_ag$TRUEPOS 
data_ag$PROP = data_ag$n_b117/data_ag$TOTAL
data_ag = data_ag[data_ag$TOTAL!=0,]
data_ag$obs = factor(1:nrow(data_ag))
sum(data_ag$TOTAL) == nrow(testdata_onlypos) # TRUE - check
head(data_ag)


# aggregated counts by date over all Laboratories
data_ag_byday = as.data.frame(table(testdata_onlypos$date, testdata_onlypos$Sdropout), check.names=F)
colnames(data_ag_byday) = c("collection_date", "S_DROPOUT", "COUNT")
data_ag_byday_sum = aggregate(COUNT ~ collection_date, data=data_ag_byday, sum)
data_ag_byday$TOTAL = data_ag_byday_sum$COUNT[match(data_ag_byday$collection_date,
                                                    data_ag_byday_sum$collection_date)]
data_ag_byday$collection_date = as.Date(data_ag_byday$collection_date)
data_ag_byday$S_DROPOUT = factor(data_ag_byday$S_DROPOUT, levels=c(FALSE,TRUE))
data_ag_byday = data_ag_byday[data_ag_byday$S_DROPOUT==TRUE,]
data_ag_byday$S_DROPOUT = NULL
colnames(data_ag_byday)[which(colnames(data_ag_byday)=="COUNT")] = "S_DROPOUT"
data_ag_byday$collection_date_num = as.numeric(data_ag_byday$collection_date)
# calculate prop of S dropout that is actually B.1.1.7 / 501Y.V1 estimated from binomial GLMM:
# (using expected marginal mean calculated using emmeans, taking into account random effects)
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) 
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ collection_date_num, 
                                     at=list(collection_date_num=seq(min(data_ag_byday$collection_date_num),
                                                                 max(data_ag_byday$collection_date_num))),
                                     type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$collection_date = as.Date(fitseq_preds$collection_date_num, origin="1970-01-01")
data_ag_byday$TRUEPOS = fitseq_preds$prob[match(data_ag_byday$collection_date, fitseq_preds$collection_date)] # prob that S dropout was B.1.1.7 / 501Y.V1
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate
data_ag_byday$n_b117 = data_ag_byday$S_DROPOUT*data_ag_byday$TRUEPOS 
data_ag_byday$PROP = data_ag_byday$n_b117/data_ag_byday$TOTAL
data_ag_byday = data_ag_byday[data_ag_byday$TOTAL!=0,]
data_ag_byday$obs = factor(1:nrow(data_ag_byday))
sum(data_ag_byday$TOTAL) == nrow(testdata_onlypos) # TRUE - check
head(data_ag_byday)

sum(tail(data_ag_byday$n_b117, 14))/sum(tail(data_ag_byday$TOTAL,14)) 
# 15.4% of the samples of last 2 weeks estimated to be by British variant 
# note: this is not the same as the estimated prop of the new infections or new diagnoses today that are of the British
# variant, which are much higher, see below)


# fit common-slope and separate-slopes binomial GLM
set_sum_contrasts()
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit1)
fit1 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(collection_date_num)+LABORATORY, family=binomial(logit), 
             data=data_ag, control=glmersettings)  # common slope model, with lab coded as fixed factor
fit2 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(collection_date_num)*LABORATORY, family=binomial(logit), 
             data=data_ag, control=glmersettings) # separate slopes model, with lab coded as fixed factor
fit3 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+ns(collection_date_num,df=2)+LABORATORY, family=binomial(logit), 
             data=data_ag, control=glmersettings2)  # common slope model, with lab coded as fixed factor & using 2 df spline ifo date to allow time-varying benefit
fit4 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+ns(collection_date_num,df=2)*LABORATORY, family=binomial(logit), 
             data=data_ag, control=glmersettings) # separate slopes model, with lab coded as fixed factor & using 2 df spline ifo date to allow time-varying benefit
BIC(fit1,fit2,fit3,fit4) 
# df      BIC
# fit1  9 1008.667
# fit2 15 1021.240
# fit3 10 1013.432
# fit4 22 1050.627


# common-slope model fit1 fits best, i.e. rate at which VOC is displacing other strains constant across regions/labs

summary(fit1)


# growth rate advantage (differences in growth rate between VOC and old strains):
# results common-slope model
fit1_emtrends = as.data.frame(emtrends(fit1, revpairwise ~ 1, var="collection_date_num", mode="link", adjust="Tukey")$emtrends)
fit1_emtrends[,c(2,5,6)]
# 0.12 [0.10-0.13] 95% CLs 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*4.7) # 1.74 [1.61-1.87] 95% CLs

# tests for differences in date of introduction
emmeans(fit1,eff~LABORATORY)$contrasts # UCL, ULB, Ghent & UZA earlier than avg (FDR p<0.05), Mons later than avg (FDR p<0.0001)
# contrast                  estimate    SE  df z.ratio p.value
# Namur effect                 0.173 0.167 Inf  1.035  0.3007 
# (Saint LUC - UCL) effect     0.427 0.143 Inf  2.985  0.0099 
# ULB effect                   0.306 0.143 Inf  2.139  0.0454 
# (UMons - Jolimont) effect   -1.918 0.214 Inf -8.973  <.0001 
# UZ Gent effect               0.430 0.153 Inf  2.813  0.0114 
# UZ leuven effect             0.258 0.151 Inf  1.707  0.1024 
# UZA effect                   0.323 0.143 Inf  2.260  0.0417 
# 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: fdr method for 7 tests 


# results spline model fit3 for growth & transmission advantage on the 1st of Febr & the 1st of Jan:
# the growth & transmission advantage as measured today on the 1st of February is probably most
# representative, as for the period between 1-14th of Jan there was quite a bit of active surveillance being done,
# whereas now testing is done more randomly :
fit3_emtrends = as.data.frame(emtrends(fit3, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=as.numeric(as.Date("2021-02-01"))),
                                       mode="link", adjust="Tukey")$emtrends)
fit3_emtrends[,c(2,5,6)]
# 0.10 [0.06-0.14] 95% CLs 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit3_emtrends[,c(2,5,6)]*4.7) # 1.62 [1.33-1.97] 95% CLs

fit3_emtrends = as.data.frame(emtrends(fit3, revpairwise ~ 1, var="collection_date_num", 
                                       at=list(collection_date_num=as.numeric(as.Date("2021-01-01"))),
                                       mode="link", adjust="Tukey")$emtrends)
fit3_emtrends[,c(2,5,6)]
# 0.14 [0.08-0.19] 95% CLs 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit3_emtrends[,c(2,5,6)]*4.7) # 1.89 [1.48-2.42] 95% CLs



# results separate-slopes model fit2:                         
fit2_emtrends = emtrends(fit2, revpairwise ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$emtrends
fit2_emtrends
# LABORATORY       collection_date_num.trend     SE  df asymp.LCL asymp.UCL
# Namur                           0.0854 0.0222 Inf    0.0418     0.129
# Saint LUC - UCL                 0.0818 0.0172 Inf    0.0480     0.116
# ULB                             0.1005 0.0175 Inf    0.0662     0.135
# UMons - Jolimont                0.0755 0.0297 Inf    0.0173     0.134
# UZ Gent                         0.1803 0.0229 Inf    0.1353     0.225
# UZ leuven                       0.1421 0.0208 Inf    0.1014     0.183
# UZA                             0.1389 0.0187 Inf    0.1022     0.175

# only Ghent has a sign abover-average rate of spread, but this could be linked to that lab's heavy focus on active surveillance,
# and so could be due to a bias:
fit2_contrasts = emtrends(fit2, eff ~ LABORATORY, var="collection_date_num", mode="link", adjust="Tukey")$contrasts
fit2_contrasts
# contrast                  estimate     SE  df z.ratio p.value
# Namur effect               -0.0295 0.0204 Inf -1.447  0.6736 
# (Saint LUC - UCL) effect   -0.0331 0.0167 Inf -1.981  0.2889 
# ULB effect                 -0.0144 0.0170 Inf -0.849  0.9707 
# (UMons - Jolimont) effect  -0.0394 0.0264 Inf -1.496  0.6364 
# UZ Gent effect              0.0654 0.0210 Inf  3.118  0.0127 
# UZ leuven effect            0.0272 0.0193 Inf  1.408  0.7027 
# UZA effect                  0.0239 0.0178 Inf  1.345  0.7475 
# 
# P value adjustment: sidak method for 7 tests 



# PLOT MODEL FIT

# for best fitting common slope model fit1
date.to = as.numeric(as.Date("2021-05-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit1_preds = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                                         # by="LABORATORY", 
                                         at=list(collection_date_num=seq(as.numeric(min(data_ag_byday$collection_date)),
                                                                     date.to)), 
                                         type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_preds$collection_date = as.Date(fit1_preds$collection_date_num, origin="1970-01-01")


total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1))$sdcor, function (x) x^2))) 
fit1_preds_bylab = as.data.frame(emmeans(fit1, ~ collection_date_num, 
                                   by="LABORATORY", 
                                   at=list(collection_date_num=seq(as.numeric(min(data_ag_byday$collection_date)),
                                                               date.to)), 
                                    type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_preds_bylab$collection_date = as.Date(fit1_preds_bylab$collection_date_num, origin="1970-01-01")
# order labs by estimated date of introduction (intercepts)
dfemmeanslabs = as.data.frame(emmeans(fit1,~LABORATORY))
levels_BE = as.character(dfemmeanslabs$LABORATORY[order(dfemmeanslabs$emmean,decreasing=T)])
fit1_preds_bylab$LABORATORY = factor(fit1_preds_bylab$LABORATORY, 
                                     levels=levels_BE)





# estimated share of VOC among currently diagnosed infections based on fit1
fit1_preds[fit1_preds$collection_date==as.Date("2021-02-08"),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 27           18659 0.2810843 0.02412842 Inf 0.2360364 0.3303836  2021-02-01
# estimated share of VOC among new infections (assuming time between infection & diagnosis of 7 days)
fit1_preds[fit1_preds$collection_date==(as.Date("2021-02-08")+7),]
#    collection_date_num     prob         SE  df asymp.LCL asymp.UCL collection_date
# 34           18666 0.4519059 0.03992995 Inf 0.3751331 0.5306089  2021-02-08

# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50%, 75% 90% by VOC:
# (really broad confidence intervals though)
fit1_preds$collection_date[fit1_preds[,"prob"]>=0.5][1]-7 # >50% by 3d of February [31 Jan - 7 Febr] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.5][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.5][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.75][1]-7 # >75% by 14th of February [10 Febr - 19 Febr] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.75][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.75][1]-7

fit1_preds$collection_date[fit1_preds[,"prob"]>=0.9][1]-7 # >90% by 23d of Febr [18 Febr - 2 March] 95% CLs
fit1_preds$collection_date[fit1_preds[,"asymp.UCL"]>=0.9][1]-7
fit1_preds$collection_date[fit1_preds[,"asymp.LCL"]>=0.9][1]-7




# PLOT MODEL FIT common-slope model fit1
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
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC 202012/01 BY NHS REGION") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(data_ag$collection_date), as.Date("2021-03-01")-1), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.01,0.99), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag, # data_ag_byday, 
             aes(x=collection_date, y=PROP, size=TOTAL,
                 colour=LABORATORY
                 ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 4), limits=c(10,10^round(log10(max(data_ag_byday$TOTAL)),0)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90))
plot_fit1


saveRDS(plot_fit1, file = paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium by lab.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium by lab.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium by lab.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium by lab.pdf"), width=8, height=6)




# same on response scale (avg over the whole of Belgium):
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
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC 202012/01 BY NHS REGION") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(data_ag_byday$collection_date), as.Date("2021-03-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_byday, 
             aes(x=collection_date, y=100*PROP, size=TOTAL,
                 # colour=LABORATORY
             ), 
             colour=I("steelblue"), 
             alpha=I(1)) +
  scale_size_continuous("number of\npositive tests", trans="sqrt", 
                        range=c(1, 4), limits=c(10,10^round(log10(max(data_ag_byday$TOTAL)),0)), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date")
plot_fit1_response



saveRDS(plot_fit1, file = paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium_response scale.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium_response scale.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\fit1_binomGLMM_VOC_Belgium_response scale.pdf"), width=8, height=6)




# 4. JOINT ANALYSIS OF BELGIAN SGTF DATA WITH UK S GENE DROPOUT (PILLAR 2 SGTF) DATA ####

sgtfdata_uk = read.csv(".//data//uk//sgtf_pillar2_UK-2021-01-25.csv") # Pillar 2 S gene targeted failure data (SGTF) (S dropout)
sgtfdata_uk$other = sgtfdata_uk$other+sgtfdata_uk$sgtf
colnames(sgtfdata_uk) = c("collection_date","REGION","SGTF","TOTAL")
sgtfdata_uk_truepos = read.csv(".//data//uk//sgtf_pillar2_UK-2021-01-25_nick davies_modelled true pos rate sgtfv.csv") # modelled proportion of S dropout that was actually the VOC
# PS this could also be estimated from the COG-UK data based on the presence of deletion 69/70, which is S dropout
sgtfdata_uk$TRUEPOS = sgtfdata_uk_truepos$sgtfv[match(interaction(sgtfdata_uk$REGION, sgtfdata_uk$collection_date),
                                                      interaction(sgtfdata_uk_truepos$group, sgtfdata_uk_truepos$date))] # modelled proportion of S dropout samples that were actually the VOC
sgtfdata_uk$n_b117 = sgtfdata_uk$SGTF*sgtfdata_uk$TRUEPOS
sgtfdata_uk$COUNTRY = "UK"
sgtfdata_uk = sgtfdata_uk[,c("collection_date","COUNTRY","REGION","VOC","TOTAL")]
range(sgtfdata_uk$collection_date) # "2020-10-01" "2021-01-24"
head(sgtfdata_uk)

data_be = data_ag_byday
data_be$REGION = "Belgium"
data_be$COUNTRY = "Belgium"
data_be = data_be[,c("collection_date","COUNTRY","REGION","VOC","TOTAL")]

# joined Belgian S dropout & COG-UK data
data_be_uk2 = rbind(data_be, sgtfdata_uk)
data_be_uk2$COUNTRY = factor(data_be_uk2$COUNTRY)
data_be_uk2$collection_date_num = as.numeric(data_be_uk2$collection_date)
data_be_uk2$PROP = data_be_uk2$n_b117/data_be_uk2$TOTAL
data_be_uk2$obs = factor(1:nrow(data_be_uk2)) # for observation-level random effect, to take into account overdispersion
data_be_uk2$REGION = factor(data_be_uk2$REGION, levels=c(c("Belgium","South East","London","East of England",
                                                           "South West","Midlands","North East and Yorkshire",
                                                           "Scotland","North West","Wales")))
head(data_be_uk2)

set_sum_contrasts()
glmersettings = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e5)) # bobyqa, PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B"))
glmersettings3 = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb"))
glmersettings4 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
fit_be_uk2_1 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(collection_date_num)+REGION, family=binomial(logit), 
                     data=data_be_uk2, control=glmersettings)  # common slope model for country
fit_be_uk2_2 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(collection_date_num)*REGION, family=binomial(logit), 
                     data=data_be_uk2, control=glmersettings) # separate slopes model for country
fit_be_uk2_3 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+ns(collection_date_num,df=3)+REGION, family=binomial(logit), 
                     data=data_be_uk2, control=glmersettings) # with additive spline term
fit_be_uk2_4 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+ns(collection_date_num,df=3)*REGION, family=binomial(logit), 
                     data=data_be_uk2, control=glmersettings3) # with spline term in interaction with region
BIC(fit_be_uk2_1, fit_be_uk2_2, fit_be_uk2_3, fit_be_uk2_4) 
# separate-slopes model very slightly better
# df      BIC
# fit_be_uk2_1 10 7281.447
# fit_be_uk2_2 17 6807.933
# fit_be_uk2_3 12 7178.812
# fit_be_uk2_4 32 6090.350


# model fit_be_uk2_4 best

summary(fit_be_uk2_4)

# PLOT MODEL PREDICTIONS fit_be_uk2_4
# growth rate advantage for BE (differences in growth rate between VOC and old strains):
# results model, with growth rate advantage evaluated today (1/2/2021):
fit_be_uk2_4_emtrends = as.data.frame(emtrends(fit_be_uk2_4, revpairwise ~ 1, 
                                               var="collection_date_num", 
                                               at=list(REGION="Belgium",
                                                       collection_date_num=as.numeric(as.Date("2021-02-01"))),
                                               mode="link", adjust="Tukey")$emtrends)
fit_be_uk2_4_emtrends[,c(2,5,6)]
# 0.096 [0.074-0.12] 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_be_uk2_4_emtrends[,c(2,5,6)]*4.7) 
# 1.57 [1.42-1.73] 95% CLs

# growth & transmission advantage evaluated one month ago (1/1/2021):
fit_be_uk2_4_emtrends = as.data.frame(emtrends(fit_be_uk2_4, revpairwise ~ 1, 
                                               var="collection_date_num", 
                                               at=list(REGION="Belgium",
                                                       collection_date_num=as.numeric(as.Date("2021-01-01"))),
                                               mode="link", adjust="Tukey")$emtrends)
fit_be_uk2_4_emtrends[,c(2,5,6)]
# 0.17 [0.11-0.24] 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_be_uk2_4_emtrends[,c(2,5,6)]*4.7) 
# 2.27 [1.66-3.12] 95% CLs


# for comparison: growth rate advantage for South East UK mid-November (differences in growth rate between VOC and old strains):
# results model :
fit_be_uk2_4_emtrends = as.data.frame(emtrends(fit_be_uk2_4, revpairwise ~ 1, 
                                               var="collection_date_num", 
                                               at=list(REGION="South East",
                                                       collection_date_num=as.numeric(as.Date("2020-11-14"))),
                                               mode="link", adjust="Tukey")$emtrends)
fit_be_uk2_4_emtrends[,c(2,5,6)]
# 0.086 [0.084-0.089] 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_be_uk2_4_emtrends[,c(2,5,6)]*4.7) 
# 1.50 [1.49-1.52] 95% CLs



# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50% or 90% by VOC
# using the joint UK+Belgium model
date.to = as.numeric(as.Date("2021-05-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_be_uk2_4))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_be_uk2_4_preds = as.data.frame(emmeans(fit_be_uk2_4, ~ collection_date_num, 
                                           by=c("REGION"), 
                                           at=list(collection_date_num=seq(as.numeric(min(data_be$collection_date)),
                                                                       date.to),
                                                   COUNTRY="Belgium"), 
                                           type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_be_uk2_4_preds$collection_date = as.Date(fit_be_uk2_4_preds$collection_date_num, origin="1970-01-01")
# fit_be_uk2_4_preds$COUNTRY = factor(fit_be_uk2_4_preds$COUNTRY)

# estimated dates at which new infections with UK variant will reach 50, 75 or 90% (at time of infection, assumed 7 days before diagnosis):
(fit_be_uk2_4_preds$collection_date[fit_be_uk2_4_preds[,"prob"]>=0.5]-7)[1] # >50% by 4th of February [1 Febr - 9 Febr] 95% CLs
(fit_be_uk2_4_preds$collection_date[fit_be_uk2_4_preds[,"asymp.UCL"]>=0.5]-7)[1]
(fit_be_uk2_4_preds$collection_date[fit_be_uk2_4_preds[,"asymp.LCL"]>=0.5]-7)[1]

(fit_be_uk2_4_preds$collection_date[fit_be_uk2_4_preds[,"prob"]>=0.9]-7)[1] # >90% by 27th of February [20th Febr - 11th March] 95% CLs
(fit_be_uk2_4_preds$collection_date[fit_be_uk2_4_preds[,"asymp.UCL"]>=0.9]-7)[1]
(fit_be_uk2_4_preds$collection_date[fit_be_uk2_4_preds[,"asymp.LCL"]>=0.9]-7)[1]


# PLOT MODEL FIT

# spline model fit_be_uk2_4
date.to = as.numeric(as.Date("2021-04-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_be_uk2_4))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_be_uk2_4_preds = as.data.frame(emmeans(fit_be_uk2_4, ~ collection_date_num, 
                                           by=c("REGION"), 
                                           at=list(collection_date_num=seq(as.numeric(min(data_be_uk2$collection_date)),
                                                                       date.to)), 
                                           type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_be_uk2_4_preds$collection_date = as.Date(fit_be_uk2_4_preds$collection_date_num, origin="1970-01-01")
# fit_be_uk2_2_preds$COUNTRY = factor(fit_be_uk2_2_preds$COUNTRY)

n = length(levels(fit_be_uk2_4_preds$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
reg_cols[2:n] = rev(reg_cols[2:n])

levels_UKregions = c("South East","London","East of England",
                     "South West","Midlands","North East and Yorkshire",
                     "Scotland","North West","Wales")

fit_be_uk2_4_preds$REGION = factor(fit_be_uk2_4_preds$REGION, levels=c("Belgium", levels_UKregions))
data_be_uk2$REGION = factor(data_be_uk2$REGION, levels=c("Belgium", levels_UKregions))

# on response scale:
plot_fit_be_uk2_4_response = qplot(data=fit_be_uk2_4_preds, x=collection_date, y=prob*100, geom="blank") +
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
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF 501Y.V1 IN BELGIUM & THE UK") +
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
  geom_point(data=data_be_uk2, 
             aes(x=collection_date, y=PROP*100, size=TOTAL,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 4), limits=c(100,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date")
plot_fit_be_uk2_4_response

saveRDS(plot_fit_be_uk2_4_response, file = paste0(".\\plots\\",dat,"\\plot_fit_be_uk2_4_S dropout data BE and UK_binomial spline GLMM.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\plot_fit_be_uk2_4_S dropout data BE and UK_binomial spline GLMM.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\plot_fit_be_uk2_4_S dropout data BE and UK_binomial spline GLMM.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\plot_fit_be_uk2_4_S dropout data BE and UK_binomial spline GLMM.pdf"), width=8, height=6)




# results separate-slopes per region model:                         
fit_be_uk2_2_emtrends = emtrends(fit_be_uk2_2, ~ REGION, 
                                 var="collection_date_num", 
                                 mode="link")
fit_be_uk2_2_emtrends
# REGION                   collection_date_num.trend       SE  df asymp.LCL asymp.UCL
# Belgium                                 0.1256 0.005374 Inf    0.1150    0.1361
# South East                              0.0888 0.001166 Inf    0.0865    0.0910
# London                                  0.0876 0.000992 Inf    0.0857    0.0896
# East of England                         0.1044 0.001251 Inf    0.1019    0.1068
# South West                              0.0913 0.001097 Inf    0.0891    0.0934
# Midlands                                0.1084 0.001341 Inf    0.1058    0.1110
# North East and Yorkshire                0.0716 0.000950 Inf    0.0698    0.0735
# North West                              0.0914 0.001588 Inf    0.0883    0.0946
# 
# Confidence level used: 0.95 

# significance of differences in slope in Belgium vs in different regions in the UK:
fit_be_uk2_2_contrasts = emtrends(fit_be_uk2_2, trt.vs.ctrl ~ REGION, var="collection_date_num", mode="link", reverse=TRUE)$contrasts
fit_be_uk2_2_contrasts
# contrast                           estimate      SE  df z.ratio p.value
# Belgium - South East                0.02714 0.00800 Inf 3.394   0.0045 
# Belgium - London                    0.02825 0.00797 Inf 3.543   0.0026 
# Belgium - East of England           0.01152 0.00801 Inf 1.438   0.5373 
# Belgium - South West                0.02463 0.00799 Inf 3.083   0.0129 
# Belgium - Midlands                  0.00751 0.00803 Inf 0.935   0.8372 
# Belgium - North East and Yorkshire  0.04428 0.00797 Inf 5.556   <.0001 
# Belgium - North West                0.02445 0.00807 Inf 3.029   0.0153 
# 
# P value adjustment: dunnettx method for 7 tests 




# 5. JOINT ANALYSIS OF BELGIAN S DROPOUT DATA WITH COG-UK SEQUENCING DATA ####
# (NOT INCLUDED IN REPORT)

data_uk = read.csv(".//data//uk//COGUKdata_agbydayregion.csv") 
data_uk = data_uk[data_uk$variant=="VOC 202012/01",]
# COG-UK sequencing data, aggregated by NHS region, from https://github.com/nicholasdavies/newcovid/tree/master/multinomial_logistic_fits/data
head(data_uk)
data_be = data_ag
colnames(data_be)[2] = "REGION"
data_be$COUNTRY = "Belgium"
data_be = data_be[,c("collection_date","COUNTRY","REGION","VOC","TOTAL")]
data_uk$COUNTRY = "UK"
data_uk = data_uk[,c("collection_date","COUNTRY","nhs_name","count","total")]
colnames(data_uk) = c("collection_date","COUNTRY","REGION","VOC","TOTAL")

# joined Belgian S dropout & COG-UK data
data_be_uk = rbind(data_be, data_uk)
data_be_uk$COUNTRY = factor(data_be_uk$COUNTRY)
data_be_uk$collection_date_num = as.numeric(data_be_uk$collection_date)
data_be_uk$PROP = data_be_uk$n_b117/data_be_uk$TOTAL
data_be_uk = data_be_uk[data_be_uk$collection_date>as.Date("2020-08-01"),]
data_be_uk$obs = factor(1:nrow(data_be_uk)) # for observation-level random effect, to take into account overdispersion
head(data_be_uk)

set_sum_contrasts()
glmersettings = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)) # PS : to try all optimizer run all_fit(fit1)
fit_be_uk1 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(collection_date_num)+COUNTRY+REGION, family=binomial(logit), 
             data=data_be_uk, control=glmersettings)  # common slope model for country
fit_be_uk2 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(collection_date_num)*COUNTRY+REGION, family=binomial(logit), 
             data=data_be_uk, control=glmersettings) # separate slopes model for country
fit_be_uk3 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+ns(collection_date_num,df=2)+COUNTRY+REGION, family=binomial(logit), 
                   data=data_be_uk, control=glmersettings)  # with additive 2 df spline
fit_be_uk4 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+ns(collection_date_num,df=2)*COUNTRY+REGION, family=binomial(logit), 
                   data=data_be_uk, control=glmersettings) # with 2 df spline in interaction with country
BIC(fit_be_uk1,fit_be_uk2,fit_be_uk3,fit_be_uk4) 
# common-slope model fits best, i.e. no evidence for the rate of the VOC displacing other variants being different in Belgium vs in the UK
#       df      BIC
# fit_be_uk1 18 2750.357
# fit_be_uk2 19 2742.109
# fit_be_uk3 19 2754.407
# fit_be_uk4 21 2744.210

summary(fit_be_uk1)
summary(fit_be_uk2)

# growth rate advantage (differences in growth rate between VOC and old strains):
# results common-slope model:
fit_be_uk1_emtrends = as.data.frame(emtrends(fit_be_uk4, revpairwise ~ 1, 
                                             var="collection_date_num", 
                                             at=list(collection_date_num=as.numeric(as.Date("2021-02-01"))),
                                             mode="link", adjust="Tukey")$emtrends)
fit_be_uk1_emtrends[,c(2,5,6)]
# 0.09 [0.07-0.11] 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_be_uk1_emtrends[,c(2,5,6)]*4.7) 
# 1.54 [1.39-1.71] 95% CLs

# results separate-slopes per country model:                         
# although one might think there are some slight differences in the growth rate advantage across the UK & Belgium:
fit_be_uk2_emtrends = emtrends(fit_be_uk4, revpairwise ~ COUNTRY, 
                               var="collection_date_num", 
                               at=list(collection_date_num=as.numeric(as.Date("2021-02-01"))),
                               mode="link")$emtrends
fit_be_uk2_emtrends
# COUNTRY collection_date_num.trend     SE  df asymp.LCL asymp.UCL
# Belgium                0.1135 0.0177 Inf    0.0787    0.1482
# UK                     0.0712 0.0132 Inf    0.0453    0.0971
# 
# Confidence level used: 0.95 

# these differences in slope are not actually significant:
fit_be_uk2_contrasts = emtrends(fit_be_uk4, pairwise ~ COUNTRY, var="collection_date_num", mode="link")$contrasts
fit_be_uk2_contrasts
# contrast     estimate      SE  df z.ratio p.value
# Belgium - UK    0.141 0.762 Inf 0.186   0.8528 




# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50% or 90% by VOC
# using the joint UK+Belgium
date.to = as.numeric(as.Date("2021-05-30")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_be_uk1))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_be_uk1_preds = as.data.frame(emmeans(fit_be_uk1, ~ collection_date_num, 
                                         by=c("COUNTRY","REGION"), 
                                         at=list(collection_date_num=seq(min(data_be_uk$collection_date_num),
                                                                     date.to),
                                                 COUNTRY="Belgium"), 
                                         type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_be_uk1_preds$collection_date = as.Date(fit_be_uk1_preds$collection_date_num, origin="1970-01-01")
fit_be_uk1_preds$COUNTRY = factor(fit_be_uk1_preds$COUNTRY)

# estimated dates:
(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"prob"]>=0.5]-7)[1] # >50% by 4th of February [31 Jan - 7 Febr] 95% CLs
(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"asymp.UCL"]>=0.5]-7)[1]
(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"asymp.LCL"]>=0.5]-7)[1]

(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"prob"]>=0.75]-7)[1] # >75% by 15th of February [11 Febr - 19 Febr 95% CLs
(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"asymp.UCL"]>=0.75]-7)[1]
(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"asymp.LCL"]>=0.75]-7)[1]

(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"prob"]>=0.9]-7)[1] # >90% by 26th of February [22 Febr - 2 March] 95% CLs
(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"asymp.UCL"]>=0.9]-7)[1]
(fit_be_uk1_preds$collection_date[fit_be_uk1_preds[,"asymp.LCL"]>=0.9]-7)[1]




# PLOT MODEL FIT

# separate slopes across countries model fit_be_uk2
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_be_uk2))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_be_uk2_preds = as.data.frame(emmeans(fit_be_uk2, ~ collection_date_num, 
                                   by=c("COUNTRY","REGION"), 
                                   at=list(collection_date_num=seq(min(data_be_uk$collection_date_num),
                                                               date.to)), 
                                   type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_be_uk2_preds$collection_date = as.Date(fit_be_uk2_preds$collection_date_num, origin="1970-01-01")
fit_be_uk2_preds$COUNTRY = factor(fit_be_uk2_preds$COUNTRY)

n = length(levels(fit_be_uk2_preds$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
reg_cols[8:n] = rev(reg_cols[8:n])

levels_UKregions = c("South East","London","East of England",
                    "South West","Midlands","North East and Yorkshire",
                    "Scotland","North West","Wales")

fit_be_uk2_preds$REGION = factor(fit_be_uk2_preds$REGION, levels=c(levels_BE, levels_UKregions))
data_be_uk$REGION = factor(data_be_uk$REGION, levels=c(levels_BE, levels_UKregions))

# on response scale:
plot_fit_be_uk2_response = qplot(data=fit_be_uk2_preds, x=collection_date, y=prob*100, geom="blank") +
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
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF 501Y.V1 IN BELGIUM & THE UK") +
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
  geom_point(data=data_be_uk, 
             aes(x=collection_date, y=PROP*100, size=TOTAL,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="sqrt", 
                        range=c(1, 4), limits=c(1,1000), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date")
plot_fit_be_uk2_response


saveRDS(plot_fit_be_uk2_response, file = paste0(".\\plots\\",dat,"\\plot_fit_S dropout BE COGUK sequencing data_binomial spline GLMM.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\plot_fit_S dropout BE COGUK sequencing data_binomial spline GLMM.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\plot_fit_S dropout BE COGUK sequencing data_binomial spline GLMM.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\plot_fit_S dropout BE COGUK sequencing data_binomial spline GLMM.pdf"), width=8, height=6)


# 6. SOME INTERNATIONAL COMPARISONS ####

# 6.1. DATA SWITZERLAND ####

# data from https://ispmbern.github.io/covid-19/variants/

data_geneva = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_GE.csv")
data_geneva$date = as.Date(data_geneva$date)
data_geneva$date_num = as.numeric(data_geneva$date)
data_geneva$obs = factor(1:nrow(data_geneva))
data_zurich = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_ZH.csv")
data_zurich$date = as.Date(data_zurich$date)
data_zurich$date_num = as.numeric(data_zurich$date)
data_zurich$obs = factor(1:nrow(data_zurich))
data_switzerland = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_CH.csv")
data_switzerland$date = as.Date(data_switzerland$date)
data_switzerland$date_num = as.numeric(data_switzerland$date)
data_switzerland$obs = factor(1:nrow(data_switzerland))
sum(data_switzerland[,"B117"]) # 302
sum(data_switzerland[,"total"]) # 12805

fit_geneva = glm(cbind(N501Y,total-N501Y)~date_num, family=quasibinomial(logit), data=data_geneva)
summary(fit_geneva)
cbind(coef(fit_geneva),confint(fit_geneva))[2,] # 0.11 [0.07-0.16] 95% CLs
exp(4.7*cbind(coef(fit_geneva),confint(fit_geneva))[2,]) # 1.67 [1.36-2.08] 
as.data.frame(emmeans(fit_geneva, ~date_num, at=list(date_num=as.numeric(as.Date("2021-02-07"))), type="response"))[,c(2,5,6)] # 71% [55-83%]
fit_zurich = glm(cbind(N501Y,total-N501Y)~date_num, family=quasibinomial(logit), data=data_zurich)
summary(fit_zurich)
cbind(coef(fit_zurich),confint(fit_zurich))[2,] # 0.10 [0.07-0.14] 95% CLs
exp(4.7*cbind(coef(fit_zurich),confint(fit_zurich))[2,]) # 1.61 [1.37-1.92] 
as.data.frame(emmeans(fit_zurich, ~date_num, at=list(date_num=as.numeric(as.Date("2021-02-07"))), type="response"))[,c(2,5,6)] # 40% [27-55%]
fit_switzerland = glm(cbind(B117,total-B117)~date_num, family=quasibinomial(logit), data=data_switzerland)
summary(fit_switzerland)
cbind(coef(fit_switzerland),confint(fit_switzerland))[2,] # 0.11 [0.095-0.14] 95% CLs
exp(4.7*cbind(coef(fit_switzerland),confint(fit_switzerland))[2,]) # 1.71 [1.56-1.88] 
as.data.frame(emmeans(fit_switzerland, ~date_num, at=list(date_num=as.numeric(as.Date("2021-02-07"))), type="response"))[,c(2,5,6)] # 42% [31-56%]
# PS: note that the effect on Rt is in https://ispmbern.github.io/covid-19/variants/
# (1) assumed to be additive as opposed to multiplicative (not quite correct -
# with gamma distributed gen time and small k Rt=exp(r*GT) and multiplicative
# is the logical choice and mult effect on Rt = exp(logistic regression slope*GT);
# (for derivation see https://cmmid.github.io/topics/covid19/uk-novel-variant.html)
# with exponentially distributed GT Rt=r*GT and multiplicative effect on Rt = 1+logistic slope*GT, 
# (2) that unlike in our analysis overdispersion is ignored and that (3) a generation time 
# of 5.2 days instead of 4.7 days is used.


# from randomly selected variants reported by Viollier lab (Basel) & Risch lab 
# (branches throughout Switzerland, e.g. in Zurich & Bern, https://www.risch.ch/de/locations)
# PS : this is the same data as above, just more recent (1 week more for Rish lab data) & split up by lab
data_switzerland2 = read.csv("https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv")
data_switzerland2[is.na(data_switzerland2)] = 0
data_switzerland2$date = as.Date(NA)
data_switzerland2$date[data_switzerland2$week>=51] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_switzerland2$week[data_switzerland2$week>=51] - 1 ) + 1
data_switzerland2$date[data_switzerland2$week<51] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_switzerland2$week[data_switzerland2$week<51] - 1 ) + 6 # PS dates were made to match the ones given in https://ispmbern.github.io/covid-19/variants/data/variants_CH.csv
data_switzerland2$date_num = as.numeric(data_switzerland2$date)
data_switzerland2$obs = factor(1:nrow(data_switzerland2))
data_switzerland2$prop = data_switzerland2$b117/data_switzerland2$n
data_switzerland2$lab = factor(data_switzerland2$lab)
sum(data_switzerland2[-nrow(data_switzerland2),"b117"]) # 302
sum(data_switzerland2[-nrow(data_switzerland2),"n"]) # 12805

summary(fit_switzerland2)
as.data.frame(emtrends(fit_switzerland2, ~ 1, "date_num"))[,c(2,5,6)]
#   date_num.trend  asymp.LCL asymp.UCL
# 1      0.1209993 0.09777233 0.1442262
exp(4.7*as.data.frame(emtrends(fit_switzerland2, ~ 1, "date_num"))[,c(2,5,6)])
#   date_num.trend asymp.LCL asymp.UCL
# 1       1.765964   1.58333  1.969664
as.data.frame(emmeans(fit_switzerland2, ~date_num, 
                      at=list(date_num=as.numeric(as.Date("2021-02-07"))), 
                      type="response"))[,c(2,5,6)] 
# 45% [29-62%]


# 6.2. DATA DENMARK ####

# analysis of data from Denmark, split by region
# from https://www.covid19genomics.dk/statistics
data_denmark = read.csv(".//data/dk//B117_denmark_20210207.csv", sep=";", dec=",")
data_denmark = data_denmark[data_denmark$Region!="Whole Denmark",]
data_denmark$percent = NULL
data_denmark$Region = gsub("SjÃ¦lland","Sjælland",data_denmark$Region)
data_denmark$WEEK = sapply(data_denmark$Week, function(s) as.numeric(strsplit(s, "W")[[1]][[2]]))
data_denmark$date = as.Date(NA)
data_denmark$date[data_denmark$WEEK>=42] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK>=42] - 1 ) + 1
data_denmark$date[data_denmark$WEEK<42] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK<42] - 1 ) + 6 
data_denmark$date_num = as.numeric(data_denmark$date)
data_denmark$obs = factor(1:nrow(data_denmark))
data_denmark$Region = factor(data_denmark$Region)
head(data_denmark)
fit_denmark = glmer(cbind(yes,total-yes) ~ (1|obs) + Region + date_num, family=binomial(logit), data=data_denmark)
summary(fit_denmark)
as.data.frame(emtrends(fit_denmark, ~ 1, "date_num"))[,c(2,5,6)]
#   date_num.trend  asymp.LCL  asymp.UCL
# 1     0.07919027 0.06722759 0.09115295
exp(4.7*as.data.frame(emtrends(fit_denmark, ~ 1, "date_num"))[,c(2,5,6)])
#   date_num.trend asymp.LCL asymp.UCL
# 1       1.450915  1.371589  1.534829
as.data.frame(emmeans(fit_denmark, ~date_num, 
                      at=list(date_num=as.numeric(as.Date("2021-02-07"))), 
                      type="response"))[,c(2,5,6)] 
# 32% [22-44%]

# PS note that https://sites.google.com/site/peterreinhardhansen/research-papers/howcontagiousisthebritishvariantofsars-cov-2
# arrives at a slightly lower growth rate advantage because it uses aggregated counts over the whole of Denmark, and that this
# pools data from regions where the variant may have been introduced at slightly different times
# other analysis for DK is presented at https://ispmbern.github.io/covid-19/variants/

# 6.3. DATA US ####

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
# and estimated increased infectiousness of
exp(4.7*as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.52649  1.328423   1.75409


# FIT FOR WHOLE US + PLOT ####

fitted_truepos = predict(fit_us_propB117amongSGTF, newdat=helix_sgtf, type="response") 
# fitted true positive rate, ie prop of S dropout samples that are B.1.1.7 for dates & states in helix_sgtf

helix_sgtf$estB117 = helix_sgtf$n_sgtf*fitted_truepos # estimated nr of B.1.1.7 samples
helix_sgtf$propB117 = helix_sgtf$estB117/helix_sgtf$n 
fit_us = glmer(cbind(estB117, n-estB117) ~ (1|state/obs)+scale(collection_date_num), 
               family=binomial(logit), data=helix_sgtf) # random intercepts by state
fit_us2 = glmer(cbind(estB117, n-estB117) ~ (collection_date_num||state/obs)+scale(collection_date_num), 
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

# PLOT MODEL FIT
plot_fitus = qplot(data=fit_us_preds2, x=collection_date, y=prob, geom="blank") +
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
  ggtitle("SPREAD OF B.1.1.7 IN THE US") +
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
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 4), limits=c(10,10^(round(log10(max(helix_sgtf$n)),0)+1)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
plot_fitus

saveRDS(plot_fitus, file = paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.pptx"), width=8, height=6)
ggsave(plot_fitus, file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.png"), width=8, height=6)
ggsave(plot_fitus, file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117.pdf"), width=8, height=6)


# PLOT MODEL FIT (response scale)
plot_fitus_resp = qplot(data=fit_us_preds2, x=collection_date, y=prob*100, geom="blank") +
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
  ggtitle("SPREAD OF B.1.1.7 IN THE US") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_calfl2_preds$collection_date), as.Date("2021-04-01")), 
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
                        range=c(1, 4), limits=c(10,10^(round(log10(max(helix_sgtf_subs2$n)),0)+1)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
plot_fitus_resp

saveRDS(plot_fitus_resp, file = paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117_resp.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117_resp.pptx"), width=8, height=6)
ggsave(plot_fitus_resp, file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117_resp.png"), width=8, height=6)
ggsave(plot_fitus_resp, file=paste0(".\\plots\\",dat,"\\fit_us_binomGLMM_B117_resp.pdf"), width=8, height=6)



# FIT JUST USING DATA FROM FLORIDA+CALIFORNIA + PLOTS ####
sel_states = c("FL","CA")
helix_sgtf_subs = helix_sgtf[helix_sgtf$state %in% sel_states,]

fit_us_propB117amongSGTF_calfl1 = glmer(cbind(n_b117, n_sgtf_seq-n_b117) ~ (1|obs)+state+scale(collection_date_num), 
                                 family=binomial(logit), data=helix_b117, subset=helix_b117$state %in% sel_states)
fit_us_propB117amongSGTF_calfl2 = glmer(cbind(n_b117, n_sgtf_seq-n_b117) ~ (1|obs)+state*scale(collection_date_num), 
                                        family=binomial(logit), data=helix_b117, subset=helix_b117$state %in% sel_states)
BIC(fit_us_propB117amongSGTF_calfl1,fit_us_propB117amongSGTF_calfl2) # fit_us_propB117amongSGTF_calfl1 fits best
# implied growth rate advantage of B.1.1.7 over other earlier strains showing S dropout:
as.data.frame(emtrends(fit_us_propB117amongSGTF_calfl1, ~ 1, var="collection_date_num"))[,c(2,5,6)]
#   collection_date_num.trend  asymp.LCL asymp.UCL
# 1                0.07797705 0.04500484 0.1109493

# with a generation time of 4.7 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of
exp(4.7*as.data.frame(emtrends(fit_us_propB117amongSGTF_calfl1, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.442665  1.235558  1.684488

# with a generation time of 5 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of
exp(5*as.data.frame(emtrends(fit_us_propB117amongSGTF_calfl1, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.476811  1.252353  1.741499

# with a generation time of 6.5 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of
exp(6.5*as.data.frame(emtrends(fit_us_propB117amongSGTF_calfl1, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.660055  1.339815  2.056839


fitted_truepos_calfl = predict(fit_us_propB117amongSGTF, newdat=helix_sgtf_subs, type="response") 
helix_sgtf_subs$estB117 = helix_sgtf_subs$n_sgtf*fitted_truepos_calfl # estimated nr of B.1.1.7 samples
helix_sgtf_subs$propB117 = helix_sgtf_subs$estB117/helix_sgtf_subs$n 

fit_calfl1 = glmer(cbind(estB117, n-estB117) ~ (1|obs)+state+scale(collection_date_num), 
               family=binomial(logit), data=helix_sgtf_subs)
fit_calfl2 = glmer(cbind(estB117, n-estB117) ~ (1|obs)+state*scale(collection_date_num), 
                   family=binomial(logit), data=helix_sgtf_subs)
BIC(fit_calfl1,fit_calfl2) # model with different slopes per state fits better
summary(fit_calfl2)

# logistic growth rates (growth rate advantage)
as.data.frame(emtrends(fit_calfl2, ~ state, "collection_date_num"))[,c(1,2,5,6)]
#   state collection_date_num.trend  asymp.LCL  asymp.UCL
# 1    CA                0.06684163 0.05852582 0.07515745
# 2    FL                0.09031168 0.08209082 0.09853253
# CA+FL:
as.data.frame(emtrends(fit_calfl2, ~ 1, "collection_date_num"))[,c(1,2,5,6)]
# 1         collection_date_num.trend  asymp.LCL  asymp.UCL
# 1 overall                0.07857665 0.07272997 0.08442334

# the lower growth rate advantage in CA vs FL is actually significant
# (though it could be caused by competition from other highly contagious strains,
# to test that theory one would have to use a multinomial model or multinomial mixed model
# as in https://cmmid.github.io/topics/covid19/uk-novel-variant.html)
contrast(emtrends(fit_calfl2, ~ state, "collection_date_num"), method="pairwise")
# contrast estimate      SE  df z.ratio p.value
# CA - FL   -0.0235 0.00597 Inf -3.934  0.0001

# increased infectiousness for GT=4.7 days: 
# CA: 37% more infectious [32-42%] 95% CLs
# FL: 53% more infectious [47-59%] 95% CLs
# CA+FL: 45% more infectious [41-49%] 95% CLs
data.frame(state=as.data.frame(emtrends(fit_calfl2, ~ state, "collection_date_num"))[,1],
      exp(4.7*as.data.frame(emtrends(fit_calfl2, ~ state, "collection_date_num"))[,c(2,5,6)]))
# state collection_date_num.trend asymp.LCL asymp.UCL
# 1    CA                  1.369103  1.316625  1.423673
# 2    FL                  1.528772  1.470830  1.588997
exp(4.7*as.data.frame(emtrends(fit_calfl2, ~ 1, "collection_date_num"))[,c(2,5,6)])
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.446736  1.407522  1.487043

# increased infectiousness for GT=5 days: 
# CA: 40% more infectious [34-46%] 95% CLs
# FL: 57% more infectious [51-64%] 95% CLs
# CA+FL: 48% more infectious [44-53%] 95% CLs
data.frame(state=as.data.frame(emtrends(fit_calfl2, ~ state, "collection_date_num"))[,1],
           exp(5*as.data.frame(emtrends(fit_calfl2, ~ state, "collection_date_num"))[,c(2,5,6)]))
# state collection_date_num.trend asymp.LCL asymp.UCL
# 1    CA                  1.396834  1.339946  1.456137
# 2    FL                  1.570758  1.507502  1.636668
exp(5*as.data.frame(emtrends(fit_calfl2, ~ 1, "collection_date_num"))[,c(2,5,6)])
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.481245   1.43857  1.525186

# increased infectiousness for GT=6.5days: 
# CA: 54% more infectious [46-63%] 95% CLs
# FL: 80% more infectious [71-90%] 95% CLs
# CA+FL: 67% more infectious [60-73%] 95% CLs
data.frame(state=as.data.frame(emtrends(fit_calfl2, ~ state, "collection_date_num"))[,1],
           exp(6.5*as.data.frame(emtrends(fit_calfl2, ~ state, "collection_date_num"))[,c(2,5,6)]))
# state collection_date_num.trend asymp.LCL asymp.UCL
# 1    CA                  1.544145  1.462896  1.629908
# 2    FL                  1.798631  1.705043  1.897356
exp(6.5*as.data.frame(emtrends(fit_calfl2, ~ 1, "collection_date_num"))[,c(2,5,6)])
# collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.666538  1.604392  1.731091


# plot model fit fit_calfl2

date.to = as.numeric(as.Date("2021-06-01"))
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_calfl2))$sdcor, function (x) x^2))) 
fit_calfl2_preds = as.data.frame(emmeans(fit_calfl2, ~ collection_date_num, 
                                         by="state", 
                                         at=list(collection_date_num=seq(min(helix_sgtf_subs$collection_date_num),
                                                                     date.to)), 
                                         type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_calfl2_preds$collection_date = as.Date(fit_calfl2_preds$collection_date_num, origin="1970-01-01")
fit_calfl2_preds$state = factor(fit_calfl2_preds$state, 
                                     levels=c("FL","CA"), labels=c("Florida","California"))

# estimated share of B.1.1.7 among currently diagnosed infections based on fit fit_calfl2
fit_calfl2_preds[fit_calfl2_preds$collection_date==as.Date("2021-02-08"),]
#     collection_date_num      state       prob          SE  df  asymp.LCL  asymp.UCL collection_date
# 157               18666 California 0.05116515 0.006268388 Inf 0.04018506 0.06494242      2021-02-08
# 335               18666    Florida 0.16243593 0.013850901 Inf 0.13708060 0.19144070      2021-02-08

# estimated share of B.1.1.7 among new infections (assuming time between infection & diagnosis of 7 days)
fit_calfl2_preds[fit_calfl2_preds$collection_date==(as.Date("2021-02-08")+7),]
#     collection_date_num      state       prob         SE  df  asymp.LCL asymp.UCL collection_date
# 164               18673 California 0.07927163 0.01140055 Inf 0.05961243 0.1046923      2021-02-15
# 342               18673    Florida 0.26736510 0.02522735 Inf 0.22089567 0.3196001      2021-02-15

# taking into account time from infection to diagnosis of ca 7 days this is 
# the time at which new infections would be by more then 50%, 75% 90% by VOC:
# in Florida:
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","prob"]>=0.5][1]-7 # >50% by 20th of February [16 Febr - 24 Febr] 95% CLs
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","asymp.UCL"]>=0.5][1]-7
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","asymp.LCL"]>=0.5][1]-7

fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","prob"]>=0.75][1]-7 # >75% by 4th of March [27 Febr - 9 March] 95% CLs
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","asymp.UCL"]>=0.75][1]-7
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","asymp.LCL"]>=0.75][1]-7

fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","prob"]>=0.9][1]-7 # >90% by 16th of March [11 March - 23 March] 95% CLs
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","asymp.UCL"]>=0.9][1]-7
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="Florida","asymp.LCL"]>=0.9][1]-7

# in California:
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","prob"]>=0.5][1]-7 # >50% by 17th of March [9 March - 27 March] 95% CLs
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","asymp.UCL"]>=0.5][1]-7
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","asymp.LCL"]>=0.5][1]-7

fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","prob"]>=0.75][1]-7 # >75% by 3d of April [24 March - 15 April] 95% CLs
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","asymp.UCL"]>=0.75][1]-7
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","asymp.LCL"]>=0.75][1]-7

fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","prob"]>=0.9][1]-7 # >90% by 19th of April [7 April - 4 May] 95% CLs
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","asymp.UCL"]>=0.9][1]-7
fit_calfl2_preds$collection_date[fit_calfl2_preds[fit_calfl2_preds$state=="California","asymp.LCL"]>=0.9][1]-7


helix_sgtf_subs2 = helix_sgtf_subs
helix_sgtf_subs2$state = factor(helix_sgtf_subs2$state, levels=c("FL","CA"), labels=c("Florida","California"))
helix_sgtf_subs2$collection_date

# PLOT MODEL FIT
plot_fitcafl2 = qplot(data=fit_calfl2_preds, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~state) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
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
  ggtitle("SPREAD OF B.1.1.7 IN THE US") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_calfl2_preds$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.01,0.99), expand=c(0,0)) +
  scale_color_manual("", values=c("red","blue")) +
  scale_fill_manual("", values=c("red","blue")) +
  geom_point(data=helix_sgtf_subs2,  
             aes(x=collection_date, y=propB117, size=n,
                 colour=state
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.3)) +
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 4), limits=c(10,10^(round(log10(max(helix_sgtf_subs2$n)),0)+1)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
plot_fitcafl2

saveRDS(plot_fitcafl2, file = paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117.pptx"), width=8, height=6)
ggsave(plot_fitcafl2, file=paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117.png"), width=8, height=6)
ggsave(plot_fitcafl2, file=paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117.pdf"), width=8, height=6)


# PLOT MODEL FIT (response scale)
plot_fitcafl2_resp = qplot(data=fit_calfl2_preds, x=collection_date, y=prob*100, geom="blank") +
  facet_wrap(~state) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
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
  ggtitle("SPREAD OF B.1.1.7 IN THE US") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_calfl2_preds$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=c("red","blue")) +
  scale_fill_manual("", values=c("red","blue")) +
  geom_point(data=helix_sgtf_subs2,  
             aes(x=collection_date, y=propB117*100, size=n,
                 colour=state
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.3)) +
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 4), limits=c(10,10^(round(log10(max(helix_sgtf_subs2$n)),0)+1)), breaks=c(10,100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
plot_fitcafl2_resp

saveRDS(plot_fitcafl2_resp, file = paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117_resp.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117_resp.pptx"), width=8, height=6)
ggsave(plot_fitcafl2_resp, file=paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117_resp.png"), width=8, height=6)
ggsave(plot_fitcafl2_resp, file=paste0(".\\plots\\",dat,"\\fit_us_cafl_binomGLMM_B117_resp.pdf"), width=8, height=6)

