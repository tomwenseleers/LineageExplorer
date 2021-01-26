# ANALYSIS OF S-GENE TARGET FAILURE DATA FROM BELGIUM TO INFER CONTAGIOUSNESS OF NEW VARIANT OF CONCERN B.1.1.7 / 501Y.V1 ####
# T. Wenseleers, 25 JAN. 2021

library(lme4)
library(emmeans)
library(ggplot2)
library(ggthemes)
# install from https://github.com/tomwenseleers/export
# library(devtools)
# devtools::install_github("tomwenseleers/export")
library(export) 
library(afex)


# 1. ESTIMATE PROPORTION OF S DROPOUT SAMPLES THAT ARE 501Y.V1 IN FUNCTION OF TIME BASED ON SEQUENCING DATA ####
# SEQUENCING DATA FROM EMMANUEL ANDRÉ 25 JAN. 2021

dat_seq = read.csv(".//data//sequencing_Sdropout_Emmanuel Andre.csv", check.names=F)
dat_seq$SAMPLE_DATE = as.Date(dat_seq$SAMPLE_DATE)
dat_seq$SAMPLE_DATE_NUM = as.numeric(dat_seq$SAMPLE_DATE)
dat_seq$PROP_501YV1 = dat_seq$VOC/dat_seq$TOTAL_SDROPOUT_SEQUENCED
dat_seq$obs = factor(1:nrow(dat_seq))
dat_seq

fit_seq = glmer(cbind(VOC,TOTAL_SDROPOUT_SEQUENCED-VOC) ~ (1|obs)+scale(SAMPLE_DATE_NUM), family=binomial(logit), data=dat_seq)
summary(fit_seq)

# PLOT MODEL FIT
extrapolate = 20 # nr of days to extrapolate fit into the future
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_seq))$sdcor, function (x) x^2))) # zero, so not really necessary here 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fitseq_preds = as.data.frame(emmeans(fit_seq, ~ SAMPLE_DATE_NUM, 
                                   at=list(SAMPLE_DATE_NUM=seq(as.numeric(min(dat_seq$SAMPLE_DATE)),
                                                               as.numeric(max(dat_seq$SAMPLE_DATE))+extrapolate)), 
                                   type="response"), bias.adjust = TRUE, sigma = total.SD)
fitseq_preds$SAMPLE_DATE = as.Date(fitseq_preds$SAMPLE_DATE_NUM, origin="1970-01-01")

# prop of S dropout samples among newly diagnosed infections that are now estimated to be 501Y.V1
fitseq_preds[fitseq_preds$SAMPLE_DATE==as.Date("2021-01-26"),]
#    SAMPLE_DATE_NUM      prob         SE  df asymp.LCL asymp.UCL SAMPLE_DATE
# 56           18653 0.9700003 0.02011941 Inf 0.8933488 0.9921018  2021-01-26

# prop of S dropout samples among new infections that are now estimated to be 501Y.V1 (using 7 days for time from infection to diagnosis)
fitseq_preds[fitseq_preds$SAMPLE_DATE==(as.Date("2021-01-26")+7),]
#    SAMPLE_DATE_NUM      prob         SE  df asymp.LCL asymp.UCL SAMPLE_DATE
# 63           18660 0.9878455 0.01042039 Inf 0.9370443 0.9977622  2021-02-02


# on logit scale:
plot_fitseq = qplot(data=fitseq_preds, x=SAMPLE_DATE, y=prob, geom="blank") +
  # facet_wrap(~laboratory) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), fill=I("steelblue"), alpha=I(0.3)) +
  geom_line(aes(y=prob), colour=I("steelblue"), alpha=I(0.8)) +
  ylab("S dropout samples that are 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("REPRESENTATION OF 501Y.V1 AMONG S DROPOUT SAMPLES") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-12-01"),as.Date("2021-02-01")), 
    ylim=c(0.01,0.99001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=dat_seq, 
             aes(x=SAMPLE_DATE, y=PROP_501YV1, size=TOTAL_SDROPOUT_SEQUENCED), colour=I("steelblue"), alpha=I(0.5)) +
  scale_size_continuous("number of S dropout\nsamples sequenced", trans="sqrt", 
                        range=c(0.01, 4), limits=c(5,100), breaks=c(10,50,100)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right") + xlab("Sampling date")
plot_fitseq

saveRDS(plot_fitseq, file = ".\\plots\\representation VOC among S dropout samples_binomial GLMM.rds")
graph2ppt(file=".\\plots\\representation VOC among S dropout samples_binomial GLMM.pptx", width=8, height=6)
ggsave(file=".\\plots\\representation VOC among S dropout samples_binomial GLMM.png", width=8, height=6)
ggsave(file=".\\plots\\representation VOC among S dropout samples_binomial GLMM.pdf", width=8, height=6)



# same on response scale:
plot_fitseq_response = qplot(data=fitseq_preds, x=SAMPLE_DATE, y=100*prob, geom="blank") +
  # facet_wrap(~laboratory) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL), fill=I("steelblue"), alpha=I(0.3)) +
  geom_line(aes(y=100*prob), colour=I("steelblue"), alpha=I(0.8)) +
  ylab("S dropout samples that are 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("REPRESENTATION OF 501Y.V1 AMONG S DROPOUT SAMPLES") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-12-01"),as.Date("2021-02-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=dat_seq, 
             aes(x=SAMPLE_DATE, y=100*PROP_501YV1, size=TOTAL_SDROPOUT_SEQUENCED), colour=I("steelblue"), alpha=I(0.5)) +
  scale_size_continuous("number of S dropout\nsamples sequenced", trans="sqrt", 
                        range=c(0.01, 4), limits=c(5,100), breaks=c(10,50,100)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right") + xlab("Sampling date")
plot_fitseq_response

saveRDS(plot_fitseq_response, file = ".\\plots\\representation VOC among S dropout samples_binomial GLMM_response.rds")
graph2ppt(file=".\\plots\\representation VOC among S dropout samples_binomial GLMM_response.pptx", width=8, height=6)
ggsave(file=".\\plots\\representation VOC among S dropout samples_binomial GLMM_response.png", width=8, height=6)
ggsave(file=".\\plots\\representation VOC among S dropout samples_binomial GLMM_response.pdf", width=8, height=6)




# 2. ESTIMATE GROWTH RATE AND TRANSMISSION ADVANTAGE OF VOC BASED ON S-GENE TARGET FAILURE DATA ####

# Read in testdata
testdata=read.csv(".//data//all valid PCRs since 1st of january.csv")
head(testdata)
testdata$date=as.Date(testdata$Analysis.created.at..UTC.)-1 # sampling date = analysis date-1

# Read in sdropdata
sdropdata=read.csv(".//data//S dropouts raw data for Niel Hens.csv")
sdropdata=sdropdata[,-ncol(sdropdata)]
head(sdropdata)
sdropdata$date=as.Date(sdropdata$Analysis.created.at..UTC.)-1 # sampling date = analysis date-1

# select data from 1st of January onwards
date.from = as.Date("2021-01-01")
date.to = min(max(sdropdata$date),max(testdata$date))
sdropdata = subset(sdropdata,(date<=date.to)&(date>=date.from))
testdata = subset(testdata,(date<=date.to)&(date>=date.from))

# Check sdrop being part of testdata
table(sdropdata$Sample.ID %in% testdata$Sample.ID)

testdata$Sdropout = testdata$Sample.ID %in% sdropdata$Sample.ID
testdata$Laboratory = factor(testdata$Laboratory)
testdata_onlypos = testdata[!(testdata$Outcome=="Not detected"|testdata$Outcome=="Negative"),] # subset with only the positive samples

# we exclude UZ Gent & UZA because they have been heavily involved in pro-active targeted sequencing of suspect infection clusters
# and we remove ULG - FF 3.x because of low sample size
# Note: travellers were screened by ULB & Namur, maybe look only at last 14 days
testdata_onlypos = testdata_onlypos[!testdata_onlypos$Laboratory %in% c("UZ Gent","UZA","ULG - FF 3.x"),]

# aggregated counts by date (sample date) and Laboratory
data_ag = as.data.frame(table(testdata_onlypos$date, testdata_onlypos$Laboratory, testdata_onlypos$Sdropout), check.names=F)
colnames(data_ag) = c("SAMPLE_DATE", "LABORATORY", "S_DROPOUT", "COUNT")
data_ag_sum = aggregate(COUNT ~ SAMPLE_DATE + LABORATORY, data=data_ag, sum)
data_ag$TOTAL = data_ag_sum$COUNT[match(interaction(data_ag$SAMPLE_DATE,data_ag$LABORATORY),
                                        interaction(data_ag_sum$SAMPLE_DATE,data_ag_sum$LABORATORY))]
data_ag$SAMPLE_DATE = as.Date(data_ag$SAMPLE_DATE)
data_ag$S_DROPOUT = factor(data_ag$S_DROPOUT, levels=c(FALSE,TRUE))
data_ag = data_ag[data_ag$S_DROPOUT==TRUE,]
data_ag$S_DROPOUT = NULL
colnames(data_ag)[which(colnames(data_ag)=="COUNT")] = "S_DROPOUT"
data_ag$LABORATORY = factor(data_ag$LABORATORY)
data_ag$SAMPLE_DATE_NUM = as.numeric(data_ag$SAMPLE_DATE)
# prop of S dropout that is actually 501Y.V1 estimated from binomial GLMM:
data_ag$TRUEPOS = predict(fit_seq, newdata=data.frame(SAMPLE_DATE_NUM=data_ag$SAMPLE_DATE_NUM), type="response", re.form=NA) 
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate:
data_ag$VOC = data_ag$S_DROPOUT*data_ag$TRUEPOS 
data_ag$PROP = data_ag$VOC/data_ag$TOTAL
data_ag = data_ag[data_ag$TOTAL!=0,]
data_ag$obs = factor(1:nrow(data_ag))
sum(data_ag$TOTAL) == nrow(testdata_onlypos) # TRUE - check
head(data_ag)

# aggregated counts by date
data_ag_byday = as.data.frame(table(testdata_onlypos$date, testdata_onlypos$Sdropout), check.names=F)
colnames(data_ag_byday) = c("SAMPLE_DATE", "S_DROPOUT", "COUNT")
data_ag_byday_sum = aggregate(COUNT ~ SAMPLE_DATE, data=data_ag_byday, sum)
data_ag_byday$TOTAL = data_ag_byday_sum$COUNT[match(data_ag_byday$SAMPLE_DATE,
                                              data_ag_byday_sum$SAMPLE_DATE)]
data_ag_byday$SAMPLE_DATE = as.Date(data_ag_byday$SAMPLE_DATE)
data_ag_byday$S_DROPOUT = factor(data_ag_byday$S_DROPOUT, levels=c(FALSE,TRUE))
data_ag_byday = data_ag_byday[data_ag_byday$S_DROPOUT==TRUE,]
data_ag_byday$S_DROPOUT = NULL
colnames(data_ag_byday)[which(colnames(data_ag_byday)=="COUNT")] = "S_DROPOUT"
data_ag_byday$SAMPLE_DATE_NUM = as.numeric(data_ag_byday$SAMPLE_DATE)
# prop of S dropout that is actually 501Y.V1 estimated from binomial GLMM:
data_ag_byday$TRUEPOS = predict(fit_seq, newdata=data.frame(SAMPLE_DATE_NUM=data_ag_byday$SAMPLE_DATE_NUM), type="response", re.form=NA) 
# estimated count of 501Y.V1, we adjust numerator of binomial GLMM to take into account true positive rate:
data_ag_byday$VOC = data_ag_byday$S_DROPOUT*data_ag_byday$TRUEPOS 
data_ag_byday$PROP = data_ag_byday$VOC/data_ag_byday$TOTAL
data_ag_byday = data_ag_byday[data_ag_byday$TOTAL!=0,]
data_ag_byday$obs = factor(1:nrow(data_ag_byday))
sum(data_ag_byday$TOTAL) == nrow(testdata_onlypos) # TRUE - check
head(data_ag_byday)


# fit common-slope and separate-slopes binomial GLM
set_sum_contrasts()
fit1 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(SAMPLE_DATE_NUM)+LABORATORY, family=binomial(logit), 
             data=data_ag)  # common slope model, with lab coded as fixed factor
fit2 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(SAMPLE_DATE_NUM)*LABORATORY, family=binomial(logit), 
             data=data_ag) # separate slopes model, with lab coded as fixed factor
BIC(fit1,fit2) # common-slope model fits best, i.e. rate at which VOC is displacing other strains is everywhere equally fast
#      df      BIC
# fit1  7 471.4720
# fit2 11 485.7201

summary(fit1)

# growth rate advantage (differences in growth rate between VOC and old strains):
# results common-slope model:
fit1_emtrends = as.data.frame(emtrends(fit1, revpairwise ~ 1, var="SAMPLE_DATE_NUM", mode="link", adjust="Tukey")$emtrends)
fit1_emtrends[,c(2,5,6)]
# 0.124 [0.0926-0.156] 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit1_emtrends[,c(2,5,6)]*4.7) # 1.79 [1.55-2.08] 95% CLs

# results separate-slopes model:                         
# although one might think there are some slight differences in the growth rate advantage in different regions:
fit2_emtrends = emtrends(fit2, revpairwise ~ LABORATORY, var="SAMPLE_DATE_NUM", mode="link", adjust="Tukey")$emtrends
fit2_emtrends
# LABORATORY      SAMPLE_DATE_NUM.trend     SE  df asymp.LCL asymp.UCL
# Namur                           0.1373 0.0393 Inf    0.0603     0.214
# Saint LUC - UCL                 0.1062 0.0297 Inf    0.0480     0.164
# ULB                             0.1064 0.0296 Inf    0.0483     0.164
# UMons - Jolimont                0.0675 0.0584 Inf   -0.0470     0.182
# UZ leuven                       0.1873 0.0385 Inf    0.1119     0.263

# these differences in slope are not actually significant:
fit2_contrasts = emtrends(fit2, revpairwise ~ LABORATORY, var="SAMPLE_DATE_NUM", mode="link", adjust="Tukey")$contrasts
fit2_contrasts
# contrast                                estimate     SE  df z.ratio p.value
# (Saint LUC - UCL) - Namur              -0.031074 0.0491 Inf -0.632  0.9699 
# ULB - Namur                            -0.030948 0.0494 Inf -0.626  0.9709 
# ULB - (Saint LUC - UCL)                 0.000126 0.0421 Inf  0.003  1.0000 
# (UMons - Jolimont) - Namur             -0.069779 0.0706 Inf -0.989  0.8606 
# (UMons - Jolimont) - (Saint LUC - UCL) -0.038705 0.0656 Inf -0.590  0.9767 
# (UMons - Jolimont) - ULB               -0.038831 0.0654 Inf -0.594  0.9761 
# UZ leuven - Namur                       0.050016 0.0550 Inf  0.910  0.8932 
# UZ leuven - (Saint LUC - UCL)           0.081090 0.0486 Inf  1.669  0.4532 
# UZ leuven - ULB                         0.080964 0.0486 Inf  1.666  0.4553 
# UZ leuven - (UMons - Jolimont)          0.119795 0.0700 Inf  1.711  0.4268 
# 
# P value adjustment: tukey method for comparing a family of 5 estimates 


# PLOT MODEL FIT

# common slope model fit1
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit1))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit1_preds = as.data.frame(emmeans(fit1, ~ SAMPLE_DATE_NUM, 
                                   # by="LABORATORY", 
                                   at=list(SAMPLE_DATE_NUM=seq(as.numeric(min(data$SAMPLE_DATE)),
                                                               date.to)), 
                                    type="response"), bias.adjust = TRUE, sigma = total.SD)
fit1_preds$SAMPLE_DATE = as.Date(fit1_preds$SAMPLE_DATE_NUM, origin="1970-01-01")
# fit1_preds$LABORATORY = factor(fit1_preds$LABORATORY)

# estimated share of VOC among currently diagnosed infections
fit1_preds[fit1_preds$SAMPLE_DATE==as.Date("2021-01-26"),]
#    SAMPLE_DATE_NUM     prob         SE  df asymp.LCL asymp.UCL SAMPLE_DATE
# 27           18653 0.168873 0.03029865 Inf 0.1170844  0.235995  2021-01-26
# estimated share of VOC among new infections (here shifted by one week)
fit1_preds[fit1_preds$SAMPLE_DATE==(as.Date("2021-01-26")+7),]
#    SAMPLE_DATE_NUM     prob         SE  df asymp.LCL asymp.UCL SAMPLE_DATE
# 34           18660 0.313647 0.06469475 Inf 0.2006735 0.4495972  2021-02-02


plot_fit1 = qplot(data=fit1_preds, x=SAMPLE_DATE, y=prob, geom="blank") +
  # facet_wrap(~LABORATORY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  # fill=LABORATORY
                  ), 
              fill=I("steelblue"), 
              alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                # colour=LABORATORY
                ), 
            colour=I("steelblue"), 
            alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC 202012/01 BY NHS REGION") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_byday, 
             aes(x=SAMPLE_DATE, y=PROP, size=TOTAL,
                 # colour=LABORATORY
                 ), 
             colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="log2", 
                        range=c(0.01, 6), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Sampling date")
plot_fit1

saveRDS(plot_fit1, file = ".\\plots\\fit1_binomGLMM_VOC_Belgium.rds")
graph2ppt(file=".\\plots\\fit1_binomGLMM_VOC_Belgium.pptx", width=8, height=6)
ggsave(file=".\\plots\\fit1_binomGLMM_VOC_Belgium.png", width=8, height=6)
ggsave(file=".\\plots\\fit1_binomGLMM_VOC_Belgium.pdf", width=8, height=6)


# same on response scale:
plot_fit1_response = qplot(data=fit1_preds, x=SAMPLE_DATE, y=100*prob, geom="blank") +
  # facet_wrap(~LABORATORY) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL, 
                  # fill=LABORATORY
  ), 
  fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=100*prob, 
                # colour=LABORATORY
  ), 
  colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC 202012/01 BY NHS REGION") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_ag_byday, 
             aes(x=SAMPLE_DATE, y=100*PROP, size=TOTAL,
                 # colour=LABORATORY
             ), 
             colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="log2", 
                        range=c(0.01, 6), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Sampling date")
plot_fit1_response

saveRDS(plot_fit1_response, file = ".\\plots\\fit1_binomGLMM_VOC_Belgium_response scale.rds")
graph2ppt(file=".\\plots\\fit1_binomGLMM_VOC_Belgium_response scale.pptx", width=8, height=6)
ggsave(file=".\\plots\\fit1_binomGLMM_VOC_Belgium_response scale.png", width=8, height=6)
ggsave(file=".\\plots\\fit1_binomGLMM_VOC_Belgium_response scale.pdf", width=8, height=6)




# 3. JOINT ANALYSIS OF BELGIAN SGTF DATA WITH COG-UK SEQUENCING DATA ####

data_uk = read.csv(".//data//COGUKdata_agbydayregion.csv") 
data_uk = data_uk[data_uk$variant=="VOC 202012/01",]
# COG-UK sequencing data, aggregated by NHS region, from https://github.com/nicholasdavies/newcovid/tree/master/multinomial_logistic_fits/data
head(data_uk)
data_be = data_ag
colnames(data_be)[2] = "REGION"
data_be$COUNTRY = "Belgium"
data_be = data_be[,c("SAMPLE_DATE","COUNTRY","REGION","VOC","TOTAL")]
data_uk$COUNTRY = "UK"
data_uk = data_uk[,c("sample_date","COUNTRY","nhs_name","count","total")]
colnames(data_uk) = c("SAMPLE_DATE","COUNTRY","REGION","VOC","TOTAL")

# joined Belgian S dropout & COG-UK data
data_be_uk = rbind(data_be, data_uk)
data_be_uk$COUNTRY = factor(data_be_uk$COUNTRY)
data_be_uk$SAMPLE_DATE_NUM = as.numeric(data_be_uk$SAMPLE_DATE)
data_be_uk$PROP = data_be_uk$VOC/data_be_uk$TOTAL
data_be_uk = data_be_uk[data_be_uk$SAMPLE_DATE>as.Date("2020-08-01"),]
data_be_uk$obs = factor(1:nrow(data_be_uk)) # for observation-level random effect, to take into account overdispersion
head(data_be_uk)

set_sum_contrasts()
fit_be_uk1 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(SAMPLE_DATE_NUM)+COUNTRY+REGION, family=binomial(logit), 
             data=data_be_uk)  # common slope model for country
fit_be_uk2 = glmer(cbind(VOC, TOTAL-VOC) ~ (1|obs)+scale(SAMPLE_DATE_NUM)*COUNTRY+REGION, family=binomial(logit), 
             data=data_be_uk) # separate slopes model for country
BIC(fit_be_uk1,fit_be_uk2) 
# common-slope model fits best, i.e. no evidence for the rate of the VOC displacing other variants being different in Belgium vs in the UK
#            df       BIC
# fit_be_uk1  5  2080.019
# fit_be_uk2  17 2085.841

summary(fit_be_uk1)
summary(fit_be_uk2)

# growth rate advantage (differences in growth rate between VOC and old strains):
# results common-slope model:
fit_be_uk1_emtrends = as.data.frame(emtrends(fit_be_uk1, revpairwise ~ 1, var="SAMPLE_DATE_NUM", mode="link", adjust="Tukey")$emtrends)
fit_be_uk1_emtrends[,c(2,5,6)]
# 0.106 [0.10-0.11] 95% CLs
# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_be_uk1_emtrends[,c(2,5,6)]*4.7) 
# 1.65 [1.60-1.70] 95% CLs

# results separate-slopes per country model:                         
# although one might think there are some slight differences in the growth rate advantage across the UK & Belgium:
fit_be_uk2_emtrends = emtrends(fit_be_uk2, revpairwise ~ COUNTRY, var="SAMPLE_DATE_NUM", mode="link")$emtrends
fit_be_uk2_emtrends
# COUNTRY SAMPLE_DATE_NUM.trend      SE  df asymp.LCL asymp.UCL
# Belgium                 0.125 0.01544 Inf    0.0946     0.155
# UK                      0.106 0.00324 Inf    0.0992     0.112
# 
# Confidence level used: 0.95 

# these differences in slope are not actually significant:
fit_be_uk2_contrasts = emtrends(fit_be_uk2, pairwise ~ COUNTRY, var="SAMPLE_DATE_NUM", mode="link")$contrasts
fit_be_uk2_contrasts
# contrast     estimate      SE  df z.ratio p.value
# Belgium - UK   0.0193  0.0158 Inf 1.225   0.2205 


# PLOT MODEL FIT

# separate slopes across countries model fit_be_uk2
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_be_uk2))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_be_uk2_preds = as.data.frame(emmeans(fit_be_uk2, ~ SAMPLE_DATE_NUM, 
                                   by=c("COUNTRY","REGION"), 
                                   at=list(SAMPLE_DATE_NUM=seq(as.numeric(min(data_be_uk$SAMPLE_DATE)),
                                                               date.to)), 
                                   type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_be_uk2_preds$SAMPLE_DATE = as.Date(fit_be_uk2_preds$SAMPLE_DATE_NUM, origin="1970-01-01")
fit_be_uk2_preds$COUNTRY = factor(fit_be_uk2_preds$COUNTRY)

n = length(levels(fit_be_uk2_preds$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
reg_cols[6:n] = rev(reg_cols[6:n])

levels_UKregions = c("South East","London","East of England",
                    "South West","Midlands","North East and Yorkshire",
                    "Scotland","North West","Wales")
levels_BE = rev(c("UMons - Jolimont","Namur","UZ leuven","ULB","Saint LUC - UCL"))

fit_be_uk2_preds$REGION = factor(fit_be_uk2_preds$REGION, levels=c(levels_BE, levels_UKregions))
data_be_uk$REGION = factor(data_be_uk$REGION, levels=c(levels_BE, levels_UKregions))

plot_fit_be_uk2 = qplot(data=fit_be_uk2_preds, x=SAMPLE_DATE, y=prob, geom="blank") +
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
  ylab("Relative abundance of 501Y.V1 (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF 501Y.V1 IN BELGIUM & THE UK") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0.001,0.999001), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_be_uk, 
             aes(x=SAMPLE_DATE, y=PROP, size=TOTAL,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="log2", 
                        range=c(0.01, 6), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Sampling date")
plot_fit_be_uk2


saveRDS(plot_fit_be_uk2, file = ".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium.rds")
graph2ppt(file=".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium.pptx", width=8, height=6)
ggsave(file=".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium.png", width=8, height=6)
ggsave(file=".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium.pdf", width=8, height=6)


# same on response scale:
plot_fit_be_uk2_response = qplot(data=fit_be_uk2_preds, x=SAMPLE_DATE, y=prob*100, geom="blank") +
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
             aes(x=SAMPLE_DATE, y=PROP*100, size=TOTAL,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests", trans="log2", 
                        range=c(0.01, 6), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Sampling date")
plot_fit_be_uk2_response


saveRDS(plot_fit_be_uk2, file = ".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium_response.rds")
graph2ppt(file=".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium_response.pptx", width=8, height=6)
ggsave(file=".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium_response.png", width=8, height=6)
ggsave(file=".\\plots\\fit_be_uk2_binomGLMM_VOC_Belgium_response.pdf", width=8, height=6)
