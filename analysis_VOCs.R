# ANALYSIS OF GROWTH ADVANTAGE OF VOCs B.1.1.7, B.1.351, P.1, B.1.214.2 COMPARED TO OTHER LINEAGES IN BELGIUM BASED ON BASELINE SURVEILLANCE SEQUENCING RESULTS ####

# last update 28 MARCH 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE)
library(emmeans)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-03-28")
today_num = as.numeric(today)
today # "2021-03-28"
plotdir = "VOCs"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))

vocs = read.csv(".//data//be_VOCs//VOCs_byweek.csv")
head(vocs)

# calculate proportions
voc_props = vocs
for (r in 1:length(unique(voc_props$PROVINCE))) {
  rws = (r-1)*6+(1:6)
  voc_props[rws,3:ncol(voc_props)] = sweep(matrix(unlist(voc_props[rws,3:ncol(voc_props)]), ncol=ncol(vocs)-2), 2, 
                                           as.numeric(as.vector(voc_props[rws[1],3:ncol(vocs)])), '/')
}
voc_props[is.na(voc_props)] = 0
voc_props

# convert wide to long
library(tidyr)
vocs_long = gather(vocs, WEEK, COUNT, Week53:Week10, factor_key=FALSE)
vocs_long$WEEK = as.numeric(gsub("Week","",vocs_long$WEEK))
vocs_long$DATE = as.Date(NA)
vocs_long$DATE[vocs_long$WEEK>=42] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( vocs_long$WEEK[vocs_long$WEEK>=42] - 1 ) + 1
vocs_long$DATE[vocs_long$WEEK<42] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( vocs_long$WEEK[vocs_long$WEEK<42] - 1 ) + 6 
vocs_long$DATE_NUM = as.numeric(vocs_long$DATE)
voc_props_long = gather(voc_props, WEEK, PROP, Week53:Week10, factor_key=FALSE)
vocs_long$PROP = voc_props_long$PROP
vocs_long = vocs_long[vocs_long$LINEAGE!="*",]
vocs_long$LINEAGE = factor(vocs_long$LINEAGE, levels=c("Other","B.1.1.7","B.1.351","P.1","B.1.214.2"))
vocs_long = vocs_long[vocs_long$DATE>=as.Date("2021-02-07"),] # we just use data from Feb 7 onward as before there was too much active surveillance
vocs_long$obs = factor(1:nrow(vocs_long))
head(vocs_long)

ggplot(data=vocs_long, aes(x=DATE, 
                           y=COUNT, fill=LINEAGE, group=LINEAGE)) +
  facet_wrap(~PROVINCE) +
  geom_area(aes(fill=LINEAGE), position = position_fill(reverse = FALSE)) +
  theme_hc() +
  scale_fill_manual("LINEAGE", values=c("grey75","red","blue","green3","orange"), 
                    labels=c("wild type","B.1.1.7","B.1.351","P.1", "B.1.214.2")) +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1), 
  #                   expand=c(0,0)) + # limits=as.Date(c("2020-12-31","2021-03-01")
  ylab("Share among newly diagnosed infections") +
  xlab("Date") +
  ggtitle("Spread of the B.1.1.7, B.1.351, P.1 & B.1.214.2\nSARS-CoV2 variants in Belgium (baseline surveillance)") +
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position = "right")
ggsave(file=paste0(".\\plots\\",plotdir,"\\VOCs_muller_plot_raw data.png"), width=10, height=7)


# multinomial fit of share of each variant

set.seed(1)
voc_mfit1 = nnet::multinom(LINEAGE ~ DATE_NUM+PROVINCE, weights=COUNT, data=vocs_long, 
                                     subset=vocs_long$PROVINCE!="All",
                                     maxit=1000) 
voc_mfit2 = nnet::multinom(LINEAGE ~ DATE_NUM*PROVINCE, weights=COUNT, data=vocs_long, 
                           subset=vocs_long$PROVINCE!="All",
                           maxit=1000) 
voc_mfit3 = nnet::multinom(LINEAGE ~ ns(DATE_NUM, df=2)+PROVINCE, weights=COUNT, data=vocs_long, 
                           subset=vocs_long$PROVINCE!="All",
                           maxit=1000) 
voc_mfit4 = nnet::multinom(LINEAGE ~ ns(DATE_NUM, df=2)*PROVINCE, weights=COUNT, data=vocs_long, 
                           subset=vocs_long$PROVINCE!="All",
                           maxit=1000) 
BIC(voc_mfit1,voc_mfit2,voc_mfit3,voc_mfit4) # voc_mfit1 has best BIC
# df      BIC
# voc_mfit1  48 9158.823
# voc_mfit2  88 9740.092
# voc_mfit3  52 9186.459
# voc_mfit4 132 9680.397

# we could also use mblogit to do a multinomial fit taking into account overdispersion
voc_mbfit1 = mblogit(LINEAGE ~ scale(DATE_NUM, center=TRUE, scale=FALSE)+PROVINCE,
                            # random = ~ 1|obs,
                            weights = COUNT, data = vocs_long,
                            subset = vocs_long$PROVINCE!="All",
                            dispersion = TRUE,
                            from.table = TRUE)
voc_mbfit2 = mblogit(LINEAGE ~ scale(DATE_NUM, center=TRUE, scale=FALSE)*PROVINCE,
                     # random = ~ 1|obs,
                     weights = COUNT, data = vocs_long,
                     subset = vocs_long$PROVINCE!="All",
                     dispersion = TRUE,
                     from.table = TRUE)
BIC(voc_mbfit1, voc_mbfit2) # voc_mbfit1 fits best
#       df      BIC
# voc_mbfit1 48 673.8215
# voc_mbfit2 88 902.7102
dispersion(voc_mbfit1, method="Afroz") # dispersion coefficient  = 1.55



# avg growth rate advantage compared to wild type (evaluated today, constant with model mbfit1) ####
delta_r_VOCs = data.frame(confint(emtrends(voc_mbfit1, trt.vs.ctrl ~ LINEAGE|1, 
                                                          var="DATE_NUM",  mode="latent",
                                                          at=list(DATE_NUM=today_num)), 
                                                 adjust="none", df=NA)$contrasts)[,-c(3,4)]
rownames(delta_r_VOCs) = delta_r_VOCs[,"contrast"]
delta_r_VOCs = delta_r_VOCs[,-1]
delta_r_VOCs
#                     estimate  asymp.LCL  asymp.UCL
# B.1.1.7 - Other   0.06468375 0.05600085 0.07336665
# B.1.351 - Other   0.06600375 0.05174988 0.08025763
# P.1 - Other       0.07060070 0.04951780 0.09168361
# B.1.214.2 - Other 0.05875080 0.04219825 0.07530335

# pairwise contrasts in growth rate (here with Tukey correction)
emtrends(voc_mbfit1, revpairwise ~ LINEAGE|1, 
         var="DATE_NUM",  mode="latent",
         at=list(DATE_NUM=today_num), 
         df=NA)$contrasts
# contrast            estimate      SE df z.ratio p.value
# B.1.1.7 - Other      0.06468 0.00443 NA 14.601  <.0001 
# B.1.351 - Other      0.06600 0.00727 NA  9.076  <.0001 
# B.1.351 - B.1.1.7    0.00132 0.00659 NA  0.200  0.9996 
# P.1 - Other          0.07060 0.01076 NA  6.563  <.0001 
# P.1 - B.1.1.7        0.00592 0.01028 NA  0.576  0.9786 
# P.1 - B.1.351        0.00460 0.01194 NA  0.385  0.9954 
# B.1.214.2 - Other    0.05875 0.00845 NA  6.957  <.0001 
# B.1.214.2 - B.1.1.7 -0.00593 0.00782 NA -0.759  0.9423 
# B.1.214.2 - B.1.351 -0.00725 0.00985 NA -0.737  0.9479 
# B.1.214.2 - P.1     -0.01185 0.01244 NA -0.952  0.8762 
# 
# Results are averaged over the levels of: PROVINCE 
# Degrees-of-freedom method: user-specified 
# P value adjustment: tukey method for comparing a family of 5 estimates 


# implied transmission advantage (assuming no immune evasion advantage of B.1.351, P.1 or B.1.214.2
# if there is such an advantage, transm advantage would be less)
exp(delta_r_VOCs*4.7) 
#                   estimate asymp.LCL asymp.UCL
# B.1.1.7 - Other   1.355288  1.301092  1.411740
# B.1.351 - Other   1.363722  1.275355  1.458212
# P.1 - Other       1.393507  1.262045  1.538662
# B.1.214.2 - Other 1.318018  1.219367  1.424649


# plot multinomial model fit

extrapolate = 60
date.from = as.numeric(as.Date("2021-01-01")) # min(vocs_long$DATE_NUM)
date.to =  max(vocs_long$DATE_NUM)+extrapolate

voc_mbfit1_preds = data.frame(emmeans(voc_mbfit1, ~ LINEAGE+PROVINCE+DATE_NUM, 
                                        at=list(DATE_NUM=seq(date.from, date.to)), mode="prob", df=NA))
voc_mbfit1_preds$DATE = as.Date(voc_mbfit1_preds$DATE_NUM, origin="1970-01-01")
voc_mbfit1_preds$LINEAGE = factor(voc_mbfit1_preds$LINEAGE, levels=c("Other","B.1.1.7","B.1.351","P.1","B.1.214.2"),
                                    labels=c("wild type","B.1.1.7","B.1.351","P.1","B.1.214.2"))

vocs_long2 = vocs_long[vocs_long$LINEAGE!="Other",]
vocs_long2$LINEAGE = droplevels(vocs_long2$LINEAGE) 
vocs_long2$LINEAGE = factor(vocs_long2$LINEAGE, levels=c("B.1.1.7","B.1.351","P.1","B.1.214.2"),
                                     labels=c("B.1.1.7","B.1.351","P.1","B.1.214.2"))

ggplot(data=voc_mbfit1_preds, aes(x=DATE, y=prob, group=LINEAGE)) + 
  facet_wrap(~ PROVINCE) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE), position="stack") +
  annotate("rect", xmin=max(vocs_long$DATE)+1, 
                   xmax=as.Date("2021-04-30"), ymin=0, ymax=1, alpha=0.35, fill="white") + # extrapolated part
  scale_fill_manual("LINEAGE", values=c("grey75","red","blue","green3","orange"), 
                    labels=c("wild type","B.1.1.7","B.1.351","P.1", "B.1.214.2")) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01"))),1,1),
                     limits=c(as.Date("2021-01-01"), as.Date("2021-04-30")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share among newly diagnosed infections") +
  ggtitle("Spread of the B.1.1.7, B.1.351, P.1 & B.1.214.2\nSARS-CoV2 variants in Belgium (baseline surveillance)")

ggsave(file=paste0(".\\plots\\",plotdir,"\\VOCs_muller_plot_multinomial fit.png"), width=10, height=7)


# PLOT MODEL FIT WITH DATA & CONFIDENCE INTERVALS

# on response scale:

qplot(data=voc_mbfit1_preds[voc_mbfit1_preds$LINEAGE!="wild type",], x=DATE, y=100*prob, geom="blank") +
  facet_wrap(~PROVINCE) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE
  ), alpha=I(1)) +
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the B.1.1.7, B.1.351, P.1 & B.1.214.2\nSARS-CoV2 variants in Belgium (baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(voc_mbfit1_preds$DATE), as.Date("2021-05-01")),
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")),
                  ylim=c(0,100), expand=c(0,0)) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("LINEAGE", values=c("red","blue","green3","orange")) +
  scale_colour_manual("LINEAGE", values=c("red","blue","green3","orange")) +
  geom_point(data=vocs_long2[vocs_long2$PROVINCE!="All",],
             aes(x=DATE, y=100*PROP, size=COUNT,
                 colour=LINEAGE),
             alpha=I(1)) +
  scale_size_continuous("number of\nsequences", trans="sqrt",
                        range=c(1, 4), limits=c(1,10^3), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Date")

ggsave(file=paste0(".\\plots\\",plotdir,"\\VOCs_muller_plot_multinomial fit_response scale.png"), width=10, height=7)

# saveRDS(plot_multinom_501YV1_501YV2_501YV3_response, file = paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2_501YV3_multinomial fit_model preds_response.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds_response.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds_response.pdf"), width=8, height=6)


# on logit scale:

be_seq_mfit0_preds3 = be_seq_mfit0_preds2
ymin = 0.001
ymax = 0.990001
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
  ylab("Share among newly diagnosed infections (%)") +
  theme_hc() + xlab("") +
  ggtitle("Spread of the British, South African & Brazilian\nSARS-CoV2 variants in Belgium (baseline surveillance)") +
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  scale_fill_manual("variant", values=c("red","blue","green3","black")) +
  scale_colour_manual("variant", values=c("red","blue","green3","black")) +
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
  xlab("Date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-04-01")), ylim=c(0.001, 0.9900001), expand=c(0,0))
plot_multinom_501YV1_501YV2_501YV3


# saveRDS(plot_multinom_501YV1_501YV2_501YV3, file = paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.rds"))
# graph2ppt(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.pptx"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.png"), width=8, height=6)
# ggsave(file=paste0(".\\plots\\",dat,"\\baseline_surveillance_501YV1 501YV2 501YV3_multinomial fit_model preds.pdf"), width=8, height=6)