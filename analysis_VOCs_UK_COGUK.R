# ANALYSIS OF GROWTH ADVANTAGE OF DIFFERENT VOCs IN THE UK BASED ON COG-UK SEQUENCING DATA ####
# (this excludes data from individuals with known travel history & surge testing/active surveillance)

# last update 4 JUNE 2021

library(nnet)
# devtools::install_github("melff/mclogit",subdir="pkg") # install latest development version of mclogit, to add emmeans support
library(mclogit)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, force = TRUE) # also use latest development version of emmeans for mclogit support
library(emmeans)
library(readr)
library(ggplot2)
library(ggthemes)

today = as.Date(Sys.time()) # we use the file date version as our definition of "today"
today = as.Date("2021-06-04")
today_num = as.numeric(today)
today # "2021-06-04"
plotdir = "VOCs_COGUK"
suppressWarnings(dir.create(paste0(".//plots//",plotdir)))


coguk = read.csv("https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv")
levels_REGION = c("England","Scotland","Wales","Northern Ireland")
coguk$REGION = factor(coguk$adm1, levels=c("UK-ENG","UK-SCT","UK-WLS","UK-NIR"), labels=levels_REGION)
coguk = coguk[!is.na(coguk$REGION),]
levels(coguk$REGION)
coguk$date = as.Date(coguk$sample_date)
head(coguk)

# use data from jan. 1 onwards & Pillar 2 (community as opposed to hospitalised) only
cogukp1 = coguk[coguk$date>=as.Date("2021-01-01")&coguk$is_pillar_2=="N",]
cogukp2 = coguk[coguk$date>=as.Date("2021-01-01")&coguk$is_pillar_2=="Y",]

# B.1.617+ cases before Apr 14 are likely mostly imported cases, so we remove those
cogukp2 = cogukp2[-which(grepl("B.1.617", cogukp2$lineage, fixed=TRUE)&cogukp2$date<=as.Date("2021-04-14")),]  

nrow(cogukp2) # 
range(cogukp2$date) # "2021-01-01" "2021-05-30"
cogukp2$Week = lubridate::week(cogukp2$date)
cogukp2$DATE_NUM = as.numeric(cogukp2$date) 
colnames(cogukp2)

unique(cogukp2$lineage)
sel_target_VOC = "B.1.617"
table(cogukp2$lineage[grepl("B.1.617",cogukp2$lineage, fixed=TRUE)])

unique(cogukp2$lineage[grepl("B.1.617",cogukp2$lineage, fixed=TRUE)])
# "B.1.617.2" "B.1.617.1" "B.1.617.3"
sum(grepl("B.1.617",cogukp2$lineage, fixed=TRUE)) #  B.1.617+
cogukp2$LINEAGE1 = cogukp2$lineage
cogukp2$LINEAGE2 = cogukp2$lineage
cogukp2[grepl(sel_target_VOC, cogukp2$LINEAGE1, fixed=T),"LINEAGE1"] = paste0(sel_target_VOC,"+")
sel_target_VOC = paste0(sel_target_VOC, "+")
cogukp2[grepl("B.1.177", cogukp2$LINEAGE1, fixed=T),"LINEAGE1"] = "B.1.177+"
cogukp2[grepl("B.1.177", cogukp2$LINEAGE2, fixed=T),"LINEAGE2"] = "B.1.177+"

sum("B.1.1.7"==cogukp2$LINEAGE1) #  B.1.1.7


table_lineage = as.data.frame(table(cogukp2$LINEAGE1))
table_lineage$Freq = table_lineage$Freq / sum(table(cogukp2$LINEAGE1))
colnames(table_lineage) = c("lineage","Prop")
table_lineage[table_lineage$Prop>0.001,]
# lineage        Prop
# 74   B.1.1.7 0.882927846
# 85  B.1.177+ 0.027505266
# 111  B.1.351 0.002023195
# 160  B.1.525 0.001810227
# 175 B.1.617+ 0.072756314

sel_ref_lineage = "B.1.1.7"

# sel_lineages = as.character(table_lineage[table_lineage$Prop>0.01,"lineage"][order(table_lineage[table_lineage$Prop>0.01,"Prop"], decreasing=TRUE)])
# sel_lineages = unique(c(sel_lineages, sel_target_VOC, sel_ref_lineage))
sel_lineages = c("B.1.1.7","B.1.177+","B.1.351", "B.1.525","B.1.617+","B.1.617.1","B.1.617.2","P.1")

cogukp2$LINEAGE1[!(cogukp2$LINEAGE1 %in% sel_lineages)] = "other"
cogukp2$LINEAGE2[!(cogukp2$LINEAGE2 %in% sel_lineages)] = "other"
# cogukp2 = cogukp2[cogukp2$LINEAGE1 %in% sel_lineages, ]
sum(table(cogukp2$LINEAGE1))
table(cogukp2$LINEAGE1)
# B.1.1.7 B.1.177+  B.1.351  B.1.525 B.1.617+    other      P.1 
# 190708     5941      437      391    15715     2686      117 
table(cogukp2$LINEAGE2)
# B.1.1.7  B.1.177+   B.1.351   B.1.525 B.1.617.1 B.1.617.2     other       P.1 
# 190708      5941       437       391       345     15363      2693       117 


levels_LINEAGE1 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617+","P.1")
levels_LINEAGE2 = c("other","B.1.1.7","B.1.177+","B.1.351","B.1.525","B.1.617.1","B.1.617.2","P.1")
cogukp2$LINEAGE1 = factor(cogukp2$LINEAGE1, levels=levels_LINEAGE1)
cogukp2$LINEAGE2 = factor(cogukp2$LINEAGE2, levels=levels_LINEAGE2)

cogukp2_lastmonth = cogukp2[cogukp2$date>=(max(cogukp2$date)-30),]
table(cogukp2_lastmonth$lineage, cogukp2_lastmonth$REGION)

str(cogukp2)

# MULTINOMIAL MODEL FIT

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages for whole of UK
data_agbyweek1 = as.data.frame(table(cogukp2$Week, cogukp2$LINEAGE2))
colnames(data_agbyweek1) = c("Week", "LINEAGE2", "count")
data_agbyweek1_sum = aggregate(count ~ Week, data=data_agbyweek1, sum)
data_agbyweek1$total = data_agbyweek1_sum$count[match(data_agbyweek1$Week, data_agbyweek1_sum$Week)]
sum(data_agbyweek1[data_agbyweek1$LINEAGE2=="B.1.617.1","total"]) == nrow(cogukp2) # TRUE
data_agbyweek1$Week = as.numeric(as.character(data_agbyweek1$Week))
data_agbyweek1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweek1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweek1$collection_date_num = as.numeric(data_agbyweek1$collection_date)
data_agbyweek1$prop = data_agbyweek1$count/data_agbyweek1$total

# aggregated by week and region
data_agbyweekregion1 = as.data.frame(table(cogukp2$Week, cogukp2$REGION, cogukp2$LINEAGE2))
colnames(data_agbyweekregion1) = c("Week", "REGION", "LINEAGE2", "count")
data_agbyweekregion1_sum = aggregate(count ~ Week + REGION, data=data_agbyweekregion1, sum)
data_agbyweekregion1$total = data_agbyweekregion1_sum$count[match(interaction(data_agbyweekregion1$Week,data_agbyweekregion1$REGION), 
                                                                  interaction(data_agbyweekregion1_sum$Week,data_agbyweekregion1_sum$REGION))]
sum(data_agbyweekregion1[data_agbyweekregion1$LINEAGE2=="B.1.617.1","total"]) == nrow(cogukp2) # TRUE
data_agbyweekregion1$Week = as.numeric(as.character(data_agbyweekregion1$Week))
data_agbyweekregion1$collection_date = lubridate::ymd( "2021-01-01" ) + lubridate::weeks( data_agbyweekregion1$Week - 1 ) - 3.5 # we use the week midpoint
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)
data_agbyweekregion1$collection_date_num = as.numeric(data_agbyweekregion1$collection_date)
data_agbyweekregion1$prop = data_agbyweekregion1$count/data_agbyweekregion1$total
data_agbyweekregion1 = data_agbyweekregion1[data_agbyweekregion1$total!=0,]

data_agbyweekregion1[data_agbyweekregion1$collection_date==max(data_agbyweekregion1$collection_date),]

# MULLER PLOT (RAW DATA)
unique(cogukp2$LINEAGE2)
levels_LINEAGE2_plot = rev(c("B.1.1.7","B.1.617.2","B.1.617.1","P.1","B.1.351","B.1.525","B.1.177+","other")) # "B.1.617.1","B.1.617.2","B.1.617.3"

library(scales)
n1 = length(levels_LINEAGE2_plot)
lineage_cols1 = c(hcl(h = seq(0, 260, length = n1-1), l = 60, c = 200), muted("dodgerblue", l=55, c=120)) # grey70
# lineage_cols1 = c(hcl(h = seq(0, 265, length = n1), l = 60, c = 200)) # grey70
# lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617+")] = "magenta" # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.1")] = muted("magenta") # muted("magenta",l=50,c=100)
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=50,c=100)

# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.1")] = "magenta" # muted("magenta",l=50,c=100)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.2")] = "magenta" # muted("magenta",l=80,c=255)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.617.3")] = "magenta" # muted("magenta",l=100,c=200)
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.1.7")] = "grey60"  
lineage_cols1[which(levels_LINEAGE2_plot=="other")] = "grey70"  
lineage_cols1[which(levels_LINEAGE2_plot=="B.1.177+")] = "grey55"  
# lineage_cols1[which(levels_LINEAGE1_plot=="B.1.351")] = muted("cyan")  
lineage_cols1[which(levels_LINEAGE2_plot=="P.1")] = "cyan3"  

data_agbyweek1$LINEAGE2 = factor(data_agbyweek1$LINEAGE2, levels=levels_LINEAGE2_plot)
data_agbyweekregion1$LINEAGE2 = factor(data_agbyweekregion1$LINEAGE2, levels=levels_LINEAGE2_plot)
data_agbyweekregion1$REGION = factor(data_agbyweekregion1$REGION, levels=levels_REGION)

library(ggplot2)
library(ggthemes)
muller_cogukp2_raw1 = ggplot(data=data_agbyweek1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  # facet_wrap(~ REGION, ncol=1) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(COG-UK data)") 
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_cogukp2_raw1


muller_cogukp2byregion_raw1 = ggplot(data=data_agbyweekregion1, aes(x=collection_date, y=count, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  # geom_col(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE1), width=1, position="fill") +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="fill") +
  scale_fill_manual("", values=lineage_cols1) +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",NA)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() +
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Share") + 
  theme(legend.position="right",  
        axis.title.x=element_blank()) +
  labs(title = "SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(COG-UK data)") +
  coord_cartesian(xlim=c(as.Date("2021-01-01"),NA)) # as.Date("2021-04-30")
# +
# coord_cartesian(xlim=c(1,max(GISAID_india_bystate$Week)))
muller_cogukp2byregion_raw1

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_raw data.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_raw data.pdf"), width=8, height=6)



# multinomial fits
library(nnet)
library(splines)
data_agbyweekregion1$LINEAGE2 = relevel(data_agbyweekregion1$LINEAGE2, ref="B.1.1.7") # we take B.1.1.7 as baseline / reference level
data_agbyweekregion1$DATE_NUM = data_agbyweekregion1$collection_date_num
set.seed(1)
fit1_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION + DATE_NUM, weights=count, data=data_agbyweekregion1, maxit=1000)
fit2_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION * DATE_NUM, weights=count, data=data_agbyweekregion1, maxit=1000)
fit3_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION + ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion1, maxit=1000)
fit4_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), weights=count, data=data_agbyweekregion1, maxit=1000)
BIC(fit1_cogukp2_multi, fit2_cogukp2_multi, fit3_cogukp2_multi, fit4_cogukp2_multi) 
# fit3_cogukp2_multi fits best (lowest BIC)

# growth rate advantages of different VOCs compared to UK type B.1.1.7 (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2 = emtrends(fit3_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2 = data.frame(confint(emtrcogukp2, 
                                   adjust="none", df=NA)$contrasts, 
                            p.value=as.data.frame(emtrcogukp2$contrasts)$p.value)
delta_r_cogukp2
#               contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1      other - B.1.1.7  0.03744192 0.001586278 NA  0.034332876  0.0405509732 9.059420e-13
# 2 (B.1.177+) - B.1.1.7 -0.12858772 0.010668178 NA -0.149496965 -0.1076784755 9.373613e-13
# 3    B.1.525 - B.1.1.7 -0.00907890 0.004875890 NA -0.018635469  0.0004776692 3.018295e-01
# 4    B.1.351 - B.1.1.7  0.01214762 0.003931767 NA  0.004441499  0.0198537443 2.106755e-02
# 5        P.1 - B.1.1.7  0.01693265 0.008003507 NA  0.001246066  0.0326192372 1.930187e-01
# 6  B.1.617.1 - B.1.1.7 -0.17784291 0.022690257 NA -0.222314994 -0.1333708197 6.544218e-09
# 7  B.1.617.2 - B.1.1.7  0.09615173 0.002478344 NA  0.091294261  0.1010091926 9.059420e-13

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_cogukp2 =  sign(delta_r_cogukp2[,c(2,5,6)])*100*(exp(abs(delta_r_cogukp2[,c(2,5,6)])*4.7)-1)
transmadv_cogukp2 =  data.frame(contrast=delta_r_cogukp2$contrast, transmadv_cogukp2)
transmadv_cogukp2
#            contrast       estimate      asymp.LCL      asymp.UCL
# 1      other - B.1.1.7   19.241069   17.5113242  20.9962750
# 2 (B.1.177+) - B.1.1.7  -83.008477 -101.9067429 -65.8790687
# 3    B.1.525 - B.1.1.7   -4.359432   -9.1536902   0.2247567
# 4    B.1.351 - B.1.1.7    5.875514    2.1094454   9.7804853
# 5        P.1 - B.1.1.7    8.283593    0.5873694  16.5686769
# 6  B.1.617.1 - B.1.1.7 -130.680088 -184.3058681 -87.1692033
# 7  B.1.617.2 - B.1.1.7   57.131542   53.5848545  60.7601321
# so this would estimate that B.1.617.2 had a 57% transmission advantage over B.1.1.7 [53%-61%] 95% CLs


# growth rate & transmission advantages of different VOCs compared to UK type B.1.1.7 by REGION (difference in growth rate per day) 
# PS p values are not Tukey corrected, you can do that with argument adjust="tukey"
emtrcogukp2_region = emtrends(fit4_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2,  
                     var="DATE_NUM",  by=c("REGION"), mode="latent",
                     at=list(DATE_NUM=max(cogukp2$DATE_NUM)))
delta_r_cogukp2_region = data.frame(confint(emtrcogukp2_region, 
                                   adjust="none", df=NA)$contrasts, 
                           p.value=as.data.frame(emtrcogukp2_region$contrasts)$p.value)
delta_r_cogukp2_region
#               contrast     estimate          SE df    asymp.LCL   asymp.UCL      p.value
# 1       other - B.1.1.7          England  0.038024062 1.688448e-03 NA  0.0347147645  0.041333360 1.356012e-10
# 2  (B.1.177+) - B.1.1.7          England -0.110892538 1.170850e-02 NA -0.1338407677 -0.087944308 1.356820e-10
# 3     B.1.525 - B.1.1.7          England -0.006019372 5.148063e-03 NA -0.0161093900  0.004070646 7.107763e-01
# 4     B.1.351 - B.1.1.7          England  0.008983786 4.222485e-03 NA  0.0007078667  0.017259705 1.788600e-01
# 5         P.1 - B.1.1.7          England  0.018842170 8.505776e-03 NA  0.0021711556  0.035513184 1.494189e-01
# 6   B.1.617.1 - B.1.1.7          England -0.167419415 2.315771e-02 NA -0.2128076951 -0.122031136 1.575568e-09
# 7   B.1.617.2 - B.1.1.7          England  0.094398747 2.725433e-03 NA  0.0890569959  0.099740498 1.356012e-10
# 8       other - B.1.1.7         Scotland  0.044229796 5.391076e-03 NA  0.0336634823  0.054796110 1.521652e-10
# 9  (B.1.177+) - B.1.1.7         Scotland -0.174732176 2.771460e-02 NA -0.2290517912 -0.120412561 8.941474e-08
# 10    B.1.525 - B.1.1.7         Scotland -0.046754939 2.121183e-02 NA -0.0883293643 -0.005180513 1.529084e-01
# 11    B.1.351 - B.1.1.7         Scotland  0.041003014 1.411968e-02 NA  0.0133289582  0.068677069 2.802554e-02
# 12        P.1 - B.1.1.7         Scotland  0.012366814 2.450354e-02 NA -0.0356592414  0.060392870 9.734125e-01
# 13  B.1.617.1 - B.1.1.7         Scotland -0.520566206 5.955395e-02 NA -0.6372898000 -0.403842612 1.370150e-10
# 14  B.1.617.2 - B.1.1.7         Scotland  0.097964726 6.746547e-03 NA  0.0847417374  0.111187714 1.356013e-10
# 15      other - B.1.1.7            Wales  0.032555191 1.172131e-02 NA  0.0095818397  0.055528543 3.935412e-02
# 16 (B.1.177+) - B.1.1.7            Wales -0.089139979 5.973362e-02 NA -0.2062157310  0.027935772 5.069965e-01
# 17    B.1.525 - B.1.1.7            Wales -0.072456492 7.681803e-02 NA -0.2230170631  0.078104079 8.335085e-01
# 18    B.1.351 - B.1.1.7            Wales  0.018416995 4.080004e-02 NA -0.0615496046  0.098383594 9.812053e-01
# 19        P.1 - B.1.1.7            Wales -0.196016011 9.502737e-02 NA -0.3822662271 -0.009765794 2.032570e-01
# 20  B.1.617.1 - B.1.1.7            Wales -0.280684720 7.193974e-02 NA -0.4216840210 -0.139685419 1.266795e-03
# 21  B.1.617.2 - B.1.1.7            Wales  0.061463782 2.688081e-02 NA  0.0087783558  0.114149208 1.283115e-01
# 22      other - B.1.1.7 Northern Ireland  0.038942185 2.108355e-02 NA -0.0023808228  0.080265194 3.005923e-01
# 23 (B.1.177+) - B.1.1.7 Northern Ireland -1.900333985 7.922923e-01 NA -3.4531983046 -0.347469665 9.997300e-02
# 24    B.1.525 - B.1.1.7 Northern Ireland -0.057848270 8.897840e-02 NA -0.2322427392  0.116546199 9.425708e-01
# 25    B.1.351 - B.1.1.7 Northern Ireland  0.077403110 3.620245e-10 NA  0.0774031090  0.077403110 1.356012e-10
# 26        P.1 - B.1.1.7 Northern Ireland  0.077403110 4.726242e-10 NA  0.0774031088  0.077403111 1.356012e-10
# 27  B.1.617.1 - B.1.1.7 Northern Ireland  0.077403110          NaN NA           NaN          NaN          NaN
# 28  B.1.617.2 - B.1.1.7 Northern Ireland  0.012256673 3.526433e-02 NA -0.0568601434  0.081373489 9.917618e-01







# CALCULATION OF TRANSMISSION ADVANTAGE THROUGH TIME BY REGION ####

gentime = 4.7 # put the generation time here that you would like to use (e.g. the one typically used to report Re values in the UK)

# growth rate advantages of B.1.617.2 compared to UK type B.1.1.7 through time (difference in growth rate per day) 
# as we would like to get a region & time-varying estimate here we will use model 
# fit4_cogukp2_multi = nnet::multinom(LINEAGE2 ~ REGION * ns(DATE_NUM, df=2), data=cogukp2, maxit=1000) here
emtrcogukp24 = emtrends(fit4_cogukp2_multi, trt.vs.ctrl1 ~ LINEAGE2, by=c("DATE_NUM", "REGION"), 
                      var="DATE_NUM",  mode="latent",
                      at=list(DATE_NUM=seq(as.numeric(as.Date("2021-05-01")),as.numeric(as.Date("2021-06-14")))))
delta_r_cogukp24 = data.frame(confint(emtrcogukp24, 
                                    adjust="none", df=NA)$contrasts)
delta_r_cogukp24 = delta_r_cogukp24[delta_r_cogukp24$contrast=="B.1.617.2 - B.1.1.7",]
delta_r_cogukp24

# If we take the exponent of the product of these growth rate advantages/disadvantages and the generation time (e.g. 4.7 days, Nishiura et al. 2020)
# we get the transmission advantage/disadvantage (here expressed in percent) :

transmadv_cogukp24 = delta_r_cogukp24
transmadv_cogukp24$estimate =  100*(exp(transmadv_cogukp24$estimate*4.7)-1)
transmadv_cogukp24$asymp.LCL =  100*(exp(transmadv_cogukp24$asymp.LCL*4.7)-1)
transmadv_cogukp24$asymp.UCL =  100*(exp(transmadv_cogukp24$asymp.UCL*4.7)-1)
transmadv_cogukp24$collection_date = as.Date(transmadv_cogukp24$DATE_NUM, origin="1970-01-01")
transmadv_cogukp24$REGION = factor(transmadv_cogukp24$REGION, levels=levels_REGION)

plot_cogukp2_mfit4_transmadv = qplot(data=transmadv_cogukp24, 
                               x=collection_date, y=estimate, ymin=asymp.LCL, ymax=asymp.UCL, colour=REGION, fill=REGION, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(fill=REGION, colour=NULL), alpha=I(0.3)) +
  geom_line(aes(colour=REGION, fill=NULL)) +
  ylab("Transmission advantage of B.1.617.2 over B.1.1.7 (%)") +
  theme_hc() + xlab("") +
  ggtitle("TRANSMISSION ADVANTAGE OF B.1.617.2 OVER B.1.1.7 BY REGION\n(COG-UK data, multinomial fit)") +
#  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
#                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
#                     limits=as.Date(c("2021-01-01","2021-05-31")), 
#                     expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  # scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  # scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) # +
  # coord_cartesian(xlim=c(as.Date("2021-05-01"),as.Date("2021-05-31")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_cogukp2_mfit4_transmadv

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit4_transm advantage.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit4_transm advantage.pdf"), width=8, height=6)




# PLOT MULTINOMIAL FIT ####

# extrapolate = 60
date.from = as.numeric(as.Date("2021-01-01"))
date.to = as.numeric(as.Date("2021-06-14")) # max(cogukp2$DATE_NUM)+extrapolate

fit_cogukp2_multi_predsbyregion = data.frame(emmeans(fit3_cogukp2_multi, 
                                                  ~ LINEAGE2,
                                                  by=c("DATE_NUM", "REGION"),
                                                  at=list(DATE_NUM=seq(date.from, date.to)), 
                                                  mode="prob", df=NA))
fit_cogukp2_multi_predsbyregion$collection_date = as.Date(fit_cogukp2_multi_predsbyregion$DATE_NUM, origin="1970-01-01")
fit_cogukp2_multi_predsbyregion$LINEAGE2 = factor(fit_cogukp2_multi_predsbyregion$LINEAGE2, levels=levels_LINEAGE2_plot) 
fit_cogukp2_multi_predsbyregion$REGION = factor(fit_cogukp2_multi_predsbyregion$REGION, levels=levels_REGION) 
# predicted incidence in different parts of the UK today
fit_cogukp2_multi_predsbyregion[fit_cogukp2_multi_predsbyregion$collection_date==today,]

fit_cogukp2_multi_preds = data.frame(emmeans(fit3_cogukp2_multi, 
                                           ~ LINEAGE2,
                                           by=c("DATE_NUM"),
                                           at=list(DATE_NUM=seq(date.from, date.to)), 
                                           mode="prob", df=NA))
fit_cogukp2_multi_preds$collection_date = as.Date(fit_cogukp2_multi_preds$DATE_NUM, origin="1970-01-01")
fit_cogukp2_multi_preds$LINEAGE2 = factor(fit_cogukp2_multi_preds$LINEAGE2, levels=levels_LINEAGE2_plot) 

muller_cogukp2_mfit = ggplot(data=fit_cogukp2_multi_preds, 
                           aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  # facet_wrap(~ STATE, ncol=1) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(cogukp2$DATE_NUM)+1, 
           xmax=as.Date(date.to, origin="1970-01-01"), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01",date.to)), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS B.1.617.1 & B.1.617.2 IN THE UK\n(COG-UK data)")
muller_cogukp2_mfit

library(ggpubr)
ggarrange(muller_cogukp2_raw1+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_cogukp2_mfit+ggtitle("Multinomial fit"), ncol=1)


muller_cogukp2byregion_mfit = ggplot(data=fit_cogukp2_multi_predsbyregion, 
                                  aes(x=collection_date, y=prob, group=LINEAGE2)) + 
  facet_wrap(~ REGION) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=LINEAGE2, group=LINEAGE2), position="stack") +
  scale_fill_manual("", values=lineage_cols1) +
  annotate("rect", xmin=max(cogukp2$DATE_NUM)+1, 
           xmax=as.Date(c("2021-06-14")), ymin=0, ymax=1, alpha=0.3, fill="white") + # extrapolated part
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-14")), expand=c(0,0)) +
  # guides(color = guide_legend(reverse=F, nrow=1, byrow=T), fill = guide_legend(reverse=F, nrow=1, byrow=T)) +
  theme_hc() + theme(legend.position="right", 
                     axis.title.x=element_blank()) + 
  ylab("Share") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit)")
muller_cogukp2byregion_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit.pdf"), width=8, height=6)


ggarrange(muller_cogukp2byregion_raw1+coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date(date.to, origin="1970-01-01")))+
            theme(legend.background = element_rect(fill = alpha("white", 0)),
                  legend.key = element_rect(fill = alpha("white", 0)),
                  legend.text=element_text(color = "white")) +
            guides(colour = guide_legend(override.aes = list(alpha = 0)),
                   fill = guide_legend(override.aes = list(alpha = 0))), 
          muller_cogukp2byregion_mfit+ggtitle("Multinomial fit"), ncol=1)

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_muller plots by region_multinom fit multipanel.pdf"), width=8, height=6)



# PLOT MODEL FIT WITH DATA & 95% CONFIDENCE INTERVALS


# on logit scale:

fit_cogukp2_multi_preds2 = fit_cogukp2_multi_preds
ymin = 0.001
ymax = 0.998
fit_cogukp2_multi_preds2$asymp.LCL[fit_cogukp2_multi_preds2$asymp.LCL<ymin] = ymin
fit_cogukp2_multi_preds2$asymp.UCL[fit_cogukp2_multi_preds2$asymp.UCL<ymin] = ymin
fit_cogukp2_multi_preds2$asymp.UCL[fit_cogukp2_multi_preds2$asymp.UCL>ymax] = ymax
fit_cogukp2_multi_preds2$prob[fit_cogukp2_multi_preds2$prob<ymin] = ymin

fit_cogukp2_multi_predsbyregion2 = fit_cogukp2_multi_predsbyregion
fit_cogukp2_multi_predsbyregion2$asymp.LCL[fit_cogukp2_multi_predsbyregion2$asymp.LCL<ymin] = ymin
fit_cogukp2_multi_predsbyregion2$asymp.UCL[fit_cogukp2_multi_predsbyregion2$asymp.UCL<ymin] = ymin
fit_cogukp2_multi_predsbyregion2$asymp.UCL[fit_cogukp2_multi_predsbyregion2$asymp.UCL>ymax] = ymax
fit_cogukp2_multi_predsbyregion2$prob[fit_cogukp2_multi_predsbyregion2$prob<ymin] = ymin

# plots of multinomial fit to cogukp2 Inst data 

plot_cogukp2_mfit_logit = qplot(data=fit_cogukp2_multi_predsbyregion2, 
                               x=collection_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-14")), 
                     expand=c(0,0)) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(1), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")+
  coord_cartesian(xlim=c(as.Date("2021-01-01"),as.Date("2021-06-14")), ylim=c(0.001, 0.998), expand=c(0,0))
plot_cogukp2_mfit_logit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.pdf"), width=8, height=6)
library(svglite)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_logit scale.svg"), width=8, height=6)

# on response scale:
plot_cogukp2_mfit = qplot(data=fit_cogukp2_multi_predsbyregion2, 
                         x=collection_date, y=100*prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=100*prob, ymin=100*asymp.LCL, ymax=100*asymp.UCL, colour=NULL,
                  fill=LINEAGE2
  ), alpha=I(0.3)) +
  geom_line(aes(y=100*prob,
                colour=LINEAGE2
  ), alpha=I(1)) +
  ylab("Share (%)") +
  theme_hc() + xlab("") +
  ggtitle("SPREAD OF SARS-CoV2 VARIANTS OF CONCERN IN THE UK\n(COG-UK data, multinomial fit)") +
  scale_x_continuous(breaks=as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01")),
                     labels=substring(months(as.Date(c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01","2021-06-01"))),1,1),
                     limits=as.Date(c("2021-01-01","2021-06-14")), 
                     expand=c(0,0)) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                     labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=as.Date(c("2021-01-01",date.to-1)),
                  ylim=c(0,100), expand=c(0,0)) +
  scale_fill_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  scale_colour_manual("variant", values=c(lineage_cols1,"darkgreen")) +
  geom_point(data=data_agbyweekregion1,
             aes(x=collection_date, y=prop*100, size=total,
                 colour=LINEAGE2, 
             ),
             alpha=I(1), pch=I(16)) +
  scale_size_continuous("total n", trans="sqrt",
                        range=c(1, 6), limits=c(1,10^5), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) +
  # guides(colour=FALSE) +
  theme(legend.position = "right") +
  xlab("Collection date")
plot_cogukp2_mfit

ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_response scale.png"), width=8, height=6)
ggsave(file=paste0(".\\plots\\",plotdir,"\\cogukp2_multinom fit_response scale.pdf"), width=8, height=6)



