# REPRODUCIBLE EXAMPLE nnet:multinom, mclogit::mblogit, VGAM::vglm, emmeans & marginaleffects
# dat = data_agbyweek1
# dat$variant2 = NULL
# write.csv(dat, "C://Users//bherr//Dropbox//temp//dat.csv", row.names=F)

# data=SARS-CoV2 coronavirus variants (variant) through time (collection_date_num)
# in India, count=actual count (nr of sequenced genomes)

dat = read.csv("https://www.dropbox.com/s/u27cn44p5srievq/dat.csv?dl=1")
dat$variant = factor(dat$variant)
# set.seed(1)
# dat$count = dat$count+runif(nrow(dat), 0, 1E-3)

# 1. multinom::net multinomial fit ####
library(nnet)
library(splines)
set.seed(1)
fit_nnet = nnet::multinom(variant ~ ns(collection_date_num, df=2), 
                          weights=count, data=dat)
summary(fit_nnet)
# Coefficients:
#   (Intercept) ns(collection_date_num, df = 2)1
# Beta                 3.798225                        4.9204017
# Delta              -10.043973                       49.4947796
# Omicron (BA.1)    -265.410841                      487.1385690
# Omicron (BA.2)    -467.537309                      831.1140917
# Omicron (BA.2.74)   -4.616715                       15.4357149
# Omicron (BA.2.75)   11.365092                      -38.7658668
# Omicron (BA.2.76)    5.437985                        0.3403524
# Omicron (BA.3)     -15.924342                       42.5896877
# Omicron (BA.4)      -7.440555                       23.8097760
# Omicron (BA.5)      -2.940342                       15.3461333
# Other               19.322779                       -8.4319599
# ns(collection_date_num, df = 2)2
# Beta                                      35.17926
# Delta                                     55.21257
# Omicron (BA.1)                           195.27695
# Omicron (BA.2)                           312.90170
# Omicron (BA.2.74)                         75.53290
# Omicron (BA.2.75)                         79.80443
# Omicron (BA.2.76)                         70.73509
# Omicron (BA.3)                            75.60153
# Omicron (BA.4)                            73.26847
# Omicron (BA.5)                            74.45010
# Other                                     54.14069
# 
# Std. Errors:
#   (Intercept) ns(collection_date_num, df = 2)1
# Beta                1.8548783                        1.9750676
# Delta               0.8914081                        0.7919794
# Omicron (BA.1)      8.9955175                       15.6207548
# Omicron (BA.2)      7.0504753                       12.1069603
# Omicron (BA.2.74)   1.1513095                        0.7721499
# Omicron (BA.2.75)   2.5900436                        5.7317343
# Omicron (BA.2.76)  15.9802663                       27.2148132
# Omicron (BA.3)     27.1955306                       46.8288796
# Omicron (BA.4)      1.0868226                        0.6370319
# Omicron (BA.5)      1.3766617                        1.5280076
# Other               0.8118766                        0.6288743
# ns(collection_date_num, df = 2)2
# Beta                                      5.998582
# Delta                                     3.133279
# Omicron (BA.1)                            5.388977
# Omicron (BA.2)                            4.888832
# Omicron (BA.2.74)                         3.177189
# Omicron (BA.2.75)                         3.862689
# Omicron (BA.2.76)                         9.112225
# Omicron (BA.3)                           14.632091
# Omicron (BA.4)                            3.116746
# Omicron (BA.5)                            3.081921
# Other                                     3.128327
# 
# Residual Deviance: 128392.2 
# AIC: 128458.2 


# emmeans & emtrends using emmeans package: works, but slow
# emmeans = current proportion of different variants
library(emmeans)
multinom_emmeans = emmeans(fit_nnet, ~ variant,  
                       mode = "prob",
                       at=list(collection_date_num = 
                                 max(data_agbyweek1$collection_date_num)))
multinom_emmeans
# variant               prob       SE df lower.CL upper.CL
# Alpha             0.00e+00 0.00e+00 33 0.00e+00 0.00e+00
# Beta              0.00e+00 0.00e+00 33 0.00e+00 0.00e+00
# Delta             7.73e-06 1.17e-06 33 5.34e-06 1.01e-05
# Omicron (BA.1)    1.82e-04 6.42e-05 33 5.14e-05 3.13e-04
# Omicron (BA.2)    1.76e-01 7.45e-03 33 1.61e-01 1.91e-01
# Omicron (BA.2.74) 9.03e-02 7.98e-03 33 7.41e-02 1.07e-01
# Omicron (BA.2.75) 1.68e-01 1.90e-02 33 1.30e-01 2.07e-01
# Omicron (BA.2.76) 2.89e-01 1.35e-02 33 2.62e-01 3.16e-01
# Omicron (BA.3)    1.34e-02 2.10e-03 33 9.10e-03 1.76e-02
# Omicron (BA.4)    1.67e-02 2.47e-03 33 1.17e-02 2.17e-02
# Omicron (BA.5)    2.03e-01 1.08e-02 33 1.81e-01 2.25e-01
# Other             4.23e-02 3.15e-03 33 3.59e-02 4.87e-02
# 
# Confidence level used: 0.95 

# pairwise contrasts in marginal trends on link scale
# = current differences in growth rate among variants
multinom_emtrends = confint(pairs(emtrends(fit_nnet, ~ variant,  
                                       var = "collection_date_num",  
                                       mode = "latent",
                                       at = list(collection_date_num = 
                                                 max(data_agbyweek1$collection_date_num))),
                              reverse=TRUE))
multinom_emtrends

# equivalent using margineffects package 
# does not work as there is no predict.multinom method with type="link"
library(marginaleffects)
multinom_preds_marginaleffects = predictions(fit_nnet,
                                         newdata = datagrid(collection_date_num = 
                                                              max(data_agbyweek1$collection_date_num)),
                                         type="link", # not supported by predict.multinom
                                         transform_post = insight::link_inverse(fit_nnet))

# Error: The `type` argument for models of class `multinom` must be an element of: probs
multinom_preds_marginaleffects 

# NOTE: would it be possible perhaps to redefine the predict.multinom function to use the
# predict.mblogit method instead, which might work after small modification,
# and which does support type="link":
# https://stackoverflow.com/questions/73010776/redefining-rs-nnetmultinom-predict-multinom-predict-method-to-support-type-l

# marginal effects
# same problem for marginal effects
multinom_marginaleffects = marginaleffects(fit_nnet, newdata = datagrid(collection_date_num = 
                                             max(data_agbyweek1$collection_date_num)),
                variables="collection_date_num",
                # type="probs"
                type="link" # not supported by predict.multinom
                )
# Error: The `type` argument for models of class `multinom` must be an element of: probs

# NOTE: I was not sure if pairwise contrasts in marginal effects 
# as in emtrends example above are supported by marginaleffects & of syntax,
# comparisons(multinom_marginaleffects) e.g. does not work


# 2. mclogit::mblogit multinomial fit ####
# this one does have a predict.mblogit method with type="link",
# https://github.com/melff/mclogit/blob/master/pkg/R/mblogit.R
library(mclogit)
fit_mblogit = mblogit(variant ~ ns(collection_date_num, df=2),
                      weight=count,
                      data=dat,
                      from.table=TRUE, dispersion=FALSE,
                      control=mclogit.control(maxit=27))
summary(fit_mblogit)
# Equation for Beta vs Alpha:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                      -2.140e+01  1.505e+06       0        1
# ns(collection_date_num, df = 2)1 -9.047e+00  3.062e+06       0        1
# ns(collection_date_num, df = 2)2  2.608e+01  1.379e+06       0        1
# 
# Equation for Delta vs Alpha:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                      -2.140e+01  1.505e+06       0        1
# ns(collection_date_num, df = 2)1 -9.047e+00  3.062e+06       0        1
# ns(collection_date_num, df = 2)2  2.608e+01  1.379e+06       0        1
# 
# Equation for Omicron (BA.1) vs Alpha:
#   Estimate Std. Error z value
# Omicron (BA.1)~(Intercept)                      -2.140e+01  1.505e+06       0
# Omicron (BA.1)~ns(collection_date_num, df = 2)1 -9.047e+00  3.062e+06       0
# Omicron (BA.1)~ns(collection_date_num, df = 2)2  2.608e+01  1.379e+06       0
# Pr(>|z|)
# Omicron (BA.1)~(Intercept)                             1
# Omicron (BA.1)~ns(collection_date_num, df = 2)1        1
# Omicron (BA.1)~ns(collection_date_num, df = 2)2        1
# 
# Equation for Omicron (BA.2) vs Alpha:
#   Estimate Std. Error z value
# Omicron (BA.2)~(Intercept)                      -6.051e+00  5.785e+05       0
# Omicron (BA.2)~ns(collection_date_num, df = 2)1 -3.120e+01  1.257e+06       0
# Omicron (BA.2)~ns(collection_date_num, df = 2)2  5.841e+01  3.135e+05       0
# Pr(>|z|)
# Omicron (BA.2)~(Intercept)                             1
# Omicron (BA.2)~ns(collection_date_num, df = 2)1        1
# Omicron (BA.2)~ns(collection_date_num, df = 2)2        1
# 
# Equation for Omicron (BA.2.74) vs Alpha:
#   Estimate Std. Error
# Omicron (BA.2.74)~(Intercept)                      -7.065e+00  5.791e+05
# Omicron (BA.2.74)~ns(collection_date_num, df = 2)1 -2.952e+01  1.258e+06
# Omicron (BA.2.74)~ns(collection_date_num, df = 2)2  5.573e+01  3.116e+05
# z value Pr(>|z|)
# Omicron (BA.2.74)~(Intercept)                            0        1
# Omicron (BA.2.74)~ns(collection_date_num, df = 2)1       0        1
# Omicron (BA.2.74)~ns(collection_date_num, df = 2)2       0        1
# 
# Equation for Omicron (BA.2.75) vs Alpha:
#   Estimate Std. Error
# Omicron (BA.2.75)~(Intercept)                      -6.311e+00  5.785e+05
# Omicron (BA.2.75)~ns(collection_date_num, df = 2)1 -3.077e+01  1.257e+06
# Omicron (BA.2.75)~ns(collection_date_num, df = 2)2  5.771e+01  3.131e+05
# z value Pr(>|z|)
# Omicron (BA.2.75)~(Intercept)                            0        1
# Omicron (BA.2.75)~ns(collection_date_num, df = 2)1       0        1
# Omicron (BA.2.75)~ns(collection_date_num, df = 2)2       0        1
# 
# Equation for Omicron (BA.2.76) vs Alpha:
#   Estimate Std. Error
# Omicron (BA.2.76)~(Intercept)                      -6.072e+00  5.784e+05
# Omicron (BA.2.76)~ns(collection_date_num, df = 2)1 -3.116e+01  1.257e+06
# Omicron (BA.2.76)~ns(collection_date_num, df = 2)2  5.835e+01  3.137e+05
# z value Pr(>|z|)
# Omicron (BA.2.76)~(Intercept)                            0        1
# Omicron (BA.2.76)~ns(collection_date_num, df = 2)1       0        1
# Omicron (BA.2.76)~ns(collection_date_num, df = 2)2       0        1
# 
# Equation for Omicron (BA.3) vs Alpha:
#   Estimate Std. Error z value
# Omicron (BA.3)~(Intercept)                      -2.140e+01  1.505e+06       0
# Omicron (BA.3)~ns(collection_date_num, df = 2)1 -9.047e+00  3.062e+06       0
# Omicron (BA.3)~ns(collection_date_num, df = 2)2  2.608e+01  1.379e+06       0
# Pr(>|z|)
# Omicron (BA.3)~(Intercept)                             1
# Omicron (BA.3)~ns(collection_date_num, df = 2)1        1
# Omicron (BA.3)~ns(collection_date_num, df = 2)2        1
# 
# Equation for Omicron (BA.4) vs Alpha:
#   Estimate Std. Error z value
# Omicron (BA.4)~(Intercept)                      -7.442e+00  5.800e+05       0
# Omicron (BA.4)~ns(collection_date_num, df = 2)1 -2.890e+01  1.260e+06       0
# Omicron (BA.4)~ns(collection_date_num, df = 2)2  5.475e+01  3.110e+05       0
# Pr(>|z|)
# Omicron (BA.4)~(Intercept)                             1
# Omicron (BA.4)~ns(collection_date_num, df = 2)1        1
# Omicron (BA.4)~ns(collection_date_num, df = 2)2        1
# 
# Equation for Omicron (BA.5) vs Alpha:
#   Estimate Std. Error z value
# Omicron (BA.5)~(Intercept)                      -6.825e+00  5.785e+05       0
# Omicron (BA.5)~ns(collection_date_num, df = 2)1 -2.992e+01  1.257e+06       0
# Omicron (BA.5)~ns(collection_date_num, df = 2)2  5.635e+01  3.120e+05       0
# Pr(>|z|)
# Omicron (BA.5)~(Intercept)                             1
# Omicron (BA.5)~ns(collection_date_num, df = 2)1        1
# Omicron (BA.5)~ns(collection_date_num, df = 2)2        1
# 
# Equation for Other vs Alpha:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                      -2.140e+01  1.505e+06       0        1
# ns(collection_date_num, df = 2)1 -9.047e+00  3.062e+06       0        1
# ns(collection_date_num, df = 2)2  2.608e+01  1.379e+06       0        1
# 
# Null Deviance:     13540 
# Residual Deviance: 7.354e-10 
# Number of Fisher Scoring iterations:  27 
# Number of observations:  43160 
# 
# 
# Note: Algorithm did not converge.

# PS increasing maxit to >30 results in 
# Error in solve.default(Information) : 
#  system is computationally singular: reciprocal condition number = 1.66419e-16
# maybe due to complete separation / Hauck-Donner effect or bug in mblogit??
# Right now SEs are very large, probably bug in mblogit or due to H-D effect

# emmeans & emtrends using emmeans package: works, but slow
# emmeans = current proportion of different variants
mblogit_emmeans = emmeans(fit_mblogit, ~ variant,  
                           mode = "prob",
                           at=list(collection_date_num = 
                                     max(data_agbyweek1$collection_date_num)))
mblogit_emmeans
# variant             prob       SE  df asymp.LCL asymp.UCL
# Alpha             0.0000 9.00e-08 Inf -1.80e-07  2.00e-07
# Beta              0.0000 1.00e-08 Inf -2.00e-08  0.00e+00
# Delta             0.0000 1.00e-08 Inf -2.00e-08  0.00e+00
# Omicron (BA.1)    0.0000 1.00e-08 Inf -2.00e-08  0.00e+00
# Omicron (BA.2)    0.3653 3.73e-02 Inf  2.92e-01  4.38e-01
# Omicron (BA.2.74) 0.0299 1.32e-02 Inf  4.09e-03  5.58e-02
# Omicron (BA.2.75) 0.1916 3.05e-02 Inf  1.32e-01  2.51e-01
# Omicron (BA.2.76) 0.3473 3.68e-02 Inf  2.75e-01  4.20e-01
# Omicron (BA.3)    0.0000 1.00e-08 Inf -2.00e-08  0.00e+00
# Omicron (BA.4)    0.0120 8.42e-03 Inf -4.52e-03  2.85e-02
# Omicron (BA.5)    0.0539 1.75e-02 Inf  1.96e-02  8.81e-02
# Other             0.0000 1.00e-08 Inf -2.00e-08  0.00e+00
# 
# Confidence level used: 0.95
# NOTE: does not quite match emmeans nnet::multinom output
# above, but maybe because I couldn't run mblogit until convergence

# pairwise contrasts in marginal trends on link scale
# = current differences in growth rate among variants
# confidence intervals are here too wide, probably bug in mblogit or H-D effect
mblogit_emtrends = confint(pairs(emtrends(fit_mblogit, ~ variant,  
                                           var = "collection_date_num",  
                                           mode = "latent",
                                           at = list(collection_date_num = 
                                                       max(data_agbyweek1$collection_date_num))),
                                  reverse=TRUE))
mblogit_emtrends


# equivalent using margineffects package 
mblogit_preds_marginaleffects = predictions(fit_mblogit,
                                             newdata = datagrid(collection_date_num = 
                                                                  max(data_agbyweek1$collection_date_num)),
                                             type="link", # supposedly supported by predict.mblogit
                                             transform_post = insight::link_inverse(fit_nnet))
                                             #transform_post = insight::link_inverse(fit_mblogit)) # function (eta) eta is not correct, inverse link for multinomial should be softmax function
# returns Error in transform_post(draws) : 
# REAL() can only be applied to a 'numeric', not a 'NULL'
# so something going wrong with transform_post
# leaving that out gives predictions on link scale (I think)
# conf intervals too wide though - prob bug in mblogit or H-D effect
mblogit_preds_marginaleffects


# 3. VGAM::vglm multinomial fit ####
# this one supposedly does have a predict.vglm method with type="link",
# https://github.com/cran/VGAM/blob/master/R/predict.vglm.q
library(VGAM)
fit_vglm = vglm(variant ~ ns(collection_date_num, df=2),
                      weight=count+1E-5, # I added 1E-5 to prevent an error
                      family=multinomial,
                      data=dat)
summary(fit_vglm) # reports problems with Hauck-Donner effects
vglm_preds_marginaleffects = predictions(fit_vglm,
                                             newdata = datagrid(collection_date_num = 
                                                                  max(data_agbyweek1$collection_date_num)),
                                             type="response" # here = "probs" 
                                             # type="link", 
                                             # transform_post = insight::link_inverse(fit_vglm))
) 
# results in error Error in `*tmp*`[["coefficients"]] : this S4 class is not subsettable
vglm_preds_marginaleffects 



