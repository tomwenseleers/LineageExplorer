---
title: "LineageExplorer: Estimate growth rate advantage of SARS-CoV2 variants of concern based on international genomic surveillance data and multinomial spline fits"
author: "Tom Wenseleers"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This repository uses multinomial spline fits to international SARS-CoV2 genomic surveillance to estimate the growth rate advantage of variants of concern or particular lineages of interest. The fitted model is a single, large `nnet` multinomial model of the form `variant ~ ns(date_num, df=2)+ns(date_num, df=2):region+division`, where `ns()` is a low 2 degree of freedom natural cubic spline and region is typically continent and division is country, or - for larger countries like the USA or India - state or province. This model is similar to some of the multinomial models used in [Davies *et al. Science* (2021)](https://www.science.org/doi/abs/10.1126/science.abg3055) to estimate the transmission advantage of the Alpha variant. In contrast to the model by [Harald Vöhringer and Moritz Gerstung and colleagues](https://www.nature.com/articles/s41586-021-04069-y) (github code [here](https://github.com/gerstung-lab/SARS-CoV-2-International)) or that of [Obermeyer et al. (2022)](https://www.science.org/doi/10.1126/science.abm1208) (github code [here](https://github.com/broadinstitute/pyro-cov)), that use multinomial models that assume the pairwise growth rate differences between variants remain constant through time, this multinomial model allows growth rate advantages to change through time, which for immune escape variants under patterns of waning immunity is a more appropriate assumption. In other aspects, the code of Moritz Gerstung or Fritz Obermeyer is superior - the code being much more professionally written for sure (I haven't had the time yet to clean up my code) and at a technical level, through the use of hierarchical Bayesian methods and the scale of analysis possible via an efficient [NumPyro implementation](https://github.com/broadinstitute/pyro-cov). My model uses a custom-written Rcpp implementation to rapidly calculate the variance-covariance matrix of the multinomial fit, as well as a fork of the `marginaleffects` package with some adaptations to be able to better calculate marginal means and marginal effects (slopes) of multinomial fits.

The current code pulls aggregated data from [cov-spectrum.org](https://github.com/gerstung-lab/SARS-CoV-2-International/blob/main/cov-spectrum.org); the underlying data is proved by [GISAID](https://www.gisaid.org/) (but there is also the option to use open NCBI data). The script uses NextcladePangolin lineage definition, as the regular Pangolin lineages provided by GISAID itself are now lagging too far behind for new variants to be of any use and also misclassify many of the lineages.

## Use of analysis pipeline

The different steps in the analysis are described in the file `global analysis.R`. To use the GISAID data a private CovSPECTRUM key has to be provided in the environment variable `Sys.setenv(REACT_APP_LAPIS_ACCESS_KEY = "XXX")`

Without a key you can still use the open NCBI data by setting `source="NCBI"` instead of using `source="GISAID"` in the beginning of the script.

The output is the estimated current daily growth rate benefit of selected variant lineages relative to a chosen baseline (typically the current dominant lineage) as well as lineage frequencies through time in different countries and/or states.

Here current growth rate advantages of different variants relative to currently predominant type BQ.1\* :

![](plots/GISAID/growth%20rate%20advantage%20VOCs_by%20continent.png)

Here fitted variant frequencies (share of sequenced genomes) through time for some selected countries or states :

![](plots/GISAID/predicted%20lineage%20freqs_New%20York_logit%20scale.png)

![](plots/GISAID/predicted%20lineage%20freqs_Israel_logit%20scale.png)

![](plots/GISAID/predicted%20lineage%20freqs_Denmark_logit%20scale.png)

![](plots/GISAID/predicted%20lineage%20freqs_logit%20scale_selected%20states%20countries.png)

I had some code before also showing new confirmed Covid cases & hospitalisations and estimated infections shown as stacked area charts & plots of Rt values by variant, but this code I still need to update.
