## LineageExplorer: Estimate growth rate advantage of SARS-CoV2 variants of concern based on international genomic surveillance data and multinomial spline fits

**Tom Wenseleers**

Last update: 16 January 2023

This repository uses multinomial spline fits to international SARS-CoV2 genomic surveillance data to estimate the growth rate advantage of variants of concern or particular lineages of interest. The fitted model is a single, large `nnet` multinomial model of the form `variant ~ ns(date_num, df=2)+ns(date_num, df=2):region+division`, where `ns()` is a low 2 degree of freedom natural cubic spline and region is typically continent and division is country, or - for a large country with decent genomic surveillance like the USA - state. This model is similar to some of the multinomial models that were used in [Davies *et al. Science* (2021)](https://www.science.org/doi/abs/10.1126/science.abg3055) to estimate the transmission advantage of the Alpha variant. The advantage of multinomial spline models is that they allow growth rate advantages to change through time, which for immune escape variants under patterns of waning immunity is the most appropriate assumption. This is the key difference compared to the plain multinomial models used by [Harald Vöhringer and Moritz Gerstung and colleagues](https://www.nature.com/articles/s41586-021-04069-y) (github code [here](https://github.com/gerstung-lab/SARS-CoV-2-International)) and [Obermeyer et al. (2022)](https://www.science.org/doi/10.1126/science.abm1208) (github code [here](https://github.com/broadinstitute/pyro-cov)), which make the implicit assumptions that pairwise growth rate differences between variants remain constant through time, which in practice is not the case. In other aspects, the code of Moritz Gerstung or Fritz Obermeyer is superior - the code being much more professionally written for sure (I haven't had the time yet to clean up my code) and at a technical level, through the use of hierarchical Bayesian methods and the scale of analysis possible via an efficient [NumPyro implementation](https://github.com/broadinstitute/pyro-cov). My model uses a custom-written Rcpp implementation to rapidly calculate the variance-covariance matrix of the multinomial fit, as well as a fork of the `marginaleffects` package with some adaptations to be able to better calculate marginal means and marginal effects (slopes) of multinomial fits.

The current code pulls aggregated data from [cov-spectrum.org](https://github.com/gerstung-lab/SARS-CoV-2-International/blob/main/cov-spectrum.org); the underlying data is provided by [GISAID](https://www.gisaid.org/) (but there is also the option to use open [NCBI](https://www.ncbi.nlm.nih.gov/) data). The script uses NextcladePangolin lineage definitions, as the regular Pangolin lineages provided by GISAID itself are now lagging too far behind for new variants to be of any use and also misclassify many of the lineages.

## Use of analysis pipeline

The different steps in the analysis are described in the file `global analysis.R`. To use the GISAID data a private CovSPECTRUM key has to be provided in the environment variable `Sys.setenv(REACT_APP_LAPIS_ACCESS_KEY = "XXX")`

Without a key you can still use the open NCBI data by setting `source="NCBI"` instead of using `source="GISAID"` in the beginning of the script.

The output is the estimated current daily growth rate benefit of selected variant lineages relative to a chosen baseline (typically the current dominant lineage) as well as lineage frequencies through time in different countries and/or states.

## Estimated current growth rate advantages

Here current growth rate advantages of different variants relative to currently predominant type BQ.1\* (continuous growth rate advantage per day in %) :

![](plots/GISAID/growth%20rate%20advantage%20VOCs_by%20continent.png)

## Estimated lineage frequencies

Here fitted variant frequencies (share of sequenced genomes) through time for some selected countries or states :

![](plots/GISAID/predicted%20lineage%20freqs_logit%20scale.png)

## Limitations

-   Code is not yet completely finished and not cleaned up. I had some code before also showing new confirmed Covid cases & hospitalisations and estimated infections by variant shown as stacked area charts & plots of Rt values by variant, but this code I still need to update.

-   Planned improvements (if I find the time): clean up code & wrap everything in proper functions, define lineages in a more generic way like [here](https://nbviewer.org/github/gerstung-lab/SARS-CoV-2-International/blob/main/genomicsurveillance-int.ipynb#Check-some-fast-growing-lineages), use an `mclogit::mblogit` multinomial fit instead of a `nnet::multinom` fit to take into account overdispersion, add plots of confirmed Covid cases, hospitalisations & deaths (working from estimated excess mortality, splitting contributions of different Covid & influenza variants [from FluNet] & temperature anomalies) and estimated infections (from IHME estimates) by variant & plots of Rt values by variant, add plots of growth rate advantage of each variant at invasion (invasion fitness), add an explicitly spatial multinomial tensor spline fit in function of latitude, longitude & time (requires geocoding all the GISAID data and correcting the many mistakes - has the advantage that it would also return interpolated lineage frequencies for countries with little or no genomic surveillance), eventually package all the code in proper functions and a proper R package.

## Acknowledgements

This repository uses the [CoV-Spectrum](https://cov-spectrum.org/) [LAPIS API](https://github.com/cevo-public/LAPIS) to import aggregated NextCladePangolin lineage counts from [GISAID](https://gisaid.org/) (if you have a private key and a GISAID login) or [NCBI](https://www.ncbi.nlm.nih.gov/) (open, see documention [here](https://lapis-docs.readthedocs.io/)). I would like to acknowledge [Cov-Spectrum](https://cov-spectrum.org/) for making this resource available, [GISAID](https://gisaid.org/) and [NCBI](https://www.ncbi.nlm.nih.gov/) for providing repositories and all the genomic surveillance data contributors and variant trackers for their work on SARS-Cov2 genomic surveillance. Thanks to Moritz Gerstung and Cornelius Römer for many useful discussions and coding advice.
