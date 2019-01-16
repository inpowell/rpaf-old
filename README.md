
<!-- README.md is generated from README.Rmd. Please edit that file -->
rpaf
====

The goal of rpaf is to allow for simple calculation of population attributable fractions (PAFs) for either disease or mortality using data from cohort studies.

It is intended as an R equivalent to the SAS macros developed by [Laaksonen et al. (2011)](#references).

Installation
------------

You **cannot** currently install the released version of rpaf from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rpaf")
```

However, you can install the development version of rpaf from [GitHub](https://github.com/inpowell/rpaf) using:

``` r
devtools::install_github("https://github.com/inpowell/rpaf")
```

Example
-------

Basic calculation of PAFs can be achieved using the `mpaf_summary` and `dpaf_summary` functions for mortality and disease PAFs respectively. For instance, in the case of the former, we may wish to determine the PAFs for mortality if everone who currently smokes in the Mini-Finland Health Survey is modified to have never smoked:

``` r
library(rpaf)
smoke_data <- gen_data(
  indata = minifhs, 
  id_var = "ID",
  ft_breaks = c(0,5,10,15,17),
  death_ind = "DEATH", death_time = "DEATH_FT",
  variables = c("B_COHORT", "SEX", "SMOKE")
)
smoke_summary <- mpaf_summary(
  survival::Surv(f_end, DEATH) ~ B_COHORT * f_period + SEX + SMOKE,
  mpaf_data = smoke_data, 
  modifications = list(SMOKE = c("Never", "<30/day", ">=30/day")),
  covar_model = c("SEX", "SMOKE")
)

knitr::kable(rbind(smoke_summary$paf, smoke_summary$paf0), digits = 4)
```

|          |     PAF|   2.5 %|  97.5 %|
|----------|-------:|-------:|-------:|
| (0,5\]   |  0.1486|  0.1134|  0.1825|
| (5,10\]  |  0.1169|  0.0895|  0.1434|
| (10,15\] |  0.0871|  0.0665|  0.1073|
| (15,17\] |  0.0787|  0.0537|  0.1030|
| (0,5\]   |  0.1486|  0.1134|  0.1825|
| (0,10\]  |  0.1302|  0.1000|  0.1593|
| (0,15\]  |  0.1117|  0.0862|  0.1365|
| (0,17\]  |  0.1066|  0.0823|  0.1304|

The process to calculate a disease PAF is similar, though the function calls vary slightly. Here we want to see the PAF for diabetes incidence if everyone in the study who has a high BMI (that is, above 25.0 kg/m^2) is instead modified to have a normal or low BMI:

``` r
bmi_data <- gen_data(
  indata = minifhs,
  id_var = "ID", ft_breaks = c(0,5,10,15,17), 
  death_ind = "DEATH", death_time = "DEATH_FT",
  disease_ind = "DIAB", disease_time = "DIAB_FT", 
  variables = c("B_COHORT", "SEX", "BMI_2")
)

bmi_summary <- dpaf_summary(
  disease_resp = survival::Surv(f_end, DIAB) ~ .,
  death_resp = survival::Surv(f_end, DEATH) ~ .,
  predictors = ~ B_COHORT * f_period + SEX + BMI_2,
  dpaf_data = bmi_data,
  modifications = list(BMI_2 = "<25.0"), 
  covar_model = c("SEX", "BMI_2")
)

knitr::kable(rbind(bmi_summary$paf, bmi_summary$paf0), digits = 4)
```

|          |     PAF|   2.5 %|  97.5 %|
|----------|-------:|-------:|-------:|
| (0,5\]   |  0.6811|  0.5490|  0.7745|
| (5,10\]  |  0.6801|  0.5491|  0.7731|
| (10,15\] |  0.6802|  0.5497|  0.7728|
| (15,17\] |  0.6795|  0.5476|  0.7729|
| (0,5\]   |  0.6811|  0.5490|  0.7745|
| (0,10\]  |  0.6804|  0.5492|  0.7735|
| (0,15\]  |  0.6803|  0.5495|  0.7732|
| (0,17\]  |  0.6803|  0.5496|  0.7730|

Note that in both these cases, the period and follow-up time variables have been automatically named as `f_period` and `f_end` respectively. These names can be chosen in the call to `gen_data`.

References
----------

-   Laaksonen, M.A., Virtala, E., Knekt, P., Oja, H. and Härkänen, T., 2011. SAS macros for calculation of population attributable fraction in a cohort study design. *J Stat Softw*, 43(7), pp.1-25.
