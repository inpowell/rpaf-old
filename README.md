
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
smoke_data <- gen_data(minifhs, "ID", c(0,5,10,15,17), "DEATH", "DEATH_FT", 
                       variables = c("B_COHORT", "SEX", "SMOKE"))
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

References
----------

-   Laaksonen, M.A., Virtala, E., Knekt, P., Oja, H. and Härkänen, T., 2011. SAS macros for calculation of population attributable fraction in a cohort study design. *J Stat Softw*, 43(7), pp.1-25.
