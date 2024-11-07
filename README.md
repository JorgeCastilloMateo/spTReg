spTReg <img src="inst/img/logospTReg.png" width="175px" align="right" />
======================

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spTReg)](https://CRAN.R-project.org/package=spTReg)
[![cran_checks](https://badges.cranchecks.info/worst/spTReg.svg)](https://cran.r-project.org/web/checks/check_results_spTReg.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/spTReg)](https://CRAN.R-project.org/package=spTReg)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/spTReg?color=red)](https://CRAN.R-project.org/package=spTReg)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
  
The R package *spTReg* provides functions to fit, spatially predict, and temporarily forecast space-time dependent data with space-time varying slope coefficients, autoregressive coefficients and scale parameters using: (1) Mean regression using Gaussian processes and Gaussian errors, (2) Quantile regression using Gaussian processes and asymmetric Laplace errors, (3) Probit and logit binary regression using Gaussian processes.

This is a development package. The currently available functionalities are limited to fitting models with iid errors and slope coefficients varying in space-time for points (1) and (2).


## Installation
You can install the **stable** version **(not available, we are working here)** from
[CRAN](https://CRAN.R-project.org/package=spTReg).

```s
install.packages("spTReg")
```

You can install the **development** version from
[GitHub](https://github.com/JorgeCastilloMateo/spTReg)

```s
if (!require("remotes")) install.packages("remotes")
remotes::install_github("JorgeCastilloMateo/spTReg")
```


## References
Castillo-Mateo J, Asín J, Cebrián AC, Gelfand AE, Abaurrea J (2023)
Spatial quantile autoregression for season within year daily maximum temperature data. 
*Annals of Applied Statistics*, **17**(3), 2305--2325.
<doi:10.1214/22-AOAS1719>.

Castillo-Mateo J, Gelfand AE, Gracia-Tabuenca Z, Asín J, Cebrián AC (*in press*).
Spatio-temporal modeling for record-breaking temperature events in Spain. 
*Journal of the American Statistical Association*.

Castillo-Mateo J, Lafuente M, Asín J, Cebrián AC, Gelfand AE, Abaurrea J (2022). 
Spatial modeling of day-within-year temperature time series: an examination of daily maximum temperatures in Aragón, Spain. 
*Journal of Agricultural, Biological and Environmental Statistics*, **27**(3), 487--505. 
<doi:10.1007/s13253-022-00493-3>.
