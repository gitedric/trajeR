<p align="center">
  <img src="logo/logotrajeR.png" height="150" />
</p>


**UNDER DEVELOPMENT**

 A R package that fit **regression mixture model - group-based trajectory modeling (GBTM)**. 

Longitudinal studies are often employed on several disciplines like finance, econometrics, psychology, sociology, biology,  ...

The aims of this package is to give tools to work with this situation. On the one hand by divide data in some clusters and on the other hand by fitting longitudinal trajectories in each clusters.

## Package Functionality


TrajeR support several distributions 

- Censored (or regular) Normal distribution ;
- Zero Inflated Poisson distribution ;
- Bernouilli distribution.

The trajectories of each clusters are modelize by poylnomials. We have the possibility to use Non Linear models too.

## Installation

To get the current released version from CRAN:

```{r}
# Comming soon
```

To get the current development version from Github:


```{r}
## install dev version of trajeR from github
devtools::install_github("gitedric/trajeR")

## load trajeR package
library(trajeR)
```

# Usage
The main function of the package is 

```{r}
trajeR(Y, A, X = NULL, TCOV = NULL, ng, degre, degre.nu = 0, 
       Model, Method = "L", 
       ssigma = FALSE, ymax = max(Y) + 1, ymin = min(Y) - 1,
       hessian = TRUE, itermax = 100, paraminit = NULL, 
       EMIRLS = TRUE, refgr = 1,
       fct = NULL, diffct = NULL, nbvar = NULL, nls.lmiter)
```

For more details about the usage of this package you can read the vignette documentation  [here](/vignettes/trajeR_vignette.pdf).

# Contact

Cédric NOEL

cedric.noel@univ-lorraine.fr
