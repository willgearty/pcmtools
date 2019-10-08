[![Build Status](https://travis-ci.com/willgearty/pcmtools.svg?branch=master)](https://travis-ci.com/willgearty/pcmtools)
[![DOI](https://zenodo.org/badge/213449197.svg)](https://zenodo.org/badge/latestdoi/213449197)

# pcmtools
Various tools to help with performing phylogenetic comparative methods and curating/visualizing the results.

## To install
```r
library(devtools)
install_github("willgearty/pcmtools")
```

## To use

### Load package
```r
library(pcmtools)
```

### OUwie Post-Analysis Tools
```r
library(OUwie)
data(tworegime)
ou.results <- list()
ou.results[[1]] <- OUwie(tree,trait,model=c("BM1"))
ou.results[[2]] <- OUwie(tree,trait,model=c("BMS"), root.station = FALSE)
ou.results[[3]] <- OUwie(tree,trait,model=c("OUM"))
ou.results[[4]] <- OUwie(tree,trait,model=c("OUMV"))

#Both regimes have same parameters for BM1 model. Both regimes have different parameters for other models.
regime.mat <- data.frame(BM1 = c(1, 1), BMS = c(1,2), OUM = c(1,2), OUMV = c(1,2), row.names = c(1,2))

#Extract and map parameters to regimes
ou.parameters <- OUwieParSumm(ou.results, regime.mat)

#Summarize model fit
ou.aic <- OUwieAICSumm(ou.parameters)

#Model average parameters
ou.avg <- OUwieModelAvg(ou.parameters)
```

### Coming Soon...
- Plotting phenotypes through time
- Calculating phenotype statistics through time
- Merging and splitting regimes in simmaps
- Plotting OU models
- And more!
