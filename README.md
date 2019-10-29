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

###Rename, merge, or split mapped states
```r
library(phytools)
#Simulate a mapped tree
set.seed(4)
Q <- matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
rownames(Q) <- colnames(Q) <- letters[1:3]
tree <- sim.history(pbtree(n=100,scale=1),Q)
cols <- setNames(c("blue","red","green","orange"),letters[1:4])

#Plot the mapping
plot(tree, cols, ftype="i", fsize=0.7)

#Split state c to state d within subclade
tree2 <- mergeMappedStates2(tree, "c", "d", 173)

#Plot new mapping
plot(tree2, cols, ftype="i", fsize=0.7)
```

### Extract posterior values from densityMaps and contMaps
```r
library(phytools)
#Simulate tree
tree <- pbtree(n=70,scale=1)

#Simulate discrete trait
Q <- matrix(c(-1,1,1,-1),2,2)
rownames(Q) <- colnames(Q)<-c(0,1)
x1 <- sim.history(tree,Q)$states

#Generate stochastic maps and density map
mtrees <- make.simmap(tree,x1,nsim=100)
map1 <- densityMap(mtrees)

#Simulate continuous trait
x2 <- fastBM(tree,sig2=0.1)

#Generate cont map
map2 <- contMap(tree, x2)

#Extract posterior densities
pp_data <- extractSimmapDensity(map1, map2)

#Plot to see how traits have evolved with respect to one another
plot(pp_data$map1, pp_data$map2)
```

### Coming Soon...
- Plotting phenotypes through time
- Calculating phenotype statistics through time
- Merging and splitting regimes in simmaps
- Plotting OU models
- And more!
