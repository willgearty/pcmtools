% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OUwieAvg.R
\name{OUwieParSumm}
\alias{OUwieParSumm}
\title{Extract parameters from OUwie results}
\usage{
OUwieParSumm(ou.results, regime.mat, params = c("Alpha", "Sigma.sq",
  "Theta", "Theta.se", "AICc"))
}
\arguments{
\item{ou.results}{A list of lists (or just a list) of unmodified results from an OUwie analysis}

\item{regime.mat}{A data frame mapping regimes to total regime options for each model (see details)}

\item{params}{A vector specifying which parameter should be calculated/returned (see details)}
}
\value{
An \code{ouwiepars} object. Basically a 4-D array that can be passed to other functions for further analysis/visualization.
}
\description{
Extract various parameter values from the results of (many) OUwie analyses and maps the parameter values to different regimes based on the inputted regime map.
Returns the parameters in a 4-D array, where the first dimension is the models, the second dimension is the parameters (including AICc), the third dimension is the regimes, and the fourth dimension is the replicates.
}
\details{
The \code{regime.mat} is the most important component, as it indicates which parameters should be mapped to which regime for each model.
For example, in an OU1 model, the user would likely want the parameters mapped to all regimes, whereas in an OUM model, the user would likely want the parameters for each regime mapped exclusively to that regime.
In more complex scenarios, the user may have multiple OUM models in which regimes are split or combined in different manners, such that the parameters for one regime in one OUM model may map to multiple regimes in the overall dataset.
The \code{rownames} of this matrix should identify names for the regimes and the \code{colnames} should identify the models.
It is assumed that the order of the models/\code{colnames} in \code{regime.mat} matches the order of the models in \code{ou.results}.

Valid options for \code{params} are "Alpha", "Sigma.sq", "Theta", "Theta.se", "Halflife" (phylogenetic half-life), "Stat.var" (stationary variance), "AIC", "AICc", and "BIC".
}
\examples{
\dontrun{
library(OUwie)
data(tworegime)
ou.results <- list()
ou.results[[1]] <- OUwie(tree,trait,model=c("BM1"))
ou.results[[2]] <- OUwie(tree,trait,model=c("BMS"), root.station = FALSE)
ou.results[[3]] <- OUwie(tree,trait,model=c("OUM"))
ou.results[[4]] <- OUwie(tree,trait,model=c("OUMV"))

#Both regimes have same parameters for BM1 model. Both regimes have different parameters for other models.
regime.mat <- data.frame(BM1 = c(1, 1), BMS = c(1,2), OUM = c(1,2), OUMV = c(1,2), row.names = c(1,2))}

OUwieParSumm(ou.results, regime.mat)
}
