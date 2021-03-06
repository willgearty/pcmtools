% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phenotypeThroughTime.R
\name{traitStatsThroughTime}
\alias{traitStatsThroughTime}
\title{Calculate phenotypic statistics through time}
\usage{
traitStatsThroughTime(tree, x, funs = c("mean", "sd"), rnd = 0.1)
}
\arguments{
\item{tree}{an object of class "phylo" or "simmap" containing a phylogenetic tree (optionally with a mapped discrete character)}

\item{x}{a vector of tip values for species; names should correspond to the respective species names in \code{tree$tip.labels}}

\item{rnd}{time bin boundary rounding factor}

\item{fun}{a vector of function names as characters}
}
\value{
A list with two elements:
\item{stats}{A data.frame containing the statistics through time}
\item{xy}{A data.frame containing the edge coordinates for the tree}
}
\description{
This function reconstructs a continuous trait through time and uses these reconstructions to estimate statistics through time.
}
\details{
For each time bin, each of the specified functions in \code{fun} will be carried out across the branches that pass through that time bin. The boundaries of the time bins are defined by the node ages rounded by \code{rnd}. If the tree has a mapped discrete character, the functions will be carried out separately for each state of the character.
}
\examples{
library(phytools)
# simulate tree
tree <- pbtree(n=70,scale=1)

# simulate discrete trait
Q <- matrix(c(-1,1,1,-1),2,2)
rownames(Q) <- colnames(Q)<-c(0,1)
tree <- sim.history(tree,Q)

# simulate continuous trait
x2 <- fastBM(tree,sig2=0.1)

# calculate stats through time
tstt <- traitStatsThroughTime(tree, x2)
coords <- tstt$xy
stats <- tstt$stats

# plot phenogram
ys <- c(coords$y1, coords$y2)
plot(NULL, xlim=c(0,1), ylim=c(min(ys), max(ys)), ylab="trait", xlab="time")
segments(coords$x1, coords$y1, coords$x2, coords$y2)

# plot stats
plot(stats$x, stats$mean)
}
