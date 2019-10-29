#' Extract posterior values from densityMap and/or contMap objects
#'
#' This function takes a list of densityMap and/or contMap objects and extracts and aggregates the posterior values across them. The age and edge index of each posterior value are also included.
#'
#' Assumes that the gradient resolution is identical for all objects (specified with \code{res} in \code{densityMap} and \code{contMap}). If the tree has a \code{root.time} element, this will be used as the max age, otherwise \code{max(nodeHeights(tree))} will be used.
#'
#' @param ... A vector of values.
#' @return A data.frame where rows represent the mapped traits through time and columns represent each of the supplied objects in \code{...}. The column names are set to the names of objects. Also includes a column for the edge index and a column for the age (in the units of the tree).
#' @importFrom zoo rollmean
#' @importFrom phytools nodeHeights
#' @export
#' @examples
#' library(phytools)
#' # simulate tree
#' tree <- pbtree(n=70,scale=1)
#'
#' # simulate discrete trait
#' Q <- matrix(c(-1,1,1,-1),2,2)
#' rownames(Q) <- colnames(Q)<-c(0,1)
#' x1 <- sim.history(tree,Q)$states
#'
#' # generate stochastic maps and density map
#' mtrees <- make.simmap(tree,x1,nsim=100)
#' map1 <- densityMap(mtrees)
#'
#' # simulate continuous trait
#' x2 <- fastBM(tree,sig2=0.1)
#'
#' # generate cont map
#' map2 <- contMap(tree, x2)
#'
#' # extract posterior densities
#' pp_data <- extractSimmapDensity(map1, map2)
#'
#' # plot to see how traits have evolved with respect to one another
#' plot(pp_data$map1, pp_data$map2)
extractSimmapDensity <- function(...) {
  objs <- list(...)
  obj_names <- sapply(substitute({...})[-1], deparse)
  pp_data <- data.frame(matrix(NA, ncol = length(objs) + 2, nrow = sum(sapply(objs[[1]]$tree$maps, length)), dimnames = list(NULL,c(obj_names, "edge", "age"))))
  breaks <- length(objs[[1]]$cols) - 1
  heights <- nodeHeights(objs[[1]]$tree)
  maxheight <- objs[[1]]$tree$root.time
  if(is.null(maxheight)) maxheight <- max(heights)
  n <- 1
  for(i in 1:length(objs[[1]]$tree$maps)) {
    sgmnts <- length(objs[[1]]$tree$maps[[i]])
    for(j in 1:length(objs)) {
      vals <- as.numeric(names(objs[[j]]$tree$maps[[i]]))/breaks
      if(is(objs[[j]], "contMap")) {
        vals <- vals * diff(objs[[j]]$lims) + objs[[j]]$lims[1]
      }
      pp_data[n:(n + sgmnts - 1), obj_names[j]] <- vals
    }
    pp_data$edge[n:(n + sgmnts - 1)] <- i
    pp_data$age[n:(n + sgmnts - 1)] <- maxheight - rollmean(c(0, cumsum(objs[[j]]$tree$maps[[i]])) + heights[i,1], 2)
    n <- n + sgmnts
  }
  return(pp_data)
}


