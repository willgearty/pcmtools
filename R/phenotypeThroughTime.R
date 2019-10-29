#calculates the x, y coordinates (where x is time and y is a quantitative character) of each edge
#tree should be an object of class phylo,
#x should be a vector of tip values for species; names(x) should be the species names
#' @importFrom phytools fastAnc
#' @importFrom ape Ntip
edgeCoord <- function(tree, x){
  try(if(!all(names(x) %in% tree$tip.label)) stop("vector names should match tip labels"))
  try(if(!all(tree$tip.label %in% names(x))) stop("vector should have measurement for every tip"))
  x <- x[tree$tip.label]
  edge <- tree$edge
  states <- setNames(c(x, fastAnc(tree, x)), seq(1, Ntip(tree) + Nnode(tree)))
  heights <- nodeHeights(tree)
  xy <- as.data.frame(cbind(heights[,1], states[edge[,1]], heights[,2], states[edge[,2]]))
  names(xy) <- c("x1", "y1", "x2", "y2")
  xy$slope <- (xy$y2 - xy$y1)/(xy$x2 - xy$x1)
  row.names(xy) <- seq(1,nrow(xy))
  return(xy)
}

calcStats <- function(xy, funs = c("mean", "sd"), rnd = .1) {
  xs <- sort(unique(round(c(xy$x1,xy$x2)/rnd)*rnd))
  states <- unique(xy$regime)
  n <- length(states)
  stats <- data.frame(x = sort(rep(xs, n)), regime = rep(NA, length(xs)*n), n = rep(NA, length(xs)*n))
  stats[funs] <- NA
  for(i in 1:length(xs)){
    for(j in 1:n){
      xy_sub <- subset(xy, x2 > xs[i] & x1 <= xs[i] & regime == states[j])
      xy_sub$y <- xy_sub$y1 + xy_sub$slope*(xs[i]-xy_sub$x1)
      stats$regime[(i-1)*n + j] <- states[j]
      stats$n[(i-1)*n + j] <- num <- nrow(xy_sub)
      for(fun in funs){
        stats[(i-1)*n + j, fun] <- get(fun)(xy_sub$y)
      }
      # stats$mean[(i-1)*n + j] <- mean(xy_sub$y)
      # stats$sd[(i-1)*n + j] <- sd(xy_sub$y)
      # stats[(i-1)*n + j, c("fifth","ninety_fifth")] <- quantile(xy_sub$y,  probs = c(5, 95)/100)
      # stats$skewness[(i-1)*n + j] <- skewness(xy_sub$y)
      # stats$ses[(i-1)*n + j] <- sqrt(6*num*(num-1)/((num-2)*(num+1)*(num+3)))
    }
  }
  return(stats)
}


#' Calculate phenotypic statistics through time
#'
#' This function reconstructs a continuous trait through time and uses these reconstructions to estimate statistics through time.
#'
#' For each time bin, each of the specified functions in \code{fun} will be carried out across the branches that pass through that time bin. The boundaries of the time bins are defined by the node ages rounded by \code{rnd}. If the tree has a mapped discrete character, the functions will be carried out separately for each state of the character.
#'
#' @param tree an object of class "phylo" or "simmap" containing a phylogenetic tree (optionally with a mapped discrete character)
#' @param x a vector of tip values for species; names should correspond to the respective species names in \code{tree$tip.labels}
#' @param fun a vector of function names as characters
#' @param rnd time bin boundary rounding factor
#' @return A list with two elements:
#' \item{stats}{A data.frame containing the statistics through time}
#' \item{xy}{A data.frame containing the edge coordinates for the tree}
#' @export
#' @examples
#' library(phytools)
#' # simulate tree
#' tree <- pbtree(n=70,scale=1)
#'
#' # simulate discrete trait
#' Q <- matrix(c(-1,1,1,-1),2,2)
#' rownames(Q) <- colnames(Q)<-c(0,1)
#' tree <- sim.history(tree,Q)
#'
#' # simulate continuous trait
#' x2 <- fastBM(tree,sig2=0.1)
#'
#' # calculate stats through time
#' tstt <- traitStatsThroughTime(tree, x2)
#' coords <- tstt$xy
#' stats <- tstt$stats
#'
#' # plot phenogram
#' ys <- c(coords$y1, coords$y2)
#' plot(NULL, xlim=c(0,1), ylim=c(min(ys), max(ys)), ylab="trait", xlab="time")
#' segments(coords$x1, coords$y1, coords$x2, coords$y2)
#'
#' # plot stats
#' plot(stats$x, stats$mean)
traitStatsThroughTime <- function(tree, x, funs = c("mean", "sd"), rnd = .1){
  xy <- edgeCoord(tree, x)
  if(exists("maps", where = tree)){
    n_states <- vector("numeric", length(tree$maps))
    for(i in 1:length(tree$maps)){
      n_states[i] <- length(tree$maps[[i]])
    }
    tot_states <- sum(n_states)
    temp <- data.frame(x1 = rep(NA, length(tot_states)), y1 = rep(NA, length(tot_states)), x2 = rep(NA, length(tot_states)), y2 = rep(NA, length(tot_states)), slope = rep(NA, length(tot_states)), regime = rep(NA, length(tot_states)))
    k <- 1
    for(i in 1:length(tree$maps)){
      x1 <- xy$x1[i]
      y1 <- xy$y1[i]
      slope <- xy$slope[i]
      for(j in 1:n_states[i]){
        dx <- tree$maps[[i]][j]
        x2 <- x1 + dx
        dy <- slope*dx
        y2 <- y1 + dy
        temp[k,1:5] <- c(x1, y1, x2, y2, slope)
        temp$regime[k] <- names(tree$maps[[i]])[j]
        x1 <- x2
        y1 <- y2
        k <- k + 1
      }
    }
    xy <- temp
  }else{
    xy$regime <- 1
  }
  stats <- calcStats(xy, funs, rnd)

  return(list(stats = stats, xy = xy))
}

