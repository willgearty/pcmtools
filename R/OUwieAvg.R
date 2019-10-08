#' Calculate AIC weights
#'
#' This function takes a vector of AIC (Akaike Information Criterion) values and returns a vector of AIC weights using the formula from Burnham and Anderson (2002).
#'
#' If \code{na.rm = FALSE} and any values in \code{AIC} are \code{NA}, all returned values will be \code{NA}.
#'
#' @param AIC A vector of values.
#' @param na.rm Whether to remove NA values.
#' @return A named vector of weights with names inherited from \code{AIC}.
#' @export
#' @examples
#' AIC <- c(NA, 5, 10, 20, 25)
#' #ignore NAs
#' AICweights(AIC)
#'
#' #should return all NAs
#' AICweights(AIC, na.rm = FALSE)
AICweights <- function(AIC, na.rm = TRUE){
  deltAIC <- deltaAIC(AIC, na.rm = na.rm)
  expAIC <- exp(-.5 * deltAIC)
  weights <- expAIC / sum(expAIC, na.rm = na.rm)
  names(weights) <- names(AIC)
  return(weights)
}

#' Calculate deltaAIC
#'
#' Calculate deltaAIC (Akaike Information Criterion), the absolute difference between the lowest AIC value and the other AIC values.
#'
#' If \code{na.rm = FALSE} and any values in \code{AIC} are \code{NA}, all returned values will be \code{NA}.
#'
#' @param AIC A vector of values.
#' @param na.rm Whether to remove NA values.
#' @return A vector of weights.
#' @export
#' @examples
#' AIC <- c(NA, 5, 10, 20, 25)
#' deltaAIC(AIC)
#'
#' #should return all NAs
#' deltaAIC(AIC, na.rm = FALSE)
deltaAIC <- function(AIC, na.rm = TRUE){
  minAIC <- min(AIC, na.rm = na.rm)
  deltAIC <- AIC - minAIC
  return(deltAIC)
}

#' Extract parameters from OUwie results
#'
#' Extract various parameter values from the results of (many) OUwie analyses and maps the parameter values to different regimes based on the inputted regime map.
#' Returns the parameters in a 4-D array, where the first dimension is the models, the second dimension is the parameters (including AICc), the third dimension is the regimes, and the fourth dimension is the replicates.
#'
#' The \code{regime.mat} is the most important component, as it indicates which parameters should be mapped to which regime for each model.
#' For example, in an OU1 model, the user would likely want the parameters mapped to all regimes, whereas in an OUM model, the user would likely want the parameters for each regime mapped exclusively to that regime.
#' In more complex scenarios, the user may have multiple OUM models in which regimes are split or combined in different manners, such that the parameters for one regime in one OUM model may map to multiple regimes in the overall dataset.
#' The \code{rownames} of this matrix should identify names for the regimes and the \code{colnames} should identify the models.
#' It is assumed that the order of the models/\code{colnames} in \code{regime.mat} matches the order of the models in \code{ou.results}.
#'
#' Valid options for \code{params} are "Alpha", "Sigma.sq", "Theta", "Theta.se", "Halflife" (phylogenetic half-life), "Stat.var" (stationary variance), "AIC", "AICc", and "BIC".
#'
#' @param ou.results A list of lists (or just a list) of unmodified results from an OUwie analysis
#' @param regime.mat A data frame mapping regimes to total regime options for each model (see details)
#' @param params A vector specifying which parameter should be calculated/returned (see details)
#' @return An \code{ouwiepars} object. Basically a 4-D array that can be passed to other functions for further analysis/visualization.
#' @export
#' @examples
#' \dontrun{
#' library(OUwie)
#' data(tworegime)
#' ou.results <- list()
#' ou.results[[1]] <- OUwie(tree,trait,model=c("BM1"))
#' ou.results[[2]] <- OUwie(tree,trait,model=c("BMS"), root.station = FALSE)
#' ou.results[[3]] <- OUwie(tree,trait,model=c("OUM"))
#' ou.results[[4]] <- OUwie(tree,trait,model=c("OUMV"))
#'
#' #Both regimes have same parameters for BM1 model. Both regimes have different parameters for other models.
#' regime.mat <- data.frame(BM1 = c(1, 1), BMS = c(1,2), OUM = c(1,2), OUMV = c(1,2), row.names = c(1,2))}
#'
#' OUwieParSumm(ou.results, regime.mat)
OUwieParSumm <- function(ou.results, regime.mat, params = c("Alpha","Sigma.sq","Theta","Theta.se", "AICc")){
  if(is(ou.results[[1]], "OUwie")) ou.results <- list(ou.results)
  nruns <- length(ou.results)
  regimes <- rownames(regime.mat)
  nregs <- length(regimes)
  mods <- colnames(regime.mat)
  nmods <- length(mods)
  nparams <- length(params)
  ou.parameters <- array(NA, dim=c(nmods, nparams, nregs, nruns), dimnames=list(mods, params, regimes))

  for(j in 1:nruns){
    for(i in 1:nmods){
      reg.temp <- colnames(ou.results[[j]][[i]]$solution)
      nreg.temp <- length(reg.temp)
      #loop through parameters
      for(param in params){
        #Record AICc values
        if(param == "AIC"){
          ou.parameters[i, "AIC", , j] <- ou.results[[j]][[i]]$AIC
        } else if(param == "AICc"){
          ou.parameters[i, "AICc", , j] <- ou.results[[j]][[i]]$AICc
        } else if(param == "BIC"){
          ou.parameters[i, "BIC", , j] <- ou.results[[j]][[i]]$BIC
        } else if(param == "Alpha"){
          if(ou.results[[j]][[i]]$model=="BMS" | ou.results[[j]][[i]]$model=="BM1"){
            ou.parameters[i, param, , j] <- 0
          } else {
            for(k in 1:nreg.temp){
              ou.parameters[i, param, which(regime.mat[,i] == reg.temp[k], useNames = FALSE), j] <- ou.results[[j]][[i]]$solution[1,k]
            }
          }
        } else if(param == "Sigma.sq"){
          for(k in 1:nreg.temp){
            ou.parameters[i, param, which(regime.mat[,i] == reg.temp[k], useNames = FALSE), j] <- ou.results[[j]][[i]]$solution[2,k]
          }
        } else if(param == "Theta"){
          if(ou.results[[j]][[i]]$model=="BM1" | ou.results[[j]][[i]]$model=="OU1" | ou.results[[j]][[i]]$model=="BMS"){
            ou.parameters[i, param, , j] <- ou.results[[j]][[i]]$theta[1,1]
          } else {
            for(k in 1:nreg.temp){
              ou.parameters[i, param, which(regime.mat[,i] == reg.temp[k], useNames = FALSE), j] <- ou.results[[j]][[i]]$theta[k,1]
            }
          }
        } else if(param == "Theta.se"){
          if(ou.results[[j]][[i]]$model=="BM1" | ou.results[[j]][[i]]$model=="OU1" | ou.results[[j]][[i]]$model=="BMS"){
            ou.parameters[i, param, , j] <- ou.results[[j]][[i]]$theta[1,2]
          } else {
            for(k in 1:nreg.temp){
              ou.parameters[i, param, which(regime.mat[,i] == reg.temp[k], useNames = FALSE), j] <- ou.results[[j]][[i]]$theta[k,2]
            }
          }
        } else if(param == "Halflife"){
          if(ou.results[[j]][[i]]$model=="BMS" | ou.results[[j]][[i]]$model=="BM1"){
            ou.parameters[i, param, , j] <- NA
          } else {
            for(k in 1:nreg.temp){
              ou.parameters[i, param, which(regime.mat[,i] == reg.temp[k], useNames = FALSE), j] <- log(2)/ou.results[[j]][[i]]$solution[1,k]
            }
          }
        } else if(param == "Stat.var"){
          if(ou.results[[j]][[i]]$model=="BMS" | ou.results[[j]][[i]]$model=="BM1"){
            ou.parameters[i, param, , j] <- NA
          } else {
            for(k in 1:nreg.temp){
              ou.parameters[i, param, which(regime.mat[,i] == reg.temp[k], useNames = FALSE), j] <- ou.results[[j]][[i]]$solution[2,k]/(2 * ou.results[[j]][[i]]$solution[1,k])
            }
          }
        }
      }
    }
  }
  class(ou.parameters) <- "ouwiepars"
  return(ou.parameters)
}

#' Model average the parameters across (many) OUwie results using AICc
#'
#' Internally calculates AICc weights and uses them to model average the parameter values output from \code{OUwieParSumm}.
#'
#' \code{na.rm = TRUE} will remove any models where \code{AICc == NA}, but will use other models for model averaging;
#' \code{na.rm = FALSE} will remove any replicates (rows of the output) with any models where \code{AICc == NA}.
#'
#' @param ou.parameters An object of class \code{ouwiepars} as output from \code{OUwieParSumm} or \code{CleanOUwieParameters}
#' @param OU.only Whether any Brownian motion models should be dropped before model averaging
#' @param na.rm Whether to ignore model results with \code{NA} AICc values
#' @return A list with two elements:
#' \item{Weights}{A data.frame of the AICc weights for each model across the replicates}
#' \item{Counts}{A named vector giving the number of replicates in which each model has the highest AICc weight}
#' @export
#' @examples
#' ou.parameters <- OUwieParSumm(ou.results, regime.mat)
#' OUwieModelAvg(ou.parameters)
OUwieModelAvg  <- function(ou.parameters, OU.only = FALSE, na.rm = TRUE){
  nmods <- dim(ou.parameters)[1]
  mods <- dimnames(ou.parameters)[[1]]
  params <- dimnames(ou.parameters)[[2]][!(dimnames(ou.parameters)[[2]] == "AICc")]
  nparams <- length(params)
  nregs <- dim(ou.parameters)[3]
  regs <- dimnames(ou.parameters)[[3]]
  nruns <- dim(ou.parameters)[4]
  ou.avg <- array(NA, dim=c(nparams, nregs, nruns), dimnames=list(params, regs))
  use <- vector(mode = "logical", length = nmods)
  if(OU.only) use <- !(grepl("BM", mods)) else use <- rep(TRUE, nmods)
  for(i in 1:nruns){
    aicc.weights <- AICweights(ou.parameters[use, "AICc", 1, i], na.rm = na.rm)
    #loop through parameters
    for (param in params){
      ou.avg[param, , i] <- colSums(aicc.weights * ou.parameters[use, param, , i], na.rm = na.rm)
    }

    #need to calculate model-averaged theta.se values differently
    #from page 162 of Burnham and Anderson 2002
    if("Theta.se" %in% params & "Theta" %in% params){
      for (j in 1:nregs){
        variance <- (ou.parameters[use, "Theta.se", j, i])^2
        theta.avg <- ou.avg["Theta", j, i]
        ou.avg["Theta.se", j, i] <- sum(aicc.weights * sqrt(variance + (ou.parameters[use, "Theta", j, i] - theta.avg)^2), na.rm = na.rm)
      }
    }
  }
  return(ou.avg)
}

#' Clean extracted parameters from OUwie results
#'
#' Cleans the parameters that are extracted by \code{OUwieParamSum} based on user-specified upper and lower bounds.
#'
#' Parameter estimates outside of these bounds will result in the AICc being changed to NA,
#' which will affect downstream model averaging (see \code{OUwieModelAvg}).
#'
#' @param ou.parameters An object of class \code{ouwiepars} as output from \code{OUwieParSumm}
#' @param lower A list of lower bounds for model parameters
#' @param upper A list of upper bounds for model parameters
#' @return An \code{ouwiepars} object
#' @export
#' @examples
#' ou.parameters <- OUwieParSumm(ou.results, regime.mat)
#'
#' #Sets the AICc for the BM1 model to NA, so it wouldn't be included in downstream model averaging
#' OUwieCleanPar(ou.parameters, upper = list("AICc" = 45))
OUwieCleanPar <- function(ou.parameters, lower = list(), upper = list()){
  nruns <- dim(ou.parameters)[4]
  params <- dimnames(ou.parameters)[[2]]
  try(if(!all(names(lower) %in% params))
    stop(paste("Incorrect lower bound parameter(s) specified: ",paste0(names(lower)[!(names(lower) %in% params)],collapse=", "))))
  try(if(!all(names(upper) %in% params))
    stop(paste("Incorrect upper bound parameter(s) specified: ",paste0(names(upper)[!(names(upper) %in% params)],collapse=", "))))

  for(i in 1:nruns){
    for(param in names(lower)){
      ou.parameters[which(rowSums(ou.parameters[ , param, , i] < lower[[param]]) > 0), "AICc", , i] <- NA
    }
    for(param in names(upper)){
      ou.parameters[which(rowSums(ou.parameters[ , param, , i] > upper[[param]]) > 0), "AICc", , i] <- NA
    }
  }
  return(ou.parameters)
}

#' Summarize the model fit across (many) OUwie results using AICc
#'
#' Returns AICc weights and counts for OUwie results that have been processed using \code{OUwieParSumm}.
#'
#' \code{na.rm = TRUE} will remove any models where \code{AICc == NA}, but will use other models for AIC weight calculation;
#' \code{na.rm = FALSE} will remove any replicates (rows of the output) with any models where \code{AICc == NA}.
#'
#' @param ou.parameters An object of class \code{ouwiepars} as output from \code{OUwieParSumm} or \code{CleanOUwieParameters}
#' @param na.rm Whether to ignore model results with \code{NA} AICc values
#' @return A list with two elements:
#' \item{Weights}{A data.frame of the AICc weights for each model across the replicates}
#' \item{Counts}{A named vector giving the number of replicates in which each model has the highest AICc weight}
#' @export
#' @examples
#' ou.parameters <- OUwieParSumm(ou.results, regime.mat)
#' OUwieAICSumm(ou.parameters)
OUwieAICSumm <- function(ou.parameters, na.rm = TRUE){
  nruns <- dim(ou.parameters)[4]
  mods <- dimnames(ou.parameters)[[1]]
  weights <- as.data.frame(array(NA, dim=c(nruns,length(mods)), dimnames=list(seq(1,nruns),mods)))
  for(i in 1:nruns){
    weights[i,] <- AICweights(ou.parameters[ , "AICc", 1, i], na.rm = na.rm)
  }
  counts <- table(factor(unlist(apply(weights,1,which.max)),levels = seq(1:length(mods))))
  names(counts) <- mods
  AICsumm <- list(weights,counts)
  names(AICsumm) <- c("Weights","Counts")
  return(AICsumm)
}


