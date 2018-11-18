maxEMIter <- 200
bestHyperpars <- list()

getBestFits.MLM <- function(passedFn, train.data, par.init, hyperpar.init, use.actual, all.data, sigma.values, ols.fits, strata.index)
{
  pars <- buildPars.MLM(passedFn, train.data, par.init)
  if(!is.null(hyperpar.init))
  {
    if(length(bestHyperpars) >= strata.index && identical(names(bestHyperpars[[strata.index]]), names(hyperpar.init)))
    {
      hyperpar <- bestHyperpars[[strata.index]]
    }
    else
    {
      hyperpar <- initHyperpars.MLM(hyperpar.init, ols.fits, par.init)
    }
  }
  else
  {
    hyperpar <- NULL
  }
  
  parNames <- names(pars)[-1]
  
  breakFlag <- FALSE
  iter <- 0
  while(!is.null(par.init) && !breakFlag && iter < maxEMIter)
  {
    #print(hyperpar)
    
    iter <- iter + 1
    breakFlag <- TRUE
    #M step
    for(parsRowIndex in 1:nrow(pars))
    {
      parRow <- pars[parsRowIndex, ]
      id <- parRow$id
      id.data <- train.data[which(train.data$Subject.ID == id),]
      par <- as.double(parRow[-1])
      names(par) <- parNames
      
      oldPar <- par
      
      upper <- c(a = Inf, log.b = Inf, c = Inf)
      lower <- c(a = 0, log.b = -Inf, c = 0)
      if(identical(passedFn, ROEDirSupp.pred))
      {
        #must limit the c value here - large c will cause negative results
        max.c = getMaxC(all.data[which(all.data$Subject.ID == id),])
        upper <- list(c = max.c)
        lower <- list(c = 0)
        par["c"] <- .001
      }
      
      fit <- optim(par = par, 
                   fn = fnWrapper.MLM,
                   lower = lower,
                   upper = upper,
                   method = "L-BFGS-B",
                   passedFn = passedFn,
                   id.data = id.data, 
                   use.actual = use.actual, 
                   hyperpar = hyperpar,
                   sigma.values = sigma.values)
      
      newPar <- fit$par
      
      if(breakFlag)
      {
        diffFound <- compareDifference(oldPar, newPar)
        breakFlag <- !diffFound
      }
      
      for(parNameIndex in 1: length(names(par)))
      {
        parName <- parNames[parNameIndex]
        pars[parsRowIndex, parName] <- fit$par[parName]
      }
    }
    
    #E step
    for(parName in parNames)
    {
      parMeanName <- paste(parName, "mean", sep=".")
      
      meanVal <- mean(pars[,parName])
      
      hyperpar[parMeanName] <- meanVal
    }
    
    cat(".")
    
    if(iter == maxEMIter)
    {
      cat("#")
    }
  }
  
  pars$std.err <- 0
  for(parsRowIndex in 1:nrow(pars))
  {
    parRow <- pars[parsRowIndex, ]
    id <- parRow$id
    id.data <- train.data[which(train.data$Subject.ID == id),]
    par <- as.double(parRow[-1])
    names(par) <- parNames
    
    pars[parsRowIndex, "std.err"] <- get.std.err(par, passedFn, id.data, use.actual)
  }
  
  bestHyperpars[[strata.index]] <<- hyperpar
  return (pars)
}

tol <- .0001
compareDifference <- function(oldPar, newPar)
{
  for(i in 1:length(oldPar))
  {
    if(abs(oldPar[i] - newPar[i]) > tol) 
    {
      return (TRUE)
    }
  }

  return (FALSE)
}

initHyperpars.MLM <- function(hyperpars, ols.fits, par.init)
{
  if(length(par.init > 0))
  {
    means <- colMeans(ols.fits)
    for(parIndex in 1:length(par.init))
    {
      parName <- names(par.init[parIndex])
      hyperparName <- paste(parName, "mean", sep=".")
      hyperparValue <- means[parName]
      
      hyperpars[hyperparName] <- hyperparValue
    }
  }
  
  return (hyperpars)
}

buildPars.MLM <- function(passedFn, file.data, par.init, hyperpar.init)
{
  ids <- unique(file.data$Subject.ID)
  pars <- data.frame(matrix(ncol = 1 + length(par.init), nrow = length(ids)))
  
  names(pars) <- c("id", names(par.init))
  pars$id <- unique(file.data$Subject.ID)
  
  if(!is.null(par.init))
  {
    for(i in 1:length(par.init))
    {
      parName <- names(par.init)[i]
      parValue <- par.init[[i]]
      pars[parName] <- parValue
    }
    
    if(identical(passedFn, ROEDirSupp.pred))
    {
      for(id in unique(file.data$Subject.ID))
      {
        id.data <- file.data[which(file.data$Subject.ID == id), ]
        #must limit the c value here - large c will cause negative results
        max.c <- getMaxC(id.data)
        pars[which(pars$id == id), "c"] = max.c/2
      }
    }
  }
  
  return(pars)
}

fnWrapper.MLM <- function(par, passedFn, id.data, use.actual, hyperpar, sigma.values)
{
  SSE <- getSSE(par, passedFn, id.data, use.actual)
  sum <- 0
  
  paramNames <- getParamNames(passedFn)
  if(!is.null(paramNames))
  {
    for(paramIndex in 1:length(paramNames))
    {
      paramName <- paramNames[paramIndex]
      paramMeanName <- paste(paramName, "mean", sep=".")
      paramCName <- paste("C", paramName, sep=".")
      
      paramMean <- hyperpar[paramMeanName]
      sigma <- sigma.values[paramCName]
      paramC <- 1/(sigma*sigma)
      
      paramVal <- par[paramName]
      resid <- paramVal - paramMean
      sum <- sum + paramC*resid*resid
    }
  }
  
  if(!is.finite(sum+SSE))
  {
    browser()
  }
  
  return (sum + SSE)
}

get.std.err <- function(par, passedFn, id.data, use.actual)
{
  SSE <- getSSE(par, passedFn, id.data, use.actual)
  std.err <- sqrt(SSE/(nrow(id.data) - getNumParams(passedFn)))
  
  return(std.err)
}

getCV.MLM <- function(data, fn, defaultParams, defaultHyperParams, use.actual, sigma.values, ols.fits)
{
  sum <- 0
  
  for(strataIndex in 1:numStrata)
  {
    data.train <- data[which(data$strata != strataIndex),]
    holdoutRows <- data[which(data$strata == strataIndex),]
    fits <- getBestFits.MLM(fn, data.train, defaultParams, defaultHyperParams, use.actual, data, sigma.values, ols.fits, strataIndex)
    
    for(holdoutIndex in 1:nrow(holdoutRows))
    {
      holdout <- holdoutRows[holdoutIndex,]
      id <- holdout$Subject.ID
      id.vals <- fits[which(fits$id == id), ]
      std.err <- id.vals$std.err
      predicted <- as.numeric(fn(id.vals, holdout))
      
      val <- 0
      if(use.actual) val <- log(holdout$B1/holdout$B2)
      else val <- holdout$generated
      
      loglik <- dnorm(val, mean = predicted, sd = std.err, log = TRUE)
      
      sum <- sum + loglik
    }
    
    cat(strataIndex)
  }
  
  return (sum)
}