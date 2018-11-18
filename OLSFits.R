#NOTE - all.data is only used to limit c for the ROE Direct Suppression function
doFit <- function(passedFn, train.data, par.init, use.actual, all.data, fnWrapper)
{
  numParams <- getNumParams(passedFn)
  if(numParams == 1)
  {
    if(identical(passedFn, ROEDirSupp.pred))
    {
      maxC <- getMaxC(all.data)    
    }
    else
    {
      maxC <- 100
    }
    
    fit <- optimize(fnWrapper,
                    interval = c(0, maxC),
                    passedFn = passedFn,
                    train.data = train.data,
                    use.actual = use.actual,
                    maximum = TRUE)
    ret <- list(par = c(c=fit$maximum), value = fit$objective)
    
    return (ret)
  }
  else
  {
    fit <- optim(par = par.init,
                 fn = fnWrapper,
                 method = "Nelder-Mead",
                 passedFn = passedFn,
                 train.data = train.data,
                 use.actual = use.actual)
    
    return(fit)
  }
}

getBestWeightedFit.OLS <- function(passedFn, train.data, par.init, use.actual, all.data)
{
  fit <- doFit(passedFn, train.data, par.init, use.actual, all.data, fnWrapper.WLS)
  
  weight.sum <- sum(train.data$weight)
  weight.sq.sum <- sum(train.data$weight^2)
  std.err <- sqrt(fit$value/(weight.sum - length(par.init)*weight.sq.sum/weight.sum))
  ret <- c(fit$par, std.err)
  names(ret) <- c(names(par.init), "std.err")
  return(ret)
}

getBestGroupFit.OLS <- function(passedFn, train.data, par.init, use.actual, all.data)
{
  fit <- doFit(passedFn, train.data, par.init, use.actual, all.data, fnWrapper.OLS)
  
  std.err <- sqrt(fit$value/(nrow(train.data) - length(par.init)))
  ret <- c(fit$par, std.err)
  names(ret) <- c(names(par.init), "std.err")
  return(ret)
}

getBestFits.OLS <- function(fn, train.data, par.init, use.actual, all.data)
{
  ids <- train.data[,"Subject.ID"]
  ret <- data.frame()
  
  for(id in unique(ids))
  {
    id.data <- train.data[which(train.data$Subject.ID == id), ]
    id.pars <- fitOneID(fn, id.data, par.init, use.actual, all.data, id)
    id.pars <- append(id, id.pars)
    ret <- rbind(ret, id.pars)
  }
  
  colnames(ret) <- c("id", names(par.init), "std.err")
  return(ret)
}

fitOneID <- function(passedFn, id.data, par.init, use.actual, all.data, id)
{
  fit <- doFit(passedFn, id.data, par.init, use.actual, all.data[which(all.data$Subject.ID == id),], fnWrapper.OLS)
  
  std.err <- sqrt(fit$value/(nrow(id.data) - length(par.init)))
  ret <- c(fit$par, std.err = std.err)
  return (ret)
}

fnWrapper.WLS <- function(par, passedFn, train.data, use.actual)
{
  WSSE <- getSSE(par, passedFn, train.data, use.actual, use.weights = TRUE)
  return (WSSE)
}

fnWrapper.OLS <- function(par, passedFn, train.data, use.actual)
{
  SSE <- getSSE(par, passedFn, train.data, use.actual, use.weights = FALSE)
  
  if(!is.finite(SSE))
  {
    browser()
  }
  
  return (SSE)
}

getCV.OLS <- function(data, fn, defaultParams, use.actual)
{
  sum <- 0
  
  for(strataIndex in unique(data$strata))
  {
    holdoutRows <- data[which(data$strata == strataIndex),]
    keptRows <- data[which(data$strata != strataIndex),]
    for(id in unique(holdoutRows$Subject.ID))
    {
      data.train <- keptRows[which(keptRows$Subject.ID == id),]
      fits <- getBestFits.OLS(fn, data.train, defaultParams, use.actual, data)
      
      for(holdoutIndex in 1:nrow(holdoutRows))
      {
        holdout <- holdoutRows[holdoutIndex,]
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
      cat(" ")
    }
  }

  return (sum)
}