fnList <- list("OML", 
             "GML",
             "Competitive Suppression",
             "Direct Suppression",
             "Generalized Competitive Suppression",
             "Generalized Direct Suppression",
             "ROE Competitive Suppression",
             "ROE Direct Suppression")

numFunctions <- length(fnList)

addPredictions <- function(file.data, fn, baseVals)
{
  ret <- file.data
  
  ret$predicted <- 0
  ret$actual <- 0
  
  for(rowIndex in 1:nrow(file.data))
  {
    row <- file.data[rowIndex,]
    id <- as.numeric(row["Subject.ID"])
    id.vals <- baseVals[which(baseVals$id == id), ]
    predicted <- as.numeric(fn(id.vals, row))
    ret[rowIndex, "predicted"] <- predicted
    ret[rowIndex, "actual"] <- log(row$B1/row$B2)
  }
  
  return (ret)
}

generateDataSet <- function(baseFunction, baseVals, stdErrMultiplier, file.data)
{
  ret <- file.data
  ret$generated <- 0
  
  for(rowIndex in 1:nrow(file.data))
  {
    row <- file.data[rowIndex,]
    id <- as.numeric(row["Subject.ID"])
    id.vals <- baseVals[which(baseVals$id == id), ]
    predicted <- as.numeric(row["predicted"])
    std.err.gen <- as.numeric(id.vals["std.err"]) * stdErrMultiplier
    generated <- rnorm(1, predicted, std.err.gen)
    ret[rowIndex, "generated"] <- generated
  }
  
  return (ret)
}

get.aic <- function(fFitIndex, data, use.actual, use.AICc)
{
  aicSum <- 0
  for(id in unique(data$Subject.ID))
  {
    aicSum <- aicSum + get.subject.aic(fFitIndex, data, use.actual, use.AICc, id)
  }
  return (aicSum)
}

get.subject.aic <- function(fFitIndex, all.data, use.actual, use.AICc, id)
{
  logliksum <- 0
  
  fnFit <- getFunction(fFitIndex)
  id.data <- all.data[which(all.data$Subject.ID == id),]
  best.fits <- getBestFits.OLS(fnFit, id.data, getDefaultParams(fFitIndex), FALSE, all.data)
  
  for(rowIndex in 1:nrow(id.data))
  {
    row <- id.data[rowIndex,]
    id.vals <- best.fits[which(best.fits$id == id), ]
    std.err <- id.vals$std.err
    predicted <- as.numeric(fnFit(id.vals, row))
    
    val <- 0
    if(use.actual) val <- log(row$B1/row$B2)
    else val <- row$generated
    
    loglik <- dnorm(val, mean = predicted, sd = std.err, log = TRUE)
    logliksum <- logliksum + loglik
  }
  
  k <- getNumParams(fnFit)
  aic <- 2*k - 2*logliksum
  
  if(use.AICc)
  {
    n <- nrow(id.data)
    AICc <- aic + 2*k*(k+1)/(n - k - 1)
    return(AICc)
  }
  else
  {
    return (aic)
  }
}

getNumParams <- function(fn)
{
  if (identical(fn, OML.pred) || identical(fn, CompSupp.pred) || identical(fn, DirSupp.pred)) return (0) 
  else if (identical(fn, GML.pred) || identical(fn, GenCompSupp.pred) || identical(fn, GenDirSupp.pred)) return (2)
  else return (1)
}

getParamNames <- function(fn)
{
  if (identical(fn, OML.pred) || identical(fn, CompSupp.pred) || identical(fn, DirSupp.pred)) return (c()) 
  else if (identical(fn, GML.pred) || identical(fn, GenCompSupp.pred) || identical(fn, GenDirSupp.pred)) return (c("a", "log.b"))
  else return (c("c"))
}

getUpper <- function(fn, id.data)
{
  if (identical(fn, OML.pred) || identical(fn, CompSupp.pred) || identical(fn, DirSupp.pred)) return (c()) 
  else if (identical(fn, GML.pred) || identical(fn, GenCompSupp.pred) || identical(fn, GenDirSupp.pred)) return (c(a = Inf, log.b = Inf))
  else if(identical(fn, ROECompSupp.pred)) return (c(c = Inf))
  else
  {
    return(c(c=getMaxC(id.data)))
  }
}

getLower <- function(fn)
{
  if (identical(fn, OML.pred) || identical(fn, CompSupp.pred) || identical(fn, DirSupp.pred)) return (c()) 
  else if (identical(fn, GML.pred) || identical(fn, GenCompSupp.pred) || identical(fn, GenDirSupp.pred)) return (c(a = -Inf, log.b = -Inf))
  else return(c(c = 0))
}

getDefaultParams <- function(index)
{
  if(index == 1 || index == 3 || index == 4) return(c())
  else if(index == 2 || index == 5 || index == 6) return(c(a = 1, log.b = 0))
  else return(c(c = 1))
}

getDefaultSigmas <- function(index)
{
  if(index == 1 || index == 3 || index == 4) return(c())
  else if(index == 2 || index == 5 || index == 6) return(c(C.a = 1, C.log.b = 1))
  else return(c(C.c = 1))
}

getDefaultHyperParams <- function(index)
{
  if(index == 1 || index == 3 || index == 4) return(c())
  else if(index == 2 || index == 5 || index == 6) return(c(a.mean = 1, a.sd = 1, log.b.mean = 0, log.b.sd = 1))
  else return(c(c.mean = 1, c.sd = 1))
}

getSSE <- function(par, passedFn, data, use.actual, use.weights)
{
  sumSq <- 0
  for(i in 1:nrow(data)) 
  {
    row <- data[i,]
    pred <- passedFn(par, row)
    val <- 0
    if(use.actual) val <- log(row$B1/row$B2)
    else val <- row$generated
    resid.sq <- (pred - val)^2
  
    if(use.weights)
    {
      sumSq <- sumSq + row$weight*resid.sq   
    }
    else
    {
      sumSq <- sumSq + resid.sq
    }
  }
  
  return (sumSq)
}

getMax <- function(parName)
{
  return (Inf)
}

getMin <- function(parName)
{
  if(parName == "c") return (0)
  return (-Inf)
}

getMaxC <- function(file.data)
{
  ret <- 10000
  for(i in 1:nrow(file.data))
  {
    row <- file.data[i,]
    
    upperMax <- row$r1/row$p1
    if(row$p1 > 0 && upperMax < ret)
    {
      ret <- upperMax
    }
    
    lowerMax <- row$r2/row$p2
    if(row$p2 > 0 && lowerMax < ret)
    {
      ret <- lowerMax
    }
  }
  return (ret)
}

getFunction <- function(index)
{
  if (index == 1) OML.pred 
  else if (index == 2) GML.pred
  else if (index == 3) CompSupp.pred
  else if (index == 4) DirSupp.pred
  else if (index == 5) GenCompSupp.pred
  else if (index == 6) GenDirSupp.pred
  else if (index == 7) ROECompSupp.pred
  else if (index == 8) ROEDirSupp.pred
}

OML.pred <- function(par, row)
{
  log.R.ratio <- log(row$r1/row$r2)
  return(log.R.ratio)
}

GML.pred <- function(par, row)
{
  log.R.ratio <- log(row$r1/row$r2)
  ret <- par["log.b"] + par["a"] * log.R.ratio
  return (ret)
}

CompSupp.pred <- function(par, row)
{
  num <- row$r1 + row$p2
  den <- row$r2 + row$p1
  return(log(num/den))
}

DirSupp.pred <- function(par, row)
{
  num <- row$r1 - row$p1
  den <- row$r2 - row$p2
  return(log(num/den))
}

GenCompSupp.pred <- function(par, row)
{
  num <- row$r1 + row$p2
  den <- row$r2 + row$p1
  log.ratio <- log(num/den)
  return(par["log.b"] + par["a"] * log.ratio)
}

GenDirSupp.pred <- function(par, row)
{
  num <- row$r1 - row$p1
  den <- row$r2 - row$p2
  log.ratio <- log(num/den)
  return(par["log.b"] + par["a"] * log.ratio)
}

ROECompSupp.pred <- function(par, row)
{
  if(length(par) > 1)
  {
    browser()
    par <- par["c"]
  }
  
  num <- row$r1 + par * row$p2
  den <- row$r2 + par * row$p1
  ret <- as.numeric(log(num/den))
  if(!is.finite(ret))
  {
    browser()
  }
  return(ret)
}

ROEDirSupp.pred <- function(par, row)
{
  if(length(par) > 1)
  {
    browser()
    par <- par["c"]
  }
  
  num <- row$r1 - par * row$p1
  den <- row$r2 - par * row$p2
  ret <- as.numeric(log(num/den))
  if(!is.finite(ret))
  {
    browser()
  }
  return(ret)
}