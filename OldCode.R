CBestSimulation <- function(file.data, OLS.fits)
{
  nDataSets <- 50
  headerWritten <- FALSE
  outfile <- "outputsCBest.csv"
  for(dataSetIndex in 1:nDataSets)
  {
    for(fGenIndex in 1:numFunctions)
    {
      baseVals <- OLS.fits[[fGenIndex]]
      fnGen <- getFunction(fGenIndex)
      fnGenName <- fnList[[fGenIndex]]
      dataWithPredictions <- addPredictions(file.data, fnGen, baseVals)
      
      print(paste("generating data set", dataSetIndex, "of", nDataSets, "from", fnGenName))
      dataGen <- generateDataSet(fnGen, baseVals, stdErrMultiplier = 1.0, dataWithPredictions)
      dataGen <- setStrata(dataGen)
      
      outputLine <- list(originator = fnGenName)
      for(fFitIndex in 1:numFunctions)
      {
        print(paste("fitting", fnList[[fFitIndex]]))
        fnFit <- getFunction(fFitIndex)
        
        cv.ols <- getCV.OLS(dataGen, getFunction(fFitIndex), getDefaultParams(fFitIndex), use.actual = FALSE)
        print(paste(fnList[[fFitIndex]], "cv.ols =", cv.ols))
        
        ols.fits <- getBestFits.OLS(fnFit, dataGen, getDefaultParams(fFitIndex), FALSE, dataGen)
        
        method <- "Nelder-Mead"
        if(length(getDefaultParams(fFitIndex) == 1))
        {
          method <- "Brent" #1-D optimization, Brent is recommended over Nelder-Mead
        }
        
        mlmFit <- optim(par = getDefaultSigmas(fFitIndex),
                        fn = CBestWrapper,
                        method = method,
                        control = list(fnscale = -1, reltol = .00001),
                        fnIndex = fFitIndex,
                        dataGen = dataGen,
                        ols.fits = ols.fits)
        cv.mlm <- mlmFit$value
        print(paste(fnList[[fFitIndex]], "cv.mlm =", cv.mlm))
        
        olsName <- paste(fnList[[fFitIndex]], "ols", sep=".")
        mlmName <- paste(fnList[[fFitIndex]], "mlm", sep=".")
        outList <- list(cv.ols, cv.mlm)
        names(outList) <- list(olsName, mlmName)
        
        outputLine <- append(outputLine, outList)
      }
      
      if(!headerWritten)
      {
        header <- paste(names(outputLine), collapse = ",")
        write(header, file=outfile, append=TRUE)
        headerWritten <- TRUE
      }
      
      output <- paste(outputLine, collapse = ",")
      write(output, file=outfile, append = TRUE)
    }
  }
}

CBestWrapper <- function(par, fnIndex, dataGen, ols.fits)
{
  if(any(par <= 0))
  {
    return(NA)
  }
  
  print(paste(par, collapse=" "))
  ret <- getCV.MLM(dataGen, 
                   getFunction(fnIndex), 
                   getDefaultParams(fnIndex), 
                   getDefaultHyperParams(fnIndex), 
                   use.actual = FALSE, 
                   sigma.values = par,
                   ols.fits = ols.fits)
  cat(paste(" CVLL:", ret, "\n"))
  
  return (ret)
}

CVariationSimulation <- function(file.data, OLS.fits)
{
  dataSetsPerC <- 50
  nCs <- 41
  maxC <- 4
  cChunk <- maxC/(nCs - 1)
  
  headerWritten <- FALSE
  outfile <- "outputsCVariation_GCS.csv"
  for(fIndex in 5:5)
  {
    baseVals <- OLS.fits[[fIndex]]
    fn <- getFunction(fIndex)
    dataWithPredictions <- addPredictions(file.data, fn, baseVals)
    for(cIndex in 23:nCs)
    {
      C <- cChunk * (cIndex - 1)
      for(dataSetIndex in 1:dataSetsPerC)
      {
        print("")
        print(paste("generating data set", dataSetIndex, "of", dataSetsPerC, "from", fnList[[fIndex]], "fitting with C =", C))
        data <- generateDataSet(fn, baseVals, stdErrMultiplier = 1.0, dataWithPredictions)
        CV.vals <- fitAllFunctions(data, C, fit.OLS = FALSE)
        
        line <- append(list(originator=fnList[[fIndex]], C=C), CV.vals)
        
        if(!headerWritten)
        {
          header <- paste(names(line), collapse = ",")
          write(header, file=outfile, append=TRUE)
          headerWritten <- TRUE
        }
        
        output <- paste(line, collapse = ",")
        write(output, file=outfile, append = TRUE)
      }
    }
  }
}

SEVariationSimulation <- function(file.data, OLS.fits)
{
  dataSetsPerStdErr <- 10
  nStdErrs <- 21
  maxStdErrMultiplier <- 2
  stdErrChunk <- maxStdErrMultiplier/(nStdErrs - 1)
  
  headerWritten <- FALSE
  for(fIndex in 1:numFunctions)
  {
    baseVals <- OLS.fits[[fIndex]]
    fn <- getFunction(fIndex)
    dataWithPredictions <- addPredictions(file.data, fn, baseVals)
    for(stdErrIndex in 6:nStdErrs)
    {
      stdErrMultiplier <- stdErrChunk * (stdErrIndex - 1)
      for(dataSetIndex in 1:dataSetsPerStdErr)
      {
        print(paste("generating from", fnList[[fIndex]], "stErrMultiplier =", stdErrMultiplier))
        data <- generateDataSet(fn, baseVals, stdErrMultiplier, dataWithPredictions)
        CV.vals <- fitAllFunctions(data, C = 1.0, fit.OLS = TRUE)
        line <- append(list(originator=fnList[[fIndex]], stdErrMultiple=stdErrMultiplier), CV.vals)
        
        if(!headerWritten)
        {
          header <- paste(names(line), collapse = ",")
          write(header, file="outputs.csv", append=TRUE)
          headerWritten <- TRUE
        }
        
        output <- paste(line, collapse = ",")
        write(output, file="outputs.csv", append = TRUE)
      }
    }
  }
}

getBestFit2.SP <- function(all.data, fn, defaultParams, use.actual)
{
  ind.fits.by.strata <- list()
  group.fits.by.strata <- list()
  for(strataIndex in 1:numStrata)
  {
    cat(strataIndex)
    data.train <- all.data[which(all.data$strata != strataIndex), ]
    ind.fits.by.strata[[strataIndex]] <- getBestFits.OLS(fn, data.train, defaultParams, use.actual, all.data)
    group.fits.by.strata[[strataIndex]] <- getBestGroupFit.OLS(fn, data.train, defaultParams, use.actual, all.data)
  }
  
  best.fit <- optimize(fnWrapper,
                       fn = fn,
                       ind.fits.by.strata = ind.fits.by.strata,
                       group.fits.by.strata = group.fits.by.strata,
                       defaultParams = defaultParams,
                       all.data = all.data,
                       use.actual = use.actual,
                       interval = c(0, 1),
                       maximum = TRUE)
  
  return (c(s = best.fit$maximum, cv = best.fit$objective))
}

fnWrapper2.SP <- function(s, fn, ind.fits.by.strata, group.fits.by.strata, defaultParams, all.data, use.actual)
{
  sum <- 0
  for(strataIndex in 1:numStrata)
  {
    group.fit <- group.fits.by.strata[[strataIndex]]
    ind.fits <- ind.fits.by.strata[[strataIndex]]
    data.test <- all.data[which(all.data$strata == strataIndex),]
    
    for(subjectIndex in 1:nrow(ind.fits))
    {
      id <- ind.fits[subjectIndex, "id"]
      
      blended.params <- defaultParams
      if(!is.null(defaultParams))
      {
        for(paramIndex in 1:length(defaultParams))
        {
          paramName <- names(defaultParams)[paramIndex]
          ind.param.value <- ind.fits[subjectIndex, paramName]
          group.param.value <- group.fit[paramName]
          
          blended.value <- s*group.param.value + (1 - s)*ind.param.value
          blended.params[paramName] <- blended.value
        }
      }
      
      ind.std.err <- ind.fits[subjectIndex, "std.err"]
      group.std.err <- group.fit["std.err"]
      blended.std.err <- sqrt(s*group.std.err^2 + (1 - s)*ind.std.err^2)
      
      id.test.data <- data.test[which(data.test$Subject.ID == id),]
      for(rowIndex in 1:nrow(id.test.data))
      {
        row <- id.test.data[rowIndex,]
        prediction <- fn(blended.params, row)
        
        if(use.actual)
        {
          output <- row$actual
        }
        else
        {
          output <- row$generated
        }
        
        loglik <- dnorm(output, prediction, blended.std.err, log = TRUE)
        sum <- sum + loglik
      }
    }
  }
  
  return (sum)
}