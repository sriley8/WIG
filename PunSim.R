source("Functions.R")
source("OLSFits.R")
source("MLMFits.R")
source("CV.R")
source("SPFits.R")
library(parallel)

sanityCheck <- function(file.data, OLS.fits)
{
  for(fGenIndex in 1:numFunctions)
  {
    baseVals <- OLS.fits[[fGenIndex]]
    fnGen <- getFunction(fGenIndex)
    fnGenName <- fnList[[fGenIndex]]
    dataWithPredictions <- addPredictions(file.data, fnGen, baseVals)
    
    dataGen <- generateDataSet(fnGen, baseVals, stdErrMultiplier = 1.0, dataWithPredictions)
    dataGen <- setStrata(dataGen)
    
    for(fFitIndex in 1:numFunctions)
    {
      print(paste("fitting", fnList[[fFitIndex]]))
      fnFit <- getFunction(fFitIndex)
      
      cv.ols <- getCV.OLS(dataGen, fnFit, getDefaultParams(fFitIndex), FALSE)
      print(paste("cv.ols =", cv.ols))
      
      fit.sp <- fnWrapper.SP(0, fnFit, getDefaultParams(fFitIndex), dataGen, FALSE)
      print(paste("cv.sp =", fit.sp["cv"], "s =", fit.sp["s"]))
    }
  }
}

AICSimulation <- function(file.data, OLS.fits)
{
  nDataSets <- 1000
  headerWritten <- FALSE
  outfile <- "outputsAIC.csv"
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
      
      outputLine <- list(originator = fnGenName)
      
      for(fFitIndex in 1:numFunctions)
      {
        aic <- get.aic(fFitIndex, dataGen, use.actual = FALSE, use.AICc = TRUE)
        print(paste(fnList[[fFitIndex]], "AICc =", aic))
        
        aicName <- paste(fnList[[fFitIndex]], "aic", sep=".")
        outList <- list(aic)
        names(outList) <- list(aicName)
        
        outputLine <- append(outputLine, outList)
      }
      
      cat("\n")
      
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

SP.worker.fn <- function(thread.id, file.data, OLS.fits, num.threads)
{
  source("Functions.R")
  source("OLSFits.R")
  source("MLMFits.R")
  source("CV.R")
  source("SPFits.R")
  set.seed(thread.id)
  sink(paste0("log", thread.id, ".txt"), append=FALSE)
  
  nDataSets <- ceiling(1200/(num.threads * 8))
  headerWritten <- FALSE
  outfile <- paste0("outputsBIG", thread.id, ".csv")
  for(dataSetIndex in 1:nDataSets)
  {
    for(fGenIndex in 1:numFunctions)
    {
      baseVals <- OLS.fits[[fGenIndex]]
      fnGen <- getFunction(fGenIndex)
      fnGenName <- fnList[[fGenIndex]]
      dataWithPredictions <- addPredictions(file.data, fnGen, baseVals)
      
      print(paste("generating data set", dataSetIndex, "of", nDataSets, "from", fnGenName, "on thread", thread.id))
      dataGen <- generateDataSet(fnGen, baseVals, stdErrMultiplier = 1.0, dataWithPredictions)
      dataGen <- setStrata(dataGen)
      
      outputLine <- list(originator = fnGenName)
      
      for(fFitIndex in 1:numFunctions)
      {
        print(paste("fitting", fnList[[fFitIndex]]))
        fnFit <- getFunction(fFitIndex)
        
        aic <- get.aic(fFitIndex, dataGen, use.actual = FALSE, use.AICc = TRUE)
        print(paste("aic =", aic))
        
        cv.ols <- getCV.OLS(dataGen, fnFit, getDefaultParams(fFitIndex), FALSE)
        print(paste("cv.ols =", cv.ols))
        
        fit.sp <- getBestFit.SP(dataGen, fnFit, getDefaultParams(fFitIndex), FALSE)
        print(paste("cv.sp =", fit.sp["cv"], "s =", fit.sp["s"]))
        
        aicName <- paste(fnList[[fFitIndex]], "aic", sep=".")
        olsName <- paste(fnList[[fFitIndex]], "ols", sep=".")
        spName <- paste(fnList[[fFitIndex]], "sp", sep=".")
        outList <- list(aic, cv.ols, fit.sp["cv"])
        names(outList) <- list(aicName, olsName, spName)
        
        outputLine <- append(outputLine, outList)
      }
      
      cat("\n")
      
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

SPSimulation <- function(file.data, OLS.fits)
{
  num.threads <- 16
  cluster <- makeCluster(num.threads, outfile="log.txt")
  parLapply(cluster, 1:num.threads, SP.worker.fn, 
            file.data = file.data, 
            OLS.fits = OLS.fits, 
            num.threads = num.threads)
  stopCluster(cluster)
}

SPAnalysis <- function(file.data)
{
  data.with.strata <- setStrata(file.data)
  
  outputLine <- list(originator = "?")
  
  for(fFitIndex in 1:numFunctions)
  {
    print(paste("fitting", fnList[[fFitIndex]]))
    fnFit <- getFunction(fFitIndex)
    
    cv.ols <- getCV.OLS(data.with.strata, fnFit, getDefaultParams(fFitIndex), TRUE)
    print(paste("cv.ols =", cv.ols))
    
    fit.sp <- getBestFit.SP(data.with.strata, fnFit, getDefaultParams(fFitIndex), TRUE)
    print(paste("cv.sp =", fit.sp["cv"], "s =", fit.sp["s"]))
    
    olsName <- paste(fnList[[fFitIndex]], "ols", sep=".")
    spName <- paste(fnList[[fFitIndex]], "sp", sep=".")
    outList <- list(cv.ols, fit.sp["cv"])
    names(outList) <- list(olsName, spName)
    
    outputLine <- append(outputLine, outList)
  }
  
  cat("\n")
  outfile <- "SPAnalysisActual.csv"
  
  header <- paste(names(outputLine), collapse = ",")
  write(header, file=outfile, append=TRUE)
  
  output <- paste(outputLine, collapse = ",")
  write(output, file=outfile, append = TRUE)
}

file.data <- read.csv("PunishmentDataCritchfield.csv")
#SPAnalysis(file.data)

OLS.fits <- vector("list", numFunctions)
for(i in 1:numFunctions)
{
  OLS.fits[[i]] <- getBestFits.OLS(getFunction(i), file.data, getDefaultParams(i), TRUE, file.data)
}

#AICSimulation(file.data, OLS.fits)
SPSimulation(file.data, OLS.fits)

#sanityCheck(file.data, OLS.fits)

