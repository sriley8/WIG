getBestFit.SP <- function(all.data, fn, defaultParams, use.actual)
{
  best.fit <- optimize(fnWrapper.SP, 
                       fn = fn,
                       defaultParams = defaultParams,
                       all.data = all.data,
                       use.actual = use.actual,
                       interval = c(0, 1),
                       maximum = TRUE,
                       tol = .0001)
  
  return (c(s = best.fit$maximum, cv = best.fit$objective))
}

fnWrapper.SP <- function(s, fn, defaultParams, all.data, use.actual)
{
  sum <- 0
  for(strataIndex in unique(all.data$strata))
  {
    data.train <- all.data[which(all.data$strata != strataIndex),]
    
    for(id in unique(data.train$Subject.ID))
    {
      data.train[which(data.train$Subject.ID == id),"weight"] <- 1
      data.train[which(data.train$Subject.ID != id),"weight"] <- s
      data.test <- all.data[which(all.data$strata == strataIndex & all.data$Subject.ID == id),]
      
      if(nrow(data.test) > 0)
      {
        fit <- getBestWeightedFit.OLS(fn, data.train, defaultParams, use.actual, all.data)
        cat(".")
        
        for(rowIndex in 1:nrow(data.test))
        {
          row <- data.test[rowIndex,]
          prediction <- fn(fit, row)
          
          if(use.actual)
          {
            output <- log(row$B1/row$B2)
          }
          else
          {
            output <- row$generated
          }
          
          loglik <- dnorm(output, prediction, fit["std.err"], log = TRUE)
          sum <- sum + loglik
        }
      }
    }
  }
  
  print(paste("s = ", s, "cv = ", sum))
  return(sum)
}