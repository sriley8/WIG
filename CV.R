num.strata <- -1

setStrata <- function(data)
{
  data$strata <- 0
  
  if(num.strata > 0)
  {
    #Stratified CV
    set.seed(0)
    ids <- unique(data$Subject.ID)
    ret <- data.frame(matrix(nrow = 0, ncol = ncol(data)))
    colnames(ret) <- colnames(data)
    
    for(id in ids)
    {
      rows.for.id <- data[which(data$Subject.ID == id), ]
      
      currentStrataIndex <- num.strata + 1
      strata.in.random.order <- 0
      rows.in.random.order <- sample(nrow(rows.for.id))
      for(rowIndex in rows.in.random.order)
      {
        if(currentStrataIndex > num.strata)
        {
          strata.in.random.order <- sample(num.strata)
          currentStrataIndex <- 1
        }
        
        rows.for.id[rowIndex, "strata"] <- strata.in.random.order[currentStrataIndex]
        currentStrataIndex <- currentStrataIndex + 1
      }
      
      ret <- rbind(ret, rows.for.id)
    }
    
    set.seed(Sys.time())
    
    return (ret)
  }
  else
  {
    #LOO CV
    for(rowIndex in 1:nrow(data))
    {
      data[rowIndex, "strata"] <- rowIndex
    }
    
    return (data)
  }
}