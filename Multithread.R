library(parallel)

test <- function(thread.id, passedData)
{
  output <- passedData
  write(output, file=paste0("testOutput", thread.id, ".txt"))
}

n.threads <- 16
cluster <- makeCluster(n.threads, outfile="log.txt")
parLapply(cluster, 1:n.threads, test, passedData="Hoihoi")
stopCluster(cluster)