#file1 <- read.csv("OutputsSP.csv")
#file2 <- read.csv("OutputsSP2.csv")
#file3 <- read.csv("OutputsSP3.csv")

#data <- rbind(rbind(file1, file2), file3)

data <- read.csv("outputsBIG1.csv")#[1:8,]
for(csvIndex in 2:16)
{
  file <- read.csv(paste0("outputsBIG", csvIndex, ".csv"))#[1:8,]
  data <- rbind(data, file)
}

columns <- names(data)
columns.ols <- columns[grepl(".ols", columns)]
columns.sp <- columns[grepl(".sp", columns)]
columns.aic <- columns[grepl(".aic", columns)]

data["ref.ols"] <- 0
data["ref.sp"] <- 0
data["ref.aic"] <- 0
data["ref.improve"] <- 0

for(rowIndex in 1:nrow(data))
{
  originator.name <- gsub(" ", ".", data[rowIndex, "originator"])
  originator.name.ols <- paste0(originator.name, ".ols")
  originator.name.sp <- paste0(originator.name, ".sp")
  originator.name.aic <- paste0(originator.name, ".aic")
  
  originator.value.ols <- as.numeric(data[rowIndex, originator.name.ols])
  originator.value.sp <- as.numeric(data[rowIndex, originator.name.sp])
  originator.value.aic <- as.numeric(data[rowIndex, originator.name.aic])
  
  values.ols <- as.numeric(data[rowIndex, columns.ols])
  values.sp <- as.numeric(data[rowIndex, columns.sp])
  values.aic <- as.numeric(data[rowIndex, columns.aic])
  
  max.ols <- max(values.ols)
  max.sp <- max(values.sp)
  min.aic <- min(values.aic)
  
  values.ols <- values.ols - max.ols
  values.sp <- values.sp - max.sp
  values.aic <- min.aic - values.aic
  
  imp.originator <- originator.value.sp - originator.value.ols 
  imp.avg <- mean(values.sp) - mean(values.ols)
  
  originator.value.ols <- originator.value.ols - max.ols
  originator.value.sp <- originator.value.sp - max.sp
  originator.value.aic <- min.aic - originator.value.aic
  
  ref.ols <- exp(originator.value.ols)/sum(exp(values.ols))
  ref.sp <- exp(originator.value.sp)/sum(exp(values.sp))
  ref.aic <- exp(originator.value.aic/2)/sum(exp(values.aic/2))
  
  data[rowIndex, "ref.ols"] <- ref.ols
  data[rowIndex, "ref.sp"] <- ref.sp
  data[rowIndex, "ref.aic"] <- ref.aic
  data[rowIndex, "ref.improve"] <- (ref.sp - ref.ols)
  data[rowIndex, "imp.originator"] <- imp.originator
  data[rowIndex, "imp.avg"] <- imp.avg
}

originators <- unique(data[,"originator"])

aic.means <- numeric(length(originators))
names(aic.means) <- originators

aic.sds <- numeric(length(originators))
names(aic.sds) <- originators

ols.means <- numeric(length(originators))
names(ols.means) <- originators

ols.sds <- numeric(length(originators))
names(ols.sds) <- originators

sp.means <- numeric(length(originators))
names(sp.means) <- originators

sp.sds <- numeric(length(originators))
names(sp.sds) <- numeric(length(originators))

imps.means <- numeric(length(originators))
names(imps.means) <- originators

imps.sds <- numeric(length(originators))
names(imps.sds) <- originators

for(originatorIndex in 1:length(originators))
{
  originator <- originators[originatorIndex]
  rows <- data[data[,"originator"] == originator,]
  
  ols.refs <- rows$ref.ols
  sp.refs <- rows$ref.sp
  aic.refs <- rows$ref.aic
  imps.refs <- rows$ref.improve
  
  aic.means[originator] <- mean(aic.refs)
  aic.sds[originator] <- sd(aic.refs)
  
  ols.means[originator] <- mean(ols.refs)
  ols.sds[originator] <- sd(ols.refs)
  
  sp.means[originator] <- mean(sp.refs)
  sp.sds[originator] <- sd(sp.refs)
  
  imps.means[originator] <- mean(imps.refs)
  imps.sds[originator] <- sd(imps.refs)
}
print("aic means")
print(aic.means)
print("ols means")
print(ols.means)
print("sp means")
print(sp.means)
print("improvement means")
print(imps.means)
print("improvement sds")
print(imps.sds)
#par(mar=c(1,1,1,1))
print(paste("mean ols ref:", mean(data$ref.ols)))
print(paste("mean sp ref:", mean(data$ref.sp)))
print(paste("mean aic ref:", mean(data$ref.aic)))
print(t.test(data$ref.improve))
print(paste("mean sp vs. ols ref improvement:", mean(data$ref.improve)))
print(paste("std err sp vs ols ref improvement:", sd(data$ref.improve)/sqrt(nrow(data))))
print(paste("std err aic ref:", sd(data$ref.aic)/sqrt(nrow(data))))
#print(nrow(data[data$ref.improve < 0,])/nrow(data))
#print(median(data$ref.improve))
#hist(data$ref.improve)

#hist(data$imp.originator)
#hist(data$imp.avg)
print(t.test(data$imp.originator - data$imp.avg))
#print(nrow(data[which(data$imp.originator > 0),]))
#print(nrow(data))
#print(t.test(data$imp.originator))
#print(sd(data$imp.originator)/sqrt(nrow(data)))