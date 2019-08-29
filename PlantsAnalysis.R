getMeasures <- function(data){
  mean <- c()
  med <- c()
  firstQ <- c()
  thirdQ <- c()
  IQR <- c()
  sd <- c()
  var <- c()
  adv <- c()
  na <- c()
  
  # Calculate statistical measures for each attribute in the data frame
  for (i in 1:ncol(data)) {
    mean <- c(mean, mean(data[,i], na.rm = TRUE))
    med <- c(med, median(data[,i], na.rm = TRUE))
    firstQ <- c(firstQ, quantile(data[,i], 1/4, na.rm = TRUE))
    thirdQ <- c(thirdQ, quantile(data[,i], 3/4, na.rm = TRUE))
    IQR <- c(IQR, IQR(data[,i], na.rm = TRUE))
    sd <- c(sd, sd(data[,i], na.rm = TRUE))
    var <- c(var, var(data[,i], na.rm = TRUE))
    adv <- c(adv, mad(data[,i], na.rm = TRUE))
    na <- c(na, sum(is.na(data[,i])))
  }
  
  measures <- NULL
  # Create dataframe
  measures <- as.data.frame(rbind(measures, mean, med, firstQ, thirdQ, IQR, sd, var, adv, na))
  # Name the rows and columns
  colnames(measures) <- colnames(data)[1:ncol(data)]
  rownames(measures) <- c("Mean", "Median", "First Quantile", "Third Quantile", "Interquartile Range", "Standard Deviation", "Variance", "Median Absolute Deviation", "Number of NAs")
  return(measures)
}

# Save attribute histograms to PDF
saveNumericAttributeHistograms <- function(data, path){
  pdf(path, onefile = TRUE)
  for (i in 1:ncol(data)) {
    p <- ggplot(data, aes(x=data[,i])) + 
      geom_histogram(bins = 17) + labs(x=colnames(data)[i], y="Frequency") + ggtitle(paste(colnames(data)[i], "Histogram"))
    print(p)
  }
  dev.off()
}

# Save class histograms to PDF
saveClassAttributeHistogram <- function(data, path){
  pdf(path, onefile = TRUE)
  p <- ggplot(data, aes(x=Class)) + geom_bar() + labs(x="Class", y="Frequency") + ggtitle("Class Histogram")
  print(p)
  dev.off()
}

calculateMaxMinCorrelation <- function(data){
  # Get correlation matrix
  cor <- cor(data, use = "na.or.complete")
  # Replace 1s with NAs
  cor[cor == 1] <- NA
  # Get the index of the min and max values
  maxIndex <- which(abs(cor) == max(abs(cor), na.rm = TRUE), arr.ind = TRUE)
  minIndex <- which(abs(cor) == min(abs(cor), na.rm = TRUE), arr.ind = TRUE)
  
  # Print the correlation values
  print(paste0("Max cor = ", colnames(data)[maxIndex[1]], " ", colnames(data)[maxIndex[2]], " cor: ", cor[maxIndex[1],maxIndex[2]]))
  print(paste0("Min cor = ", colnames(data)[minIndex[1]], " ", colnames(data)[minIndex[2]], " cor: ", cor[minIndex[1],minIndex[2]]))
}

replaceMissingValues <- function(data, method){
  # For each column
  for (i in 1:(ncol(data)-1)) {
    # Replace with 0
    if(method == "zero"){
      data[which(is.na(data[,i])),i] <- 0
    }else{
      for (class in unique(data$Class)) {
        # Replace with mean taken from that class
        if(method == "mean"){
          data[which(is.na(data[,i]) & data$Class == class),i] <- mean(data[which(!is.na(data[,i]) & data$Class == class),i])
        }else{ # Replace with median taken from that class
          data[which(is.na(data[,i]) & data$Class == class),i] <- median(data[which(!is.na(data[,i]) & data$Class == class),i])
        }
      }
    }
  }
  return(data)
}

# Center the data by subtracting the mean
meanCenterData <- function(data){
  for (i in 1:(ncol(data)-1)) {
    data[, i] <- data[, i] - mean(data[, i])
  }
  return(data)
}

# Normalise the data
normaliseData <- function(data){
  for (i in 1:(ncol(data)-1)) {
    data[, i] <- (data[, i] - min(data[, i])) / (max(data[, i]) - min(data[, i]))
  }
  return(data)
}

# Standardise the data
standardiseData <- function(data){
  for (i in 1:(ncol(data)-1)) {
    data[, i] <- (data[, i] - mean(data[, i])) / sd(data[, i])
  }
  return(data)
}

removeRowsMultipleMissingValues <- function(data){
  rem2 <- c()
  rem3 <- c()
  # Remove rows with more than 2 NAs
  for (i in 1:nrow(data)) {
    if(sum(is.na(data[i,])) == 2){
      print(paste0("Count ", sum(is.na(data[i,])), " at ", i))
      rem2 <- c(rem2, i)
    }else if(sum(is.na(data[i,])) == 3){
      print(paste0("Count ", sum(is.na(data[i,])), " at ", i))
      rem3 <- c(rem3, i)
    }
  }
  # Do the removal
  data <- data[-rem2,]
  data <- data[-rem3,]
  return(data)
}

highCorrelationFilter <- function(data, threshold){
  cor <- abs(cor(data, use = "na.or.complete"))
  cor[cor == 1] <- NA
  cor <- as.data.frame(cor)
  # Calculate row and column mean
  cor <- rbind(cor, colMeans(cor, na.rm = TRUE))
  cor <- cbind(cor, rowMeans(cor, na.rm = TRUE))
  # Sort the matrix by the mean
  cor <- cor[order(cor[,ncol(cor)], decreasing = TRUE), order(cor[nrow(cor),], decreasing = TRUE)]
  # Remove the mean columns
  cor <- cor[ , apply(cor, 2, function(x) any(is.na(x)))]
  cor <- cor[apply(cor, 1, function(x) any(is.na(x))),]
  
  iCanCont <- TRUE
  jCanCont <- TRUE
  i <- 1
  j <- 1
  removing <- c()
  while (iCanCont == TRUE) {
    while (jCanCont == TRUE) {
      # If we can ignore the value
      if(is.na(cor[i, j]) == TRUE | cor[i, j] <= threshold){
        if(j < ncol(cor)){
          # go to the next iteration
          j <- j + 1
        }else{
          # Stop j and move onto the next i
          jCanCont <- FALSE
        }
        next
      }else{ # correlation value above the threshold
        # Calculate the mean for row i
        rowMean <- rowMeans(cor, na.rm=TRUE)
        # Calculate the mean for all the rows except j
        restOfMean <- (sum(colMeans(cor, na.rm = TRUE))/ncol(cor)) - rowMean[j]
        if(rowMean[i] > restOfMean){
          print(paste("Removing variable:", colnames(cor)[i]))
          removing <- c(removing, colnames(cor)[i])
          # remove the variable from the martix
          cor <- cor[-c(i),-c(i)]
        }
        # start over
        j <- 1
        i <- 1
      }
    }
    if(i < nrow(cor)){
      i <- i + 1
      j <- 1
      jCanCont <- TRUE
    }else{
      iCanCont <- FALSE
    }
  }
  # remove the attribues from the dataframe
  data <- data[,!(names(data) %in% removing)]
  return(data)
}

# Sort the confusion matrix by maximising the diagonal
createSortedConfusionMatrix <- function(data, method) {
  bestPerm <- c()
  bestDiag <- 0
  matrix <- table(data$Class, data[,method])
  
  # Try every permutation
  for (i in 1:5) {
    for (j in 1:5) {
      for (k in 1:5) {
        for (l in 1:5) {
          for (m in 1:5) {
            if(length(unique(c(i,j,k,l,m))) != 5){
              next
            }
            permutation <- matrix[c(i,j,k,l,m),]
            # If its better than the current best, take it
            if(sum(diag(permutation)) > bestDiag){
              bestPerm <- c(i,j,k,l,m)
              bestDiag <- sum(diag(permutation))
            }
          }
        }
      }
    }
  }
  return(matrix[bestPerm,])
}

# Get the accuracy of the clustering
getAccuracy <- function(matrix){
  sum(diag(matrix))/sum(matrix)
}

# Perform Hierarchical clustering
clusterHC <- function(data, params){
  if(!is.null(params)){
    hc <- hclust(dist(data, method = params[2]), method = params[1])
  }else{
    hc <- hclust(dist(data))
  }
  return(cutree(hc,5))
}

# Perform k-means clustering
clusterKM <- function(data, params){
  set.seed(20)
  if(!is.null(params)){
    km5 = kmeans(data,5, 30, nstart = params[1], algorithm = params[2])
  }else{
    km5 = kmeans(data,5)
  }
  return(km5$cluster)
}

# Perform PAM clustering
clusterPAM <- function(data, params){
  if(!is.null(params)){
    return(pam(data, 5, metric = params[1], stand = params[2], cluster.only = TRUE, do.swap = params[3]))
  }else{
    return(pam(data, 5, cluster.only = T))
  }
}

# Print the sorted confusion matrix and accuracy
printMatrixAndAccuracy <- function(data, method){
  cMatrix <- createSortedConfusionMatrix(data, method)
  print("Confusion matrix:")
  print(cMatrix)
  print(paste("Accuracy:", getAccuracy(cMatrix)))
}

# Tune Hierarchical clustering
tuneHC <- function(data){
  bestParam <- NA
  bestAcc <- -1
  
  # for each option get the accuracy
  for(param in c("ward.D", "ward.D2", "single", "complete", "ave", "mcq", "med", "cent")){
    for(distParam in c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
      data$HC5 <- clusterHC(data[,2:19], c(param, distParam))
      cMatrix <- createSortedConfusionMatrix(data, "HC5")
      accuracy <- getAccuracy(cMatrix)
      
      if(accuracy >= bestAcc){
        bestAcc <- accuracy
        bestParam <- c(param, distParam)
      }
    }
  }
  return(bestParam)
}

# Tune k-means
tuneKM <- function(data){
  bestParam <- NA
  bestAcc <- -1
  val <- 1
  
  # for each option get the accuracy and hill climb for nstart
  for (algorithm in c("Hart", "Lloyd", "Forgy", "MacQueen")) {
    canCont <- TRUE
    while(canCont){
      set.seed(20)
      data$KM5 = kmeans(data[2:19],5, iter.max = 50, nstart = val, algorithm = algorithm)$cluster
      cMatrix <- createSortedConfusionMatrix(data, "KM5")
      accuracy <- getAccuracy(cMatrix)
      
      if(accuracy > bestAcc){
        bestAcc <- accuracy
        bestParam <- c(val, algorithm)
        val <- val + 1
      }else{
        canCont <- FALSE
      }
    }
  }
  return(bestParam)
}

# Tune PAM
tunePAM <- function(data){
  bestAcc <- -1
  
  #for each option get the accuracy
  for (algorithm in c("euclidean", "manhattan")) {
    for (stand in c(FALSE, TRUE)) {
      for (swap in c(FALSE, TRUE)) {
        data$PAM5 = pam(data[2:19], 5, metric = algorithm, stand = stand, cluster.only = TRUE, do.swap = swap)
        cMatrix <- createSortedConfusionMatrix(data, "PAM5")
        accuracy <- getAccuracy(cMatrix)
        
        if(accuracy > bestAcc){
          bestAcc <- accuracy
          bestParam <- c(algorithm, stand, swap)
        }
      }
    }
  }
  return(bestParam)
}

# Get the raw data
plants <- read.csv("dataset.csv")

##### Question 1a #####
getMeasures(plants[,2:19])

library("ggplot2")
saveNumericAttributeHistograms(plants[,2:19], "originalDataAttHist.pdf")
saveClassAttributeHistogram(plants, "originalDataClassHist.pdf")

##### Question 1b #####
cor(plants$Orientation1, plants$Orientation7, use = "na.or.complete")
cor(plants$Mass, plants$Orientation0, use = "na.or.complete")
cor(plants$Orientation7, plants$Orientation8, use = "na.or.complete")

calculateMaxMinCorrelation(plants[,2:19])

p <- ggplot(plants, aes(x=Class, y=Orientation2)) + geom_point() + ggtitle("Orientation2 / Class Scatter Plot")
print(p)
p <- ggplot(plants, aes(x=Class, y=Depth)) + geom_point() + ggtitle("Depth / Class Scatter Plot")
print(p)
p <- ggplot(plants, aes(x=Class, y=LeafArea)) + geom_point() + ggtitle("Leaf Area / Class Scatter Plot")
print(p)

##### Question 1d #####
plants_0 <- replaceMissingValues(plants, "zero")
saveNumericAttributeHistograms(plants_0[,2:19], "0FillDataAttHist.pdf")
plants_mean <- replaceMissingValues(plants, "mean")
saveNumericAttributeHistograms(plants_mean[,2:19], "MeanFillDataAttHist.pdf")
plants_med <- replaceMissingValues(plants, "med")
saveNumericAttributeHistograms(plants_med[,2:19], "MedFillDataAttHist.pdf")

##### Question 1e #####
plants_0_centered <- meanCenterData(plants_0[,2:20])
plants_mean_centered <- meanCenterData(plants_mean[,2:20])
plants_med_centered <- meanCenterData(plants_med[,2:20])

plants_0_normalised <- normaliseData(plants_0[,2:20])
plants_mean_normalised <- normaliseData(plants_mean[,2:20])
plants_med_normalised <- normaliseData(plants_med[,2:20])

plants_0_standardised <- standardiseData(plants_0[,2:20])
plants_mean_standardised <- standardiseData(plants_mean[,2:20])
plants_med_standardised <- standardiseData(plants_med[,2:20])

print(summary(plants_med_centered))
print(summary(plants_med_normalised))
print(summary(plants_med_standardised))

##### Question 1fi #####
plants_fi <- plants
plants_fi$Sample_ID <- NULL
plants_fi$Leaf.weight <- NULL
print(paste("There are", nrow(plants_fi[rowSums(is.na(plants_fi)) > 0,]), "rows with a missing value"))
print(nrow(plants_fi[rowSums(is.na(plants_fi)) > 0,]) / nrow(plants_fi))
plants_fi <- removeRowsMultipleMissingValues(plants_fi)
print(paste("There are", nrow(plants_fi[rowSums(is.na(plants_fi)) > 0,]), "rows with a missing value"))
print(nrow(plants_fi[rowSums(is.na(plants_fi)) > 0,]) / nrow(plants_fi))
getMeasures(plants_fi[,1:17])

##### Question 1fii #####
plants_fii <- highCorrelationFilter(plants[2:19], 0.75)
plants_fii$Class <- plants$Class
plants_fii$Leaf.weight <- NULL
print(head(plants_fii))
calculateMaxMinCorrelation(plants_fii[,1:9])
print(paste("There are", nrow(plants_fii[rowSums(is.na(plants_fii)) > 0,]), "rows with a missing value"))
print(nrow(plants_fii[rowSums(is.na(plants_fii)) > 0,]) / nrow(plants_fii))

##### Question 1fiii #####
plants_cleaned <- replaceMissingValues(plants_fi, "mean")
plants_stand <- as.data.frame(scale(plants_cleaned[,1:17]))
pca <- prcomp(plants_stand,scale=T, center = T)

print(summary(pca))
screeplot(pca, type="lines",col=3, main="Variance explained by PC")
title(xlab="Principal Components")

print(pca$rotation)
biplot(pca, cex = 0.7, xlabs=rep(".", nrow(plants_stand)))
abline(h = 0, v = 0, lty = 2, col = 8)

plants_fiii <- as.data.frame(predict(pca,plants_stand)[,1:7])
plants_fiii$Class <- plants_cleaned$Class

##### Question 2a #####
print("Hierarchical clustering:")
plants_mean$HC5 <- clusterHC(plants_mean[,2:19], NULL)
cMatrix <- createSortedConfusionMatrix(plants_mean, "HC5")
print("Confusion matrix for HC clustering")
print(cMatrix)
print(paste("Accuracy:", getAccuracy(cMatrix)))

print("K-means clustering:")
plants_mean$KM5 <- clusterKM(plants_mean[,2:19], NULL)
cMatrix <- createSortedConfusionMatrix(plants_mean, "KM5")
print("Confusion matrix for KM clustering")
print(cMatrix)
print(paste("Accuracy:", getAccuracy(cMatrix)))

library("cluster")
print("PAM clustering:")
plants_mean$PAM5 <- clusterPAM(plants_mean[,2:19], NULL)
cMatrix <- createSortedConfusionMatrix(plants_mean, "PAM5")
print("Confusion matrix for PAM clustering")
print(cMatrix)
print(paste("Accuracy:", getAccuracy(cMatrix)))

##### Question 2b #####
# Hierarchical clustering tuning
parametersHC <- tuneHC(plants_mean)
print(paste("The best parameter for HC is", parametersHC[1], "and the best parameter for dist is", parametersHC[2]))
plants_mean$HC5 <- clusterHC(plants_mean[,2:19], parametersHC)
cMatrix <- createSortedConfusionMatrix(plants_mean, "HC5")
print("Confusion matrix for tuned HC clustering")
print(cMatrix)
print(paste("Accuracy:", getAccuracy(cMatrix)))

# K-means clustering tuning
parametersKM <- tuneKM(plants_mean)
print(paste("The best parameter for nstart is", parametersKM[1], "and the best parameter for algorithm is", parametersKM[2]))
plants_mean$KM5 <- clusterKM(plants_mean[,2:19], parametersKM)
cMatrix <- createSortedConfusionMatrix(plants_mean, "KM5")
print("Confusion matrix for tuned KM clustering")
print(cMatrix)
print(paste("Accuracy:", getAccuracy(cMatrix)))

# PAM clustering tuning
parametersPAM <- tunePAM(plants_mean)
print(paste("The best parameters are; metric", parametersPAM[1], "stand", parametersPAM[2], "swap", parametersPAM[3]))
plants_mean$PAM5 <- clusterPAM(plants_mean[,2:19], parametersPAM)
cMatrix <- createSortedConfusionMatrix(plants_mean, "PAM5")
print("Confusion matrix for tuned PAM clustering")
print(cMatrix)
print(paste("Accuracy:", getAccuracy(cMatrix)))

##### Question 2c #####
print("K-means clustering with PCA dataset:")
plants_pca_10 <- as.data.frame(predict(pca,plants_stand)[,1:10])
plants_pca_10$Class <- plants_cleaned$Class
plants_pca_10$KM5 <- clusterKM(plants_pca_10[,1:10], parametersKM)
printMatrixAndAccuracy(plants_pca_10, "KM5")

print("K-means clustering with dataset after deletion:")
plants_2ci <- na.omit(plants_fi)
plants_2ci$KM5 <- clusterKM(plants_2ci[,1:17], parametersKM)
printMatrixAndAccuracy(plants_2ci, "KM5")

print("K-means clustering with dataset filled with 0:")
plants_0$KM5 <- clusterKM(plants_0[,2:19], parametersKM)
printMatrixAndAccuracy(plants_0, "KM5")

print("K-means clustering with dataset filled with mean:")
plants_mean$KM5 <- clusterKM(plants_mean[,2:19], parametersKM)
printMatrixAndAccuracy(plants_mean, "KM5")

print("K-means clustering with dataset filled with median:")
plants_med$KM5 <- clusterKM(plants_med[,2:19], parametersKM)
printMatrixAndAccuracy(plants_med, "KM5")
