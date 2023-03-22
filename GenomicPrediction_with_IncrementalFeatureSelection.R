#Loading required libraries####
library(ranger)
library(data.table)
library(ggplot2)

#set.seed(42) #Uncomment to reproduce the published results

#Parse args####
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  print("Usage: Rscript GenomicPrediction_with_IncrementalFeatureSelection plinkBinaryPrefix threadCount")
  q()
}
plinkPrefix = args[1]
threadCount = as.numeric(args[2])

#Additional parameters####
chrCount = 40 # Required for PLINK
#Cross validation
repeatNumber = 10
foldNumber = 5
#Incremental feature selection
stepLength = 5
stepLength2 = 10
stepLength3 = 100
stepLength4 = 500
stepLength5 = 1000
stepLength6 = 5000
#Random forest
numberTrees = 500

#Preparing and reading data####
#Create temporary raw file using PLINK
system(command = paste0("plink --bfile ", plinkPrefix, " --allow-no-sex --recode A --out tempRawFile --chr-set ", chrCount), show.output.on.console = F)
rawData = fread("tempRawFile.raw", header = T, data.table = F)
snpData = subset(rawData, select = -c(FID,IID,PAT,MAT,SEX))
#Remove allele suffix from SNP IDs
colnames(snpData)[2:ncol(snpData)] = sub("_[^_]+$", "", colnames(snpData)[2:ncol(snpData)])
snpCount = ncol(rawData) - 6
sampleCount = nrow(rawData)
sampleIndices = 1:sampleCount
foldsIndices <- cut(sampleIndices,breaks=foldNumber,labels=FALSE)

#Determine number of steps####
if(snpCount < 100){
  steps = c(1:snpCount)
} else{
  if(snpCount < 500){
    steps = c(c(1:100), seq(from = 101, to = snpCount, by = stepLength), snpCount)
  } else{
    if(snpCount < 1000){
      steps = c(c(1:100), seq(from = 101, to = 500, by = stepLength), seq(from = 501, to = snpCount, by = stepLength2), snpCount)
    } else{
      if(snpCount < 5000){
        steps = c(c(1:100), seq(from = 101, to = 500, by = stepLength), seq(from = 501, to = 1000, by = stepLength2),  seq(from = 1001, to = snpCount, by = stepLength3), snpCount)
      } else{
        if(snpCount < 10000){
          steps = c(c(1:100), seq(from = 101, to = 500, by = stepLength), seq(from = 501, to = 1000, by = stepLength2),  seq(from = 1001, to = 5000, by = stepLength3), seq(from = 5001, to = snpCount, by = stepLength4), snpCount)
        } else{
          if(snpCount < 50000){
            steps = c(c(1:100), seq(from = 101, to = 500, by = stepLength), seq(from = 501, to = 1000, by = stepLength2),  seq(from = 1001, to = 5000, by = stepLength3), seq(from = 5001, to = 10000, by = stepLength4), seq(from = 10001, to = snpCount, by = stepLength5), snpCount)
          } else{
            steps = c(c(1:100), seq(from = 101, to = 500, by = stepLength), seq(from = 501, to = 1000, by = stepLength2),  seq(from = 1001, to = 5000, by = stepLength3), seq(from = 5001, to = 10000, by = stepLength4), seq(from = 10001, to = 50000, by = stepLength5), seq(from = 50001, to = snpCount, by = stepLength6), snpCount)
          }
        }
      }
    }
  } 
}
# Create dataframes to store results####
allSNPsResultsDF = data.frame(matrix(nrow = 1, ncol = repeatNumber*2))
resultsDF = data.frame(matrix(nrow = length(steps), ncol=repeatNumber*2+1))
tmpCounter = 1
for(i in 1:repeatNumber){
  colnames(resultsDF)[tmpCounter] = paste0("Repetition_",i,"_R2")
  colnames(allSNPsResultsDF)[tmpCounter] = paste0("Repetition_",i,"_R2")
  tmpCounter = tmpCounter+1
  colnames(resultsDF)[tmpCounter] = paste0("Repetition_",i,"_SE")
  colnames(allSNPsResultsDF)[tmpCounter] = paste0("Repetition_",i,"_SE")
  tmpCounter = tmpCounter+1
}
colnames(resultsDF)[tmpCounter] = "SNP_Count"
resultsDF$SNP_Count = steps

#Running genomic prediction with cross-validation####
for(i in 1:repeatNumber){
  cat("Repetition: ", i,"\n")
  repeatIndices = sample(sampleIndices)
  rawData = rawData[repeatIndices,]
  snpData = snpData[repeatIndices,]
  resultAllSNPMatrix = matrix(nrow = 1, ncol = sampleCount, dimnames = list("AllSNPs", rawData$IID))
  resultMatrix = matrix(nrow = length(steps), ncol = sampleCount, dimnames = list(steps, rawData$IID))
  for(j in 1:foldNumber){
    cat("Fold: ", j,"\n")
    testIndices = repeatIndices[which(foldsIndices == j, arr.ind = T)]
    trainIndices = setdiff(repeatIndices, testIndices)
    #Writing sample IDs for training in file for PLINK
    write.table(rawData[trainIndices,c("FID","IID")], file = "tempTrainingSetSampleIDs.txt", quote = F, row.names = F, col.names = F)
    #Running PLINK GWAS
    system(command = paste0("plink --bfile ", plinkPrefix, " --allow-no-sex --assoc --out tempTrainingResult --keep tempTrainingSetSampleIDs.txt --chr-set ", chrCount), show.output.on.console = F)
    #Read in PLINK result
    plinkResult = fread("tempTrainingResult.qassoc", header = T)
    orderedSNPs = plinkResult[order(plinkResult$P, decreasing = F),]$SNP
    #Prepare training and test data with sorted SNPs
    trainData = snpData[trainIndices, c("PHENOTYPE",orderedSNPs)]
    testData = snpData[testIndices, c("PHENOTYPE",orderedSNPs)]
    #Run genomic prediction with all SNPs
    rf = ranger(x = subset(trainData, select =  -get("PHENOTYPE")), y = trainData[,1], num.threads = threadCount, num.trees = numberTrees)
    predValues = predict(rf, testData)$predictions
    resultAllSNPMatrix[1,testIndices] = predValues
    #Run genomic prediction with incrementally added SNPs
    entryCounter = 1
    for(k in steps){
      cat("Number of SNPs: ", k,"\r")
      currentTrainData = trainData[,1:(k+1)]
      currentTestData = testData[,1:(k+1)]
      rf = ranger(x = subset(currentTrainData, select =  -get("PHENOTYPE")), y = currentTrainData[,1], num.threads = threadCount, num.trees = numberTrees)
      predValues = predict(rf, currentTestData)$predictions
      resultMatrix[entryCounter,testIndices] = predValues
      entryCounter = entryCounter+1
    }
  }
  #Combine results from all folds to one vector and calculate rsquared and standard error####
  rSquaredVec = rep(-1,nrow(resultMatrix))
  seVec = rep(-1, nrow(resultMatrix))
  for(j in 1:nrow(resultMatrix)){
    predValues = resultMatrix[j,]
    rSquaredVec[j] = 1 - sum(((predValues - snpData$PHENOTYPE)^2)) / sum(((mean(snpData$PHENOTYPE) - snpData$PHENOTYPE)^2))
    foldRSquaredVec = rep(-1, foldNumber)
    for(k in 1:foldNumber){
      testIndices = repeatIndices[which(foldsIndices == k, arr.ind = T)]
      predValues = resultMatrix[j,testIndices]
      foldRSquaredVec[k] = 1 - sum(((predValues - snpData$PHENOTYPE[testIndices])^2)) / sum(((mean(snpData$PHENOTYPE[testIndices]) - snpData$PHENOTYPE[testIndices])^2))
    }
    seVec[j] = sd(foldRSquaredVec) / sqrt(foldNumber)
  }
  resultsDF[,paste0("Repetition_",i,"_R2")] = rSquaredVec
  resultsDF[,paste0("Repetition_",i,"_SE")] = seVec
  #Evaluating model with all SNPs####
  predValues = resultAllSNPMatrix[1,]
  rSquaredAllSNPs = 1 - sum(((predValues - snpData$PHENOTYPE)^2)) / sum(((mean(snpData$PHENOTYPE) - snpData$PHENOTYPE)^2))
  foldRSquaredVec = rep(-1, foldNumber)
  for(j in 1:foldNumber){
    testIndices = repeatIndices[which(foldsIndices == j, arr.ind = T)]
    predValues = resultAllSNPMatrix[1,testIndices]
    foldRSquaredVec[j] = 1 - sum(((predValues - snpData$PHENOTYPE[testIndices])^2)) / sum(((mean(snpData$PHENOTYPE[testIndices]) - snpData$PHENOTYPE[testIndices])^2))
  }
  seAllSNPs = sd(foldRSquaredVec) / sqrt(foldNumber)
  allSNPsResultsDF[1,paste0("Repetition_",i,"_R2")] = rSquaredAllSNPs
  allSNPsResultsDF[1,paste0("Repetition_",i,"_SE")] = seAllSNPs
}

#Combine results from all repetitions to calculate mean and standard error####
resultsDF$Mean_R2 = c()
resultsDF$Mean_SE = c()
for(i in 1:length(steps)){
  meanVec = c()
  for(j in 1:repeatNumber){
    meanVec = c(meanVec,resultsDF[i,paste0("Repetition_",j,"_R2")])
  }
  resultsDF$Mean_R2[i] = mean(meanVec)
  resultsDF$Mean_SE[i] = sd(meanVec) / sqrt(repeatNumber)
}
allSNPsResultsDF$SNP_Count = "AllSNPs"
allSNPsResultsDF$Mean_R2 = c()
allSNPsResultsDF$Mean_SE = c()
meanVec = c()
for(i in 1:repeatNumber){
  meanVec = c(meanVec,allSNPsResultsDF[1,paste0("Repetition_",i,"_R2")])
}
allSNPsResultsDF$Mean_R2[1] = mean(meanVec)
allSNPsResultsDF$Mean_SE[1] = sd(meanVec) / sqrt(repeatNumber)

allResults = rbind(resultsDF,allSNPsResultsDF)

write.table(allResults,paste0(plinkPrefix,"_ResultDF"))

#Create plot####
#Log transformation of SNP counts for better visualization
resultsDF$logSNP_Count = log10(resultsDF$SNP_Count)
#Apply Friedman's super smoother to determine the trend curve and its maximum
fit_supersmooth = supsmu(resultsDF$logSNP_Count, resultsDF$Mean_R2)
loessMaxValue = max(fit_supersmooth$y)
loessMaxIndex = fit_supersmooth$x[which(fit_supersmooth$y == loessMaxValue)]
if( round(loessMaxValue - allSNPsResultsDF$Mean_R2,3) <= 0 | round((10^(loessMaxIndex) / max(resultsDF$SNP_Count) * 100),2) == 100){
  cat(paste0("No increase of R2 possible\n"))
}else{
  cat(paste0("Increase of R2 by ", round(loessMaxValue-allSNPsResultsDF$Mean_R2,3), " percentage points (from ", round(allSNPsResultsDF$Mean_R2,3), " to ", round(loessMaxValue,3),") using ", round((10^(loessMaxIndex) / max(resultsDF$SNP_Count) * 100),2), "% of all SNPs (",round(10^(loessMaxIndex),0), " out of ", max(resultsDF$SNP_Count), " SNPs)\n"))
}
resultsDF$x = fit_supersmooth$x
resultsDF$y = fit_supersmooth$y
resultsDF$maxIdx = loessMaxIndex
p = ggplot(resultsDF, aes(x = logSNP_Count, y= Mean_R2)) + geom_point() + 
  geom_line(aes(x = x, y = y)) + geom_vline(aes(xintercept = maxIdx), colour = "green", size = 1) +
  geom_hline(aes(yintercept = allSNPsResultsDF$Mean_R2), colour = "black")  + geom_ribbon(aes(ymin = allSNPsResultsDF$Mean_R2 - allSNPsResultsDF$Mean_SE, ymax = allSNPsResultsDF$Mean_R2 + allSNPsResultsDF$Mean_SE), colour = "black", alpha = 0.3) +
  labs(title = paste0(plinkPrefix),y = expression("R"^2), x = expression("log"["10"]*"(#SNPs)")) +
  theme(text = element_text(size = 20))
ggsave(paste0(plinkPrefix,"_ResultPlot.png"), plot = p, device = "png", width = 2560/100, height = 1380/100, dpi = 100)
