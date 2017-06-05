required_packages <- c('outliers', 'DMwR', 'caret', 'MASS', 'entropy', 'shiny', 'DT', 'e1071')

for(i in 1:length(required_packages)){
  if(required_packages[i] %in% installed.packages()){
    library(required_packages[i], character.only = TRUE)
  }
  else{
    install.packages(required_packages[i], dependencies = TRUE)
    library(required_packages[i], character.only = TRUE)
  }
}


source('outlierDetectionMethods.R')

#datasetpath<-paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/waveform/Random/")
#datasetname<-paste("waveform-random-size500-sim100.csv")
#maxOut=10
#n=10
detectOutliers <- function(datasetpath,datasetname, maxOut, k, n){
 # maxOut=10
 # n=10
  #k=5
 
  #datasetPath <- sprintf('data/sources/%s.data', dataset)
  #dat <- read.csv(datasetPath, header = FALSE)
  #Take non categorical columns
  print('Detecting Outliers')
  dat <- read.csv(paste(datasetpath,datasetname, sep=""), header = FALSE) #original detection
 #   dat<-datasetpath #for secondary detection from samples
 dat <- sapply(dat,function(x) as.numeric(levels(x))[x])
 sizef<-dim(dat)[2]-7

 data_frame <- dat[-1,1:sizef]
 maxOut<-maxOut*dim(data_frame)[1]/100
  
 tempdf <- data.frame(matrix(0, nrow = nrow(data_frame), ncol = 0))

  for(i in 1:dim(data_frame)[2]){
    data_frame[,i] <- as.numeric(data_frame[,i])
    if(sum(data_frame[,i] != 0)){
      tempdf[,colnames(data_frame)[i]] <- data_frame[,i]
    }
  }
  
  data_frame <- tempdf
  #Create an empty list to record number of outlier for every method
  num_outliers <- list()
  
  #Creating a data frame to store output of every outlier detection method
  df_outlier <- data.frame(matrix(0, nrow = nrow(data_frame), ncol = 7))
  colnames(df_outlier) <- c('LOF','Mahalanobis','kMeans', 'ChiSq', 'BoxPlot', 'MAD', 'threeSigma')
 
  #' For a given dataset, implement outlier detection methods and fill values in the df_outlier matrix
  #' Try/Catch block used to check if NO outliers detected or if there's an ERROR in anyone of the methods
  #' Check the terminal for output of the error! 
  
 #LOF
 LOF <- NA
 out <- tryCatch({
   LOF <- lof(data_frame, maxOut, k)
   num_outliers$LOF <- maxOut
 },
 error = function(x){
   num_outliers$LOF <- NA
 },
 finally = {
   if(!is.na(LOF[1])){
     df_outlier$LOF <- replace(df_outlier$LOF, LOF, 1)  
     df_outlier$LOF <- replace(df_outlier$LOF, !LOF, 0)
   }
   else{
     df_outlier$LOF <- NA
   }
 }
 )
 
 #Mahalanobis
 Mahalanobis <- NA
 out <- tryCatch({
   Mahalanobis <- mahal(data_frame, maxOut)
   num_outliers$Mahalanobis <- maxOut
 },
 error = function(x){
   print(x)
   num_outliers$Mahalanobis <- NA},
finally = {
  if(length(Mahalanobis) > 0){
    if(!is.na(Mahalanobis[1])){
      df_outlier$Mahalanobis <- replace(df_outlier$Mahalanobis, Mahalanobis, 1)  
      df_outlier$Mahalanobis <- replace(df_outlier$Mahalanobis, !Mahalanobis, 0)
    }
    else{
      df_outlier$Mahalanobis <- NA
    } 
  }
}
 )

#Using KMeans
kMeans <- NA
out <- tryCatch({
  kMeans <- km(data_frame, maxOut, n)
  num_outliers$kMeans <- maxOut
},
error = function(x){
  print(x)
  num_outliers$kMeans <- NA
},
finally = {
  if(!is.na(kMeans[1])){
    df_outlier$kMeans <- replace(df_outlier$kMeans, kMeans, 1)  
    df_outlier$kMeans <- replace(df_outlier$kMeans, !kMeans, 0)  
  }
  else{
    df_outlier$kMeans <- NA
  }
}
)

#ChiSq
ChiSq <- NA
out <- tryCatch({
  ChiSq <- chisq(data_frame)
},
error = function(x){
  print(x)
  num_outliers$ChiSq <- NA
},
finally = {
  if(!is.na(ChiSq[1])){
    num_outliers$ChiSq <- length(ChiSq)
    df_outlier$ChiSq <- replace(df_outlier$ChiSq, ChiSq, 1) 
    df_outlier$ChiSq <- replace(df_outlier$ChiSq, !ChiSq, 0)
  }
  else{
    df_outlier$ChiSq <- NA
  }
}
)

#Boxplot
BoxPlot <- NA
out <- tryCatch({
  BoxPlot <- box(data_frame)
},
error = function(x){
  print(x)
  num_outliers$BoxPlot <- NA
},
finally = {
  if(length(BoxPlot)){
    if(!is.na(BoxPlot[1])){
      num_outliers$BoxPlot <- length(BoxPlot)
      df_outlier$BoxPlot <- replace(df_outlier$BoxPlot, BoxPlot, 1)
      df_outlier$BoxPlot <- replace(df_outlier$BoxPlot, !BoxPlot, 0)
    }
    else{
      df_outlier$BoxPlot <- NA
    } 
  }
}
)

#MAD
MAD <- NA
out <- tryCatch({
  MAD <- MADO(data_frame)
},
error = function(x){
  print(x)
  num_outliers$MAD <- NA
},
finally = {
  if(length(MAD) > 0){
    if(!is.na(MAD[1])){
      num_outliers$MAD <- length(MAD)
      df_outlier$MAD <- replace(df_outlier$MAD, MAD, 1)
      df_outlier$MAD <- replace(df_outlier$MAD, !MAD, 0)
    }
    else{
      df_outlier$MAD <- NA
    } 
  }
}
)

#ThreeSigma
threeSigma <- NA
out <- tryCatch({
  threeSigma <- threeSig(data_frame)
},
error = function(x){
  print(x)
  num_outliers$threeSigma <- NA
},
finally = {
  if(length(threeSigma) > 0){
    if(!is.na(threeSigma[1])){
      num_outliers$threeSigma <- length(threeSigma)
      df_outlier$threeSigma <- replace(df_outlier$threeSigma, threeSigma, 1)
      df_outlier$threeSigma <- replace(df_outlier$threeSigma, !threeSigma, 0)
    }
    else{
      df_outlier$threeSigma <- NA
    } 
  }
}
)
  #If all of them agree it's an outlier
  #for(i in 1:nrow(df_outlier)){
  #if(df_outlier[i,1] == df_outlier[i,2] && df_outlier[i,2] == df_outlier[i,3] && df_outlier[i,3] == df_outlier[i,4]
  #   && df_outlier[i,5] == df_outlier[i, 4] && df_outlier[i,5] == df_outlier[i, 6] && df_outlier[i,1] == 1){
  #    training_mat[i,5] = 1
  #  }
  #}
  
  #Calculate Disagreement for every column
  #disagreement <- vote_entropy(df_outlier)
  #Order disagreement and fetch only top values according to budged
  #disagreement <- order(disagreement, decreasing = TRUE)[1:maxOut]
 colnames(df_outlier) <- c('sLOF','sMAHA','skMEAN', 'sCHISQ', 'sBoxPlot', 'sMAD', 's3SIGMA')
 dataset <- sub("^([^.]*).*", "\\1", datasetname) 
  print(dataset)
  df_outlier <- cbind(dat[-1,], df_outlier)
  #filename = sprintf("%s_outliers.rds", dataset)
  #saveRDS(df_outlier, filename)
  filename = sprintf("%s_outliers.csv", dataset)
  write.csv(df_outlier, filename, row.names = FALSE)
}

#BLOCK WAVEFORM

Blk.num <- c(10, 25, 50, 10, 25, 50)
Blk.sizes <- c(25, 10, 5, 50, 20, 10)
Nsim <- 100
sample.sizes <- c(250, 500)
foldernames<-c("waveform")#,"ionosphere","australian","arrhythmia","iris",
             "pima-indians-diabetes","spambase","covtype","isolet5","Skin_NonSkin")
for (i in 1:length(foldernames)){
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/",foldernames[i],"/Block/", sep="")

for (i in 1:length(Blk.num)){
  for (n in 1:Nsim){
    if (n <10){
      filename <- paste("waveform-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("waveform-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim0", n, ".csv", sep="")
      } else  filename <- paste("waveform-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim", n, ".csv", sep="")
    }  
detectOutliers(pth, filename, 10,5,10) 
 }
}
}


#RANDOM WAVEFORM
foldernames<-c("waveform")#,"ionosphere","australian","arrhythmia","iris",
"pima-indians-diabetes","spambase","covtype","isolet5","Skin_NonSkin")
for (i in 1:length(foldernames)){
  pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/",foldernames[i],"/Random/", sep="")
for (m in 1:length(sample.sizes)){
  for (n in 1:Nsim){  
    if (n <10){
      filename <- paste("waveform-random-size", sample.sizes[m], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("waveform-random-size", sample.sizes[m], "-sim0", n, ".csv", sep="")
      } else filename <- paste("waveform-random-size", sample.sizes[m], "-sim", n, ".csv", sep="")
    }
    
    detectOutliers(pth, filename, 10,5,10) 
  }
}
}


#Block SPAM
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/spam/Block/")
Blk.num <- c(10, 25, 50, 10, 25, 50)
Blk.sizes <- c(20, 8, 4, 40, 16, 8)
Nsim <- 100
sample.sizes <- c(200, 400)


for (i in 1:length(Blk.num)){
  for (n in 1:Nsim){
    
    if (n <10){
      filename <- paste("spam-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("spam-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim0", n, ".csv", sep="")
      } else  filename <- paste("spam-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim", n, ".csv", sep="")
    }
    
    detectOutliers(pth, filename, 10,5,10) 
  }
}


#Random SPAM
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/spam/Random/")
for (m in 1:length(sample.sizes)){
  for (n in 1:Nsim){
    sample.data <- random.sample(spam, sample.sizes[m])
    if (n <10){
      filename <- paste("spam-random-size", sample.sizes[m], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("spam-random-size", sample.sizes[m], "-sim0", n, ".csv", sep="")
      } else filename <- paste("spam-random-size", sample.sizes[m], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10) 
  }
}


#Block Diabetes



Blk.num <- c(10, 25, 50, 10, 25, 50)
Blk.sizes <- c(5, 2, 1, 10, 4, 2)
Nsim <- 100
sample.sizes <- c(50, 100)
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/diabetes/Block/")
for (i in 1:length(Blk.num)){
  for (n in 1:Nsim){
    block.data <- block.sample(diabetes, Blk.num[i], Blk.sizes[i])
    if (n <10){
      filename <- paste("diabetes-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("diabetes-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim0", n, ".csv", sep="")
      } else  filename <- paste("diabetes-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10) 
  }
}

#Random Diabetes
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/diabetes/Random/")
for (m in 1:length(sample.sizes)){
  for (n in 1:Nsim){
    sample.data <- random.sample(diabetes, sample.sizes[m])
    if (n <10){
      filename <- paste("diabetes-random-size", sample.sizes[m], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("diabetes-random-size", sample.sizes[m], "-sim0", n, ".csv", sep="")
      } else filename <- paste("diabetes-random-size", sample.sizes[m], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10) 
  }
  
}

#Block Iris PROBLEM KMean+3sigma
Blk.num <- c(1, 5, 10, 1, 5, 10)
Blk.sizes <- c(10, 2, 1, 20, 4, 2)
Nsim <- 100
sample.sizes <- c(10, 20)
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/iris/Block/")


for (i in 1:length(Blk.num)){
  for (n in 1:Nsim){
   
    if (n <10){
      filename <- paste("iris-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("iris-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim0", n, ".csv", sep="")
      } else  filename <- paste("iris-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,3) 
  }
}

#Random Iris
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/iris/Random/")
for (m in 1:length(sample.sizes)){
  for (n in 1:Nsim){
   
    if (n <10){
      filename <- paste("iris-random-size", sample.sizes[m], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("iris-random-size", sample.sizes[m], "-sim0", n, ".csv", sep="")
      } else filename <- paste("iris-random-size", sample.sizes[m], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,3) 
  }
  
}


#Block ionosphere
Blk.num <- c(1, 5, 10, 1, 5, 10)
Blk.sizes <- c(20, 4, 2, 40, 8, 4)
Nsim <- 100
sample.sizes <- c(20, 40)

pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/ionosphere/Block/")


for (i in 1:length(Blk.num)){
  for (n in 1:Nsim){
    block.data <- block.sample(ionosphere, Blk.num[i], Blk.sizes[i])
    if (n <10){
      filename <- paste("ionosphere-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("ionosphere-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim0", n, ".csv", sep="")
      } else  filename <- paste("ionosphere-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10) 
  }
}
#Random ionosphere
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/ionosphere/Random/")
for (m in 1:length(sample.sizes)){
  for (n in 1:Nsim){
   
    if (n <10){
      filename <- paste("ionosphere-random-size", sample.sizes[m], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("ionosphere-random-size", sample.sizes[m], "-sim0", n, ".csv", sep="")
      } else filename <- paste("ionosphere-random-size", sample.sizes[m], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10) 
  }
  
}


#Block australian
Blk.num <- c(1, 5, 10, 1, 5, 10)
Blk.sizes <- c(30, 6, 3, 60, 12, 6)
Nsim <- 100
sample.sizes <- c(30, 60)
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/australian/Block/")

for (i in 1:length(Blk.num)){
  for (n in 1:Nsim){
 
    if (n <10){
      filename <- paste("australian-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("australian-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim0", n, ".csv", sep="")
      } else  filename <- paste("australian-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10)  }
}

#Random australian
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/australian/Random/")

for (m in 1:length(sample.sizes)){
  for (n in 1:Nsim){

    if (n <10){
      filename <- paste("australian-random-size", sample.sizes[m], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("australian-random-size", sample.sizes[m], "-sim0", n, ".csv", sep="")
      } else filename <- paste("australian-random-size", sample.sizes[m], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10)  }
  
}




#Random arrythmia
Blk.num <- c(1, 5, 10, 1, 5, 10)
Blk.sizes <- c(20, 4, 2, 50, 10, 5)
Nsim <- 100
sample.sizes <- c(20, 50)
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/arrhythmia/Block/")
for (i in 1:length(Blk.num)){
  for (n in 1:Nsim){
    block.data <- block.sample(arrhythmia, Blk.num[i], Blk.sizes[i])
    if (n <10){
      filename <- paste("arrhythmia-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim00", n, ".csv", sep="")
    } else{
      if (n <100){
        filename <- paste("arrhythmia-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim0", n, ".csv", sep="")
      } else  filename <- paste("arrhythmia-block-num", Blk.num[i], "-size", Blk.sizes[i], "-sim", n, ".csv", sep="")
    }
    detectOutliers(pth, filename, 10,5,10)
  }
}
#Random arrythmia
pth <- paste("/Users/Laure/Desktop/home/Aberti/Experiments/Datasets/JM-Expts/RWdata/arrhythmia/Random/")

for (m in 1:length(sample.sizes)){
  for (n in 1:9){
    for (n in 1:Nsim){
    print(n)
 #  if (n <10){
      filename <- paste("arrhythmia-random-size", sample.sizes[m], "-sim00", n, ".csv", sep="")
  #  } 
   #else{
   #   if (n <100){
  #      filename <- paste("arrhythmia-random-size", sample.sizes[m], "-sim0", n, ".csv", sep="")
   #   } else filename <- paste("arrhythmia-random-size", sample.sizes[m], "-sim", n, ".csv", sep="")

      detectOutliers(pth, filename, 10,5,10)
    }
  }  
}

