################################################################
### Developped by Ji Meng Loh 
### Functions to generate synthetic data         
### with controlled distributions and percentage of outliers 
################################################################



F.generateData <- function(N, pois.mu, norm.mu, norm.sd, norm2.mu, norm2.cov, str.length){
  #### to generate a poisson rv, a univariate normal, a trivariate normal and a random string
  #### return as a data.frame
#  tri.norm <- mvrnorm(N, mu=norm3.mu, Sigma=norm3.cov)
  bi.norm <- rmvnorm(N, mean=norm2.mu, sigma=norm2.cov)
  rand.norm <- rnorm(N, norm.mu, norm.sd)
  rand.pois <- rpois(N, pois.mu)
  
  rand.str <- c(1:str.length)  ## initialize
  for (i in 1:Nrow){
    rand.str[i] <- paste(sample(c(0:9, letters, LETTERS), str.length, replace=TRUE), collapse="")
  }
  mydata <- data.frame(str=rand.str, pois=rand.pois, norm=rand.norm, binorm1=bi.norm[,1], binorm2=bi.norm[,2])
  return(mydata)
}


### functions to add glitches (addMissing and addDuplicate functions are not used in the current paper)

F.addMissing <- function(good.data, prop.missing=0.1, which.cols=NULL){
  ### add missing to dataset, specifying how much missing to add and to which cols
  ### which.cols should NOT have repeated numbers
  bad.data <- good.data
  if (is.null(which.cols)){  ## i.e. any col
    total <- nrow(good.data) * ncol(good.data)
    num.missing <- max(round(prop.missing * total), 1)
    where.missing <- sample(c(1:total), num.missing, replace=FALSE)
    
    row.locations <- (where.missing+(ncol(bad.data)-1)) %/% ncol(bad.data)
    col.locations <- where.missing %% ncol(bad.data)
    col.locations[col.locations==0] <- ncol(bad.data)
#    cat(row.locations, "\n\n")
#    cat(col.locations, "\n\n")
    for (i in 1:length(row.locations)){
      bad.data[row.locations[i], col.locations[i]] <- NA
    }
  } else{ ## only in columns specified in which.cols
    ncols <- length(which.cols)
    total <- nrow(good.data) * ncols
    num.missing <- max(round(prop.missing * total), 1)
    where.missing <- sample(c(1:total), num.missing, replace=FALSE)
    row.locations <- (where.missing+(ncols-1)) %/% ncols
    col.locations <- where.missing %% ncols
    col.locations[col.locations==0] <- ncols
#    cat(row.locations, "\n\n")
#    cat(col.locations, "\n\n")
    for (i in 1:length(row.locations)){
      bad.data[row.locations[i], which.cols[col.locations[i]]] <- NA
    }
  }
  return(bad.data)
}


F.addDuplicates <- function(good.data, duplicate.dist, duplicate.prop){
  ### NOTE: THESE ARE EXACT DUPLICATES - NEED TO ADD SMALL CHANGES TO THE DUPLICATED RECORDS SO THAT THEY ARE NOT EXACT DUPLICATES
  
  ### duplicate records
  ### duplicate.dist lists the probability distribution that a record is duplicated 1, 2, ... times
  num.dup <- max(round(duplicate.prop * nrow(good.data)), 1)
  records.dup <- sample(1:nrow(good.data), num.dup, replace=FALSE)
  dup.times <- sample(1:length(duplicate.dist), num.dup, replace=TRUE, prob=duplicate.dist)
  bad.data <- good.data
  
  for (i in 1:num.dup){  ## for each record to be duplicated, which is in row records.dup[i]
    for (j in 1:dup.times[i]){ ### take the number of times it is duplicated
      bad.data <- rbind(bad.data, bad.data[records.dup[i],])  ### duplicate it
    }
  }
  return(bad.data)
}

F.addOutliers.norm <- function(good.data, out.mean, out.sd, out.prop, clustered=FALSE){
  out.num <- max(round(out.prop * nrow(good.data)), 1)  ### number of outliers
  out.values <- rnorm(out.num, out.mean, out.sd)
  if (clustered){ ### add code to generate clustered locations
  } else out.where <- sample(1:nrow(good.data), out.num, replace=FALSE)  ### where they will occur
  
  bad.data <- good.data
  bad.data$norm[out.where] <- out.values
  return(bad.data)
}

F.addOutliers.pois <- function(good.data, out.mean, out.prop, clustered=FALSE){
  out.num <- max(round(out.prop * nrow(good.data)), 1)  ### number of outliers
  out.values <- rpois(out.num, out.mean)
  if (clustered){ ### add code to generate clustered locations
  } else out.where <- sample(1:nrow(good.data), out.num, replace=FALSE)  ### where they will occur
  
  bad.data <- good.data
  bad.data$pois[out.where] <- out.values
  return(bad.data)
}

F.addOutliers.trinorm <- function(good.data, out.mean, out.sd, out.prop, clustered=FALSE, Ncols=3,
                                  col.names=c("trinorm1", "trinorm2", "trinorm3")){
  if (length(col.names) != Ncols){
    cat("Number of col.names is not equal to", Ncols, "\n")
    return()
  }
  out.num <- max(round(out.prop * nrow(good.data)*Ncols), 1)  ### number of outliers
  out.values <- rnorm(out.num, out.mean, out.sd)
  bad.data <- good.data
  if (clustered){ ### add code to generate clustered locations
  } else{
    
    out.where <- sample(1:(Ncols*nrow(good.data)), out.num, replace=FALSE)  ### where they will occur
    row.locations <- (out.where+(Ncols-1)) %/% Ncols
    col.locations <- out.where %% Ncols
    col.locations[col.locations==0] <- Ncols
    tri.col <- list()
    for (i in 1:Ncols){
      tri.col[[i]] <- row.locations[col.locations==i]
    }
    start <- 1
    for (i in 1:length(tri.col)){
      bad.data[[col.names[i]]][tri.col[[i]]] <- out.values[start:(start+length(tri.col[[i]])-1)]
      start <- start + length(tri.col[[i]])            
    }        
  }
  return(bad.data)
}

F.addOutliers.binorm <- function(good.data, out.mean, out.cov, out.prop, clustered=FALSE, Ncols=2,
                                  col.names=c("binorm1", "binorm2")){
  if (length(col.names) != Ncols){
    cat("Number of col.names is not equal to", Ncols, "\n")
    return()
  }
  out.num <- max(round(out.prop * nrow(good.data)), 1)  ### number of outliers
#  out.values <- rnorm(out.num, out.mean, out.sd)
  out.values <- rmvnorm(out.num, mean=out.mean, sigma=out.cov)
  bad.data <- good.data
  if (clustered){ ### add code to generate clustered locations
  } else{
    
    out.where <- sample(1:nrow(good.data), out.num, replace=FALSE)  ### where they will occur
    for (i in 1:out.num){
      bad.data[[col.names[1]]][out.where[i]] <- out.values[i, 1]
      bad.data[[col.names[2]]][out.where[i]] <- out.values[i, 2]
    }        
  }
  return(bad.data)
}

