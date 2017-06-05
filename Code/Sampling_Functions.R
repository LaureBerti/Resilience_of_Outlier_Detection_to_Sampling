
#### R functions developped by Ji Meng Loh to do sampling and computing of glitch summaries

# Stratified sampling function (Parni's code)
stratified.sample <- function(y,k) {
	#y is matrix from which sampling
	# M=desired size of bootstrap sample, N=size of original sample/data
	# k column on which to stratify for sampling
	sp <- split(y, y[,k]) # splitting on combination of variables list(d$age, d$lc))
	samples <- lapply(sp, function(x,M,N) x[sample(1:nrow(x), M*nrow(x)/N, T),],M=nrow(y),N=nrow(y))
	boot <- do.call(rbind, samples)
}

### Block sampling
block.sample <- function(y, Blk.num, Blk.size){
    ## y is the data
    ## Blk.num is the number of blocks
    ## Blk.size is the size of the blocks (measured by number of consecutive rows)
    start.pts <- c(1:(nrow(y)-Blk.size+1))
    blk.start <- sample(start.pts, Blk.num, replace=T)
    samples <- lapply(blk.start, function(x, n) y[x:(x+n-1), ], n=Blk.size)
    boot <- do.call(rbind, samples)
    return(boot)
}

### simple random sampling
random.sample <- function(y, sample.size){
    ## y is the data
    ## sample.size is the size of the bootstrap sample
    samples <- sample(c(1:nrow(y)), sample.size, replace=T)
    boot <- y[samples,]
    return(boot)
}

#### subsetting
subset.sample <- function(y, numsets, myname){
    ### break y up into subsets based on the proportion
    y.size <- nrow(y)
    size <- round(y.size / numsets)
    start <- 1
    for (N in 1:numsets){
        end <- start + size - 1
        if (N != numsets){
          y.subset <- y[(start:end), ]
        } else y.subset <- y[start:y.size,]
        start <- start+size
        if (N < 10){
           filename <- paste(myname, "-subset", numsets, "-00", N, ".txt", sep="")
        } else {
             if (N < 100){
                filename <- paste(myname, "-subset", numsets, "-0", N, ".txt", sep="")
             } else filename <- paste(myname, "-subset", numsets, "-", N, ".txt", sep="")
        }
        write.table(y.subset, filename, sep="|", row.names=F, col.names=T, append=F, quote=F)
    }
}


#### Next functions to compute the glitch summary statistics of a given dataset (actual or bootstrapped)
### 1. totals or proportions of each type of glitch
### 2. glitch correlations - number of glitches that are 1, 2, 3, etc rows from any glitch

### y is the data, plus glitch matrix in columns c1 to c2
### d is the range that we will compute correlations to
F.Glitch <- function(y, c1, c2, d){
    glitch.mat <- y[, c1:c2]
    glitch.prop <- apply(glitch.mat, 2, mean) #### proportion of each glitch type in y
    glitch.cor <- sapply(glitch.mat, F.glitch.cor, d=d)  ### glitch.cor should be d x (c2-c1+1), i.e. each column is a pair correlation function up to distance d
    Glitch.stat <- list(proportion=glitch.prop, correlation=glitch.cor)
    return(Glitch.stat)
}


F.glitch.cor <- function(y, d){ ### computes the pair correlation function for the vector y of 1's and 0's
    ### to distance separation d
    y[is.na(y)] <- 0
    n <- length(y)
    ng <- sum(y) # number of glitches
    ans <- rep(0, d)
    glitch.locations <- which(y==1)
    for (i in glitch.locations){
        v1 <- (i + c(1:d))
        v2 <- (i - c(1:d))
        num.obs <- (v1 <= n) + (v2 > 0) ## vector of TRUE and FALSE as to whether points in v1 or v2 are outside obs region
        wt <- 1/num.obs ### edge correction
        ans <- ans + wt * (c(y[v1[v1<=n]], rep(0, sum(v1>n))) + c(rep(0, sum(v2<1)), y[v2[v2>0]]))
    }
    ans <- n * ans / (ng^2)
    if (ng==0) ans <- rep(0,d)
    return(ans)
}


F.processGlitchSummaries <- function(blockdata, randomdata, subsetdata){
    ### uses Blk.num, Nsim, cor.length, sample.sizes, and mysubsets variables in global workspace
    ans <- list()
    ans$Block <- list()
    for (m in 1:length(Blk.num)){
        temp  <- list()
        temp1 <- NULL
        temp2 <- NULL
        for (i in 1:Nsim){
            temp1 <- rbind(temp1, blockdata[[m]][[i]]$proportion)
            temp2 <- rbind(temp2, blockdata[[m]][[i]]$correlation)
        }
        temp$prop.mean.sd <- rbind(apply(temp1, 2, mean), apply(temp1, 2, sd))

        cor.mean <- NULL
        cor.sd <- NULL
        for (j in 1:ncol(temp2)){
            temp3 <- matrix(temp2[,j], nrow=cor.length, byrow=F)
            cor.mean <- cbind(cor.mean, apply(temp3, 1, mean))
            cor.sd <- cbind(cor.sd, apply(temp3, 1, sd))
        }
        temp$cor.mean <- cor.mean
        temp$cor.sd <- cor.sd
        ans$Block[[m]] <- temp
    }
    
    ans$Random <- list()
    for (m in 1:length(sample.sizes)){
        temp  <- list()
        temp1 <- NULL
        temp2 <- NULL
        for (i in 1:Nsim){
            temp1 <- rbind(temp1, randomdata[[m]][[i]]$proportion)
            temp2 <- rbind(temp2, randomdata[[m]][[i]]$correlation)
        }
        temp$prop.mean.sd <- rbind(apply(temp1, 2, mean), apply(temp1, 2, sd))

        cor.mean <- NULL
        cor.sd <- NULL
        for (j in 1:ncol(temp2)){
            temp3 <- matrix(temp2[,j], nrow=cor.length, byrow=F)
            cor.mean <- cbind(cor.mean, apply(temp3, 1, mean))
            cor.sd <- cbind(cor.sd, apply(temp3, 1, sd))
        }
        temp$cor.mean <- cor.mean
        temp$cor.sd <- cor.sd
        ans$Random[[m]] <- temp
    }
    
    ans$Subset <- list()
    for (m in 1:length(mysubsets)){
        temp  <- list()
        temp1 <- NULL
        temp2 <- NULL
        for (i in 1:mysubsets[m]){
            temp1 <- rbind(temp1, subsetdata[[m]][[i]]$proportion)
            temp2 <- rbind(temp2, subsetdata[[m]][[i]]$correlation)
        }
        temp$prop.mean.sd <- rbind(apply(temp1, 2, mean), apply(temp1, 2, sd))
        
        cor.mean <- NULL
        cor.sd <- NULL
        for (j in 1:ncol(temp2)){
            temp3 <- matrix(temp2[,j], nrow=cor.length, byrow=F)
            cor.mean <- cbind(cor.mean, apply(temp3, 1, mean))
            cor.sd <- cbind(cor.sd, apply(temp3, 1, sd))
        }
        temp$cor.mean <- cor.mean
        temp$cor.sd <- cor.sd
        ans$Subset[[m]] <- temp
    }
    return(ans)
}


F.processGlitchSummaries1 <- function(blockdata, randomdata, subsetdata){
    ### uses Blk.num, Nsim, cor.length, sample.sizes, and mysubsets variables in global workspace
    ans <- list()
    ans$Block <- list()
    for (m in 1:length(Blk.num)){
        temp  <- list()
        temp1 <- NULL
        temp2 <- NULL
        for (i in 1:Nsim){
            temp1 <- rbind(temp1, blockdata[[m]][[i]]$proportion)
            temp2 <- rbind(temp2, blockdata[[m]][[i]]$correlation)
        }
        temp$prop.mean.sd <- rbind(apply(temp1, 2, mean), apply(temp1, 2, sd))

        cor.mean <- NULL
        cor.sd <- NULL
        for (j in 1:ncol(temp2)){
            temp3 <- matrix(temp2[,j], nrow=cor.length, byrow=F)
            cor.mean <- cbind(cor.mean, apply(temp3, 1, mean))
            cor.sd <- cbind(cor.sd, apply(temp3, 1, sd))
        }
        temp$cor.mean <- cor.mean
        temp$cor.sd <- cor.sd
        ans$Block[[m]] <- temp
    }
    
    ans$Random <- list()
    for (m in 1:length(sample.sizes)){
        temp  <- list()
        temp1 <- NULL
        temp2 <- NULL
        for (i in 1:Nsim){
            temp1 <- rbind(temp1, randomdata[[m]][[i]]$proportion)
            temp2 <- rbind(temp2, randomdata[[m]][[i]]$correlation)
        }
        temp$prop.mean.sd <- rbind(apply(temp1, 2, mean), apply(temp1, 2, sd))

        cor.mean <- NULL
        cor.sd <- NULL
        for (j in 1:ncol(temp2)){
            temp3 <- matrix(temp2[,j], nrow=cor.length, byrow=F)
            cor.mean <- cbind(cor.mean, apply(temp3, 1, mean))
            cor.sd <- cbind(cor.sd, apply(temp3, 1, sd))
        }
        temp$cor.mean <- cor.mean
        temp$cor.sd <- cor.sd
        ans$Random[[m]] <- temp
    }
    
    ans$Subset <- list()
    for (m in 1:length(mysubsets)){
        temp  <- list()
        temp1 <- NULL
        temp2 <- NULL
        for (i in 1:mysubsets[m]){
            if (i==20) break
            temp1 <- rbind(temp1, subsetdata[[m]][[i]]$proportion)
            temp2 <- rbind(temp2, subsetdata[[m]][[i]]$correlation)
        }
        temp$prop.mean.sd <- rbind(apply(temp1, 2, mean), apply(temp1, 2, sd))
        
        cor.mean <- NULL
        cor.sd <- NULL
        for (j in 1:ncol(temp2)){
            temp3 <- matrix(temp2[,j], nrow=cor.length, byrow=F)
            cor.mean <- cbind(cor.mean, apply(temp3, 1, mean))
            cor.sd <- cbind(cor.sd, apply(temp3, 1, sd))
        }
        temp$cor.mean <- cor.mean
        temp$cor.sd <- cor.sd
        ans$Subset[[m]] <- temp
    }
    return(ans)
}
