# Original error estimation function
loessErrfun <- function (trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = tot)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-07
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4, 
  ], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ], 1 - 
    colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12, 
    ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), 
                          "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}


#enforce monotonicity by changing the parameters (span=2 and the log-transformed totals as weights) of the loess function
loessErrfun_mod1 <- function (trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A",
                                                "C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot,
                         rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-07
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4,
  ], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ], 1 -
    colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12,
    ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),
                          "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}

# alter loess arguments (weights and span) & enforce monotonicity
loessErrfun_mod2 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# Just enforce monotonicity
loessErrfun_mod3 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# Loess weights only
loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        # only change the weights
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# Fit erorr models

library(dada2)
library(tidyverse)
set.seed(666) # set random seed for reproducability
samdf <- read.csv("sample_data/Sample_info.csv", stringsAsFactors = FALSE)
runs <- unique(samdf$fcid)

dir.create("output/error_models", recursive=TRUE)

i=1
for (i in 1:length(runs)){
  
  run_data <- samdf %>%
    filter(fcid == runs[i])
  
  filtpath <- paste0("data/", runs[i], "/02_filtered" )
  
  # Load forward reads
  filtFs <- list.files(filtpath, pattern="R1_001.*", full.names = TRUE)
  filtFs <- filtFs[!str_detect(filtFs, "Undetermined")]
  
  # Load reverse reads
  filtRs <- list.files(filtpath, pattern="R2_001.*", full.names = TRUE)
  filtRs <- filtRs[!str_detect(filtRs, "Undetermined")]
  
  # Loop over the error est functions
  err_funs <- grep("loessErrfun",names(.GlobalEnv),value=TRUE)
  
  for (e in 1:length(err_funs)){
    
    # Learn error rates from a subset of the samples and reads (rather than running self-consist with full dataset)
    errF <- learnErrors(filtFs,
                        multithread = 8, 
                        nbases = 1e+9, # set to 1e+9
                        randomize = TRUE, 
                        qualityType = "FastqQuality", 
                        errorEstimationFunction = get(err_funs[e]),
                        verbose=TRUE)
    
    errR <- learnErrors(filtFs,
                        multithread = 8, 
                        nbases = 1e+9, # set to 1e+9
                        randomize = TRUE, 
                        qualityType = "FastqQuality", 
                        errorEstimationFunction = get(err_funs[e]),
                        verbose=TRUE)
    
    # Write out fit error modes
    saveRDS(errF, paste0(runs[i],"_",err_funs[e], "_errF.rds"))
    saveRDS(errR, paste0(runs[i],"_",err_funs[e], "_errR.rds"))
    
    #output error plots to see how well the algorithm modelled the errors in the different runs
    p1 <- plotErrors(errF, nominalQ = TRUE) + ggtitle(paste0(runs[i], " Forward Reads, errorfun:", err_funs[e]))
    p2 <- plotErrors(errR, nominalQ = TRUE) + ggtitle(paste0(runs[i], " Reverse Reads, errorfun:", err_funs[e]))
    pdf(paste0("output/error_models/",runs[i],"_",err_funs[e],"_errormodel.pdf"), width = 11, height = 8 , paper="a4r")
    plot(p1)
    plot(p2)
    try(dev.off(), silent=TRUE)
  }
  
}

