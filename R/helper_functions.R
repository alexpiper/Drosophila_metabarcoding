## Parralel slidenucdiag
slideNucDiag_para <- function (DNAbin, sppVector, width, interval = 1, cores=1) {
  
  # Define parralel nucdiag function
  nucDiag_para <- function (DNAbin, sppVector, cores=1) {
    DNAbin <- as.matrix(DNAbin)
    inform <- seg.sites(DNAbin)
    sppSeqs <- lapply(unique(sppVector), function(x) which(sppVector == x))
    # Define sitecheck
    siteCheck <- function(DNAbin, spp, inform) {
      res <- vector("logical", length=length(inform))
      for (j in 1:length(inform)) {
        site <- inform[j]
        res[j] <- as.character(DNAbin[spp, site]) %in% as.character(DNAbin[-spp, site])
        res[j] <- as.logical(sum(as.numeric(res[j])))
      }
      return(res)
    }
    li <- vector("list", length=length(sppSeqs))
    if (cores > 1){
      cl <- parallel::makeCluster(cores)
      registerDoParallel(cl)
      li <- foreach(i=1:length(sppSeqs)) %dopar% {
        siteCheck(DNAbin, sppSeqs[[i]], inform=inform)
      } 
      parallel::stopCluster(cl)
      
    } else {
      for (i in 1:length(sppSeqs)) {
        li[[i]] <- NA
        for (j in 1:length(inform)) {
          li[[i]][j] <- siteCheck(sppSeqs[[i]], inform[j])
        }
      }
    }
    out <- lapply(li, function(x) inform[which(!x)])
    names(out) <- unique(sppVector)
    return(out)
  }
  
  nd <- nucDiag_para(DNAbin, sppVector, cores=cores)
  if (interval == "codons") 
    interval <- 3
  win <- seq(1, dim(DNAbin)[2] - width, by = interval) #Get all possible windows
  mat <- matrix(NA, nrow = length(nd), ncol = length(win))
  
  #setup parallel backend to use many processors
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    registerDoParallel(cl)
    mat <- foreach(i=1:length(win), .combine = cbind) %dopar% {
      sapply(nd, function(x) length(which(x %in% win[i]:(win[i] + width))))
    }
    parallel::stopCluster(cl)
  } else{
    for (i in 1:length(win)) {
      mat[, i] <- sapply(nd, function(x) length(which(x %in% win[i]:(win[i] + width))))
    }
  }
  dimnames(mat)[[1]] <- unique(sppVector)
  dimnames(mat)[[2]] <- NULL
  return(mat)
}