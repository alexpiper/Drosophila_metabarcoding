lca_probs <- function(x, method="mbed",  k=5, nstart = 20, ranks=c("kingdom", "phylum", "class", "order", "family", "genus", "species"), delim=";"){
  # Convert to DNAbin
  if (!is(x, "DNAbin")) {
    x <- ape::as.DNAbin(x)
  }
  
  # Replace any trailing delimiters
  names(x) <- names(x) %>%
    stringr::str_remove(";$")
  # Count number of remaining delimiters
  ndelim <- stringr::str_count(names(x)[1], ";")
  
  #check if delims match ranks
  if(any(ndelim <=1)) {
    stop("Error: Needs to be in hierarchial format first, run get_lineage or get_ott_lineage")
  } else if (!ndelim==length(ranks)){
    stop("Error: Number of delimiters does not match number of ranks")
  }
  
  # Create distance matrix - could probably do this using random samples
  if(method == "kdist"){
    dist <- as.matrix(kmer::kdistance(x, k=k, method="edgar", residues="DNA"))
  } else if(method == "mbed"){
    dist <- as.matrix(kmer::mbed(x, k=k, residues="DNA")[,])
  }
  
  #convert to pairwise distances
  xy <- t(combn(colnames(dist), 2))
  pw <- data.frame(xy, dist=(100 - round(dist[xy] * 100)))
  
  # subset to the values in loop
  
  sim <- sort(unique(pw$dist), decreasing = TRUE)
  simlist <- vector("list", length=length(sim))
  s=1
  for (s in 1:length(sim)) {
    
    subsets <- pw %>%
      filter(dist==sim[s])
    
    df1 <- subsets %>%
      tidyr::separate(X1, into=c("Acc",ranks), sep=";")%>%
      select(rev(ranks))
    
    df2 <- subsets %>%
      tidyr::separate(X2, into=c("Acc",ranks), sep=";")%>%
      select(rev(ranks))
    
    #Get all shared ranks
    logidf <- as.data.frame(df1 == df2)  #
    
    #Get lowest common rank
    keepvec <- apply(logidf, 1, which.max)
    rows <- seq(logidf[,1])
    #cols <- seq(logidf[1,])
    #unselect <- matrix(ncol=2,c(rep(rows, length(cols)), sort(rep(cols, length(rows)))))
    select <- matrix(ncol=2,c(rows, keepvec))
    
    logidf[select] <- "KEEP"
    logidf[!logidf=="KEEP"] <- 0
    logidf[logidf=="KEEP"] <- 1
    logidf <- logidf  %>%
      mutate_all(as.numeric) %>%
      colSums() / length(rows)
    
    simlist[[s]]  <- data.frame(rank=names(logidf), prob=logidf, sim=sim[s])
    
  }
  names(simlist) <- sim
  out <- bind_rows(simlist) %>%
    group_by(rank, sim) %>%
    summarise(prob = mean(prob))
  return(out)
}
