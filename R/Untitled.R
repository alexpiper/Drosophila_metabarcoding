
get_primer_statistics <- function(x, metrics = "all", disambiguate = TRUE){
  #Get metrics
  availmetrics <- c("base_freq","clamps",
                   "tm", "homopolymer",
                   "degeneracy", "primer_length")
  if (length (metrics)== 1 && metrics =="all"){
    metrics <- availmetrics
  }
  if(any(!metrics %in% availmetrics )){
    stop("Error, invalid metric: ", metrics[!metrics %in% availmetrics])
  }

  # Replace inosines
  x <- x %>% stringr::str_replace_all("I", "N")
  if(disambiguate){
    query <- unlist(DECIPHER::Disambiguate(DNAStringSet(x)))
  } else {
    query <- DNAString(x)
  }
  if("base_freq" %in% metrics){
    base_freq<- letterFrequency(query,letters="ACGT", OR=0) %>%
    prop.table() %>%
    colSums() %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(`GC%` = (G+C))
  }
  if ("clamps" %in% metrics){ # Should ambiguous bases count?
    clamps <- data.frame(
      GC_last2 = x %>% str_extract("..$") %>% str_count("G|C"),
      GC_last5 = x %>% str_extract(".....$") %>% str_count("G|C")
    )
  }
  if ("tm" %in% metrics){
      tm <- data.frame(tm= TmCalculator::Tm_NN(x, ambiguous = disambiguate))
  }  
  if("homopolymer" %in% metrics){
    homopolymer <- data.frame(
        `poly_A` = query %>% as.character() %>% purrr::map_chr(longestConsecutive, "A") %>% max(),
        `poly_T` = query %>% as.character() %>%purrr::map_chr(longestConsecutive, "T") %>% max(),
        `poly_G` = query %>% as.character() %>% purrr::map_chr(longestConsecutive, "G") %>% max(),
        `poly_C` = query %>% as.character() %>% purrr::map_chr(longestConsecutive, "C") %>% max(),
        stringsAsFactors = FALSE
      )
  }
  if("degeneracy" %in% metrics){
    degeneracy <- data.frame(degeneracy = query %>% DECIPHER::Disambiguate() %>% length())
  } 
  if("primer_length" %in% metrics){
    primer_length <- data.frame(length = x %>% nchar())
  } 

  out <- metrics %>% 
    purrr::map(get) %>%
    bind_cols
}

