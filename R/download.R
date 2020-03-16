## Load packages
## Load Necessary packages
sapply(c("rentrez", "bold", "taxize","taxizedb", "usethis", "tidyverse", "spider", "insect", "ape", "DECIPHER", "ggpubr", "RColorBrewer", "plotly", "ggforce", "seqinr", "shortread", "patchwork", "viridis","ggridges","UpSetR"), require, character.only = TRUE)


library(taxreturn)

## Download data for all insecta

## Fetch sequences from GenBank 
#genbank <- fetchSeqs("Insecta", database="genbank", downstream="Family", quiet=FALSE, marker="COI OR COI OR COX1 OR COXI", output = "gb-binom",compress=FALSE, cores=1)


## Fetch missing
#
#  taxlist <- dplyr::bind_rows(taxize::downstream("Insecta", db = "ncbi", downto = "Family")) %>% 
#    dplyr::filter(rank =="family") %>%
#    pull(childtaxa_name)
#  
#message(length(taxlist))
#message("taxlist done")
#
##already downloaded files
#  dl <- list.files("genbank/") %>%
#    str_split_fixed("_", n=2) %>%
#    as.data.frame() %>%
#    pull(V1) %>%
#    as.character()
#message(length(dl))
#message("dl done")
#
#  
#  notdl <- taxlist[which(!taxlist %in% dl)]
#message(length(notdl))
#genbank <- fetchSeqs(notdl, database="genbank", downstream=FALSE, quiet=FALSE, marker="COI OR COI OR COX1 OR COXI", output = "gb-binom", compress=TRUE, cores=1)
  
## Fetch sequences from BOLD
bold <- fetchSeqs("Insecta", database="bold", downstream="Family",quiet=FALSE, marker="COI-5P", output = "gb-binom",compress=TRUE, cores=1)

## Fetch mitochondrial genomes from genbank
fetchSeqs("Insecta", database="genbank", quiet=FALSE, marker="mitochondria", output = "gb-binom", compress=TRUE, cores=1)


## Download data for Arachnida

## Fetch sequences from GenBank 
fetchSeqs("Arachnida", database="genbank", downstream="Order", quiet=FALSE, output = "gb-binom", compress=TRUE, cores=1)

## Fetch sequences from BOLD
fetchSeqs("Arachnida", database="bold", downstream="Order", quiet=FALSE, marker="COI-5P", output = "gb-binom",compress=FALSE, cores=1)


# Get some outgroups using random subsampling of search

outgroup_classes <- dplyr::bind_rows(taxize::downstream("Arthropoda", db="itis", downto="Class")) %>%
  filter(!taxonname=="Insecta") %>%
  pull(taxonname)
outgroup_phyla <- dplyr::bind_rows(taxize::upstream("Insecta", db="itis", upto="Phylum")) %>%
  pull(taxonname)
outgroup_kingdoms <- c("Bacteria", "Fungi")

outgroups <- c(outgroup_classes, outgroup_phyla, outgroup_kingdoms)

fetchSeqs(outgroups, database="genbank", out.dir="outgroups", quiet=FALSE, output = "gb-binom", subsample = 100, compress=TRUE, cores=1)

