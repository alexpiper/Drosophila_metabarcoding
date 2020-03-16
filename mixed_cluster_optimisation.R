# Input seqs
input <- "reference/subset.fa.gz"
message(input)


#Load libraries
library(taxreturn)
library(aphid)
library(insect)
library(Biostrings)
library(ape)
library(pbapply)
library(parallel)
library(tidyverse)


seqs <-  insect::readFASTA(input)
seqs <- insect::subset.DNAbin(seqs, subset = !duplicated(names(seqs)))

cat("seqs read \n")

thresholds <- rev(seq(0.9, 1, 0.01))
cat(paste0("thresholds:", thresholds," \n"))

threshlist <- vector("list", length=length(thresholds))

db <- get_ncbi_lineage()
for (i in 1:length(thresholds)){
  threshlist[[i]] <- taxreturn::get_mixed_clusters(x = seqs, db=db, rank = c("species","genus","family"), threshold = thresholds[i], confidence=0.6, quiet = FALSE) 
}

names(threshlist) <- thresholds 
results <- bind_rows(threshlist)
saveRDS(results, str_replace(input, ".fa.gz", "_mixedclusters.RDS"))
write.csv(results, str_replace(input, ".fa.gz", "_mixedclusters.csv"))

## Read in results - check if the resolve_twaxoomy function improves or makes things worse
subset_mixed <- read_csv("reference/optimisation/subset_mixedclusters.csv")

resolved_mixed <- read_csv("reference/optimisation/resolved_mixedclusters.csv")

## Count at each rank

gg.resolved <- resolved_mixed %>%
  group_by(threshold, rank) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=threshold, y=n, fill=rank, group=rank)) + 
  geom_bar(stat="identity", position="dodge") +
  ggtitle("Resolved")

gg.subset <- subset_mixed %>%
  group_by(threshold, rank) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=threshold, y=n, fill=rank, group=rank)) + 
  geom_bar(stat="identity", position="dodge") +
  ggtitle("Subset")

# Need to check the resolving code. Does it use the taxonomy code or the name? if its the taxonomy code should probably purge before resolving names?

x <- seqs[1:2000]
## need to do a subset_mixed > resolved


