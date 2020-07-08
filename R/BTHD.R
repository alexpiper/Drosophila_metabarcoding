
# Load libraries
library(tidyverse)
library(insect)
library(taxreturn)


# top hit distribution ----------------------------------------------------

seqtab <- readRDS("run4_HLVKYDMXX_seqtab.rds")

seqs <- insect::char2dna(colnames(seqtab))
  

# make blast Databases
makeblastdb("merged_final_bftrimmed.fa.gz")
makeblastdb("pruned.fa.gz")
makeblastdb("mergedseqs.fa.gz")
makeblastdb("codon_filt.fa.gz")
makeblastdb("purged.fa.gz")
makeblastdb("resolved.fa.gz")
makeblastdb("subset.fa.gz")
makeblastdb("filtered.fa.gz")
makeblastdb("uniqSeqs.fa.gz")

# Top hit distribution

merged_final_bftrimmed <- blast_top_hit(query=seqs, db="merged_final_bftrimmed.fa", threshold=50)
write_csv(merged_final_bftrimmed, "merged_final_bftrimmed_tophit.csv")

pruned <- blast_top_hit(query=seqs, db="pruned.fa", threshold=50)
write_csv(pruned, "pruned_tophit.csv")

mergedseqs <- blast_top_hit(query=seqs, db="mergedseqs.fa", threshold=50)
write_csv(mergedseqs, "mergedseqs.csv")

codon_filt <- blast_top_hit(query=seqs, db="codon_filt.fa", threshold=50)
write_csv(codon_filt, "codon_filt.csv")

purged <- blast_top_hit(query=seqs, db="purged.fa", threshold=50)
write_csv(purged, "purged.csv")

resolved <- blast_top_hit(query=seqs, db="resolved.fa", threshold=50)
write_csv(resolved, "resolved.csv")

subset <- blast_top_hit(query=seqs, db="subset.fa", threshold=50)
write_csv(subset, "subset.csv")

filtered <- blast_top_hit(query=seqs, db="filtered.fa", threshold=50)
write_csv(filtered, "filtered.csv")
uniqSeqs <- blast_top_hit(query=seqs, db="uniqSeqs.fa", threshold=50)

