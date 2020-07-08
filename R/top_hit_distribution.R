


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

merged_final_bftrimmed <- blast_top_hit(query=seqs, db="merged_final_bftrimmed.fa")
pruned <- blast_top_hit(query=seqs, db="pruned.fa")
mergedseqs <- blast_top_hit(query=seqs, db="mergedseqs.fa")
codon_filt <- blast_top_hit(query=seqs, db="codon_filt.fa")
purged <- blast_top_hit(query=seqs, db="purged.fa")
resolved <- blast_top_hit(query=seqs, db="resolved.fa")
subset <- blast_top_hit(query=seqs, db="subset.fa")
filtered <- blast_top_hit(query=seqs, db="filtered.fa")
uniqSeqs <- blast_top_hit(query=seqs, db="uniqSeqs.fa")



#  Plotting ---------------------------------------------------------------


# Top hit dist

# 01 - merged
merged <- read_csv("reference/old/mergedseqs.csv")

gg.merged <- merged %>%
  ggplot(aes(x=pident))+ 
  geom_histogram(colour="black") + 
  ggtitle("mergedSeqs tophit")

# 02- unique

# 03 - filtered

# 04 - subset

# 05 - resolved

# 06 - Purged
purged <- read_csv("reference/old/purged.csv")

gg.purged <- purged %>%
  ggplot(aes(x=pident))+ 
  geom_histogram(colour="black") + 
  ggtitle("purged tophit")

# 07 - Codon filt
codon_filt <- read_csv("reference/old/codon_filt.csv")

gg.codon <- codon_filt %>%
  ggplot(aes(x=pident))+ 
  geom_histogram(colour="black") + 
  ggtitle("codon_filt tophit")

# 08 - Pruned
pruned <- read_csv("reference/old/pruned_tophit.csv")

gg.pruned <- pruned %>%
  ggplot(aes(x=pident))+ 
  geom_histogram(colour="black") + 
  ggtitle("Pruned tophit")

# 09 - Final
final <- read_csv("reference/old/merged_final_bftrimmed_tophit.csv")

gg.final <- final %>%
  ggplot(aes(x=pident))+ 
  geom_histogram(colour="black") + 
  ggtitle("final tophit")


library(patchwork)

gg.merged / gg.purged / gg.codon / gg.pruned / gg.final

