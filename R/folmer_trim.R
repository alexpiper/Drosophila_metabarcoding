sapply(c("rentrez", "bold", "taxize","taxizedb", "usethis", "tidyverse", "spider", "insect", "ape", "DECIPHER", "ggpubr", "RColorBrewer", "plotly", "ggforce", "seqinr", "shortread", "patchwork", "viridis","ggridges","UpSetR"), require, character.only = TRUE)

#build PHMM from midori longest - sequences need to be same length
midori <-  Biostrings::readDNAStringSet("reference/MIDORI_LONGEST_20180221_COI.fasta")
insecta_midori <- as.DNAbin(midori[str_detect(names(midori),pattern=";Insecta;"),])
folmer <- insect::virtualPCR(insecta_midori, up = "TITCIACIAAYCAYAARGAYATTGG",down= "TAIACYTCIGGRTGICCRAARAAYCA",cores=2, rcdown = TRUE, trimprimers = FALSE)
insect::writeFASTA(folmer, "folmer_full.fa")