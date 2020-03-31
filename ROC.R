
## ROC CURVE
df <- data.frame(value = c(1,2,3,5,8,4,6,7,9,10),
                 truth = c(rep(0,5), rep(1,5)))

# Value could be threshold value
# Truth: 0 = false, 1=true

library(pROC)
roc.demo <- roc(truth ~ value, data = df)
par(pty = "s") # make it square
plot(roc.demo) # plot ROC curve

roc.demo <- roc(truth ~ value, data = df, percent = T)

roc.demo$thresholds[roc.demo$sensitivities == 0.8]

quantile(df$value[df$truth == 1], 
         probs = c(0.00, 0.10, 0.20, 0.30), type = 1) # percentile giving the closest number


## FUNCTION - create ROC

# inputs:
# a PS object
# A expected table
# a range of filtering thresholds to try default (seq 0,0.1,0.001)
# option to use all taxa or specific taxa
# Return a ROC object


# Read in expected table
exp <- read_csv("sample_data/expected_quant.csv") %>%
  gather(Species, Actual, -X1) %>%
  mutate(Species = str_replace(Species, pattern=" ",replacement="_")) %>%
  #filter(str_detect(X1,"D100M|D250M|D500M|D1000M|DLarv")) %>%
  drop_na() %>%
  set_colnames(c("Sample","Species","Actual")) %>%
  mutate(Species = str_replace(Species, pattern="Drosophila_simulans", replacement="Drosophila_mauritiana/simulans"),
         Sample = str_replace(Sample, pattern="D100M", replacement="DM"))

# PS needs replacing 
#mutate(Species = str_replace(Species, pattern="Drosophila_mauritiana/simulans", replacement= "Drosophila_simulans")) %>%

physeq <- ps.merged

create_roc <- function(physeq, exp, taxrank=NULL, thresholds = seq(0,0.1,0.001)){
  
  # Agglomerate taxa
  if (stringr::str_to_sentence(taxrank) %in% colnames(tax_table(physeq))){
    physeq <-  physeq %>%
      speedyseq::tax_glom(taxrank=taxrank)
  }
  
  # Start filter for loop here
  roc <- vector("list", length=length(thresholds))
  for (i in 1:length(thresholds)){
    
    # Filter
    ps_filt <- physeq %>% 
      transform_sample_counts(function (x) x/sum(x)) %>%
      speedyseq::psmelt() %>%
      filter(Abundance > thresholds[i])
  
    # Convert to pres/abs
    ps_qual <- ps_filt %>%
      mutate(Abundance = case_when(
        Abundance > 0 ~ 1,
        Abundance == 0 ~  0
      )) %>%
      dplyr::select(Sample, Abundance, Species)

    # Get joint
    joint <- ps_qual  %>%
      distinct() %>%
      right_join(exp, by=c("Sample", "Species")) %>%
      mutate(Actual = case_when(
        Actual > 0 ~ 1,
        Actual == 0 ~  0
      )) %>%
      mutate(Abundance = replace_na(Abundance, 0)) %>% 
      mutate(confusion = case_when(
        Abundance > 0 & Actual == 0 ~ "FP", 
        Abundance == 0 & Actual > 0 ~ "FN",
        Abundance > 0 & Actual > 0 ~ "TP",
        Abundance == 0 & Actual == 0 ~ "TN"
      ))
    
    
    roc[[i]] <- joint$confusion
    
    

   #roc[[i]] <- tibble(
   #  thresh = thresholds[i],
   #  tpr = sum(joint$confusion == "TP") / (sum(joint$confusion == "TP") + sum(joint$confusion == "FN")),
   #  fpr = sum(joint$confusion == "FP") / (sum(joint$confusion == "FP") + sum(joint$confusion == "TN"))
   #)
    
  }
  
  test <- bind_cols(roc) %>% 
    set_colnames(thresholds) %>%
    pivot_longer(cols=everything(),
                 names_to = "value",
                 values_to = "truth") %>%
    filter(truth %in% c("TP", "FP")) %>%
    mutate(truth = case_when(
      truth == "TP" ~ 1,
      truth == "FP" ~ 0
    )) %>%
    mutate(value = as.numeric(value))
  library(pROC)
  roc.demo <- roc(truth ~ value, data = test)
  par(pty = "s") # make it square
  plot(roc.demo) # plot ROC curve
  


}

# Need to make a confusion matrix instead

## Plot true and false carsonellas by run - this shows may need to filter seperately by run, as seqrun2 has much higher than others.
gg.switch <- sam %>% 
  group_by(switched) %>%
  filter(!is.na(switched)) %>%
  filter(Abundance > 0) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x=seqrun, y=freq, group=switched, fill=switched)) +
  geom_bar(stat="identity")


#to check if index switching could explain the multiple OTU's check if top=FALSE carsonella exist as top in another sample in the same run

thresholds <- seq(0,0.1,0.0001)

df <- data.frame(thresh = thresholds, TP= thresholds, FP = thresholds)
for(i in 1:length(thresholds)){
  filt <- sam %>%
    filter(Abundance > thresholds[i]) %>%
    group_by(switched) %>%
    summarise(sum=n())
  df$FP[i] <- filt$sum[2]
  df$TP[i] <- filt$sum[1]
}

gg.filt <- df %>%
  bind_rows() %>%
  gather(key=type, value=n, -thresh) %>%
  ggplot(aes(x=thresh, y=n, fill=type)) + 
  geom_density(stat="identity", alpha=0.5) + 
  #scale_x_continuous(breaks=seq(0,0.05,0.001)) +
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  geom_vline(xintercept = 0.0035)+ 
  geom_vline(xintercept = 0.001)

print(gg.filt)

## Filter run 2 
run2_ps <- ps2 %>% 
  subset_samples(seqrun==2)

run2_pass <- run2_ps %>%
  transform_sample_counts(function (x) x/sum(x)) %>%  # Convert to proportions
  transform_sample_counts(function (x) (x > 0.001) * 1)

newotu1 <- otu_table(run2_ps) * otu_table(run2_pass)

run2_newps <- run2_ps
otu_table(run2_newps) <- otu_table(newotu1, taxa_are_rows = FALSE) 

#Filter run 1 and 3 

run13_ps <- ps2 %>%
  subset_samples(!seqrun == 2)

run13_pass <- run13_ps %>%
  transform_sample_counts(function (x) x/sum(x)) %>%  # Convert to proportions
  transform_sample_counts(function (x) (x > 0.0005) * 1)

newotu2  <- otu_table(run13_ps) * otu_table(run13_pass)

run13_newps <- run13_ps
otu_table(run13_newps) <- otu_table(newotu2, taxa_are_rows = FALSE) 

# Create new phyloseq and drop missing taxa
ps3 <- merge_phyloseq(run2_newps, run13_newps) %>%
  filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
ps3 <- prune_samples(sample_sums(ps3) >0 , ps3) # Drop empty samples

#Count number of overall taxa pre and post filtering
print(paste(ntaxa(ps2) - ntaxa(ps3), " taxa Dropped when using filtering threshold of: ", 0.0035, " for run 2 and ", 0.001, "for run 1 and 3"))

# Count number of carsonella OTU's per sample pre and post filtering
n_carson <- speedyseq::psmelt(ps2) %>% 
  filter(Genus == "Candidatus_Carsonella") %>%
  filter(Abundance > 0) %>%
  dplyr::select(Sample.Name, OTU) %>%
  group_by(Sample.Name) %>%
  add_tally() %>%
  dplyr::select(-OTU) %>%
  unique() %>%
  mutate(type = "pre") %>%
  bind_rows(.,speedyseq::psmelt(ps3) %>% 
              filter(Genus == "Candidatus_Carsonella") %>%
              filter(Abundance > 0) %>%
              dplyr::select(Sample.Name, OTU) %>%
              group_by(Sample.Name) %>%
              add_tally() %>%
              dplyr::select(-OTU) %>%
              unique() %>%
              mutate(type = "post"))

gg.ncarson <- ggplot(n_carson, aes(x=reorder(Sample.Name, -n), y=n)) + 
  geom_bar(stat="identity") +
  facet_grid(~type) + 
  xlab("Sample.Name") +
  ylab("Number of carsonella OTU's per sample") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0))

print(gg.ncarson)

#check if any samples dont have carsonella
table(!n_carson$Sample.Name %in% speedyseq::psmelt(ps2)$Sample.Name)

## Get the name of taxa that dont have carsonella after filtering
