ps_traps <- readRDS("ps_traps.rds")

library(mvabund)
library(phyloseq)
library(tidyverse)

# Get species occuracne table
otutab <- otu_table(ps_traps) %>%
  as.matrix() %>%
  as.data.frame()

# Get metadata 
metadata <- sample_data(ps_traps) %>%
  as("data.frame")%>%
  mutate(orchard_type = paste0(orchard, "_", type)) %>%
  dplyr::select(sample_id, orchard, type, orchard_type)


spp_occur <- mvabund(otutab)

# Mean - variance relationship
meanvar.plot(spp_occur)

model_res <- vector("list", length=4)

# Model 1 - trap type
model_res[[1]] <- manyglm(spp_occur ~ metadata$type, family="negative.binomial")

# Model 2 - orchard
model_res[[2]] <- manyglm(spp_occur ~ metadata$orchard, family="negative.binomial")

# Model 3 - orchard + type
model_res[[3]] <- manyglm(spp_occur ~ metadata$orchard+metadata$type, family="negative.binomial")

# Model 4 - orchard * type
model_res[[4]] <- manyglm(spp_occur ~ metadata$orchard*metadata$type, family="negative.binomial")

saveRDS(model_res, "model_res.rds")


# Run multivariate anovas
anova_res <- vector("list", length=length(model_res))

anova_res[[1]] <- anova(model_res[[1]])
anova_res[[2]] <- anova(model_res[[2]])
anova_res[[3]] <- anova(model_res[[3]])
anova_res[[4]] <- anova(model_res[[4]])

saveRDS(anova_res, "anova_res.rds")

# run univariate tests 
univa_res <- vector("list", length=length(model_res))

univa_res[[1]] <- anova(model_res[[1]], p.uni="adjusted")
univa_res[[2]] <- anova(model_res[[2]], p.uni="adjusted")
univa_res[[3]] <- anova(model_res[[3]], p.uni="adjusted")
univa_res[[4]] <- anova(model_res[[4]], p.uni="adjusted")

saveRDS(univa_res, "univa_res.rds")
