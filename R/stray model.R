# Stray modelling
devtools::install_github("twbattaglia/MicrobeDS")

library(MicrobeDS) #Microbiome example datasets
library(phyloseq)
library(dplyr)

devtools::install_github("jsilve24/stray")
library(stray)
set.seed(899)

data("RISK_CCFA")
# drop low abundant taxa and samples
dat <- RISK_CCFA %>% 
  subset_samples(disease_stat!="missing", 
                 immunosup!="missing") %>% 
  subset_samples(diagnosis %in% c("no", "CD")) %>% 
  subset_samples(steroids=="false") %>% 
  subset_samples(antibiotics=="false") %>% 
  subset_samples(biologics=="false") %>% 
  subset_samples(biopsy_location=="Terminal ileum") %>% 
  speedyseq::tax_glom("Family") %>% 
  prune_samples(sample_sums(.) >= 5000,.) %>%
  filter_taxa(function(x) sum(x > 3) > 0.10*length(x), TRUE)

# Create desing matrix
sample_dat <- as.data.frame(as(sample_data(dat),"matrix"), stringsAsFactors=TRUE) %>% 
  mutate(age = as.numeric(as.character(age)),
         diagnosis = relevel(diagnosis, ref="no"), 
         disease_stat = relevel(disease_stat, ref="non-inflamed"))
X <- t(model.matrix(~diagnosis + disease_stat+age, data=sample_dat))
Y <- otu_table(dat)

# Specify priors
upsilon <- ntaxa(dat)+3 
Omega <- diag(ntaxa(dat))
G <- cbind(diag(ntaxa(dat)-1), -1)
Xi <- (upsilon-ntaxa(dat))*G%*%Omega%*%t(G)

Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
Gamma <- diag(nrow(X))                             
            
priors <- stray::pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
print(priors)

# Convert into CLR coordinate system
priors <- to_clr(priors)  
summary(priors, pars="Lambda")  

# Plot priors
names_covariates(priors) <- rownames(X)
plot(priors, par="Lambda") + ggplot2::xlim(c(-10, 10))  

# Fit the data
priors$Y <- Y # remember pibblefit objects are just lists so can add stuff - Y is the OTUtable!
posterior <- refit(priors, optim_method="adam")

# Add the taxon names
tax <- tax_table(dat)[,c("Class", "Family")]
tax <- apply(tax, 1, paste, collapse="_")
names_categories(posterior) <- tax

#plot the prior predictive distribution
ppc(posterior) + ggplot2::coord_cartesian(ylim=c(0, 30000))
ppc_summary(posterior)

# plot again using fromscratch
ppc(posterior, from_scratch=TRUE) +ggplot2::coord_cartesian(ylim=c(0, 30000))
ppc_summary(posterior, from_scratch=TRUE)

# look at the posterior
posterior_summary <- summary(posterior, pars="Lambda")$Lambda
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
focus <- unique(focus$coord)
plot(posterior, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[2:4])


# Plot just the coordinates with non-zero effect for diagnosis CD
posterior_summary <- filter(posterior_summary, covariate=="diagnosisCD") 
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
focus <- unique(focus$coord)

tax_table(dat)[taxa_names(dat)[which(names_coords(posterior) %in% focus)]]
plot(posterior, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[2])
