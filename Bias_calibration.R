devtools::install_github("mikemc/2019-bias-manuscript")
library(BiasManuscript)
library(tidyverse)

st0 <- read_csv("output/csv/filtered/all_genglom_filt.csv")

stcoi <- st0 %>%
  filter(biome=="Laboratory") %>%
  filter(feature !="Pool-0250") %>%
  filter(Abundance > 0) %>%
  filter(Genus !="NA") %>%
  filter(loci =="COI-Eukaryota")

stcoi  %>%
  arrange(Abundance)

exp <- read_csv("sample_data/expected/Test_expected_quant.csv") %>%
  gather(Genus,Abundance,-X1)

colnames(exp) <- c("Sample","Genus","Abundance")

stm <- stcoi %>%
  left_join(exp, by = c("Sample","Genus"))%>%
  mutate_by(Sample, Actual = close_elts(Abundance.y))%>%
  rename(Observed = Abundance.x) %>%
  rename(Taxon = Genus)


# Estimate bias

bias <- stm %>%
  nest %>%
  mutate(Bias = map(data, estimate_bias, method = "rss"))%>%
  unnest(Bias) 

#COmpute predicted proportions

stm <- stm %>%
  left_join(bias, by = c("Taxon")) %>%
  mutate_by(Sample, Predicted = close_elts(Actual * Bias_est))

stm <- stm %>% mutate_by(Sample, Proportion=close_elts(Observed))

#Plotting

library(ggpubr)
library(patchwork)
gg.act <- ggplot(stm, aes(x= Sample, y=Actual,fill= Taxon)) + 
  geom_bar(stat = "identity", position = "stack", color = "NA")  + 
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title=element_text(hjust = 0.5),
        legend.position = "none") + 
  ggtitle(paste0("Expected COI")) + 
  scale_fill_manual(values=c("#0c4687","#ae0707","#fa6e24","#3a9e82","#95cf77"))+
  coord_flip()

gg.obs <- ggplot(stm, aes(x= Sample, y=Proportion,fill= Taxon)) + 
  geom_bar(stat = "identity", position = "stack", color = "NA")  + 
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title=element_text(hjust = 0.5),
        legend.position = "none") + 
  ggtitle(paste0("Observed COI")) + 
  scale_fill_manual(values=c("#0c4687","#ae0707","#fa6e24","#3a9e82","#95cf77"))+
  coord_flip()

gg.pre <- ggplot(stm, aes(x= Sample, y=Predicted,fill= Taxon)) + 
  geom_bar(stat = "identity", position = "stack", color = "NA")  + 
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title=element_text(hjust = 0.5),
        legend.position = "none") + 
  ggtitle(paste0("Calibrated COI")) + 
  scale_fill_manual(values=c("#0c4687","#ae0707","#fa6e24","#3a9e82","#95cf77"))+
  coord_flip()

#Multiplot

#Create final figure 1 by stitching all the subfigures together using patchwork
gg.act + gg.obs + gg.pre + plot_layout(ncol = 3)



#Plot expected vs observed

correction <- stm %>% select(Sample,Taxon,Observed,Proportion,Actual,Predicted)

correction <- correction %>%
  gather(Type,Abundance,-Sample,-Taxon,-Actual,-Proportion)

g.cor <- ggplot(correction, aes(x=Actual,y=Abundance)) + geom_point(aes(color=Type)) + geom_abline(slope=1, intercept = 0)+
  stat_cor(aes(color=Type), label.x = 0.1)  + xlim(0,1) + ylim(0,1) + scale_colour_manual(values=c("#A9A9A9","#ae0707")) + theme_pubr()

#Might be worth trying this again with all genes in together?
