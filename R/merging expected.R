
dros_exp <- read_csv("sample_data/expected_quant2.csv") %>%
  dplyr::rename(Sample = X1)%>%
  pivot_longer(-Sample,
               names_to= "Species",
               values_to= "Abundance") %>%
  mutate(Species = str_replace(Species, pattern=" ",replacement="_")) 


carp_exp <- read_csv("sample_data/expected_quant_carpophilus.csv")%>%
  dplyr::rename(Sample = X1) %>%
  mutate(Sample = Sample %>%
           str_remove("-ex[0-9]") %>%
           str_remove("fwhF2-fwhR2n-|fwhF2-HexCOIR4-|SauronS878-HexCOIR4-|BF1-BR1-")) %>%
  distinct()%>%
  pivot_longer(-Sample,
               names_to= "Species",
               values_to= "Abundance") %>%
  mutate(Species = str_replace(Species, pattern=" ",replacement="_"))

merged <- bind_rows(dros_exp, carp_exp) %>%
  distinct() %>%
  pivot_wider(names_from = Species,
              values_from = Abundance,
              values_fill = list(Species=0)) %>%
  filter(!str_detect(Sample, "-s$"),
         !is.na(Sample)) %>%
  mutate_at(c(2:ncol(.)), ~replace(., is.na(.), 0))

write_csv(merged, "sample_data/expected_quant_merged.csv")  
