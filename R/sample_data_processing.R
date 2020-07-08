
# Samdf processing --------------------------------------------------------

samdf <- read.csv("sample_data/sample_info.csv", header=TRUE) %>% 
  janitor::clean_names() %>%
  mutate(sample_id = case_when(
    fcid=="HLVKYDMXX" ~ paste0(sample_name, "_", replicate),
    !fcid=="HLVKYDMXX" ~ paste0(fcid, "_", sample_name, "_", replicate)
  )) %>%
  mutate(extract_id = sample_name) %>%
  mutate(sample_name = str_remove(sample_name, "-exp*.$")) %>%
  filter(!(index=="ATCGATCG" & index2=="ATCACACG"), #CT11-ex1 duplicated
         !(index=="TCGCTGTT" & index2=="ACTCCATC") # CT12-ex1 duplicated
  ) %>%
  filter(!fcid=="CK3HD") %>%
  mutate(type = case_when(
    str_detect(sample_id, "D[0-9][0-9][0-9]M|D[0-9][0-9][0-9][0-9]M|DM[0-9]")  ~ "DrosMock",
    str_detect(sample_id, "SPD")  ~ "SPD",
    str_detect(sample_id, "ACV")  ~ "ACV",
    str_detect(sample_id, "DC")  ~ "DC",
    str_detect(sample_id, "Sach")  ~ "Sachet",
    str_detect(sample_id, "FF")  ~ "FF",
    str_detect(sample_id, "NTC")  ~ "NTC",
    str_detect(sample_id, "DLarv")  ~ "DrosLarv",
    str_detect(sample_id, "POS|SynMock")  ~ "POS",
    str_detect(sample_id, "extblank|BLANK")  ~ "Extblank",
    str_detect(sample_id, "pcrblank")  ~ "PCRblank",
    str_detect(sample_id, "CT")  ~ "CarpTrap",
    str_detect(sample_id, "CM[0-9]")  ~ "CarpMock",
    str_detect(sample_id, "CML[0-9]")  ~ "CarpLarval"
  )) %>%
  mutate(target_subfragment = case_when(
    str_detect(fprimer, "GGDACWGGWTGAACWGTWTAYCCHCC") & str_detect(rprimer, "GTRATWGCHCCDGCTARWACWGG") ~ "fwhF2-fwhR2n",
    str_detect(fprimer, "ACWGGWTGRACWGTNTAYCC") & str_detect(rprimer, "ARYATDGTRATDGCHCCDGC") ~ "BF1-BR1",
    str_detect(fprimer, "GGDRCWGGWTGAACWGTWTAYCCNCC") & str_detect(rprimer, "TATDGTRATDGCHCCNGC") ~ "SauronS878-HexCOIR4",
    str_detect(fprimer, "GGDACWGGWTGAACWGTWTAYCCHCC") & str_detect(rprimer, "TATDGTRATDGCHCCNGC") ~ "fwhF2-HexCOIR4",
  )) %>%
  magrittr::set_rownames(.$sample_id)

write_csv(samdf, "sample_data/sample_info2.csv")
