
process_errors <- function(dq){
  nti <- c("A", "C", "G", "T")
  ntj <- c("A", "C", "G", "T")
  dq <- getErrors(dq, detailed = TRUE, enforce = FALSE)
  
  transdf = reshape2::melt(dq$trans, factorsAsStrings = TRUE)
  colnames(transdf) <- c("Transition", "Qual", "count")
  
  transdf$from <- substr(transdf$Transition, 1, 1)
  transdf$to <- substr(transdf$Transition, 3, 3)
  
  tot.count <- tapply(transdf$count, list(transdf$from, transdf$Qual), sum)
  
  transdf$tot <- mapply(function(x, y) tot.count[x, y], 
                        transdf$from, as.character(transdf$Qual))
  transdf$Observed <- transdf$count/transdf$tot
  
  transdf$Estimated <- mapply(function(x, y) dq$err_out[x, y], transdf$Transition, as.character(transdf$Qual))
  
  ei <- dq$err_in
  ei <- ei[[1]]
  transdf$Input <- mapply(function(x, y) ei[x, y], transdf$Transition, 
                          as.character(transdf$Qual))
  
  transdf$Nominal <- (1/3) * 10^-(transdf$Qual/10)
  transdf$Nominal[transdf$Transition %in% 
                    c("A2A", "C2C", "G2G", "T2T")] <- 1 - 10^-(
                      transdf$Qual[transdf$Transition %in% c("A2A", "C2C", "G2G", "T2T")]/10)
  out <- transdf[transdf$from %in% nti & transdf$to %in% ntj, ]
  return(out)
}

error <- c(
  "output/error_models/CB3DR_loessErrfun_errF.rds",
  "output/error_models/CB3DR_loessErrfun_mod4_errF.rds",
  "output/error_models/HLVKYDMXX_loessErrfun_errF.rds",
  "output/error_models/HLVKYDMXX_loessErrfun_mod4_errF.rds"
            ) %>%
  purrr::set_names() %>%
  purrr::map(readRDS) %>%
  purrr::map(process_errors) %>%
  bind_rows(.id = "source") %>%
  dplyr::mutate(
    fcid = source %>% str_remove("^.*/") %>% str_remove("_(.*?)$"),
    model = source %>% str_remove("^.*/") %>% str_replace("loessErrfun_mod4", "mod4") %>% str_remove("^(.*?)_")%>% str_remove("_.*$"),
    dir = source %>% str_remove("^.*_err") %>% str_remove(".rds")
                ) %>%
  mutate(fcid = fcid %>% str_replace("CB3DR", "MiSeq") %>% str_replace("HLVKYDMXX", "NovaSeq"),
         model = model %>% str_replace("loessErrfun", "Original") %>% str_replace("mod4", "Modified"),
         model = paste0(fcid, " ", model)
         ) %>%
  mutate(model = factor(model, levels = c("MiSeq Original", "MiSeq Modified", "NovaSeq Original", "NovaSeq Modified")))

# Map a function to extract error information

# Plot with Read direction as different colours
gg.errors <- ggplot(data = error, aes(x = Qual, group=Transition, colour=fcid)) +
  geom_point(aes(y = Observed), na.rm = TRUE) + 
  geom_line(aes(y = Estimated), colour="black") +
  scale_y_log10()+
  scale_color_brewer(palette="Set1", direction = -1) +
  facet_grid(Transition~model) +
  xlab("Consensus quality score") + 
  ylab("Error frequency (log10)")+
  base_theme+
  theme(axis.text.x = element_text(angle=0,hjust = 0.5))
  #theme_bw()

gg.errors


pdf(file="fig/supplementary/error_models.pdf", width = 8, height = 11, paper="a4")
  plot(gg.errors)
try(dev.off(), silent=TRUE)
  
  