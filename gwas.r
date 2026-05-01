# Title: Publication Composite Figure 2 (GWAS)
# Description: Generates the 5 plots with CORRECTED phenotype mappings.
#              Uses Blue/Yellow alternating chromosome colors.

library(data.table)
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)
library(stringr)
library(patchwork)

ASSOC_ROOT <- "/Users/an/Desktop/Discovery data/"
FILE_PATTERN <- "chr{chr}_pheno{pheno}_lmm.assoc.txt"
OUTPUT_DIR <- "/Users/an/Desktop/PLOTS_PAPERSTYLE/"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# THE CORRECTED PHENOTYPE MAPPING
# ==============================================================================
PHENO_MAP <- c(
  "1" = "Nicotine Composite Score",
  "2" = "Alcohol Consumption Composite Score",
  "3" = "Alcohol Dependence Composite Score",
  "4" = "Illicit Drug Composite Score",
  "5" = "Behavioral Disinhibition Composite Score"
)

# CHANGED: Alternating Blue and Yellow hex codes
PLOT_COLORS <- c("#1f78b4", "#ffc125") 
SIGNIFICANCE_THRESHOLD <- 6

# ==============================================================================
# THE MASTER DICTIONARY (Re-aligned to match the correct traits)
# ==============================================================================
annotation_mapping <- tribble(
  ~PhenotypeID, ~CoordinateString, ~GeneSymbol,
  
  # Pheno 1 (Nicotine) 
  "1", "3:29195780",  "EOMES",       
  
  # Pheno 2 (Alcohol Consumption)
  "2", "7:153433844", "DPP6",       
  "2", "22:50725237", "TOP6BL",     
  "2", "22:50724593", "ADAM32",     
  
  # Pheno 3 (Alcohol Dependence)
  "3", "2:180554440", "TNS1", 
  "3", "3:139986727", "NMNAT3", 
  "3", "5:148871502", "HTR4",       
  "3", "13:48927587", "RCBTB2",
  "3", "22:50725237", "ADAM32",   
  "3", "22:50724593", "TOP6BL",   
  
  # Pheno 4 (Illicit Drug) 
  "4", "2:169756930", "LRP1B",      
  "4", "3:54112290",  "ERC2",       
  "4", "7:153433844", "DPP6",       
  "4", "13:48927587", "RCBTB2",     
  "4", "17:10702703", "ABCA8",
  
  # Pheno 5 (Behavioral Disinhibition)
  "5", "2:180554440", "TNS1", 
  "5", "3:139986727", "NMNAT3", 
  "5", "5:148871502", "HTR4",       
  "5", "13:48927587", "RCBTB2"
)

create_manhattan <- function(pheno_id, title_text, show_x_axis = TRUE) {
  message(sprintf("  -> Processing Phenotype %s: %s...", pheno_id, title_text))
  gwas_raw_all <- data.frame()
  
  for (c in 1:22) {
    filename <- str_glue(FILE_PATTERN, chr = c, pheno = pheno_id)
    full_path <- file.path(ASSOC_ROOT, filename)
    
    if (file.exists(full_path)) {
      this_chr <- fread(full_path, data.table = FALSE)
      actual_cols <- colnames(this_chr)
      
      if ("chr" %in% actual_cols) names(this_chr)[names(this_chr) == "chr"] <- "CHR"
      if ("ps" %in% actual_cols) names(this_chr)[names(this_chr) == "ps"] <- "BP"
      if ("pos" %in% actual_cols) names(this_chr)[names(this_chr) == "pos"] <- "BP"
      
      if ("rs" %in% actual_cols) names(this_chr)[names(this_chr) == "rs"] <- "SNP"
      else if (!"SNP" %in% actual_cols) this_chr$SNP <- NA 
      
      if ("p_lrt" %in% actual_cols) names(this_chr)[names(this_chr) == "p_lrt"] <- "P"
      else if ("p_wald" %in% actual_cols) names(this_chr)[names(this_chr) == "p_wald"] <- "P"
      
      this_chr <- this_chr %>%
        mutate(CHR = as.integer(CHR), BP = as.numeric(BP), P = as.numeric(P), SNP = as.character(SNP)) %>%
        mutate(CoordKey = paste0(CHR, ":", sprintf("%d", as.integer(BP)))) %>%
        mutate(DefaultLabel = ifelse(is.na(SNP) | SNP == "" | SNP == ".", CoordKey, SNP)) %>%
        select(CHR, BP, P, CoordKey, DefaultLabel) 
      
      gwas_raw_all <- bind_rows(gwas_raw_all, this_chr)
    }
  }
  
  gwas_raw_all <- gwas_raw_all %>% filter(!is.na(P) & P > 0)
  gwas_dat <- gwas_raw_all %>% mutate(logP = -log10(P))
  nCHR <- length(unique(gwas_dat$CHR))
  
  gwas_dat_cumulative <- gwas_dat %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas_dat, ., by = "CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot)
  
  axis_set <- gwas_dat_cumulative %>% group_by(CHR) %>% summarize(center = mean(range(BPcum))) 
  
  current_dict <- annotation_mapping %>% filter(PhenotypeID == pheno_id)
  sig_snps <- gwas_dat_cumulative %>% filter(logP > SIGNIFICANCE_THRESHOLD) %>% left_join(current_dict, by = c("CoordKey" = "CoordinateString"))
  
  dict_matches <- sig_snps %>% filter(!is.na(GeneSymbol)) %>% mutate(FinalLabel = GeneSymbol)
  chrs_with_dict <- unique(dict_matches$CHR)
  fallback_snps <- sig_snps %>% filter(! CHR %in% chrs_with_dict) %>% group_by(CHR) %>% slice_max(order_by = logP, n = 1, with_ties = FALSE) %>% ungroup() %>% mutate(FinalLabel = DefaultLabel)
  
  top_snps <- bind_rows(dict_matches, fallback_snps)
  
  # Brute force label cleanup (Universal)
  top_snps <- top_snps %>% mutate(FinalLabel = case_when(
    FinalLabel == "22:50725237" ~ "ADAM32", FinalLabel == "22:50724593" ~ "TOP6BL",
    FinalLabel == "7:153433844" ~ "DPP6", FinalLabel == "3:29195780"  ~ "EOMES",
    FinalLabel == "2:180554440" ~ "TNS1", FinalLabel == "3:139986727" ~ "NMNAT3", 
    FinalLabel == "5:148871502" ~ "HTR4", FinalLabel == "13:48927587" ~ "RCBTB2",
    FinalLabel == "2:169756930" ~ "LRP1B", FinalLabel == "3:54112290"  ~ "ERC2",             
    FinalLabel == "17:10702703" ~ "ABCA8", TRUE ~ FinalLabel
  ))
  
  man_plot <- ggplot(gwas_dat_cumulative, aes(x = BPcum, y = logP)) +
    theme_minimal(base_size = 14) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(color = "lightgrey", linetype = "solid"), panel.grid.minor.y = element_blank(), axis.text.x = element_text(size = 10, vjust = 0.5), plot.title = element_text(hjust = 0.5, face = "bold", size = 16), axis.title.y = element_text(face="bold")) +
    labs(title = title_text, x = "Chromosome", y = expression(bold(paste("-log"[10], "(P)")))) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(PLOT_COLORS, nCHR)) +
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center, expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(gwas_dat_cumulative$logP) + 1.5)) +
    geom_hline(yintercept = SIGNIFICANCE_THRESHOLD, color = "red", linetype = "dashed", size = 0.8) +
    theme(legend.position = "none") +
    geom_label_repel(data = top_snps, aes(label = FinalLabel), fontface = "bold", color = "black", fill = "white", box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"), segment.color = "grey50", size = 4, label.padding = unit(0.2, "lines"), label.r = unit(0.1, "lines"), label.size = 0.25, min.segment.length = 0, max.overlaps = Inf)
  
  if (!show_x_axis) man_plot <- man_plot + theme(axis.title.x = element_blank())
  return(man_plot)
}

message("\nGenerating CORRECTED GWAS Composite (Blue/Yellow)...")
plot_list <- list()
for (i in 1:5) { pheno_id <- as.character(i); plot_list[[i]] <- create_manhattan(pheno_id, PHENO_MAP[pheno_id], show_x_axis = (i == 5)) }

composite_figure <- (plot_list[[1]] / plot_list[[2]] / plot_list[[3]] / plot_list[[4]] / plot_list[[5]]) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face = "bold")) 
output_file <- file.path(OUTPUT_DIR, "Figure_2_Composite_Manhattan_CORRECTED_MAPPING.png")
ggsave(filename = output_file, plot = composite_figure, width = 14, height = 25, units = "in", dpi = 300, limitsize = FALSE)
message(sprintf("Success! Saved to: %s", output_file))