# POINT MUTATION & CNV MULTI-PANEL FIGURES

#Importing Necessary Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(scales)
library(patchwork)

setwd("PathTofiles")

# GLOBAL PALETTES
TYPE_COLORS <- c(SNP = "#4E79A7", DEL = "#E15759", DNP = "#F28E2B", TNP = "#59A14F")
CLASS_COLORS <- c(IGR = "#4E79A7", INTRON = "#F28E2B", THREE_PRIME_UTR = "#E15759", RNA = "#76B7B2", FIVE_PRIME_FLANK = "#59A14F")
SNV_COLORS <- c("C>A" = "#00BCD4", "C>G" = "#212121", "C>T" = "#F44336", "T>A" = "#B0BEC5", "T>C" = "#4CAF50", "T>G" = "#FF80AB")
arm_cols <- c("p-arm" = "#E07B54", "q-arm" = "#4472C4")

#For point Mutation
chr_order <- c("chr1","chr3","chr4","chr5","chr6","chr7","chr10","chr15","chr16","chr19","chr20","chr21","chr22")

# publication theme
theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(colour = "black", linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.4),
    axis.title         = element_text(face = "bold", size = 10),
    axis.text          = element_text(size = 8.5, colour = "black"),
    legend.title       = element_text(face = "bold", size = 9),
    legend.text        = element_text(size = 8),
    plot.title         = element_blank(), # Enforce NO titles inside panels
    plot.subtitle      = element_blank(),
    plot.caption       = element_blank()
  )

# ── 2. DATA PROCESSING ────────────────────────────────────────────────────────
# Mutation Burden Data
df <- read.delim("ALL_SAMPLES_PASS_final.tsv", stringsAsFactors = FALSE)
df$Sample <- as.character(df$Sample)
df$CHROM  <- factor(df$CHROM, levels = chr_order)
df$Variant_Type  <- factor(df$Variant_Type,  levels = names(TYPE_COLORS))
df$Variant_Class <- factor(df$Variant_Class, levels = names(CLASS_COLORS))
df$Sample <- factor(df$Sample, levels = as.character(sort(as.integer(unique(df$Sample)))))

# CNV Raw Data
df_raw <- read_tsv("Copy Number Variants.tsv", show_col_types = FALSE)
df_cnv <- df_raw %>%
  mutate(
    cnv_size_mb = (cnv_end - cnv_start) / 1e6,
    Patient     = str_extract(Samples, "(?<=Tumor-)P-\\d+"),
    chr_num     = factor(str_replace(cnv_chr, "chr", ""), levels = c(as.character(1:22), "X","M"), ordered = TRUE)
  )

df_cb <- df_cnv %>%
  mutate(band_list = str_split(Cytoband, ",\\s*")) %>%
  unnest(band_list) %>%
  rename(band = band_list) %>%
  mutate(
    band        = str_trim(band),
    chr_band    = paste0(cnv_chr, " ", band),
    arm         = if_else(str_starts(band, "p"), "p-arm", "q-arm"),
    chr_label_f = factor(str_remove(cnv_chr, "chr"), levels = c(as.character(1:22), "X","M"), ordered = TRUE)
  )

# ── 3. GENERATE INDIVIDUAL PLOTS (p1 - p4, cp1 - cp3b) ────────────────────────

# PLOT 1 – Chromosome Burden
p1 <- df %>%
  count(CHROM, Variant_Type) %>%
  ggplot(aes(x = CHROM, y = n, fill = Variant_Type)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5), size = 2.2, colour = "white", fontface = "bold") +
  scale_x_discrete(labels = function(x) sub("chr", "", x)) +
  scale_fill_manual(values = TYPE_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", y = "Mutation Count", fill = "Variant Type") +
  theme_pub + theme(legend.position = "top")

# PLOT 2 – Per-Sample Mutation Burden
sample_totals <- df %>% count(Sample, name = "total")
p2 <- df %>%
  count(Sample, Variant_Type) %>%
  left_join(sample_totals, by = "Sample") %>%
  ggplot(aes(x = Sample, y = n, fill = Variant_Type)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5), size = 2.2, colour = "white", fontface = "bold") +
  geom_text(data = sample_totals, aes(x = Sample, y = total, label = total, fill = NULL), vjust = -0.4, size = 2.5, fontface = "bold", colour = "#333333") +
  scale_x_discrete(labels = function(x) paste0("P", x)) +
  scale_fill_manual(values = TYPE_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Sample ID", y = "Mutation Count") +
  theme_pub + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# PLOT 3 – SNV Mutation Spectrum
complement <- c(A = "T", T = "A", G = "C", C = "G")
canonical_sub <- function(ref, alt) { if (ref %in% c("C", "T")) paste0(ref, ">", alt) else paste0(complement[ref], ">", complement[alt]) }
sub_order <- c("C>A","C>G","C>T","T>A","T>C","T>G")

snp_counts <- df %>%
  filter(Variant_Type == "SNP") %>%
  rowwise() %>% mutate(Sub = canonical_sub(REF, ALT)) %>% ungroup() %>%
  filter(!is.na(Sub)) %>% count(Sub) %>%
  complete(Sub = sub_order, fill = list(n = 0)) %>% mutate(Sub = factor(Sub, levels = sub_order))

p3 <- snp_counts %>%
  ggplot(aes(x = Sub, y = n, fill = Sub)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n), vjust = -0.4, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = SNV_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Substitution Class", y = "SNV Count") +
  theme_pub + theme(legend.position = "none")

# PLOT 4 – Distribution of Variant Class
p4 <- df %>%
  count(Variant_Class, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1), combined_label = paste0(" ", Variant_Class, " [", n, " (", pct, "%)]"), Variant_Class = fct_reorder(Variant_Class, n)) %>%
  ggplot(aes(x = n, y = Variant_Class, fill = Variant_Class)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = combined_label), x = 0, hjust = 0, size = 2.4, fontface = "bold", colour = "#222222") +
  scale_fill_manual(values = CLASS_COLORS) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Variant Count", y = NULL) +
  theme_pub + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

# PLOT 5 (cp1) — CNV Frequency per Chromosome
chr_freq <- df_cnv %>% count(chr_num, name = "n_CNVs") %>% arrange(chr_num)
cp1 <- chr_freq %>%
  ggplot(aes(x = chr_num, y = n_CNVs)) +
  geom_col(colour = "white", fill = "#4472C4", width = 0.7) +
  geom_text(aes(label = n_CNVs), vjust = -0.4, size = 2.5, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks()) +
  labs(x = "Chromosome", y = "CNV Segment Count") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# PLOT 6 (cp2) — CNV Count per Chromosome Band
cyto_counts <- df_cb %>% count(chr_label_f, band, arm, name = "n_cnvs") %>% filter(n_cnvs > 3) %>% arrange(chr_label_f, arm, band) %>% mutate(chr_band_label = paste0("chr", chr_label_f, " ", band))
cyto_order <- unique(cyto_counts$chr_band_label)

cp2 <- cyto_counts %>%
  mutate(chr_band_label = factor(chr_band_label, levels = cyto_order)) %>%
  ggplot(aes(x = chr_band_label, y = n_cnvs, fill = arm)) +
  geom_col(colour = "white", width = 0.72) +
  geom_text(aes(label = n_cnvs), vjust = -0.45, size = 2.4, fontface = "bold", colour = "black") +
  geom_vline(xintercept = {
    cyto_counts %>% mutate(chr_band_label = factor(chr_band_label, levels = cyto_order)) %>%
      group_by(chr_label_f) %>% summarise(last_pos = max(as.numeric(chr_band_label)), .groups = "drop") %>%
      filter(chr_label_f != last(chr_label_f)) %>% pull(last_pos) + 0.5
  }, colour = "grey70", linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = arm_cols, name = "Arm") +
  scale_y_continuous(breaks = pretty_breaks(n = 4), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Chromosome Band", y = "CNV Segment Count") +
  theme_pub + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))

# PLOT 7 & 8 (cp3a / cp3b) — Patient CNV Summary Profiles
patient_burden <- df_cnv %>% group_by(Patient) %>% summarise(n_cnvs = n(), total_mb = sum(cnv_size_mb), .groups = "drop")

cp3a <- patient_burden %>%
  mutate(Patient = fct_reorder(Patient, total_mb)) %>%
  ggplot(aes(x = Patient, y = total_mb)) +
  geom_col(colour = "white", width = 0.65, fill = "#F28E2B") +
  coord_flip() + labs(x = "Patient ID", y = "Total Altered Region (Mb)") +
  theme_pub + theme(axis.text.y = element_text(size = 7))

cp3b <- patient_burden %>%
  mutate(Patient = fct_reorder(Patient, n_cnvs)) %>%
  ggplot(aes(x = Patient, y = n_cnvs)) +
  geom_col(colour = "white", width = 0.65, fill = "#F28E2B") +
  coord_flip() + labs(x = "Patient ID", y = "CNV Segment Count") +
  theme_pub + theme(axis.text.y = element_text(size = 7))


# ── 4. STACKED STRUCTURAL MATRIX (PATCHWORK ASSEMBLY) ─────────────────────────

row1 <- p1 + p2   + plot_layout(widths = c(1.2, 1)) 
row2 <- p4 + p3   + plot_layout(widths = c(1, 1.1))
row3 <- cp1 + cp2 + plot_layout(widths = c(1, 2.3)) # More width given to complex bands
row4 <- cp3a + cp3b + plot_layout(widths = c(1, 1))

#Subset plots - SNV

final_publication_plot <- ((
  row1 / 
    row2
) + 
  plot_annotation(tag_levels = "A"))& 
  theme(
    plot_tag   = element_text(face = "bold", size = 12, family = "Arial"),
    axis.title = element_text(size = 10, face = "bold", family = "Arial"),
    axis.text  = element_text(size = 8, family = "Arial"),
  )

ggsave(
  filename = "SNV_figure.tiff", 
  plot     = final_publication_plot, 
  width    = 10.0,     
  height   = 10.0,    
  device = "tiff",
  dpi      = 300, 
  bg       = "white",
  compression = "lzw"
)


#Subset plots - CNV

final_publication_plot <- ((
  row3 / 
    row4
) + 
  plot_annotation(tag_levels = "A"))& 
  theme(
    plot_tag   = element_text(face = "bold", size = 12, family = "Arial"),
    axis.title = element_text(size = 10, face = "bold", family = "Arial"),
    axis.text  = element_text(size = 8, family = "Arial"),
  )

ggsave(
  filename = "CNV_figure_vertical.tiff", 
  plot     = final_publication_plot, 
  width    = 10.0,     
  height   = 10.0,    
  device = "tiff",
  dpi      = 300, 
  bg       = "white",
  compression = "lzw"
)

# GENE WORD CLOUD FOR COPY NUMBER GAINS

library(tidyverse)
library(ggwordcloud)


cnv_file <- "Copy Number Variants.tsv" # Update path if needed
df_raw <- read_tsv(cnv_file, show_col_types = FALSE)

# (Data processing steps remain identical to your previous script)
gene_counts <- df_raw %>%
  select(Samples, Genes) %>%
  filter(!is.na(Genes) & Genes != "" & Genes != "-") %>%
  mutate(Gene_Split = str_split(Genes, ";\\s*")) %>%
  unnest(Gene_Split) %>%
  mutate(Gene_Clean = str_trim(Gene_Split)) %>%
  count(Gene_Clean, sort = TRUE, name = "Frequency") %>%
  slice_max(Frequency, n = 52, with_ties = FALSE)

# ── 3. Plot with Zero Margins ─────────────────────────────────────────────────
set.seed(42) 

wordcloud_plot <- ggplot(gene_counts, 
                         aes(label = Gene_Clean, 
                             size = Frequency, 
                             color = Frequency)) +
  # area_scale optimization forces the words to pack much more tightly
  geom_text_wordcloud_area(
    rm_outside = TRUE, 
    shape = "circle",
    eccentricity = 1.0  # Keeps the aspect ratio close to 1:1 circular packing
  ) +
  scale_size_area(max_size = 14) +
  scale_color_viridis_c(option = "plasma", end = 0.85) +
  # theme_void() strips all axes, gridlines, AND default padding margins
  theme_void() + 
  theme(
    plot.margin      = margin(0, 0, 0, 0, "pt"), # Complete zero-out of canvas padding
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# ── 4. Save with tight dimensions ─────────────────────────────────────────────
# Adjusting to a strict 5x5 square layout forces the canvas bounding box 
# to match the native circular shape of the cloud.
ggsave(
  filename = "gene_gains_wordcloud_tight.png",
  plot     = wordcloud_plot,
  width    = 5, 
  height   = 5,
  dpi      = 600,
  bg       = "white"
)

