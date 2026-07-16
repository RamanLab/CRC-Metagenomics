# Importing all necessary libraries
setwd("~/Downloads/NII") #Setting the working directory

library(phyloseq)
library(biomformat)
library(microbiome) 
library(RColorBrewer) 
library(ggpubr) 
library(DT) 
library(data.table) 
library(dplyr) 
library(tidyverse) 
library(viridis)
library(igraph)
library(vegan)
library(Rtsne)
library(qiime2R)
library(venn)
library(Boruta)
library(gapminder)
library(modelr)
library(ggtext)
library(patchwork)
library(randomcoloR)

# GLOBAL FONT / THEME CONSTANTS
BASE_SIZE   <- 10   # base_size passed to theme_classic() for all plots
AXIS_TEXT   <- 8    # size for all axis tick labels
AXIS_TITLE  <- 9    # size for all axis titles
STRIP_TEXT  <- 8    # size for all facet strip labels
LEGEND_TEXT <- 7    # size for legend item labels
LEGEND_TTL  <- 8    # size for legend title

# DATA LOADING AND PHYLOSEQ OBJECT CREATION
otu_data  <- read.delim2("NII_Abundance_custom_point1.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
taxa_data <- read.delim2("NII_Taxa_custom_point1.tsv", sep = "\t", header = TRUE, row.names = 1) %>% as.matrix()
meta_data <- read.delim2("NII_Metadata.tsv", sep = "\t", header = TRUE, row.names = 1)

my_otu  <- otu_table(otu_data, taxa_are_rows = TRUE)
my_taxa <- tax_table(taxa_data)
my_meta <- sample_data(meta_data)

ps_NII  <- phyloseq(my_otu, my_taxa, my_meta)
PS_NII  <- subset_taxa(ps_NII, Domain == "Bacteria" | Domain == "Archaea")

# Transform to relative abundance
PS_NII_relabundance <- transform_sample_counts(PS_NII, function(x) x / sum(x))

set.seed(123)  
custom_colors <- distinctColorPalette(50)  

# ANALYSIS 1: BORUTA & WILCOXON (TUMOR VS TUMOR-ADJACENT)
my_taxa_names <- as.data.frame(tax_table(PS_NII_relabundance))
taxa_names(PS_NII_relabundance) <- my_taxa_names$Species

pseq_df_type <- as.data.frame(t(phyloseq::otu_table(PS_NII_relabundance)))
pseq_df_type$type <- phyloseq::sample_data(PS_NII_relabundance)$type %>% as.factor()

set.seed(12345)
sig_taxa_type <- Boruta(type ~ ., data = pseq_df_type, pValue = 0.05, mcAdj = TRUE, maxRuns = 999)
boruta_fixed_type <- TentativeRoughFix(sig_taxa_type)
selected_biomarkers <- getSelectedAttributes(boruta_fixed_type)

# Wilcoxon Modeling
wilcox_model <- function(df) { wilcox.test(abund ~ type, data = df) }
wilcox_pval  <- function(df) { wilcox.test(abund ~ type, data = df)$p.value }

wilcox_results <- pseq_df_type %>%
  gather(key = OTU, value = abund, -type) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value     = map(data, wilcox_pval))  

wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest(cols = c(p_value))

taxa_info <- data.frame(tax_table(PS_NII_relabundance)) %>% rownames_to_column(var = "OTU")

wilcox_results <- wilcox_results %>%
  full_join(taxa_info, by = "OTU") %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())

# PLOT A: BIOMARKERS BOXPLOT
target_physeq_type <- prune_taxa(selected_biomarkers, PS_NII_relabundance)
richness_melted_type <- phyloseq::psmelt(target_physeq_type)

richness_melted_type <- richness_melted_type %>%
  mutate(OTU_wrap = stringr::str_wrap(OTU, width = 12)) %>%
  mutate(type = factor(type, levels = c("Tumor-Adjacent", "Tumor")))

p1_biomarkers <- ggplot(richness_melted_type, aes(x = type, y = Abundance, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, color = "black", width = 0.55) +
  geom_jitter(color = "#1a1a1a", width = 0.1, size = 1.2, alpha = 0.8) +
  facet_wrap(~ OTU_wrap, scales = "free_y") +
  stat_compare_means(
    method      = "wilcox.test",
    comparisons = list(c("Tumor-Adjacent", "Tumor")),
    label       = "p.format",
    tip.length  = 0.01,
    size        = 3
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("Tumor" = "#1a3a6b", "Tumor-Adjacent" = "#a8c8e8")) +
  labs(x = NULL, y = "Relative Abundance") +
  theme_classic(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    #  Axis labels & titles 
    axis.text.x      = element_text(size = AXIS_TEXT, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
    axis.text.y      = element_text(size = AXIS_TEXT, color = "black"),
    axis.title.y     = element_text(size = AXIS_TITLE, face = "bold"),
    #  Strip text
    strip.text       = element_text(size = STRIP_TEXT, face = "italic", hjust = 0.5, lineheight = 0.9),
    strip.background = element_blank(),
    legend.position  = "none",
    plot.background  = element_blank(),
    plot.margin      = margin(6, 6, 6, 6)
  )

# ANALYSIS 2: BORUTA - DUMBBELL (LEFT VS RIGHT)
pseq_df_side <- as.data.frame(t(phyloseq::otu_table(PS_NII_relabundance)))
pseq_df_side$Side <- phyloseq::sample_data(PS_NII_relabundance)$Side %>% as.factor()

set.seed(12345)
sig_taxa_side <- Boruta(Side ~ ., data = pseq_df_side, pValue = 0.05, mcAdj = TRUE, maxRuns = 999)
boruta_fixed_side <- TentativeRoughFix(sig_taxa_side)
selected_sides <- getSelectedAttributes(boruta_fixed_side)

target_physeq_side <- prune_taxa(selected_sides, PS_NII_relabundance)
richness_melted_side <- phyloseq::psmelt(target_physeq_side)
richness_melted_side$Side <- factor(richness_melted_side$Side, levels = c("Left", "Right"))

dumbbell_df <- richness_melted_side %>%
  group_by(OTU, Side) %>%
  summarise(mean_ab = mean(Abundance), .groups = "drop")

pvals <- richness_melted_side %>%
  group_by(OTU) %>%
  summarise(p = wilcox.test(Abundance ~ Side)$p.value,
            sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
            .groups = "drop")

dumbbell_wide <- dumbbell_df %>%
  tidyr::pivot_wider(names_from = Side, values_from = mean_ab) %>%
  left_join(pvals, by = "OTU") %>%
  mutate(max_val = pmax(Left, Right),
         OTU = factor(OTU, levels = rev(unique(OTU))))

# PLOT B: DUMBBELL PLOT 
p_dumbbell <- ggplot(dumbbell_wide) +
  geom_segment(aes(x = Right, xend = Left, y = OTU, yend = OTU), color = "grey70", linewidth = 0.6) +
  geom_point(aes(x = Right, y = OTU), color = "#a8c8e8", size = 3, shape = 19) +
  geom_point(aes(x = Left,  y = OTU), color = "#1a3a6b", size = 3, shape = 19) +
  geom_text(aes(x = max_val * 1.5, y = OTU, label = sig), size = 2.8, hjust = 0, color = "black") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.05))  
  )+
  labs(
    x     = "Mean Relative Abundance",
    y     = NULL,
    title = "<span style='color:#1a3a6b'>&#9679; Left</span> &nbsp; <span style='color:#a8c8e8'>&#9679; Right</span>"
  ) +
  theme_classic(base_size = BASE_SIZE, base_family = "Arial") +          
  theme(
    plot.title         = element_markdown(size = STRIP_TEXT, margin = margin(b = 3)),
    # Y-axis (species names) 
    axis.text.y        = element_text(face = "italic", size = AXIS_TEXT, color = "black"),
    axis.title.x       = element_text(size = AXIS_TITLE, face = "bold", margin = margin(t = 6)),
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.line.x = element_line(color = "black"),   # make sure x axis line is drawn
    axis.ticks.length = unit(0.15, "cm"),           # reduce tick length pushing label down
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    legend.position    = "none",
    plot.background    = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),   # match the tick length C uses implicitly
    plot.margin = margin(t = 8, r = 40, b = 0, l = 60)  # ← b = 0 removes bottom phantom space
  )

p_dumbbell
# ANALYSIS 3: BRAY-CURTIS DISTANCE STRIP PLOT
sample_data(PS_NII) <- sample_data(meta_data)
dist_matrix <- distance(PS_NII, method = "bray")

patient_groups <- table(sample_data(PS_NII)$Patient, sample_data(PS_NII)$type)
patients_with_pairs <- rownames(patient_groups)[rowSums(patient_groups[, c("Tumor", "Tumor-Adjacent")] > 0) == 2]

patient_distances <- data.frame()
for(patient in patients_with_pairs) {
  patient_samples        <- sample_data(PS_NII)$Patient == patient
  tumor_samples          <- patient_samples & sample_data(PS_NII)$type == "Tumor"
  tumor_adjacent_samples <- patient_samples & sample_data(PS_NII)$type == "Tumor-Adjacent"
  
  if(sum(tumor_samples) > 0 & sum(tumor_adjacent_samples) > 0) {
    tumor_dist <- as.matrix(dist_matrix)[tumor_samples, tumor_adjacent_samples]
    pairs_df   <- expand.grid(tumor_sample    = sample_names(PS_NII)[tumor_samples],
                              nontumor_sample = sample_names(PS_NII)[tumor_adjacent_samples],
                              stringsAsFactors = FALSE)
    pairs_df$bray_distance <- as.vector(tumor_dist)
    pairs_df$Patient       <- patient
    patient_distances      <- rbind(patient_distances, pairs_df)
  }
}

n_rows_c     <- 4
row_heights_c <- seq(-0.15, 0.15, length.out = n_rows_c)

patient_distances_fixed <- patient_distances %>%
  arrange(bray_distance) %>%
  mutate(swarm_y = rep(row_heights_c, length.out = n()))

# PLOT C: STRIP PLOT 
p_strip <- ggplot(patient_distances_fixed, aes(x = bray_distance, y = swarm_y)) +
  geom_point(size = 9, color = "#D85A30", fill = "#F0997B", shape = 21, stroke = 1, alpha = 0.9) +
  geom_text(aes(label = Patient), size = 2.2, color = "#4A1B0C", fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(-0.25, 0.25)) +
  labs(x = "Bray-Curtis dissimilarity", y = NULL) +
  theme_classic(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    axis.text.x         = element_text(size = AXIS_TEXT, color = "black"),
    axis.title.x        = element_text(size = AXIS_TITLE, face = "bold", margin = margin(t = 6)),
    axis.text.y         = element_blank(),
    axis.title.y        = element_blank(),
    axis.ticks.y        = element_blank(),
    axis.line.y         = element_blank(),
    panel.grid          = element_blank(),
    panel.grid.major.x  = element_line(color = "grey92", linewidth = 0.4),
    plot.background     = element_blank(),
    plot.margin         = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = NA, color = NA),
    plot.title.position = "panel"
  )
  

# ANALYSIS 4: TOP 15 TAXA COMPOSITION PLOT
ps_NII_rel_abund <- microbiome::transform(PS_NII, "compositional")
top15_taxa <- names(sort(taxa_sums(ps_NII_rel_abund), decreasing = TRUE))[1:15]
top15 <- prune_taxa(top15_taxa, ps_NII_rel_abund)

species_palette_15 <- c(
  "#C10020", "#004949", "#009E73", "#56B4E9", "#006CD1",
  "#4B0092", "#B66DFF", "#B6DBFF", "#E69F00", "#D55E00",
  "#CC79A7", "#FF6DB6", "#F0E442", "#924900", "#888888"
)
names(species_palette_15) <- top15_taxa

#  PLOT D: TOP 15 BARPLOT 
p_top15 <- plot_bar(top15, x = "Patient", fill = "Species") +
  facet_wrap(~ type, scales = "free_x") +
  labs(x = NULL, y = "Relative abundance", fill = "Species") +
  scale_fill_manual(values = species_palette_15) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    axis.text.y      = element_text(size = AXIS_TEXT,  color = "black"),
    axis.text.x      = element_text(size = AXIS_TEXT,  color = "black"),
    axis.title.x     = element_text(size = AXIS_TITLE, face = "bold"),
    axis.title.y     = element_text(size = AXIS_TITLE, face = "bold"),
    strip.text       = element_text(size = STRIP_TEXT, face = "bold", color = "black"),
    strip.background = element_blank(),
    legend.text      = element_text(size = LEGEND_TEXT, face = "italic"),
    legend.title     = element_text(size = LEGEND_TTL,  face = "bold"),
    legend.key.size  = unit(0.35, "cm"),
    legend.position  = "right",
    panel.spacing    = unit(0.8, "lines"),
    plot.background  = element_blank(),
    plot.margin      = margin(4, 4, 4, 4)
  ) +
  coord_flip()

# 
# PATCHWORK ASSEMBLY

# ── Nest A and B first so they align only with each other ────────────────────
AB_row <- (p1_biomarkers | p_dumbbell) +
  plot_layout(widths = c(1, 1)) &   # equal width for A and B
  theme(plot.background = element_rect(fill = "white", color = NA))

# ── Stack AB / C / D with CCCC and DDDD spanning full width ──────────────────
combined_figure <- (AB_row / p_strip / p_top15) +
  plot_layout(
    heights = c(2.5, 0.7, 2)   # top row taller, C slim, D tall
  ) &
  theme(plot.background = element_rect(fill = "white", color = NA))

# ── Annotation: tags, border, outer margin ───────────────────────────────────
combined_figure <- combined_figure +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag        = element_text(size = 12, face = "bold"),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      plot.margin     = margin(12, 12, 12, 12)
    )
  )

print(combined_figure)

ggsave(
  "combined_figure_final.tiff",
  combined_figure,
  device = "tiff",
  width  = 13,
  height = 12,
  units  = "in",
  dpi    = 300,
  compression = "lzw"
)

