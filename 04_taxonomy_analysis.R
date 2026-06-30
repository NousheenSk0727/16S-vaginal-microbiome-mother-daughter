library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

theme_set(theme_bw())

# Load saved objects from pipeline script
# Run 01_dada2_pipeline.R first, or load saved objects below
seqtab.nochim <- readRDS("/Users/nousheenjahanshaik/Desktop/seqtab.nochim.rds")
samdf         <- readRDS("/Users/nousheenjahanshaik/Desktop/srr.rds")
taxa          <- read.csv("~/Desktop/taxonomy_assignments.csv", row.names = 1)

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(as.matrix(taxa))
)
# TAXONOMIC COMPOSITION ANALYSIS

# Glom to Genus level - combines all ASVs belonging to the same genus
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus

ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Top 10 genera overall
top10_names <- names(sort(taxa_sums(ps_genus), decreasing = TRUE))[1:10]
ps_top10    <- prune_taxa(top10_names, ps_genus_rel)

top10_df <- psmelt(ps_top10)
top10_df$Genus[is.na(top10_df$Genus)] <- "Unknown"

ggplot(top10_df, aes(x = reorder(Genus, -Abundance), y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "none"
  ) +
  labs(
    title = "Top 10 Genera - Overall Relative Abundance",
    x     = "Genus",
    y     = "Mean Relative Abundance"
  )

# Top 10 genera by Group (Mother vs Daughter)
top10_group <- top10_df %>%
  group_by(Group, Genus) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop")

ggplot(top10_group, aes(x = reorder(Genus, -Mean_Abundance), y = Mean_Abundance, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Mother" = "#6C8CD5", "Daughter" = "#F4A6A6")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  labs(
    title = "Top 10 Genera - Mother vs Daughter",
    x     = "Genus",
    y     = "Mean Relative Abundance",
    fill  = "Group"
  )

# Community composition per sample
top10_df$Sample <- factor(top10_df$Sample)

ggplot(top10_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.text     = element_text(face = "italic"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Community Composition per Sample",
    x     = "Samples",
    y     = "Relative Abundance",
    fill  = "Genus"
  )

# LACTOBACILLUS DOMINANCE
# which() instead of boolean indexing on the Genus column - the direct
# comparison was returning NAs for taxa with missing genus assignments
genus_col  <- as.character(tax_table(ps_genus_rel)[, "Genus"])
lacto_idx  <- which(genus_col == "Lactobacillus")
lacto_taxa <- taxa_names(ps_genus_rel)[lacto_idx]

lacto_abund <- rowSums(otu_table(ps_genus_rel)[, lacto_taxa, drop = FALSE])

lacto_df <- data.frame(
  SampleID      = sample_names(ps_genus_rel),
  Lactobacillus = as.numeric(lacto_abund),
  Group         = sample_data(ps_genus_rel)$Group,
  AGE           = sample_data(ps_genus_rel)$AGE
)

head(lacto_df)
summary(lacto_df$Lactobacillus)

wilcox.test(Lactobacillus ~ Group, data = lacto_df)

ggplot(lacto_df, aes(x = Group, y = Lactobacillus, fill = Group)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(aes(color = AGE), width = 0.2, size = 3) +
  scale_fill_manual(values  = c("Mother" = "#6C8CD5", "Daughter" = "#F4A6A6")) +
  scale_color_gradient(low = "#FFD166", high = "#D62828", na.value = "grey50") +
  stat_compare_means(
    method      = "wilcox.test",
    label       = "p.format",
    comparisons = list(c("Mother", "Daughter"))
  ) +
  labs(
    title = "Lactobacillus Relative Abundance by Group",
    x     = "Group",
    y     = "Relative Abundance",
    color = "Age"
  )

# BV-ASSOCIATED GENERA
# Gardnerella, Atopobium, Sneathia
bv_genera <- c("Gardnerella", "Atopobium", "Sneathia")

bv_df <- data.frame(SampleID = sample_names(ps_genus_rel))
for (g in bv_genera) {
  g_idx  <- which(genus_col == g)
  g_taxa <- taxa_names(ps_genus_rel)[g_idx]
  bv_df[[g]] <- rowSums(otu_table(ps_genus_rel)[, g_taxa, drop = FALSE])
}
bv_df$Group    <- sample_data(ps_genus_rel)$Group
bv_df$AGE      <- sample_data(ps_genus_rel)$AGE
bv_df$BV_Total <- rowSums(bv_df[, bv_genera], na.rm = TRUE)

wilcox.test(Gardnerella ~ Group, data = bv_df)
wilcox.test(Atopobium   ~ Group, data = bv_df)
wilcox.test(Sneathia    ~ Group, data = bv_df)
wilcox.test(BV_Total    ~ Group, data = bv_df)

bv_long <- bv_df %>%
  pivot_longer(cols = all_of(bv_genera), names_to = "Genus", values_to = "Abundance")

ggplot(bv_long, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_fill_manual(values = c("Mother" = "#6C8CD5", "Daughter" = "#F4A6A6")) +
  stat_compare_means(
    method      = "wilcox.test",
    label       = "p.format",
    comparisons = list(c("Mother", "Daughter"))
  ) +
  theme(strip.text = element_text(face = "italic"), legend.position = "none") +
  labs(
    title = "BV-Associated Genera by Group",
    x     = "Group",
    y     = "Relative Abundance"
  )

# DOMINANT GENUS PER SAMPLE
genus_mat   <- as.data.frame(otu_table(ps_genus_rel))
genus_names <- tax_table(ps_genus_rel)[, "Genus"]

dominant_genus <- apply(genus_mat, 1, function(x) as.character(genus_names[which.max(x)]))

dominant_df <- data.frame(
  SampleID       = names(dominant_genus),
  Dominant_Genus = dominant_genus,
  Group          = sample_data(ps_genus_rel)$Group,
  familyid       = sample_data(ps_genus_rel)$familyid
)

table(dominant_df$Group, dominant_df$Dominant_Genus)

ggplot(dominant_df, aes(x = Group, fill = Dominant_Genus)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.text = element_text(face = "italic")) +
  labs(
    title = "Dominant Genus per Sample by Group",
    x     = "Group",
    y     = "Proportion of Samples",
    fill  = "Dominant Genus"
  )

# Do mother-daughter pairs share the same dominant genus?
mothers_dom   <- dominant_df %>%
  filter(Group == "Mother") %>%
  select(familyid, Dominant_M = Dominant_Genus)
daughters_dom <- dominant_df %>%
  filter(Group == "Daughter") %>%
  select(familyid, Dominant_D = Dominant_Genus)

paired_dominant <- inner_join(mothers_dom, daughters_dom, by = "familyid") %>%
  mutate(Same_Dominant = Dominant_M == Dominant_D)

table(paired_dominant$Same_Dominant)
prop.table(table(paired_dominant$Same_Dominant))
table(paired_dominant$Dominant_M, paired_dominant$Dominant_D)

paired_dominant_long <- paired_dominant %>%
  pivot_longer(
    cols      = c(Dominant_M, Dominant_D),
    names_to  = "Role",
    values_to = "Dominant_Genus"
  ) %>%
  mutate(Role = ifelse(Role == "Dominant_M", "Mother", "Daughter"))

ggplot(paired_dominant_long, aes(x = Role, fill = Dominant_Genus)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.text = element_text(face = "italic")) +
  labs(
    title = "Dominant Genus in Paired Mother-Daughter Samples",
    x     = NULL,
    y     = "Proportion",
    fill  = "Dominant Genus"
  )

# HEATMAP - top 10 genera across paired families
paired_samples <- inner_join(
  data.frame(sample_data(ps_genus_rel)) %>%
    filter(Group == "Mother") %>%
    select(SampleID_M = Run, familyid),
  data.frame(sample_data(ps_genus_rel)) %>%
    filter(Group == "Daughter") %>%
    select(SampleID_D = Run, familyid),
  by = "familyid"
)

all_paired_samples <- c(paired_samples$SampleID_M, paired_samples$SampleID_D)
ps_paired <- prune_samples(all_paired_samples, ps_genus_rel)

top10_paired <- prune_taxa(top10_names, ps_paired)
heatmap_df   <- psmelt(top10_paired)
heatmap_df$Genus[is.na(heatmap_df$Genus)] <- "Unknown"

ggplot(heatmap_df, aes(x = Sample, y = Genus, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#D62828") +
  facet_wrap(~ Group, scales = "free_x") +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(face = "italic")
  ) +
  labs(
    title = "Top 10 Genera Heatmap - Paired Samples",
    x     = "Samples",
    y     = "Genus",
    fill  = "Relative\nAbundance"
  )
