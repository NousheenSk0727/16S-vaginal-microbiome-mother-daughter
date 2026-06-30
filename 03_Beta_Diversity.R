library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)

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
)# BETA DIVERSITY ANALYSIS
# Jaccard = presence/absence (which species are shared)
# Bray-Curtis = abundance-weighted (how much of each species is shared)

# JACCARD

# Jaccard is usually calculated on presence/absence data,
# so first convert counts to 1 (present) / 0 (absent)
ps.pa <- transform_sample_counts(ps, function(x) ifelse(x > 0, 1, 0))

ord.nmds.jaccard <- ordinate(ps.pa, method = "NMDS", distance = "jaccard")

plot_ordination(ps.pa, ord.nmds.jaccard, color = "Group") +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Mother" = "#6C8CD5", "Daughter" = "#F4A6A6")) +
  stat_ellipse(type = "t", linetype = 2) +
  labs(title = "Jaccard NMDS - Mother vs Daughter")

# Same paired vs unpaired framework used for alpha diversity
jaccard_dist_nr <- distance(ps, method = "jaccard", binary = TRUE)
jaccard_mat_nr  <- 1 - as.matrix(jaccard_dist_nr)  # convert distance to similarity

meta_nr <- data.frame(sample_data(ps))
meta_nr$SampleID <- rownames(meta_nr)

mothers_nr   <- meta_nr %>% filter(Group == "Mother")
daughters_nr <- meta_nr %>% filter(Group == "Daughter")

# Paired - true mother-daughter pairs
direct_nr <- inner_join(
  mothers_nr   %>% select(SampleID_M = SampleID, familyid),
  daughters_nr %>% select(SampleID_D = SampleID, familyid),
  by = "familyid"
) %>%
  rowwise() %>%
  mutate(
    Jaccard = jaccard_mat_nr[SampleID_M, SampleID_D],
    Comparison = "Paired"
  ) %>%
  ungroup()

# Unpaired - every cross-family combination
partial_nr <- inner_join(
  mothers_nr   %>% select(SampleID_M = SampleID, FamilyID_M = familyid),
  daughters_nr %>% select(SampleID_D = SampleID, FamilyID_D = familyid),
  by = character()
) %>%
  filter(FamilyID_M != FamilyID_D) %>%
  rowwise() %>%
  mutate(
    Jaccard = jaccard_mat_nr[SampleID_M, SampleID_D],
    Comparison = "Unpaired"
  ) %>%
  ungroup()

plot_df_j_nr <- bind_rows(
  direct_nr  %>% select(Jaccard, Comparison),
  partial_nr %>% select(Jaccard, Comparison)
)
plot_df_j_nr$Comparison <- factor(plot_df_j_nr$Comparison, levels = c("Paired", "Unpaired"))

wilcox_jaccard_nr <- wilcox.test(
  plot_df_j_nr$Jaccard[plot_df_j_nr$Comparison == "Paired"],
  plot_df_j_nr$Jaccard[plot_df_j_nr$Comparison == "Unpaired"]
)
print(wilcox_jaccard_nr)

ggplot(plot_df_j_nr, aes(x = Comparison, y = Jaccard, fill = Comparison)) +
  geom_jitter(aes(color = Comparison), width = 0.15, size = 1.5, alpha = 0.4) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +
  stat_compare_means(
    method      = "wilcox.test",
    label       = "p.format",
    comparisons = list(c("Paired", "Unpaired"))
  ) +
  scale_fill_manual(values  = c("Paired" = "#534AB7", "Unpaired" = "#D85A30")) +
  scale_color_manual(values = c("Paired" = "#534AB7", "Unpaired" = "#D85A30")) +
  theme(legend.position = "none") +
  labs(
    title = "Jaccard Similarity - Paired vs Unpaired (Non-Rarefied)",
    x     = NULL,
    y     = "Jaccard Similarity"
  )

# BRAY-CURTIS

# Bray-Curtis uses actual abundance counts (not presence/absence)
# so ps is used directly without the binary transformation
bc_dist <- phyloseq::distance(ps, method = "bray")
bc_mat  <- 1 - as.matrix(bc_dist)

meta_bc <- data.frame(sample_data(ps))
meta_bc$SampleID <- rownames(meta_bc)

mothers_bc   <- meta_bc %>% filter(Group == "Mother")
daughters_bc <- meta_bc %>% filter(Group == "Daughter")

# Paired
direct_bc <- inner_join(
  mothers_bc   %>% select(SampleID_M = SampleID, familyid),
  daughters_bc %>% select(SampleID_D = SampleID, familyid),
  by = "familyid"
) %>%
  rowwise() %>%
  mutate(
    BrayCurtis = bc_mat[SampleID_M, SampleID_D],
    Comparison = "Paired"
  ) %>%
  ungroup()

# Unpaired
partial_bc <- inner_join(
  mothers_bc   %>% select(SampleID_M = SampleID, FamilyID_M = familyid),
  daughters_bc %>% select(SampleID_D = SampleID, FamilyID_D = familyid),
  by = character()
) %>%
  filter(FamilyID_M != FamilyID_D) %>%
  rowwise() %>%
  mutate(
    BrayCurtis = bc_mat[SampleID_M, SampleID_D],
    Comparison = "Unpaired"
  ) %>%
  ungroup()

plot_df_bc <- bind_rows(
  direct_bc  %>% select(BrayCurtis, Comparison),
  partial_bc %>% select(BrayCurtis, Comparison)
)
plot_df_bc$Comparison <- factor(plot_df_bc$Comparison, levels = c("Paired", "Unpaired"))

wilcox_bc <- wilcox.test(
  plot_df_bc$BrayCurtis[plot_df_bc$Comparison == "Paired"],
  plot_df_bc$BrayCurtis[plot_df_bc$Comparison == "Unpaired"]
)
print(wilcox_bc)

ggplot(plot_df_bc, aes(x = Comparison, y = BrayCurtis, fill = Comparison)) +
  geom_jitter(aes(color = Comparison), width = 0.15, size = 1.5, alpha = 0.4) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +
  stat_compare_means(
    method      = "wilcox.test",
    label       = "p.format",
    comparisons = list(c("Paired", "Unpaired"))
  ) +
  scale_fill_manual(values  = c("Paired" = "#534AB7", "Unpaired" = "#D85A30")) +
  scale_color_manual(values = c("Paired" = "#534AB7", "Unpaired" = "#D85A30")) +
  theme(legend.position = "none") +
  labs(
    title = "Bray-Curtis Similarity - Paired vs Unpaired (Non-Rarefied)",
    x     = NULL,
    y     = "Bray-Curtis Similarity"
  )

# Density plot - distribution comparison
ggplot(plot_df_bc, aes(x = BrayCurtis, fill = Comparison, color = Comparison)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values  = c("Paired" = "#534AB7", "Unpaired" = "#D85A30")) +
  scale_color_manual(values = c("Paired" = "#534AB7", "Unpaired" = "#D85A30")) +
  theme(legend.position = "top") +
  labs(
    title = "Bray-Curtis Similarity Distribution",
    x     = "Bray-Curtis Similarity",
    y     = "Density"
  )

ord_bc <- ordinate(ps, method = "NMDS", distance = "bray")

plot_ordination(ps, ord_bc, color = "Group") +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Mother" = "#6C8CD5", "Daughter" = "#F4A6A6")) +
  stat_ellipse(type = "t", linetype = 2) +
  labs(title = "Bray-Curtis NMDS - Mother vs Daughter", color = "Group")

# PERMANOVA - tests whether Group and Age explain community variation
# set.seed fixes the random permutation order so the p-value
# is exactly reproducible on re-run
meta_perm_bc <- data.frame(sample_data(ps))
meta_perm_bc_clean <- meta_perm_bc %>% filter(!is.na(AGE))

bc_dist_clean <- phyloseq::distance(
  prune_samples(rownames(meta_perm_bc_clean), ps), method = "bray"
)

set.seed(123)
permanova_bc <- adonis2(bc_dist_clean ~ Group + AGE, data = meta_perm_bc_clean, permutations = 999)
print(permanova_bc)

# Age correlations with Bray-Curtis paired similarity
age_bc <- inner_join(
  direct_bc %>% select(familyid, BrayCurtis),
  data.frame(sample_data(ps)) %>% filter(Group == "Mother") %>% select(familyid, Age_M = AGE),
  by = "familyid"
) %>%
  inner_join(
    data.frame(sample_data(ps)) %>% filter(Group == "Daughter") %>% select(familyid, Age_D = AGE),
    by = "familyid"
  ) %>%
  mutate(Age_Gap = abs(Age_M - Age_D))

cor_agegap_bc <- cor.test(age_bc$Age_Gap, age_bc$BrayCurtis, method = "spearman")
print(cor_agegap_bc)

cor_aged_bc <- cor.test(age_bc$Age_D, age_bc$BrayCurtis, method = "spearman")
print(cor_aged_bc)

ggplot(age_bc, aes(x = Age_Gap, y = BrayCurtis)) +
  geom_point(size = 3, color = "#534AB7", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate(
    "text",
    x     = max(age_bc$Age_Gap, na.rm = TRUE) * 0.7,
    y     = max(age_bc$BrayCurtis, na.rm = TRUE) * 0.95,
    label = paste0("rho = ", round(cor_agegap_bc$estimate, 3),
                   "\np = ", round(cor_agegap_bc$p.value, 3)),
    size = 4
  ) +
  labs(
    title = "Age Gap vs Bray-Curtis Similarity",
    x     = "Age Gap (years)",
    y     = "Bray-Curtis Similarity"
  )

