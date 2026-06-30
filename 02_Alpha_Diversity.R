library(phyloseq)
library(ggplot2)
library(dplyr)
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

# ALPHA DIVERSITY ANALYSIS
# Shannon captures richness + evenness

# Step 1: Calculate Shannon diversity values for each sample
shannon_df <- estimate_richness(ps, measures = "Shannon")
head(shannon_df)

# Step 2: Add grouping information (Mother vs Daughter)
shannon_df$Group <- sample_data(ps)$Group
head(shannon_df)

# Step 3: Summary statistics for reporting
aggregate(Shannon ~ Group, data = shannon_df, mean)
aggregate(Shannon ~ Group, data = shannon_df, median)

# Step 4: Wilcoxon test - non-parametric, recommended for microbiome data
wilcox.test(Shannon ~ Group, data = shannon_df)

# Boxplot with significance
ggplot(shannon_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("Mother" = "#6C8CD5", "Daughter" = "#F4A6A6")) +
  stat_compare_means(
    method      = "wilcox.test",
    label       = "p.format",
    comparisons = list(c("Mother", "Daughter"))
  ) +
  labs(
    title = "Shannon Diversity by Group (Non-Rarefied)",
    x     = "Group",
    y     = "Shannon Index"
  )

# PAIRED vs UNPAIRED FRAMEWORK
# Question: are true mother-daughter pairs more similar to each other
# than random unrelated mother-daughter combinations?
meta_nr <- data.frame(sample_data(ps))
meta_nr$SampleID <- rownames(meta_nr)

alpha_div_nr <- estimate_richness(ps, measures = "Shannon")
alpha_div_nr$SampleID <- rownames(alpha_div_nr)

df_nr <- left_join(
  alpha_div_nr,
  meta_nr[, c("SampleID", "Group", "familyid")],
  by = "SampleID"
)

mothers_nr <- df_nr %>%
  filter(Group == "Mother") %>%
  select(familyid, Shannon_M = Shannon)
daughters_nr <- df_nr %>%
  filter(Group == "Daughter") %>%
  select(familyid, Shannon_D = Shannon)

paired_families_nr <- inner_join(mothers_nr, daughters_nr, by = "familyid")

# Paired - true mother-daughter pairs
shannon_direct_nr <- paired_families_nr %>%
  mutate(
    Shannon_Diff = abs(Shannon_M - Shannon_D),
    Comparison = "Paired"
  )

# Unpaired - mother crossed with every unrelated daughter (null model)
cross_nr <- expand.grid(
  FamilyID_M = paired_families_nr$familyid,
  FamilyID_D = paired_families_nr$familyid,
  stringsAsFactors = FALSE
) %>%
  filter(FamilyID_M != FamilyID_D) %>%
  left_join(paired_families_nr %>% select(familyid, Shannon_M), by = c("FamilyID_M" = "familyid")) %>%
  left_join(paired_families_nr %>% select(familyid, Shannon_D), by = c("FamilyID_D" = "familyid"))

shannon_partial_nr <- cross_nr %>%
  mutate(
    Shannon_Diff = abs(Shannon_M - Shannon_D),
    Comparison = "Unpaired"
  )

shannon_diff_nr <- bind_rows(
  shannon_direct_nr %>% select(Shannon_Diff, Comparison),
  shannon_partial_nr %>% select(Shannon_Diff, Comparison)
)
shannon_diff_nr$Comparison <- factor(shannon_diff_nr$Comparison, levels = c("Paired", "Unpaired"))

wilcox_shannon_nr <- wilcox.test(
  shannon_diff_nr$Shannon_Diff[shannon_diff_nr$Comparison == "Paired"],
  shannon_diff_nr$Shannon_Diff[shannon_diff_nr$Comparison == "Unpaired"]
)
print(wilcox_shannon_nr)

ggplot(shannon_diff_nr, aes(x = Comparison, y = Shannon_Diff, fill = Comparison)) +
  geom_jitter(aes(color = Comparison), width = 0.15, size = 2, alpha = 0.5) +
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
    title = "Shannon Difference - Paired vs Unpaired (Non-Rarefied)",
    x     = NULL,
    y     = "|Shannon Difference|"
  )

# SIMPSON DIVERSITY
# Simpson = 1 - dominance. Values closer to 1 = more diverse,
# closer to 0 = one species dominates. Useful here since the
# vaginal microbiome is often dominated by a single Lactobacillus species.

simpson_df <- estimate_richness(ps, measures = "Simpson")
simpson_df$Group    <- sample_data(ps)$Group
simpson_df$familyid <- sample_data(ps)$familyid
simpson_df$AGE      <- sample_data(ps)$AGE
simpson_df$SampleID <- rownames(simpson_df)
head(simpson_df)

aggregate(Simpson ~ Group, data = simpson_df, mean)
aggregate(Simpson ~ Group, data = simpson_df, median)

wilcox_simpson <- wilcox.test(Simpson ~ Group, data = simpson_df)
print(wilcox_simpson)

ggplot(simpson_df, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(aes(color = AGE), width = 0.2, size = 3) +
  scale_fill_manual(values = c("Daughter" = "#F4A6A6", "Mother" = "#6C8CD5")) +
  scale_color_gradient(low = "#FFD166", high = "#D62828", na.value = "grey50") +
  stat_compare_means(
    method = "wilcox.test",
    label  = "p.format",
    comparisons = list(c("Mother", "Daughter"))
  ) +
  labs(
    title = "Simpson Diversity by Group (Non-Rarefied)",
    x     = "Group",
    y     = "Simpson Index",
    color = "Age"
  )

# Same paired vs unpaired framework as Shannon above
mothers_s <- simpson_df %>%
  filter(Group == "Mother") %>%
  select(familyid, Simpson_M = Simpson)
daughters_s <- simpson_df %>%
  filter(Group == "Daughter") %>%
  select(familyid, Simpson_D = Simpson)

paired_simpson <- inner_join(mothers_s, daughters_s, by = "familyid") %>%
  mutate(
    Simpson_Diff = abs(Simpson_M - Simpson_D),
    Comparison   = "Paired"
  )

cross_s <- expand.grid(
  FamilyID_M = paired_simpson$familyid,
  FamilyID_D = paired_simpson$familyid,
  stringsAsFactors = FALSE
) %>%
  filter(FamilyID_M != FamilyID_D) %>%
  left_join(paired_simpson %>% select(familyid, Simpson_M), by = c("FamilyID_M" = "familyid")) %>%
  left_join(paired_simpson %>% select(familyid, Simpson_D), by = c("FamilyID_D" = "familyid")) %>%
  mutate(
    Simpson_Diff = abs(Simpson_M - Simpson_D),
    Comparison   = "Unpaired"
  )

simpson_diff <- bind_rows(
  paired_simpson %>% select(Simpson_Diff, Comparison),
  cross_s        %>% select(Simpson_Diff, Comparison)
)
simpson_diff$Comparison <- factor(simpson_diff$Comparison, levels = c("Paired", "Unpaired"))

wilcox_simpson_paired <- wilcox.test(
  simpson_diff$Simpson_Diff[simpson_diff$Comparison == "Paired"],
  simpson_diff$Simpson_Diff[simpson_diff$Comparison == "Unpaired"]
)
print(wilcox_simpson_paired)

ggplot(simpson_diff, aes(x = Comparison, y = Simpson_Diff, fill = Comparison)) +
  geom_jitter(aes(color = Comparison), width = 0.15, size = 2, alpha = 0.5) +
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
    title = "Simpson Difference - Paired vs Unpaired (Non-Rarefied)",
    x     = NULL,
    y     = "|Simpson Difference|"
  )

# Age correlations with Simpson
cor.test(simpson_df$AGE, simpson_df$Simpson, method = "spearman")

cor.test(simpson_df$AGE[simpson_df$Group == "Mother"],
         simpson_df$Simpson[simpson_df$Group == "Mother"],
         method = "spearman")

cor.test(simpson_df$AGE[simpson_df$Group == "Daughter"],
         simpson_df$Simpson[simpson_df$Group == "Daughter"],
         method = "spearman")

ggplot(simpson_df, aes(x = AGE, y = Simpson, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Mother" = "#6C8CD5", "Daughter" = "#F4A6A6")) +
  labs(
    title = "Age vs Simpson Diversity by Group",
    x     = "Age",
    y     = "Simpson Index"
  )

