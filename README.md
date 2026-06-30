# 16S-vaginal-microbiome-mother-daughter

## Overview
This project investigates whether mother-daughter pairs share more similar vaginal microbiome communities than unrelated women, using publicly available 16S rRNA amplicon sequencing data from 87 samples across 33 paired families (NCBI SRA: PRJNA779415).

## Data
- **Source:** NCBI SRA (PRJNA779415)
- **Samples:** 87 vaginal microbiome samples - 33 mother-daughter pairs
- **Sequencing:** 16S rRNA amplicon sequencing (paired-end, Illumina)
- **Reference database:** SILVA v138.1

## Pipeline
Raw FASTQ → Quality Control → Filtering & Trimming → Error Learning → Denoising → ASV Table → Chimera Removal → Taxonomy Assignment → phyloseq → Diversity & Taxonomy Analysis

## Scripts
| Script | Description |
|--------|-------------|
| `01_dada2_pipeline.R` | Raw FASTQ processing through phyloseq object construction |
| `02_alpha_diversity.R` | Shannon and Simpson diversity - group comparison and paired vs unpaired framework |
| `03_beta_diversity.R` | Jaccard and Bray-Curtis - NMDS, paired vs unpaired, PERMANOVA, age correlations |
| `04_taxonomy_analysis.R` | Genus-level composition, Lactobacillus dominance, BV-associated genera, dominant genus sharing across pairs |

## Key Findings
## Key Findings

Alpha Diversity
- Shannon diversity did not differ significantly between mothers and daughters (p = 0.46)
- Simpson diversity did not differ significantly between groups (p = 0.30)
- No significant difference in Shannon or Simpson when comparing paired vs unrelated combinations
- Alpha diversity is similar between groups - both mothers and daughters carry communities of comparable richness and evenness

Beta Diversity
- True mother-daughter pairs shared significantly more similar microbial communities than unrelated pairs
- Jaccard similarity (presence/absence): paired vs unpaired p < 0.00001
- Bray-Curtis similarity (abundance-weighted): paired vs unpaired p = 0.0006
- PERMANOVA showed Group and Age together explained ~3.8% of community variation (p = 0.054, trending)
- No significant age gap or daughter age correlation with Bray-Curtis similarity

Taxonomy
- *Lactobacillus* and *Gardnerella* identified as dominant genera across all samples
- Daughters showed significantly higher *Lactobacillus* abundance (p = 0.024)
- Mothers carried significantly more BV-associated bacteria - *Atopobium* (p = 0.033), total BV burden (p = 0.038)
- 30% of mother-daughter pairs shared the same dominant genus
- Most common combination: Mother *Gardnerella* dominant + Daughter *Lactobacillus* dominant (13 pairs)

Limitations & Biological Context
- The original research paper on this dataset discusses vertical transfer of vaginal microorganisms from mother to daughter
- Our findings show that mother-daughter pairs are significantly more similar to each other than unrelated pairs in terms of beta diversity - however this alone is not sufficient evidence to confirm vertical transmission
- The observed similarity could reflect shared genetics, shared household environment, similar dietary habits, or other familial factors rather than direct microbial transfer
- Longitudinal data tracking microbiome composition from birth through adolescence would be needed to establish vertical transmission more definitively
- This analysis provides a cross-sectional diversity framework and should be interpreted as exploratory rather than mechanistic

## Tools & Packages
R 4.6.0 · DADA2 · phyloseq · vegan · ggplot2 · ggpubr · dplyr · tidyr · Biostrings

## Author
Nousheen Jahan Shaik
M.S. Bioinformatics & Data Science, University of Delaware
Independent Research Project, 2026 | Guided by Dr Ryan Moore
