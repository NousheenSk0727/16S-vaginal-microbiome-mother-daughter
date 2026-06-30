# Workflow
# 1. Quality assessment of trimmed FASTQ files
# 2. Filtering and trimming of reads
# 3. Error learning
# 4. Dereplication and sample inference
# 5. Paired-end merging
# 6. Sequence table construction
# 7. Chimera removal
# 8. Read tracking through the pipeline
# 9. Taxonomy assignment
# 10. Import into phyloseq

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("dada2", "ShortRead", "phyloseq"))

library(dada2)

# dada2 is used for denoising, merging, chimera removal,
# and taxonomy assignment of amplicon sequencing data.
library(ShortRead)
library(ggplot2)
# ggplot2 is used for plotting quality profiles and
# downstream alpha diversity visualizations.
library(phyloseq)
library(Biostrings)
# Biostrings is used to store DNA sequences inside the phyloseq object.

theme_set(theme_bw())

# SET INPUT PATH AND IDENTIFY FORWARD / REVERSE READ FILES
path <- "/Users/nousheenjahanshaik/Documents/trimmed_fastq"

list.files(path)

# Forward and reverse FASTQ filenames follow the format:
fnFs <- sort(list.files(path, pattern = "_1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2_trimmed.fastq.gz", full.names = TRUE))

# Extract sample names from filenames by splitting at "_".
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# QUALITY CONTROL (QC) OF INPUT READS
# Plot quality profiles for the first two forward-read files.
# This helps inspect where sequence quality starts to decline,
# which guides decisions about truncation length.
plotQualityProfile(fnFs[1:2])

# Check that the number of forward and reverse files is equal.
# This confirms that the dataset is properly paired-end.
length(fnFs) == length(fnRs)

# Save combined QC plots for all samples into one PDF file.
# Each page contains forward + reverse read quality profile for a single sample.
pdf("All_QC_profiles_combined.pdf", width = 12, height = 6)
for (i in seq_along(fnFs)) {
  p1 <- plotQualityProfile(fnFs[i]) + ggtitle("Forward")
  p2 <- plotQualityProfile(fnRs[i]) + ggtitle("Reverse")
  gridExtra::grid.arrange(p1, p2, ncol = 2, top = basename(fnFs[i]))
}
dev.off ()

# FILTERING AND TRIMMING
dir.create(file.path(path, "filtered"), showWarnings = FALSE)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# truncLen = c(240, 220) trims forward reads to 240 bp and reverse reads to 220 bp
# maxN = 0 removes reads containing ambiguous bases (N)
# maxEE = c(2, 5) sets maximum expected errors for forward and reverse reads
# truncQ = 2 truncates reads at the first quality score <= 2
# rm.phix = TRUE removes PhiX contamination
# Ran in batches of 10 samples - running all 87 at once kept hanging on my machine
batch_ranges <- list(1:10, 11:20, 21:30, 31:40, 41:50,
                     51:60, 61:70, 71:80, 81:87)

out_list <- lapply(batch_ranges, function(idx) {
  filterAndTrim(
    fnFs[idx], filtFs[idx], fnRs[idx], filtRs[idx],
    truncLen = c(240, 220),
    maxN     = 0,
    maxEE    = c(2, 5),
    truncQ   = 2,
    rm.phix  = TRUE,
    compress = TRUE,
    multithread = 2
  )
})

out <- do.call(rbind, out_list)
gc()
head(out)

# LEARN ERROR RATES
# DADA2 uses this learned error model to distinguish sequencing errors
# from real biological sequence variation.
errF <- learnErrors(filtFs, multithread = 2)

filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
errR <- learnErrors(filtRs, multithread = 2)

# The observed error rates should roughly follow the fitted black line.
plotErrors(errF, nominalQ = TRUE)

# DEREPLICATION AND SAMPLE INFERENCE
# Dereplication collapses identical reads into unique sequences with abundances.
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

sample.names <- sample.names[file.exists(filtFs)]
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# This step identifies true sequence variants (ASVs) while correcting errors.
dadaFs <- dada(derepFs, err = errF, multithread = 2)
dadaRs <- dada(derepRs, err = errR, multithread = 2)

dadaFs[[1]]

# MERGE PAIRED-END READS
# Reconstructs the full amplicon sequence by overlapping the paired reads.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers[[1]])

# CONSTRUCT SEQUENCE TABLE
# Rows = samples, columns = unique sequence variants, values = counts
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# REMOVE CHIMERIC SEQUENCES
# Chimeras are artificial sequences formed during PCR and should be removed
# before downstream ecological analysis.
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 2, verbose = TRUE)
dim(seqtab.nochim)

# Proportion of reads retained after chimera removal - high = good data quality
sum(seqtab.nochim) / sum(seqtab)

# TRACK READS THROUGH THE PIPELINE
getN <- function(x) sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# TAXONOMY ASSIGNMENT
# The SILVA reference file provides taxonomic labels for bacterial 16S sequences.
silva_path <- "~/Downloads/silva_nr99_v138.1_train_set.fa.gz"
file.exists(silva_path)

taxa <- assignTaxonomy(seqtab.nochim, silva_path, multithread = 2)
taxa[1:5, ]

write.csv(taxa, "~/Desktop/taxonomy_assignments.csv", quote = FALSE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# LOAD AND PREPARE SAMPLE METADATA
metadata <- read.csv("/Users/nousheenjahanshaik/Downloads/SraRunTable.csv")

samples.out <- rownames(seqtab.nochim)
samdf <- metadata[metadata$Run %in% samples.out, ]
rownames(samdf) <- samdf$Run

# M = Mother, D = Daughter
samdf$Group <- ifelse(samdf$MotherDaughter == "M", "Mother", "Daughter")

# CREATE PHYLOSEQ OBJECT
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)

# Add DNA sequences for each ASV into the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)

# SAVE KEY OBJECTS
saveRDS(seqtab.nochim, file = "/Users/nousheenjahanshaik/Desktop/seqtab.nochim.rds")
saveRDS(samdf, file = "/Users/nousheenjahanshaik/Desktop/srr.rds")
