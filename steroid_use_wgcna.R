# Install/Load ----
suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(magrittr)
  library(tidyverse)   # ggplot2, readr, tidyr, etc.
  library(WGCNA)
  library(DESeq2)
  library(gridExtra)
  library(psych)
  library(CorLevelPlot)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(broom)       # for tidy() on lm fits
})

cor <- WGCNA::cor

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

# Load counts and pre=processing
# Counts ----
data <- read.csv("all_counts.csv", stringsAsFactors = FALSE) %>%
  na.omit() %>%
  group_by(Gene = .[[1]]) %>%          # collapse duplicate gene rows if any
  summarise(across(-1, ~ round(sum(.)))) %>%
  as.data.frame() %>%
  { rownames(.) <- .$Gene; .[, -1, drop = FALSE] }

# Metadata (samples in rows) ----
phenoData <- read.csv("metadata.csv", row.names = 1, stringsAsFactors = FALSE)

# Basic alignment checks (samples = columns of counts)
stopifnot(all(rownames(phenoData) %in% colnames(data)))
phenoData <- phenoData[colnames(data), , drop = FALSE]


# Quality control
# WGCNA expects samples in rows, genes in columns -> transpose for diagnostic only
gsg <- goodSamplesGenes(t(data))
if (!gsg$allOK) {
  data <- data[gsg$goodGenes, , drop = FALSE]
}

# (Optional) outlier inspection ----
htree <- hclust(dist(t(data)), method = "average")
plot(htree, main = "Sample clustering (raw)")

pca <- prcomp(t(data))
pca.dat <- as.data.frame(pca$x)
pvar <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
ggplot(pca.dat, aes(PC1, PC2, label = rownames(pca.dat))) +
  geom_point() + geom_text(nudge_y = 0.5, size = 3) +
  labs(x = paste0("PC1: ", pvar[1], "%"), y = paste0("PC2: ", pvar[2], "%")) +
  theme_classic()

# Remove specific outliers if needed (empty here)
samples.to.be.excluded <- character(0)
data.subset <- data[, !(colnames(data) %in% samples.to.be.excluded), drop = FALSE]
phenoData.subset <- phenoData[!(rownames(phenoData) %in% samples.to.be.excluded), , drop = FALSE]

# Sanity
stopifnot(identical(rownames(phenoData.subset), colnames(data.subset)))


# Preprocess
# Build DESeq2 object ----
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData   = phenoData.subset,
                              design    = ~ 1)

# Filter low counts: keep genes with >= 15 counts in at least 75% samples
min_count <- 15
keep_n    <- ceiling(ncol(dds) * 0.75)
dds75 <- dds[rowSums(counts(dds) >= min_count) >= keep_n, ]
message("Genes kept after filter: ", nrow(dds75))

# Variance stabilizing transform -> WGCNA-friendly expression
dds_norm <- vst(dds75, blind = TRUE)

# WGCNA format: samples in rows, genes in columns
norm.counts <- t(assay(dds_norm))

# For safety, everything numeric
storage.mode(norm.counts) <- "double"

# Choose power and build network
# Soft-threshold selection (signed network) ----
powers <- c(1:10, 12, 14, 16, 18, 20)
sft <- pickSoftThreshold(norm.counts, powerVector = powers,
                         networkType = "signed", verbose = 5)
fi <- sft$fitIndices
soft_power <- if (any(fi$SFT.R.sq >= 0.8)) {
  fi$Power[min(which(fi$SFT.R.sq >= 0.8))]
} else fi$Power[which.max(fi$SFT.R.sq)]
message("Chosen soft power = ", soft_power)

# Visual check (optional)
a1 <- ggplot(fi, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() + geom_text(nudge_y = 0.06, size = 3) +
  geom_hline(yintercept = 0.8, color = "red") + theme_classic() +
  labs(y = "Scale-free Topology (R^2)")
a2 <- ggplot(fi, aes(Power, mean.k., label = Power)) +
  geom_point() + geom_text(vjust = -0.8, size = 3) +
  theme_classic() + labs(y = "Mean connectivity")
gridExtra::grid.arrange(a1, a2, nrow = 2)

# Network ----
bwnet <- blockwiseModules(
  norm.counts,
  power             = soft_power,
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = TRUE,
  maxBlockSize      = 20000,
  verbose           = 3,
  corType           = "pearson",
  corFnc            = "WGCNA::cor",          # <- force WGCNA's cor()
  corOptions        = "use = 'p'"            # pairwise complete obs
)

module_colors <- bwnet$colors
module_eigengenes <- bwnet$MEs          # rows = samples, cols = MEs

# Interaction model
# Design: Subtype(High=1, Low=0) and Steroid(ACS=1) ----
design <- phenoData.subset %>%
  mutate(
    Subtype = case_when(T2_Status == "High" ~ 1L,
                        T2_Status == "Low"  ~ 0L,
                        TRUE ~ NA_integer_),      # drop NA/Healthy for this part
    Steroid = as.integer(grepl("ACS", patient_type))
  ) %>%
  select(Subtype, Steroid) %>%
  filter(!is.na(Subtype))

# Align to MEs
common <- intersect(rownames(design), rownames(module_eigengenes))
design <- design[common, , drop = FALSE]
MEs_for_lm <- module_eigengenes[common, , drop = FALSE]
stopifnot(identical(rownames(design), rownames(MEs_for_lm)))

# Fit ME ~ Subtype * Steroid for each module
results_list <- lapply(colnames(MEs_for_lm), function(me) {
  df <- cbind(ME = MEs_for_lm[, me], design)
  fit <- lm(ME ~ Subtype * Steroid, data = df)
  broom::tidy(fit) %>% mutate(Module = me)
})
results_df <- bind_rows(results_list)

# Heatmap of coefficients
coef_mat <- results_df %>%
  select(Module, term, estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate)

p_interaction <- ggplot(coef_mat %>%
                          tidyr::pivot_longer(-Module, names_to = "Effect", values_to = "Beta"),
                        aes(Effect, Module, fill = Beta)) +
  geom_tile(color = "black") +
  geom_text(aes(label = sprintf("%.2f", Beta)), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#113F67", mid = "white", high = "#DC143C", midpoint = 0) +
  labs(title = "Subtype × Steroid effects on Module Eigengenes", x = "", y = "", fill = "β") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_interaction)


# 4 column subgroup heatmap
# Helper: lock ME signs to module gene means (prevents arbitrary flips) ----
align_ME_sign <- function(MEs, datExpr, colors){
  out <- MEs
  for (ME in colnames(MEs)) {
    col <- sub("^ME","", ME)
    genes <- names(colors)[colors == col]
    genes <- intersect(genes, colnames(datExpr))
    if (length(genes) >= 2) {
      anchor <- rowMeans(datExpr[rownames(MEs), genes, drop = FALSE])
      if (suppressWarnings(cor(MEs[[ME]], anchor, use = "pairwise.complete.obs")) < 0)
        out[[ME]] <- -out[[ME]]
    }
  }
  WGCNA::orderMEs(out)
}

# Meta restricted to High/Low, plus Steroid flag ----
meta <- phenoData.subset %>%
  transmute(
    T2_Status2 = ifelse(T2_Status %in% c("High","Low"), T2_Status, NA),
    Steroid    = ifelse(grepl("ACS", patient_type), "Steroid", "NoSteroid")
  )

# Align meta & MEs and drop NA (Healthy) ----
common2 <- intersect(rownames(meta), rownames(module_eigengenes))
meta    <- meta[common2, , drop = FALSE]
MEs_all <- module_eigengenes[common2, , drop = FALSE]
meta    <- meta[!is.na(meta$T2_Status2), , drop = FALSE]
MEs_all <- MEs_all[rownames(meta), , drop = FALSE]

# OPTIONAL: lock ME signs (recommended for consistent plots)
MEs_all <- align_ME_sign(MEs_all, norm.counts[rownames(MEs_all), , drop = FALSE], module_colors)

# Build 4-group indicator with FIXED order ----
meta$grp <- factor(paste(meta$T2_Status2, meta$Steroid, sep = "_"),
                   levels = c("High_Steroid","High_NoSteroid","Low_Steroid","Low_NoSteroid"))
G <- model.matrix(~ 0 + grp, data = meta)
colnames(G) <- c("T2High_Steroid","T2High_NoSteroid","T2Low_Steroid","T2Low_NoSteroid")

# Sanity checks (avoid subtle bugs)
stopifnot(identical(rownames(MEs_all), rownames(G)))

# Correlate MEs (rows=samples) vs subgroup indicators (rows=samples)
corr_mat <- WGCNA::bicor(MEs_all, G, use = "pairwise.complete.obs")   # modules × 4 groups
p_mat    <- WGCNA::corPvalueStudent(corr_mat, nrow(MEs_all))



# Plot function ----
p_subgroups <- {
  df <- reshape2::melt(corr_mat, varnames = c("Module","Group"), value.name = "Correlation") %>%
    dplyr::mutate(Module = as.character(Module)) %>%
    dplyr::group_by(Module) %>% dplyr::mutate(.ord = max(abs(Correlation), na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::mutate(Module = forcats::fct_reorder(Module, .ord, .desc = TRUE)) %>%
    dplyr::select(-.ord)
  
  ggplot(df, aes(Group, Module, fill = Correlation)) +
    geom_tile(color = "black", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3, fontface = "bold") +
    scale_fill_gradient2(
      low = "#113F67", mid = "white", high = "#DC143C",
      midpoint = 0,
      limits = c(-0.3, 0.3),     # <<< legend & color scale from -0.3 to 0.3
      oob = scales::squish
    ) +
    labs(title = "Modules vs T2-High/Low × Steroid/NoSteroid", x = "", y = "", fill = "r") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
print(p_subgroups)




# ==== Save genes grouped into ME-color columns ====

# Ensure module_colors are named by gene IDs
genes <- colnames(norm.counts)
mods  <- module_colors
if (is.null(names(mods))) names(mods) <- genes

library(dplyr)
library(tidyr)

# Build a table (Module, Gene)
mod_gene_tbl <- tibble(
  Module = as.character(mods[genes]),
  Gene   = genes
)

# Split into wide format: one column per Module
mod_gene_wide <- mod_gene_tbl %>%
  group_by(Module) %>%
  mutate(Row = row_number()) %>%         # create row index within each module
  ungroup() %>%
  pivot_wider(
    id_cols   = Row,
    names_from = Module,
    values_from = Gene
  ) %>%
  select(-Row)

# Save to CSV
out_path <- "/Users/gourianil/Library/CloudStorage/Dropbox-UniversityofMichigan/Gouri Anil/Kozik Lab - Dry Lab/Kozik Lab - Gouri/Consolidated/WGCNA/wgcna_module_genes_by_column.csv"
write.csv(mod_gene_wide, out_path, row.names = FALSE)

message("Saved wide-format file: ", out_path)




