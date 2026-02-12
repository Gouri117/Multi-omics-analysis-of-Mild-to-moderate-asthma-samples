# Load libraries
library(edgeR)
library(dplyr)
library(readr)
library(forcats)
library(ggplot2)

# Step 1: Load Data ----------------------------------------------------

# KO abundance matrix (KO terms as rows, patients as columns)
ko_counts <- read.csv("kodat_table_MAinput.csv", row.names = 1, check.names = FALSE)

# Metadata
meta <- read.csv("kodat_meta_MAinput.csv")
colnames(meta)[colnames(meta) == "X.NAME"] <- "SampleID"

# Step 2: Ensure Sample Matching ---------------------------------------

# Match and reorder metadata to match count columns
common_samples <- intersect(colnames(ko_counts), meta$SampleID)
ko_counts <- ko_counts[, common_samples]
meta <- meta[match(common_samples, meta$SampleID), ]

# Sanity check
stopifnot(all(colnames(ko_counts) == meta$SampleID))

# Step 3: Create DGEList ------------------------------------------------

dge <- DGEList(counts = ko_counts)
dge <- calcNormFactors(dge)  # TMM normalization

# Step 4: Differential Expression Analysis ------------------------------

# 1. High vs Low (bleosstatus3)
meta$bleosstatus3 <- factor(meta$bleosstatus3, levels = c("Low", "High"))
design1 <- model.matrix(~ bleosstatus3, data = meta)
dge1 <- estimateDisp(dge, design1)
fit1 <- glmQLFit(dge1, design1)
res1 <- glmQLFTest(fit1, coef = 2)
top_tags_high_vs_low <- topTags(res1, n = Inf)

# Save results
write.csv(top_tags_high_vs_low$table, "High_vs_Low_KO_differential.csv")

# Extract and prepare results
res_df <- top_tags_high_vs_low$table
res_df$KO <- rownames(res_df)

# Get top 15 upregulated (highest logFC) and top 15 downregulated (lowest logFC)
top_up <- res_df %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 10)

top_down <- res_df %>%
  arrange(logFC) %>%
  slice_head(n = 10)

# Combine and label
top_combined <- bind_rows(top_up, top_down) %>%
  mutate(
    Direction = ifelse(logFC > 0, "Upregulated", "Downregulated"),
    KO = fct_reorder(KO, logFC)  # reorders for better plotting
  )

# Step 1: Download KEGG pathway table (mapXXXXX)
kegg_pathways <- read_tsv("https://rest.kegg.jp/list/pathway", col_names = FALSE)
colnames(kegg_pathways) <- c("KO", "Pathway_Name")  # This 'KO' column has 'mapXXXXX'

# Step 2: Convert mapXXXXX → koXXXXX
kegg_pathways <- kegg_pathways %>%
  mutate(KO = gsub("^map", "ko", KO))
# Merge pathway names
top_combined_annotated <- left_join(top_combined, kegg_pathways, by = "KO")

# Readable labels and order
top_combined_annotated <- top_combined_annotated %>%
  mutate(
    Label = ifelse(is.na(Pathway_Name), KO, Pathway_Name),
    Label = fct_reorder(Label, logFC)
  )

# Add stars directly into the label
top_combined_annotated <- top_combined_annotated %>%
  mutate(
    sig_label = ifelse(PValue < 0.1, "*", ""),
    Label = ifelse(is.na(Pathway_Name), KO, Pathway_Name),
    Label = paste0(sig_label, " ", Label),   # put * before term
    Label = fct_reorder(Label, logFC),
    Direction = ifelse(logFC > 0, "Enriched in High", "Enriched in Low")
  )



ggplot(top_combined_annotated, aes(x = Label, y = logFC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(
    name = "Difference (relative to T2 Low)",
    values = c("Enriched in High" = "#BB6653", "Enriched in Low" = "#640D5F"),
    labels = c("Enriched in High" = "Higher in T2 High",
               "Enriched in Low"  = "Lower in T2 High")
  ) +
  labs(
    title = "",
    x = "KEGG Pathway",
    y = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.margin = margin(10, 30, 10, 10)
  ) +
  scale_y_continuous(
    limits = c(-8.5, 6.5),
    breaks = seq(-8, 6, by = 2),
    expand = expansion(mult = c(0.02, 0.12))
  )





meta_high <- meta %>% filter(bleosstatus3 == "High")

# Match KO counts to High group samples
ko_high <- ko_counts[, colnames(ko_counts) %in% meta_high$SampleID]

# Reorder metadata to match counts
meta_high <- meta_high[match(colnames(ko_high), meta_high$SampleID), ]
stopifnot(all(colnames(ko_high) == meta_high$SampleID))

group <- factor(meta_high$Patient_Type, levels = c("ATH", "ACS"))

dge <- DGEList(counts = ko_high)
dge <- calcNormFactors(dge)

design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
res <- glmQLFTest(fit, coef = 2)  # ACS vs ATH

# Extract table
res_df <- topTags(res, n = Inf)$table
res_df$KO <- rownames(res_df)

top_up <- res_df %>% arrange(desc(logFC)) %>% slice_head(n = 10)
top_down <- res_df %>% arrange(logFC) %>% slice_head(n = 10)

top_combined <- bind_rows(top_up, top_down) %>%
  mutate(Direction = ifelse(logFC > 0, "Enriched in ACS", "Enriched in ATH"))

# Download KEGG pathway list (mapXXXXX)
kegg_pathways <- read_tsv("https://rest.kegg.jp/list/pathway", col_names = FALSE)
colnames(kegg_pathways) <- c("KO", "Pathway_Name")

# Convert mapXXXXX → koXXXXX
kegg_pathways <- kegg_pathways %>%
  mutate(KO = gsub("^map", "ko", KO))

# Merge pathway names
top_combined_annotated <- left_join(top_combined, kegg_pathways, by = "KO")

# Label for plotting
top_combined_annotated <- top_combined_annotated %>%
  mutate(
    Label = ifelse(is.na(Pathway_Name), KO, Pathway_Name),
    Label = fct_reorder(Label, logFC)
  )

# Add stars directly into the label
# Label for plotting (+ prepend * for PValue < 0.1)
top_combined_annotated <- top_combined_annotated %>%
  mutate(
    sig_label = ifelse(PValue < 0.1, "*", ""),
    Label = ifelse(is.na(Pathway_Name), KO, Pathway_Name),
    Label = paste0(sig_label, " ", Label),   # * to the LEFT of term
    Label = fct_reorder(Label, logFC)
    # IMPORTANT: do NOT overwrite Direction here
  )

ggplot(top_combined_annotated, aes(x = Label, y = logFC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(
    name   = "Difference (relative to Asthma + No steroids)",
    values = c("Enriched in ACS" = "#468A9A", "Enriched in ATH" = "#541212"),
    labels = c(
      "Enriched in ACS" = "Higher in Asthma + Inhaled steroids",
      "Enriched in ATH" = "Lower in Asthma + Inhaled steroids"
    )
  ) +
  labs(
    title = "",
    x = "KEGG Pathway",
    y = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  scale_y_continuous(
    limits = c(-8, 2),
    breaks = seq(-8, 2, by = 2),
    expand = c(0, 0)
  )




# Filter metadata for T2-Low group
meta_low <- meta %>% filter(bleosstatus3 == "Low")

# Subset KO count matrix to only T2-Low sample columns
ko_low <- ko_counts[, colnames(ko_counts) %in% meta_low$SampleID]

# Match order of metadata to count matrix columns
meta_low <- meta_low[match(colnames(ko_low), meta_low$SampleID), ]
stopifnot(all(colnames(ko_low) == meta_low$SampleID))

group <- factor(meta_low$Patient_Type, levels = c("ATH", "ACS"))

dge <- DGEList(counts = ko_low)
dge <- calcNormFactors(dge)

design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
res <- glmQLFTest(fit, coef = 2)

# Extract result table and add KO IDs as column
res_df <- topTags(res, n = Inf)$table
res_df$KO <- rownames(res_df)

top_up <- res_df %>% arrange(desc(logFC)) %>% slice_head(n = 10)
top_down <- res_df %>% arrange(logFC) %>% slice_head(n = 10)

top_combined <- bind_rows(top_up, top_down) %>%
  mutate(Direction = ifelse(logFC > 0, "Enriched in ACS", "Enriched in ATH"))

# Download KEGG pathway list (mapXXXXX)
kegg_pathways <- read_tsv("https://rest.kegg.jp/list/pathway", col_names = FALSE)
colnames(kegg_pathways) <- c("KO", "Pathway_Name")

# Convert 'mapXXXXX' → 'koXXXXX'
kegg_pathways <- kegg_pathways %>%
  mutate(KO = gsub("^map", "ko", KO))

# Merge with result table
top_combined_annotated <- left_join(top_combined, kegg_pathways, by = "KO")

# Prepare labels for plotting
top_combined_annotated <- top_combined_annotated %>%
  mutate(
    Label = ifelse(is.na(Pathway_Name), KO, Pathway_Name),
    Label = fct_reorder(Label, logFC)
  )

top_combined_annotated <- top_combined_annotated %>%
  mutate(
    sig_label = ifelse(PValue < 0.1, "*", ""),
    Label = ifelse(is.na(Pathway_Name), KO, Pathway_Name),
    Label = paste0(sig_label, " ", Label),   # * to the LEFT of term
    Label = fct_reorder(Label, logFC)
    # IMPORTANT: do NOT overwrite Direction here
  )


ggplot(top_combined_annotated, aes(x = Label, y = logFC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(name = "Difference (relative to Asthma + No steroids)",
                    values = c("Enriched in ACS" = "#468A9A", "Enriched in ATH" = "#541212"),
                    labels = c(
                      "Enriched in ACS" = "Higher in Asthma + Inhaled steroids",
                      "Enriched in ATH" = "Low in Asthma + Inhaled steroids"
                    )) +
  labs(
    title = "",
    x = "KEGG Pathway",
    y = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  scale_y_continuous(
    limits = c(-6, 6),
    breaks = seq(-6, 6, by = 2),
    expand = c(0, 0)
  )

library(edgeR)
library(dplyr)
library(readr)
library(ggplot2)
library(forcats)

# Step 1: Load KO counts and metadata
ko_counts <- read.csv("kodat_full_values.csv", row.names = 1, check.names = FALSE)
ko_counts <- t(ko_counts)

meta <- read.csv("eos300_wgs_clinical_meta.csv")
colnames(meta)[colnames(meta) == "X.NAME"] <- "SampleID"

# Step 2: Match sample IDs between KO data and metadata
common_samples <- intersect(colnames(ko_counts), meta$SampleID)
ko_counts <- ko_counts[, common_samples]
meta <- meta[match(common_samples, meta$SampleID), ]
stopifnot(all(colnames(ko_counts) == meta$SampleID))

# Step 3: Create DGEList and normalize
dge <- DGEList(counts = ko_counts)
dge <- calcNormFactors(dge)

# Step 4: Differential expression analysis - Asthma vs Healthy
meta$status <- factor(meta$status, levels = c("Healthy", "Asthma"))
design1 <- model.matrix(~ status, data = meta)
dge1 <- estimateDisp(dge, design1)
fit1 <- glmQLFit(dge1, design1)
res1 <- glmQLFTest(fit1, coef = 2)
res_df <- topTags(res1, n = Inf)$table
res_df$KO <- rownames(res_df)

# Step 5: Select top 10 enriched in each direction
top_up <- res_df %>% arrange(desc(logFC)) %>% slice_head(n = 10)
top_down <- res_df %>% arrange(logFC) %>% slice_head(n = 10)

top_combined <- bind_rows(top_up, top_down) %>%
  mutate(Direction = ifelse(logFC > 0, "Enriched in Asthma", "Enriched in Healthy"))

# Step 6: Map KO IDs to KEGG Pathway Names
kegg_pathways <- read_tsv("https://rest.kegg.jp/list/pathway", col_names = FALSE)
colnames(kegg_pathways) <- c("KO", "Pathway_Name")
kegg_pathways <- kegg_pathways %>%
  mutate(KO = gsub("^map", "ko", KO))  # map00120 → ko00120

top_combined_annotated <- left_join(top_combined, kegg_pathways, by = "KO") %>%
  mutate(
    Label = ifelse(is.na(Pathway_Name), KO, Pathway_Name),
    Label = fct_reorder(Label, logFC)
  )

# Step 7: Plot
ggplot(top_combined_annotated, aes(x = Label, y = logFC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(name = "Difference (relative to Healthy)",
                    values = c("Enriched in Asthma" = "#154D71", "Enriched in Healthy" = "#568F87"),
                    labels = c(
                      "Enriched in Asthma" = "Higher in Asthma",
                      "Enriched in Healthy" = "Lower in Asthma"
                    )) +
  labs(
    title = "",
    x = "KEGG Pathway",
    y = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  scale_y_continuous(
    limits = c(-4, 8),
    breaks = seq(-4, 8, by = 2),
    expand = c(0, 0)
  )

