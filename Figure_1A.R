# Figure 1A
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Input is in /Users/gourianil/Library/CloudStorage/Dropbox-UniversityofMichigan/Gouri Anil/Kozik Lab - Dry Lab/Kozik Lab - Gouri/Consolidated/Relative abundance plot/Asthma_vs_Healthy_taxonomy_Fig1AB.csv

# Load your taxa abundance file
taxa_abundance <- read.csv("Asthma_vs_Healthy_taxonomy_Fig1AB.csv", header = TRUE)

# Replace underscores
taxa_abundance$X <- gsub("_", " ", taxa_abundance$X)

# Convert to long format
taxa_long <- taxa_abundance %>%
  pivot_longer(cols = c(Asthma, Healthy),
               names_to = "Group",
               values_to = "Abundance")

# Get top 30 species
top_taxa <- taxa_long %>%
  group_by(X) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  arrange(desc(MeanAbundance)) %>%
  slice(1:30) %>%
  pull(X)

# Label non-top as "Others"
taxa_long <- taxa_long %>%
  mutate(Species = ifelse(X %in% top_taxa, X, "Others"))

# Aggregate and normalize
taxa_long <- taxa_long %>%
  group_by(Group, Species) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(RelativeAbundance = TotalAbundance / sum(TotalAbundance) * 100)

# Order species within each group (Others first, then highest abundance)
taxa_long <- taxa_long %>%
  arrange(Group, Species == "Others", desc(RelativeAbundance))

# Define darker color palette
palette_dark <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(taxa_long$Species)))

# Factor Species for correct stacking order
taxa_long$Species <- factor(taxa_long$Species, levels = unique(taxa_long$Species))

# Vertical stacked bar plot
# Vertical stacked bar plot with smaller legend
ggplot(taxa_long, aes(x = Group, y = RelativeAbundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  labs(
    title = "",
    x = "Group",
    y = "Relative Abundance (%)",
    fill = "Taxon (Species)"  # Legend title
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10, face = "italic"),   # italic legend labels
    legend.title = element_text(size = 15)   
  ) +
  scale_fill_manual(values = palette_dark)

