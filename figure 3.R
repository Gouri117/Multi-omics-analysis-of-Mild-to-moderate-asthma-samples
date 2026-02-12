# Load required libraries
library(tidyverse)
library(ggpubr)

# Read data
data <- read.csv("asthma_metadata_with_sampleid_status.csv")

# Exclude specific Subject IDs
excluded_ids <- c(6113032, 6123016, 6181101)
data <- data[!data$Subject %in% excluded_ids, ]

# Ensure Status is a factor
data$Status <- factor(data$Status, levels = c("Low", "High"))

### ===== IL4_ccdt2 =====
il4_data <- data %>% filter(!is.na(IL4_ccdt2))
il4_test <- wilcox.test(IL4_ccdt2 ~ Status, data = il4_data, exact = FALSE)
cat("IL4 Wilcoxon p-value:", il4_test$p.value, "\n")

p_il4 <- ggboxplot(il4_data, x = "Status", y = "IL4_ccdt2",
                   fill = "Status", palette = "jco", add = "jitter") +
  scale_y_log10() +
  labs(
    title = "IL-4 levels by T2 Status",
    subtitle = paste("p-value:", signif(il4_test$p.value, 3)),
    x = "T2 Status",
    y = "IL-4 concentration (log10 scale)"
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


print(p_il4)


### ===== IL5_ccdt2 =====
il5_data <- data %>% filter(!is.na(IL5_ccdt2))
il5_test <- wilcox.test(IL5_ccdt2 ~ Status, data = il5_data, exact = FALSE)
cat("IL5 Wilcoxon p-value:", il5_test$p.value, "\n")

p_il5 <- ggboxplot(il5_data, x = "Status", y = "IL5_ccdt2",
                   fill = "Status", palette = "jco", add = "jitter") +
  # Removed stat_compare_means to hide p-value box
  scale_y_log10() +
  labs(
    title = "IL-5 levels by T2 Status",
    subtitle = paste("p-value:", signif(il5_test$p.value, 3)),
    x = "T2 Status",
    y = "IL-5 concentration (log10 scale)"
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


print(p_il5)


### ===== IL13_ccdt2 =====
il13_data <- data %>% filter(!is.na(IL13_ccdt2))
il13_test <- wilcox.test(IL13_ccdt2 ~ Status, data = il13_data, exact = FALSE)
cat("IL13 Wilcoxon p-value:", il13_test$p.value, "\n")

p_il13 <- ggboxplot(il13_data, x = "Status", y = "IL13_ccdt2",
                    fill = "Status", palette = "jco", add = "jitter") +
  # Removed stat_compare_means to hide p-value box
  scale_y_log10() +
  labs(
    title = "IL-13 levels by T2 Status",
    subtitle = paste("p-value:", signif(il13_test$p.value, 3)),
    x = "T2 Status",
    y = "IL-13 concentration (log10 scale)"
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

print(p_il13)


