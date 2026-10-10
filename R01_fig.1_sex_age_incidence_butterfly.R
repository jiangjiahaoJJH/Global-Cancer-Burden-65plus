# =============================================================================
# Figure 1: Sex- and age-specific distribution of cancer incidence (65+)
# -----------------------------------------------------------------------------
# Description:
#   This script generates a butterfly (mirrored stacked bar) plot showing the
#   proportional distribution of cancer types by age group and sex among
#   adults aged 65 years and older, using GLOBOCAN 2024 incidence data.
#   It also exports supplementary tables for proportions, counts, and rates.
#
# Input:
#   - Processed_65plus_global_age_ASR.csv
#     (contains columns: sex_label, age_label, type_label, cancer_label,
#      total, rate)
#
# Outputs:
#   - Fig1_sex_age_incidence_proportion.csv
#   - Table_S_incidence_proportion.csv
#   - Table_S_incidence_count.csv
#   - Table_S_incidence_rate.csv
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
setwd("path/to/working/directory")   # <- update to your local path

library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)

# -----------------------------------------------------------------------------
# 1. Load and filter data
# -----------------------------------------------------------------------------
IBD_china <- fread("Processed_65plus_global_age_ASR.csv")

data1 <- subset(
  IBD_china,
  sex_label %in% c("Male", "Female") &
    age_label %in% c("65-69", "70-74", "75-79", "80-84", "85+") &
    type_label == "Incidence"
)

# Remove rows with total = 0
data1 <- subset(data1, total != 0)

# Recode cancer labels
data1 <- data1 %>%
  mutate(cancer_label = recode(
    cancer_label,
    "Brain, central nervous system" = "Brain & CNS",
    "Liver and intrahepatic bile ducts" = "Liver & IBD"
  ))

# -----------------------------------------------------------------------------
# 2. Compute proportions using "All cancers" as denominator
# -----------------------------------------------------------------------------
all_cancer_total <- data1 %>%
  filter(cancer_label == "All cancers excl. non-melanoma skin cancer") %>%
  select(age_label, sex_label, total) %>%
  rename(all_cancer_total = total)

# Remove "All cancers" from the plotting data
data1 <- subset(data1, cancer_label != "All cancers excl. non-melanoma skin cancer")

result <- data1 %>%
  left_join(all_cancer_total, by = c("age_label", "sex_label")) %>%
  mutate(val_ratio = round(total / all_cancer_total * 100, 2))

# Male proportions are plotted on the negative side
result <- result %>%
  mutate(val_ratio = ifelse(sex_label == "Male", -val_ratio, val_ratio))

# Age group order (top to bottom)
result$age_label <- factor(
  result$age_label,
  levels = rev(c("65-69", "70-74", "75-79", "80-84", "85+"))
)

# Cancer type order
cause_levels <- c(
  "Trachea, bronchus and lung", "Colorectum", "Prostate", "Breast",
  "Stomach", "Liver & IBD", "Bladder", "Pancreas", "Oesophagus",
  "Non-Hodgkin lymphoma", "Leukaemia", "Kidney", "Melanoma of skin",
  "Corpus uteri", "Lip, oral cavity", "Cervix uteri", "Thyroid",
  "Multiple myeloma", "Brain & CNS", "Ovary", "Larynx",
  "Gallbladder", "Oropharynx", "Hypopharynx", "Nasopharynx",
  "Vulva", "Salivary glands", "Mesothelioma", "Penis",
  "Hodgkin lymphoma", "Vagina", "Kaposi sarcoma", "Testis"
)

result$cancer_label <- factor(result$cancer_label, levels = cause_levels)

# Colour palette
colors <- c(
  "#f89c59", "#dcdcdc", "#20b2aa", "#FF6B6B",
  "#5c3c92", "#1D395F", "#2F4F2F", "#008000", "#add8e6",
  "#9a3e5a", "#8b008b", "#ffb6c1", "#FFD700",
  "#6a5acd", "#3C1E40", "#bc8f8f", "#e84393",
  "#2a9d8f", "#4169e1", "#ff4500", "#63ADEE",
  "#FFFACD", "#e22627", "#0d4e92", "#e9f5f3",
  "#c07020", "#3d5a80", "#1a1a2e", "#4682B4",
  "#7B3F00", "#8fbc8f", "#b31f39", "#5F9EA0"
)

# -----------------------------------------------------------------------------
# 3. Butterfly plot
# -----------------------------------------------------------------------------
p_result <- ggplot(
  data = result,
  mapping = aes(x = val_ratio, y = age_label, fill = cancer_label)
) +
  geom_bar(
    stat = "identity", position = "stack", width = 0.9,
    alpha = 0.9, color = "black", linewidth = 0.15
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  scale_x_continuous(
    name = "Proportion of new cases (%)",
    breaks = seq(-100, 100, by = 25),
    labels = abs(seq(-100, 100, by = 25)),
    limits = c(-100.1, 100.1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = colors, name = "") +
  labs(y = "", title = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, size = 7, color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.25),
    axis.line = element_line(color = "black", size = 0.25),
    axis.ticks = element_line(color = "black", size = 0.25),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  ) +
  guides(fill = guide_legend(keysize = 0.4))

print(p_result)

# Recommended export size: 3 x 9 inches (male left, female right)

# -----------------------------------------------------------------------------
# 4. Supplementary table 1: proportions (val_ratio)
# -----------------------------------------------------------------------------
newdata <- result[, c("age_label", "sex_label", "cancer_label", "val_ratio")]
newdata <- newdata %>%
  mutate(cancer_label = recode(
    cancer_label,
    "Trachea, bronchus and lung" = "Lung",
    "Non-Hodgkin lymphoma" = "NHL"
  ))

newdata <- newdata %>%
  pivot_wider(names_from = cancer_label, values_from = val_ratio)

desired_order <- c(
  "Lung", "Colorectum", "Prostate", "Breast",
  "Stomach", "Liver & IBD", "Bladder", "Pancreas", "Oesophagus",
  "NHL", "Leukaemia", "Kidney", "Melanoma of skin",
  "Corpus uteri", "Lip, oral cavity", "Cervix uteri", "Thyroid",
  "Multiple myeloma", "Brain & CNS", "Ovary", "Larynx",
  "Gallbladder", "Oropharynx", "Hypopharynx", "Nasopharynx",
  "Vulva", "Salivary glands", "Mesothelioma", "Penis",
  "Hodgkin lymphoma", "Vagina", "Kaposi sarcoma", "Testis"
)

newdata <- newdata %>%
  select(age_label, sex_label, all_of(desired_order)) %>%
  arrange(age_label, sex_label)

write.csv(newdata, "Table_S_incidence_proportion.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Supplementary table 2: incidence counts (total)
# -----------------------------------------------------------------------------
newdata2 <- result[, c("age_label", "sex_label", "cancer_label", "total")]
newdata2 <- newdata2 %>%
  mutate(cancer_label = recode(
    cancer_label,
    "Trachea, bronchus and lung" = "Lung",
    "Non-Hodgkin lymphoma" = "NHL"
  ))

newdata2 <- newdata2 %>%
  pivot_wider(names_from = cancer_label, values_from = total)

newdata2 <- newdata2 %>%
  select(age_label, sex_label, all_of(desired_order)) %>%
  arrange(age_label, sex_label)

write.csv(newdata2, "Table_S_incidence_count.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Supplementary table 3: incidence rates (rate)
# -----------------------------------------------------------------------------
newdata3 <- result[, c("age_label", "sex_label", "cancer_label", "rate")]
newdata3 <- newdata3 %>%
  mutate(cancer_label = recode(
    cancer_label,
    "Trachea, bronchus and lung" = "Lung",
    "Non-Hodgkin lymphoma" = "NHL"
  ))

newdata3 <- newdata3 %>%
  pivot_wider(names_from = cancer_label, values_from = rate)

newdata3 <- newdata3 %>%
  select(age_label, sex_label, all_of(desired_order)) %>%
  arrange(age_label, sex_label)

write.csv(newdata3, "Table_S_incidence_rate.csv", row.names = FALSE)