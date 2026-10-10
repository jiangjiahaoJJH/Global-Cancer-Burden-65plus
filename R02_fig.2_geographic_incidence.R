# =============================================================================
# Figure 2: Geographic variation in cancer incidence among adults aged 65+
# -----------------------------------------------------------------------------
# Description:
#   This script generates a stacked bar chart showing the proportional
#   distribution of cancer types across 21 UN regions and 6 WHO regions,
#   using GLOBOCAN 2024 age-standardised incidence rates (ASR) for adults
#   aged 65 years and older. It also exports supplementary tables of
#   proportions and counts.
#
# Input:
#   - GLOBOCAN_65plus_ASR_incidence_by_region.csv
#     (contains columns: Sex, Country label, Cancer label, Number)
#
# Outputs:
#   - Table_S_incidence_proportion_by_region.csv
#   - Table_S_incidence_count_by_region.csv
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
setwd("path/to/working/directory")   # <- update to your local path

library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(data.table)
library(tidyr)
library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
data <- fread("GLOBOCAN_65plus_ASR_incidence_by_region.csv")
cat("Original data dimensions:", dim(data), "\n")
cat("Column names:", colnames(data), "\n\n")

# -----------------------------------------------------------------------------
# 2. Data cleaning
# -----------------------------------------------------------------------------

# Keep Both sexes only
data1 <- subset(data, Sex %in% c("Both"))

# Recode cancer labels
data1 <- data1 %>%
  mutate(`Cancer label` = recode(
    `Cancer label`,
    "Brain, central nervous system" = "Brain & CNS",
    "Liver and intrahepatic bile ducts" = "Liver & IBD"
  ))

# Recode region labels
data1 <- data1 %>%
  mutate(`Country label` = recode(
    `Country label`,
    "WHO Africa region (AFRO)"             = "WHO Africa region",
    "WHO Americas region (PAHO)"           = "WHO Americas region",
    "WHO East Mediterranean region (EMRO)" = "WHO East Mediterranean region",
    "WHO Europe region (EURO)"             = "WHO Europe region",
    "WHO South-East Asia region (SEARO)"   = "WHO South-East Asia region",
    "WHO Western Pacific region (WPRO)"    = "WHO Western Pacific region"
  ))

cat("Number of cancer types:", length(unique(data1$`Cancer label`)), "\n")
cat("Number of regions:", length(unique(data1$`Country label`)), "\n\n")

# -----------------------------------------------------------------------------
# 3. Compute proportions using "All cancers" as denominator
# -----------------------------------------------------------------------------

# Extract "All cancers" totals as denominator
all_cancer_total <- data1 %>%
  filter(`Cancer label` == "All cancers excluded NMSC") %>%
  group_by(`Country label`) %>%
  summarise(all_cancer_total = sum(Number, na.rm = TRUE), .groups = "drop")

# Remove "All cancers" from plotting data
data1 <- data1 %>%
  filter(`Cancer label` != "All cancers excluded NMSC")

# Aggregate incidence by region and cancer
data_agg <- data1 %>%
  group_by(`Country label`, `Cancer label`) %>%
  summarise(Number = sum(Number, na.rm = TRUE), .groups = "drop")

# Merge denominator and compute proportions
result <- data_agg %>%
  left_join(all_cancer_total, by = "Country label") %>%
  mutate(val_ratio = round((Number / all_cancer_total) * 100, 2))

cat("Processed data dimensions:", dim(result), "\n\n")

# Verify: check that proportions sum to ~100 within each region
verification <- result %>%
  group_by(`Country label`) %>%
  summarise(sum_ratio = sum(val_ratio, na.rm = TRUE), .groups = "drop")
cat("Verification of proportion sums by region:\n")
print(verification)
cat("\n")

# -----------------------------------------------------------------------------
# 4. Set factor order
# -----------------------------------------------------------------------------

# Region order (top to bottom)
region_levels <- c(
  # WHO regions
  "WHO Africa region",
  "WHO Americas region",
  "WHO East Mediterranean region",
  "WHO Europe region",
  "WHO South-East Asia region",
  "WHO Western Pacific region",
  # UN subregions
  # Africa
  "Northern Africa",
  "Eastern Africa",
  "Middle Africa",
  "Southern Africa",
  "Western Africa",
  # Americas
  "Northern America",
  "Caribbean",
  "Central America",
  "South America",
  # Eastern Mediterranean
  "Western Asia",
  # Europe
  "Western Europe",
  "Northern Europe",
  "Southern Europe",
  "Eastern Europe",
  # South-East Asia
  "South-Eastern Asia",
  "South Central Asia",
  # Western Pacific / Oceania
  "Eastern Asia",
  "Australia-New Zealand",
  "Melanesia",
  "Micronesia",
  "Polynesia"
)

# Keep only regions present in the data
region_levels <- intersect(region_levels, unique(result$`Country label`))
cat("Number of region factor levels:", length(region_levels), "\n")

result$`Country label` <- factor(result$`Country label`, levels = region_levels)

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

# Keep only cancers present in the data
cause_levels <- intersect(cause_levels, unique(result$`Cancer label`))
cat("Number of cancer factor levels:", length(cause_levels), "\n\n")

result$`Cancer label` <- factor(result$`Cancer label`, levels = cause_levels)

# -----------------------------------------------------------------------------
# 5. Colour palette
# -----------------------------------------------------------------------------
colors <- c(
  "#f89c59",  # 1  Trachea, bronchus and lung
  "#dcdcdc",  # 2  Colorectum
  "#20b2aa",  # 3  Prostate
  "#FF6B6B",  # 4  Breast
  "#5c3c92",  # 5  Stomach
  "#1D395F",  # 6  Liver & IBD
  "#2F4F2F",  # 7  Bladder
  "#008000",  # 8  Pancreas
  "#add8e6",  # 9  Oesophagus
  "#9a3e5a",  # 10 Non-Hodgkin lymphoma
  "#8b008b",  # 11 Leukaemia
  "#ffb6c1",  # 12 Kidney
  "#FFD700",  # 13 Melanoma of skin
  "#6a5acd",  # 14 Corpus uteri
  "#3C1E40",  # 15 Lip, oral cavity
  "#bc8f8f",  # 16 Cervix uteri
  "#e84393",  # 17 Thyroid
  "#2a9d8f",  # 18 Multiple myeloma
  "#4169e1",  # 19 Brain & CNS
  "#ff4500",  # 20 Ovary
  "#63ADEE",  # 21 Larynx
  "#FFFACD",  # 22 Gallbladder
  "#e22627",  # 23 Oropharynx
  "#0d4e92",  # 24 Hypopharynx
  "#e9f5f3",  # 25 Nasopharynx
  "#c07020",  # 26 Vulva
  "#3d5a80",  # 27 Salivary glands
  "#1a1a2e",  # 28 Mesothelioma
  "#4682b4",  # 29 Penis
  "#7B3F00",  # 30 Hodgkin lymphoma
  "#8fbc8f",  # 31 Vagina
  "#b31f39",  # 32 Kaposi sarcoma
  "#5F9EA0"   # 33 Testis
)

# Truncate to match the number of cancer types
colors <- colors[1:length(cause_levels)]

# Check for duplicates
if (any(duplicated(colors))) {
  warning("Duplicate colours detected! Please check.")
  print(colors[duplicated(colors)])
} else {
  cat("Colour check passed: all colours are unique.\n")
}

cat("Number of colours:", length(colors), "\n\n")

# -----------------------------------------------------------------------------
# 6. Stacked bar chart
# -----------------------------------------------------------------------------
p_result <- ggplot(
  data = result,
  mapping = aes(x = val_ratio, y = `Country label`, fill = `Cancer label`)
) +
  geom_bar(
    stat = "identity", position = "stack", width = 0.75,
    alpha = 0.9, color = "black", linewidth = 0.15
  ) +
  scale_x_continuous(
    name = "Proportion of new cases (%)",
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(limits = rev) +
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
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.25),
    axis.line = element_line(color = "black", size = 0.25),
    axis.ticks = element_line(color = "black", size = 0.25),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5)
  ) +
  guides(fill = guide_legend(keysize = 0.4))

print(p_result)

# Recommended export size: 3.5 x 9 inches

# -----------------------------------------------------------------------------
# 7. Export supplementary tables
# -----------------------------------------------------------------------------

# Table 1: proportion of new cases (%)
result_wide <- result %>%
  select(`Country label`, `Cancer label`, val_ratio) %>%
  pivot_wider(names_from = `Cancer label`, values_from = val_ratio)

result_wide <- result_wide %>%
  mutate(`Country label` = factor(`Country label`, levels = region_levels)) %>%
  arrange(`Country label`) %>%
  mutate(`Country label` = as.character(`Country label`))

result_wide <- result_wide %>%
  select(`Country label`, all_of(cause_levels[cause_levels %in% names(result_wide)]))

write.csv(result_wide, "Table_S_incidence_proportion_by_region.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Table 2: incidence counts
result_wide <- result %>%
  select(`Country label`, `Cancer label`, Number) %>%
  pivot_wider(names_from = `Cancer label`, values_from = Number)

result_wide <- result_wide %>%
  mutate(`Country label` = factor(`Country label`, levels = region_levels)) %>%
  arrange(`Country label`) %>%
  mutate(`Country label` = as.character(`Country label`))

result_wide <- result_wide %>%
  select(`Country label`, all_of(cause_levels[cause_levels %in% names(result_wide)]))

write.csv(result_wide, "Table_S_incidence_count_by_region.csv",
          row.names = FALSE, fileEncoding = "UTF-8")