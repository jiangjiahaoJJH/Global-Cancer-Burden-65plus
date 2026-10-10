# =============================================================================
# Figure 5 / Supplementary: Joinpoint regression analysis of cancer trends
# -----------------------------------------------------------------------------
# Description:
#   This script prepares input data for Joinpoint regression analysis of
#   age-standardised incidence rates (ASIR) and mortality rates (ASMR) for
#   adults aged 65+ years, using GBD 2000-2023 data. It then processes the
#   AAPC output from Joinpoint and exports a supplementary table of annual
#   percentage changes (APC) by cancer type.
#
# Inputs:
#   - GBD_65plus_ASR_incidence_APC_export.txt   (Joinpoint export for incidence)
#   - GBD_65plus_ASR_mortality_APC_export.txt   (Joinpoint export for mortality)
#   - cause_id_lookup.csv                       (maps cause_id to cause_name)
#
# Outputs:
#   - Joinpoint_input_ASIR.csv
#   - Joinpoint_input_ASMR.csv
#   - Table_S_global_cancer_incidence_APC.csv
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
setwd("path/to/working/directory")   # <- update to your local path

library(data.table)
library(dplyr)

# =============================================================================
# Part 1: Prepare Joinpoint input data (incidence and mortality)
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Load GBD data for incidence
# -----------------------------------------------------------------------------
IBD_china <- fread("GBD_65plus_ASR_incidence_APC_export.txt")

# Inspect
unique(IBD_china$measure_name)
unique(IBD_china$location_name)
unique(IBD_china$cause_name)

# Filter: incidence, age-standardised, rate, Both sexes, 2000-2023
ASIR <- subset(
  IBD_china,
  measure_name == "Incidence" &
    age_name == "Age-standardized" &
    metric_name == "Rate" &
    cause_name != "Non-melanoma skin cancer" &
    sex_name == "Both" &
    year >= 2000 & year <= 2023
)

# Calculate standard error from 95% UI
ASIR$SE <- (ASIR$upper - ASIR$lower) / (1.96 * 2)

# Keep required columns
ASIR <- ASIR[, c("cause_id", "year", "val", "SE")]

# Sort by cause_id and year (required by Joinpoint)
ASIR <- ASIR[order(ASIR$cause_id, ASIR$year), ]

# Export
write.csv(ASIR, "Joinpoint_input_ASIR.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 1.2 Load GBD data for mortality
# -----------------------------------------------------------------------------
# Re-use IBD_china but filter for deaths
ASMR <- subset(
  IBD_china,
  measure_name == "Deaths" &
    age_name == "Age-standardized" &
    metric_name == "Rate" &
    cause_name != "Non-melanoma skin cancer" &
    sex_name == "Both" &
    year >= 2000 & year <= 2023
)

# Calculate standard error
ASMR$SE <- (ASMR$upper - ASMR$lower) / (1.96 * 2)

# Keep required columns
ASMR <- ASMR[, c("cause_id", "year", "val", "SE")]

# Sort
ASMR <- ASMR[order(ASMR$cause_id, ASMR$year), ]

# Export
write.csv(ASMR, "Joinpoint_input_ASMR.csv", row.names = FALSE)

# =============================================================================
# Part 2: Process Joinpoint AAPC output and export supplementary table
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 Read AAPC output from Joinpoint (incidence)
# -----------------------------------------------------------------------------
AAPC <- read.table("GBD_65plus_ASR_incidence_APC_export.txt", header = TRUE)

# Rename columns 6-8 to val, lower, upper
names(AAPC)[6:8] <- c("val", "lower", "upper")

# Format to two decimal places
AAPC$val   <- sprintf("%0.2f", AAPC$val)
AAPC$lower <- sprintf("%0.2f", AAPC$lower)
AAPC$upper <- sprintf("%0.2f", AAPC$upper)

# Combine into a single string: "val (lower, upper)"
AAPC$AAPC <- paste0(AAPC$val, " (", AAPC$lower, ", ", AAPC$upper, ")")

# Keep cause_id and AAPC
newdata <- AAPC[, c(1, 6, 12)]

# -----------------------------------------------------------------------------
# 2.2 Merge with cause_id lookup
# -----------------------------------------------------------------------------
data1 <- fread("cause_id_lookup.csv")

merged_data <- newdata %>%
  left_join(data1, by = "cause_id")

# Set factor order for cause_name (as in the manuscript)
merged_data$cause_name <- factor(
  merged_data$cause_name,
  levels = c(
    "Lung", "Colorectal", "Prostate", "Breast", "Stomach",
    "Liver", "Bladder", "Pancreas", "Esophagus", "Non-Hodgkin lymphoma",
    "Total cancers", "Lip and oral cavity", "Nasopharynx",
    "Other pharynx", "Larynx", "Biliary tract", "Kidney",
    "Testis", "Uterus", "Cervix", "Ovary", "Leukemia", "Hodgkin lymphoma",
    "Multiple myeloma", "Brain and CNS", "NB and OPNCT", "Eye",
    "Thyroid", "Melanoma", "Bone and cartilage", "Soft tissue",
    "Mesothelioma", "Others"
  )
)

# Sort by cause_name
merged_data <- merged_data %>%
  arrange(cause_name)

# Export
write.csv(merged_data, "Table_S_global_cancer_incidence_APC.csv", row.names = FALSE)

# =============================================================================
# Notes:
# - The Joinpoint software must be run separately using the exported CSV files.
# - The AAPC output is then read back into R for formatting and table export.
# - The same workflow can be applied to mortality data by replacing the input
#   file and output table name.
# =============================================================================