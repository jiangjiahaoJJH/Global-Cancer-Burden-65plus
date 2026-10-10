# =============================================================================
# Figure 4: Association between HDI and cancer mortality among adults aged 65+
# -----------------------------------------------------------------------------
# Description:
#   This script examines the correlation between the Human Development Index
#   (HDI) and age-standardised mortality rates (ASMR) for all cancers
#   (excluding non-melanoma skin cancer) among adults aged 65 years and older.
#   Data are from GLOBOCAN 2024 and the UNDP HDI 2022.
#   A Spearman correlation test is performed, and a scatter plot with a LOESS
#   smooth curve is generated.
#
# Input:
#   - GLOBOCAN_65plus_ASR_mortality_186countries.csv
#     (columns: Cancer label, Country label, Cancer code, Country code,
#      Sex, type_label, age_label, ASR)
#   - HDI_2022.csv
#     (columns: Country, HDI)
#
# Outputs:
#   - Total_cancers_mortality_HDI_correlation.csv
#   - Figure4_HDI_ASMR.png (recommended size: 3.5 x 4 inches)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
setwd("path/to/working/directory")   # <- update to your local path

library(reshape)
library(ggplot2)
library(ggrepel)
library(data.table)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
IBD_china <- fread("GLOBOCAN_65plus_ASR_mortality_186countries.csv",
                   header = TRUE, stringsAsFactors = FALSE)

# Rename columns for convenience
colnames(IBD_china)[colnames(IBD_china) == "Cancer label"] <- "Cancer_label"
colnames(IBD_china)[colnames(IBD_china) == "Country label"] <- "Country"
colnames(IBD_china)[colnames(IBD_china) == "Cancer code"] <- "Cancer_code"
colnames(IBD_china)[colnames(IBD_china) == "Country code"] <- "Country_code"

# Inspect data
head(IBD_china)
unique(IBD_china$Cancer_label)
unique(IBD_china$Sex)
unique(IBD_china$type_label)
unique(IBD_china$age_label)

# -----------------------------------------------------------------------------
# 2. Extract all-cancer mortality data (Both sexes, age-standardised)
# -----------------------------------------------------------------------------
data1 <- subset(
  IBD_china,
  age_label == "Age-standardized" &
    type_label == "Mortality" &
    Sex == "Both" &
    Cancer_label == "All cancers excl. non-melanoma skin cancer"
)

# Load HDI data
HDI_data <- read.csv("HDI_2022.csv", header = TRUE)

# -----------------------------------------------------------------------------
# 3. Merge ASMR with HDI
# -----------------------------------------------------------------------------
IBD_ASMR <- data1
IBD_ASMR <- IBD_ASMR[, .(Country, ASR)]

# Fix country name encoding
name_map <- c(
  "C?te d'Ivoire" = "Côte d'Ivoire",
  "Cura?ao" = "Curaçao",
  "France, La R\xa8\xa6union" = "France, La Réunion",
  "T\xa8\xb9rkiye" = "Türkiye"
)

for (old_name in names(name_map)) {
  IBD_ASMR$Country[IBD_ASMR$Country == old_name] <- name_map[old_name]
}

# Merge with HDI
IBD_ASMR_HDI <- merge(IBD_ASMR, HDI_data, by = "Country")

# Check unmatched countries
not_matched <- IBD_ASMR[!Country %in% HDI_data$Country, ]
print(not_matched)

# Rename ASR to val
names(IBD_ASMR_HDI)[names(IBD_ASMR_HDI) == "ASR"] <- "val"

# Check structure
str(IBD_ASMR_HDI$val)
str(IBD_ASMR_HDI$HDI)

# -----------------------------------------------------------------------------
# 4. Spearman correlation test
# -----------------------------------------------------------------------------
result <- cor.test(IBD_ASMR_HDI$val, IBD_ASMR_HDI$HDI, method = "spearman")

correlation_coefficient <- round(result$estimate, 3)
p_value <- ifelse(result$p.value < 0.001, "<0.001", round(result$p.value, 3))

# -----------------------------------------------------------------------------
# 5. Prepare annotation
# -----------------------------------------------------------------------------
y_range <- c(0, 1600)
y_text1 <- y_range[1] + 0.9 * diff(y_range)
y_text2 <- y_range[1] + 0.7 * diff(y_range)

anno_df <- data.frame(
  x = 0.3,
  y = c(y_text1, y_text2),
  lab = c(
    paste0("rs: ", correlation_coefficient),
    paste0("p: ", p_value)
  )
)

# -----------------------------------------------------------------------------
# 6. Scatter plot with LOESS smooth
# -----------------------------------------------------------------------------
P <- ggplot(IBD_ASMR_HDI, aes(x = HDI, y = val)) +
  geom_point(color = "#0A7BC1", size = 1.5) +
  geom_smooth(colour = "#C90B0B", method = "loess", se = FALSE,
              span = 0.9, size = 1) +
  labs(
    x = "HDI",
    y = "ASMR (per 100,000)",
    title = "All cancers excluded NMSC"
  ) +
  scale_x_continuous(
    limits = c(0.25, 1),
    breaks = seq(0.3, 1.0, by = 0.1),
    labels = seq(0.3, 1.0, by = 0.1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 1600), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_line(color = "black", size = 0.25),
    axis.ticks = element_line(color = "black", size = 0.25),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5)
  ) +
  geom_text(data = anno_df, aes(x = x, y = y, label = lab),
            size = 4, hjust = 0)

print(P)

# Recommended export size: 3.5 x 4 inches
ggsave("Figure4_HDI_ASMR.png", plot = P, width = 3.5, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# 7. Export correlation results
# -----------------------------------------------------------------------------
output <- data.frame(
  correlation_coefficient = result$estimate,
  p_value = result$p.value
)
write.csv(output, "Total_cancers_mortality_HDI_correlation.csv", row.names = FALSE)

# =============================================================================
# Notes:
# - The same script can be adapted for other cancer types by changing the
#   Cancer_label filter (e.g., Esophageal cancer, Gastric cancer,
#   Colorectal cancer, Liver cancer, Biliary tract cancer, Pancreatic cancer).
# =============================================================================