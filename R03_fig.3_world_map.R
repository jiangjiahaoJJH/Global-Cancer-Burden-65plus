# =============================================================================
# Figure 3: World map of cancer burden among adults aged 65+ (2024)
# -----------------------------------------------------------------------------
# Description:
#   This script generates a choropleth world map showing the age-standardised
#   mortality rate (ASMR) or incidence rate (ASIR) for all cancers (excluding
#   non-melanoma skin cancer) among adults aged 65 years and older in 2024.
#   Data are from GLOBOCAN 2024.
#
# Input:
#   - GLOBOCAN_65plus_allcancers_deaths.csv
#     (contains columns: Country label, ASR)
#
# Output:
#   - World_map_ASMR_2024.png (saved at 5 x 8 inches)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
setwd("path/to/working/directory")   # <- update to your local path

library(reshape)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggmap)
library(maps)
library(data.table)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
a <- fread("GLOBOCAN_65plus_allcancers_deaths.csv")

# -----------------------------------------------------------------------------
# 2. Prepare ASR data
# -----------------------------------------------------------------------------
ASIR <- a[, c("Country label", "ASR")]

# Manually add Taiwan using China's value (if needed)
ASIR <- rbind(ASIR, setNames(
  data.frame("Taiwan", ASIR[ASIR$`Country label` == "China", "ASR"][1]),
  c("Country label", "ASR")
))

# -----------------------------------------------------------------------------
# 3. Prepare world map data
# -----------------------------------------------------------------------------
worldData <- map_data("world")

# Rename for joining
colnames(ASIR) <- c("region", "val")

# Join map data with ASR data
total <- full_join(worldData, ASIR, by = "region")

# Check for unmatched country names
world_countries <- unique(worldData$region)
unmatched <- ASIR %>%
  filter(!region %in% world_countries) %>%
  distinct(region)

cat("Number of unmatched regions:", nrow(unmatched), "\n")
print(unmatched)

# -----------------------------------------------------------------------------
# 4. Sort to prevent map distortion
# -----------------------------------------------------------------------------
total <- total %>% arrange(group, order)

# -----------------------------------------------------------------------------
# 5. Create quantile-based categories
# -----------------------------------------------------------------------------
val_non_na <- total$val[!is.na(total$val)]

if (length(val_non_na) > 0) {
  quantiles <- quantile(val_non_na, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
  labels <- paste0(round(quantiles[-length(quantiles)], 2), " to ", round(quantiles[-1], 2))
  total <- total %>%
    mutate(val2 = cut(val, breaks = quantiles, labels = labels, include.lowest = TRUE))
} else {
  stop("The 'val' column is entirely NA. Please check the data merge.")
}

# -----------------------------------------------------------------------------
# 6. Custom colour palette
# -----------------------------------------------------------------------------
custom_colors <- setNames(
  c("#0d4e92", "#0095c0", "#aed9ea", "#e9f5f3", "#fdeb97", "#f89c59", "#e22627", "#b31f39"),
  labels
)

# -----------------------------------------------------------------------------
# 7. Plot world map
# -----------------------------------------------------------------------------
p <- ggplot() +
  geom_polygon(data = total, aes(x = long, y = lat, group = group, fill = val2),
               colour = "black", size = 0.12) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  labs(fill = "ASMR (per 100,000)") +   # change legend title as needed
  theme(
    legend.position = c(0.05, 0.22),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  ggtitle("All cancers excluding NMSC (2024)")

print(p)

# Save figure (5 x 8 inches)
ggsave("World_map_ASMR_2024.png", plot = p, width = 5, height = 8, dpi = 300)

# =============================================================================
# Notes:
# - The manual addition of Taiwan using China's value may be adjusted based on
#   data availability.
# - Colour palette can be modified to suit journal style.
# - Legend title and plot title should be updated according to the measure
#   (ASIR vs ASMR) and year.
# =============================================================================