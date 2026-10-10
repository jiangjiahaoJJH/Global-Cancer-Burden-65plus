# =============================================================================
# Das Gupta decomposition of cancer burden changes, 2000-2023
# -----------------------------------------------------------------------------
# Description:
#   This script decomposes the change in cancer incidence (or mortality) among
#   adults aged 65+ into three components: population ageing (a_effect),
#   population growth (p_effect), and epidemiological change (r_effect).
#   The decomposition is performed for Both, Male, and Female sexes at the
#   global level, using GBD 2000-2023 data.
#
# Inputs:
#   - data/GBD_older_adults_global_cancer_summary.csv
#       (columns: location_name, sex_name, year, age_name, metric_name,
#        measure_name, cause_name, val)
#   - data/GBD_global_population_2000_2023_age65plus.csv
#       (columns: location_name, sex_name, year, age_name, val)
#
# Outputs:
#   - output/incidence_decomposition_by_sex_global_65plus.csv
#   - figures/incidence_decomposition_65plus.png
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
setwd("path/to/project/root")   # <- update to your local project root

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

data_dir   <- "data"
output_dir <- "output"
fig_dir    <- "figures"

if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(fig_dir))    dir.create(fig_dir)

# -----------------------------------------------------------------------------
# 1. Read disease burden data
# -----------------------------------------------------------------------------
disease_file <- file.path(data_dir, "GBD_older_adults_global_cancer_summary.csv")
cancer_data  <- fread(disease_file)

# Filter years and cancer type
cancer_data <- cancer_data %>%
  filter(
    year >= 2000 & year <= 2023,
    cause_name == "Total Cancers excluding Non-melanoma skin cancer"
  )

# Older age groups
age_groups <- c("65-69 years", "70-74 years", "75-79 years",
                "80-84 years", "85+ years")

cancer_older <- cancer_data %>%
  filter(age_name %in% age_groups) %>%
  mutate(age_name = factor(gsub(" years", "", age_name),
                           levels = c("65-69", "70-74", "75-79", "80-84", "85+")))

# -----------------------------------------------------------------------------
# 2. Read population data
# -----------------------------------------------------------------------------
pop_file   <- file.path(data_dir, "GBD_global_population_2000_2023_age65plus.csv")
population <- fread(pop_file)

population <- population %>%
  select(location_name, sex_name, year, age_name, val) %>%
  filter(age_name %in% age_groups, year >= 2000 & year <= 2023) %>%
  mutate(age_name = factor(gsub(" years", "", age_name),
                           levels = c("65-69", "70-74", "75-79", "80-84", "85+")))

# -----------------------------------------------------------------------------
# 3. Decomposition setup
# -----------------------------------------------------------------------------
decomposition_name <- c(
  "location_name", "sex_name", "overall_difference",
  "a_effect", "p_effect", "r_effect",
  "a_percent", "p_percent", "r_percent"
)
decomposition_data <- as.data.frame(
  matrix(nrow = 0, ncol = length(decomposition_name))
)
names(decomposition_data) <- decomposition_name

target_location <- "Global"
sex_categories  <- c("Both", "Male", "Female")

# -----------------------------------------------------------------------------
# 4. Loop over sexes
# -----------------------------------------------------------------------------
for (sex in sex_categories) {
  cat("Processing sex:", sex, "\n")
  
  # --- 4.1 Population 2000 ---
  pop_2000 <- population %>%
    filter(location_name == target_location, year == 2000, sex_name == sex) %>%
    arrange(age_name)
  
  if (nrow(pop_2000) == 0) {
    cat("Warning: no population data for", sex, "in 2000. Skipping.\n")
    next
  }
  
  P_2000        <- sum(pop_2000$val)
  pop_2000$percent <- pop_2000$val / P_2000
  
  # --- 4.2 Population 2023 ---
  pop_2023 <- population %>%
    filter(location_name == target_location, year == 2023, sex_name == sex) %>%
    arrange(age_name)
  
  if (nrow(pop_2023) == 0) {
    cat("Warning: no population data for", sex, "in 2023. Skipping.\n")
    next
  }
  
  P_2023        <- sum(pop_2023$val)
  pop_2023$percent <- pop_2023$val / P_2023
  
  # --- 4.3 Age proportions and total population ---
  a_2000 <- pop_2000$percent
  a_2023 <- pop_2023$percent
  p_2000 <- P_2000
  p_2023 <- P_2023
  
  # --- 4.4 Incidence rates ---
  has_metric <- all(c("metric_name", "measure_name") %in% names(cancer_older))
  
  if (!has_metric) {
    case_2000 <- cancer_older %>%
      filter(year == 2000, location_name == target_location, sex_name == sex) %>%
      arrange(age_name)
    case_2023 <- cancer_older %>%
      filter(year == 2023, location_name == target_location, sex_name == sex) %>%
      arrange(age_name)
  } else {
    case_2000 <- cancer_older %>%
      filter(
        year == 2000, location_name == target_location, sex_name == sex,
        metric_name == "Rate", measure_name == "Incidence",
        cause_name == "Total Cancers excluding Non-melanoma skin cancer"
      ) %>%
      arrange(age_name)
    case_2023 <- cancer_older %>%
      filter(
        year == 2023, location_name == target_location, sex_name == sex,
        metric_name == "Rate", measure_name == "Incidence",
        cause_name == "Total Cancers excluding Non-melanoma skin cancer"
      ) %>%
      arrange(age_name)
  }
  
  if (nrow(case_2000) == 0 | nrow(case_2023) == 0) {
    cat("Warning: incomplete rate data for", sex, ". Skipping.\n")
    next
  }
  
  r_2000 <- as.numeric(case_2000$val) / 1e5
  r_2023 <- as.numeric(case_2023$val) / 1e5
  
  # --- 4.5 Das Gupta effects ---
  a_effect <- round(
    (sum(a_2023 * p_2000 * r_2000) + sum(a_2023 * p_2023 * r_2023)) / 3 +
      (sum(a_2023 * p_2000 * r_2023) + sum(a_2023 * p_2023 * r_2000)) / 6 -
      (sum(a_2000 * p_2000 * r_2000) + sum(a_2000 * p_2023 * r_2023)) / 3 -
      (sum(a_2000 * p_2000 * r_2023) + sum(a_2000 * p_2023 * r_2000)) / 6,
    3
  )
  
  p_effect <- round(
    (sum(a_2000 * p_2023 * r_2000) + sum(a_2023 * p_2023 * r_2023)) / 3 +
      (sum(a_2000 * p_2023 * r_2023) + sum(a_2023 * p_2023 * r_2000)) / 6 -
      (sum(a_2000 * p_2000 * r_2000) + sum(a_2023 * p_2000 * r_2023)) / 3 -
      (sum(a_2000 * p_2000 * r_2023) + sum(a_2023 * p_2000 * r_2000)) / 6,
    3
  )
  
  r_effect <- round(
    (sum(a_2000 * p_2000 * r_2023) + sum(a_2023 * p_2023 * r_2023)) / 3 +
      (sum(a_2000 * p_2023 * r_2023) + sum(a_2023 * p_2000 * r_2023)) / 6 -
      (sum(a_2000 * p_2000 * r_2000) + sum(a_2023 * p_2023 * r_2000)) / 3 -
      (sum(a_2000 * p_2023 * r_2000) + sum(a_2023 * p_2000 * r_2000)) / 6,
    3
  )
  
  overall_difference <- round(a_effect + p_effect + r_effect, 2)
  
  # --- 4.6 Percentage contributions ---
  if (overall_difference != 0) {
    a_percent <- round(a_effect / overall_difference * 100, 2)
    p_percent <- round(p_effect / overall_difference * 100, 2)
    r_percent <- round(r_effect / overall_difference * 100, 2)
  } else {
    a_percent <- p_percent <- r_percent <- 0
  }
  
  # --- 4.7 Append results (preserve numeric types) ---
  temp <- data.frame(
    location_name      = target_location,
    sex_name           = sex,
    overall_difference = overall_difference,
    a_effect           = a_effect,
    p_effect           = p_effect,
    r_effect           = r_effect,
    a_percent          = a_percent,
    p_percent          = p_percent,
    r_percent          = r_percent,
    stringsAsFactors   = FALSE
  )
  decomposition_data <- rbind(decomposition_data, temp)
  
  cat("Finished sex:", sex, "\n")
}

# -----------------------------------------------------------------------------
# 5. Total number change (for validation)
# -----------------------------------------------------------------------------
burden_data <- fread(disease_file)

has_metric_burden <- all(c("metric_name", "measure_name") %in% names(burden_data))

if (!has_metric_burden) {
  num_2000 <- burden_data %>%
    filter(
      age_name == "65+ years", location_name == target_location,
      cause_name == "Total Cancers excluding Non-melanoma skin cancer",
      year == 2000
    ) %>%
    select(location_name, sex_name, val) %>%
    rename(val_2000 = val)
  
  num_2023 <- burden_data %>%
    filter(
      age_name == "65+ years", location_name == target_location,
      cause_name == "Total Cancers excluding Non-melanoma skin cancer",
      year == 2023
    ) %>%
    select(location_name, sex_name, val) %>%
    rename(val_2023 = val)
} else {
  num_2000 <- burden_data %>%
    filter(
      age_name == "65+ years", location_name == target_location,
      cause_name == "Total Cancers excluding Non-melanoma skin cancer",
      metric_name == "Number", measure_name == "Incidence",
      year == 2000
    ) %>%
    select(location_name, sex_name, val) %>%
    rename(val_2000 = val)
  
  num_2023 <- burden_data %>%
    filter(
      age_name == "65+ years", location_name == target_location,
      cause_name == "Total Cancers excluding Non-melanoma skin cancer",
      metric_name == "Number", measure_name == "Incidence",
      year == 2023
    ) %>%
    select(location_name, sex_name, val) %>%
    rename(val_2023 = val)
}

decomposition_data <- decomposition_data %>%
  left_join(num_2000, by = c("location_name", "sex_name")) %>%
  left_join(num_2023, by = c("location_name", "sex_name")) %>%
  mutate(diff1 = as.numeric(val_2023) - as.numeric(val_2000))

# Convert all numeric columns safely
num_cols <- c(
  "overall_difference", "a_effect", "p_effect", "r_effect",
  "a_percent", "p_percent", "r_percent",
  "val_2000", "val_2023", "diff1"
)
decomposition_data[num_cols] <- lapply(decomposition_data[num_cols], as.numeric)

print(decomposition_data)

# Save results
write.csv(
  decomposition_data,
  file.path(output_dir, "incidence_decomposition_by_sex_global_65plus.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Plot
# -----------------------------------------------------------------------------
decomposition_data <- decomposition_data %>%
  rename(
    `Overall difference`    = overall_difference,
    Aging                   = a_effect,
    Population              = p_effect,
    `Epidemiological change` = r_effect
  )

decomposition_plot <- decomposition_data[, 1:6] %>%
  pivot_longer(
    cols = c(Aging, Population, `Epidemiological change`),
    names_to  = "varname",
    values_to = "value"
  ) %>%
  mutate(
    value    = as.numeric(value),
    sex_name = factor(sex_name, levels = c("Both", "Male", "Female"))
  )

my_colors <- c("#C90B0B", "#0A7BC1", "#b0a090")

p2 <- ggplot(decomposition_plot, aes(x = sex_name, y = value, fill = varname)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(
    values = my_colors,
    name   = "Effect",
    labels = c("Aging", "Population", "Epidemiological change")
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 8, color = "black"),
    axis.text.x  = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.text  = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line    = element_line(color = "black", linewidth = 0.25),
    axis.ticks   = element_line(color = "black", linewidth = 0.25),
    plot.title   = element_text(size = 12, color = "black", hjust = 0.5)
  ) +
  labs(x = "", y = "Incidence") +
  ggtitle("Incidence of total cancer (65+ years)")

decomposition_data$sex_name <- factor(decomposition_data$sex_name,
                                      levels = c("Both", "Male", "Female"))

plot2 <- p2 +
  geom_point(
    data = decomposition_data,
    mapping = aes(x = sex_name, y = `Overall difference`),
    fill = "black", color = "black", size = 2
  )

print(plot2)

ggsave(
  file.path(fig_dir, "incidence_decomposition_65plus.png"),
  plot2, width = 4, height = 5, dpi = 300
)

# -----------------------------------------------------------------------------
# 7. Session info
# -----------------------------------------------------------------------------
sessionInfo()