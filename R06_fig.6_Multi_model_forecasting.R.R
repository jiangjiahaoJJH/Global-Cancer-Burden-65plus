# =============================================================================
# Forecasting cancer burden among adults aged 65+ (GBD 1990-2023)
# -----------------------------------------------------------------------------
# Workflow:
#   1. Load GBD data and construct three independent series
#      (Male, Female, Both) for incidence and deaths.
#   2. Evaluate five candidate time-series models via internal time-series
#      holdout validation (train 1990-2018, test 2019-2023).
#   3. Rank models using a normalised composite score and select the best.
#   4. Refit the selected model on the full period (1990-2023) and project
#      to 2050 with prediction intervals.
#   5. Rolling-origin validation and sensitivity analyses
#      (training window, COVID-19, GAM smoothing parameter k).
#
# Outputs (written to the working directory):
#   - FiveModel_holdout_by_series.csv
#   - FiveModel_composite_score.csv
#   - BestModel_main_forecast_Incidence.csv
#   - BestModel_main_forecast_Deaths.csv
#   - Supplementary_Table_S18_wide_CI.csv
#   - Supplementary_Table_S2_best_model_params.csv
#   - Supplementary_Table_S21_rolling_origin.csv
#   - Supplementary_Table_S21_summary.csv
#   - Supplementary_Table_S22_training_window.csv
#   - Supplementary_Table_S22_window_range.csv
#   - Supplementary_Table_S23_covid_sensitivity.csv
#   - Supplementary_Table_S23_covid_pct.csv
#   - Supplementary_Table_S24_k_sensitivity.csv
#   - Supplementary_Table_S24_k_range.csv
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
setwd("D:/path/to/working/directory")   # <- set your own path

library(forecast)
library(mgcv)
library(data.table)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# 1. Load GBD data
# -----------------------------------------------------------------------------
IBD_china <- fread("GBD-input-65plus-global.csv")
IBD_china <- IBD_china[year >= 1990 & year <= 2023, ]

# Extract a single series (one sex, one measure)
get_series_df <- function(df, sex, measure = "Incidence") {
  df %>%
    filter(
      age_name      == "65+ years",
      sex_name      == sex,
      location_name == "Global",
      metric_name   == "Number",
      measure_name  == measure,
      cause_name    == "Total Cancers excluding Non-melanoma skin cancer"
    ) %>%
    arrange(year) %>%
    select(year, val)
}

# Build three independent series: Male, Female, Both
get_three_series <- function(measure) {
  df_m <- get_series_df(IBD_china, "Male",   measure)
  df_f <- get_series_df(IBD_china, "Female", measure)
  df_b <- df_m
  df_b$val <- df_m$val + df_f$val
  list(Male = df_m, Female = df_f, Both = df_b)
}

# -----------------------------------------------------------------------------
# 2. Candidate models
# -----------------------------------------------------------------------------
# Each fit_* function returns list(pred = numeric vector, fit = fitted model,
#                                  ok = logical flag)

fit_arima <- function(train_ts, h) {
  use_log <- all(train_ts > 0)
  fit <- tryCatch(
    auto.arima(train_ts,
               max.p = 3, max.d = 1, max.q = 3,
               lambda = if (use_log) 0 else NULL,
               stepwise = FALSE, approximation = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(pred = NA, ok = FALSE))
  fc <- forecast(fit, h = h)
  list(pred = as.numeric(fc$mean), fit = fit, ok = TRUE)
}

fit_poisson <- function(train_df, h) {
  fit <- tryCatch(
    glm(val ~ year, family = poisson(link = "log"), data = train_df),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(pred = NA, ok = FALSE))
  nd <- data.frame(year = (max(train_df$year) + 1):(max(train_df$year) + h))
  pred <- tryCatch(predict(fit, newdata = nd, type = "response"),
                   error = function(e) NA)
  list(pred = as.numeric(pred), fit = fit, ok = TRUE)
}

fit_gam <- function(train_df, h) {
  fit <- tryCatch(
    gam(val ~ s(year, k = 10), family = gaussian(), method = "REML",
        data = train_df),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(pred = NA, ok = FALSE))
  nd <- data.frame(year = (max(train_df$year) + 1):(max(train_df$year) + h))
  pred <- tryCatch(predict(fit, newdata = nd, type = "response"),
                   error = function(e) NA)
  list(pred = as.numeric(pred), fit = fit, ok = TRUE)
}

fit_loglinear <- function(train_df, h) {
  fit <- tryCatch(lm(val ~ log(year), data = train_df),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(pred = NA, ok = FALSE))
  nd <- data.frame(year = (max(train_df$year) + 1):(max(train_df$year) + h))
  pred <- tryCatch(predict(fit, newdata = nd), error = function(e) NA)
  list(pred = as.numeric(pred), fit = fit, ok = TRUE)
}

fit_linear <- function(train_df, h) {
  fit <- tryCatch(lm(val ~ year, data = train_df),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(pred = NA, ok = FALSE))
  nd <- data.frame(year = (max(train_df$year) + 1):(max(train_df$year) + h))
  pred <- tryCatch(predict(fit, newdata = nd), error = function(e) NA)
  list(pred = as.numeric(pred), fit = fit, ok = TRUE)
}

# Unified dispatch
predict_best <- function(train_df, h, model) {
  train_ts <- ts(train_df$val, start = min(train_df$year), frequency = 1)
  switch(model,
         ARIMA     = fit_arima(train_ts, h),
         Poisson   = fit_poisson(train_df, h),
         GAM       = fit_gam(train_df, h),
         LogLinear = fit_loglinear(train_df, h),
         Linear    = fit_linear(train_df, h)
  )
}

# -----------------------------------------------------------------------------
# 3. Internal time-series holdout validation and composite score
# -----------------------------------------------------------------------------
compute_metrics <- function(obs, pred) {
  rmse <- sqrt(mean((obs - pred)^2))
  mae  <- mean(abs(obs - pred))
  mape <- mean(abs((obs - pred) / obs)) * 100
  data.frame(RMSE = rmse, MAE = mae, MAPE = mape)
}

run_holdout <- function(measure, train_end = 2018,
                        test_start = 2019, test_end = 2023) {
  results <- list()
  series_list <- get_three_series(measure)
  
  for (sex in c("Male", "Female", "Both")) {
    df       <- series_list[[sex]]
    train_df <- df %>% filter(year <= train_end)
    test_df  <- df %>% filter(year >= test_start, year <= test_end)
    if (nrow(test_df) == 0) next
    
    h        <- nrow(test_df)
    obs      <- test_df$val
    train_ts <- ts(train_df$val, start = min(train_df$year), frequency = 1)
    
    models <- list(
      ARIMA     = fit_arima(train_ts, h),
      Poisson   = fit_poisson(train_df, h),
      GAM       = fit_gam(train_df, h),
      LogLinear = fit_loglinear(train_df, h),
      Linear    = fit_linear(train_df, h)
    )
    
    for (m in names(models)) {
      res <- models[[m]]
      if (!isTRUE(res$ok) || all(is.na(res$pred))) {
        results[[length(results) + 1]] <- data.frame(
          Measure = measure, Sex = sex, Model = m,
          RMSE = NA, MAE = NA, MAPE = NA, Fail = 1
        )
      } else {
        met <- compute_metrics(obs, res$pred)
        results[[length(results) + 1]] <- data.frame(
          Measure = measure, Sex = sex, Model = m,
          met, Fail = 0
        )
      }
    }
  }
  do.call(rbind, results)
}

val_inc   <- run_holdout("Incidence")
val_death <- run_holdout("Deaths")
val_all   <- rbind(val_inc, val_death)

fwrite(val_all, "FiveModel_holdout_by_series.csv")
print(val_all)

# Min-max normalisation, with a safeguard for constant vectors
normalize <- function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (rng == 0) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

model_summary <- val_all %>%
  group_by(Model) %>%
  summarise(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    SD_RMSE   = sd(RMSE, na.rm = TRUE),
    Mean_MAPE = mean(MAPE, na.rm = TRUE),
    FailRate  = mean(Fail, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    nRMSE     = normalize(Mean_RMSE),
    nSD       = normalize(SD_RMSE),
    nMAPE     = normalize(Mean_MAPE),
    nFail     = normalize(FailRate),
    Composite = 0.40 * nRMSE + 0.25 * nSD + 0.20 * nMAPE + 0.15 * nFail
  ) %>%
  arrange(Composite)

print(model_summary)
fwrite(model_summary, "FiveModel_composite_score.csv")

best_model <- model_summary$Model[1]
cat("\nBest model:", best_model, "\n")

# -----------------------------------------------------------------------------
# 4. Prediction intervals
# -----------------------------------------------------------------------------
# - ARIMA: intervals from forecast()
# - GAM: parameter uncertainty (se.fit) + residual uncertainty
#        + dampened growth factor + minimum width floor
# - Other models: residual uncertainty + dampened growth factor
# -----------------------------------------------------------------------------
get_pi <- function(res, train_df, h, model) {
  years <- (max(train_df$year) + 1):(max(train_df$year) + h)
  nd    <- data.frame(year = years)
  
  # ARIMA
  if (model == "ARIMA") {
    fc <- forecast(res$fit, h = h, level = c(80, 95))
    return(list(lo80 = as.numeric(fc$lower[, "80%"]),
                hi80 = as.numeric(fc$upper[, "80%"]),
                lo95 = as.numeric(fc$lower[, "95%"]),
                hi95 = as.numeric(fc$upper[, "95%"])))
  }
  
  # Dampened growth factor
  D      <- 0.11
  growth <- sqrt(1 + (0:(h - 1)) * D)
  z80    <- qnorm(0.90)
  z95    <- qnorm(0.975)
  
  if (model == "GAM") {
    pred_obj  <- predict(res$fit, newdata = nd, type = "response", se.fit = TRUE)
    base_pred <- as.numeric(pred_obj$fit)
    se_param  <- as.numeric(pred_obj$se.fit)
    se_resid  <- summary(res$fit)$scale^0.5
    se_total  <- sqrt(se_param^2 + se_resid^2)
  } else {
    base_pred <- as.numeric(predict(res$fit, newdata = nd, type = "response"))
    se_resid  <- sd(residuals(res$fit))
    se_total  <- rep(se_resid, h)
  }
  
  se_final <- se_total * growth
  
  # Minimum width floor: at least 5% of the first-year point estimate
  min_w    <- 0.05 * base_pred[1]
  se_final <- pmax(se_final, min_w / z95)
  
  list(lo80 = base_pred - z80 * se_final,
       hi80 = base_pred + z80 * se_final,
       lo95 = base_pred - z95 * se_final,
       hi95 = base_pred + z95 * se_final)
}

# -----------------------------------------------------------------------------
# 5. Final projections (2024-2050)
# -----------------------------------------------------------------------------
# Male, Female, and Both series are modelled independently.
# Their projections are not constrained to satisfy Both = Male + Female.
# -----------------------------------------------------------------------------
run_main_forecast <- function(measure, model, h = 27) {
  series_list <- get_three_series(measure)
  out_list    <- list()
  
  for (sex in c("Male", "Female", "Both")) {
    df  <- series_list[[sex]]
    res <- predict_best(df, h, model)
    if (!isTRUE(res$ok)) stop(paste(sex, "forecast failed"))
    
    pred  <- res$pred
    pi    <- get_pi(res, df, h, model)
    years <- (max(df$year) + 1):(max(df$year) + h)
    
    out_list[[length(out_list) + 1]] <- data.frame(
      Measure = measure, Sex = sex, Year = years,
      Value = pred,
      lo80 = pi$lo80, hi80 = pi$hi80,
      lo95 = pi$lo95, hi95 = pi$hi95,
      Type = "Forecast"
    )
    
    out_list[[length(out_list) + 1]] <- data.frame(
      Measure = measure, Sex = sex, Year = df$year,
      Value = df$val,
      lo80 = NA, hi80 = NA, lo95 = NA, hi95 = NA,
      Type = "Actual"
    )
  }
  bind_rows(out_list)
}

main_inc   <- run_main_forecast("Incidence", best_model)
main_death <- run_main_forecast("Deaths",   best_model)

fwrite(main_inc,   "BestModel_main_forecast_Incidence.csv")
fwrite(main_death, "BestModel_main_forecast_Deaths.csv")

# -----------------------------------------------------------------------------
# 6. Supplementary Table S18: wide format with 95% CI (integers)
# -----------------------------------------------------------------------------
s18_long <- bind_rows(main_inc, main_death) %>%
  filter(Type == "Forecast") %>%
  select(Measure, Sex, Year, Value, lo95, hi95) %>%
  mutate(
    Value = round(Value, 0),
    lo95  = round(lo95, 0),
    hi95  = round(hi95, 0),
    Cell  = paste0(Value, " (", lo95, ", ", hi95, ")")
  )

s18_wide <- s18_long %>%
  select(Year, Measure, Sex, Cell) %>%
  pivot_wider(
    id_cols     = Year,
    names_from  = c(Measure, Sex),
    values_from = Cell
  )

setnames(s18_wide,
         old = c("Incidence_Both", "Incidence_Male", "Incidence_Female",
                 "Deaths_Both", "Deaths_Male", "Deaths_Female"),
         new = c("Incidence Both", "Incidence Male", "Incidence Female",
                 "Death Both", "Death Male", "Death Female"),
         skip_absent = TRUE)

s18_wide <- s18_wide %>%
  select(Year,
         `Incidence Both`, `Incidence Male`, `Incidence Female`,
         `Death Both`, `Death Male`, `Death Female`)

print(s18_wide, n = Inf)
fwrite(s18_wide, "Supplementary_Table_S18_wide_CI.csv")

# Check interval width in 2050
s18_check <- bind_rows(main_inc, main_death) %>%
  filter(Type == "Forecast", Year == 2050) %>%
  mutate(
    width95 = hi95 - lo95,
    pct95   = 100 * width95 / Value
  )
cat("\n2050 95% interval width (% of point estimate):\n")
print(s18_check %>% select(Measure, Sex, Value, lo95, hi95, pct95))

# -----------------------------------------------------------------------------
# 7. Supplementary Table S2: model parameters of the selected model
# -----------------------------------------------------------------------------
if (best_model == "ARIMA") {
  params_list <- list()
  for (measure in c("Incidence", "Deaths")) {
    series_list <- get_three_series(measure)
    for (sex in c("Male", "Female", "Both")) {
      df      <- series_list[[sex]]
      tsx     <- ts(df$val, start = min(df$year), frequency = 1)
      use_log <- all(tsx > 0)
      fit <- auto.arima(tsx, max.p = 3, max.d = 1, max.q = 3,
                        lambda = if (use_log) 0 else NULL,
                        stepwise = FALSE, approximation = FALSE)
      n    <- length(tsx)
      k    <- length(fit$coef) + 1
      aicc <- -2 * fit$loglik + 2 * k + 2 * k * (k + 1) / (n - k - 1)
      lb   <- Box.test(residuals(fit), lag = 8, type = "Ljung-Box",
                       fitdf = length(fit$coef))
      params_list[[length(params_list) + 1]] <- data.frame(
        Measure = measure, Sex = sex,
        Order = paste(arimaorder(fit), collapse = ","),
        AICc = round(aicc, 2),
        LjungBox_Q = round(lb$statistic, 2),
        LjungBox_df = length(fit$coef),
        LjungBox_p = round(lb$p.value, 4),
        LogTransform = ifelse(use_log, "Yes", "No")
      )
    }
  }
  s2 <- do.call(rbind, params_list)
  fwrite(s2, "Supplementary_Table_S2_best_model_params.csv")
  print(s2)
}

if (best_model == "GAM") {
  params_list <- list()
  for (measure in c("Incidence", "Deaths")) {
    series_list <- get_three_series(measure)
    for (sex in c("Male", "Female", "Both")) {
      df  <- series_list[[sex]]
      fit <- gam(val ~ s(year, k = 10), family = gaussian(),
                 method = "REML", data = df)
      s   <- summary(fit)
      params_list[[length(params_list) + 1]] <- data.frame(
        Measure = measure, Sex = sex,
        Edf     = round(s$s.table[1, "edf"], 2),
        F       = round(s$s.table[1, "F"], 2),
        P       = signif(s$s.table[1, "p-value"], 4),
        DevExpl = round(s$dev.expl * 100, 2)
      )
    }
  }
  s2 <- do.call(rbind, params_list)
  fwrite(s2, "Supplementary_Table_S2_best_model_params.csv")
  print(s2)
}

# =============================================================================
# 8. Rolling-origin validation
# =============================================================================
# Training origins: 2005, 2010, 2015, 2018; each followed by 5 test years.
# All five models are evaluated at every origin.
# =============================================================================
run_rolling_origin <- function(measure,
                               origins = c(2005, 2010, 2015, 2018),
                               horizon = 5) {
  results     <- list()
  series_list <- get_three_series(measure)
  
  for (origin in origins) {
    for (sex in c("Male", "Female", "Both")) {
      df       <- series_list[[sex]]
      train_df <- df %>% filter(year <= origin)
      test_df  <- df %>% filter(year > origin, year <= origin + horizon)
      if (nrow(test_df) == 0) next
      
      h        <- nrow(test_df)
      obs      <- test_df$val
      train_ts <- ts(train_df$val, start = min(train_df$year), frequency = 1)
      scale    <- mean(abs(diff(train_df$val)))  # naive scale for MASE
      
      models <- list(
        ARIMA     = fit_arima(train_ts, h),
        Poisson   = fit_poisson(train_df, h),
        GAM       = fit_gam(train_df, h),
        LogLinear = fit_loglinear(train_df, h),
        Linear    = fit_linear(train_df, h)
      )
      
      for (m in names(models)) {
        res <- models[[m]]
        if (!isTRUE(res$ok) || all(is.na(res$pred))) {
          results[[length(results) + 1]] <- data.frame(
            Origin = origin, Measure = measure, Sex = sex, Model = m,
            RMSE = NA, MAE = NA, MAPE = NA, MASE = NA, Fail = 1
          )
        } else {
          rmse <- sqrt(mean((obs - res$pred)^2))
          mae  <- mean(abs(obs - res$pred))
          mape <- mean(abs((obs - res$pred) / obs)) * 100
          mase <- mae / scale
          results[[length(results) + 1]] <- data.frame(
            Origin = origin, Measure = measure, Sex = sex, Model = m,
            RMSE = rmse, MAE = mae, MAPE = mape, MASE = mase, Fail = 0
          )
        }
      }
    }
  }
  do.call(rbind, results)
}

rolling_inc   <- run_rolling_origin("Incidence")
rolling_death <- run_rolling_origin("Deaths")
rolling_all   <- rbind(rolling_inc, rolling_death)

print(rolling_all)
fwrite(rolling_all, "Supplementary_Table_S21_rolling_origin.csv")

rolling_summary <- rolling_all %>%
  group_by(Model, Origin) %>%
  summarise(
    Mean_MASE = mean(MASE, na.rm = TRUE),
    Mean_MAPE = mean(MAPE, na.rm = TRUE),
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(Model, Origin)

print(rolling_summary)
fwrite(rolling_summary, "Supplementary_Table_S21_summary.csv")

# =============================================================================
# 9. Sensitivity analysis 1: training window
# =============================================================================
# Train GAM on 1990-2015, 1990-2018, 1990-2023; project to 2050.
# The forecast horizon h is computed per window so that every scenario
# produces a 2050 projection.
# =============================================================================
run_training_window_sensitivity <- function(measure,
                                            train_ends = c(2015, 2018, 2023)) {
  results     <- list()
  series_list <- get_three_series(measure)
  
  for (train_end in train_ends) {
    for (sex in c("Male", "Female", "Both")) {
      df       <- series_list[[sex]]
      train_df <- df %>% filter(year <= train_end)
      h        <- 2050 - max(train_df$year)
      res      <- fit_gam(train_df, h)
      if (!isTRUE(res$ok)) next
      
      years     <- (max(train_df$year) + 1):(max(train_df$year) + h)
      pred_2050 <- res$pred[which(years == 2050)]
      
      results[[length(results) + 1]] <- data.frame(
        Measure = measure, Sex = sex, TrainEnd = train_end,
        Pred2050 = round(pred_2050, 0)
      )
    }
  }
  do.call(rbind, results)
}

sens_window <- rbind(
  run_training_window_sensitivity("Incidence"),
  run_training_window_sensitivity("Deaths")
)
print(sens_window)
fwrite(sens_window, "Supplementary_Table_S22_training_window.csv")

sens_window_wide <- sens_window %>%
  pivot_wider(
    id_cols     = c(Measure, Sex),
    names_from  = TrainEnd,
    values_from = Pred2050,
    names_prefix = "TrainEnd_"
  ) %>%
  mutate(
    Range    = apply(select(., starts_with("TrainEnd_")), 1,
                     function(x) max(x) - min(x)),
    PctRange = 100 * Range / TrainEnd_2023
  )

print(sens_window_wide)
fwrite(sens_window_wide, "Supplementary_Table_S22_window_range.csv")

# =============================================================================
# 10. Sensitivity analysis 2: COVID-19
# =============================================================================
# Scenarios:
#   Base           : full training period (1990-2023)
#   Truncate2019   : training ends in 2019 (forecast 2020-2050)
#   Exclude2020_21 : drop 2020-2021 from training (forecast 2024-2050)
# =============================================================================
run_covid_sensitivity <- function(measure) {
  results     <- list()
  series_list <- get_three_series(measure)
  
  for (sex in c("Male", "Female", "Both")) {
    df         <- series_list[[sex]]
    train_base <- df
    train_a    <- df %>% filter(year <= 2019)
    train_b    <- df %>% filter(!year %in% c(2020, 2021))
    
    for (tag in c("Base", "Truncate2019", "Exclude2020_21")) {
      train_df <- switch(tag,
                         Base           = train_base,
                         Truncate2019   = train_a,
                         Exclude2020_21 = train_b)
      
      h   <- 2050 - max(train_df$year)
      res <- fit_gam(train_df, h)
      if (!isTRUE(res$ok)) next
      
      years     <- (max(train_df$year) + 1):(max(train_df$year) + h)
      pred_2050 <- res$pred[which(years == 2050)]
      
      results[[length(results) + 1]] <- data.frame(
        Measure = measure, Sex = sex, Scenario = tag,
        Pred2050 = round(pred_2050, 0)
      )
    }
  }
  do.call(rbind, results)
}

sens_covid <- rbind(
  run_covid_sensitivity("Incidence"),
  run_covid_sensitivity("Deaths")
)
print(sens_covid)
fwrite(sens_covid, "Supplementary_Table_S23_covid_sensitivity.csv")

sens_covid_wide <- sens_covid %>%
  pivot_wider(
    id_cols     = c(Measure, Sex),
    names_from  = Scenario,
    values_from = Pred2050
  ) %>%
  mutate(
    PctChange_Truncate = 100 * (Truncate2019   - Base) / Base,
    PctChange_Exclude  = 100 * (Exclude2020_21 - Base) / Base
  )

print(sens_covid_wide)
fwrite(sens_covid_wide, "Supplementary_Table_S23_covid_pct.csv")

# =============================================================================
# 11. Sensitivity analysis 3: GAM smoothing parameter k
# =============================================================================
run_k_sensitivity <- function(measure, k_values = c(5, 10, 15), h = 27) {
  results     <- list()
  series_list <- get_three_series(measure)
  
  for (sex in c("Male", "Female", "Both")) {
    df <- series_list[[sex]]
    for (k in k_values) {
      fit <- tryCatch(
        gam(val ~ s(year, k = k), family = gaussian(), method = "REML",
            data = df),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      nd   <- data.frame(year = (max(df$year) + 1):(max(df$year) + h))
      pred <- tryCatch(predict(fit, newdata = nd, type = "response"),
                       error = function(e) NA)
      pred_2050 <- pred[length(pred)]
      
      results[[length(results) + 1]] <- data.frame(
        Measure = measure, Sex = sex, k = k,
        Pred2050 = round(pred_2050, 0)
      )
    }
  }
  do.call(rbind, results)
}

sens_k <- rbind(
  run_k_sensitivity("Incidence"),
  run_k_sensitivity("Deaths")
)
print(sens_k)
fwrite(sens_k, "Supplementary_Table_S24_k_sensitivity.csv")

sens_k_wide <- sens_k %>%
  pivot_wider(
    id_cols     = c(Measure, Sex),
    names_from  = k,
    values_from = Pred2050,
    names_prefix = "k"
  ) %>%
  mutate(
    Range    = apply(select(., starts_with("k")), 1,
                     function(x) max(x) - min(x)),
    PctRange = 100 * Range / k10
  )

print(sens_k_wide)
fwrite(sens_k_wide, "Supplementary_Table_S24_k_range.csv")

cat("\nAll analyses completed.\n")