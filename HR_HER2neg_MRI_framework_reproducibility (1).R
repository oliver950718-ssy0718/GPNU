suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rsample)
  library(glmnet)
  library(xgboost)
  library(pROC)
  library(survival)
  library(survminer)
  library(boot)
  library(broom)
  library(DescTools)
})

# -----------------------------------------------------------------------------
# 2. Assemble radiomics and externally extracted DL features
# -----------------------------------------------------------------------------

radiomics <- read_required(path$radiomics)
dl_features <- read_required(path$dl)

feature_data <- radiomics %>%
  inner_join(dl_features, by = "patient_id")


rad_cols <- names(feature_data)[str_detect(names(feature_data), "^RAD_")]
dl_cols <- names(feature_data)[str_detect(names(feature_data), "^DL_")]

# -----------------------------------------------------------------------------
# 3. Generic feature filtering and scaling
# -----------------------------------------------------------------------------

remove_near_zero_variance <- function(x) {
  keep <- vapply(x, function(z) {
    z <- z[is.finite(z)]
    length(unique(z)) > 1 && stats::sd(z) > 0
  }, logical(1))
  x[, keep, drop = FALSE]
}

remove_correlated_features <- function(x, cutoff = 0.90) {
  cm <- abs(cor(x, use = "pairwise.complete.obs", method = "spearman"))
  diag(cm) <- 0
  keep <- rep(TRUE, ncol(x))
  names(keep) <- colnames(x)

  while (any(cm[keep, keep, drop = FALSE] >= cutoff)) {
    current <- cm[keep, keep, drop = FALSE]
    pair <- which(current == max(current), arr.ind = TRUE)[1, ]
    candidate <- colnames(current)[pair]
    mean_cor <- colMeans(current[, candidate, drop = FALSE])
    keep[candidate[which.max(mean_cor)]] <- FALSE
  }
  x[, keep, drop = FALSE]
}

fit_scaler <- function(x) {
  list(
    center = vapply(x, mean, numeric(1), na.rm = TRUE),
    scale = vapply(x, sd, numeric(1), na.rm = TRUE)
  )
}

apply_scaler <- function(x, scaler) {
  scaler$scale[!is.finite(scaler$scale) | scaler$scale == 0] <- 1
  out <- sweep(as.matrix(x), 2, scaler$center, "-")
  sweep(out, 2, scaler$scale, "/")
}

lasso_select <- function(x, y, folds) {
  cvfit <- cv.glmnet(
    x = as.matrix(x), y = y, family = "binomial",
    foldid = folds, standardize = FALSE
  )
  beta <- coef(cvfit, s = "lambda.1se")
  rownames(beta)[as.vector(beta != 0) & rownames(beta) != "(Intercept)"]
}


development <- read_required(path$development)
train_meta <- development %>% filter(split == "Training")
test_meta <- development %>% filter(split == "Testing")

train_features <- feature_data %>%
  semi_join(train_meta, by = "patient_id") %>%
  arrange(match(patient_id, train_meta$patient_id))

test_features <- feature_data %>%
  semi_join(test_meta, by = "patient_id") %>%
  arrange(match(patient_id, test_meta$patient_id))

x_train_all <- train_features %>% select(all_of(c(rad_cols, dl_cols)))
x_train_all <- remove_near_zero_variance(x_train_all)
x_train_all <- remove_correlated_features(x_train_all)
scaler <- fit_scaler(x_train_all)
x_train_scaled <- apply_scaler(x_train_all, scaler)

fold_id <- integer(nrow(train_meta))
for (outcome_level in sort(unique(train_meta$non_pcr))) {
  idx <- which(train_meta$non_pcr == outcome_level)
  fold_id[idx] <- sample(rep(seq_len(5), length.out = length(idx)))
}

selected_features <- lasso_select(
  x_train_scaled, train_meta$non_pcr, folds = fold_id
)

# -----------------------------------------------------------------------------
# 4. Stage 1 XGBoost model: externalized study-specific configuration
# -----------------------------------------------------------------------------

load_stage1_configuration <- function(file) {
  if (!file.exists(file)) {
    stop(
      "The study-specific model configuration is not distributed in this ",
      "script. Provide the locked parameter file to refit the original model."
    )
  }
  readRDS(file)
}

fit_stage1_xgboost <- function(x, y, configuration) {
  y_numeric <- if (is.factor(y)) as.numeric(as.character(y)) else as.numeric(y)
  dtrain <- xgb.DMatrix(data = as.matrix(x), label = y_numeric)
  xgb.train(
    params = configuration$params,
    data = dtrain,
    nrounds = configuration$nrounds,
    verbose = 0
  )
}

predict_stage1_xgboost <- function(model, x) {
  as.numeric(predict(model, xgb.DMatrix(as.matrix(x))))
}

make_stratified_folds <- function(y, v = 5, seed = 2026) {
  set.seed(seed)
  fold <- integer(length(y))
  for (level in sort(unique(y))) {
    idx <- which(y == level)
    fold[idx] <- sample(rep(seq_len(v), length.out = length(idx)))
  }
  fold
}

# All preprocessing objects used for a validation fold are estimated from the
# corresponding fold-training data only.
fit_fold_preprocessor <- function(x, y, seed) {
  x <- as.data.frame(x, check.names = FALSE)
  x <- remove_near_zero_variance(x)
  x <- remove_correlated_features(x)
  fold_scaler <- fit_scaler(x)
  x_scaled <- apply_scaler(x, fold_scaler)
  inner_fold <- make_stratified_folds(y, v = 5, seed = seed)
  fold_features <- lasso_select(x_scaled, y, folds = inner_fold)
  if (length(fold_features) == 0) {
    fold_features <- colnames(x_scaled)
  }
  list(
    input_features = colnames(x),
    scaler = fold_scaler,
    selected_features = fold_features
  )
}

apply_fold_preprocessor <- function(x, preprocessor) {
  missing_features <- setdiff(preprocessor$input_features, colnames(x))
  if (length(missing_features) > 0) {
    stop("Missing Stage 1 input features: ", paste(missing_features, collapse = ", "))
  }
  x <- as.data.frame(x[, preprocessor$input_features, drop = FALSE])
  x_scaled <- apply_scaler(x, preprocessor$scaler)
  x_scaled[, preprocessor$selected_features, drop = FALSE]
}

# SMOTE is applied only after a fold has been split and only to its training
# portion; validation patients are never resampled.
smote_training_fold <- function(x, y, seed) {
  set.seed(seed)
  training_data <- data.frame(
    non_pcr = factor(y, levels = c(0, 1)),
    as.data.frame(x, check.names = FALSE),
    check.names = FALSE
  )
  smote_recipe <- recipes::recipe(non_pcr ~ ., data = training_data) %>%
    themis::step_smote(non_pcr)
  balanced <- recipes::prep(
    smote_recipe, training = training_data, retain = TRUE
  ) %>%
    recipes::juice()
  list(
    x = balanced %>% select(-non_pcr) %>% as.matrix(),
    y = as.numeric(as.character(balanced$non_pcr))
  )
}

fit_stage1_pipeline <- function(x, y, configuration, seed) {
  preprocessor <- fit_fold_preprocessor(x, y, seed = seed)
  x_processed <- apply_fold_preprocessor(x, preprocessor)
  balanced <- smote_training_fold(x_processed, y, seed = seed)
  model <- fit_stage1_xgboost(
    balanced$x, balanced$y, configuration = configuration
  )
  list(preprocessor = preprocessor, model = model)
}

predict_stage1_pipeline <- function(pipeline, x) {
  x_processed <- apply_fold_preprocessor(x, pipeline$preprocessor)
  predict_stage1_xgboost(pipeline$model, x_processed)
}

# Study-specific XGBoost settings are stored externally. This keeps the public
# script executable without disclosing the final tuned parameter values.
stage1_configuration <- load_stage1_configuration(path$stage1_configuration)
model_feature_sets <- list(
  RAD = intersect(rad_cols, colnames(x_train_all)),
  DL = intersect(dl_cols, colnames(x_train_all)),
  DR = intersect(c(rad_cols, dl_cols), colnames(x_train_all))
)

# Generate pooled out-of-fold probabilities for each Stage 1 model.
stage1_oof <- imap_dfr(model_feature_sets, function(feature_names, model_name) {
  oof_probability <- rep(NA_real_, nrow(train_meta))

  for (fold in sort(unique(fold_id))) {
    analysis_idx <- which(fold_id != fold)
    assessment_idx <- which(fold_id == fold)
    fold_pipeline <- fit_stage1_pipeline(
      x = x_train_all[analysis_idx, feature_names, drop = FALSE],
      y = train_meta$non_pcr[analysis_idx],
      configuration = stage1_configuration$models[[model_name]],
      seed = 2026 + fold
    )
    oof_probability[assessment_idx] <- predict_stage1_pipeline(
      fold_pipeline,
      x_train_all[assessment_idx, feature_names, drop = FALSE]
    )
  }

  tibble(
    patient_id = train_meta$patient_id,
    cohort = "Training",
    non_pcr = train_meta$non_pcr,
    model = model_name,
    prediction_source = "OOF",
    p_non_pcr = clip_probability(oof_probability)
  )
})

# Select the operating threshold from pooled DR OOF probabilities using the
# sensitivity-prioritized rule specified in the external configuration.
select_operating_threshold <- function(truth, probability, minimum_sensitivity) {
  candidates <- sort(unique(probability))
  performance <- map_dfr(candidates, function(threshold) {
    predicted <- as.integer(probability >= threshold)
    tp <- sum(predicted == 1 & truth == 1)
    tn <- sum(predicted == 0 & truth == 0)
    fp <- sum(predicted == 1 & truth == 0)
    fn <- sum(predicted == 0 & truth == 1)
    tibble(
      threshold = threshold,
      sensitivity = tp / (tp + fn),
      specificity = tn / (tn + fp)
    )
  })

  eligible <- performance %>% filter(sensitivity >= minimum_sensitivity)
  if (nrow(eligible) == 0) {
    stop("No OOF threshold satisfied the sensitivity-prioritized criterion.")
  }
  eligible %>%
    arrange(desc(specificity), desc(sensitivity), desc(threshold)) %>%
    slice(1) %>%
    pull(threshold)
}

operating_threshold <- stage1_oof %>%
  filter(model == "DR") %>%
  summarise(
    threshold = select_operating_threshold(
      truth = non_pcr,
      probability = p_non_pcr,
      minimum_sensitivity = stage1_configuration$minimum_sensitivity
    )
  ) %>%
  pull(threshold)

# Refit each preprocessing pipeline and XGBoost classifier on the complete
# training set. These objects are then locked for all subsequent predictions.
locked_stage1_models <- imap(model_feature_sets, function(feature_names, model_name) {
  fit_stage1_pipeline(
    x = x_train_all[, feature_names, drop = FALSE],
    y = train_meta$non_pcr,
    configuration = stage1_configuration$models[[model_name]],
    seed = 2026
  )
})
locked_stage1_model <- locked_stage1_models$DR

# Assemble patient identifiers and outcomes only; probabilities are generated
# below by the locked model rather than read from a precomputed file.
prediction_index <- bind_rows(
  test_meta %>% transmute(patient_id, cohort = "Testing", non_pcr),
  read_required(path$fuscc) %>%
    transmute(patient_id, cohort = "FUSCC", non_pcr),
  read_required(path$yncc) %>%
    transmute(patient_id, cohort = "YNCC", non_pcr),
  read_required(path$fuscc_nat) %>%
    transmute(
      patient_id, cohort = "FUSCC NAT",
      non_pcr = 1L - as.integer(as.character(pcr))
    ),
  read_required(path$ispy2_nat) %>%
    transmute(
      patient_id, cohort = "I-SPY2 NAT",
      non_pcr = 1L - as.integer(as.character(pcr))
    )
) %>%
  inner_join(feature_data, by = "patient_id")

locked_predictions <- imap_dfr(model_feature_sets, function(feature_names, model_name) {
  tibble(
    patient_id = prediction_index$patient_id,
    cohort = prediction_index$cohort,
    non_pcr = prediction_index$non_pcr,
    model = model_name,
    prediction_source = "Locked model",
    p_non_pcr = predict_stage1_pipeline(
      locked_stage1_models[[model_name]],
      prediction_index[, feature_names, drop = FALSE]
    )
  )
})

# Downstream analyses use the DR probabilities. Training values are pooled OOF
# predictions; all other values come from the locked full-training model.
predictions <- bind_rows(
  stage1_oof %>% filter(model == "DR"),
  locked_predictions %>% filter(model == "DR")
) %>%
  mutate(
    p_non_pcr = clip_probability(p_non_pcr),
    p_pcr_nac = 1 - p_non_pcr,
    predicted_non_pcr = as.integer(p_non_pcr >= operating_threshold)
  )

# -----------------------------------------------------------------------------
# 5. Stage 1 discrimination, threshold metrics, and calibration
# -----------------------------------------------------------------------------

auc_with_ci <- function(data) {
  roc_obj <- roc(
    response = data$non_pcr,
    predictor = data$p_non_pcr,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  ci_obj <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)
  tibble(
    auc = as.numeric(auc(roc_obj)),
    lower = as.numeric(ci_obj[1]),
    upper = as.numeric(ci_obj[3])
  )
}

threshold_metrics <- function(data, threshold = operating_threshold) {
  truth <- as.integer(as.character(data$non_pcr))
  estimate <- as.integer(data$p_non_pcr >= threshold)
  tp <- sum(estimate == 1 & truth == 1, na.rm = TRUE)
  tn <- sum(estimate == 0 & truth == 0, na.rm = TRUE)
  fp <- sum(estimate == 1 & truth == 0, na.rm = TRUE)
  fn <- sum(estimate == 0 & truth == 1, na.rm = TRUE)

  binom_ci <- function(success, total) {
    if (total == 0) {
      return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    ci <- binom.test(success, total)$conf.int
    c(estimate = success / total, lower = ci[1], upper = ci[2])
  }

  bind_rows(
    sensitivity = binom_ci(tp, tp + fn),
    specificity = binom_ci(tn, tn + fp),
    ppv = binom_ci(tp, tp + fp),
    npv = binom_ci(tn, tn + fn),
    .id = "metric"
  ) %>%
    mutate(
      n = nrow(data),
      actual_non_pcr = sum(truth == 1, na.rm = TRUE),
      predicted_non_pcr = sum(estimate == 1, na.rm = TRUE)
    )
}

calibration_statistics <- function(data) {
  p <- clip_probability(data$p_non_pcr)
  fit <- glm(data$non_pcr ~ qlogis(p), family = binomial())
  tibble(
    intercept = unname(coef(fit)[1]),
    slope = unname(coef(fit)[2]),
    brier = mean((data$non_pcr - p)^2, na.rm = TRUE)
  )
}

bootstrap_statistic <- function(data, indices, statistic_function) {
  as.numeric(statistic_function(data[indices, , drop = FALSE]))
}

calibration_bootstrap <- function(data, repetitions = 2000) {
  statistic_function <- function(d) {
    unlist(calibration_statistics(d)[1, c("intercept", "slope", "brier")])
  }
  b <- boot(
    data = data,
    statistic = function(d, i) bootstrap_statistic(d, i, statistic_function),
    R = repetitions
  )
  point <- statistic_function(data)
  ci <- apply(b$t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  tibble(
    metric = c("intercept", "slope", "brier"),
    estimate = point,
    lower = ci[1, ],
    upper = ci[2, ]
  )
}

stage1_auc <- predictions %>%
  group_by(cohort) %>%
  group_modify(~ auc_with_ci(.x)) %>%
  ungroup()

stage1_threshold <- predictions %>%
  group_by(cohort) %>%
  group_modify(~ threshold_metrics(.x, threshold = operating_threshold)) %>%
  ungroup()

stage1_calibration <- predictions %>%
  group_by(cohort) %>%
  group_modify(~ calibration_bootstrap(.x, repetitions = 2000)) %>%
  ungroup()


calibration_curve_data <- predictions %>%
  group_by(cohort) %>%
  group_modify(~ {
    cal_df <- data.frame(
      non_pcr = as.numeric(as.character(.x$non_pcr)),
      p = clip_probability(.x$p_non_pcr)
    )
    fit <- glm(non_pcr ~ qlogis(p), data = cal_df, family = binomial())
    grid <- tibble(p_non_pcr = seq(0.01, 0.99, length.out = 200))
    pr <- predict(fit, newdata = data.frame(p = grid$p_non_pcr),
                  type = "link", se.fit = TRUE)
    grid %>%
      mutate(
        observed = plogis(pr$fit),
        lower = plogis(pr$fit - 1.96 * pr$se.fit),
        upper = plogis(pr$fit + 1.96 * pr$se.fit)
      )
  }) %>%
  ungroup()

ggplot(calibration_curve_data,
       aes(x = p_non_pcr, y = observed, color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability of non-pCR",
       y = "Observed probability of non-pCR") +
  theme_classic()

# -----------------------------------------------------------------------------
# 6. Import frozen Stage 2 assignments and evaluate Stage 1-Stage 2 association
# -----------------------------------------------------------------------------

subtype_data <- read_required(path$subtype) %>%
  inner_join(
    predictions %>% select(patient_id, cohort, non_pcr, p_non_pcr),
    by = c("patient_id", "cohort")
  ) %>%
  filter(non_pcr == 1) %>%
  mutate(subtype_binary = as.integer(subtype == "Subtype B"))

stage1_stage2_association <- subtype_data %>%
  group_by(cohort) %>%
  group_modify(~ {
    wilcox_result <- wilcox.test(p_non_pcr ~ subtype, data = .x, exact = FALSE)
    correlation <- cor.test(
      .x$p_non_pcr, .x$subtype_binary,
      method = "pearson"
    )
    roc_obj <- roc(
      .x$subtype_binary, .x$p_non_pcr,
      levels = c(0, 1), direction = "<", quiet = TRUE
    )
    roc_ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)

    tibble(
      median_a = median(.x$p_non_pcr[.x$subtype == "Subtype A"]),
      median_b = median(.x$p_non_pcr[.x$subtype == "Subtype B"]),
      median_difference = median_b - median_a,
      wilcoxon_p = wilcox_result$p.value,
      r_pb = unname(correlation$estimate),
      correlation_p = correlation$p.value,
      subtype_auc = as.numeric(auc(roc_obj)),
      auc_lower = as.numeric(roc_ci[1]),
      auc_upper = as.numeric(roc_ci[3])
    )
  }) %>%
  ungroup()

habitat_patterns <- read_required(path$habitat_pattern) %>%
  inner_join(subtype_data %>% select(patient_id, cohort, subtype),
             by = c("patient_id", "cohort"))

habitat_pattern_tests <- habitat_patterns %>%
  group_by(cohort) %>%
  group_modify(~ {
    tab <- table(.x$subtype, .x$spatial_pattern)
    if (all(dim(tab) >= 2) && all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
      tibble(p_value = chisq.test(tab)$p.value)
    } else {
      tibble(p_value = NA_real_)
    }
  }) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 7. Clinicopathologic associations and DFS
# -----------------------------------------------------------------------------

external_clinical <- bind_rows(
  read_required(path$fuscc),
  read_required(path$yncc)
) %>%
  inner_join(
    subtype_data %>% select(patient_id, cohort, subtype, p_non_pcr),
    by = c("patient_id", "cohort")
  ) %>%
  mutate(subtype = relevel(factor(subtype), ref = "Subtype A"))

binary_outcomes <- c(
  "grade3", "n2_3", "lvi_positive", "luminal_b",
  "her2_low", "ki67_high", "mp_grade4_5"
)

fit_binary_association <- function(data, outcome) {
  available <- data %>% filter(!is.na(.data[[outcome]]))
  if (nrow(available) == 0 || length(unique(available[[outcome]])) < 2) {
    return(tibble())
  }
  unadjusted <- glm(
    reformulate("subtype", response = outcome),
    data = available, family = binomial()
  )
  adjusted <- glm(
    reformulate(c("subtype", "p_non_pcr"), response = outcome),
    data = available, family = binomial()
  )
  bind_rows(
    tidy(unadjusted, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(model = "Unadjusted"),
    tidy(adjusted, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(model = "Adjusted")
  ) %>%
    filter(term == "subtypeSubtype B") %>%
    mutate(outcome = outcome)
}

clinicopathologic_models <- external_clinical %>%
  group_by(cohort) %>%
  group_modify(~ {
    cohort_data <- .x
    map_dfr(
      binary_outcomes,
      function(outcome) fit_binary_association(cohort_data, outcome)
    )
  }) %>%
  ungroup()

dfs_data <- external_clinical %>%
  filter(
    (cohort == "FUSCC" & potential_followup_years >= 5) |
      cohort == "YNCC"
  )

dfs_models <- dfs_data %>%
  group_by(cohort) %>%
  group_modify(~ {
    unadjusted <- coxph(Surv(dfs_time, dfs_event) ~ subtype, data = .x)
    adjusted <- coxph(
      Surv(dfs_time, dfs_event) ~ subtype + p_non_pcr,
      data = .x
    )
    bind_rows(
      tidy(unadjusted, exponentiate = TRUE, conf.int = TRUE) %>%
        mutate(model = "Unadjusted"),
      tidy(adjusted, exponentiate = TRUE, conf.int = TRUE) %>%
        mutate(model = "Adjusted")
    ) %>%
      filter(term == "subtypeSubtype B")
  }) %>%
  ungroup()

logrank_results <- dfs_data %>%
  group_by(cohort) %>%
  group_modify(~ {
    test <- survdiff(Surv(dfs_time, dfs_event) ~ subtype, data = .x)
    tibble(logrank_p = pchisq(test$chisq, df = 1, lower.tail = FALSE))
  }) %>%
  ungroup()

km_fit <- survfit(Surv(dfs_time, dfs_event) ~ subtype + cohort, data = dfs_data)
ggsurvplot_facet(
  km_fit,
  data = dfs_data,
  facet.by = "cohort",
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("#D62728", "#2CA02C"),
  xlab = "Time since surgery (years)",
  ylab = "Disease-free survival"
)

# -----------------------------------------------------------------------------
# 8. SNF, mutation, RNA, and RPPA analyses
# -----------------------------------------------------------------------------

snf_data <- read_required(path$snf) %>%
  inner_join(subtype_data %>% select(patient_id, subtype), by = "patient_id")

snf_distribution <- snf_data %>%
  count(subtype, snf_subtype) %>%
  group_by(subtype) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

mutation_data <- read_required(path$mutation) %>%
  inner_join(subtype_data %>% select(patient_id, subtype), by = "patient_id")

mutation_genes <- setdiff(
  names(mutation_data),
  c("patient_id", "cohort", "subtype")
)

mutation_results <- map_dfr(mutation_genes, function(gene) {
  tab <- table(mutation_data$subtype, mutation_data[[gene]])
  if (all(dim(tab) >= 2) && all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
    test <- fisher.test(tab)
    tibble(gene = gene, p_value = test$p.value)
  } else {
    tibble(gene = gene, p_value = NA_real_)
  }
}) %>%
  mutate(
    adjusted_p = p.adjust(p_value, method = "BH"),
    significance = case_when(
      adjusted_p < .001 ~ "***",
      adjusted_p < .01 ~ "**",
      adjusted_p < .05 ~ "*",
      TRUE ~ "ns"
    )
  )


rna_results <- read_required(path$rna) %>%
  mutate(adjusted_p = p.adjust(p_value, method = "BH"))

rppa <- read_required(path$rppa) %>%
  inner_join(subtype_data %>% select(patient_id, subtype), by = "patient_id")

rppa_pathways <- c(
  "PI3K", "MAPK", "RTK", "Proliferation",
  "Immune_Inflammation", "ER_Luminal"
)

rppa_results <- map_dfr(rppa_pathways, function(pathway) {
  test <- wilcox.test(rppa[[pathway]] ~ rppa$subtype, exact = FALSE)
  tibble(
    pathway = pathway,
    subtype_a_median = median(
      rppa[[pathway]][rppa$subtype == "Subtype A"], na.rm = TRUE
    ),
    subtype_b_median = median(
      rppa[[pathway]][rppa$subtype == "Subtype B"], na.rm = TRUE
    ),
    p_value = test$p.value
  )
})

# -----------------------------------------------------------------------------
# 9. End-to-end Stage 1-4 application in FUSCC NAT and I-SPY2 NAT cohorts
# -----------------------------------------------------------------------------

nat_clinical <- bind_rows(
  read_required(path$fuscc_nat),
  read_required(path$ispy2_nat)
)

nat_analysis <- nat_clinical %>%
  inner_join(
    predictions %>%
      select(patient_id, cohort, p_non_pcr, predicted_non_pcr),
    by = c("patient_id", "cohort")
  ) %>%
  left_join(
    read_required(path$subtype) %>% select(patient_id, cohort, subtype),
    by = c("patient_id", "cohort")
  ) %>%
  mutate(
    pcr = as.integer(as.character(pcr)),
    p_pcr_nac = 1 - p_non_pcr,
    delta = pcr - p_pcr_nac
  )

end_to_end_counts <- nat_analysis %>%
  group_by(cohort) %>%
  summarise(
    total_nat = n(),
    predicted_non_pcr = sum(predicted_non_pcr == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    nat_analysis %>%
      filter(predicted_non_pcr == 1) %>%
      count(cohort, subtype, pcr) %>%
      pivot_wider(
        names_from = c(subtype, pcr), values_from = n,
        values_fill = 0
      ),
    by = "cohort"
  )

stage4_data <- nat_analysis %>%
  filter(predicted_non_pcr == 1, !is.na(subtype))

bootstrap_cell <- function(data, repetitions = 2000) {
  statistic <- function(d, i) mean(d$delta[i], na.rm = TRUE)
  b <- boot(data, statistic = statistic, R = repetitions)
  ci <- quantile(b$t, c(0.025, 0.975), na.rm = TRUE)
  tibble(
    n = nrow(data),
    observed_pcr_n = sum(data$pcr, na.rm = TRUE),
    observed_pcr_rate = mean(data$pcr, na.rm = TRUE),
    expected_pcr = mean(data$p_pcr_nac, na.rm = TRUE),
    mean_delta = mean(data$delta, na.rm = TRUE),
    lower = ci[1],
    upper = ci[2],
    oer = mean(data$pcr, na.rm = TRUE) / mean(data$p_pcr_nac, na.rm = TRUE)
  )
}

opportunity_matrix <- stage4_data %>%
  group_by(cohort, subtype, regimen_group) %>%
  group_modify(~ bootstrap_cell(.x, repetitions = 2000)) %>%
  ungroup()

# Omics-concordant regimens are defined a priori from Stage 3 annotations.
stage4_data <- stage4_data %>%
  mutate(
    omics_concordant = case_when(
      subtype == "Subtype A" &
        regimen_group %in% c("Endocrine", "Immune") ~ TRUE,
      subtype == "Subtype B" &
        regimen_group %in% c("PARP/DDR", "RTK/MAPK-Angio", "Ganetespib") ~ TRUE,
      TRUE ~ FALSE
    )
  )

bootstrap_concordance <- function(data, repetitions = 2000) {
  # Guard against cases where one concordance group is entirely missing
  if (length(unique(data$omics_concordant)) < 2) {
    return(tibble(
      difference_in_delta = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      p_value = NA_real_
    ))
  }

  statistic <- function(d, i) {
    sampled <- d[i, , drop = FALSE]
    concordant <- sampled$delta[sampled$omics_concordant]
    discordant <- sampled$delta[!sampled$omics_concordant]
    if (length(concordant) == 0 || length(discordant) == 0) {
      return(NA_real_)
    }
    mean(concordant, na.rm = TRUE) - mean(discordant, na.rm = TRUE)
  }
  b <- boot(data, statistic = statistic, R = repetitions,
            strata = data$omics_concordant)
  observed <- statistic(data, seq_len(nrow(data)))
  ci <- quantile(b$t, c(0.025, 0.975), na.rm = TRUE)
  p_value <- min(
    1,
    2 * min(mean(b$t <= 0, na.rm = TRUE), mean(b$t >= 0, na.rm = TRUE))
  )
  tibble(
    difference_in_delta = observed,
    lower = ci[1],
    upper = ci[2],
    p_value = p_value
  )
}

omics_concordance <- stage4_data %>%
  group_by(cohort) %>%
  group_modify(~ bootstrap_concordance(.x, repetitions = 2000)) %>%
  ungroup()
