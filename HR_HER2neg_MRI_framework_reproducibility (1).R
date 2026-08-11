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
# 1. File locations and required variables
# -----------------------------------------------------------------------------

path <- list(
  development = "data/development_cohort.csv",
  fuscc = "data/fuscc_validation_cohort.csv",
  yncc = "data/yncc_validation_cohort.csv",
  fuscc_nat = "data/fuscc_nat_cohort.csv",
  ispy2_nat = "data/ispy2_nat_cohort.csv",
  radiomics = "data/radiomics_features.csv",
  dl = "data/dl_features_from_public_pipeline.csv",
  subtype = "data/frozen_stage2_subtype_assignments.csv",
  habitat_pattern = "data/spatial_pattern_from_public_pipeline.csv",
  locked_predictions = "data/locked_stage1_predictions.csv",
  locked_stage1_model = "model/locked_stage1_xgboost.json",
  stage1_parameters = "config/stage1_model_parameters.rds",
  rna = "data/rna_pathway_results.csv",
  rppa = "data/rppa_pathway_scores.csv",
  mutation = "data/fuscc_mutation_matrix.csv",
  snf = "data/fuscc_snf_subtypes.csv"
)

# Minimum patient-level fields expected in the relevant files:
# patient_id, cohort, split, center, pcr, non_pcr, subtype,
# dfs_time, dfs_event, regimen_group, and the clinicopathologic variables.
# pcr is coded 1 for pCR and 0 for non-pCR; non_pcr = 1 - pcr.

read_required <- function(file) {
  if (!file.exists(file)) stop("Required file not found: ", file)
  fread(file, data.table = FALSE)
}

clip_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

# -----------------------------------------------------------------------------
# 2. Assemble radiomics and externally extracted DL features
# -----------------------------------------------------------------------------

radiomics <- read_required(path$radiomics)
dl_features <- read_required(path$dl)

feature_data <- radiomics %>%
  inner_join(dl_features, by = "patient_id")

# DL feature extraction is not repeated here. The imported DL table should
# contain the patient-level 192-dimensional ViT-Tiny representation generated
# using the cited public implementation.

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

# Reproducibility analyses should use the prespecified training/testing split.
development <- read_required(path$development)
train_meta <- development %>% filter(split == "Training")
test_meta <- development %>% filter(split == "Testing")

train_features <- feature_data %>%
  semi_join(train_meta, by = "patient_id") %>%
  arrange(match(patient_id, train_meta$patient_id))

test_features <- feature_data %>%
  semi_join(test_meta, by = "patient_id") %>%
  arrange(match(patient_id, test_meta$patient_id))

# All filtering and fitted transformations must be estimated from training data
# only and transferred unchanged to held-out or external data.
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

# Exact paper-level downstream results are reproduced from locked, patient-level
# out-of-fold/testing/external probabilities. To refit the locked model, use:
# stage1_configuration <- load_stage1_configuration(path$stage1_parameters)
# stage1_model <- fit_stage1_xgboost(
#   x_train_scaled[, selected_features, drop = FALSE],
#   train_meta$non_pcr,
#   stage1_configuration
# )
# xgb.save(stage1_model, path$locked_stage1_model)

predictions <- read_required(path$locked_predictions) %>%
  mutate(
    p_non_pcr = clip_probability(p_non_pcr),
    p_pcr_nac = 1 - p_non_pcr,
    predicted_non_pcr = as.integer(p_non_pcr >= 0.58)
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

threshold_metrics <- function(data, threshold = 0.58) {
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
  group_modify(~ threshold_metrics(.x, threshold = 0.58)) %>%
  ungroup()

stage1_calibration <- predictions %>%
  group_by(cohort) %>%
  group_modify(~ calibration_bootstrap(.x, repetitions = 2000)) %>%
  ungroup()

# Logistic calibration curves are used for evaluation only; fitted values from
# these curves are not used to replace or modify locked model probabilities.
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

# Consensus clustering, frozen nearest-centroid assignment, and spatial-pattern
# classification are not reimplemented here because those procedures use the
# cited public code. Their patient-level outputs are imported above and below.

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

# FUSCC DFS is restricted to the prespecified subset with at least 5 years of
# potential follow-up; all eligible YNCC non-pCR patients are retained.
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

# RNA pathway enrichment is performed with the cited pathway-analysis workflow;
# its patient- or pathway-level result table is imported for summarization.
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

# -----------------------------------------------------------------------------
# 10. Export analysis tables and software environment
# -----------------------------------------------------------------------------

dir.create("results", showWarnings = FALSE, recursive = TRUE)

write_csv(stage1_auc, "results/stage1_auc.csv")
write_csv(stage1_threshold, "results/stage1_threshold_metrics.csv")
write_csv(stage1_calibration, "results/stage1_calibration.csv")
write_csv(stage1_stage2_association, "results/stage1_stage2_association.csv")
write_csv(habitat_pattern_tests, "results/habitat_pattern_tests.csv")
write_csv(clinicopathologic_models, "results/clinicopathologic_models.csv")
write_csv(dfs_models, "results/dfs_cox_models.csv")
write_csv(logrank_results, "results/dfs_logrank_tests.csv")
write_csv(snf_distribution, "results/snf_distribution.csv")
write_csv(mutation_results, "results/mutation_results.csv")
write_csv(rna_results, "results/rna_pathway_results.csv")
write_csv(rppa_results, "results/rppa_pathway_results.csv")
write_csv(end_to_end_counts, "results/end_to_end_counts.csv")
write_csv(opportunity_matrix, "results/opportunity_matrix.csv")
write_csv(omics_concordance, "results/omics_concordance.csv")

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
