# Run from the repository root: Rscript scripts/q4/q4_analysis.R
# Required packages: glmnet, randomForest, rsample; ggplot2 for the report.
# The execution block saves report tables and figures.

# 1. Load the measurements and group labels
# We keep the original data reading and response encoding. Raw measurements let us
# estimate preprocessing within training folds; the clean data reproduce class.

load_q4_data <- function(project_root = ".") {
  raw_path <- file.path(project_root, "data", "biomarker-raw.csv")
  headers <- read.csv(raw_path, header = FALSE, nrows = 2, check.names = FALSE)
  raw <- read.csv(raw_path, header = FALSE, skip = 2,
                  na.strings = c("-", ""), check.names = FALSE)
  raw <- raw[!is.na(raw[[1]]), ]
  protein_columns <- 3:(ncol(raw) - 1)
  protein_names <- as.character(unlist(headers[2, protein_columns]))
  protein_levels <- as.matrix(raw[, protein_columns])
  colnames(protein_levels) <- protein_names
  clean_environment <- new.env()
  load(file.path(project_root, "data", "biomarker-clean.RData"), clean_environment)
  clean <- clean_environment$biomarker_clean
  stopifnot(identical(as.character(raw[[1]]), as.character(clean$group)),
            identical(protein_names, names(clean)[-(1:2)]),
            all(is.finite(protein_levels)), all(protein_levels > 0))
  list(log_levels = log10(protein_levels), outcome = as.integer(raw[[1]] == "ASD"),
       clean_levels = as.matrix(clean[, -(1:2)]),
       protein_names = protein_names,
       protein_full_names = as.character(unlist(headers[1, protein_columns])))
}

# 2. Prepare training and validation measurements
# Log values are centered, scaled, and clipped using training parameters only.
# This preserves the professor's transformation without using validation subjects.
preprocess_partition <- function(training_levels, validation_levels) {
  training_mean <- colMeans(training_levels)
  training_sd <- apply(training_levels, 2, sd)
  transform_levels <- function(levels) {
    scaled <- sweep(sweep(levels, 2, training_mean), 2, training_sd, "/")
    scaled[scaled > 3] <- 3
    scaled[scaled < -3] <- -3
    scaled
  }
  list(training = transform_levels(training_levels),
       validation = transform_levels(validation_levels))
}

# 3. Create balanced validation folds
# Stratification keeps both groups represented. The same outer folds are used
# for every model so that performance differences are paired by subject.
stratified_folds <- function(outcome, folds = 5) {
  fold_assignment <- integer(length(outcome))
  for (group in c(0L, 1L)) {
    group_indices <- which(outcome == group)
    fold_assignment[group_indices] <- sample(rep(seq_len(folds), length.out = length(group_indices)))
  }
  fold_assignment
}

# 4. Rank proteins and reproduce the in-class baseline
# Welch tests measure group differences. The baseline retains the professor's
# adjusted ranking and intersects its top ten with random forest's top ten.
welch_p_values <- function(protein_levels, outcome) {
  apply(protein_levels, 2, function(levels)
    t.test(levels[outcome == 1], levels[outcome == 0], var.equal = FALSE)$p.value)
}

select_baseline_panel <- function(protein_levels, outcome) {
  p_values <- welch_p_values(protein_levels, outcome)
  protein_count <- length(p_values)
  p_order <- order(p_values)
  # Reproduce the professor's rankwise adjustment, including its lack of cummin.
  adjusted_p <- protein_count * (log(protein_count) + 1 / (2 * protein_count) -
                                 digamma(1)) * p_values[p_order] / seq_len(protein_count)
  test_panel <- names(adjusted_p)[adjusted_p <= sort(adjusted_p)[10]]
  forest <- randomForest::randomForest(
    x = protein_levels, y = factor(outcome, levels = c(1, 0), labels = c("ASD", "TD")),
    ntree = 1000, importance = TRUE)
  importance <- forest$importance[, "MeanDecreaseGini"]
  forest_panel <- names(importance)[importance >= sort(importance, decreasing = TRUE)[10]]
  list(proteins = intersect(test_panel, forest_panel),
       test_panel = test_panel, forest_panel = forest_panel)
}

# 5. Fit logistic regression for a selected panel
# Both the baseline and the small-panel procedure use this classifier. Sharing
# the fitting step keeps their differences focused on how proteins are selected.
fit_logistic_panel <- function(training_levels, outcome, validation_levels, proteins) {
  training_data <- data.frame(class = outcome, training_levels[, proteins, drop = FALSE],
                              check.names = FALSE)
  model <- glm(class ~ ., data = training_data, family = binomial())
  list(probabilities = as.numeric(predict(model,
         newdata = as.data.frame(validation_levels[, proteins, drop = FALSE]), type = "response")),
       coefficients = setNames(unname(coef(model)[-1]), proteins))
}

# 6. Preserve the existing LASSO comparison
# A training-specific penalty path permits comparison of relative penalties.
# Inner-fold deviance chooses the one-standard-error penalty, already used in
# the original Q4 implementation. The minimum-deviance size is retained only
# to explain the effect of stronger regularization, not to choose a winner.
fit_lasso_path <- function(protein_levels, outcome, penalty_ratios) {
  maximum_penalty <- glmnet::glmnet(protein_levels, outcome, family = "binomial",
                                     standardize = FALSE)$lambda[1]
  glmnet::glmnet(protein_levels, outcome, family = "binomial", standardize = FALSE,
                 lambda = maximum_penalty * penalty_ratios)
}

select_lasso_panel <- function(log_levels, outcome, validation_log_levels) {
  penalty_ratios <- exp(seq(0, log(0.001), length.out = 50))
  inner_folds <- stratified_folds(outcome)
  fold_deviance <- matrix(NA_real_, 5, length(penalty_ratios))
  for (fold in 1:5) {
    validation_indices <- which(inner_folds == fold)
    training_indices <- which(inner_folds != fold)
    partition <- preprocess_partition(log_levels[training_indices, , drop = FALSE],
                                      log_levels[validation_indices, , drop = FALSE])
    model <- fit_lasso_path(partition$training, outcome[training_indices], penalty_ratios)
    probabilities <- predict(model, newx = partition$validation, type = "response")
    probabilities <- pmin(pmax(probabilities, 1e-15), 1 - 1e-15)
    fold_deviance[fold, ] <- colMeans(-2 * (outcome[validation_indices] * log(probabilities) +
                                           (1 - outcome[validation_indices]) * log1p(-probabilities)))
  }
  mean_deviance <- colMeans(fold_deviance)
  best_index <- which.min(mean_deviance)
  standard_error <- sd(fold_deviance[, best_index]) / sqrt(5)
  selected_index <- which(mean_deviance <= mean_deviance[best_index] + standard_error)[1]
  partition <- preprocess_partition(log_levels, validation_log_levels)
  model <- fit_lasso_path(partition$training, outcome, penalty_ratios)
  coefficients <- as.matrix(coef(model))[, selected_index]
  proteins <- names(coefficients)[coefficients != 0 & names(coefficients) != "(Intercept)"]
  list(probabilities = as.numeric(predict(model, newx = partition$validation,
                                         type = "response")[, selected_index]),
       proteins = proteins, coefficients = coefficients[proteins],
       penalty_ratio = penalty_ratios[selected_index],
       minimum_deviance_size = sum(as.matrix(coef(model))[-1, best_index] != 0))
}

# 7. Select a deliberately small panel using training data only
# We restrict the candidates to the top one, two, or three Welch-ranked proteins.
# In each inner fold we redo preprocessing, ranking, and logistic fitting.
# The smallest size within one standard error of the best mean deviance is used.
# The cap of three is fixed before evaluating this new procedure on outer folds.
select_small_panel <- function(log_levels, outcome, validation_log_levels) {
  candidate_sizes <- 1:3
  inner_folds <- stratified_folds(outcome)
  fold_deviance <- matrix(NA_real_, 5, length(candidate_sizes))
  for (fold in 1:5) {
    training_indices <- which(inner_folds != fold)
    validation_indices <- which(inner_folds == fold)
    partition <- preprocess_partition(log_levels[training_indices, ], log_levels[validation_indices, ])
    ranking <- names(sort(welch_p_values(partition$training, outcome[training_indices])))
    for (panel_size in candidate_sizes) {
      fitted <- fit_logistic_panel(partition$training, outcome[training_indices],
                                   partition$validation, ranking[seq_len(panel_size)])
      probabilities <- pmin(pmax(fitted$probabilities, 1e-15), 1 - 1e-15)
      observed <- outcome[validation_indices]
      fold_deviance[fold, panel_size] <- mean(-2 * (observed * log(probabilities) +
                                                    (1 - observed) * log1p(-probabilities)))
    }
  }
  mean_deviance <- colMeans(fold_deviance)
  best_size <- which.min(mean_deviance)
  standard_error <- sd(fold_deviance[, best_size]) / sqrt(5)
  selected_size <- which(mean_deviance <= mean_deviance[best_size] + standard_error)[1]
  partition <- preprocess_partition(log_levels, validation_log_levels)
  ranking <- names(sort(welch_p_values(partition$training, outcome)))
  proteins <- ranking[seq_len(selected_size)]
  fitted <- fit_logistic_panel(partition$training, outcome, partition$validation, proteins)
  list(probabilities = fitted$probabilities, proteins = proteins, coefficients = fitted$coefficients)
}

# 8. Measure classification performance
# ASD is the positive class. Accuracy uses the same 0.5 threshold for all models;
# AUROC uses probability ranks and gives tied predictions their average rank.
classification_metrics <- function(outcome, probabilities) {
  predicted_asd <- probabilities > 0.5
  positive_count <- sum(outcome == 1)
  negative_count <- sum(outcome == 0)
  auc <- (sum(rank(probabilities)[outcome == 1]) - positive_count * (positive_count + 1) / 2) /
    (positive_count * negative_count)
  c(Accuracy = mean(predicted_asd == outcome), AUROC = auc,
    Sensitivity = mean(predicted_asd[outcome == 1]),
    Specificity = mean(!predicted_asd[outcome == 0]))
}

# 9. Compare the procedures and identify the final candidate panel
# First reproduce the historical estimate, then evaluate all procedures on the
# original repeated outer folds. Inner selection never sees outer test subjects.
# The small-panel procedure is the primary new analysis; LASSO is a retained
# reference, not a competing procedure selected using the outer results.
run_q4_analysis <- function(project_root = ".", save_tables = FALSE) {
  data <- load_q4_data(project_root)
  outcome <- data$outcome
  # The historical result intentionally retains the original full-data selection.
  # It is a separate reference, not a leakage-free estimate for model comparison.
  set.seed(101422)
  historical_panel <- select_baseline_panel(data$clean_levels, outcome)
  set.seed(101422)
  historical_split <- rsample::initial_split(data.frame(subject = seq_along(outcome)), prop = 0.8)
  historical_training <- rsample::training(historical_split)$subject
  historical_testing <- rsample::testing(historical_split)$subject
  historical_probabilities <- fit_logistic_panel(
    data$clean_levels[historical_training, ], outcome[historical_training],
    data$clean_levels[historical_testing, ], historical_panel$proteins)$probabilities
  historical_metrics <- classification_metrics(outcome[historical_testing], historical_probabilities)

  # Generate the original outer partitions before fitting any models. Separate
  # seeds for the new procedure preserve the earlier baseline and LASSO fits.
  set.seed(197104)
  fold_assignments <- replicate(5, stratified_folds(outcome))
  baseline_predictions <- alternative_predictions <- matrix(NA_real_, length(outcome), 5)
  small_predictions <- matrix(NA_real_, length(outcome), 5)
  baseline_sizes <- alternative_sizes <- numeric(25)
  small_sizes <- minimum_deviance_sizes <- numeric(25)
  fold_differences <- numeric(25)
  small_panels <- vector("list", 25)
  small_fold_differences <- numeric(25)
  for (repeat_index in 1:5) {
    for (fold in 1:5) {
      partition_index <- (repeat_index - 1) * 5 + fold
      testing_indices <- which(fold_assignments[, repeat_index] == fold)
      training_indices <- which(fold_assignments[, repeat_index] != fold)
      stopifnot(length(intersect(training_indices, testing_indices)) == 0)
      partition <- preprocess_partition(data$log_levels[training_indices, ],
                                        data$log_levels[testing_indices, ])
      set.seed(197104 + partition_index)
      baseline_panel <- select_baseline_panel(partition$training, outcome[training_indices])$proteins
      baseline_probabilities <- fit_logistic_panel(partition$training, outcome[training_indices],
                                                   partition$validation, baseline_panel)$probabilities
      alternative <- select_lasso_panel(data$log_levels[training_indices, ], outcome[training_indices],
                                         data$log_levels[testing_indices, ])
      set.seed(197500 + partition_index)
      small <- select_small_panel(data$log_levels[training_indices, ], outcome[training_indices],
                                  data$log_levels[testing_indices, ])
      baseline_predictions[testing_indices, repeat_index] <- baseline_probabilities
      alternative_predictions[testing_indices, repeat_index] <- alternative$probabilities
      baseline_sizes[partition_index] <- length(baseline_panel)
      alternative_sizes[partition_index] <- length(alternative$proteins)
      minimum_deviance_sizes[partition_index] <- alternative$minimum_deviance_size
      small_predictions[testing_indices, repeat_index] <- small$probabilities
      small_sizes[partition_index] <- length(small$proteins)
      small_panels[[partition_index]] <- small$proteins
      small_fold_differences[partition_index] <- mean((small$probabilities > 0.5) == outcome[testing_indices]) -
        mean((baseline_probabilities > 0.5) == outcome[testing_indices])
      fold_differences[partition_index] <- mean((alternative$probabilities > 0.5) == outcome[testing_indices]) -
        mean((baseline_probabilities > 0.5) == outcome[testing_indices])
    }
    message("Completed validation repeat ", repeat_index, " of 5")
  }
  baseline_metrics <- apply(baseline_predictions, 2, function(p) classification_metrics(outcome, p))
  alternative_metrics <- apply(alternative_predictions, 2, function(p) classification_metrics(outcome, p))
  small_metrics <- apply(small_predictions, 2, function(p) classification_metrics(outcome, p))
  # Reapply the same training-only rules to all subjects to name final panels.
  # These fits do not supply additional performance estimates.
  set.seed(197999)
  final_lasso <- select_lasso_panel(data$log_levels, outcome, data$log_levels)
  set.seed(198000)
  final_panel <- select_small_panel(data$log_levels, outcome, data$log_levels)
  # Weight subjects equally within each repeat, then average the five repeats.
  # Equivalently, weight each paired fold difference by its test-fold size.
  accuracy_difference <- mean(alternative_metrics["Accuracy", ] - baseline_metrics["Accuracy", ])
  small_accuracy_difference <- mean(small_metrics["Accuracy", ] - baseline_metrics["Accuracy", ])
  # Center intervals on those same differences. Retain the corrected paired-fold
  # variance as an approximation for these nearly balanced folds (30-32 subjects).
  # The 1/4 test/training ratio and overlap correction are approximate, not exact.
  corrected_se <- sqrt((1 / 25 + 1 / 4) * var(fold_differences))
  difference_interval <- accuracy_difference + c(-1, 1) * qt(0.975, df = 24) * corrected_se
  small_corrected_se <- sqrt((1 / 25 + 1 / 4) * var(small_fold_differences))
  small_difference_interval <- small_accuracy_difference +
    c(-1, 1) * qt(0.975, df = 24) * small_corrected_se
  comparison <- data.frame(
    Model = c("Historical in-class", "Validated in-class procedure", "Validated LASSO procedure",
              "Validated small-panel procedure"),
    rbind(historical_metrics, rowMeans(baseline_metrics), rowMeans(alternative_metrics), rowMeans(small_metrics)),
    Mean_proteins = c(length(historical_panel$proteins), mean(baseline_sizes), mean(alternative_sizes), mean(small_sizes)),
    Min_proteins = c(length(historical_panel$proteins), min(baseline_sizes), min(alternative_sizes), min(small_sizes)),
    Max_proteins = c(length(historical_panel$proteins), max(baseline_sizes), max(alternative_sizes), max(small_sizes)),
    row.names = NULL, check.names = FALSE)
  selected_proteins <- data.frame(
    Protein = final_panel$proteins,
    Full_name = data$protein_full_names[match(final_panel$proteins, data$protein_names)],
    Coefficient = unname(final_panel$coefficients),
    Selection_frequency = vapply(final_panel$proteins, function(protein)
      mean(vapply(small_panels, function(panel) protein %in% panel, logical(1))), numeric(1)),
    row.names = NULL)
  selection_counts <- sort(table(unlist(small_panels)), decreasing = TRUE)
  selection_frequency <- data.frame(Protein = names(selection_counts),
                                     Frequency = as.numeric(selection_counts) / 25)
  results <- list(comparison = comparison, selected_proteins = selected_proteins,
                  historical_panel = historical_panel, historical_testing = historical_testing,
                  baseline_sizes = baseline_sizes, alternative_sizes = alternative_sizes,
                  accuracy_difference = accuracy_difference,
                  difference_interval = difference_interval,
                  baseline_metrics = baseline_metrics, alternative_metrics = alternative_metrics,
                  final_penalty_ratio = final_lasso$penalty_ratio,
                  final_lasso_proteins = final_lasso$proteins,
                  minimum_deviance_sizes = minimum_deviance_sizes,
                  small_sizes = small_sizes, small_metrics = small_metrics,
                  small_accuracy_difference = small_accuracy_difference,
                  small_difference_interval = small_difference_interval,
                  selection_frequency = selection_frequency)
  if (save_tables) {
    output_paths <- file.path(project_root, "results", "tables",
                              c("q4_model_comparison.csv", "q4_selected_proteins.csv"))
    # These two existing Q4 tables are explicitly authorized for updating.
    write.csv(comparison, output_paths[1], row.names = FALSE)
    write.csv(selected_proteins, output_paths[2], row.names = FALSE)
  }
  results
}

# 10. Display selection stability and paired accuracy differences
# Frequencies show how consistently the small procedure selects each protein.
# Repeat-level differences show partition variability without treating repeated
# predictions as independent subjects. The report embeds these plots directly.
plot_q4_results <- function(results) {
  frequency <- head(results$selection_frequency, 10)
  selection_plot <- ggplot2::ggplot(frequency,
      ggplot2::aes(x = Frequency, y = reorder(Protein, Frequency))) +
    ggplot2::geom_col(fill = "#236B8E", width = 0.65) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                               labels = function(value) paste0(100 * value, "%")) +
    ggplot2::labs(title = "Which proteins recur in the small panels?",
                  subtitle = "Training-only selection across 25 outer fits",
                  x = "Selection frequency (denominator: 25 fits)", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  differences <- data.frame(
    Repeat = factor(1:5),
    Difference = 100 * (results$small_metrics["Accuracy", ] - results$baseline_metrics["Accuracy", ]))
  difference_plot <- ggplot2::ggplot(differences, ggplot2::aes(x = Difference, y = Repeat)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey50") +
    ggplot2::geom_vline(xintercept = -5, color = "#AC503C", linetype = "dashed") +
    ggplot2::geom_point(size = 3, color = "#236B8E") +
    ggplot2::expand_limits(x = c(-5, 0)) +
    ggplot2::labs(title = "Small panel versus baseline in each repeat",
                  subtitle = "Dashed line: exploratory 5-percentage-point loss margin",
                  x = "Accuracy difference (small panel minus baseline, percentage points)",
                  y = "Validation repeat") +
    ggplot2::theme_minimal(base_size = 11)
  list(selection = selection_plot, difference = difference_plot)
}

# Run analysis ------------------------------------------------

q4_results <- run_q4_analysis(
  project_root = ".",
  save_tables = TRUE
)

# Model comparison
print(q4_results$comparison)

# Selected proteins
print(q4_results$selected_proteins)

# Figures
q4_plots <- plot_q4_results(q4_results)

print(q4_plots$selection)
print(q4_plots$difference)

ggplot2::ggsave("results/graphs/q4_fig-selection-1.png", q4_plots$selection,
                width = 7, height = 4.5, dpi = 300)
ggplot2::ggsave("results/graphs/q4_fig-difference-1.png", q4_plots$difference,
                width = 7, height = 4.5, dpi = 300)
