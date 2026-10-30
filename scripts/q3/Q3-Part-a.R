# Question 3(a): Training/test split and protein selection

library(tidyverse) 
library(infer) 
library(randomForest) 
library(tidymodels) 
library(modelr) 
library(yardstick)

source("scripts/q3/preprocess-split.R", local = TRUE)
q3_partition <- prepare_q3_split()
train_data <- q3_partition$train
test_data <- q3_partition$test

# MULTIPLE TESTING (only on training)
# function to compute tests
test_fn <- function(.df) {
  t_test(
    .df,
    formula = level ~ group,
    order = c("ASD", "TD"),
    alternative = "two-sided",
    var.equal = FALSE
  )
}

# Run a t-test for each protein
ttests_out <- train_data %>%
  select(-ados) %>%
  pivot_longer(
    -group,
    names_to = "protein",
    values_to = "level"
  ) %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  mutate(p.adj = p.adjust(p_value, method = "BY")) %>%
  arrange(p.adj, p_value, protein)

# select top 10 significant proteins
proteins_s1 <- ttests_out %>%
  slice_head(n = 10) %>%
  pull(protein)


# RANDOM FOREST (only on training)

# store predictors and response separately using train_data
predictors <- train_data %>%
  select(-c(group, ados))

response <- train_data %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  arrange(desc(MeanDecreaseGini), protein) %>%
  slice_head(n = 10) %>%
  pull(protein)

# LOGISTIC REG & EVAL
# Find the overlapping proteins selected using training data only
proteins_sstar <- intersect(proteins_s1, proteins_s2)

# Check how many proteins were selected
length(proteins_s1)
length(proteins_s2)
length(proteins_sstar)

proteins_sstar

# Prepare training data
biomarker_sstar_train <- train_data %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = factor(group == "ASD",
                        levels = c(FALSE, TRUE))) %>%
  select(-group)

# Fit logistic regression model
fit <- glm(
  class ~ .,
  data = biomarker_sstar_train,
  family = binomial()
)


# Define classification metrics
class_metrics <- metric_set(
  sensitivity,
  specificity,
  accuracy
)

# Prepare testing data
biomarker_sstar_test <- test_data %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(
    class = factor(group == "ASD",
                   levels = c(FALSE, TRUE))
  ) %>%
  select(-group)

# Predict probabilities on the test set
biomarker_sstar_test <- biomarker_sstar_test %>%
  mutate(
    pred = predict(
      fit,
      newdata = biomarker_sstar_test,
      type = "response"
    ),
    estimate = factor(
      pred >= 0.5,
      levels = c(FALSE, TRUE)
    )
  )

# Sensitivity, specificity, and accuracy
classification_results <- class_metrics(
  data = biomarker_sstar_test,
  truth = class,
  estimate = estimate,
  event_level = "second"
)

# ROC AUC uses predicted probabilities, not class predictions
auc_result <- roc_auc(
  data = biomarker_sstar_test,
  truth = class,
  pred,
  event_level = "second"
)

classification_results
auc_result


# Calculate test-set metrics
cat("\n--- Q3(a): Test Set Classification Results ---\n")
print(classification_results)

cat("\n--- Q3(a): ROC AUC ---\n")
print(auc_result)



write.csv(bind_rows(classification_results, auc_result),
          "results/tables/q3a_metrics.csv", row.names = FALSE)
write.csv(data.frame(protein = proteins_sstar),
          "results/tables/q3a_selected_proteins.csv", row.names = FALSE)
