## ----setup, message=FALSE----------------------------------------------------
library(tidyverse)
library(infer)
library(randomForest)
library(rsample)
library(yardstick)

# Set the project root before the data-loading chunk runs
knitr::opts_knit$set(root.dir = normalizePath("../.."))


## ----------------------------------------------------------------------------
load("data/biomarker-clean.RData")

# Use the same split as Part (a)
set.seed(101422)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8, strata = group)
train_data <- training(biomarker_split)
test_data <- testing(biomarker_split)

table(train_data$group)
table(test_data$group)


## ----------------------------------------------------------------------------
# Rank proteins by their adjusted t-test p-values
test_fn <- function(df) {
  t_test(df, level ~ group, order = c("ASD", "TD"),
         alternative = "two-sided", var.equal = FALSE)
}

ttests_out <- train_data %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  mutate(p.adj = p.adjust(p_value, method = "BY")) %>%
  arrange(p.adj, p_value, protein)

# Break adjusted-p-value ties using raw p-values, then protein names
ttest_rank <- ttests_out$protein

# Fit random forest once so both comparisons use the same ranking
set.seed(101422)
rf_out <- randomForest(
  x = train_data %>% select(-group, -ados),
  y = factor(train_data$group),
  ntree = 1000, importance = TRUE
)

rf_rank <- rf_out$importance %>%
  as_tibble(rownames = "protein") %>%
  arrange(desc(MeanDecreaseGini), protein) %>%
  pull(protein)


## ----------------------------------------------------------------------------
# Repeat the same model fitting for top 10 and top 20
comparison <- tibble()
selected_proteins <- tibble()

for (k in c(10, 20)) {
  proteins_s1 <- head(ttest_rank, k)
  proteins_s2 <- head(rf_rank, k)
  panel <- intersect(proteins_s1, proteins_s2)

  selected_proteins <- bind_rows(selected_proteins,
                                tibble(top_k = k, protein = panel))

  model_train <- train_data %>%
    select(group, all_of(panel)) %>%
    mutate(class = as.integer(group == "ASD")) %>%
    select(-group)

  fit <- glm(class ~ ., data = model_train, family = binomial())
  prob <- predict(fit, newdata = test_data, type = "response")
  truth <- factor(test_data$group, levels = c("TD", "ASD"))
  predicted <- factor(ifelse(prob >= 0.5, "ASD", "TD"),
                      levels = c("TD", "ASD"))

  comparison <- bind_rows(comparison, tibble(
    top_k = k,
    panel_size = length(panel),
    accuracy = accuracy_vec(truth, predicted),
    sensitivity = sens_vec(truth, predicted, event_level = "second"),
    specificity = spec_vec(truth, predicted, event_level = "second"),
    AUROC = roc_auc_vec(truth, prob, event_level = "second")
  ))
}

knitr::kable(comparison, digits = 3)
knitr::kable(selected_proteins)

# Save the results for the group report
write.csv(comparison, "results/tables/q3b_comparison.csv", row.names = FALSE)
write.csv(selected_proteins, "results/tables/q3b_selected_proteins.csv",
          row.names = FALSE)


## ----------------------------------------------------------------------------


