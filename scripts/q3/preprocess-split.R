# Shared Q3 split: estimate preprocessing using training subjects only.
prepare_q3_split <- function() {
  headers <- readr::read_csv("data/biomarker-raw.csv", col_names = FALSE,
                            n_max = 2, show_col_types = FALSE)
  protein_names <- as.character(unlist(headers[2, 3:(ncol(headers) - 1)]))
  raw <- readr::read_csv("data/biomarker-raw.csv", skip = 2,
    col_names = c("group", "empty", protein_names, "ados"),
    na = c("-", ""), show_col_types = FALSE) %>%
    filter(!is.na(group)) %>% select(-empty)
  stopifnot(all(vapply(raw[protein_names], function(x)
    all(is.finite(x) & x > 0), logical(1))))
  set.seed(101422)
  split <- rsample::initial_split(raw, prop = 0.8, strata = group)
  train <- rsample::training(split)
  test <- rsample::testing(split)
  log_train <- log10(as.matrix(train[protein_names]))
  centers <- colMeans(log_train)
  scales <- apply(log_train, 2, sd)
  stopifnot(all(is.finite(scales) & scales > 0))
  transform <- function(data) {
    z <- sweep(sweep(log10(as.matrix(data[protein_names])), 2, centers),
               2, scales, "/")
    z[z > 3] <- 3
    z[z < -3] <- -3
    data[protein_names] <- as.data.frame(z)
    data
  }
  list(train = transform(train), test = transform(test))
}
