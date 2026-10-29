## ----message=FALSE, warning=FALSE---------------
library(tidyverse)
load("data/biomarker-clean.RData") 

#ls()
#dim(biomarker_clean)


## ----message=FALSE, warning=FALSE---------------
biomarker_raw <- read_csv("data/biomarker-raw.csv")

#dim(biomarker_raw)


## ----message=FALSE, warning=FALSE---------------
# get names
var_names <- read_csv("data/biomarker-raw.csv", 
                     col_names = F, 
                     n_max = 2, 
                     col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()


# read in data
biomarker_clean2 <- read_csv("data/biomarker-raw.csv", 
         skip = 2,
         col_select = -2L,
         col_names = c('group', 
                       'empty',
                       pull(var_names, abbreviation),
                       'ados'),
         na = c('-', '')) %>%
  filter(!is.na(group))
  # Removed log transform, center and scale, and trim


## -----------------------------------------------
# need to remove non-protein columns.
# These columns are the first and the last column that are not proteins.
protein_raw <- select(biomarker_clean2, -group, -ados)


set.seed(888)
protein_sample <- sample(protein_raw, 6)
#colnames(protein_sample)
#protein_sample

# level represents the amount of the specific protein in the person.
# protein_sample_raw and protein_sample_log were created for reusability.
protein_sample_raw <- protein_sample %>%
  pivot_longer(everything(), names_to = "protein", values_to = "level") %>%
  group_by(protein)

summary_raw <- protein_sample_raw %>%
  summarize(min = min(level),
            max = max(level),
            mean = mean(level),
            median = median(level),
            standard_deviation = sd(level),
            skewness = ((3*(mean(level) - median(level))) / sd(level))
            ) %>%
  arrange(desc(skewness))
print(summary_raw)


protein_sample_log <- protein_sample %>%
  pivot_longer(everything(), names_to = "protein", values_to = "level") %>%
  #log10() was used here because in preprocessing.r, log10() was used.
  mutate(level = log10(level)) %>%
  group_by(protein)

summary_log <- protein_sample_log %>%
    summarize(min = min(level),
              max = max(level),
              mean = mean(level),
              median = median(level),
              standard_deviation = sd(level),
              skewness = ((3*(mean(level) - median(level))) / sd(level))
              ) %>%
    arrange(desc(skewness))
print(summary_log)


# saves tables as csv files
write_csv(summary_raw, "results/tables/q1_summary_raw.csv")
write_csv(summary_log, "results/tables/q1_summary_log.csv")



## ----q1-raw-hist--------------------------------
library(ggplot2)
hist1 <- protein_sample_raw %>%
  ggplot(aes(x = level)) +
  geom_histogram(bins = 35, fill = 'skyblue', color = 'black') +
  # facet_wrap is used here instead of creating 6 different slides for protein_sample_raw
  # vars(protein), protein is a data frame, vars() helps find protein column
  # scales = free allows each histogram to have its own x and y-axis ranges.
  facet_wrap(vars(protein), scales = "free") +
  labs(title = "Raw Protein Levels",
       x = "Protein Level", y = "Count")
print(hist1)


## ----q1-log-hist--------------------------------
hist2 <- protein_sample_log %>%
  ggplot(aes(x = level)) +
  geom_histogram(bins = 35, fill = "skyblue", color = "black") +
  facet_wrap(vars(protein), scales = "free") +
  labs(title = "Log Protein Levels",
       x = "Protein Level", y = "Count")
print(hist2)

## ----message=FALSE, warning=FALSE---------------
# read in data
biomarker_clean3 <- read_csv("data/biomarker-raw.csv", 
         skip = 2,
         col_select = -2L,
         col_names = c('group', 
                       'empty',
                       pull(var_names, abbreviation),
                       'ados'),
         na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(.cols = -c(group, ados), 
                ~ (scale(log10(.x))[, 1]))) %>%
# I decided to add a column called subject, which gives a number to each individual to distinguish them from one another. This way we can sort them by who had the most outliers.
  mutate(subject = row_number()) %>%
  select(subject, group, everything(), -ados)
# removed ados, don't need for question 2


## -----------------------------------------------
biomarker3_long <- biomarker_clean3 %>%
  pivot_longer(cols = -c(subject, group), names_to = "protein", 
               values_to = "z")

outliers_subject <- biomarker3_long %>%
  group_by(subject, group) %>%
  summarize(
    num_z_outliers = sum(abs(z) > 3, na.rm = TRUE),
    num_z_outliers_perc = ((num_z_outliers / 1317) * 100)) %>%
  arrange(desc(num_z_outliers))
print(head(outliers_subject, 10))

print(median(outliers_subject$num_z_outliers))
print(sum(outliers_subject$num_z_outliers > 8.5))

# saves tables as csv files
write_csv(outliers_subject, "results/tables/q2_outliers_per_subject.csv")


## -----------------------------------------------
group_summary <- outliers_subject %>%
  group_by(group) %>%
  summarize(
    n_subjects = n(),
    mean_outliers = mean(num_z_outliers),
    median_outliers = median(num_z_outliers),
    sd_outliers = sd(num_z_outliers),
    # the max amount of outliers a single person had out of the group.
    max_outliers = max(num_z_outliers) )

print(group_summary)

write_csv(group_summary, "results/tables/q2_group_summary.csv")


## ----q2-boxplot---------------------------------
box1 <- outliers_subject %>%
  ggplot(aes(x = group, y = num_z_outliers, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = .1, alpha = .5, size = 1) +
  labs(title = "Outliers (|z| > 3) per Subject by Group",
    x = "Group", y = "Number of outliers per subject")
print(box1)
