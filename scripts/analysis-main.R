# PSTAT 197A - Module 1 Biomarker Analysis
# Team 9

# Setup 

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet)

# Run from the repository root. Create report output directories first.
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/graphs", recursive = TRUE, showWarnings = FALSE)

# Data preprocessing 

load("data/biomarker-clean.RData")


# Questions 1 & 2: Log transformation and outliers 

source(
  "scripts/q1 & q2/Christian Question 1 and 2 Analysis.R",
  local = new.env(parent = globalenv())
)


# Question 3: Methodological variations 

# (a) Training-only feature selection
source(
  "scripts/q3/Q3-Part-a.R",
  local = new.env(parent = globalenv())
)

# (b) Larger protein panels
source(
  "scripts/q3/Q3-Part-b.R",
  local = new.env(parent = globalenv())
)

# (c) Fuzzy intersection
source(
  "scripts/q3/Q3-Part-C.R",
  local = new.env(parent = globalenv())
)


# Question 4: Alternative classifier 

source(
  "scripts/q4/q4_analysis.R",
  local = new.env(parent = globalenv())
)
