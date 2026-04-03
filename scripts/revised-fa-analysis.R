######################## LIBRARIES
library(dplyr)        # data manipulation
library(tidyr)        # data reshaping
library(readr)        # read CSV files
library(car)          # leveneTest()
library(ggplot2)      # data visualization
library(ggridges)     # geom_density_ridges()



######################## USER CONFIGURATION
# File paths -- copy from your 2-inspection-constants.R:
f_fromStep9  = "data/revised-fa-fromStep9.csv"
f_hypothesis = "data/revised-fa-hypothesis.csv"

# Grouping variables -- copy group_cols from your 2-inspection-constants.R:
group_cols = c("period")    # e.g. c("period") or c("Sex", "Dimension")

# Column names -- copy from your 2-inspection-constants.R:
col_individual = "ind"  # your individual ID column (e.g. "uid" or "ind")
col_side       = "side"  # your side column (e.g. "side")
col_side_left  = "L"  # your left side value (e.g. "L")
col_side_right = "R"  # your right side value (e.g. "R")

# Column names in f_hypothesis output:
trait_col  = "trait"       # trait column name
fa_col     = "fa10_noME"   # FA index column name
pct_me_col = "pct_me"      # percent ME column name



######################## DATA PREPARATION

# build_sd_long: computes signed side differences (R-L) per individual per
# trait from the final cleaned dataset. Used for distribution visualization.
# Grouping columns are taken from group_cols in fa-constants.R.
# Arguments:
#   f_data -- path to pipeline CSV (e.g. f_fromStep9)
# Returns: long dataframe with Ind, grouping columns, Trait, and diff (R-L)
build_sd_long = function(f_data) {
  dat = read_csv(f_data, show_col_types = FALSE) |>
    filter(!is.na(Value)) |>
    group_by(.data[[col_individual]], across(all_of(group_cols)), Trait, .data[[col_side]]) |>
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

  left_df  = dat |> filter(.data[[col_side]] == col_side_left) |>
    rename(val_l = Value) |> select(-all_of(col_side))
  right_df = dat |> filter(.data[[col_side]] == col_side_right) |>
    rename(val_r = Value) |> select(-all_of(col_side))

  join_cols = c(col_individual, group_cols, "Trait")
  inner_join(left_df, right_df, by = join_cols) |>
    mutate(diff = val_r - val_l) |>
    select(-val_l, -val_r)
}

fa_distribution = build_sd_long(f_fromStep9)


#### DISTRIBUTION VISUALIZATION
# Ridge/density plot of R-L distributions per analytical group.
# Useful visual confirmation that final distributions are approximately
# centred on zero before hypothesis testing.
# Adapt the aes() interaction() call to your group_cols.
# Example with group_cols = c("Sex", "Dimension"):
#   y = interaction(Sex, Dimension), fill = interaction(Sex, Dimension)

ggplot(fa_distribution,
       aes(x    = diff,
           y    = interaction(!!!syms(group_cols)),
           fill = interaction(!!!syms(group_cols)))) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8) +
  labs(x = "R-L", y = NULL, fill = NULL) +
  theme_classic()
ggsave("revised-fig-fa-distribution.png", width = 9, height = 3, dpi = 300)

rm(fa_distribution, build_sd_long)


#### LOAD HYPOTHESIS DATA (FA10a/b index -- for Levene's approach)
# Reads f_hypothesis output from 2-inspection.R.
# Rename or derive additional columns needed for your hypotheses here.
# group_cols columns are already present; add any trait sub-structure
# columns your hypotheses require (e.g. trait class, position, region).

fa = read_csv(f_hypothesis, show_col_types = FALSE)
names(fa) = tolower(names(fa))

# Rename FA index column to fa10 for convenience.
if (fa_col %in% names(fa)) fa = fa |> rename(fa10 = all_of(fa_col))

# Cast grouping columns and trait to factor.
for (col in c(tolower(group_cols), tolower(trait_col))) {
  if (col %in% names(fa)) fa[[col]] = as.factor(fa[[col]])
}

# Add any additional factor columns derived from your trait names here.
# These will be used as grouping variables in your hypothesis tests below.
# Example: derive a broad category from the first character of the trait name:
#   fa$category = as.factor(substr(fa[[trait_col]], 1, 1))

fa = droplevels(fa)



######################## FA10 vs PCT_ME (DIAGNOSTIC PLOT)
# Shows which traits had the most measurement error relative to FA.
# High pct_me means ME was a large fraction of apparent FA for that trait.
# Useful context when interpreting hypothesis test results.
# Adapt aes(color) to your trait column name.

ggplot(fa, aes(x = pct_me, y = fa10_nome, color=period)) +
  geom_point(size = 3) +
  scale_color_viridis_d() +
  theme_minimal() +
  scale_color_viridis_d(end = .6) +
  labs(x = "% ME", y = "FA10")
ggsave("revised-fig-fa-pct-me.png", width = 8, height = 4, dpi = 300)



######################## LEVENE
leveneTest(fa10_nome ~ period, data = fa, center = mean)

# not significant; hypothesis supported
# df 1,5 value 1.73 p 0.25

######################## TIDY
rm(list = ls())
gc()

