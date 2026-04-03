# 2-inspection-functions.R
# Palmer-Strobeck Fluctuating Asymmetry Pipeline
# Shared analytical functions
#
# PURPOSE:
#   All analytical and utility functions for the Palmer-Strobeck FA protocol.
#   No study-specific values are defined here. All constants (file paths,
#   trait variable lists, grouping variable definitions, outlier removals)
#   are defined in 2-inspection-constants.R.
#
# USAGE:
#   Source this file at the start of each session before running any
#   inspection or analysis script:
#     source("2-inspection-functions.R")
#     source("2-inspection-constants.R")
#
# GROUPING VARIABLE DESIGN:
#   Functions that filter or loop over analytical groups use group_cols
#   (defined in 2-inspection-constants.R) to identify which columns define groups.
#   The group key string (e.g. "female-breadth") is the values of those
#   columns joined by hyphens in the order they appear in group_cols.
#   load_group() parses the key back into filter conditions using group_cols,
#   making all functions robust to any number of grouping variables.
#
# FUNCTION OVERVIEW:
#   Utilities:         parse_group_key
#   Data reshaping:    pivot_reps_wide, avg_reps, build_sd_wide, build_rl
#   Group loading:     load_group, null_replicate
#   ME detection:      run_me_dixon, diagnose_me_flagged        (Steps 1-2)
#   FA detection:      run_fa_dixon, diagnose_fa_flagged        (Steps 3-5)
#   Mixed model ANOVA: run_univariates, run_mixedmodel          (Step 6)
#   ME variance:       lev_to_df                                (Step 7)
#   Size dependency:   run_step8a_scatters, run_step8b_spearman (Step 8)
#   Ideal FA tests:    run_descriptives, run_idealfa,
#                      assess_da_vs_fa4a                        (Step 9)
#   Normality:         run_normality_scatters
#   PS helpers:        bonferroni_sorted, derive_M, skew_se,
#                      kurt_se, moment_p
#
# STATISTICAL APPROACH (Palmer & Strobeck 1994, 2003):
#   ME outlier detection (Steps 1-2): Dixon's test on consecutive pair
#     differences of |R-L| values per individual per trait per group.
#   FA outlier detection (Step 5): Dixon's test on signed R-L values per
#     individual per trait per group.
#   Both ME and FA detection test against the group mean AND against zero.
#   Bonferroni correction is applied within each analytical group.
#   Data I/O uses readr (read_csv / write_csv) throughout.



# ============================================================
# UTILITIES
# ============================================================

# parse_group_key --------------------------------------------------------
parse_group_key = function(grp, grp_cols = group_cols) {
  # Parses a group key string into a named list of grouping variable values.
  # Used by load_group() and all functions that need to filter or label by
  # grouping variable values.
  #
  # The key is the grouping variable values joined by hyphens in the order
  # defined by group_cols (from 2-inspection-constants.R). For example:
  #   group_cols = c("Sex", "Dimension"), key = "female-breadth"
  #   returns: list(Sex = "female", Dimension = "breadth")
  #
  # Arguments:
  #   grp      -- group key string (e.g. "female-breadth", "Female-soldier")
  #   grp_cols -- character vector of grouping column names (default: group_cols)
  #
  # Returns: named list of grouping variable values, one element per group_col.
  parts = strsplit(grp, "-")[[1]]
  if (length(parts) != length(grp_cols)) {
    stop(sprintf(paste0(
      "Group key '%s' has %d parts but group_cols has %d elements.\n",
      "  group_cols: %s\n",
      "  Key parts:  %s\n",
      "  Keys must be grouping variable values joined by hyphens in the\n",
      "  same order as group_cols."
    ), grp, length(parts), length(grp_cols),
    paste(grp_cols, collapse = ", "),
    paste(parts, collapse = ", ")))
  }
  setNames(as.list(parts), grp_cols)
}



# ============================================================
# DATA RESHAPING
# ============================================================

# pivot_reps_wide --------------------------------------------------------
pivot_reps_wide = function(df, trait_cols) {
  # Pivots replicates from long to wide.
  # Arguments:
  #   df         -- long-format dataframe with Ind, Side, Rep, Trait, Value
  #   trait_cols -- character vector of trait names
  # Returns: dataframe with Ind, Side, Trait, rep1:repN columns
  df[[col_replicate]] = as.factor(df[[col_replicate]])
  long = pivot_longer(df |> filter(Trait %in% trait_cols),
                      cols = Value, names_to = "val_col", values_to = "val")
  wide = pivot_wider(long, names_from = all_of(col_replicate), values_from = val,
                     names_prefix = "rep")
  return(wide)
}

# consec_diffs -----------------------------------------------------------
consec_diffs = function(x) {
  # Computes consecutive differences from a numeric vector.
  # For ME outlier detection: input is |R-L| values per individual per
  # trait; output is consecutive pair differences across replicates.
  # Arguments:
  #   x -- numeric vector of length n (ordered by replicate number)
  # Returns: numeric vector of length n-1
  diff(x)
}

# avg_reps ---------------------------------------------------------------
avg_reps = function(df, trait_cols) {
  # Averages replicates per trait per individual-side row.
  # Works on long-format data (Trait/Value columns).
  # Arguments:
  #   df         -- long-format dataframe with Ind, Side, Rep, Trait, Value
  #   trait_cols -- character vector of trait names
  # Returns: dataframe with one row per Ind-Side, one column per trait
  df[[col_replicate]] = as.integer(as.character(df[[col_replicate]]))
  df |>
    filter(Trait %in% trait_cols) |>
    group_by(.data[[col_individual]], .data[[col_side]], Trait) |>
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = Trait, values_from = Value)
}

# build_sd_wide ----------------------------------------------------------
build_sd_wide = function(df, trait_cols, abs_val = FALSE) {
  # Returns wide dataframe of R-L (or |R-L|) side differences per trait,
  # averaged across replicates. Used by Steps 4 and 5.
  # Arguments:
  #   df         -- long-format dataframe with Ind, Side, Rep, Trait, Value
  #   trait_cols -- character vector of trait names
  #   abs_val    -- logical; if TRUE returns |R-L|, if FALSE returns signed R-L
  # Returns: wide dataframe with Ind and <trait>_sd columns
  df_avg   = avg_reps(df, trait_cols)
  left_df  = df_avg[df_avg[[col_side]] == col_side_left, ] |> dplyr::select(-all_of(col_side))
  right_df = df_avg[df_avg[[col_side]] == col_side_right, ] |> dplyr::select(-all_of(col_side))

  long_l = pivot_longer(left_df,  cols = all_of(trait_cols),
                        names_to = "var", values_to = "val_left")
  long_r = pivot_longer(right_df, cols = all_of(trait_cols),
                        names_to = "var", values_to = "val_right")

  sd_df      = merge(long_l, long_r, by = c(col_individual, "var"))
  sd_df      = na.omit(sd_df)
  sd_df$diff = sd_df$val_right - sd_df$val_left
  if (abs_val) sd_df$diff = abs(sd_df$diff)
  sd_df$val_left  = NULL
  sd_df$val_right = NULL

  sd_wide        = pivot_wider(sd_df, names_from = var, values_from = diff)
  names(sd_wide) = ifelse(names(sd_wide) == col_individual, col_individual,
                          paste0(names(sd_wide), "_sd"))
  return(sd_wide)
}

# build_rl ---------------------------------------------------------------
build_rl = function(df, trait_cols) {
  # Returns merged left/right replicate-averaged dataframe (long).
  # Used by Steps 3 and 4 for L vs R and cross-trait scatterplot generation.
  # Arguments:
  #   df         -- long-format dataframe with Ind, Side, Rep, Trait, Value
  #   trait_cols -- character vector of trait names
  # Returns: long dataframe with val_left and val_right per Ind-var pair
  df_avg   = avg_reps(df, trait_cols)
  left_df  = df_avg[df_avg[[col_side]] == col_side_left, ]
  right_df = df_avg[df_avg[[col_side]] == col_side_right, ]
  long_l   = pivot_longer(left_df,  cols = all_of(trait_cols),
                          names_to = "var", values_to = "val_left")
  long_r   = pivot_longer(right_df, cols = all_of(trait_cols),
                          names_to = "var", values_to = "val_right")
  rl       = merge(long_l, long_r, by = c(col_individual, "var"))
  rl$var   = as.factor(rl$var)
  return(rl)
}



# ============================================================
# GROUP LOADING
# ============================================================

# load_group -------------------------------------------------------------
load_group = function(f_data, grp, trait_cols, grp_cols = group_cols) {
  # Reads the pipeline CSV and filters to one analytical group.
  # Parses the group key using group_cols to determine which column values
  # to filter on. Works with any number of grouping variables.
  # Returns the filtered dataframe with factor columns coerced, and the
  # trait vector intersected with traits actually present and non-NA.
  #
  # Arguments:
  #   f_data     -- path to the pipeline CSV
  #   grp        -- group key string (e.g. "female-breadth")
  #   trait_cols -- character vector of expected trait names for this group
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: list with df (filtered dataframe) and trait_cols (intersected)
  grp_vals = parse_group_key(grp, grp_cols)

  df = read_csv(f_data, show_col_types = FALSE)

  for (col in names(grp_vals)) {
    df = df[df[[col]] %in% grp_vals[[col]], ]
  }

  df[[col_individual]] = as.factor(df[[col_individual]])
  df[[col_side]]       = as.factor(df[[col_side]])
  df[[col_replicate]]  = as.factor(df[[col_replicate]])
  df$Trait             = as.factor(df$Trait)

  trait_cols = intersect(trait_cols, unique(as.character(df$Trait)))
  trait_cols = trait_cols[sapply(trait_cols, function(t) {
    any(!is.na(df$Value[df$Trait == t]))
  })]

  list(df = df, trait_cols = trait_cols)
}

# null_replicate ---------------------------------------------------------
null_replicate = function(df, ind_id, trait, rep) {
  # Sets a specific replicate Value to NA for a given individual and trait.
  # Used for ME and FA outlier removal where the removal unit is the specific
  # offending replicate, not the whole individual or trait.
  # Operates on the long-format CSV: matches on Ind, Rep, and Trait columns.
  #
  # Arguments:
  #   df     -- long-format dataframe with Ind, Rep, Trait, Value columns
  #   ind_id -- individual ID string
  #   trait  -- trait name
  #   rep    -- integer replicate number to null
  #
  # Returns: dataframe with the specified Value set to NA
  df$Value[df[[col_individual]] == ind_id & df[[col_replicate]] == rep & df$Trait == trait] = NA
  return(df)
}



# ============================================================
# OUTLIER DETECTION (STEPS 1-2 AND 3-5)
# ============================================================

# dixon_q ----------------------------------------------------------------
dixon_q = function(vals, ref) {
  # Computes Dixon's Q statistic for the most extreme value in a sample
  # relative to a reference value. Two-tailed: tests both the minimum and
  # maximum values and returns the larger Q statistic and its p-value.
  # Uses Dixon's Q critical value approximation via the outliers package.
  #
  # Arguments:
  #   vals -- numeric vector (NAs removed before calling)
  #   ref  -- reference value to test against (group mean or zero)
  #
  # Returns: list with Q statistic, p-value (two-tailed), direction
  #          ("high" or "low"), and the extreme value
  n = length(vals)
  if (n < 3) return(list(Q = NA_real_, p = NA_real_, direction = NA_character_,
                         extreme_val = NA_real_))

  shifted = vals - ref
  rng     = diff(range(shifted))
  if (rng == 0) return(list(Q = NA_real_, p = NA_real_, direction = NA_character_,
                            extreme_val = NA_real_))

  sorted = sort(shifted)
  Q_low  = (sorted[2]   - sorted[1])   / rng
  Q_high = (sorted[n]   - sorted[n-1]) / rng

  if (Q_high >= Q_low) {
    Q         = Q_high
    direction = "high"
    extreme   = vals[which.max(shifted)]
  } else {
    Q         = Q_low
    direction = "low"
    extreme   = vals[which.min(shifted)]
  }

  p_result = tryCatch(
    outliers::dixon.test(vals, type = 0, opposite = (direction == "low"),
                         two.sided = TRUE)$p.value,
    error = function(e) NA_real_
  )

  list(Q = Q, p = p_result, direction = direction, extreme_val = extreme)
}

# run_me_dixon -----------------------------------------------------------
run_me_dixon = function(df, trait_cols, grp, grp_cols = group_cols) {
  # Step 2 ME outlier detection using Dixon's test.
  # Detects individuals whose replicate-to-replicate consistency in side
  # differences is extreme relative to the group.
  #
  # For each individual x trait x group:
  #   1. Computes |R-L| for each replicate (one absolute side difference per rep).
  #   2. Computes consecutive differences between those |R-L| values across
  #      replicates -- measures how much the side difference varies between
  #      replicate pairs.
  #   3. Runs Dixon's two-tailed test on those consecutive replicate differences
  #      against the group mean AND against zero.
  # Bonferroni correction applied within each group separately.
  #
  # Arguments:
  #   df         -- long-format dataframe for one group (from load_group)
  #   trait_cols -- character vector of trait names for this group
  #   grp        -- group key string
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with grouping columns, trait, and test results
  cat("\n---", grp, "---\n")
  grp_vals = parse_group_key(grp, grp_cols)

  rows = list()

  for (v in trait_cols) {
    cat("\n", toupper(v), "\n")

    trait_df = df[df$Trait == v, c(col_individual, col_replicate, col_side, "Value")]
    trait_df = trait_df[!is.na(trait_df$Value), ]
    trait_df[[col_replicate]] = as.integer(as.character(trait_df[[col_replicate]]))

    inds = unique(trait_df[[col_individual]])
    ind_sd_list = list()
    for (ind in inds) {
      sub  = trait_df[trait_df[[col_individual]] == ind, ]
      reps = sort(unique(sub[[col_replicate]]))
      sd_vals = vapply(reps, function(r) {
        r_val = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_right]
        l_val = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_left]
        if (length(r_val) == 1 && length(l_val) == 1 &&
            !is.na(r_val) && !is.na(l_val)) {
          abs(r_val - l_val)
        } else {
          NA_real_
        }
      }, numeric(1))
      names(sd_vals) = reps
      ind_sd_list[[as.character(ind)]] = sd_vals
    }

    all_consec = unlist(lapply(ind_sd_list, function(x) {
      x_clean = x[!is.na(x)]
      if (length(x_clean) >= 2) consec_diffs(x_clean) else numeric(0)
    }))
    group_mean = mean(all_consec, na.rm = TRUE)
    group_sd   = sd(all_consec,   na.rm = TRUE)
    cat(sprintf("  Group mean consec diff: %.6f  SD: %.6f  n: %d\n",
                group_mean, group_sd, length(all_consec)))

    for (ind in inds) {
      sd_vals  = ind_sd_list[[as.character(ind)]]
      sd_clean = sd_vals[!is.na(sd_vals)]
      if (length(sd_clean) < 2) next
      cd = consec_diffs(sd_clean)
      if (length(cd) < 3) next

      res_mean = dixon_q(cd, ref = group_mean)
      res_zero = dixon_q(cd, ref = 0)

      row = as.data.frame(grp_vals, stringsAsFactors = FALSE)
      row$trait       = v
      row[[col_individual]] = as.character(ind)
      row$n_consec    = length(cd)
      row$group_mean  = group_mean
      row$group_sd    = group_sd
      row$extreme_val = res_mean$extreme_val
      row$direction   = res_mean$direction
      row$Q_vs_mean   = res_mean$Q
      row$p_vs_mean   = res_mean$p
      row$Q_vs_zero   = res_zero$Q
      row$p_vs_zero   = res_zero$p
      rows[[length(rows) + 1]] = row
    }
  }

  out = do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) return(data.frame())

  out = out[order(out$p_vs_mean), ]
  bonf_mean     = bonferroni_sorted(out$p_vs_mean)
  out$bonf_mean = bonf_mean$bonf
  out$rank_mean = bonf_mean$rank
  out$sig_mean  = bonf_mean$alpha

  out = out[order(out$p_vs_zero), ]
  bonf_zero     = bonferroni_sorted(out$p_vs_zero)
  out$bonf_zero = bonf_zero$bonf
  out$rank_zero = bonf_zero$rank
  out$sig_zero  = bonf_zero$alpha

  out = out[order(out[[grp_cols[1]]], out$trait, out[[col_individual]]), ]
  rownames(out) = NULL

  print_cols = c(grp_cols, "trait", col_individual, "n_consec", "group_mean",
                 "extreme_val", "Q_vs_mean", "p_vs_mean", "sig_mean",
                 "Q_vs_zero", "p_vs_zero", "sig_zero")
  print(out[, intersect(print_cols, names(out))])
  out
}

# run_fa_dixon -----------------------------------------------------------
run_fa_dixon = function(df, trait_cols, grp, grp_cols = group_cols) {
  # Step 5 FA outlier detection using Dixon's test.
  # Detects individuals whose side difference itself is extreme relative
  # to the group -- distinct from ME outliers (Steps 1-2) whose
  # replicate-to-replicate variation in side difference is extreme.
  #
  # For each individual x trait x group:
  #   Computes signed R-L for each replicate (one side difference per rep).
  #   Runs Dixon's two-tailed test on those signed R-L values against the
  #   group mean R-L AND against zero.
  # Bonferroni correction applied within each group separately.
  #
  # Arguments:
  #   df         -- long-format dataframe for one group (from load_group)
  #   trait_cols -- character vector of trait names for this group
  #   grp        -- group key string
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with grouping columns, trait, and test results
  cat("\n---", grp, "---\n")
  grp_vals = parse_group_key(grp, grp_cols)

  rows = list()

  for (v in trait_cols) {
    cat("\n", toupper(v), "\n")

    trait_df = df[df$Trait == v, c(col_individual, col_replicate, col_side, "Value")]
    trait_df = trait_df[!is.na(trait_df$Value), ]
    trait_df[[col_replicate]] = as.integer(as.character(trait_df[[col_replicate]]))

    inds = unique(trait_df[[col_individual]])
    ind_rl_list = list()
    for (ind in inds) {
      sub  = trait_df[trait_df[[col_individual]] == ind, ]
      reps = sort(unique(sub[[col_replicate]]))
      rl_vals = vapply(reps, function(r) {
        r_val = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_right]
        l_val = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_left]
        if (length(r_val) == 1 && length(l_val) == 1 &&
            !is.na(r_val) && !is.na(l_val)) {
          r_val - l_val
        } else {
          NA_real_
        }
      }, numeric(1))
      names(rl_vals) = reps
      ind_rl_list[[as.character(ind)]] = rl_vals
    }

    all_rl     = unlist(lapply(ind_rl_list, function(x) x[!is.na(x)]))
    group_mean = mean(all_rl, na.rm = TRUE)
    group_sd   = sd(all_rl,   na.rm = TRUE)
    cat(sprintf("  Group mean R-L: %.6f  SD: %.6f  n: %d\n",
                group_mean, group_sd, length(all_rl)))

    for (ind in inds) {
      rl_vals  = ind_rl_list[[as.character(ind)]]
      rl_clean = rl_vals[!is.na(rl_vals)]
      if (length(rl_clean) < 3) next

      res_mean = dixon_q(rl_clean, ref = group_mean)
      res_zero = dixon_q(rl_clean, ref = 0)

      row = as.data.frame(grp_vals, stringsAsFactors = FALSE)
      row$trait       = v
      row[[col_individual]] = as.character(ind)
      row$n_reps      = length(rl_clean)
      row$group_mean  = group_mean
      row$group_sd    = group_sd
      row$extreme_val = res_mean$extreme_val
      row$direction   = res_mean$direction
      row$Q_vs_mean   = res_mean$Q
      row$p_vs_mean   = res_mean$p
      row$Q_vs_zero   = res_zero$Q
      row$p_vs_zero   = res_zero$p
      rows[[length(rows) + 1]] = row
    }
  }

  out = do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) return(data.frame())

  out = out[order(out$p_vs_mean), ]
  bonf_mean     = bonferroni_sorted(out$p_vs_mean)
  out$bonf_mean = bonf_mean$bonf
  out$rank_mean = bonf_mean$rank
  out$sig_mean  = bonf_mean$alpha

  out = out[order(out$p_vs_zero), ]
  bonf_zero     = bonferroni_sorted(out$p_vs_zero)
  out$bonf_zero = bonf_zero$bonf
  out$rank_zero = bonf_zero$rank
  out$sig_zero  = bonf_zero$alpha

  out = out[order(out[[grp_cols[1]]], out$trait, out[[col_individual]]), ]
  rownames(out) = NULL

  print_cols = c(grp_cols, "trait", col_individual, "n_reps", "group_mean",
                 "extreme_val", "Q_vs_mean", "p_vs_mean", "sig_mean",
                 "Q_vs_zero", "p_vs_zero", "sig_zero")
  print(out[, intersect(print_cols, names(out))])
  out
}



# ============================================================
# DIAGNOSTIC
# ============================================================

# diagnose_fa_flagged ----------------------------------------------------
diagnose_fa_flagged = function(dixon_out, raw_groups, grp_cols = group_cols) {
  # For each flagged individual x trait from run_fa_dixon, computes the
  # signed R-L values per replicate and identifies which replicate produced
  # the extreme value.
  #
  # Arguments:
  #   dixon_out  -- dataframe returned by run_fa_dixon (full output, all rows)
  #   raw_groups -- named list of long-format group dataframes, keyed by group key
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with one row per replicate per flagged individual x trait
  flagged = dixon_out[dixon_out$sig_mean == "*" | dixon_out$sig_zero == "*", ]
  if (nrow(flagged) == 0) {
    cat("No flagged cases to diagnose.\n")
    return(data.frame())
  }

  rows = list()

  for (i in seq_len(nrow(flagged))) {
    grp_key = paste(sapply(grp_cols, function(col) flagged[[col]][i]),
                    collapse = "-")
    trait = flagged$trait[i]
    ind   = flagged[[col_individual]][i]

    df = raw_groups[[grp_key]]
    if (is.null(df)) {
      cat("Warning: no raw data found for group", grp_key, "\n")
      next
    }

    sub = df[df[[col_individual]] == ind & df$Trait == trait, c(col_replicate, col_side, "Value")]
    sub = sub[!is.na(sub$Value), ]
    sub[[col_replicate]] = as.integer(as.character(sub[[col_replicate]]))
    reps = sort(unique(sub[[col_replicate]]))

    rl_vals = vapply(reps, function(r) {
      rv = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_right]
      lv = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_left]
      if (length(rv) == 1 && length(lv) == 1) rv - lv else NA_real_
    }, numeric(1))
    names(rl_vals) = reps
    rl_vals = rl_vals[!is.na(rl_vals)]

    gm         = flagged$group_mean[i]
    deviations = abs(rl_vals - gm)

    for (j in seq_along(rl_vals)) {
      row = as.data.frame(
        setNames(as.list(sapply(grp_cols, function(col) flagged[[col]][i])),
                 grp_cols),
        stringsAsFactors = FALSE)
      row$trait      = trait
      row[[col_individual]] = ind
      row$rep        = as.integer(names(rl_vals)[j])
      row$rl_val     = rl_vals[j]
      row$group_mean = gm
      row$deviation  = deviations[j]
      row$is_extreme = deviations[j] == max(deviations)
      rows[[length(rows) + 1]] = row
    }
  }

  out = do.call(rbind, rows)
  out = out[order(out[[grp_cols[1]]], out$trait, out[[col_individual]]), ]
  out = out[order(!out$is_extreme), ]
  rownames(out) = NULL
  print(out)
  out
}

# diagnose_me_flagged ----------------------------------------------------
diagnose_me_flagged = function(dixon_out, raw_groups, grp_cols = group_cols) {
  # For each flagged individual x trait from run_me_dixon, computes |R-L|
  # values and consecutive pair differences with replicate pair labels,
  # and identifies which pair produced the extreme value.
  #
  # Arguments:
  #   dixon_out  -- dataframe returned by run_me_dixon (full output, all rows)
  #   raw_groups -- named list of long-format group dataframes, keyed by group key
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with one row per consecutive pair per flagged individual x trait
  flagged = dixon_out[dixon_out$sig_mean == "*" | dixon_out$sig_zero == "*", ]
  if (nrow(flagged) == 0) {
    cat("No flagged cases to diagnose.\n")
    return(data.frame())
  }

  rows = list()

  for (i in seq_len(nrow(flagged))) {
    grp_key = paste(sapply(grp_cols, function(col) flagged[[col]][i]),
                    collapse = "-")
    trait = flagged$trait[i]
    ind   = flagged[[col_individual]][i]

    df = raw_groups[[grp_key]]
    if (is.null(df)) {
      cat("Warning: no raw data found for group", grp_key, "\n")
      next
    }

    sub = df[df[[col_individual]] == ind & df$Trait == trait, c(col_replicate, col_side, "Value")]
    sub = sub[!is.na(sub$Value), ]
    sub[[col_replicate]] = as.integer(as.character(sub[[col_replicate]]))
    reps = sort(unique(sub[[col_replicate]]))

    abs_rl = vapply(reps, function(r) {
      rv = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_right]
      lv = sub$Value[sub[[col_replicate]] == r & sub[[col_side]] == col_side_left]
      if (length(rv) == 1 && length(lv) == 1) abs(rv - lv) else NA_real_
    }, numeric(1))
    names(abs_rl) = reps
    abs_rl = abs_rl[!is.na(abs_rl)]

    if (length(abs_rl) < 2) next

    cd          = consec_diffs(abs_rl)
    rep_keys    = names(abs_rl)
    pair_labels = paste0("M", rep_keys[-1], "-M", rep_keys[-length(rep_keys)])

    for (j in seq_along(cd)) {
      row = as.data.frame(
        setNames(as.list(sapply(grp_cols, function(col) flagged[[col]][i])),
                 grp_cols),
        stringsAsFactors = FALSE)
      row$trait       = trait
      row[[col_individual]] = ind
      row$pair        = pair_labels[j]
      row$abs_rl_a    = abs_rl[j]
      row$abs_rl_b    = abs_rl[j + 1]
      row$consec_diff = cd[j]
      row$is_extreme  = abs(cd[j] - flagged$group_mean[i]) ==
        max(abs(cd - flagged$group_mean[i]))
      rows[[length(rows) + 1]] = row
    }
  }

  out = do.call(rbind, rows)
  out = out[order(out[[grp_cols[1]]], out$trait, out[[col_individual]]), ]
  out = out[order(!out$is_extreme), ]
  rownames(out) = NULL
  print(out)
  out
}



# ============================================================
# PS WORKSHEET HELPERS
# ============================================================

# bonferroni_sorted ------------------------------------------------------
bonferroni_sorted = function(p_vals, alpha = 0.05) {
  # Sequential Bonferroni correction (Rice 1989).
  # Step-down procedure: sort p values ascending, threshold at
  # rank i = alpha/(n-i+1). Stop at first non-significant result;
  # all subsequent rows are also non-significant.
  # Data must be sorted by p-value ascending before calling.
  n    = length(p_vals)
  rank = seq_len(n)
  bonf = alpha / (n - rank + 1)
  sig  = logical(n)
  for (i in seq_len(n)) {
    if (!is.na(p_vals[i]) && p_vals[i] < bonf[i]) {
      sig[i] = TRUE
    } else {
      break
    }
  }
  data.frame(rank  = rank,
             bonf  = bonf,
             alpha = ifelse(sig, "*", "n.s."),
             stringsAsFactors = FALSE)
}

# derive_M ---------------------------------------------------------------
derive_M = function(df_me, df_sides, df_ind) {
  # Derives M (number of replicates) from ANOVA df values.
  # Excel equivalent: =(df_me / ((df_sides + 1) * (df_ind + 1))) + 1
  # fa10_m is always derived from df -- never hardcoded.
  (df_me / ((df_sides + 1) * (df_ind + 1))) + 1
}

# skew_se ----------------------------------------------------------------
skew_se = function(n) {
  # Exact standard error of skewness.
  # Excel: =SQRT((6*n*(n-1)) / ((n-2)*(n+1)*(n+3)))
  # Warning: NOT sqrt(6/n) -- use this exact formula.
  sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))
}

# kurt_se ----------------------------------------------------------------
kurt_se = function(n) {
  # Exact standard error of kurtosis.
  # Excel: =SQRT((24*n*(n-1)^2) / ((n-3)*(n-2)*(n+3)*(n+5)))
  # Warning: NOT sqrt(24/n) -- use this exact formula.
  sqrt((24 * n * (n - 1)^2) / ((n - 3) * (n - 2) * (n + 3) * (n + 5)))
}

# moment_p ---------------------------------------------------------------
moment_p = function(ts) {
  # P-value for skew/kurt tests using df=10000 approximation.
  # Worksheet approximation for infinity degrees of freedom.
  pt(-abs(ts), df = 10000) * 2
}

# lev_to_df --------------------------------------------------------------
lev_to_df = function(lev, term) {
  # Converts leveneTest output to a labelled dataframe for export.
  # Arguments:
  #   lev  -- output of leveneTest()
  #   term -- character label for the tested term
  # Returns: dataframe with term and row columns added
  df      = as.data.frame(lev)
  df$term = term
  df$row  = rownames(df)
  df
}



# ============================================================
# ANOVA FUNCTIONS (STEP 6)
# ============================================================

# run_univariates --------------------------------------------------------
run_univariates = function(df, trait_cols, grp, grp_cols = group_cols) {
  # Runs Sides x Individuals ANOVA per trait for one group (Step 6).
  # Pivots to wide format internally before running ANOVA.
  # This step tests whether FA variance significantly exceeds ME variance
  # per trait. Results feed into run_mixedmodel().
  #
  # Arguments:
  #   df         -- long-format dataframe for one group (from load_group)
  #   trait_cols -- character vector of trait names for this group
  #   grp        -- group key string
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with grouping columns, trait, and ANOVA MS/df values
  cat("\n\n---", grp, "---\n")
  grp_vals = parse_group_key(grp, grp_cols)

  df_wide = df |>
    filter(Trait %in% trait_cols) |>
    select(all_of(c(col_individual, col_side, col_replicate)), Trait, Value) |>
    pivot_wider(names_from = Trait, values_from = Value)

  df_wide[[col_individual]] = as.factor(df_wide[[col_individual]])
  df_wide[[col_side]]       = as.factor(df_wide[[col_side]])
  df_wide[[col_replicate]]  = as.factor(df_wide[[col_replicate]])

  trait_cols = intersect(trait_cols, names(df_wide))

  cat("\nAOV Univariates:\n")
  aov_rows = list()
  for (v in trait_cols) {
    cat("\n", toupper(v), "\n")
    tryCatch({
      fit = suppressWarnings(aov(
        as.formula(paste(v, "~", col_side, "+ ", col_side, "*", col_individual,
                         "+ Error(", col_individual, "/(", col_side, "*", col_individual, "))")),
        data = df_wide
      ))
      s = summary(fit)
      print(s)

      s_ind  = s[[paste0("Error: ", col_individual)]][[1]]
      s_side = s[[paste0("Error: ", col_individual, ":", col_side)]][[1]]
      s_with = s[["Error: Within"]][[1]]
      rownames(s_ind)  = trimws(rownames(s_ind))
      rownames(s_side) = trimws(rownames(s_side))

      ms_I  = s_ind[col_individual,                          "Mean Sq"]
      df_I  = s_ind[col_individual,                          "Df"]
      ms_S  = s_side[col_side,                               "Mean Sq"]
      df_S  = s_side[col_side,                               "Df"]
      ms_SI = s_side[paste0(col_side, ":", col_individual),  "Mean Sq"]
      df_SI = s_side[paste0(col_side, ":", col_individual),  "Df"]
      ms_ME = s_with["Residuals", "Mean Sq"]
      df_ME = s_with["Residuals", "Df"]

      row = as.data.frame(grp_vals, stringsAsFactors = FALSE)
      row$trait = v
      row$MS_S  = ms_S;  row$df_S  = df_S
      row$MS_I  = ms_I;  row$df_I  = df_I
      row$MS_SI = ms_SI; row$df_SI = df_SI
      row$MS_ME = ms_ME; row$df_ME = df_ME
      aov_rows[[v]] = row
    }, error = function(e) {
      cat("  Skipping", v, "-- insufficient data for ANOVA:", conditionMessage(e), "\n")
    })
  }
  aov_rows = Filter(Negate(is.null), aov_rows)
  if (length(aov_rows) == 0) return(data.frame())
  aov_df = do.call(rbind, aov_rows)
  rownames(aov_df) = NULL
  aov_df
}

# run_mixedmodel ---------------------------------------------------------
run_mixedmodel = function(grp, aov_inputs) {
  # Computes FA > ME test with Bonferroni correction from run_univariates
  # output (Step 6). Tests F = MS_SxI / MS_ME per trait per group.
  # Also derives FA10 (ME-corrected FA index) and pct_me.
  #
  # Arguments:
  #   grp        -- group key string (key into aov_inputs)
  #   aov_inputs -- named list of dataframes from run_univariates
  #
  # Returns: dataframe sorted by fa_me_p with Bonferroni columns,
  #          pct_me, fa10_noME, and fa10_df
  ws = aov_inputs[[grp]]
  n  = nrow(ws)

  cat(sprintf("\n=== %s  (n traits = %d) ===\n", grp, n))

  sides_f = numeric(n); sides_p = numeric(n)
  ind_f   = numeric(n); ind_p   = numeric(n)
  fa_me_f = numeric(n); fa_me_p = numeric(n)
  me_pct     = numeric(n)
  fa10_m     = numeric(n)
  fa10_value = numeric(n)
  fa10_df    = numeric(n)

  for (i in seq_len(n)) {
    ms_s  = ws$MS_S[i];  df_s  = ws$df_S[i]
    ms_i  = ws$MS_I[i];  df_i  = ws$df_I[i]
    ms_si = ws$MS_SI[i]; df_si = ws$df_SI[i]
    ms_me = ws$MS_ME[i]; df_me = ws$df_ME[i]

    sides_f[i] = ms_s  / ms_si
    sides_p[i] = pf(sides_f[i], df_s,  df_si, lower.tail = FALSE)

    ind_f[i]   = ms_i  / ms_si
    ind_p[i]   = pf(ind_f[i],   df_i,  df_si, lower.tail = FALSE)

    fa_me_f[i] = ms_si / ms_me
    fa_me_p[i] = pf(fa_me_f[i], df_si, df_me, lower.tail = FALSE)

    me_pct[i]     = 100 * ms_me / ms_si
    fa10_m[i]     = derive_M(df_me, df_s, df_i)
    fa10_value[i] = max(ms_si - ms_me, 0) / fa10_m[i]

    if (is.na(ms_si) || is.na(ms_me)) {
      fa10_df[i] = NA_real_
    } else if (ms_si > ms_me) {
      num        = (ms_si - ms_me)^2
      denom      = (ms_si^2 / df_si) + (ms_me^2 / df_me)
      fa10_df[i] = num / denom
    } else {
      fa10_df[i] = NA_real_
    }
  }

  grp_cols_present = intersect(group_cols, names(ws))
  anova_base = ws[, c(grp_cols_present, "trait")]
  anova_base$sides_ms = ws$MS_S;  anova_base$sides_df = ws$df_S
  anova_base$ind_ms   = ws$MS_I;  anova_base$ind_df   = ws$df_I
  anova_base$sxi_ms   = ws$MS_SI; anova_base$sxi_df   = ws$df_SI
  anova_base$me_ms    = ws$MS_ME; anova_base$me_df    = ws$df_ME
  anova_base$sides_f  = sides_f;  anova_base$sides_p  = sides_p
  anova_base$ind_f    = ind_f;    anova_base$ind_p    = ind_p
  anova_base$fa_me_f  = fa_me_f;  anova_base$fa_me_p  = fa_me_p

  anova_sorted = anova_base[order(anova_base$fa_me_p), ]
  bonf_anova   = bonferroni_sorted(anova_sorted$fa_me_p)
  anova_out    = cbind(anova_sorted,
                       anova_bonf  = bonf_anova$bonf,
                       anova_rank  = bonf_anova$rank,
                       anova_alpha = bonf_anova$alpha,
                       pct_me      = me_pct[order(anova_base$fa_me_p)],
                       n_reps      = fa10_m[order(anova_base$fa_me_p)],
                       fa10_noME   = fa10_value[order(anova_base$fa_me_p)],
                       fa10_df     = fa10_df[order(anova_base$fa_me_p)])

  cat("\n  ANOVA sort (by fa_me_p):\n")
  print_cols = c(grp_cols_present, "trait", "fa_me_f", "fa_me_p",
                 "anova_bonf", "anova_rank", "anova_alpha",
                 "fa10_noME", "fa10_df", "pct_me")
  print(anova_out[, intersect(print_cols, names(anova_out))])
  anova_out
}



# ============================================================
# DESCRIPTIVES AND IDEAL FA (STEP 9)
# ============================================================

# run_descriptives -------------------------------------------------------
run_descriptives = function(df, trait_cols, grp, grp_cols = group_cols) {
  # Computes skew, kurtosis, and side diff/avg descriptives on post-Step8
  # data. Input for run_idealfa() (Step 9).
  #
  # Arguments:
  #   df         -- long-format dataframe for one group (from load_group)
  #   trait_cols -- character vector of trait names for this group
  #   grp        -- group key string
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with grouping columns, trait, and descriptive stats
  cat("\n\n---", grp, "---\n")
  grp_vals = parse_group_key(grp, grp_cols)

  df_avg   = avg_reps(df, trait_cols)
  left_df  = df_avg[df_avg[[col_side]] == col_side_left, ] |> dplyr::select(-all_of(col_side))
  right_df = df_avg[df_avg[[col_side]] == col_side_right, ] |> dplyr::select(-all_of(col_side))

  long_l = pivot_longer(left_df,  cols = all_of(trait_cols),
                        names_to = "var", values_to = "left")
  long_r = pivot_longer(right_df, cols = all_of(trait_cols),
                        names_to = "var", values_to = "right")

  desc_df       = merge(long_l, long_r, by = c(col_individual, "var"))
  desc_df       = na.omit(desc_df)

  diff_df       = desc_df
  diff_df$diff  = diff_df$right - diff_df$left
  diff_df$left  = NULL
  diff_df$right = NULL
  diff_wide     = pivot_wider(diff_df, names_from = var, values_from = diff)
  names(diff_wide) = ifelse(names(diff_wide) == col_individual, col_individual,
                            paste0(names(diff_wide), "_sd"))

  av_df         = desc_df
  av_df$average = (av_df$right + av_df$left) / 2
  av_df$left    = NULL
  av_df$right   = NULL
  av_wide       = pivot_wider(av_df, names_from = var, values_from = average)

  cat("\nSkew and Kurtosis on R-L differences:\n")
  sd_cols   = grep("_sd$", names(diff_wide), value = TRUE)
  skew_rows = list()
  for (v in sd_cols) {
    trait = sub("_sd$", "", v)
    vals  = diff_wide[[v]]; vals = vals[!is.na(vals)]
    cat("\n", toupper(trait), "\n")
    sk = Skew(vals, na.rm = TRUE)
    ku = Kurt(vals, na.rm = TRUE, method = 1)
    cat("  Skew:    ", sk, "\n")
    cat("  Kurtosis:", ku, "\n")
    skew_rows[[trait]] = data.frame(trait = trait, skew = sk, kurt = ku,
                                    stringsAsFactors = FALSE)
  }
  skew_df = do.call(rbind, skew_rows)
  rownames(skew_df) = NULL

  cat("\nSide Difference Descriptives (R-L):\n")
  desc_rows = list()
  for (v in sd_cols) {
    trait   = sub("_sd$", "", v)
    vals_sd = diff_wide[[v]];   vals_sd = vals_sd[!is.na(vals_sd)]
    vals_av = av_wide[[trait]]; vals_av = vals_av[!is.na(vals_av)]
    cat("\n", toupper(trait), "\n")
    print(describe(vals_sd, type = 2))
    desc_rows[[trait]] = data.frame(
      trait        = trait,
      n            = length(vals_sd),
      mean_diff    = mean(vals_sd),
      mean_diff_se = sd(vals_sd) / sqrt(length(vals_sd)),
      mean_av      = mean(vals_av),
      mean_av_se   = sd(vals_av) / sqrt(length(vals_av)),
      stringsAsFactors = FALSE
    )
  }

  cat("\nSide Average Descriptives ((R+L)/2):\n")
  for (v in setdiff(names(av_wide), col_individual)) {
    cat("\n", toupper(v), "\n")
    print(describe(av_wide[[v]], type = 2))
  }

  desc_out = do.call(rbind, desc_rows)
  rownames(desc_out) = NULL

  result = merge(skew_df, desc_out, by = "trait", all = TRUE)
  for (col in names(grp_vals)) result[[col]] = grp_vals[[col]]
  result
}

# run_idealfa ------------------------------------------------------------
run_idealfa = function(grp, desc_inputs) {
  # Computes DA, skew, and kurtosis blocks with Bonferroni correction
  # (Step 9). Input is from run_descriptives on post-Step8 data.
  # Bonferroni correction applied within each group separately.
  #
  # Arguments:
  #   grp         -- group key string (key into desc_inputs)
  #   desc_inputs -- named list of dataframes from run_descriptives
  #
  # Returns: list with da, skew, kurt blocks (each independently Bonferroni-sorted)
  ws = desc_inputs[[grp]]
  n  = nrow(ws)

  cat(sprintf("\n=== %s  (n traits = %d) ===\n", grp, n))

  da_ts  = ws$mean_diff / ws$mean_diff_se
  da_p   = pt(-abs(da_ts), df = ws$n - 1) * 2

  se_skew = skew_se(ws$n)
  skew_ts = ws$skew / se_skew
  skew_p  = moment_p(skew_ts)

  se_kurt = kurt_se(ws$n)
  kurt_ts = ws$kurt / se_kurt
  kurt_p  = moment_p(kurt_ts)

  grp_cols_present = intersect(group_cols, names(ws))

  da_base = ws[, c(grp_cols_present, "trait", "mean_av", "mean_av_se",
                   "mean_diff", "mean_diff_se", "n", "skew", "kurt")]
  da_base$da_ts = da_ts
  da_base$da_p  = da_p
  da_sorted = da_base[order(da_base$da_p), ]
  bonf_da   = bonferroni_sorted(da_sorted$da_p)
  da_out    = cbind(da_sorted, da_bonf = bonf_da$bonf,
                    da_rank = bonf_da$rank, da_alpha = bonf_da$alpha)

  skew_base = ws[, c(grp_cols_present, "trait", "skew", "n")]
  names(skew_base)[names(skew_base) == "skew"] = "skew_obs"
  names(skew_base)[names(skew_base) == "n"]    = "skew_n"
  skew_base$skew_expected = 0
  skew_base$skew_se       = se_skew
  skew_base$skew_ts       = skew_ts
  skew_base$skew_df       = 10000
  skew_base$skew_p        = skew_p
  skew_sorted = skew_base[order(skew_base$skew_p), ]
  bonf_skew   = bonferroni_sorted(skew_sorted$skew_p)
  skew_out    = cbind(skew_sorted, skew_bonf = bonf_skew$bonf,
                      skew_rank = bonf_skew$rank, skew_alpha = bonf_skew$alpha)

  kurt_base = ws[, c(grp_cols_present, "trait", "kurt", "n")]
  names(kurt_base)[names(kurt_base) == "kurt"] = "kurt_obs"
  names(kurt_base)[names(kurt_base) == "n"]    = "kurt_n"
  kurt_base$kurt_expected = 0
  kurt_base$kurt_se       = se_kurt
  kurt_base$kurt_ts       = kurt_ts
  kurt_base$kurt_df       = 10000
  kurt_base$kurt_p        = kurt_p
  kurt_sorted = kurt_base[order(kurt_base$kurt_p), ]
  bonf_kurt   = bonferroni_sorted(kurt_sorted$kurt_p)
  kurt_out    = cbind(kurt_sorted, kurt_bonf = bonf_kurt$bonf,
                      kurt_rank = bonf_kurt$rank, kurt_alpha = bonf_kurt$alpha)

  cat("\n  DA sort (by da_p):\n")
  print_cols = c(grp_cols_present, "trait", "da_ts", "da_p",
                 "da_bonf", "da_rank", "da_alpha")
  print(da_out[, intersect(print_cols, names(da_out))])

  list(da = da_out, skew = skew_out, kurt = kurt_out)
}

# assess_da_vs_fa4a ------------------------------------------------------
assess_da_vs_fa4a = function(f_data, highDA_list, grp_cols = group_cols) {
  # Step 9 DA assessment using Palmer & Strobeck rule of thumb.
  # For each trait in step9_highDA, compares |mean(R-L)| (DA) to FA4a
  # (mean |R-L|). If DA > FA4a, the directional component dominates and
  # the trait is flagged for exclusion.
  #
  # Arguments:
  #   f_data      -- path to post-Step6 CSV
  #   highDA_list -- step9_highDA constant from 2-inspection-constants.R
  #   grp_cols    -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with grouping columns, trait, da, fa4a, exceeds (logical)
  dat  = read_csv(f_data, show_col_types = FALSE)
  rows = list()

  for (entry in highDA_list) {
    for (tr in entry$trait) {
      sub = dat
      for (col in grp_cols) {
        if (col %in% names(entry)) {
          sub = sub[sub[[col]] == entry[[col]], ]
        }
      }
      sub = sub[sub$Trait == tr & !is.na(sub$Value), ]

      avg = sub |>
        group_by(.data[[col_individual]], .data[[col_side]]) |>
        summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

      left_df  = avg |> filter(.data[[col_side]] == col_side_left) |>
        rename(val_l = Value) |> select(all_of(col_individual), val_l)
      right_df = avg |> filter(.data[[col_side]] == col_side_right) |>
        rename(val_r = Value) |> select(all_of(col_individual), val_r)
      merged   = inner_join(left_df, right_df, by = col_individual) |>
        mutate(diff = val_r - val_l)

      if (nrow(merged) == 0) next

      da   = abs(mean(merged$diff))
      fa4a = mean(abs(merged$diff))

      row_key = paste(c(sapply(grp_cols, function(col) entry[[col]]), tr),
                      collapse = "-")
      row = as.data.frame(
        setNames(as.list(sapply(grp_cols, function(col) entry[[col]])),
                 grp_cols),
        stringsAsFactors = FALSE)
      row$trait   = tr
      row$da      = da
      row$fa4a    = fa4a
      row$exceeds = da > fa4a
      rows[[row_key]] = row
    }
  }

  out = do.call(rbind, rows)
  out = out[do.call(order, out[, c(grp_cols, "trait")]), ]
  rownames(out) = NULL
  cat("\nDA vs FA4a Assessment:\n")
  print(out)
  out
}



# ============================================================
# SIZE DEPENDENCY (STEP 8)
# ============================================================

# run_step8a_scatters ----------------------------------------------------
run_step8a_scatters = function(df, trait_cols, grp) {
  # Step 8A: scatter plots of |R-L| vs (R+L)/2 per trait per group.
  # Slow step -- run once, skip on Spearman reruns.
  # Saves interactive HTML scatterplots to figures-inspection/.
  #
  # Arguments:
  #   df         -- long-format dataframe for one group (from load_group)
  #   trait_cols -- character vector of trait names for this group
  #   grp        -- group key string (used in file names and plot titles)
  cat("\n---", grp, "---\n")

  df_avg   = avg_reps(df, trait_cols)
  left_df  = df_avg[df_avg[[col_side]] == col_side_left, ] |> dplyr::select(-all_of(col_side))
  right_df = df_avg[df_avg[[col_side]] == col_side_right, ] |> dplyr::select(-all_of(col_side))

  long_l = pivot_longer(left_df,  cols = all_of(trait_cols),
                        names_to = "var", values_to = "val_left")
  long_r = pivot_longer(right_df, cols = all_of(trait_cols),
                        names_to = "var", values_to = "val_right")

  rl           = merge(long_l, long_r, by = c(col_individual, "var"))
  rl           = na.omit(rl)
  rl$abs_diff  = abs(rl$val_right - rl$val_left)
  rl$trait_avg = (rl$val_right + rl$val_left) / 2

  for (v in trait_cols) {
    sub = rl[rl$var == v, ]
    if (nrow(sub) < 3) next
    p = ggplot(sub, aes(trait_avg, abs_diff, text = .data[[col_individual]])) +
      geom_point() +
      theme_classic() +
      labs(x = "(R+L)/2", y = "|R-L|",
           title = paste("Step 8: Size dependency --", grp, "--", v))
    fname = file.path("figures-inspection",
                      paste0("step8-", grp, "-", v, ".html"))
    saveWidget(ggplotly(p, tooltip = "text"), file = fname, selfcontained = TRUE)
  }
  invisible(NULL)
}

# run_step8b_spearman ----------------------------------------------------
run_step8b_spearman = function(df, trait_cols, grp, grp_cols = group_cols) {
  # Step 8B: Spearman correlation of |R-L| vs (R+L)/2 per trait per group.
  # Fast; re-runnable independently of scatter plots.
  # Sequential Bonferroni correction applied within each group.
  # Significant positive association indicates trait size dependency;
  # negative associations are not biologically meaningful for this test.
  #
  # Arguments:
  #   df         -- long-format dataframe for one group (from load_group)
  #   trait_cols -- character vector of trait names for this group
  #   grp        -- group key string
  #   grp_cols   -- character vector of grouping column names (default: group_cols)
  #
  # Returns: dataframe with grouping columns, Spearman z, p, Bonferroni correction
  cat("\n---", grp, "---\n")
  grp_vals = parse_group_key(grp, grp_cols)

  df_avg   = avg_reps(df, trait_cols)
  left_df  = df_avg[df_avg[[col_side]] == col_side_left, ] |> dplyr::select(-all_of(col_side))
  right_df = df_avg[df_avg[[col_side]] == col_side_right, ] |> dplyr::select(-all_of(col_side))

  long_l = pivot_longer(left_df,  cols = all_of(trait_cols),
                        names_to = "var", values_to = "val_left")
  long_r = pivot_longer(right_df, cols = all_of(trait_cols),
                        names_to = "var", values_to = "val_right")

  rl           = merge(long_l, long_r, by = c(col_individual, "var"))
  rl           = na.omit(rl)
  rl$abs_diff  = abs(rl$val_right - rl$val_left)
  rl$trait_avg = (rl$val_right + rl$val_left) / 2

  rows = list()
  for (v in trait_cols) {
    sub = rl[rl$var == v, ]
    if (nrow(sub) < 3) next
    tst   = spearman_test(abs_diff ~ trait_avg, data = sub)
    zstat = as.numeric(statistic(tst))
    pval  = pvalue(tst)
    row   = as.data.frame(grp_vals, stringsAsFactors = FALSE)
    row$trait = v
    row$n     = nrow(sub)
    row$z     = zstat
    row$p     = pval
    rows[[v]] = row
  }

  out       = do.call(rbind, rows)
  out       = out[order(out$p), ]
  bonf      = bonferroni_sorted(out$p)
  out$bonf  = bonf$bonf
  out$rank  = bonf$rank
  out$alpha = bonf$alpha
  cat("\n")
  print(out)
  out
}



# ============================================================
# NORMALITY SCREENING
# ============================================================

# run_normality_scatters -------------------------------------------------
run_normality_scatters = function(f_data, prefix, grp_vars = group_trait_vars) {
  # Scatterplot matrix (QQ diagonal) per analytical group.
  # Group label is derived from the key string directly.
  # Saves one TIFF per group to figures-inspection/.
  #
  # Arguments:
  #   f_data   -- path to the pipeline CSV
  #   prefix   -- filename prefix for output TIFFs (e.g. "normal-scatters-raw")
  #   grp_vars -- named list of trait vectors per group (default: group_trait_vars)
  for (grp in names(grp_vars)) {
    g    = load_group(f_data, grp, grp_vars[[grp]])
    df   = g$df
    vars = g$trait_cols

    if (length(vars) == 0) {
      cat("Skipping", grp, "-- no traits retained\n")
      next
    }

    df_wide = df |>
      filter(Trait %in% vars) |>
      select(all_of(c(col_individual, col_replicate, col_side)), Trait, Value) |>
      pivot_wider(names_from = Trait, values_from = Value) |>
      select(all_of(vars))

    if (all(is.na(df_wide))) {
      cat("Skipping", grp, "-- all values NA\n")
      next
    }

    vars = vars[sapply(vars, function(v) any(!is.na(df_wide[[v]])))]
    if (length(vars) == 0) {
      cat("Skipping", grp, "-- no non-NA traits\n")
      next
    }
    df_wide = df_wide[, vars, drop = FALSE]

    fname = file.path("figures-inspection",
                      paste0(prefix, "-", grp, ".tiff"))
    tiff(fname, units = "in", width = 20, height = 20,
         compression = "lzw", res = 300)
    # Attempt scatterplot matrix (QQ diagonal, pairwise scatters).
    # Falls back to individual QQ plots per trait if scatterplotMatrix fails
    # -- this occurs when pairwise complete cases are insufficient for some
    # trait pairs, which is common in datasets with many missing traits
    # (e.g. dental studies where not all teeth are present in all individuals).
    suppressWarnings(tryCatch({
      scatterplotMatrix(
        df_wide[, vars, drop = FALSE],
        cex        = 0.5,
        diagonal   = list(method = "qqplot"),
        var.labels = vars,
        main       = paste("Normality --", grp)
      )
    }, error = function(e) {
      cat("  scatterplotMatrix failed for", grp, "-- falling back to individual QQ plots\n")
      ncols = min(4, length(vars))
      nrows = ceiling(length(vars) / ncols)
      par(mfrow = c(nrows, ncols))
      for (v in vars) {
        vals = df_wide[[v]]
        vals = vals[!is.na(vals)]
        if (length(vals) >= 3) {
          suppressWarnings({
            qqnorm(vals, main = paste(grp, "--", v))
            qqline(vals, col = "red")
          })
        }
      }
      par(mfrow = c(1, 1))
    }))
    dev.off()
    cat("Written:", fname, "\n")
  }
}