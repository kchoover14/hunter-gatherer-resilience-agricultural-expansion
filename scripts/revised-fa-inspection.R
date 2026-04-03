######################## LIBRARIES
library(readr)       # read_csv(), write_csv()
library(dplyr)       # data manipulation
library(tidyr)       # pivot_longer(), pivot_wider()
library(ggplot2)     # plotting
library(car)         # leveneTest(), scatterplotMatrix()
library(outliers)    # dixon.test()
library(DescTools)   # Skew(), Kurt()
library(Hmisc)       # describe()
library(coin)        # spearman_test()
library(plotly)      # ggplotly()
library(htmlwidgets) # saveWidget()



######################## SOURCE SHARED FILES
source("scripts/revised-fa-inspection-functions.R")
source("scripts/revised-fa-inspection-constants.R")



######################## NORMALITY SCREENING (PRE-SCREENING)
# Scatterplot matrix per analytical group. One TIFF per group saved to
# figures-inspection/. Diagonal panels show QQ plots; off-diagonal panels
# show pairwise scatters across traits.
#
# Purpose: visual assessment of distributional properties before any
# screening begins. Look for gross departures from normality -- heavy
# multi-modality or extreme skew that exceeds typical biological variation.
# Mild tails and moderate S-curves are expected in biological FA data and
# do not preclude parametric testing. This is a visual check only; no
# formal decision is required here.

dir.create("figures-inspection", showWarnings = FALSE)
run_normality_scatters(f_wrangled, "normal-scatters-raw")



######################## STEP 1: REPLICATE-PAIR PLOTS (VISUAL INSPECTION)
# For each individual x trait x group, plots |R-L| values across replicates
# so you can visually identify any single replicate whose side difference
# departs markedly from the others for that individual.
#
# One interactive HTML plot per group, faceted by trait. Individual IDs
# appear as tooltips. Saved to figures-inspection/.
#
# NOTE: The number of plots scales with individuals x traits x groups.
# For large datasets this step may be impractical. In that case, proceed
# directly to Step 2 (Dixon's test), which provides the statistical
# equivalent. For small or moderate datasets, run Step 1 first as it gives
# valuable visual context for interpreting Step 2 flags.

for (grp in names(group_trait_vars)) {
  g = load_group(f_wrangled, grp, group_trait_vars[[grp]])

  # Build one |R-L| value per individual per replicate per trait.
  step1_rows = list()
  for (ind in levels(g$df[[col_individual]])) {
    ind_df = g$df[g$df[[col_individual]] == ind, ]
    for (tr in g$trait_cols) {
      tr_df = ind_df[ind_df$Trait == tr, ]
      tr_df[[col_replicate]] = as.integer(as.character(tr_df[[col_replicate]]))
      for (r in sort(unique(tr_df[[col_replicate]]))) {
        rv = tr_df$Value[tr_df[[col_replicate]] == r & tr_df[[col_side]] == col_side_right]
        lv = tr_df$Value[tr_df[[col_replicate]] == r & tr_df[[col_side]] == col_side_left]
        if (length(rv) == 1 && length(lv) == 1 &&
            !is.na(rv) && !is.na(lv)) {
          step1_rows[[length(step1_rows) + 1]] = data.frame(
            Ind = ind, Trait = tr, Rep = r,
            abs_rl = abs(rv - lv),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(step1_rows) == 0) next
  step1_df = do.call(rbind, step1_rows)

  p = ggplot(step1_df, aes(x = Rep, y = abs_rl, group = Ind,
                           colour = Ind, text = Ind)) +  # Ind/Rep are internal plot column names
    geom_line(alpha = 0.5) +
    geom_point(size = 1.5) +
    facet_wrap(~Trait, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Replicate", y = "|R-L|",
         title = paste("Step 1: Replicate-pair |R-L| --", grp))

  fname = file.path("figures-inspection",
                    paste0("step1-replpair-", grp, ".html"))
  saveWidget(ggplotly(p, tooltip = "text"), file = fname,
             selfcontained = TRUE)
  cat("Written:", fname, "\n")
}


######################## STEP 2: ME OUTLIER DETECTION (DIXON'S TEST)
# Measurement error (ME) outliers are individuals whose replicate-to-replicate
# consistency in side differences is extreme relative to the group.
#
# For each individual x trait x group, the procedure is:
#   1. Compute |R-L| for each replicate -- one absolute side difference per rep.
#   2. Compute consecutive differences between those |R-L| values across
#      replicates (e.g. |rep2_sidediff - rep1_sidediff|, etc.).
#      This measures how much the side difference varies between replicate pairs.
#   3. Run Dixon's two-tailed test on those consecutive replicate differences
#      against the group mean of all such differences AND against zero.
#   Bonferroni correction applied within each group separately.
#
# Output: f_checkpoint_step12
#   One row per individual x trait per group.
#   sig_mean == "*" : flagged as outlier relative to group mean replicate diff
#   sig_zero == "*" : flagged as outlier relative to zero

me_dixon_rows = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_wrangled, grp, group_trait_vars[[grp]])
  me_dixon_rows[[grp]] = run_me_dixon(g$df, g$trait_cols, grp)
}

me_dixon_all = do.call(rbind, me_dixon_rows)
write_csv(me_dixon_all, f_checkpoint_step12)
cat("Written:", f_checkpoint_step12, "\n")


######################## CHECKPOINT: STEPS 1-2
# Review f_checkpoint_step12.
# sig_mean or sig_zero == "*" indicates a flagged individual x trait.
#
# For each flagged case:
#   1. Run diagnose_me_flagged() below to see the consecutive replicate
#      differences and identify which specific replicate pair is driving
#      the flag. abs_rl_a is the lower-numbered replicate in the pair;
#      abs_rl_b is the higher-numbered replicate.
#   2. Compare the extreme replicate's |R-L| value to the others for
#      that individual. If one replicate is clearly aberrant, it is the
#      offending measurement.
#   3. Add confirmed removals to me_dixon_removals in scripts/revised-fa-inspection-constants.R.
#      Removal unit: one specific replicate for one individual x trait.
#
# After populating me_dixon_removals, save scripts/revised-fa-inspection-constants.R and continue.

raw_groups = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_wrangled, grp, group_trait_vars[[grp]])
  raw_groups[[grp]] = g$df
}
me_diagnostic = diagnose_me_flagged(me_dixon_all, raw_groups)



######################## STEPS 3-5
# Steps 3-5 detect FA outliers -- individuals whose side difference itself
# is extreme relative to the group, as opposed to ME outliers in Steps 1-2
# whose replicate-to-replicate variation in side difference is extreme.
#
# All three steps read from f_wrangled and can be re-run independently.
#
# Step 3: L vs R scatterplots -- reveals extreme asymmetries and outlier-
#         sized individuals within traits. Visual inspection only.
# Step 4: (R-L) vs (R-L) cross-trait scatterplots -- reveals subtler FA
#         outliers visible as consistent off-diagonal patterns across traits.
#         An individual appearing extreme in many trait pairs is a stronger
#         candidate for removal than one appearing in a single pair.
#         Visual inspection only.
# Step 5: Dixon's test on signed R-L side differences -- statistical
#         confirmation of univariate FA outliers per trait.
#         Cross-trait patterns identified visually in Step 4 may not always
#         be confirmed here; visual judgment takes precedence for individuals
#         with consistent cross-trait patterns.



######################## STEP 3: L vs R SCATTERPLOTS (ALL GROUPS)
# One interactive HTML plot per group, faceted by trait.
# Individual IDs appear as tooltips on hover.
# Saved to figures-inspection/.
#
# Look for: points far from the identity line (extreme asymmetry), or
# individuals that appear consistently extreme across multiple traits.

for (grp in names(group_trait_vars)) {
  g   = load_group(f_wrangled, grp, group_trait_vars[[grp]])
  rl  = build_rl(g$df, g$trait_cols)

  rl_lims = range(c(rl$val_left, rl$val_right), na.rm = TRUE)
  rl_pad  = diff(rl_lims) * 0.05
  rl_lo   = rl_lims[1] - rl_pad
  rl_hi   = rl_lims[2] + rl_pad

  p = ggplot(rl, aes(val_left, val_right, text = .data[[col_individual]])) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey50") +
    coord_equal(xlim = c(rl_lo, rl_hi), ylim = c(rl_lo, rl_hi)) +
    facet_wrap(~var) + theme_classic() +
    labs(x = "Left", y = "Right",
         title = paste("Step 3: L vs R --", grp))

  fname = file.path("figures-inspection",
                    paste0("step3-lvr-", grp, ".html"))
  saveWidget(ggplotly(p, tooltip = "text"), file = fname,
             selfcontained = TRUE)
  cat("Written:", fname, "\n")
}



######################## STEP 4: CROSS-TRAIT (R-L) SCATTERPLOTS (ALL GROUPS)
# One interactive HTML plot per trait pair per group.
# Individual IDs appear as tooltips on hover.
# Saved to figures-inspection/.
#
# Look for: individuals appearing as consistent outliers across multiple
# trait pairs. Consistency across pairs is stronger evidence of a genuine
# FA anomaly than a single-trait flag.

for (grp in names(group_trait_vars)) {
  g        = load_group(f_wrangled, grp, group_trait_vars[[grp]])
  sd_wide  = build_sd_wide(g$df, g$trait_cols)
  sd_cols  = grep("_sd$", names(sd_wide), value = TRUE)

  if (length(sd_cols) < 2) next

  trait_pairs = combn(sd_cols, 2)
  for (i in seq_len(ncol(trait_pairs))) {
    xv = trait_pairs[1, i]; yv = trait_pairs[2, i]
    p  = ggplot(sd_wide, aes(.data[[xv]], .data[[yv]], text = .data[[col_individual]])) +
      geom_point() +
      theme_classic() +
      labs(x = xv, y = yv,
           title = paste("Step 4: Cross-trait (R-L) --", grp))
    fname = file.path("figures-inspection",
                      paste0("step4-sd-", grp, "-",
                             gsub("_sd$", "", xv), "_vs_",
                             gsub("_sd$", "", yv), ".html"))
    saveWidget(ggplotly(p, tooltip = "text"), file = fname,
               selfcontained = TRUE)
  }
  cat("Step 4 plots saved:", grp, "\n")
}



######################## STEP 5: DIXON'S TEST ON SIGNED (R-L) SIDE DIFFERENCES
# For each individual x trait x group:
#   Computes signed R-L for each replicate -- one side difference per rep.
#   Runs Dixon's two-tailed test on those signed R-L values against the
#   group mean R-L AND against zero.
#   Bonferroni correction applied within each group separately.
#
# This step is fully independent of Steps 3-4 and reads f_wrangled directly.
#
# Output: f_checkpoint_step345
#   One row per individual x trait per group.
#   sig_mean == "*" : flagged as outlier relative to group mean R-L
#   sig_zero == "*" : flagged as outlier relative to zero

fa_dixon_rows = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_wrangled, grp, group_trait_vars[[grp]])
  fa_dixon_rows[[grp]] = run_fa_dixon(g$df, g$trait_cols, grp)
}

fa_dixon_all = do.call(rbind, fa_dixon_rows)
write_csv(fa_dixon_all, f_checkpoint_step345)
cat("Written:", f_checkpoint_step345, "\n")


######################## CHECKPOINT: STEPS 3-5
# Review f_checkpoint_step345 together with the Step 3 and Step 4 plots.
# sig_mean or sig_zero == "*" indicates a flagged individual x trait.
#
# For each flagged case:
#   1. Run diagnose_fa_flagged() below to see the signed R-L values per
#      replicate and identify which replicate produced the extreme value.
#   2. Cross-reference with the Step 4 cross-trait plots for that individual.
#      Consistent off-diagonal position across multiple trait pairs strengthens
#      the case for removal.
#   3. Add confirmed removals to fa_dixon_removals in scripts/revised-fa-inspection-constants.R.
#      Removal unit: one specific replicate for one individual x trait.
#
# After populating fa_dixon_removals, save scripts/revised-fa-inspection-constants.R and continue.

raw_groups = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_wrangled, grp, group_trait_vars[[grp]])
  raw_groups[[grp]] = g$df
}
fa_diagnostic = diagnose_fa_flagged(fa_dixon_all, raw_groups)



######################## OUTLIER REMOVAL PRIOR TO STEP 6
# Applies all ME (Steps 1-2) and FA (Steps 3-5) replicate removals in one
# pass. Removal is at the level of the specific offending replicate per
# individual per trait -- not the whole individual or the whole trait.
# Re-source scripts/revised-fa-inspection-constants.R to pick up me_dixon_removals and fa_dixon_removals.

source("scripts/revised-fa-inspection-constants.R")

all_removals = c(me_dixon_removals, fa_dixon_removals)
cat("Total replicate removals loaded:", length(all_removals), "\n")

posts5 = read_csv(f_wrangled, show_col_types = FALSE)

for (entry in all_removals) {
  posts5 = null_replicate(posts5, entry$ind, entry$trait, entry$rep)
}

cat("Replicate removals applied:", length(all_removals), "\n")
write_csv(posts5, f_fromStep5)
cat("Written:", f_fromStep5, "\n")



######################## STEP 6: SIDES x INDIVIDUALS MIXED-MODEL ANOVA (FA > ME)
# Tests whether FA variance (Sides x Individuals interaction mean square)
# significantly exceeds ME variance (within-cell residual mean square) for
# each trait per group. This is the central test of the P&S protocol:
#   F = MS_SxI / MS_ME
# Only traits where FA > ME is significant after Bonferroni correction pass
# this checkpoint and are retained for downstream analysis.
#
# Output columns of note:
#   fa_me_f / fa_me_p : F statistic and p-value for FA > ME test
#   anova_alpha       : "*" = significant, "n.s." = non-significant
#   pct_me            : ME as a percentage of the FA variance component.
#                       High pct_me indicates measurement error is a large
#                       fraction of apparent FA. Useful context even for
#                       traits that pass the FA > ME test.
#   fa10_noME         : ME-corrected FA index (FA10, Palmer & Strobeck 1994
#                       Table 1). This is the FA measure used in hypothesis
#                       testing. FA10a when no ln transformation is applied;
#                       FA10b when ln_transform = TRUE (set at Step 8).
#
# Requires: f_fromStep5

aov_inputs = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_fromStep5, grp, group_trait_vars[[grp]])
  aov_inputs[[grp]] = run_univariates(g$df, g$trait_cols, grp)
}

step6_out = list()
for (grp in names(aov_inputs)) {
  step6_out[[grp]] = run_mixedmodel(grp, aov_inputs)
}

step6_all = do.call(rbind, step6_out)
write_csv(step6_all, f_step6)
cat("\nWritten:", f_step6, "\n")


######################## CHECKPOINT: STEP 6
# Review f_step6.
# anova_alpha == "n.s." indicates FA does not significantly exceed ME for
# that trait in that group after Bonferroni correction.
#
# For each non-significant trait:
#   Add it to step6_eliminations in scripts/revised-fa-inspection-constants.R.
#   Each entry must include values for all group_cols plus the trait name.
#   To remove a trait from every group simultaneously, set all group_cols
#   values to "all".
#   Example with group_cols = c("Sex", "Dimension"):
#     list(Sex = "female", Dimension = "breadth", trait = "trait1")
#     list(Sex = "all",    Dimension = "all",     trait = "trait2")
#
# Leave step6_eliminations as an empty list if all traits are significant.
# Save scripts/revised-fa-inspection-constants.R before continuing.



######################## STEP 6 REMOVAL
# Removes traits with non-significant FA > ME from the dataset.
# Removal unit: entire trait column for a group (Value set to NA).
# Re-source scripts/revised-fa-inspection-constants.R to pick up step6_eliminations.

source("scripts/revised-fa-inspection-constants.R")
cat("Step 6 eliminations loaded:", length(step6_eliminations), "entries\n")

posts6 = read_csv(f_fromStep5, show_col_types = FALSE)

for (entry in step6_eliminations) {
  mask = rep(TRUE, nrow(posts6))
  for (col in group_cols) {
    if (!identical(entry[[col]], "all")) {
      mask = mask & (posts6[[col]] == entry[[col]])
    }
  }
  posts6$Value[mask & posts6$Trait %in% entry$trait] = NA
}

write_csv(posts6, f_fromStep6)
cat("Written:", f_fromStep6, "\n")



######################## STEP 7: LEVENE'S TEST ON ME VARIANCE
# Tests whether ME variance is homogeneous across all grouping variables
# and traits (Palmer & Strobeck 2003, Step 7).
#
# ME is measured as |consecutive replicate differences| per individual per
# trait per side -- the same replicate difference logic used in Steps 1-2,
# but applied to raw measurement values rather than side differences.
#
# If ME variance is heterogeneous, this should be reported in any publication.
# Heterogeneous ME does not invalidate the analysis but affects interpretation
# of FA comparisons across groups.
#
# The Levene's formula is built dynamically from group_cols and Trait.
# Requires: f_fromStep5

me_var = read_csv(f_fromStep5, show_col_types = FALSE) |>
  filter(!is.na(Value)) |>
  arrange(.data[[col_individual]], across(all_of(group_cols)), Trait, .data[[col_side]], .data[[col_replicate]]) |>
  group_by(.data[[col_individual]], across(all_of(group_cols)), Trait, .data[[col_side]]) |>
  summarise(
    me = list(abs(diff(Value))),
    .groups = "drop"
  ) |>
  tidyr::unnest(me) |>
  filter(!is.na(me))

lev_formula = as.formula(
  paste("me ~", paste(c(group_cols, "Trait"), collapse = " * "))
)
leveneTest(lev_formula, data = me_var, center = mean)



######################## STEP 8: TRAIT SIZE DEPENDENCY
# Tests whether FA (|R-L|) is correlated with trait size [(R+L)/2] per
# trait per group. A significant positive association indicates that larger
# individuals show more asymmetry, which may reflect scaling rather than
# developmental instability.
#
# FA10a vs FA10b:
#   If no significant positive association is found, raw values are used
#   in hypothesis testing and the FA index is FA10a.
#   If any trait in any group shows a significant positive association,
#   ln transformation is applied to ALL values before hypothesis testing
#   (FA10b). Mixing transformed and untransformed values across traits
#   would make FA indices incomparable, so the transformation must be
#   applied globally if needed at all.
#   Record your decision as ln_transform = TRUE or FALSE in scripts/revised-fa-inspection-constants.R.
#
# Elimination only justified for significant AND positive associations.
# Negative associations are not biologically meaningful for this test
# (Palmer & Strobeck 2003).
# Requires: f_fromStep6


#### STEP 8A: SCATTER PLOTS (run once -- slow)
# Generates interactive |R-L| vs (R+L)/2 scatterplots per trait per group.
# Saved to figures-inspection/. Run once; skip on reruns.

for (grp in names(group_trait_vars)) {
  g = load_group(f_fromStep6, grp, group_trait_vars[[grp]])
  run_step8a_scatters(g$df, g$trait_cols, grp)
}


#### STEP 8B: SPEARMAN CORRELATION (re-runnable independently)
# Rank-based correlation, robust to non-normality.
# Sequential Bonferroni correction applied within each group.
#
# Output: f_checkpoint_step8
#   One row per trait per group.
#   alpha == "*" and z > 0 : significant positive association -- action required
#   alpha == "*" and z < 0 : significant negative association -- no action
#   alpha == "n.s."        : no significant association

step8_out = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_fromStep6, grp, group_trait_vars[[grp]])
  step8_out[[grp]] = run_step8b_spearman(g$df, g$trait_cols, grp)
}

step8_combined = do.call(rbind, step8_out)
write_csv(step8_combined, f_checkpoint_step8)
cat("\nWritten:", f_checkpoint_step8, "\n")


######################## CHECKPOINT: STEP 8
# Review f_checkpoint_step8 and the Step 8A scatterplots.
#
# If any trait shows a significant positive association (alpha == "*" and
# z > 0) after Bonferroni correction:
#   Set ln_transform = TRUE in scripts/revised-fa-inspection-constants.R.
#   The ln transformation will be applied to all values in all groups at
#   the hypothesis testing stage -- no trait removal is needed.
#
# If no significant positive associations are found:
#   Leave ln_transform = FALSE in scripts/revised-fa-inspection-constants.R.
#
# Save scripts/revised-fa-inspection-constants.R before continuing.



######################## STEP 9: IDEAL FA ASSESSMENT (DA, SKEW, KURTOSIS)
# Tests each trait per group for departures from ideal FA:
#   Directional asymmetry (DA): mean(R-L) significantly different from zero.
#   Skewness: distribution of (R-L) values departs from symmetry.
#   Kurtosis: distribution of (R-L) values departs from normality.
# Each test is independently Bonferroni-sorted within each group.
# Requires: f_fromStep6

desc_inputs = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_fromStep6, grp, group_trait_vars[[grp]])
  desc_inputs[[grp]] = run_descriptives(g$df, g$trait_cols, grp)
}

step9_out = list()
for (grp in names(group_trait_vars)) {
  step9_out[[grp]] = run_idealfa(grp, desc_inputs)
}

step9_da_all   = do.call(rbind, lapply(step9_out, function(x) x$da))
step9_skew_all = do.call(rbind, lapply(step9_out, function(x) x$skew))
step9_kurt_all = do.call(rbind, lapply(step9_out, function(x) x$kurt))

write_csv(step9_da_all,   sub("\\.csv", "-da.csv",   f_checkpoint_step9))
write_csv(step9_skew_all, sub("\\.csv", "-skew.csv", f_checkpoint_step9))
write_csv(step9_kurt_all, sub("\\.csv", "-kurt.csv", f_checkpoint_step9))
cat("\nWritten: Step 9 checkpoint files\n")


######################## CHECKPOINT: STEP 9 (PART 1 -- DA)
# Review f_checkpoint_step9-da.csv.
# da_alpha == "*" indicates significant directional asymmetry.
#
# Significant DA does not automatically disqualify a trait. See the
# step9_highDA section in scripts/revised-fa-inspection-constants.R for the Palmer & Strobeck (2003)
# rule of thumb on when DA is small enough to retain a trait.
#
# Populate step9_highDA in scripts/revised-fa-inspection-constants.R with all traits showing
# significant DA. Save scripts/revised-fa-inspection-constants.R before running the DA assessment below.

source("scripts/revised-fa-inspection-constants.R")


######################## DA ASSESSMENT USE ONLY IF DA SIG
# For each trait in step9_highDA, compares |mean(R-L)| (DA) to FA4a
# (mean |R-L|). If DA does not exceed FA4a, the trait may be retained.
# If DA exceeds FA4a, the trait must be excluded. Borderline cases should
# be excluded conservatively.
#
# Output: f_checkpoint_step9_daTest
#   exceeds == TRUE  : DA dominates, trait must be excluded
#   exceeds == FALSE : DA is within acceptable range, trait may be retained

#da_assessment = assess_da_vs_fa4a(f_fromStep6, step9_highDA)
#write_csv(da_assessment,
     #     sub("\\.csv", "-da-assessment.csv", f_checkpoint_step9_daTest))
# cat("Written: DA assessment\n")


######################## CHECKPOINT: STEP 9 (PART 2 -- SKEW, KURTOSIS, DA)
# Review f_checkpoint_step9-skew.csv and f_checkpoint_step9-kurt.csv.
# skew_alpha or kurt_alpha == "*" indicates a flagged trait.
# There is no retention rule for significant skew or kurtosis -- these
# traits must be excluded.
#
# Also review the DA assessment output. Traits where exceeds == TRUE
# must be excluded.
#
# Populate step9_eliminations in scripts/revised-fa-inspection-constants.R combining DA exceedances,
# skew flags, and kurtosis flags. Save scripts/revised-fa-inspection-constants.R before continuing.



######################## STEP 9 REMOVALS
# Removes all traits flagged for sig DA, skew, or kurtosis.
# Removal unit: entire trait column for a group (Value set to NA).

source("scripts/revised-fa-inspection-constants.R")

posts9 = read_csv(f_fromStep6, show_col_types = FALSE)

for (entry in step9_eliminations) {
  for (tr in entry$trait) {
    mask = rep(TRUE, nrow(posts9))
    for (col in group_cols) {
      if (!identical(entry[[col]], "all")) {
        mask = mask & (posts9[[col]] == entry[[col]])
      }
    }
    posts9$Value[mask & posts9$Trait == tr] = NA
  }
}

write_csv(posts9, f_fromStep9)
cat("Written:", f_fromStep9, "\n")



######################## NORMALITY SCREENING (POST-SCREENING)
# Repeat the normality scatterplot matrix on the final cleaned dataset.
# Compare to the pre-screening plots. Mild tails and moderate S-curves
# are typical for biological FA data and do not preclude parametric testing.

run_normality_scatters(f_fromStep9, "normal-scatters-final")



######################## GENERATE HYPOTHESIS TESTING DATA
# Computes the ME-corrected FA index (FA10) from the cleaned dataset for
# use in 3-analysis.R.
#
# FA10a vs FA10b (controlled by ln_transform in 2-inspection-constants.R):
#   ln_transform = FALSE : raw values used, index is FA10a
#   ln_transform = TRUE  : ln transformation applied to all values in all
#                          groups before ANOVA, index is FA10b.
#                          Set when any trait showed a significant positive
#                          size dependency in Step 8.
#
# The ln_transform flag is set at the Step 8 checkpoint in 2-inspection-constants.R.

source("scripts/revised-fa-inspection-constants.R")

aov_inputs_final = list()
for (grp in names(group_trait_vars)) {
  g = load_group(f_fromStep9, grp, group_trait_vars[[grp]])
  if (ln_transform) {
    g$df$Value = log(g$df$Value)
    cat(grp, "-- ln transformation applied (FA10b)\n")
  } else {
    cat(grp, "-- no transformation (FA10a)\n")
  }
  aov_inputs_final[[grp]] = run_univariates(g$df, g$trait_cols, grp)
}

mm_inputs_final = list()
for (grp in names(aov_inputs_final)) {
  mm_inputs_final[[grp]] = run_mixedmodel(grp, aov_inputs_final)
}

hypothesis_all = do.call(rbind, mm_inputs_final)
write_csv(hypothesis_all, f_hypothesis)
cat("Written:", f_hypothesis, "\n")



######################## TIDY
rm(list = ls())
gc()
