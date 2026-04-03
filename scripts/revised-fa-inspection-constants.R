######################## FILE PATHS
# Update these paths to match your project directory structure.
# All paths should be relative to your working directory.

f_wrangled               = "data/revised-fa-revised-wrangled.csv"
f_checkpoint_step12      = "data/revised-fa-checkpoint-step1-2-dixon.csv"
f_checkpoint_step345     = "data/revised-fa-checkpoint-step3-4-5-dixon.csv"
f_fromStep5              = "data/revised-fa-fromStep5.csv"
f_step6                  = "data/revised-fa-checkpoint-step6-mixedmodel.csv"
f_fromStep6              = "data/revised-fa-fromStep6.csv"
f_checkpoint_step7       = "data/revised-fa-checkpoint-step7-leveneME.csv"
f_checkpoint_step8       = "data/revised-fa-checkpoint-step8-spearman.csv"
f_checkpoint_step9       = "data/revised-fa-checkpoint-step9-daSkewKurtosis.csv"
f_checkpoint_step9_daTest = "data/revised-fa-checkpoint-step9-daTest.csv"
f_fromStep9              = "data/revised-fa-fromStep9.csv"
f_hypothesis             = "data/revised-fa-hypothesis.csv"



######################## GROUPING VARIABLES
group_cols = c("period")  # replace with your actual column names





######################## COLUMN NAMES
col_individual = "uid"
col_side       = "side"
col_replicate  = "rep"
col_side_right = "R"
col_side_left  = "L"





######################## TRAIT VARIABLES
trait_vars_all = c('xlm3', 'xlm2', 'xlm1', 'xlp2', 'xlp1', 'xlc', 'xli2', 'xli1', 'xbm3', 'xbm2', 'xbm1', 'xbp2', 'xbp1', 'xbc', 'xbi2', 'xbi1', 'mlm3', 'mlm2', 'mlm1', 'mlp2', 'mlp1', 'mlc', 'mli2', 'mli1', 'mbm3', 'mbm2', 'mbm1', 'mbp2', 'mbp1', 'mbc', 'mbi2', 'mbi1'
)



######################## TRAIT VARIABLE SETS PER GROUP
group_trait_vars = list(
  "Jomon" = trait_vars_all,
  "Yayoi" = trait_vars_all
)



######################## OUTLIER REMOVALS: STEPS 1-2 (ME OUTLIERS)
# Populated after reviewing f_checkpoint_step12 and the ME diagnostic output.
#
# Removal unit: one specific replicate for one individual x trait combination.
# Use diagnose_me_flagged() output to identify the offending replicate.
# sig_mean or sig_zero == "*" in the checkpoint file indicates a flagged case.
#
# Format: one list entry per removal. ind is the uid value.
#   Example with group_cols = c("period"):
#     list(period = "Jomon", ind = "Wakimisaki_1", trait = "xlm3", rep = 2)
#
# Leave as an empty list until you have reviewed the checkpoint output.

me_dixon_removals = list(
  # list(period = "...", ind = "...", trait = "...", rep = N),
)



######################## OUTLIER REMOVALS: STEPS 3-5 (FA OUTLIERS)
# Populated after reviewing f_checkpoint_step345 and the FA diagnostic output.
# Same format as me_dixon_removals above.
#
# sig_mean or sig_zero == "*" indicates a flagged individual x trait.
# Cross-trait patterns identified visually in Step 4 may not always be
# confirmed statistically in Step 5; visual judgment takes precedence.
#
# Leave as an empty list until you have reviewed the checkpoint output.

fa_dixon_removals = list(
  # list(period = "...", ind = "...", trait = "...", rep = N),
)



######################## STEP 6 ELIMINATIONS (FA > ME NON-SIGNIFICANT)
# Populated after reviewing f_step6.
#
# Removal unit: entire trait column for one group (or all groups).
# Remove any trait where FA > ME is non-significant after Bonferroni
# correction. If a trait fails in every group, use "all" for that
# grouping variable value.
#
# Format: one list entry per group, trait as a vector.
#   Example with group_cols = c("period"):
#     list(period = "Jomon", trait = c("trait1", "trait2", "trait3"))
#     list(period = "all",   trait = c("trait4"))
#
# Leave as an empty list if all traits have significant FA > ME.

step6_eliminations = list(
  list(period = "Jomon", trait = c("mlp2", "xbm1", "xlm2", "xli2", "xbi2", "mlc", "xlp2", "mbp2", "mbm1", "mbm2", "mbp1", "mlm1", "xlp1", "xlm1", "xbp2", "xbc", "mli2", "mbc", "mli1", "mlp1", "xli1", "xbm2", "xbp1", "mbi2", "mbi1", "xlc")),
  list(period = "Yayoi", trait = c("mlm1", "mbp1", "mbp2", "mli2", "xlc", "xlp1", "mbi2", "xlm2", "xbm2", "mbi1", "mlp1", "mbm1", "mlm2", "mlp2", "mbc", "xlp2", "xli1", "mli1", "mlc", "xbp2", "mbm2"))
)



######################## STEP 8: LN TRANSFORMATION (FA10a vs FA10b)
# Populated after reviewing f_checkpoint_step8.
#
# If any trait in any group shows a significant POSITIVE Spearman correlation
# between |R-L| and (R+L)/2 after sequential Bonferroni correction, a size
# dependency exists. In this case, ln transformation is applied to ALL values
# before hypothesis testing.
#
# Why all values, not just the affected trait or group?
#   Mixing raw (mm) and ln-transformed values across traits would make FA
#   indices incomparable. If any trait requires transformation, all traits
#   in all groups are transformed for consistency. The resulting index is
#   FA10b (ln-corrected) rather than FA10a (raw).
#
# Negative associations are not biologically meaningful for this test and
# do not justify transformation (Palmer & Strobeck 2003).
#
# Set to TRUE if any significant positive association was found; FALSE otherwise.
# Default is FALSE (no transformation, FA10a).

ln_transform = FALSE #one value with sig ts dependency removed skew/kurt



######################## STEP 9: HIGH DA TRAITS (FOR DA vs FA4a ASSESSMENT)
# Populated after reviewing the DA checkpoint file (f_checkpoint_step9 -da.csv).
#
# List traits with statistically significant directional asymmetry (DA) for
# the DA vs FA4a assessment. Significant DA does not automatically disqualify
# a trait. Palmer & Strobeck (2003) provide the following guidance:
#   If |mean(R-L)| (DA) is no larger than FA4a (mean |R-L|), the predisposition
#   towards one side is smaller than the average random deviation, and the trait
#   may be retained. If DA exceeds FA4a, the directional component dominates
#   and the trait should be excluded. Borderline cases should be excluded
#   conservatively. Note that in large samples even very slight DA may become
#   statistically significant; the DA vs FA4a comparison is the relevant
#   criterion, not statistical significance alone.
#
# There is no equivalent retention rule for significant skew or kurtosis --
# traits flagged for those must be excluded (see step9_eliminations below).
#
# Format: one list entry per group with high-DA traits, trait as a vector.
#   Example with group_cols = c("period"):
#     list(period = "Jomon", trait = c("trait1", "trait2"))

step9_highDA = list(
  # list(period = "...", trait = c("trait1", "trait2")),
)


#none sig for da
######################## STEP 9 ELIMINATIONS (DA, SKEW, KURTOSIS)
# Populated after reviewing all three Step 9 checkpoint files and the DA
# assessment output.
#
# Removal unit: entire trait column for one group.
# Remove traits where:
#   (a) DA exceeds FA4a (from the DA assessment -- see step9_highDA above), OR
#   (b) Skew or kurtosis is significant after sequential Bonferroni correction.
#
# Format: one list entry per group with traits to remove, trait as a vector.
#   Example with group_cols = c("period"):
#     list(period = "Jomon", trait = c("trait1", "trait2"))

# yayoi xli2 and xlm1 had sig kurt and xlm1 had sig skew

step9_eliminations = list(
  list(period = "Yayoi", trait = c("xli2", "xlm1"))
)