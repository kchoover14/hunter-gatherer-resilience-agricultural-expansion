############## LIBRARIES

library(readxl)       # read_excel()
library(dplyr)        # data manipulation
library(tidyr)        # pivot_longer()
library(readr)        # write_csv()


############## PATHS
path_input    = file.path("data/raw-nagasaki-18July13.xlsx")
path_output = file.path("data/revised-fa-wrangled.csv")



############## LOOKUP TABLE FROM PO
lookup = read_excel(path_input, sheet = "demog", col_names = TRUE)

lookup_key = lookup |>
  mutate(key = gsub("wakimiksaki", "wakimisaki", tolower(paste(site, ind, sep = "_"))),
         key = gsub("[ -]", "_", key)) |>
  select(key, period) |>
  mutate(period = ifelse(grepl("Jomon", period), "Jomon", period),
         period = ifelse(grepl("yayoi", tolower(period)), "Yayoi", period))



############## USER CONFIGURATION
# Name of the Side column (values should be "R" and "L"):
col_side = "side"

col_side_right = "R"
col_side_left  = "L"

# Name of the Replicate column in your sheets (read_sheets_with_rep), or
# the name to assign to the inferred column (read_sheets_infer_rep):
col_replicate = "rep"

# Names of all trait (measurement) to be used by enforce_min_reps() to identify paired measurements.
trait_cols = c('xlm3', 'xlm2', 'xlm1', 'xlp2', 'xlp1', 'xlc', 'xli2', 'xli1', 'xbm3', 'xbm2', 'xbm1', 'xbp2', 'xbp1', 'xbc', 'xbi2', 'xbi1', 'mlm3', 'mlm2', 'mlm1', 'mlp2', 'mlp1', 'mlc', 'mli2', 'mli1', 'mbm3', 'mbm2', 'mbm1', 'mbp2', 'mbp1', 'mbc', 'mbi2', 'mbi1'
)

# Sheet names to read as replicate trials.
rep_sheets = c("fa1", "fa2")

col_side       = "side"
col_replicate  = "rep"
col_side_right = "R"
col_side_left  = "L"





#### FUNCTIONS
# read_sheets_with_rep
read_sheets_with_rep = function(path, sheets, col_rep = col_replicate) {
  # Use when each replicate sheet already contains a Replicate column.
  # Reads the specified sheets and binds them into a single dataframe.
  # Stops with an informative error if the Replicate column is not found
  # and suggests switching to read_sheets_infer_rep().
  #
  # Arguments:
  #   path    -- path to the Excel workbook
  #   sheets  -- numeric or character vector of sheet positions or names
  #              (use rep_sheets from USER CONFIGURATION)
  #   col_rep -- name of the Replicate column (default: col_replicate)
  #
  # Returns: combined dataframe, column structure identical to input sheets.
  all_sheets = list()

  for (s in sheets) {
    sheet = read_excel(path, sheet = s, col_names = TRUE,
                       .name_repair = "minimal")

    if (!col_rep %in% names(sheet)) {
      stop(sprintf(paste0(
        "Sheet '%s': Replicate column '%s' not found.\n",
        "  Columns present: %s\n",
        "  If your sheets have no Replicate column, use read_sheets_infer_rep() instead."
      ), s, col_rep, paste(names(sheet), collapse = ", ")))
    }

    # Coerce trait columns to numeric.
    present_traits = intersect(trait_cols, names(sheet))
    sheet[present_traits] = lapply(sheet[present_traits],
                                   function(x) suppressWarnings(as.numeric(x)))

    all_sheets[[as.character(s)]] = sheet
  }

  out = bind_rows(all_sheets)
  cat(sprintf("read_sheets_with_rep: %d sheets read, %d rows total\n",
              length(sheets), nrow(out)))
  out
}



# enforce_min_reps -------------------------------------------------------
enforce_min_reps = function(df,
                            col_ind = col_individual,
                            col_s   = col_side,
                            col_rep = col_replicate) {
  # Operates on long-format data (Trait and Value columns must be present).
  # For each individual x trait x rep, checks whether both R and L are
  # non-NA (a paired rep). If fewer than 2 paired reps exist for an
  # individual x trait combination, sets Value to NA for all rows of
  # that combination.
  #
  # Palmer & Strobeck (2003) require a minimum of 2 paired replicates for
  # the Sides x Individuals ANOVA. This function enforces that requirement
  # at the wrangling stage so no downstream step receives insufficient data.
  #
  # Must be called AFTER pivot_to_long() since it requires Trait and Value
  # columns.
  #
  # Arguments:
  #   df      -- long-format dataframe with Trait and Value columns
  #   col_ind -- name of the Individual column (default: col_individual)
  #   col_s   -- name of the Side column (default: col_side)
  #   col_rep -- name of the Replicate column (default: col_replicate)
  #
  # Returns: dataframe with insufficient Individual x Trait combinations
  #          set to NA. All other columns are unchanged.

  paired_counts = df |>
    group_by(.data[[col_ind]], Trait, .data[[col_rep]]) |>
    summarise(
      has_R = any(.data[[col_s]] == col_side_right & !is.na(Value)),
      has_L = any(.data[[col_s]] == col_side_left  & !is.na(Value)),
      .groups = "drop"
    ) |>
    mutate(paired = has_R & has_L) |>
    group_by(.data[[col_ind]], Trait) |>
    summarise(n_paired = sum(paired), .groups = "drop")

  insufficient = paired_counts |>
    filter(n_paired < 2) |>
    select(all_of(col_ind), Trait)

  if (nrow(insufficient) == 0) {
    cat("enforce_min_reps: all Individual x Trait combinations have >= 2 paired replicates.\n")
    return(df)
  }

  cat(sprintf(
    "enforce_min_reps: %d Individual x Trait combination(s) have < 2 paired replicates -- setting to NA.\n",
    nrow(insufficient)
  ))

  df = df |>
    left_join(insufficient |> mutate(drop = TRUE),
              by = c(col_ind, "Trait")) |>
    mutate(Value = ifelse(!is.na(drop), NA_real_, Value)) |>
    select(-drop)

  df
}






# WRANGLE

# read all replicate sheets and combine into one dataframe.
raw = read_sheets_with_rep(path_input, sheets = rep_sheets)


# create composite individual ID (site + ind) for uniqueness across sites.
raw = raw |>
  mutate(uid = paste(site, ind, sep = "_"))
col_individual = "uid"


# recode sex to full labels (wide format).
raw = raw |>
  mutate(sex = case_when(
    tolower(sex) == "f" ~ "female",
    tolower(sex) == "m" ~ "male"
  ))

# pivot to long format for inspection pipeline.
raw = pivot_longer(raw,
                   cols      = all_of(trait_cols),
                   names_to  = "Trait",
                   values_to = "Value")

# enforce minimum paired replicates (sets insufficient combos to NA).
raw = enforce_min_reps(raw)

# remove rows where Value is NA.
raw = raw |>
  filter(!is.na(Value))

# assign period from lookup table using match -- no join, no row multiplication.
raw = raw |>
  mutate(
    key    = gsub("[ -]", "_", gsub("wakimiksaki", "wakimisaki", tolower(paste(site, ind, sep = "_")))),
    period = lookup_key$period[match(key, lookup_key$key)]
  ) |>
  select(-key) |>
  relocate(period, .after = site)

# remove all third molars (trait contains m3)
raw = raw |> filter(!grepl("m3", Trait))

# write the combined CSV.
write_csv(raw, path_output)

######################## TIDY
rm(list = ls())
gc()
