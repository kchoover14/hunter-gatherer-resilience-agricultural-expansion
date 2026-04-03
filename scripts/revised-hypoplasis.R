############## LIBRARIES
library(readxl)       # read_excel()
library(dplyr)        # data manipulation
library(tidyr)        # fill()
library(ggplot2)      # data visualization
library(cowplot)      # plot_grid() for panel figures
library(car)          # leveneTest()
library(MuMIn)        # r.squaredGLMM()


############## PATHS
path_input  = file.path("data/raw-nagasaki-18July13.xlsx")





############## LOOKUP TABLE FROM demog
lookup = read_excel(path_input, sheet = "demog", col_names = TRUE)

lookup_key = lookup |>
  mutate(
    uid = gsub("[ -]", "_", gsub("wakimiksaki", "wakimisaki", tolower(paste(site, ind, sep = "_"))))
  ) |>
  dplyr::select(uid, period, sex, age, nTeeth)



############## FUNCTIONS

# tooth number to trait name map (universal dental numbering system)
tooth_map = c(
  "1"  = "xm3r", "2"  = "xm2r", "3"  = "xm1r", "4"  = "xp2r", "5"  = "xp1r",
  "6"  = "xcr",  "7"  = "xi2r", "8"  = "xi1r", "9"  = "xi1l", "10" = "xi2l",
  "11" = "xcl",  "12" = "xp1l", "13" = "xp2l", "14" = "xm1l", "15" = "xm2l",
  "16" = "xm3l",
  "17" = "mm3l", "18" = "mm2l", "19" = "mm1l", "20" = "mp2l",
  "21" = "mp1l", "22" = "mcl",  "23" = "mi2l", "24" = "mi1l", "25" = "mi1r",
  "26" = "mi2r", "27" = "mcr",  "28" = "mp1r", "29" = "mp2r", "30" = "mm1r",
  "31" = "mm2r", "32" = "mm3r"
)

# read_hypo_sheet
read_hypo_sheet = function(path, sheet) {
  raw = read_excel(path, sheet = sheet, col_names = TRUE,
                   .name_repair = "minimal")

  # fill down site and ind within individual blocks
  raw = raw |>
    mutate(
      site = ifelse(site == "" | is.na(site), NA, site),
      ind  = ifelse(ind  == "" | is.na(ind),  NA, ind)
    ) |>
    fill(site, ind, .direction = "down")

  # coerce tooth columns to numeric (handles Pit entries -> NA)
  tooth_nums = as.character(1:32)
  present_teeth = intersect(tooth_nums, names(raw))
  raw[present_teeth] = lapply(raw[present_teeth],
                              function(x) suppressWarnings(as.numeric(x)))
  raw
}



############## WRANGLE

# read hypo sheet
hypo = read_hypo_sheet(path_input, "hypo")

# create uid
hypo = hypo |>
  mutate(
    uid = gsub("[ -]", "_", gsub("wakimiksaki", "wakimisaki", tolower(paste(site, ind, sep = "_"))))
  )

# rename tooth columns using dental map
hypo = hypo |> rename(any_of(setNames(names(tooth_map), tooth_map)))

# remove 0 values
hypo = hypo |> mutate(across(all_of(unname(tooth_map)), ~ ifelse(. == 0, NA, .)))

# create final dataset
hypo = hypo |>
  group_by(uid) |>
  summarise(
    hypoCount = sum(colSums(!is.na(across(all_of(unname(tooth_map)))))),
    .groups   = "drop"
  ) |>
  left_join(lookup_key, by = "uid") |>
  mutate(period = case_when(
    grepl("Jomon", period) ~ "Jomon",
    tolower(period) == "yayoi" ~ "Yayoi"
  )) |>
  mutate(
    hypoRate = hypoCount / nTeeth
  ) |>
  select(-hypoCount) |>
  relocate(period, .after = uid) |>
  relocate(sex, .after = period) |>
  relocate(age, .after = sex) |>
  relocate(nTeeth, .after = age)

rm(path_input, lookup, lookup_key, tooth_map, read_hypo_sheet)





############## DISTRIBUTIONS
# summary stats by period
aggregate(hypoRate ~ period, hypo, function(x) c(mean = mean(x), var = var(x)))

# hypoRate
a = ggplot(hypo, aes(x = hypoRate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 10) +
  geom_vline(aes(xintercept = mean(hypoRate, na.rm = TRUE)),
             color = "black", linetype = "dashed", linewidth = 1) +
  geom_density(alpha = .2) +
  annotate("text", x = mean(hypo$hypoRate, na.rm = TRUE), y = max(density(hypo$hypoRate, na.rm = TRUE)$y),
           label = round(mean(hypo$hypoRate, na.rm = TRUE), 2),
           color = "black", size = 4, hjust = -0.1) +
  theme_classic() + xlab("hypoRate") + ylab("Density")

b = ggplot(hypo, aes(x = hypoRate, color = period)) +
  geom_histogram(aes(y = ..density..), fill = "white", position = "dodge", bins = 10) +
  geom_density(alpha = 0.2) +
  scale_color_viridis_d(end = .6) +
  theme_classic() + xlab("hypoRate") + ylab("Density")

plot_grid(a, b, nrow = 1, labels = c("A", "B"))
ggsave("revised-fig-hypo-distribution.png", plot = last_plot(), width = 10, height = 4, dpi = 300)

# levene test for variance homogeneity by period
car::leveneTest(hypoRate ~ period, data = hypo)

# all distributions OK and levene's n.s. so can proceed with models




############## ANALYSIS
# offset
hypo$tooth_offset = log(hypo$nTeeth / 28)

# if random effect variance is zero, drop to lm
lm_rate = lm(hypoRate ~ period  + offset(tooth_offset), data = hypo)
summary(lm_rate)

# R2
r.squaredGLMM(lm_rate)

# visualize hypoRate by period
means_rate = hypo |>
  group_by(period) |>
  summarise(mean = mean(hypoRate, na.rm = TRUE))

ggplot(hypo, aes(x = period, y = hypoRate, color = period)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.6) +
  geom_point(data = means_rate, aes(y = mean), shape = 18, size = 4, color = "black") +
  scale_color_viridis_d(end = .6) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") + ylab("Hypoplasia Rate")
ggsave("revised-fig-hypoRate.png", plot = last_plot(), units = "in", height = 5, width = 4, dpi = 300)

# jomon trends lower but not sig diff


######################## TIDY
rm(list = ls())
gc()