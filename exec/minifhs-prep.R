# The original data frame we will work with has two time
# columns -- one for diabetes incidence (DIAB_FT), and one
# for death (DEATH_FT).
head(minifhs)

# we are interested in the first of diabetes incidence and death
ct_minifhs <- collapse_times(minifhs, "DIAB", "DIAB_FT", "DEATH", "DEATH_FT")
head(ct_minifhs)

# Now we can time-vary the data
ext_minifhs <- gen_data_fun(ct_minifhs, c(0,5,10,15,17), time = "time",
                            primary = "DIAB", secondary = "DEATH",
                            id_col = "ID")

# Fit survival regression models for every event, and then summarise
mods_minifhs <- fit_models(~ B_COHORT * period + SEX + BMI_2, ext_minifhs)
mods_minifhs
