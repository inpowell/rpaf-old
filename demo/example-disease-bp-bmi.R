dpaf_bpxbmi <- gen_data(
  minifhs,
  "ID",
  c(0,5,10,15,17),
  "DEATH", "DEATH_FT", "DIAB", "DIAB_FT",
  c("B_COHORT", "SEX", "BP", "BMI_2")
)

summary_bpxbmi <- dpaf_summary(
  survival::Surv(f_end, DIAB) ~ .,
  survival::Surv(f_end, DEATH) ~ .,
  ~ 0 + B_COHORT * f_period + SEX + BP + BP:BMI_2,
  dpaf_data = dpaf_bpxbmi,
  modifications = list(BMI_2 = "<25.0"),
  covar_model = c("SEX", "BP", "BP:BMI_2"),
  group_vars = "BP"
)

summary_bpxbmi
