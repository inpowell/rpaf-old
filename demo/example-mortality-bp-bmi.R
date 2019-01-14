mpaf_bpxbmi <- mpaf_gen_data(minifhs, "ID", c(0,5,10,15,17), "DEATH", "DEATH_FT",
                             c("B_COHORT", "SEX", "BP", "BMI_2"))
summary_bpxbmi <- mpaf_summary(
  survival::Surv(f_end, DEATH) ~ B_COHORT * f_period + SEX + BP + BP:BMI_2,
  mpaf_data = mpaf_bpxbmi,
  modifications = list(BMI_2 = "<25.0"),
  covar_model = c("SEX", "BP", "BP:BMI_2"),
  group_vars = "BP"
)

summary_bpxbmi
