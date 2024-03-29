mpaf_smoke <- gen_data(minifhs, "ID", c(0,5,10,15,17), "DEATH", "DEATH_FT",
                       variables = c("SEX", "SMOKE", "B_COHORT"))
summary_smoke <- mpaf_summary(
  survival::Surv(f_end, DEATH) ~ B_COHORT * f_period + SEX + SMOKE,
  paf_data = mpaf_smoke,
  modifications = list("SMOKE" = c("Never", "<30/day", ">=30/day")),
  hr_out = c("SEX", "SMOKE")
)

summary_smoke
