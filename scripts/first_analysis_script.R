# load packages
library(tidyverse)
library(gamm4)

# import data
VASCO_data_org <- readRDS(file = "data/20230609_data_VASCO_CIS.rds")

# give structure and summary
VASCO_data_org |> str()
VASCO_data_org |> summary()

# transform data
VASCO_data <- VASCO_data_org |> 
  rename_with(
    .fn = tolower
  ) |> 
  mutate(
    study_id = studyid |> 
      factor(),
    studyid = NULL,
    fatigue_score = cis_subj_fat,
    cis_subj_fat = NULL,
    variant = variant |> 
      word(start = 1) |> 
      factor()
  )

VASCO_data_omicron_delta <- VASCO_data |> 
  filter(variant %in% c("Delta", "Omicron")) |> 
  droplevels()

# fixed effects only, P-spline for time
omicron_delta_gam_ps <- gam(
  formula = fatigue_score ~ variant + s(time, bs = "ps", by = variant),
  data = VASCO_data_omicron_delta)

# fixed effects only, adaptive-spline for time
omicron_delta_gam_ad <- gam(
  formula = fatigue_score ~ variant + s(time, bs = "ad", by = variant),
  data = VASCO_data_omicron_delta)
plot(omicron_delta_gam_ps)
plot(omicron_delta_gam_ad)
# Doesn't work/converge.....
# omicron_delta_gam_re <- gam(
#   formula = fatigue_score ~ variant + s(time, bs = "ps", by = variant) + s(study_id, bs = "re"),
#   data = VASCO_data_omicron_delta)

omicron_delta_gamm4 <- gamm4(
  formula = fatigue_score ~ variant + s(time, bs = "ps", by = variant),
  random = ~(1|study_id),
  REML = FALSE,
  data = VASCO_data_omicron_delta)

plot(omicron_delta_gamm4$gam, pages = 1)

summary(omicron_delta_gamm4$gam) ## summary of gam
summary(omicron_delta_gamm4$mer) ## underlying mixed model
anova(omicron_delta_gamm4$gam) 

omicron_delta_gamm4_0 <- gamm4(
  formula = fatigue_score ~ variant + s(time, bs = "ps"),
  random = ~(1|study_id),
  REML = FALSE,
  data = VASCO_data_omicron_delta)

# test for interaction effect between spline and variant
anova(omicron_delta_gamm4$mer, omicron_delta_gamm4_0$mer)

B <- omicron_delta_gamm4$gam$gcv.ubre
A <- omicron_delta_gamm4_0$gam$gcv.ubre


df_B <- omicron_delta_gamm4$gam |> logLik() |> attributes() |> pluck("df")
df_A <- omicron_delta_gamm4_0$gam |> logLik() |> attributes() |> pluck("df")

pchisq(
  q = A-B, 
  df = df_B-df_A,
  lower.tail = FALSE)

# omicron_delta_gamm <- gamm(
#   formula = fatigue_score ~ variant + s(time, bs = "ps", by = variant),
#   random = list(study_id = ~1), 
#   data = VASCO_data_omicron_delta)

newdata <- expand_grid(
  variant = VASCO_data_omicron_delta$variant |> fct_unique(),
  time = -90:300
)

newdata <- newdata |> 
  mutate(
    pred_gamm4 = predict(
      object  = omicron_delta_gamm4$gam,
      newdata = newdata)
  )


newdata |> 
  ggplot(
    mapping = aes(x = time, color = variant)
  ) +
  geom_line(
    mapping = aes(y = pred_gamm4)
  )


# Make `lpmatrix` predictions for calculating differences and confidence bands thereof
newdata <- tibble(
  time = -90:300
)

X_lp_omicron <- predict(
  object  = omicron_delta_gamm4$gam,
  newdata = newdata |> 
    mutate(variant = "Omicron"),
  type = "lpmatrix")

X_lp_delta <- predict(
  object  = omicron_delta_gamm4$gam,
  newdata = newdata |> 
    mutate(variant = "Delta"),
  type = "lpmatrix")

X_lp_difference <- X_lp_delta - X_lp_omicron

betas <- coef(omicron_delta_gamm4$gam)
varcov <- omicron_delta_gamm4$gam$Vp

newdata <- newdata |> 
  mutate(
    # calculate prediction of difference
    pred_difference = X_lp_difference %*% betas,
    # calculate standard errors of prediction of difference
    se_pred_difference =  (X_lp_difference %*% varcov %*% t(X_lp_difference)) |> 
      diag() |> 
      sqrt(),
    # calculate confidence limits of prediction of difference
    ci_lwr_pred_difference = pred_difference - 1.96 * se_pred_difference,
    ci_upr_pred_difference = pred_difference + 1.96 * se_pred_difference)

# visualize predicted differences   
ggplot(
  data = newdata,
  mapping = aes(x = time, y = pred_difference)
) + 
  # add predictions
  geom_line(
  ) +
  # add confidence bands
  geom_ribbon(
    mapping = aes(
      ymin = ci_lwr_pred_difference,
      ymax = ci_upr_pred_difference),
    alpha = .3) +
  # add reference line at 0
  geom_hline(yintercept = 0, linetype = 2) +
  # add annotations
  labs(
    title ="predicted differences between Delta and Omicron variants",
    subtitle = "with 95% confidence limits",
    y = "fatigue_score difference"
  )

# add season effect (as cyclic spline)
VASCO_data_omicron_delta <- VASCO_data_omicron_delta |> 
  mutate(
    maand_infectie = infectie_datum |> 
      lubridate::month()
  )

omicron_delta_season_gamm4 <- gamm4(
  formula = fatigue_score ~ variant + s(time, bs = "ps", by = variant) + s(maand_infectie, bs = "cp", k = 10),
  knots = list(
    maand_infectie = c(0, 12)),
  random = ~(1|study_id),
  REML = TRUE,
  data = VASCO_data_omicron_delta)

plot(omicron_delta_season_gamm4$gam, pages = 1)


C <- omicron_delta_season_gamm4$gam$gcv.ubre
B <- omicron_delta_gamm4$gam$gcv.ubre
A <- omicron_delta_gamm4_0$gam$gcv.ubre

df_C <- omicron_delta_season_gamm4$gam |> logLik() |> attributes() |> pluck("df")
df_B <- omicron_delta_gamm4$gam |> logLik() |> attributes() |> pluck("df")
df_A <- omicron_delta_gamm4_0$gam |> logLik() |> attributes() |> pluck("df")

pchisq(
  q = B-C, 
  df = df_C-df_B,
  lower.tail = FALSE)
