library(tidyverse)
library(MixedPsy)
library(gnlm)

load("vibro_exp3.RData")

# GLM ----
model_list_vibro <- PsychModels(vibro_exp3, 
                                group_factors = c("subject", "vibration"),
                                formula = cbind(faster, slower) ~ speed)
parameters_glm <- PsychParameters(model_list_vibro)

parameters_glm <- parameters_glm %>%
  mutate(
    slope    = map_dbl(model_list_vibro, ~ .x$model$coefficients["speed"]),
    slope_se = map_dbl(model_list_vibro, ~ summary(.x$model)$coefficients["speed", "Std. Error"])
  )

## inference tests -----
params_glm_wider <- parameters_glm %>%
  select(-pse_se, -jnd_se, -slope_se)%>%
  pivot_wider(names_from = vibration, values_from = c(pse,jnd, slope))

t_glm_slope <-  t.test(params_glm_wider$slope_32, params_glm_wider$slope_0,paired = TRUE)

# GNM ---- 
source("gnlm_functions_slope.R")
parameters_gnm <- vibro_exp3 %>%
  group_split(subject, vibration) %>%
  map_dfr(process_subject_vibro)

parameters_gnm$jnd <- qnorm(0.75)/parameters_gnm$slope

params_gnm_wider <- parameters_gnm %>%
  select(pse, jnd, slope, vibration, subject)%>%
  pivot_wider(names_from = vibration, values_from = c(pse,jnd,slope))

t_gnm_slope <-  t.test(params_gnm_wider$slope_32, params_gnm_wider$slope_0,paired = TRUE)

