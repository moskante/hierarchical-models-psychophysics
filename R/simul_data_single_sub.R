library(tidyverse)
library(MixedPsy)
library(gnlm)

set.seed(123)
simul_data <- PsySimulate(ntrials = 160, nsubjects = 10, guess = TRUE, lapse = TRUE)
simul_data$Subject <- factor(simul_data$Subject)

parameters_simul <- simul_data %>%
  group_by(Subject) %>%
  summarise(across(everything(), first),
            pse = -Intercept/Slope,
            jnd = qnorm(0.75)/Slope,
            gamma = Gamma,
            lambda = Lambda) %>%
  select(Subject, pse, jnd, gamma, lambda)

PSE <- mean(parameters_simul$pse)
JND <- mean(parameters_simul$jnd)

# GLM ----
glm_list <- PsychModels(formula = cbind(Longer, Total - Longer) ~ X,
                        data = simul_data,
                        group_factors = "Subject")

parameters_glm <- PsychParameters(glm_list)

## inference tests -----
t_glm_pse <- t.test(parameters_glm$pse, mu = PSE)
t_glm_jnd <- t.test(parameters_glm$jnd, mu = JND)


# GNL ----
source("R/gnlm_functions_psejnd.R")
parameters_gnm <- simul_data %>%
  group_split(Subject) %>%
  map_dfr(process_subject)


## inference tests -----
t_gnm_pse <- t.test(parameters_gnm$pse, mu = PSE)
t_gnm_jnd <- t.test(parameters_gnm$jnd, mu = JND)



