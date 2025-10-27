library(tidyverse)
library(MixedPsy)
library(gnlm)

load("simul_data.RData") 
simul_data$Subject <- factor(simul_data$Subject)

# GLM ----
glm_list <- PsychModels(formula = cbind(Longer, Total - Longer) ~ X,
                        data = simul_data,
                        group_factors = "Subject")

parameters_glm <- PsychParameters(glm_list)

## inference tests -----
t_glm_pse <- t.test(parameters_glm$pse, mu = PSE)
t_glm_jnd <- t.test(parameters_glm$jnd, mu = JND)


# GNL ----
source("gnlm_functions_psejnd.R")
parameters_gnm <- simul_data %>%
  group_split(Subject) %>%
  map_dfr(process_subject)


## inference tests -----
t_gnm_pse <- t.test(parameters_gnm$pse, mu = PSE)
t_gnm_jnd <- t.test(parameters_gnm$jnd, mu = JND)



