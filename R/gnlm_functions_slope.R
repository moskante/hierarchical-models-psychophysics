mu <- function(p){
  beta_0 <- p[1]
  beta_1 <- p[2]
  gamma_i <- atan2(p[3],1)/pi + 0.5
  lambda_i <- (1-gamma_i)*(atan2(p[4],1)/pi + 0.5)
  gamma_i + (1 - gamma_i - lambda_i) * pnorm(beta_0 + beta_1*x)
}

pstart <- function(x,y){
  glm_start <- glm(y~x, family=binomial(link = "probit"))
  coef_start <- coef(glm_start)
  prop_correct <- y[,1]/rowSums(y)
  gamma_i_start <- min(min(prop_correct, na.rm = TRUE),0.99)
  lambda_i_start <- max(1 - max(prop_correct, na.rm = TRUE), 0.01)
  beta_0_start <- coef_start[1]
  beta_1_start <- coef_start[2]
  p3_start <- tan(pi*(gamma_i_start-0.5))
  p4_start <- tan(pi*((lambda_i_start/(1-gamma_i_start))-0.5))
  return(c(beta_0_start, beta_1_start, p3_start, p4_start))
}

PsychParametersGNM <- function(gnm, p=0.75){
  pout <- gnm$coefficients
  pse <- -pout[1]/pout[2]
  slope <- pout[2]
  gamma_i <- atan2(pout[3],1)/pi + 0.5
  lambda_i <- (1-gamma_i)*(atan2(pout[4],1)/pi + 0.5)
  return(c(pse = pse, slope = slope, gamma = gamma_i, lambda = lambda_i))
}

PsychBootGNM_vibro <- function(data,indices, p=0.75){
  data_resampled <- data[indices,]
  y <- with(data_resampled, cbind(faster, slower))
  x <- assign("x", data_resampled$speed, envir = .GlobalEnv)
  pmu <- pstart(x,y)
  gnm <- tryCatch(
    {model <- gnlr(y = y, distribution = "binomial",
                   mu = mu, pmu = pmu, iterlim = 10000)
    PsychParametersGNM(model, p)
    },error = function(e) rep(NA, 4) 
  )
  x <- assign("x", data$speed, envir = .GlobalEnv)
  return(gnm)
}

process_subject_vibro <- function(sub_data) {
  y <- with(sub_data, cbind(faster, slower))
  x <- assign("x", sub_data$speed, envir = .GlobalEnv)
  
  pmu <- pstart(x, y)
  
  fit_gnm <- gnlr(y = y,
                  distribution = "binomial",
                  mu = mu, pmu = pmu, iterlim = 10000)
  
  boot_results <- boot::boot(
    data = sub_data,
    statistic = PsychBootGNM_vibro,
    R = 500,
    p = 0.75 #change for JND
  )
  
  estimates <- PsychParametersGNM(fit_gnm, p = 0.75)
  ses <- apply(boot_results$t, 2, sd, na.rm = TRUE)
  low <- apply(boot_results$t, 2, function(x){quantile(x, 0.025, na.rm = TRUE)})
  high <- apply(boot_results$t, 2, function(x){quantile(x, 0.975, na.rm = TRUE)})
  # Turn estimates and ses into a single row tibble
  param_tibble <- as_tibble(t(estimates))
  colnames_se <- c("pse_se", "slope_se", "gamma_se", "lambda_se")
  se_tibble <- as_tibble(t(ses))
  colnames(se_tibble) <- colnames_se
  
  colnames_low <- c("pse_low", "slope_low", "gamma_low", "lambda_low")
  low_tibble <- as_tibble(t(low))
  colnames(low_tibble) <- colnames_low
  
  colnames_high <- c("pse_high", "slope_high", "gamma_high", "lambda_high")
  high_tibble <- as_tibble(t(high))
  colnames(high_tibble) <- colnames_high
  
  bind_cols(param_tibble, se_tibble, low_tibble, high_tibble, tibble(subject = unique(sub_data$subject), vibration = unique(sub_data$vibration)))
}

