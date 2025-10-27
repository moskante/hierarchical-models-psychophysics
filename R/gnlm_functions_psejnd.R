mu <- function(p){
  mu_i <- p[1]
  sigma_i <- p[2]
  gamma_i <- atan2(p[3],1)/pi + 0.5
  lambda_i <- (1-gamma_i)*(atan2(p[4],1)/pi + 0.5)
  gamma_i + (1 - gamma_i - lambda_i) * pnorm(x, mean = mu_i, sd = sigma_i)
}

pstart <- function(x,y){
  glm_start <- glm(y~x, family=binomial(link = "probit"))
  coef_start <- coef(glm_start)
  if (is.na(coef_start[2]) || coef_start[2] == 0) {
    coef_start[2] <- 1e-6
  }
  prop_correct <- y[,1]/rowSums(y)
  gamma_i_start <- min(min(prop_correct, na.rm = TRUE),0.99)
  lambda_i_start <- max(1 - max(prop_correct, na.rm = TRUE), 0.01)
  mu_i_start <- -coef_start[1]/coef_start[2]
  sigma_i_start <- 1/coef_start[2]
  p3_start <- tan(pi*(gamma_i_start-0.5))
  p4_start <- tan(pi*((lambda_i_start/(1-gamma_i_start))-0.5))
  return(c(mu_i_start, sigma_i_start, p3_start, p4_start))
}

PsychParametersGNM <- function(gnm, p=0.75){
  pout <- gnm$coefficients
  pse <- pout[1]
  jnd <- qnorm(p)*pout[2]
  gamma_i <- atan2(pout[3],1)/pi + 0.5
  lambda_i <- (1-gamma_i)*(atan2(pout[4],1)/pi + 0.5)
  return(c(pse = pse, jnd = jnd, gamma = gamma_i, lambda = lambda_i))
}

PsychBootGNM <- function(data,indices, p=0.75){
  data_resampled <- data[indices,]
  y <- with(data_resampled, cbind(Longer, Total-Longer))
  x <- assign("x", data_resampled$X, envir = .GlobalEnv)
  pmu <- pstart(x,y)
  gnm <- tryCatch(
    {model <- gnlr(y = y, distribution = "binomial",
                   mu = mu, pmu = pmu, iterlim = 10000)
    PsychParametersGNM(model, p)
    },error = function(e) rep(NA, 4) 
  )
  x <- assign("x", data$X, envir = .GlobalEnv)
  return(gnm)
}

process_subject <- function(sub_data) {
  y <- with(sub_data, cbind(Longer, Total - Longer))
  x <- assign("x", sub_data$X, envir = .GlobalEnv)
  
  pmu <- pstart(x, y)
  
  fit_gnm <- gnlr(y = y,
                  distribution = "binomial",
                  mu = mu, pmu = pmu, iterlim = 10000)
  
  boot_results <- boot::boot(
    data = sub_data,
    statistic = PsychBootGNM,
    R = 500,
    p = 0.75 #change for JND
  )
  
  estimates <- PsychParametersGNM(fit_gnm, p = 0.75)
  ses <- apply(boot_results$t, 2, sd, na.rm = TRUE)
  
  # Turn estimates and ses into a single row tibble
  param_tibble <- as_tibble(t(estimates))
  colnames_se <- c("pse_se", "jnd_se", "gamma_se", "lambda_se")
  se_tibble <- as_tibble(t(ses))
  colnames(se_tibble) <- colnames_se
  
  bind_cols(param_tibble, se_tibble, tibble(Subject = unique(sub_data$Subject)))
}

