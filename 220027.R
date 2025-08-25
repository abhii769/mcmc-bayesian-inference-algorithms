log_posterior_beta <- function(beta, y, X, lambda, tau) {
  eta <- X %*% beta
  log_lik <- sum(y * eta - exp(eta))  # Poisson log-likelihood
  prior_var <- (tau^2) * (lambda^2)
  log_prior <- sum(dnorm(beta, mean = 0, sd = sqrt(prior_var), log = TRUE))
  return(log_lik + log_prior)
}

grad_log_posterior_beta <- function(beta, y, X, lambda, tau) {
  eta <- X %*% beta
  grad_ll <- t(X) %*% (y - exp(eta))
  prior_var <- (tau^2) * (lambda^2)
  grad_prior <- -beta / prior_var
  return(as.numeric(grad_ll + grad_prior))
}

poissonHorse <- function(y, X, n_iter, beta_start, lambda_start, tau_start) {
  p <- dim(X)[2]
  beta <- matrix(0, ncol = p, nrow = n_iter)
  lambda <- matrix(0, ncol = p, nrow = n_iter)
  tau <- numeric(length = n_iter)
  
  beta[1, ] <- beta_start
  lambda[1, ] <- lambda_start
  tau[1] <- tau_start
  
  step_size <- 1e-4
  
  for (t in 2:n_iter) {
    beta_prev <- beta[t - 1, ]
    lambda_prev <- lambda[t - 1, ]
    tau_prev <- tau[t - 1]
    
    grad_lp <- grad_log_posterior_beta(beta_prev, y, X, lambda_prev, tau_prev)
    
    beta_prop <- beta_prev + (step_size^2 / 2) * grad_lp + step_size * rnorm(p)
    
    if (any(is.na(beta_prop))) {
      beta[t, ] <- beta_prev
      lambda[t, ] <- lambda_prev
      tau[t] <- tau_prev
      next
    }
    
    log_post_curr <- log_posterior_beta(beta_prev, y, X, lambda_prev, tau_prev)
    log_post_prop <- log_posterior_beta(beta_prop, y, X, lambda_prev, tau_prev)
    
    grad_prop <- grad_log_posterior_beta(beta_prop, y, X, lambda_prev, tau_prev)
    
    # Handle NaNs in gradient
    if (any(is.na(grad_lp)) || any(is.na(grad_prop))) {
      beta[t, ] <- beta_prev
      lambda[t, ] <- lambda_prev
      tau[t] <- tau_prev
      next
    }
    
    log_q_curr_to_prop <- sum(dnorm(beta_prop,
                                    mean = beta_prev + (step_size^2 / 2) * grad_lp,
                                    sd = step_size, log = TRUE))
    log_q_prop_to_curr <- sum(dnorm(beta_prev,
                                    mean = beta_prop + (step_size^2 / 2) * grad_prop,
                                    sd = step_size, log = TRUE))
    
    log_accept_ratio <- (log_post_prop + log_q_prop_to_curr) - (log_post_curr + log_q_curr_to_prop)
    
    if (is.na(log_accept_ratio)) {
      beta[t, ] <- beta_prev
    } else if (log(runif(1)) < log_accept_ratio) {
      beta[t, ] <- beta_prop
    } else {
      beta[t, ] <- beta_prev
    }
    
    # --- Sample lambda and tau ---
    v <- 1 / rgamma(p, shape = 1, rate = 1 + 1 / (lambda_prev^2))
    eta <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / (tau_prev^2))
    
    rate_lambda <- 1 / v + (beta[t, ]^2) / (2 * (tau_prev^2))
    lambda2_new <- rgamma(p, shape = 1, rate = rate_lambda)
    lambda[t, ] <- sqrt(1 / lambda2_new)
    
    rate_tau <- 1 / eta + 0.5 * sum((beta[t, ]^2) / (lambda[t, ]^2))
    tau2_new <- rgamma(1, shape = (p + 1) / 2, rate = rate_tau)
    tau[t] <- sqrt(1 / tau2_new)
  }
  
  samples <- list("beta" = beta, "lambda" = lambda, "tau" = tau)
  return(samples)
}
