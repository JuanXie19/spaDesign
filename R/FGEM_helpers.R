library(movMF)
library(matrixStats)
library(Directional)


#' @keywords internal
# L2 norm function
norm_vec <- function(v) {
  norm_val <- sqrt(sum(v^2))
  if (norm_val < 1e-10) return(1e-10)  # Avoid division by zero
  return(norm_val)
}

#' @keywords internal
#' Function to calculate log(Cd(tau))
# Cd(tau) = tau^{d/2-1}(2pi)^{-d/2}(I_{d/2-1}(tau))^{-1}
# besselI() function in base R package computes the modified Bessel function of the first kind

log_Cd <- function(tau, d) {
  # Check for invalid inputs
  if (tau < 0) stop("tau must be non-negative")
  
  # Handle the case where tau is zero
  if (tau == 0) {
    # For tau = 0, the vMF distribution is uniform on the sphere
    # The normalizing constant is C_d(0) = gamma(d/2)/(2*pi^(d/2))
    # So log(C_d(0)) = log(gamma(d/2)) - log(2) - (d/2)*log(pi)
    # lgamma computes the logarithm of the gamma function
    return(-log(2) - (d/2) * log(pi) + lgamma(d/2))
  }
  
  # Compute log_Cd(tau) for tau > 0
  log_tau <- log(tau)
  bessel_term <- besselI(tau, nu = (d/2 - 1), expon.scaled = TRUE)
  
  # Compute log_Cd(tau)
  res <- (d/2 - 1) * log_tau - log(bessel_term) - (d/2) * log(2*pi)
  return(res)
}

#' @keywords internal
# Function to calculates FG-kernel density given the parameters cntr(center), rd(radius), sigma.sq, phi(direction), tau(concentration)
# this comes from the FG paper code. Note that we 
FG_kernel <- function(x_i, c_k, r_k, sigma_sq, phi_k, tau_k) { 
  d <- length(x_i)
  
  # Check for invalid inputs
  if (sigma_sq <= 0) stop("sigma_sq must be positive")
  if (tau_k <= 0) stop("tau_k must be positive")
  
  # Compute z = x_i - c_k
  z <- x_i - c_k
  
  # Compute temp = ||tau_k * phi_k + (r_k * z) / sigma_sq||
  temp <- norm_vec(tau_k * phi_k + (r_k * z) / sigma_sq)
  
  
  # Compute the log-density
  #ln_h_x_theta <- log_Cd_tau - (d * log(2 * pi * sigma_sq) / 2) - 
  #  log_Cd_temp - ((norm_vec(z))^2 + (r_k)^2) / (2 * sigma_sq)
  
  # when expon.scaled = T, need to adjust the extra terms due to scaling
  ln_h_x_theta <- log_Cd(tau_k, d) - (d * log(2 * pi * sigma_sq) / 2) - 
    log_Cd(temp, d) - ((norm_vec(z))^2 + (r_k)^2) / (2 * sigma_sq) + temp - tau_k
  
  # Handle numerical underflow by returning exp(ln_h_x_theta)
  return(exp(ln_h_x_theta))
}

 
#' @keywords internal
# Function to compute E(y_i|x_i,z_i, Theta)
compute_y_expect <- function(x, M, phi, tau, c, r, sigma_sq) {
  n <- nrow(x)
  d <- ncol(x)
  
  E_y <- list()
  
  for (k in 1:M) {
    # Compute nu for all points in the cluster at once
    nu <- t(tau[k] * phi[k, ] + t(r[k] * (x - matrix(c[k, ], n, d, byrow = TRUE)) / sigma_sq))
    
    # Compute kappa_k and mu_k
    kappa_k <- apply(nu, 1, norm_vec)
    # Handle cases where kappa_k is zero
    kappa_k <- pmax(kappa_k, 1e-10)  # Small value to avoid division by zero
    mu_k <- nu / kappa_k
    
   
    # Compute A_d(kappa_k) using exponentially scaled Bessel functions for numerical stability
    log_Ad_k <- log(besselI(kappa_k, d/2, expon.scaled = TRUE)) - 
      log(besselI(kappa_k, d/2 - 1, expon.scaled = TRUE))
    Ad_k <- exp(log_Ad_k)
    Ad_k <- pmin(pmax(Ad_k, 1e-10), 1 - 1e-10)
    
    # Compute E[y_i | x_i, z_i = k]
    y_expect_k <- Ad_k * mu_k
    
    E_y[[k]] <- y_expect_k
  }
  
  return(E_y)
}


# Function to calculate the incomplete log-likelihood
### logSumExp calculate the log of the sum of exponentials
## eg., originally, we calculate pi_1*FG1, pi_2*FG2, etc, then we calculate log(pi_1*FG1 + pi_2*FG2)
## now we calculate t1= log(pi_1*FG1), t2= log(pi_2*FG2), etc, then logSumExp(t1, t2), which is log(exp(t1)+ exp(t2))
## which is just log(pi_1*FG1 + pi_2*FG2)


log_likelihood <- function(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq) {
  n <- nrow(x)
  M <- length(pi_hat)
  log_lik <- matrix(0, n, M)
  temp <- numeric(n)
  
  for (i in 1:n) {
    for (k in 1:M) {
      fg_value <- FG_kernel(x[i, ], c_hat[k, ], r_hat[k], sigma_sq, phi_hat[k, ], tau_hat[k])
      # Handle extreme cases robustly
      if (fg_value <= 0 || !is.finite(fg_value)) {
        if (fg_value == 0 || is.na(fg_value)) {
          warning(sprintf("FG_kernel=0 at i=%d, k=%d. Check for outliers or mis-specified M.", i, k))
        }
        fg_value <- 1e-100  # Fallback to tiny positive value
      }
      
      log_lik[i, k] <- log(pi_hat[k]) + log(fg_value)
    }
    temp[i] <- logSumExp(log_lik[i, ])  # Use logSumExp for numerical stability
    
  }
  
  rst <- sum(temp)
  return(rst)
}

## function to compute the posterior probability w_{ik}= P(z_i=k|x_i, Theta)
compute_W <- function(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq){
  n <- nrow(x)  # Number of data points
  M <- length(pi_hat)  # Number of clusters
  
  W <- matrix(0, n, M)  
  for(i in 1:n){
    for (k in 1: M){
      # compute p(x_i|z_i=k, Theta)
      p_x_given_z <- FG_kernel(x[i,], c_hat[k, ], r_hat[k], sigma_sq, phi_hat[k, ], tau_hat[k])
      W[i, k] <- pi_hat[k] * p_x_given_z
    }
  }
  return(W)
}

