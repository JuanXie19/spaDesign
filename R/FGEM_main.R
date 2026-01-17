
#' EM algorithm for simplified Fisher-Gaussian model
#'
#' Fits a Fisher-Gaussian mixture model to spatial coordinates using the EM algorithm.
#'
#' @param x Numeric matrix of coordinates (rows = observation, cols = dimensions)
#' @param M Number of mixture components(clusters) (default: 5)
#' @param iter_max Maximum number of EM iterations (default: 1000)
#' @param tol Convergence tolerance for change in log-likelihood (default: 1e-1)
#' @return A list containing estimated parameters:
#'   \item{pi}{Mixture weights}
#'   \item{c}{Cluster centers}
#'   \item{r}{Cluster radii}
#'   \item{phi}{Mean direction of each cluster}
#'   \item{tau}{Angular spread of each cluster}
#'   \item{sigma_sq}{Variance of Gaussian noise}
#'   \item{log_likelihood}{Log-likelihood trajectory over iterations}
#'   \item{n_iter}{Number of iterations performed}
#'   \item{W}{Posterior cluster assignment probabilities}
#' @examples
#' # Generate example 2D data
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol = 2)
#' fit <- FG_EM(x, M = 3, iter_max = 100, tol = 1e-2)
#' @noRd
#' 
FG_EM <- function(x, M, iter_max = 500, tol = 1e-1 ){
  
  n <- nrow(x)
  d <- ncol(x)
  
  ################### initialization
  initialPars <- initial_kmeans(x, M)
  
  log_likelihood_vals <- rep(NA, iter_max) # use NA to distinguish unconverged iterations
  pi_hat <- initialPars$pi_hat
  c_hat <- initialPars$c_hat
  r_hat <- initialPars$r_hat
  sigma_sq <-initialPars$sigma_sq
  phi_hat <- initialPars$phi_hat
  tau_hat <- initialPars$tau_hat
  
  ############ EM iterations
  for (iter in 1: iter_max){
    # E-step
    e_step_rst <- E_step(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
    W <- e_step_rst$W
    y_expect_list <- e_step_rst$y_expect_list
    
    ## log-likelihood
    log_likelihood_vals[iter] <- log_likelihood(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
    
    
    # M-step
    m_step_rst <- M_step(x, W, y_expect_list, c_hat, r_hat, phi_hat, tau_hat)
    pi_hat <- m_step_rst$pi_hat
    c_hat <- m_step_rst$c_hat
    r_hat <- m_step_rst$r_hat
    phi_hat <- m_step_rst$phi_hat
    tau_hat <- m_step_rst$tau_hat
    sigma_sq <- m_step_rst$sigma_sq
    
    # convergence check
    if (iter > 1 && abs(log_likelihood_vals[iter] - log_likelihood_vals[iter - 1]) < tol) {
      message("Converged at iteration ", iter)
      break
    }
    
  }
  
  converged_loglik <- log_likelihood_vals[1:iter]
  return(list(
    pi = pi_hat, c = c_hat, r = r_hat, 
    phi = phi_hat, tau = tau_hat, sigma_sq = sigma_sq, 
    log_likelihood = converged_loglik,
    n_iter = iter,
    W = W))
}

#' Initialize parameters for FG-EM algorithm using k-means
#' 
#' Computes initial estimates of mixture weights, cluster centers, radii, angular directions, 
#' angular spread, and residual variance for the Fisher-Gaussian EM algorithm.
#'
#' @param x Numeric matrix of coordinates
#' @param M Number of clusters
#' @return A list containing initialized parameters:
#'   \item{pi_hat}{Initial mixture weights}
#'   \item{c_hat}{Initial cluster centers}
#'   \item{r_hat}{Initial radii}
#'   \item{sigma_sq}{Initial residual variance}
#'   \item{phi_hat}{Initial mean directions}
#'   \item{tau_hat}{Initial angular spreads}
#' @noRd
initial_kmeans <- function(x, M){
  n <- nrow(x)
  pi_hat <- rep(1/M, M)
  set.seed(12345)
  init_kmeans <- kmeans(x, centers = M, nstart = 1000, algorithm = 'Lloyd')
  c_hat <- init_kmeans$centers
  r_hat <- numeric(M)
  for (k in 1: M){
    cluster_points <- x[init_kmeans$cluster == k, , drop = FALSE]
    avg_distance <- mean(apply(cluster_points, 1, function(row) sqrt(sum((row - c_hat[k,])^2))))
    r_hat[k] <- avg_distance
  }

  # initialize sigma_sq based on the squared residual
  residuals <- numeric(n)
  for (i in 1:n) {
    cluster_center <- init_kmeans$centers[init_kmeans$cluster[i], ]
    residuals[i] <- sum((x[i, ] - cluster_center)^2)
  }
  sigma_sq <- mean(residuals)

  # recover y from x
  y <- matrix(0, nrow = nrow(x), ncol = ncol(x))
	for (i in 1:nrow(x)) {
		cluster_idx <- init_kmeans$cluster[i]
		y_i <- (x[i, ] - c_hat[cluster_idx, ]) / r_hat[cluster_idx]
		y[i, ] <- y_i / sqrt(sum(y_i^2))  # Normalize y_i
	}
                               
  # initialize phi_hat (mean direction)
  phi_hat <- matrix(0, nrow = M, ncol = ncol(x))
	for (k in 1:M) {
		cluster_points <- which(init_kmeans$cluster == k)
		if (length(cluster_points) > 0) {
		# Compute mean direction phi_k
			y_cluster <- y[cluster_points, , drop = FALSE]
			phi_hat[k, ] <- colSums(y_cluster)
			phi_hat[k, ] <- phi_hat[k, ] / sqrt(sum(phi_hat[k, ]^2))
		}
	}
  # initialize tau_hat based on the angular spread of points in the cluster
  tau_hat <- numeric(M)
	for (k in 1:M) {
		# Get data points in cluster k
		cluster_points <- x[init_kmeans$cluster == k, ]
    
		# Skip if the cluster has no points or only one point
		if (nrow(cluster_points) <= 1) {
			tau_hat[k] <- 1  # Default value for small clusters
			next
		}    
		phi_k <- phi_hat[k,]   
		# Compute angular distances for all points in the cluster
		dot_products <- rowSums(cluster_points * phi_k)  # Dot product between each point and phi_k
		norms <- sqrt(rowSums(cluster_points^2))  # Norm of each point
		angular_distances <- acos(pmin(pmax(dot_products / norms, -1), 1))  # Ensure acos input is in [-1, 1]
    
		# Estimate tau_k as the inverse of the average angular distance
		# Use median for robustness to outliers
		tau_hat[k] <- 1 / median(angular_distances)
	}
  return(list(pi_hat = pi_hat, c_hat = c_hat, r_hat = r_hat, sigma_sq = sigma_sq, phi_hat = phi_hat, tau_hat = tau_hat))
}

#' E-step of the Fisher-Gaussian EM algorithm
#'
#' Computes posterior probabilities (W) of cluster assignments and expected 
#' normalized directions (y_expect) for each cluster.
#'
#' @param x Numeric matrix of coordinates (n x d)
#' @param pi_hat Numeric vector of mixture weights
#' @param c_hat Numeric matrix of cluster centers (M x d)
#' @param r_hat Numeric vector of cluster radii
#' @param phi_hat Numeric matrix of mean directions (M x d)
#' @param tau_hat Numeric vector of angular spread parameters
#' @param sigma_sq Numeric scalar of Gaussian noise
#' @return A list with:
#'   \item{W}{Posterior probabilities of cluster membership (n x M)}
#'   \item{y_expect_list}{List of expected normalized vectors for each cluster}
#' @noRd

E_step <- function(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq){
  n <- nrow(x)
  M <- length(pi_hat)
  W <- matrix(0, n, M)
  y_expect_list <- vector('list', M)
  
  W <- compute_W(x, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq)
  W <- W/rowSums(W) # normalize 
  
  
  # Compute E(y_i|x_i, z_i=k, Theta)
  y_expect_list <- compute_y_expect(x, M, phi = phi_hat, tau = tau_hat, c= c_hat, r= r_hat, sigma_sq)
    
  return(list(W = W, y_expect_list = y_expect_list))
}


#' M-step of the Fisher-Gaussian EM algorithm
#'
#' Updates mixture weights, cluster centers, radii, mean directions, angular spreads, and Gaussian noise
#' based on the current E-step results.
#'
#' @param x Numeric matrix of coordinates (n x d)
#' @param W Posterior cluster probabilities (n x M)
#' @param y_expect_list List of expected normalized vectors for each cluster
#' @param c_hat Current cluster centers (M x d)
#' @param r_hat Current cluster radii
#' @param phi_hat Current mean directions (M x d)
#' @param tau_hat Current angular spreads
#' @return A list with updated parameters:
#'   \item{pi_hat}{Updated mixture weights}
#'   \item{c_hat}{Updated cluster centers}
#'   \item{r_hat}{Updated radii}
#'   \item{phi_hat}{Updated mean directions}
#'   \item{tau_hat}{Updated angular spreads}
#'   \item{sigma_sq}{Updated Gaussian noise}
#' @noRd

M_step <- function(x, W, y_expect_list, c_hat, r_hat, phi_hat, tau_hat){
  n <- nrow(x)
  d <- ncol(x)
  M <- ncol(W)
  
  # update mixture weight pi_k
  pi_hat <- colSums(W)/n
  
  # update c_k, r_k, phi_k and tau_k for each cluster
  
  
for (k in 1: M){
    y_expect_k <- y_expect_list[[k]]
	
    # update c_k --> checked, correct
    c_hat[k, ] <- colSums(W[, k] * (x - r_hat[k] * y_expect_k)) / sum(W[, k])
    
    # update r_k --- checked, correct
    temp <- numeric(n)
    for(i in 1:n){
      temp[i] <- W[i,k]*sum((x[i,] - c_hat[k,])*y_expect_k[i,])
    }
    r_hat[k] <- sum(temp)/sum(W[,k])
    
    
    # update phi_k --- checked, correct
    phi_hat[k, ] <- colSums(W[, k] * y_expect_k)
    phi_hat[k, ] <- phi_hat[k, ] / norm_vec(phi_hat[k, ])
	}
  
  ## opt4: based on the formula(2) from the movMF package 
  for (k in 1:M){
    #rbar <- norm_vec(phi_hat[k])/pi_hat[k]  # should not use pi_hat[k], but sum(W[, k]) instead
    rbar <- norm_vec(colSums(W[, k] * y_expect_list[[k]])) / sum(W[, k])  # also not used phi_hat[k], as it is normalized, should use unnormalized phi_hat here
    temp1 <- (rbar * d - rbar^3)/(1 - rbar^2)
    tau_hat[k] <- pmax(temp1, 1e-6)

  } 
  
  # update sigma_sq --- checked, correct
    temp <- matrix(0, n, M)
    for (i in 1:n) {
      for (k in 1:M) {
        y_expect_k <- y_expect_list[[k]]  # Get E[y_i | x_i, z_i = k] for cluster k
        a <- x[i,]-c_hat[k,]
        norm_a <- sum(a^2)
        temp[i, k] <- W[i,k] * (norm_a - 2 * r_hat[k] * sum(a * y_expect_k[i,]) + (r_hat[k])^2)
      }
    }
    sigma_sq <- sum(temp) / (n * d)
 
  
  return(list(pi_hat = pi_hat, c_hat = c_hat, r_hat = r_hat,
              phi_hat = phi_hat, tau_hat = tau_hat, sigma_sq = sigma_sq))
}

#' Fit FG-EM model and compute model selection criteria
#'
#' Runs the Fisher-Gaussian EM algorithm for a given number of clusters and computes
#' AIC and BIC for model selection.
#'
#' @param x Numeric matrix of coordinates (n x d)
#' @param M Integer. Number of mixture components
#' @param iter_max Maximum EM iterations (default = 500)
#' @param tol Convergence tolerance for EM algorithm (default = 1e-1)
#' @return A list containing:
#'   \item{model}{List returned by FG_EM, including fitted parameters and log-likelihoods}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{n_parameters}{Number of estimated parameters in the model}
#'   \item{converged_iter}{Number of EM iterations until convergence}
#' @noRd
FG_EM_with_criteria <- function(x, M, iter_max, tol) {
  fit <- FG_EM(x, M, iter_max, tol)
  
  # Compute number of parameters (p)
  n <- nrow(x)
  d <- ncol(x)
  p <- (M-1) + M*d + M + M*(d-1) + M + 1  # π_k (M-1), c_k (M*d), r_k (M), φ_k (M*(d-1) due to norm=1 constrain), τ_k (M), σ² (1)
  
  # Use the LAST (converged) log-likelihood value
  final_loglik <- tail(fit$log_likelihood, 1)
  
  # Compute BIC and AIC
  BIC <- -2 * final_loglik + p * log(n)
  AIC <- -2 * final_loglik + 2 * p
  
  # Return model with BIC
  return(list(
    model = fit,
    AIC = AIC,
    BIC = BIC,
    n_parameters = p,
    converged_iter = fit$n_iter
  ))
}

#' Select the best number of clusters via BIC/AIC
#'
#' Fits FG-EM models for multiple candidate numbers of clusters and selects the best model
#' based on BIC and AIC.
#'
#' @param x Numeric matrix of coordinates (n x d)
#' @param M_candidates Integer vector of candidate cluster numbers (default = 2:5)
#' @param iter_max Maximum EM iterations (default = 1000)
#' @param tol Convergence tolerance (default = 1e-1)
#' @return A list containing:
#'   \item{all_models}{List of models fitted for each candidate M}
#'   \item{criteria_table}{Data frame summarizing AIC, BIC, convergence, and iterations for each M}
#'   \item{best_M_AIC}{Best M according to AIC}
#'   \item{best_M_BIC}{Best M according to BIC}
#'   \item{best_model_AIC}{FG-EM model corresponding to best AIC}
#'   \item{best_model_BIC}{FG-EM model corresponding to best BIC}
#' @noRd

select_best_M <- function(x, M_candidates = 2:5, iter_max = 1000, tol = 1e-1) {
  # Initialize storage
  results <- list()
  criteria_df <- data.frame(
    M = M_candidates, 
    AIC = NA_real_,
    BIC = NA_real_,
    converged = NA,
    n_iter = NA
  )
  
  # Fit models for all M candidates
  for (i in seq_along(M_candidates)) {
    M <- M_candidates[i]
    message("\nFitting M = ", M)
    result <- FG_EM_with_criteria(x, M, iter_max, tol)
    
    # Store full results (including model, AIC, BIC, etc.)
    results[[paste0("M_", M)]] <- result
    
    # Update criteria table
    criteria_df$AIC[i] <- result$AIC
    criteria_df$BIC[i] <- result$BIC
    criteria_df$converged[i] <- (result$model$n_iter < iter_max)
    criteria_df$n_iter[i] <- result$model$n_iter
  }
  
  # Find best M by AIC/BIC (optional, since we keep all results)
  best_AIC_idx <- which.min(criteria_df$AIC)
  best_BIC_idx <- which.min(criteria_df$BIC)
  
  # Return ALL models + summary table
  return(list(
    all_models = results,           # Full results for every M
    criteria_table = criteria_df,    # Summary table
    best_M_AIC = M_candidates[best_AIC_idx],
    best_M_BIC = M_candidates[best_BIC_idx],
    best_model_AIC = results[[best_AIC_idx]]$model,  # Optional
    best_model_BIC = results[[best_BIC_idx]]$model   # Optional
  ))
}


## Given the FG-mixture results and a number N, this function generates N samples from estimated density
#' @import movMF
#' @import mvtnorm
#' @import matrixStats
#' @import Directional
#' @import MASS
#' @noRd
plot_density_FG_EM <- function(N, rst){
  cntr <- rst$c
  rd <- rst$r
  phi <- rst$phi
  tau <- rst$tau
  Pi <- rst$pi
  sigma.sq <- rst$sigma_sq
  D <- ncol(cntr)
  M <- length(tau)
  
  z <- matrix(data = NA,nrow = N,ncol = D)
  
  
  kr_label <- sapply(1:N,function(u){v <- stats::rmultinom(n = 1,size = 1,prob = Pi); return(which.max(v))})
  for(k in 1:M){
      idx <- which(kr_label == k)
      if(length(idx) > 0){
          yy <- Directional::rvmf(length(idx), phi[k, ],tau[k])
          if(length(idx) == 1){
            z[idx,] <- cntr[k, ] + rd[k] * yy
          }
          else{
            z[idx,] <- t(apply(yy, 1, function(u){return(cntr[k, ] + rd[k] * u)})) 
          }
      }
      
  }
  error <- MASS::mvrnorm(N, rep(0, D), sigma.sq * diag(D))
  w <- z + error
  return(list(z,w))
}
