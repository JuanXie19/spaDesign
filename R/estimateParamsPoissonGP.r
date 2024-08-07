#' Estimate Parameters for Poisson Gaussian Process Model
#' @param spadesign A \code{spaDesign} object
#' @param n_iter Number of MCMC iterations
#' @param n_burn Number of MCMC iterations to exclude as burn-in period
#' @param cov_type Integer specifying the covariance kernel (0 for Gaussian, 1 for Exponential)
#' @param save Logical value to save the posterior samples for the parameters (default = FALSE)
#' @return Updated \code{spaDesign} object with parameter estimates and optionally posterior samples
#' @import igraph
#' @import dplyr
#' @export
#' @examples
#' ## Example usage
#' toyDATA <- createSpaDesignObject(count_matrix = toyData$toyCount, loc = toyData$loc)
#' toyDATA <- estimateParamsPoissonGP(toyDATA, n_iter = 1000, n_burn = 200, cov_type = 1, save = FALSE)
#'

estimateParamsPoissonGP <- function(spadesign, n_iter = 50000, n_burn = 30000, cov_type = 1, save = FALSE){

    loc_file <- refcolData(spadesign)[, c('x','y','domain')]
    count_matrix <- refCounts(spadesign)
    
    ## scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[,c('x','y')]),xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain

    DOMAIN <- sort(unique(loc_file$domain))
    topGenes <- topGenes(spadesign)

    if (is.null(topGenes) || length(topGenes) == 0) {
        stop("topGenes is NULL or empty. Check the spaDesign object.")
    }
    
    RST <- list()
    POST <- list()
    for (i in seq_along(topGenes)){
        d <- names(topGenes)[i] # the domain name
        GENES <- rownames(topGenes[[i]])  # the informative genes for that domain
		
		if (length(GENES) == 0) {
          next
        }
        
        idx <- which(coords_norm$domain == d) # fit the model for within domain expression
        counts.genes <- count_matrix[rownames(count_matrix) %in% GENES, , drop = FALSE]
        counts.sub <- counts.genes[, idx, drop = FALSE]   # informative genes within-domain gene expression
        coords_norm_sub <- coords_norm[idx, ]  # domain coordinates
        
        GPfit <- PoissonGP_multi(
            coords = coords_norm_sub[, c('x', 'y')],
            counts = counts.sub,
            n_iter = n_iter,
            n_burn = n_burn,
            cov_type = cov_type
        )
		
		POST[[i]] <- GPfit
        
        PAR <- list()
        ## calculate posterior mean as parameter estimates
        for (g in seq_along(GPfit)){
        
            alpha.est <- mean(GPfit[[g]]$alpha)
            rho.est <- mean(GPfit[[g]]$rho)
            mu.est <- mean(GPfit[[g]]$mu[,1])
            s.est <- colMeans(GPfit[[g]]$s)
            temp <- list(alpha.est = alpha.est, rho.est = rho.est, mu.est = mu.est, s.est = s.est)
            PAR[[g]] <- temp
        }
        names(PAR) <- names(GPfit)
        RST[[i]] <- PAR
    }
    names(RST) <- names(topGenes)
    names(POST) <- names(topGenes)
    if(save == TRUE){
        spadesign@GP.post <- POST
        spadesign@paramsPoissonGP <- RST
    }else{
        rm(POST)
        spadesign@paramsPoissonGP <- RST
    }
    return(spadesign)    
}



#' Fit Poisson Gaussian Process Model for Multiple Domain-Informatives Genes
#' @param coords Matrix of 2-dimensional locations
#' @param counts Matrix of gene expression counts (rows = genes, columns = locations)
#' @param n_iter Number of MCMC iterations
#' @param n_burn Number of MCMC iterations to exclude as burn-in period
#' @param cov_type Integer specifying the covariance kernel (0 for Gaussian, 1 for Exponential)
#' @return A list of results for each gene including MCMC samples and computation time
#' @import igraph
#' @import nimble
#' @import parallel
#' @export


PoissonGP_multi <- function(coords, counts, n_iter, n_burn, cov_type ) {
  
  library(parallel)
 
  # check input
  if (!is.numeric(n_iter) || n_iter <= 0 || n_iter != as.integer(n_iter)) {
    stop('`n_iter` must be a positive integer.')
  }
  if (!is.numeric(n_burn) || n_burn <= 0 || n_burn != as.integer(n_burn)) {
    stop('`n_burn` must be a positive integer.')
  }
  if (cov_type != 0 && cov_type != 1) {
    stop('`cov_type` must be either 0 (squared exponential) or 1 (absolute exponential).')
  }

  ## prepare the test data
  test_data <- lapply(1:nrow(counts), function(i) counts[i, ])
  names(test_data) <- rownames(counts)
  
  
  code <- nimbleCode({
    mu0 ~ dnorm(0, sd = 100)
    alpha ~ dnorm(0, sd = 0.1)
    rho ~ dinvgamma(shape = 5, scale = 5)
    
    mu[1:N] <- ones[1:N] * mu0
    cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], alpha, rho, cov_type)
    s[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
    
    #likelihood
    for (i in 1:N) {
      lambda[i] <- exp(s[i])
      y[i] ~ dpois(lambda[i])                        
    }
  })
  # The key consideration is to ensure that all NIMBLE execution, including model building, is conducted inside the parallelized code
  # https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html
  run_MCMC_allcode <- function(idx, coords, test_data, code, cov_type, n_iter, n_burn) {
    
    library(nimble)
    
    tic <- proc.time()[[3]]
    dists <- as.matrix(dist(cbind(coords[ ,'x'], coords[ ,'y'])))
    n <- nrow(coords)
    
    expcov <- nimbleFunction(     
      run = function(dists = double(2), alpha = double(0), rho = double(0),
                     cov_type = integer()) {
        returnType(double(2))
        n <- dim(dists)[1]
        result <- matrix(nrow = n, ncol = n, init = FALSE)
        alpha2 <- alpha * alpha
        for(i in 1:n){
          for(j in 1:n){
            if(cov_type == 0) {
              result[i, j] <- alpha2 * exp(-(dists[i, j]^2) / (2 * rho^2))
            } else if(cov_type == 1){
              result[i, j] <- alpha2 * exp(-dists[i, j] / rho)
            } else{
              stop('The current version only supports exponential or  Gaussian kernel. Please choose from one of them.')
            }
          }                            
        }
        for(i in 1:n)
          result[i,i] <- result[i,i] + 1e-6
        return(result)
      })
    
    
    assign('expcov', expcov, envir = .GlobalEnv)
    cExpcov <- compileNimble(expcov)
    
    
    constants <- list(N = n, dists = dists, ones = rep(1, n), cov_type = cov_type)
    inits <- list(alpha = 1, rho = 0.2, mu0 = 0)
    
    set.seed(1)
    
    ## setup initial spatially-correlated latent process values
    inits$cov <- cExpcov(dists,inits$alpha, inits$rho, cov_type)
    inits$s <-  t(chol(inits$cov)) %*% rnorm(n)
    inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix
    
    
    model <- nimbleModel(code, constants = constants, data = list(y = test_data[[idx]]), inits = inits)
    cModel <- compileNimble(model)
    conf <- configureMCMC(model, monitors = c('alpha', 'rho', 'mu', 's'))
    MCMC <- buildMCMC(conf)
    cMCMC <- compileNimble(MCMC, project = cModel)
    results <- runMCMC(cMCMC, niter = n_iter, nburnin = n_burn, setSeed = 12345)
    
    samples_alpha <- results[,grep('alpha', colnames(results))]
    samples_mu <- results[,grep('mu', colnames(results))]
    samples_rho <- results[,grep('rho', colnames(results))]
    samples_s <- results[,grep('s', colnames(results))]
    
    # create output object    
    out <- list(n_iter = n_iter, n_burn = n_burn, cov_type = cov_type, alpha = samples_alpha, mu = samples_mu, rho = samples_rho, s = samples_s)
    toc <- proc.time()[[3]]
    out$comput_time <- unname(toc - tic)
    return(out)
    }
    
	 available_cores <- parallel::detectCores()
     num_batches <- ceiling(nrow(counts) / available_cores)

     message('Running MCMC for ', nrow(counts), ' genes using ', available_cores, ' cores in ', num_batches, ' batches')

     results <- vector("list", nrow(counts))
     for (batch in seq_len(num_batches)) {
        start_idx <- (batch - 1) * available_cores + 1
        end_idx <- min(batch * available_cores, nrow(counts))
    
        this_cluster <- parallel::makeCluster(min(available_cores, end_idx - start_idx + 1))
        batch_results <- tryCatch({
          parallel::parLapply(
          cl = this_cluster,
          X = start_idx:end_idx,
          fun = run_MCMC_allcode,
          test_data = test_data,
          code = code,
          coords = coords,
          n_iter = n_iter,
          n_burn = n_burn,
          cov_type = cov_type
         )
        }, error = function(e) {
        parallel::stopCluster(this_cluster)
        stop(e)
       })
      parallel::stopCluster(this_cluster)
    
      results[start_idx:end_idx] <- batch_results
     }

    names(results) <- names(test_data)
    return(results)

}
