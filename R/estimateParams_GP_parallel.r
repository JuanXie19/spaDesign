#' Estimate the parameters for Poisson Gaussian process model
#' @param spadesign A \code{spaDesign} object
#' @param nIter Number of MCMC iterations
#' @param nBurnin Number of MCMC iterations to exclude as burn-in period
#' @param covType Integer, specifying the covariance kernel, either exponential (default, \code{covType} = 1) or Gaussian kernel (\code{covType} = 0, a.k.a squared exponential, radial basis function kernel)
#' @param save Logical value, to save the posterior samples for the parameters (\code{save} = TRUE) or not (default, \code{save} = FALSE)
#' @return Returns a \code{spaDesign} object with updated parameters
#' @import pbapply
#' @import parallel
#' @import igraph
#' @export
#' @examples
#' # Example usage:
#' # toyDATA <- createspadesignObject(count_matrix = toyData$toyCount, loc = toyData$loc)
#' # result <- estimateParams_GP(spadesign = toyDATA, nIter = 1000, nBurnin = 200, covType = 1, save = TRUE)
estimateParams_GP_parallel <- function(spadesign, nIter, nBurnin, covType, save = FALSE){
    library(pbapply)
    library(parallel)
    
    loc_file <- refcolData(spadesign)[, c('x', 'y', 'domain')]
    count_matrix <- refCounts(spadesign)
    
    # Scale the coordinates to [0,1] range
    coords_norm <- igraph::norm_coords(as.matrix(loc_file[, c('x', 'y')]), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    coords_norm <- as.data.frame(coords_norm)
    coords_norm$domain <- loc_file$domain
    
    DOMAIN <- sort(unique(loc_file$domain))
    topGenes <- topGenes(spadesign)
    
    RST <- list()
    POST <- list()
    
    # Parallel setup
    n_cores <- detectCores()
    cl <- makeCluster(n_cores)-2
    clusterExport(cl, varlist = c("coords_norm", "count_matrix", "topGenes", "PoissonGP_multi", "nIter", "nBurnin", "covType"))
    
    # Progress bar
    pb <- txtProgressBar(min = 0, max = length(DOMAIN) * length(topGenes), style = 3)
    
    process_domain <- function(d) {
        domain_data <- coords_norm[coords_norm$domain == d, ]
        counts_data <- count_matrix[, rownames(domain_data)]
        
        genes <- names(topGenes)[which(sapply(topGenes, function(g) any(g %in% rownames(counts_data))))] 
        domain_results <- list()
        
        for (g in genes) {
            GENES <- rownames(topGenes[[g]])
            counts_genes <- counts_data[which(rownames(counts_data) %in% GENES), ]
            coords_domain <- domain_data
            
            GPfit <- PoissonGP_multi(coords = coords_domain[, c('x', 'y')], counts = counts_genes, n_iter = nIter, n_burn = nBurnin, cov_type = covType)
            POST[[g]] <- GPfit
            
            PAR <- lapply(GPfit, function(fit) {
                list(
                    alpha.est = mean(fit$alpha),
                    rho.est = mean(fit$rho),
                    mu.est = mean(fit$mu[, 1]),
                    s.est = colMeans(fit$s)
                )
            })
            names(PAR) <- names(GPfit)
            domain_results[[g]] <- PAR
            
            # Update progress bar
            setTxtProgressBar(pb, which(DOMAIN == d) * length(topGenes) + which(genes == g))
        }
        
        return(domain_results)
    }
    
    # Process domains in parallel
    results <- pblapply(DOMAIN, process_domain, cl = cl)
    names(results) <- DOMAIN
    
    # Update spaDesign object
    if (save) {
        spadesign@GP.post <- POST
        spadesign@paramGP <- results
    } else {
        spadesign@paramGP <- results
    }
    
    close(pb)
    stopCluster(cl)
    
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
  
    message('Running MCMC for ', nrow(counts), ' genes')
    this_cluster <- makeCluster(nrow(counts))
    chain_output <- parLapply(cl = this_cluster, X = 1:nrow(counts), fun = run_MCMC_allcode, test_data = test_data, code = code,
                            coords = coords, n_iter = n_iter, n_burn = n_burn, cov_type = cov_type)
    names(chain_output) <- names(test_data)
    return(chain_output)  
}
