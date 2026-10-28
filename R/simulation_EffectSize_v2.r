#' Prepare fixed gene-level spatial means
#'
#' Predicts selected gene-domain blocks using existing fitted GP models. Other
#' blocks use the full reference mean for that gene in the corresponding domain.
#'
#' @param spaDesign A spaDesign object with refCounts, refcolData and paramsGP.
#' @return A numeric gene-by-spot matrix of Poisson means.
#' @details Coordinates are normalized globally to [0,1] on each axis. This must match the coordinate transformation
#'   used to fit paramsGP. Predictions are transformed using
#'   round(exp(prediction + tau.sq / 2)) - 1, preserving the supplied convention.
#' @export
prepare_genelevel_mean <- function(spaDesign) {
  counts <- refCounts(spaDesign)
  meta <- as.data.frame(refcolData(spaDesign))
  models <- paramsGP(spaDesign)
  .genelevel_reference_check(counts, meta)
  domains <- unique(as.character(meta$domain))
  if (!is.list(models) || !length(models) ||
      !.genelevel_names_ok(names(models), length(models)) ||
      any(!names(models) %in% domains)) {
    stop('paramsGP must be a nonempty list with unique reference-domain names.')
  }
  for (d in names(models)) {
    z <- models[[d]]
    if (!is.list(z) || (length(z) &&
        !.genelevel_names_ok(names(z), length(z)))) {
      stop('Each paramsGP domain must contain a list of uniquely named gene models.')
    }
  }
  genes <- unique(unlist(lapply(models, names), use.names = FALSE))
  if (!length(genes) || any(!genes %in% rownames(counts))) {
    stop('Selected genes must be nonempty and present in refCounts.')
  }
  xy <- as.matrix(meta[, c('x', 'y'), drop = FALSE])
  if (!is.numeric(xy) || any(!is.finite(xy)) ||
      any(apply(xy, 2L, function(z) diff(range(z))) <= 0)) {
    stop('Coordinates must be finite and vary on both axes.')
  }
  xy <- igraph::norm_coords(xy, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  mu <- matrix(NA_real_, length(genes), ncol(counts),
               dimnames = list(genes, colnames(counts)))
  for (d in domains) {
    idx <- which(as.character(meta$domain) == d)
    fits <- models[[d]]
    for (g in genes) {
      if (g %in% names(fits)) {
        fit <- fits[[g]]
        tau2 <- unname(fit$Theta['tau.sq'])
        if (length(tau2) != 1L || !is.finite(tau2) || tau2 < 0) {
          stop('Invalid tau.sq for gene ', g, ' in domain ', d, '.')
        }
        pred <- BRISC::BRISC_prediction(fit, xy[idx, , drop = FALSE])
        value <- round(exp(as.numeric(pred$prediction) + tau2 / 2)) - 1
        if (length(value) != length(idx) || any(!is.finite(value)) ||
            any(value < 0)) {
          stop('Invalid GP mean for gene ', g, ' in domain ', d, '.')
        }
        mu[g, idx] <- value
      } else {
        value <- as.numeric(counts[g, idx, drop = TRUE])
        if (any(!is.finite(value)) || any(value < 0)) {
          stop('Invalid reference counts for gene ', g, ' in domain ', d, '.')
        }
        mu[g, idx] <- mean(value)
      }
    }
  }
  .genelevel_mean_check(mu)
  mu
}

#' Adjust between-domain expression effects
#'
#' @param mean_matrix Numeric gene-by-spot mean matrix with unique names.
#' @param domain Domain labels in exactly the same order as matrix columns.
#'   If named, names must match the matrix columns exactly.
#' @param effect_power Finite numeric scalar > 0. One leaves means unchanged;
#'   values above one amplify domain mean ratios.
#' @return An adjusted mean matrix preserving each gene's total expected count.
#' @details For positive domain means m[d], the adjusted domain mean ratio for
#'   domains d and e is (m[d]/m[e])^effect_power. Within each domain all spots
#'   receive the same multiplier. A nonzero gene with any zero domain mean
#'   causes an error for effect_power other than one.
#'   All-zero genes remain zero. 
#' @export
adjust_genelevel_effect <- function(mean_matrix, domain, effect_power = 1) {
  .genelevel_mean_check(mean_matrix)
  .genelevel_scalar_check(effect_power, 'effect_power', allow_zero = TRUE)
  if (length(domain) != ncol(mean_matrix) || anyNA(domain) ||
      any(!nzchar(as.character(domain)))) {
    stop('domain must supply one nonmissing, nonempty label per spot.')
  }
  if (!is.null(names(domain)) &&
      !identical(names(domain), colnames(mean_matrix))) {
    stop('Names of domain must match mean_matrix columns in exactly the same order.')
  }
  if (effect_power == 1) return(mean_matrix)
  domain <- as.character(domain)
  groups <- lapply(unique(domain), function(d) which(domain == d))
  adjusted <- mean_matrix
  for (i in seq_len(nrow(mean_matrix))) {
    baseline <- mean_matrix[i, ]
    A <- sum(baseline)
    if (!is.finite(A)) stop('Nonfinite total mean for gene ', rownames(mean_matrix)[i], '.')
    if (A == 0) next
    dm <- vapply(groups, function(idx) mean(baseline[idx]), numeric(1))
    if (any(dm <= 0)) {
      stop('Gene ', rownames(mean_matrix)[i],
           ' has a zero domain mean; resolve zero means before adjusting effect.')
    }
    logw <- (effect_power - 1) * log(dm)
    if (any(!is.finite(logw))) stop('Effect scaling exceeds numeric range.')
    w <- exp(logw - max(logw))
    if (any(w == 0)) stop('Effect scaling underflow; reduce effect_power.')
    tmp <- baseline
    for (k in seq_along(groups)) {
      idx <- groups[[k]]
      tmp[idx] <- baseline[idx] * w[k]
    }
    B <- sum(tmp)
    if (!is.finite(B) || B <= 0) stop('Invalid adjusted total mean.')
    value <- (tmp / B) * A
    if (any(!is.finite(value)) || any(value[baseline > 0] <= 0)) {
      stop('Effect scaling exceeds numeric range.')
    }
    adjusted[i, ] <- value
  }
  adjusted
}

#' Simulate gene-level spatial counts with effect and depth scaling
#'
#' @param spaDesign A spaDesign object containing reference counts and metadata.
#' @param seq_depth_factor Finite numeric scalar > 0; multiplies Poisson means.
#' @param effect_size_factor Finite numeric scalar > 0; passed to effect_power
#'   in adjust_genelevel_effect. 
#' @param SEED Nonnegative integer seed.
#' @param mean_matrix Optional fixed means returned by prepare_genelevel_mean
#'   for this reference and these fitted models. Columns must match reference counts exactly; rows must be
#'   reference genes. The caller must regenerate this matrix if counts,
#'   coordinates, domains or fitted models change.
#' @return The spaDesign object with simCounts and simcolData updated.
#' @details Every call samples independent Poisson entries conditional on the
#'   adjusted means, setting the seed once per dataset. For a given gene, the
#'   expected total is seq_depth_factor times its baseline total.
#' @examples
#' \dontrun{
#' mu <- prepare_genelevel_mean(my_spadesign)
#' simulated <- simulation_effectSize2(
#'   my_spadesign, seq_depth_factor = 0.5,
#'   effect_size_factor = 3, SEED = 1, mean_matrix = mu
#' )
#' simulated@simCounts
#' }
#' @export
simulation_effectSize2 <- function(spaDesign, seq_depth_factor = 1,
                                 effect_size_factor = 1, SEED = 123,
                                 mean_matrix = NULL) {
  .genelevel_scalar_check(seq_depth_factor, 'seq_depth_factor')
  .genelevel_scalar_check(effect_size_factor, 'effect_size_factor', TRUE)
  if (!is.null(SEED) && (!is.numeric(SEED) || length(SEED) != 1L ||
      !is.finite(SEED) || SEED < 0 || SEED > .Machine$integer.max ||
      SEED != floor(SEED))) {
    stop('SEED must be NULL or a nonnegative integer within the R integer range.')
  }
  counts <- refCounts(spaDesign)
  metadata <- refcolData(spaDesign)
  meta <- as.data.frame(metadata)
  .genelevel_reference_check(counts, meta)
  if (is.null(mean_matrix)) mean_matrix <- prepare_genelevel_mean(spaDesign)
  .genelevel_mean_check(mean_matrix)
  if (!identical(colnames(mean_matrix), colnames(counts)) ||
      any(!rownames(mean_matrix) %in% rownames(counts))) {
    stop('mean_matrix must use reference genes and the exact reference spot order.')
  }
  mu <- adjust_genelevel_effect(mean_matrix, meta$domain, effect_size_factor)
  lambda <- seq_depth_factor * as.vector(mu)
  if (any(!is.finite(lambda))) stop('Depth scaling produces nonfinite Poisson means.')
  if (!is.null(SEED)) set.seed(SEED)
  values <- stats::rpois(length(lambda), lambda = lambda)
  if (any(!is.finite(values))) stop('Poisson sampling produced nonfinite counts.')
  spaDesign@simCounts <- matrix(values, nrow(mu), ncol(mu), dimnames = dimnames(mu))
  spaDesign@simcolData <- metadata
  spaDesign
}

# Internal validation helpers; no exports.
.genelevel_names_ok <- function(x, n) {
  !is.null(x) && length(x) == n && !anyNA(x) &&
    all(nzchar(x)) && !anyDuplicated(x)
}

.genelevel_scalar_check <- function(x, label, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (allow_zero) x < 0 else x <= 0)) {
    stop(label, if (allow_zero) ' must be a finite numeric scalar >= 0.'
         else ' must be a finite numeric scalar > 0.')
  }
  invisible(NULL)
}

.genelevel_mean_check <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || any(dim(x) == 0L) ||
      any(!is.finite(x)) || any(x < 0) ||
      !.genelevel_names_ok(rownames(x), nrow(x)) ||
      !.genelevel_names_ok(colnames(x), ncol(x))) {
    stop('Means must be a finite nonnegative numeric matrix with unique gene and spot names.')
  }
  invisible(NULL)
}

.genelevel_reference_check <- function(counts, meta) {
  if (length(dim(counts)) != 2L || any(dim(counts) == 0L) ||
      !.genelevel_names_ok(rownames(counts), nrow(counts)) ||
      !.genelevel_names_ok(colnames(counts), ncol(counts))) {
    stop('refCounts must have nonzero dimensions and unique gene and spot names.')
  }
  if (!all(c('x', 'y', 'domain') %in% names(meta)) ||
      !identical(colnames(counts), rownames(meta)) ||
      anyNA(meta$domain) || any(!nzchar(as.character(meta$domain)))) {
    stop('refcolData needs x, y and nonempty domain labels, with rows in exact refCounts column order.')
  }
  invisible(NULL)
}
