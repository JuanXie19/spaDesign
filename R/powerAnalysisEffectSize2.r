#' Evaluate clustering across effect sizes and sequencing depths
#'
#' Simulates counts using fitted-GP means and evaluates
#' clustering with evaluatePowerSeurat().
#'
#' @param spaDesign A spaDesign object containing reference counts, spot metadata,
#'   and fitted GP models required by prepare_genelevel_mean().
#' @param es_range Nonempty numeric vector of finite effect powers > 0.
#'   One preserves baseline domain means, zero equalizes positive domain means,
#'   and values above one amplify domain mean ratios.
#' @param seq_depth_range Nonempty numeric vector of finite depth factors > 0.
#' @param n_rep Positive integer number of replicates per scenario.
#' @param mean_matrix Optional fixed baseline mean matrix returned by
#'   prepare_genelevel_mean() for the same reference and fitted models.
#'
#' @details The baseline GP mean matrix is prepared once. Each distinct effect
#'   power is applied once before parallel simulation. 
#'   Seeds 1 through n_rep are reused across scenarios, preserving the supplied
#'   simulation convention. Reusing seeds does not implement nested read
#'   thinning. The caller's future plan is restored on exit. Evaluation uses
#'   future.seed = TRUE; for reproducible stochastic clustering as well as
#'   counts, call set.seed() before invoking this function.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{seq_depth}{Sequencing depth scaling factor}
#'     \item{effect_size}{Effect size scaling factor}
#'     \item{total_counts}{Total counts in simulated dataset}
#'     \item{NMI}{Normalized Mutual Information between predicted and true labels}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' results <- powerAnalysisEffectSize(
#'   my_spadesign, es_range = c(0.5, 1, 3),
#'   seq_depth_range = c(0.5, 1, 3, 5), n_rep = 5)
#' }
powerAnalysisEffectSize2 <- function(spaDesign, es_range, seq_depth_range,
                                   n_rep, mean_matrix = NULL) {
  if (!inherits(spaDesign, 'spaDesign')) {
    stop('spaDesign must be a spaDesign object.')
  }
  if (!is.numeric(es_range) || !length(es_range) ||
      any(!is.finite(es_range)) || any(es_range < 0)) {
    stop('es_range must contain finite numeric effect powers >= 0.')
  }
  if (!is.numeric(seq_depth_range) || !length(seq_depth_range) ||
      any(!is.finite(seq_depth_range)) || any(seq_depth_range <= 0)) {
    stop('seq_depth_range must contain finite numeric factors > 0.')
  }
  check_integer <- function(x, label) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x > .Machine$integer.max || x != floor(x)) {
      stop(label, ' must be a positive integer within the R integer range.')
    }
  }
  check_integer(n_rep, 'n_rep')

  # Validate reference alignment even when a precomputed matrix is supplied.
  counts <- refCounts(spaDesign)
  meta <- as.data.frame(refcolData(spaDesign))
  .genelevel_reference_check(counts, meta)
  # prepare the mean_matrix
  if (is.null(mean_matrix)) mean_matrix <- prepare_genelevel_mean(spaDesign)
  .genelevel_mean_check(mean_matrix)
  if (!identical(colnames(mean_matrix), colnames(counts)) ||
      any(!rownames(mean_matrix) %in% rownames(counts))) {
    stop('mean_matrix must use reference genes and the exact reference spot order.')
  }

  effects <- unique(es_range)
  # adjust the mean matrix
  adjusted_means <- lapply(effects, function(es) {
    adjust_genelevel_effect(mean_matrix, meta$domain, effect_power = es)
  })
  param_grid <- expand.grid(
    seq_depth = seq_depth_range, effect_size = es_range,
    SEED = seq_len(n_rep), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  effect_index <- match(param_grid$effect_size, effects)

  process_row <- function(i) {
    depth <- param_grid$seq_depth[[i]]
    effect <- param_grid$effect_size[[i]]
    seed <- param_grid$SEED[[i]]
    tryCatch({
      # Effect is already applied: factor 1 prevents applying it twice.
      simulated <- simulation_effectSize2(
        spaDesign = spaDesign,
        seq_depth_factor = depth,
        effect_size_factor = 1,
        SEED = seed,
        mean_matrix = adjusted_means[[effect_index[[i]]]]
      )
      evaluated <- evaluatePowerSeurat(simulated)
      nmi <- evaluated@NMI
      if (!is.numeric(nmi) || length(nmi) != 1L || !is.finite(nmi)) {
        stop('evaluatePowerSeurat() must return one finite numeric NMI.')
      }
      data.frame(seq_depth = depth, effect_size = effect, SEED = seed,
                 total_counts = sum(simulated@simCounts), NMI = nmi)
    }, error = function(e) {
      stop(sprintf('Scenario failed (depth=%s, effect=%s, SEED=%s): %s',
                   depth, effect, seed, conditionMessage(e)), call. = FALSE)
    })
  }
  results <- lapply(
    seq_len(nrow(param_grid)), process_row)
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  results
}
