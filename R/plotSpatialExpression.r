#' Plot spatial expression of one reference or simulated gene
#'
#' @param spaDesign A spaDesign object.
#' @param gene One gene name or an exact count-matrix row name.
#' @param source Either "simulated" (default) or "reference".
#' @param domain Optional model-domain name used to resolve a simulated row
#'   such as "Layer2-SNAP25". This selects a simulated gene version; it does
#'   not restrict the displayed spots. Not used for reference data.
#' @param transform Either "log1p" (natural log of count + 1) or "none".
#' @param point_size Positive numeric point size.
#' @param reverse_y Logical; reverse the y-axis for image-style coordinates.
#' @return A ggplot object. Its data contain barcode, x, y, count, and expression.
#' @details Counts and coordinates are aligned by barcode. Only the selected
#'   row is converted to a numeric vector, preserving sparse-matrix efficiency.
#'   Exact row-name matches take priority. Otherwise, known domain prefixes
#'   are matched literally. Ambiguous matches produce an error listing choices.
#'   Domain-specific simulated rows are never summed. Expression is displayed
#'   as counts or log1p(counts), without library-size normalization.
#' @importFrom rlang .data
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' plotSpatialExpression(sim_obj, "SNAP25")
#' plotSpatialExpression(sim_obj2, "SNAP25", domain = "Layer2")
#' plotSpatialExpression(DATA, "SNAP25", source = "reference")
#' }
plotSpatialExpression <- function(spaDesign, gene,
                                  source = c("simulated", "reference"),
                                  domain = NULL,
                                  transform = c("log1p", "none"),
                                  point_size = 1.5, reverse_y = TRUE) {
  source <- match.arg(source)
  transform <- match.arg(transform)
  fail <- function(...) stop(..., call. = FALSE)
  scalar_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if (!methods::is(spaDesign, "spaDesign")) fail("spaDesign must be a spaDesign object.")
  if (!scalar_string(gene)) fail("gene must be one nonempty character string.")
  if (!is.null(domain) && !scalar_string(domain)) fail("domain must be NULL or one nonempty character string.")
  if (source == "reference" && !is.null(domain)) fail("domain selects simulated model versions; omit it for reference data.")
  if (!is.numeric(point_size) || length(point_size) != 1L ||
      !is.finite(point_size) || point_size <= 0) fail("point_size must be positive and finite.")
  if (!is.logical(reverse_y) || length(reverse_y) != 1L || is.na(reverse_y)) fail("reverse_y must be TRUE or FALSE.")

  counts <- methods::slot(spaDesign, if (source == "simulated") "simCounts" else "refCounts")
  meta <- methods::slot(spaDesign, if (source == "simulated") "simcolData" else "refcolData")
  if (length(dim(counts)) != 2L || any(dim(counts) == 0L)) fail("No ", source, " counts available.")
  valid_names <- function(x) !is.null(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  if (!valid_names(rownames(counts)) || !valid_names(colnames(counts))) fail("Count rows and columns must have unique, nonempty names.")
  if (!valid_names(rownames(meta)) || !all(c("x", "y") %in% names(meta))) fail("Metadata require unique barcode row names and x/y columns.")
  idx <- match(colnames(counts), rownames(meta))
  if (anyNA(idx)) fail("Some count-matrix barcodes are missing from the spatial metadata.")
  meta <- meta[idx, , drop = FALSE]
  if (!is.numeric(meta$x) || !is.numeric(meta$y) ||
      any(!is.finite(meta$x)) || any(!is.finite(meta$y))) fail("Coordinates must be finite numeric values.")

  rows <- rownames(counts)
  selected <- gene
  if (!is.null(domain)) {
    # A domain argument resolves prefixed rows only; never guesses an exact
    # unprefixed row's model-domain identity.
    selected <- if (startsWith(gene, paste0(domain, "-")) && gene %in% rows) gene else paste0(domain, "-", gene)
  } else if (!gene %in% rows && source == "simulated") {
    domains <- unique(c(names(methods::slot(spaDesign, "paramsGP")),
                        as.character(meta$domain)))
    domains <- domains[!is.na(domains) & nzchar(domains)]
    hits <- intersect(paste0(domains, "-", gene), rows)
    if (length(hits) > 1L) fail("Multiple simulated versions match '", gene,
      "'. Select an exact row name or supply domain. Choices: ", paste(hits, collapse = ", "))
    if (length(hits) == 1L) selected <- hits
  }
  if (!selected %in% rows) fail("Gene row '", selected, "' not found in ", source, " counts. Check rownames of the count matrix.")
  values <- as.numeric(counts[selected, , drop = FALSE])
  if (any(!is.finite(values)) || any(values < 0)) fail("Selected gene counts must be finite and nonnegative.")
  plot_data <- data.frame(barcode = colnames(counts), x = meta$x, y = meta$y,
                          count = values,
                          expression = if (transform == "log1p") log1p(values) else values)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$expression)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_viridis_c(option = "H") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = paste0(selected, " (", source, ")"),
                  colour = if (transform == "log1p") "log1p(count)" else "Count", x = "x", y = "y") +
    ggplot2::theme_classic()
  if (reverse_y) p <- p + ggplot2::scale_y_reverse()
  p
}


