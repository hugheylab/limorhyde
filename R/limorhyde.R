#' Convert a periodic time variable into components usable in linear models
#'
#' `limorhyde` decomposes a periodic time variable into multiple components
#' based on either the first harmonic of a Fourier series or on a periodic
#' smoothing spline.
#'
#' @param df Data frame containing one row for each sample.
#' @param timeColname Character string indicating the column in `df` that
#' contains the time at which each sample was acquired.
#' @param period Number corresponding to the period to use for the
#' decomposition (in units of the values in `timeColname`).
#' @param sinusoid If `TRUE`, the decomposition is based on sine and cosine
#' curves. If `FALSE`, the decomposition is based on a periodic smoothing
#' spline.
#' @param nKnots Number of knots for the periodic spline. Only used if
#' `sinusoid` is `FALSE`.
#'
#' @return A data frame with a row for each sample and a column for each
#' component of the time decomposition.
#'
#' @example R/limorhyde_example.R
#'
#' @export
limorhyde = function(df, timeColname, period = 24, sinusoid = TRUE,
                     nKnots = 3) {
  if (!(timeColname %in% colnames(df))) {
    stop('timeColname must be a named column in df.')}

  if (sinusoid) {
    d = data.frame(tcos = cos(df[[timeColname]] / period * 2 * pi),
                   tsin = sin(df[[timeColname]] / period * 2 * pi))
    colnames(d) = paste0(timeColname, c('_cos', '_sin'))
  } else {
    knots = seq(0, period - period / nKnots, length = nKnots)
    d = bigsplines::ssBasis(df[[timeColname]] %% period, knots = knots,
                            xmin = 0, xmax = period, periodic = TRUE)$X
    d = as.data.frame(d)
    colnames(d) = paste0(timeColname, '_knot', 1:nKnots)}

  return(d)}
