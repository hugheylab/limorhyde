#' Convert a periodic time variable into components usable in linear models
#'
#' `limorhyde` decomposes a periodic time variable into multiple components
#' based on either the first harmonic of a Fourier series or on a periodic
#' smoothing spline.
#'
#' @param time Numeric vector of times, e.g., at which samples were acquired.
#' @param colnamePrefix Character string with which to prefix the column names
#' of the basis.
#' @param period Number corresponding to the period to use for the
#' decomposition (in the same units as `time`).
#' @param sinusoid If `TRUE`, the decomposition is based on cosinor, i.e.,
#' cosine and sine. If `FALSE`, the decomposition is based on a periodic
#' smoothing spline from the `pbs` package.
#' @param nKnots Number of internal knots for the periodic spline. Only used if
#' `sinusoid` is `FALSE`.
#' @param intercept If `TRUE`, a column of ones will be included in the basis.
#'
#' @return A matrix with a row for each sample and a column for each
#' component of the time decomposition.
#'
#' @example R/limorhyde_example.R
#'
#' @export
limorhyde = function(time, colnamePrefix = NULL, period = 24, sinusoid = TRUE,
                     nKnots = 3, intercept = FALSE) {
  if (sinusoid) {
    b = getCosinorBasis(time, period, intercept)
  } else {
    b = getSplineBasis(time, period, nKnots, intercept)}
  colnames(b) = paste0(colnamePrefix, colnames(b))
  return(b)}


#' Basis matrix for cosinor
#'
#' Generate basis matrix for cosinor regression.
#'
#' @param x Values of the predictor variable.
#' @param period Period for the predictor variable.
#' @param intercept If `TRUE`, a column of ones will be included in the basis.
#'
#' @return A matrix with a row for each value of `x` and a column for each
#' component of the decomposition.
#'
#' @examples
#' b = getCosinorBasis(seq(0, 20, 4), period = 24, intercept = FALSE)
#'
#' @export
getCosinorBasis = function(x, period, intercept) {
  b = cbind(cos(x / period * 2 * pi),
            sin(x / period * 2 * pi))
  colnames(b) = c('cos', 'sin')
  b = addIntercept(b, intercept)
  return(b)}


#' Basis matrix for periodic splines
#'
#' Generate basis matrix for a periodic B-spline using \link[pbs]{pbs}.
#'
#' @param x Values of the predictor variable.
#' @param period Period for the predictor variable.
#' @param nKnots Number of internal knots.
#' @param intercept If `TRUE`, a column of ones will be included in the basis.
#'
#' @return A matrix with a row for each value of `x` and a column for each
#' component of the decomposition.
#'
#' @examples
#' b = getSplineBasis(seq(0, 20, 4), period = 24, nKnots = 3, intercept = FALSE)
#'
#' @export
getSplineBasis = function(x, period, nKnots, intercept) {
  knots = seq(0, period - period / nKnots, length = nKnots)
  b = pbs::pbs(x %% period, knots = knots,
               Boundary.knots = c(0, period))[, , drop = FALSE]
  colnames(b) = paste0('knot', 1:nKnots)
  b = addIntercept(b, intercept)
  return(b)}


addIntercept = function(b, intercept) {
  if (intercept) {
    b = cbind(1, b)
    colnames(b)[1] = 'intercept'}
  return(b)}
