#' Complementary Log-Log Function
#'
#' \code{CLL} returns the complementary log-log function
#' \code{ICLL} returns the inverse complementary log-log function
#'
#' The complementary log-log function and its inverse are defined respectively by:
#'
#' $$\text{CLL}(x) = \log(-\log(1-x)) \quad \quad \quad \text{ICLL}(x) = 1-\exp(-\exp(x)).$$
#'
#' The complementary log-log function transforms a probability value to an extended real number and the inverse complementary log-log function
#' transforms an extended real number back to a probability value.  The functions are vectorised so that a vector of inputs is transformed to a
#' vector of outputs each transformed with the relevant function.  These functions also allow complex inputs and outputs.
#'
#' @usage \code{CLL/ICLL}
#' @param x An input value (numeric/complex scalar or vector)
#' @return The complementary log-log transformation of the input

CLL  <- function(x) { Log(-Log(1-x)) }

ICLL <- function(x) { 1-exp(-exp(x)) }

Log <- function(x, base = exp(1)) {
  LOG <- base::log(as.complex(x), base = base)
  if (all(Im(LOG) == 0)) { LOG <- Re(LOG) }
  LOG }
