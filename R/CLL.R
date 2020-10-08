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

Log <- function(x, base = NULL) {

  #Check input
  if (!is.vector(x))                            { stop('Error: input x should be a vector') }
  if ((!is.numeric(x))&(!is.complex(x)))        { stop('Error: input x should be a numeric/complex vector') }
  if (!missing(base)) {
  if (!is.vector(base))                         { stop('Error: input base should be a vector')  }
  if ((!is.numeric(base))&(!is.complex(base)))  { stop('Error: input base should be a numeric/complex vector') }
  if (length(base) != 1) {
  if (length(base) != length(x))                { stop('Error: input base must either be a scalar or a vector with the same length as x') } } }

  #Compute natural logarithm
  LOG <- complex(real = log(abs(x)), imaginary = atan2(Im(x), Re(x)))

  #Adjust for non-natural base
  if (!missing(base)) {
    LOGB <- Log(base)
    LOG  <- LOG/LOGB
    if (length(base) == 1)   { LOGB <- rep(LOGB, length(LOG)) }
    for (i in 1:length(LOG)) { if (LOGB[i] == 0) LOG[i] <- ifelse(LOG[i] == 0, 1, NA) } }

  #Set back to real vector
  if (all(Im(LOG) == 0))     { LOG <- Re(LOG) }

  LOG }





