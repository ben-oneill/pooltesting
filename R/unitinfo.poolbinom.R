#' Unit Information for the Pooled Binomial Distribution
#'
#' \code{unitinfo.poolbinom} returns the unit information for the pooled binomial distribution
#'
#' This function computes the unit information for the pooled binomial distribution.  The user gives a poolsize value \code{poolsize} for the
#' unit.  The distribution can be parameterised either by the parameter \code{prev} (the prevalence probabilities) or by the parameter
#' \code{CLLprev} (the complementary-log-log of the prevalence probabilities).  The distribution also uses an intensity parameter
#' \code{intensity} which is also a scalar parameter.  The output of the function is the unit information provided by a unit in a pool of the
#' specified pool size.
#'
#' @usage \code{unitinfo.poolbinom()}
#' @param poolsize The poolsize value (positive integer)
#' @param prev The prevalence probability
#' @param CLLprev The transformed prevalence probability
#' @param intensity A positive intensity parameter
#' @return The unit information for a unit in a pool of the specified pool-size.

unitinfo.poolbinom <- function(poolsize, prev = NULL, CLLprev = NULL, intensity = 1) {

  #Check input poolsize
  if (!is.vector(poolsize))            { stop('Error: poolsize should be a positive integer') }
  if (!is.numeric(poolsize))           { stop('Error: poolsize should be a positive integer') }
  if (as.integer(poolsize) != poolsize)   { stop('Error: poolsize should be a positive integer') }
  if (poolsize < 1)                    { stop('Error: poolsize should be a positive integer') }

  #Check inputs prev and CLLprev
  if (missing(prev) & missing(CLLprev)) { stop('Error: must specify either prev or CLLprev') }
  if (!missing(prev)) {
    if (!is.vector(prev))              { stop('Error: prev should be a single probability value') }
    if (!is.numeric(prev))             { stop('Error: prev should be a single probability value') }
    if (length(prev) != 1)             { stop('Error: prev should be a single probability value') }
    if (prev < 0)                      { stop('Error: prev should be a probability value (between zero and one)') }
    if (prev > 1)                      { stop('Error: prev should be a probability value (between zero and one)') }
    CLL <- FALSE }
  if (!missing(CLLprev)) {
    if (!is.vector(CLLprev))           { stop('Error: CLLprev should be a single numeric value') }
    if (!is.numeric(CLLprev))          { stop('Error: CLLprev should be a single numeric value') }
    if (length(CLLprev) != 1)          { stop('Error: CLLprev should be a single numeric value') }
    CLL <- TRUE }
  if (!missing(prev) & !missing(CLLprev)) {
    CLLprob2 <- log(-log(1-prev))
    mse      <- (CLLprev - CLLprob2)^2
    if (mse > 1e-10)                   { stop('Error: specify prev or CLLprev but not both') } else {
      warning('Warning: specify prev or CLLprev but not both') }
    CLL <- TRUE }

  #Check input intensity
  if (!is.vector(intensity))           { stop('Error: intensity should be a single positive value') }
  if (!is.numeric(intensity))          { stop('Error: intensity should be a single positive value') }
  if (length(intensity) != 1)          { stop('Error: intensity should be a single positive value') }
  if (intensity <= 0)                  { stop('Error: intensity should be a single positive value') }

  #Compute information matrix
  if (CLL) {
    PROB  <- 1-exp(-(poolsize^intensity)*exp(CLLprev))
    INFO  <- exp(2*CLLprev)*(poolsize^(2*intensity-1))*(1-PROB)/PROB } else {
    PROB <- 1-(1-prev)^(poolsize^intensity)
    INFO  <- (poolsize^(2*intensity-1))*(1-PROB)/(PROB*(1-prev)^2) }

  #Output information matrix
  INFO }
