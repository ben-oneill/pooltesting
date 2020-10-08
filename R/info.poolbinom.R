#' Information Matrix for the Pooled Binomial Distribution
#'
#' \code{info.poolbinom} returns the Fisher information matrix for the pooled binomial distribution
#'
#' This function computes the Fisher information matrix for the pooled binomial distribution.  The user gives a count vector \code{poolcount} and
#' a poolsize vector \code{poolsizet} that give the number of pools and their respective sizes.  (If the latter is omitted then the pool sizes are
#' taken to increase sequentially as $s=1,2,3,...$.)  The distribution can be parameterised either by the parameter \code{prev} (the prevalence
#' probabilities) or by the parameter \code{CLLprev} (the complementary-log-log of the prevalence probabilities).  The distribution also uses an
#' intensity parameter \code{intensity} which is also a scalar parameter.
#'
#' @usage \code{info.poolbinom()}
#' @param poolcount Vector of the number of pools (columns correspond to pool sizes; rows are for vectorisation of inputs)
#' @param poolsize Vector of the size of pools (columns correspond to pool sizes; rows are for vectorisation of inputs)
#' @param prev The prevalence probability
#' @param CLLprev The transformed prevalence probability
#' @param intensity A positive intensity parameter
#' @return The Fisher information matrix for the distribution.

info.poolbinom <- function(poolcount, poolsize = NULL, prev = NULL, CLLprev = NULL, intensity = 1) {

  #Check input poolcount
  if (!is.vector(poolcount))                { stop('Error: poolcount should be a vector of count values') }
  if (!is.numeric(poolcount))               { stop('Error: poolcount should be a vector of count values') }
  if (any(as.integer(poolcount) != poolcount))   { stop('Error: poolcount should contain only non-negative integer values') }
  if (min(poolcount) < 0)                   { stop('Error: poolcount should contain only non-negative integer values') }
  COLS <- length(poolcount)

  #Check input poolsize
  if (missing(poolsize)) { poolsize <- 1:COLS } else {
    if (!is.vector(poolsize))          { stop('Error: poolsize should be a vector of positive integers') }
    if (!is.numeric(poolsize))         { stop('Error: poolsize should be a vector of positive integers') }
    if (any(as.integer(poolsize) != poolsize))   { stop('Error: poolsize should contain only positive integer values') }
    if (min(poolsize) < 1)             { stop('Error: poolsize should contain only positive integer values') } }

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
    PROBS <- 1-exp(-(poolsize^intensity)*exp(CLLprev))
    RR    <- poolcount*(poolsize^(2*intensity))*(1-PROBS)/PROBS
    I11   <- exp(2*CLLprev)*sum(RR)
    I12   <- -exp(3*CLLprev - exp(CLLprev))*sum(RR*(log(poolsize)^2))
    I22   <- exp(CLLprev + exp(CLLprev))*sum(RR*log(poolsize))
    INFO  <- matrix(c(I11, I12, I12, I22), byrow = TRUE, nrow = 2, ncol = 2) } else {
    PROBS <- 1-(1-prev)^(poolsize^intensity)
    RR    <- poolcount*(poolsize^(2*intensity))*(1-PROBS)/PROBS
    I11   <- sum(RR)/(1-prev)^2
    I12   <- (log(1-prev)^2)*sum(RR*(log(poolsize)^2))
    I22   <- -log(1-prev)*sum(RR*log(poolsize))/(1-prev)
    INFO  <- matrix(c(I11, I12, I12, I22), byrow = TRUE, nrow = 2, ncol = 2) }

  #Output information matrix
  INFO }
