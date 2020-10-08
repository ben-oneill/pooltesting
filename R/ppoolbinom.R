#' Cumulative distribution function for the pooled binomial distribution
#'
#' \code{ppoolbinom} returns the sumulative distribution function for the pooled binomial distribution.
#'
#' This function computes the CDF or log-CDF for the pooled binomial distribution.  For each computation there must be a count
#' vector \code{poolcount} and a poolsize vector \code{poolsize} that give the number of pools and their respective sizes.  (If the latter
#' is omitted then the pool sizes are taken to increase sequentially as $s=1,2,3,...$.)  There must also be a count vector \code{x} that
#' gives the number of positive pools of each pool size.  For vectorisation of the function the objects \code{poolcount} and \code{x} are matrices
#' instead of vectors, with each row corresponding to one density computation.  The density can be parameterised either by the parameter
#' \code{prev} (the prevalence probabilities) or by the parameter \code{CLLprev} (the complementary-log-log of the prevalence probabilities).
#' The chosen parameter should either be a scalar value that applies to all computations, or a vector with length equal to the number of rows
#' of the \code{poolcount} and \code{x} objects.  The output of the function is a vector of CDF or log-CDF values for the distribution.
#' The density also allows an intensity parameter \code{intensity} which is also a scalar/vector parameter.
#'
#' @usage \code{dpoolbinom()}
#' @param x Vector/matrix of the number of pools (columns correspond to pool sizes; rows are for vectorisation of inputs)
#' @param poolcount Vector/matrix of the number of pools (columns correspond to pool sizes; rows are for vectorisation of inputs)
#' @param poolsize Vector/matrix of the size of pools (columns correspond to pool sizes; rows are for vectorisation of inputs)
#' @param prev A scalar/vector of prevalence probabilities (if a vector, it should have length equal to the \code{poolcount} and \code{x} inputs)
#' @param CLLprev A scalar/vector of transformed prevalence probabilities (if a vector, it should have length equal to the \code{poolcount} and \code{x} inputs)
#' @param intensity A scalar/vector of intensity parameters (if a vector, it should have length equal to the \code{poolcount} and \code{x} inputs)
#' @param log.p Logical value; if \code{TRUE} the function returns the log-CDF instead of the CDF
#' @return A vector of outputs of the CDF or log-CDF values.

ppoolbinom <- function(x, poolcount, poolsize = NULL, prev = NULL, CLLprev = NULL, intensity = 1, log.p = FALSE) {

  #Check input poolcount and put it in matrix form
  if (!is.numeric(poolcount))               { stop('Error: poolcount should be numeric') }
  if (is.array(poolcount)) {
    if (length(dim(poolcount)) > 2)         { stop('Error: poolcount should be a vector or matrix, not a larger array') } }
  if (any(as.integer(poolcount) != poolcount))   { stop('Error: poolcount should contain only non-negative integer values') }
  if (min(poolcount) < 0)                   { stop('Error: poolcount should contain only non-negative integer values') }
  if (!is.array(poolcount)) { poolcount <- matrix(poolcount, nrow = 1) }
  ROWS <- nrow(poolcount);
  COLS <- ncol(poolcount);

  #Check input x and put it in matrix form
  if (!is.numeric(x))                  { stop('Error: x should be numeric') }
  if (is.array(x)) {
    if (length(dim(x)) > 2)            { stop('Error: x should be a vector or matrix, not a larger array') } }
  if (!is.array(x)) { x <- matrix(x, nrow = 1) }
  if (nrow(x) != ROWS)                 { stop('Error: x should be the same dimensions as poolcount') }
  if (ncol(x) != COLS)                 { stop('Error: x should be the same dimensions as poolcount') }

  #Check input poolsize
  if (missing(poolsize)) {
    poolsize <- matrix(1:COLS, byrow = TRUE, nrow = ROWS) } else {
    if (!is.numeric(poolsize))         { stop('Error: poolsize must be numeric') }
    if (is.array(poolsize)) {
      if (length(dim(poolsize)) > 2)   { stop('Error: poolsize should be a vector or matrix, not a larger array') } }
    if (any(as.integer(poolsize) != poolsize))   { stop('Error: poolsize should contain only positive integer values') }
    if (min(poolsize) < 1)             { stop('Error: poolsize should contain only positive integer values') }
    if (!is.array(poolsize)) { poolsize <- matrix(poolsize, nrow = 1) } }

  #Check inputs prev and CLLprev
  if (missing(prev) & missing(CLLprev)) { stop('Error: must specify either prev or CLLprev') }
  if (!missing(prev)) {
    if (!is.numeric(prev))             { stop('Error: prev should be numeric') }
    if (is.array(prev)) {
      if (length(dim(prev)) > 1)         { stop('Error: prev should be a vector, not a larger array') } }
    if (length(prev) == 1) { prev <- rep(prev, ROWS) }
    if (length(prev) != ROWS)          { stop('Error: prev must either be a scalar or a vector with one element for each row of n') }
    if (min(prev) < 0)                 { stop('Error: prev values must be probability values (between zero and one)') }
    if (max(prev) > 1)                 { stop('Error: prev values must be probability values (between zero and one)') }
    CLL <- FALSE }
  if (!missing(CLLprev)) {
    if (!is.numeric(CLLprev))          { stop('Error: CLLprev should be numeric') }
    if (is.array(CLLprev)) {
      if (length(dim(CLLprev)) > 1)      { stop('Error: CLLprev should be a vector, not a larger array') } }
    if (length(CLLprev) == 1) { CLLprev <- rep(CLLprev, ROWS) }
    if (length(CLLprev) != ROWS)       { stop('Error: CLLprev must either be a scalar or a vector with one element for each row of n') }
    CLL <- TRUE }
  if (!missing(prev) & !missing(CLLprev)) {
    CLLprob2 <- log(-log(1-prev))
    mse      <- (CLLprev - CLLprob2)^2
    if (mse > 1e-10)                   { stop('Error: specify prev or CLLprev but not both') } else {
      warning('Warning: specify prev or CLLprev but not both') }
    CLL <- TRUE }

  #Check input intensity
  if (!is.numeric(intensity))          { stop('Error: intensity should be numeric') }
  if (is.array(intensity)) {
    if (length(dim(intensity)) > 1)    { stop('Error: intensity should be a vector, not a larger array') } }
  if (length(intensity) == 1) { intensity <- rep(intensity, ROWS) }
  if (length(intensity) != ROWS)       { stop('Error: intensity must either be a scalar or a vector with one element for each row of n') }

  #Check input log.p
  if (!is.logical(log.p))              { stop('Error: log.p must be a logical value') }
  if (is.array(log.p))                 { stop('Error: log.p must be a single logical value') }
  if (length(log.p) != 1)              { stop('Error: log.p must be a single logical value') }

  #Compute log-density values
  LOGCDF <- rep(NA, ROWS)
  for (i in 1:ROWS) {
    if (CLL) { PROBS <- 1-exp(-(poolsize[i,]^intensity)*exp(CLLprev)) } else {
               PROBS <- 1-(1-prev[i])^(poolsize[i,]^intensity) }
    LOGCDF[i] <- sum(pbinom(x[i,], size = poolcount[i,], prob = PROBS, log.p = TRUE)); }

  #Output values
  if (log.p) { LOGCDF } else { exp(LOGCDF) } }
