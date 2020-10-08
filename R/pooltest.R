#' Fit the Pooled Binomial Model
#'
#' \code{pooltest} returns the fitting pooled binomial model
#'
#' This function fits a pooled binomial model using a formula and dataset specified by the user.  The response variable for the formula
#' shouldbe the count of positive pools in the dataset and the explanatory variables can be any covariates selected by the user.  The
#' user must also specify a variable ```poolcount``` giving the pool counts and a variable ```poolsize``` giving the pool-sizes.  The default
#' model allows for variable intensity values, but the user may also fix the testing intensity so that the model does not allow excess
#' intensity.  The output is a model object containing the information for the model.
#'
#' @usage \code{pooltest}
#' @param formula A formula for the model
#' @param poolcount The name of the variable giving the pool counts (state this input without quote marks)
#' @param poolsize The name of the variable giving the pool-size numbers (state this input without quote marks)
#' @param data The data-frame for the analysis
#' @param fixed.intensity Logical value; if ```TRUE``` the model does not allow excess testing intensity
#' @param base.intensity The value for the base intensity in the test (treated as an offset in the model)
#' @return A model object of classes ```'pooltest'```, ```'glm'``` and ```'lm'```

pooltest <- function(formula, poolcount = NULL, poolsize, data, fixed.intensity = FALSE, base.intensity = 1, ...) {

  #Check inputs
  if (!("formula" %in% class(formula)))               { stop('Error: formula should be a formula') }
  if (!is.data.frame(data))                           { stop('Error: data should be a data frame') }
  if (!is.vector(fixed.intensity))                    { stop('Error: fixed.intensity should be a single logical value') }
  if (!is.logical(fixed.intensity))                   { stop('Error: fixed.intensity should be a single logical value') }
  if (length(fixed.intensity) != 1)                   { stop('Error: fixed.intensity should be a single logical value') }
  if (!is.vector(base.intensity))                     { stop('Error: base.intensity should be a single non-negative value') }
  if (!is.numeric(base.intensity))                    { stop('Error: base.intensity should be a single non-negative value') }
  if (length(base.intensity) != 1)                    { stop('Error: base.intensity should be a single non-negative value') }
  if (base.intensity < 0)                             { stop('Error: base intensity cannot be negative') }

  #Extract formula variables and check response
  dataname <- deparse(substitute(data))
  y <- as.character(formula)[2]
  x <- as.character(formula)[3]
  if (!(y %in% names(data)))                          { stop(paste0('Error: Cannot find response variable ', y, ' in data frame ', dataname)) }
  if (!all.equal(as.integer(data[[y]]), data[[y]]))   { stop(paste0('Error: Response variable ', y, ' should be a count variable (i.e., non-negative integers)')) }
  if (min(data[[y]]) < 0)                             { stop(paste0('Error: Response variable ', y, ' should be a count variable (i.e., non-negative integers)')) }

  #Extract and check input poolcount
  #If this input is missing then all poolcounts are taken to be one and the response vector must be binary
  if (missing(poolcount)) {
    data$n <- 1
    n <- 'n'
    if (max(data[[y]]) > 1)                           { stop(paste0('Error: Response variable ', y, ' cannot be larger than one if no poolcount variable is specified')) } }
  if (!missing(poolcount)) {
    n <- deparse(substitute(poolcount))
    if (!(n %in% names(data)))                        { stop(paste0('Error: Cannot find poolcount variable ', n, ' in data frame ', dataname)) }
    if (!all.equal(as.integer(data[[n]]), data[[n]])) { stop(paste0('Error: Poolcount variable ', n, ' should be a positive count variable (i.e., positive integers)')) }
    if (min(data[[n]]) <= 0)                          { stop(paste0('Error: Poolcount variable ', n, ' should be a positive count variable (i.e., positive integers)')) }
    if (min(data[[n]] - data[[y]]) < 0)               { stop(paste0('Error: Response variable ', y, ' cannot be larger than poolcount variable', n)) } }

  #Extract and check input poolsize
  s <- deparse(substitute(poolsize))
  if (!(s %in% names(data)))                          { stop(paste0('Error: Cannot find poolsize variable ', s, ' in data frame ', dataname)) }
  if (!all.equal(as.integer(data[[s]]), data[[s]]))   { stop(paste0('Error: Poolsize variable ', s, ' should be a positive count variable (i.e., positive integers)')) }
  if (min(data[[s]]) <= 0)                            { stop(paste0('Error: Poolsize variable ', s, ' should be a positive count variable (i.e., positive integers)')) }

  #Edit data frame
  data$Positive        <- data[[y]]
  data$Negative        <- data[[n]] - data[[y]]
  data$BaseIntensity   <- base.intensity*log(data[[s]])
  data$ExcessIntensity <- log(data[[s]])

  #Fit the model and give warnings if there was bad model convergence
  if (fixed.intensity) {
    FORMULA <- formula(paste0('cbind(Positive, Negative) ~ offset(BaseIntensity) + ', x)) } else {
    FORMULA <- formula(paste0('cbind(Positive, Negative) ~ offset(BaseIntensity) + ExcessIntensity + ', x)) }
  MODEL <- glm(FORMULA, data = data, family = binomial(link = cloglog), ...)
  if (!MODEL$converged) { warning(paste0('Attempted fit to glm function did not converge after ', MODEL$iter, ' iterations')) }
  if (MODEL$boundary)   { warning(paste0('Attempted fit to glm function obtained at a boundary point')) }

  #Amend model details
  if (fixed.intensity) {
    MODEL$family$family <- paste0('Pooled Binomial Model with fixed intensity (base intensity = ', base.intensity, ')')    } else {
    MODEL$family$family <- paste0('Pooled Binomial Model with variable intensity (base intensity = ', base.intensity, ')') }
  MODEL$fixed.intensity <- fixed.intensity
  MODEL$base.intensity  <- base.intensity
  MODEL$call            <- sys.call()

  #Add unconstrained model
  FORMULA.UNCONSTRAINED     <- formula(paste0('cbind(Positive, Negative) ~ factor(', s, ') + ', x))
  MODEL$unconstrained.model <- suppressWarnings(glm(FORMULA.UNCONSTRAINED, data = data, family = binomial(link = cloglog)))

  #Set model class
  class(MODEL) <- c('pooltest', class(MODEL))

  #Give output
  MODEL }
