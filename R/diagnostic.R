#' Diagnostic Analysis of Pooled Binomial Model
#'
#' \code{diagnostic} returns an ANOVA diagnostic test for a fitting pooled binomial model
#'
#' This function operates on a fitted pooled binomial model to produce an analysis of variance (ANOVA) test comparing the fitted model to the
#' unconstrained model obtained when the poolsize is treated as a factor variable.  This test operates as a diagnostic test for the stipulated
#' form of the positive pool probabilities, as a function of the pool size.  The input for the function is a model of class \code{'pooltest'} and
#' the output is an ANOVA table comparing the actual model and the unconstrained model.
#'
#' @usage \code{diagnostic}
#' @param model A pooled binomial model produced by the \code{pooltest} function (with class \code{'pooltest'})
#' @return An ANOVA table

diagnostic <- function(model) {

  #Check that input is a pooltest object
  if (!('pooltest' %in% class(model)))                 { stop('Error: Input should be a \'pooltest\' model') }

  #Create ANOVA table
  OUTPUT <- anova(model, model$unconstrained.model, test = 'Chisq')
  names(OUTPUT) <- c('Resid. Df', 'Resid. Dev', 'Df', 'Deviance', 'p-value')
  attributes(OUTPUT)$row.names  <- c('Actual Model', 'Unconstrained Model')
  attributes(OUTPUT)$heading[1] <- '\nAnalysis of Deviance Table (Actual Model vs Unconstrained Model)\n'
  if (model$fixed.intensity) {
    attributes(OUTPUT)$heading[2] <- paste0('Model: Pooled Binomial Model with fixed intensity (base intensity = ', model$base.intensity, ')\n')    } else {
    attributes(OUTPUT)$heading[2] <- paste0('Model: Pooled Binomial Model with variable intensity (base intensity = ', model$base.intensity, ')\n') }
  attributes(OUTPUT)$heading[3] <- 'Using chi-squared test of nested models\n'

  #Give output
  OUTPUT }





