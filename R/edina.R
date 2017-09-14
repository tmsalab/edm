#' Construct an EDINA object
#'
#' Constructor function to satisfy the EDINA object definition
#' @param coefs        Matrix with means of guessing and slipping coefficients.
#' @param pis          Estimated latent classes
#' @param est_q        Estimated Q Matrix
#' @param or_tested    Sample empirical odds ratio compared against simulated
#'                     odds ratio divided by the number of simulations.
#' @param sample_or    Sample empirical odds ratio based on the trial matrix.
#' @inheritParams edina
#' @param timing       Number of Seconds that have elapsed since run time.
#' @keywords internal
new_edina = function(coefs, pis, est_q, or_tested, sample_or, k, burnin, chain_length,  timing) {

    item_nums = paste0("Item", seq_len(nrow(coefs)) )

    colnames(coefs) = c("Guessing", "Slipping")
    rownames(coefs) = item_nums

    colnames(est_q) = paste0("Trait", seq_len(ncol(est_q)))
    rownames(est_q) = item_nums

    structure(list("coefficients" = coefs,
                   "pi_classes"   = pis,
                   "est_q"        = est_q,
                   "or_tested"    = or_tested,
                   "sample_or"    = sample_or,
                   "k"            = k,
                   "burnin"       = burnin,
                   "chain_length" = chain_length,
                   "timing"       = timing),
                   class = "edina")
}


#' Construct a Summary EDINA object
#'
#' Constructor function to satisfy the Summary EDINA object definition
#' @param edina      An `edina` object
#' @param model_fit  Computed model heuristic value
#' @param alpha      The region used in the computation of the heuristic.
#' @keywords internal
new_edina_summary = function(edina, model_fit, alpha) {

    edina[["model_fit"]] = model_fit
    edina[["alpha"]] = alpha

    class(edina) = c("summary_edina", "edina")

    edina
}


#' Model Heuristic used for Model Selection
#'
#' Computes the model heuristic for model selection.
#' @inheritParams summary.edina
#' @return A `double` that is the mean.
#' @export
model_heuristic = function(object, alpha = 0.05) {
    or_tested = object$or_tested

    mean(or_tested[upper.tri(or_tested)] < alpha |
         or_tested[upper.tri(or_tested)] > (1-alpha))
}

#' EDINA estimation routine
#'
#' Exploratory Determinatistic Input, Noise and Gate Model (EDINA)
#'
#' @param data         Binary responses to assessements in `matrix`
#'                     form with dimensions \eqn{N \times J}{N x J}.
#' @param k            Number of Attribute Levels as a positive `integer`.
#' @param burnin       Number of Observations to discard on the chain.
#' @param chain_length Length of the MCMC chain
#' @return An `edina` object that contains:
#'
#' @export
#' @importFrom balamuta is.whole
#' @examples
#' \dontrun{
#' library("tmsadata")
#'
#' # Load data
#' data("trial_matrix", package="tmsadata")
#'
#' # Coerce to matrix
#' trial_matrix = as.matrix(trial_matrix)
#'
#' edina_model = edina(trial_matrix, k = 2)
#' }
#'
edina = function(data, k = 3, burnin = 10000, chain_length = 20000){

    stopifnot(is.matrix(data))

    stopifnot(is.whole(k) && length(k) == 1 && k >= 1)

    stopifnot(is.whole(chain_length) && length(chain_length) == 1)

    time_info = system.time({

        edina_model = edina_Gibbs_Q(data, k,
                                    burnin = burnin,
                                    chain_length = chain_length)

    })[1:3]

    new_edina(edina_model$coefficients,
              edina_model$pis,
              edina_model$est_q,
              edina_model$or_tested,
              edina_model$sample_or,
              k,
              burnin, chain_length,
              time_info
              )
}

#' Printing out the EDINA Object
#'
#' Custom print method for computing the EDINA.
#' @param x An `edina` object
#' @export
print.edina = function(x, ...){
    cat("The EDINA model for K =", x$k,
        "took", format(convert_seconds_to_time(x$timing[3])), "\n")

    cat("\nThe estimated coefficients for the EDINA model are:\n")
    print(x$coefficients)

    cat("\nThe estimated Q matrix is:\n")
    print(x$est_q)
}

#' Summarizing the EDINA Object
#'
#' Determine whether the `edina` object is summarized appropriately.
#' @param object An `edina` object
#' @param alpha  Defining region to indicate the level of extremeness
#'               the data must before the model is problematic.
#' @export
summary.edina = function(object, alpha = 0.05, ...) {

    new_edina_summary(
        object,
        model_fit = model_heuristic(object, alpha),
        alpha = alpha
    )

}

#' Printing out the Summary EDINA Object
#'
#' @param x A `summary_edina` object
#' @export
print.summary_edina = function(x, ...) {
    # Rely upon the specification of the `edina` object in the summary class.
    # NextMethod()
    cat("The EDINA model for K =", x$k,
        "took", format(convert_seconds_to_time(x$timing[3])), "\n")
    cat("\nThe model fit under the heuristic alpha =", x$alpha ,"is:", x$model_fit, "\n")
}
