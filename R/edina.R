#' Construct an EDINA object
#'
#' Constructor function to satisfy the EDINA object definition
#' @inheritParams edina
#' @param x            `edina` object.
#' @param timing       Number of Seconds that have elapsed since run time.
#' @param sample_OR    Odds Ratio based on the trial matrix.
#' @param model_fit    How well the model fit according to the heuristic
#' @keywords internal
new_edina = function(x, k, burnin, chain_length, timing, sample_OR, model_fit) {
    structure(list("edina"        = x,
                   "k"            = k,
                   "burnin"       = burnin,
                   "chain_length" = chain_length,
                   "timing"       = timing,
                   "sample_OR"    = sample_OR,
                   "model_fit"    = model_fit),
                   class = "edina")
}


#' Construct a Summary EDINA object
#'
#' Constructor function to satisfy the Summary EDINA object definition
#' @param m_gs    Mean of the Guessing Values
#' @param m_ss    Mean of the Slipping Values
#' @param m_pi    Mean of the PIs
#' @inheritParams edina
#' @inheritParams new_edina
#' @keywords internal
new_edina_summary = function(m_gs, m_ss, m_pi,
                             k, burnin, chain_length, timing, sample_OR, model_fit) {
    structure(list("m_gs"         = m_gs,
                   "m_ss"         = m_ss,
                   "m_pi"         = m_pi,
                   "k"            = k,
                   "burnin"       = burnin,
                   "chain_length" = chain_length,
                   "timing"       = timing,
                   "sample_OR"    = sample_OR,
                   "model_fit"    = model_fit),
              class = "summary_edina")
}


#' Model Heuristic used for Model Selection
#'
#' Computes the model heuristic for model selection.
#' @param sample_OR   A `matrix` with dimensions \code{J x J}
#'                    containing the empirical odds ratio.
#' @param bayesian_OR An `array` with dimensions \code{J x J}.
#'                    containing the bayesian odds ratio.
#' @return A `double` that is the mean.
#' @export
model_heuristic = function(sample_OR, bayesian_OR) {

    obj = matrix(0, ncol(sample_OR), ncol(sample_OR))

    chain_length = dim(bayesian_OR)[3]

    for ( i in seq_len(chain_length) ) {
        obj = obj + as.matrix( (bayesian_OR[,,i]) > sample_OR)*1
    }

    d = obj / dim(bayesian_OR)[3]

    mean(d[upper.tri(d)] < 0.05 | d[upper.tri(d)] > 0.95)
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


    sample_OR = OddsRatio(nrow(data), ncol(data), data)

    new_edina(edina_model, k,
              burnin, chain_length,
              time_info,
              sample_OR,
              model_heuristic(sample_OR, edina_model$ORs)
              )
}

print.edina = function(x, ...){
    cat("The edina function for K = ", x$k, " took ", convert_seconds_to_time(x$timing[3]))

}

#' Summarizing the EDINA Object
#'
#' @param object An `edina` object
#' @export
summary.edina = function(object, ...) {

    m_gs = apply(object$edina$GS, 1, mean)
    # m_gs
    m_ss = apply(object$edina$SS, 1, mean)
    # m_ss
    m_pi = apply(object$edina$PIs, 1, mean)

    new_edina_summary(m_gs,
                      m_ss,
                      m_pi,
                      object$k,
                      object$burnin,
                      object$chain_length,
                      object$timing,
                      object$sample_OR,
                      object$model_fit)
}

#' Printing out the Summary EDINA Object
#'
#' @param x A `summary_edina` object
#' @export
print.summary_edina = function(x, ...) {
    cat("The means of the Guessing Parameters are: ", fill = TRUE)
    print(x$m_gs)
    cat("\n")

    cat("The means of the Slipping Parameters are: ",fill = TRUE)
    print(x$m_ss)
    cat("\n")

    cat("The means of the PIs are: ", fill = TRUE)
    print(x$m_pi)
    cat("\n")

    cat("The model fit under the heuristic is: ", fill = TRUE)
    print(x$model_fit)
    cat("\n")
}
