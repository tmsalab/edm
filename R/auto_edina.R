#' Auto EDINA model selection routine
#'
#' Automatically select an appropriate $K$ dimension for a $Q$ matrix
#' under the Exploratory Determinatistic Input, Noise And gate (EDINA) Model.
#'
#' @param data         Binary responses to assessements in `matrix`
#'                     form with dimensions \eqn{N \times J}{N x J}.
#' @param k            Number of Attribute Levels as a positive `integer`.
#' @param burnin       Number of Observations to discard on the chain.
#' @param chain_length Length of the MCMC chain
#' @return An `auto_edina` object
#' @export
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
#' edina_models = auto_edina(trial_matrix, k = 1:2)
#' }
auto_edina = function(data, k = 2:4,
                      burnin = 20000, chain_length = 10000,
                      save_results = FALSE, save_filename = "edina_model_data") {

    stopifnot(is.logical(save_results))

    ## Note:
    # Chain length is adjusted in edina

    # Compute the number of _K_ to estimate
    num_k = length(k)
    num_j = ncol(data)

    message("The estimated runtime for searching this space is: ", format(convert_seconds_to_time(sum(2^(k+4)))))

    if(save_results) {
        # Get the length of the string e.g. 200 => 3
        nlen = max(nchar(k))

        # Format string
        save_file_ids = sprintf(paste0("%0",nlen,"d"), k)
    }

    # Setup storage for EDINA Object
    outobj = vector('list', num_k)

    heuristics = rep(NA, num_k)

    for(i in seq_along(k)) {
        k_idx = k[i]
        message("Working on k = ", k_idx)
        message("Estimated runtime is: ", format(convert_seconds_to_time(2^(k_idx+4))))

        modeled_value = edina(data,
                              k = k_idx,
                              burnin = burnin,
                              chain_length = chain_length)

        modeled_value_summary = summary(modeled_value)

        # Launch job
        outobj[[i]] = modeled_value_summary

        heuristics[i] = outobj[[i]][["model_fit"]]

        message("Time Elapsed: ",  outobj[[i]][["timing"]][3])

        if(save_results) {
            # Create a single edina object
            ## consider doing an `assign()` to match the filename
            edina_obj = outobj[[i]]
            save(edina_obj, file = paste0(save_filename, "_", save_file_ids[i],".rda"))
        }
    }

    best_model_id = which.min(heuristics)
    # Output all EDINA objects
    structure(list("edina_models" = outobj,
                   "best_fit" = c("best_model_id" = best_model_id,
                                  "best_k" = k[best_model_id],
                                  "best_model_fit" = heuristics[best_model_id]),
                   "model_fits" = heuristics,
                   "k" = k,
                   "j" = num_j
                   )
              , class = "auto_edina" )
}

#' @export
print.auto_edina = function(x, ...){
    cat("The results of searching Q-matrices between", min(x$k), "and", max(x$k), "...\n")
    cat("Best Fit Model K =", x$best_fit['best_k'],
           "with value", round(x$best_fit['best_model_fit'], 2), "\n")
}

#' Graph the Auto EDINA Object
#'
#' Presents either the fitting of model heuristics or the evolution of parameters
#' on a graph
#' @importFrom ggplot2 autoplot ggplot geom_line geom_point geom_vline facet_wrap labs aes theme_bw
#' @param object An `auto_edina` object.
#' @param type   Kind of graph to display. Valid types: `"selection"`, `"guessing"`, or `"slipping"`.
#' @export "autoplot.auto_edina"
autoplot.auto_edina = function(object, type = "selection", ...) {

    switch(type,
           "selection" = model_selection_graph(object, ...),
           "guessing"  = parameter_evolution_graph(object, type = type, ...),
           "slipping"  =  parameter_evolution_graph(object, type = type, ...),
            stop('Only the following types are valid: `"selection"`, `"guessing"`, or `"slipping"`')
          )

}

#' @export
model_selection_graph = function(x, ...){
    UseMethod("model_selection_graph", x)
}

#' @export
model_selection_graph.auto_edina = function(x, ...) {
    df = data.frame(k = x$k, Heuristic = x$model_fits)

    ggplot(df, aes(x = k, y = Heuristic)) +
        geom_line() +
        geom_point(colour="blue") +
        geom_point(data = df[x$best_fit['best_model_id'], ],
                   colour="red", size = 3) +
        labs(title = "Auto EDINA Model Selection",
             subtitle = paste0("Best Fit Model K=", x$best_fit['best_k'],
                               " with value ", round(x$best_fit['best_model_fit'], 2)  )) +
        theme_bw()
}

#' @export
model_selection_graph.default = function(x, ...){
    stop("Please supply an `auto_edina` object.")
}

#' @export
parameter_evolution_graph = function(x, type = c("guessing", "slipping"), ...) {
    UseMethod("parameter_evolution_graph", x)
}

#' @export
parameter_evolution_graph.auto_edina = function(x, type = c("guessing", "slipping"), ...) {

    J = x$j
    nmodels = length(x$edina_models)

    # Get the length of the string e.g. 200 => 3
    nlen = nchar(J)


    if(type == "guessing") {
        param_name = "g"
        param_estimates = sapply(x$edina_models, `[[`, 1)[, 1]
    } else {
        param_name = "s"
        param_estimates = sapply(x$edina_models, `[[`, 1)[, 2]
    }

    # Potentially add pis class? unlist(m_pi))


    o = data.frame(k        = rep(seq_len(nmodels), each = J),
                   method   = rep(paste0(param_name,
                                         sprintf(paste0("%0", nlen, "d"), seq_len(J))),
                                  nmodels),
                   estimate = param_estimates
    )

    ggplot(o, aes(x = k, y = estimate)) +
        geom_point() + geom_line() +
        geom_vline(xintercept = x$best_fit["best_k"], color = "red") +
        facet_wrap( ~ method) +
        labs(
            title = paste0("Evolution of the ", type, " Parameters"),
            subtitle = "Over Changes in Q Matrix's K Dimension",
            y = paste0("Parameter Estimate of the ", type, " Parameter"),
            x = "Q Matrix of a given K Dimension"
        ) +
        theme_bw()
}

#' @export
parameter_evolution_graph.default = function(x, type = c("guessing", "slipping"), ...){
    stop("Please supply an `auto_edina` object.")
}
