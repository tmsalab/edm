#' Extract Q Matrix
#'
#' Given a modeling object, extract the Q Matrix
#'
#' @param x    An `edina`, `dina`, `errum`, or `rrum` object
#' @param ...  Additional parameters
#'
#' @return A `matrix` that is either dichotomous or estimated.
#'
#' @rdname extract_q
#' @export
extract_q_matrix = function(x, ...) {
    UseMethod("extract_q_matrix")
}

#' @param binary_q   Classified Q matrix or a rounded Q matrix.
#' @rdname extract_q
#' @export
extract_q_matrix.edina = function(x, binary_q = FALSE) {
    stopifnot(inherits(x, "edina"))
    pull_est_q_matrix(x, binary_q)
}

#' @rdname extract_q
#' @export
extract_q_matrix.errum = function(x, binary_q = FALSE, ...) {
    stopifnot(inherits(x, "errum"))

    pull_est_q_matrix(x, binary_q)
}

#' @rdname extract_q
#' @export
extract_q_matrix.default = function(x, ...) {
    stop("`x` must be a supported type.")
}


pull_est_q_matrix = function(x, binary_q = FALSE) {
    if(binary_q) {
        x$est_q
    } else {
        x$avg_q
    }
}

format_q_matrix = function(x) {
    colnames(x) = paste0("Trait", seq_len(ncol(x)))
    rownames(x) = paste0("Item", seq_len(nrow(x)) )
    x
}
