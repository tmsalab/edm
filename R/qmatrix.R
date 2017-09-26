#' Q Matrix
#'
#' Extract a Q matrix
#' @param x An `edina` object
#' @export
extract_q_matrix = function(x, type = c("whole", "avg")) {
    stopifnot(inherits(x, "edina"))

    if(type == "whole") {
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
