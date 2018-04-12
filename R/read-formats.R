#' Create a Generic Matrix Read Function
#'
#' Imports data from a flat file as a matrix with a predefined
#' naming scheme.
#'
#' @param colname_prefix Name to attach to columns read in. Default is `NULL`
#'                       meaning names are not filled in.
#' @param rowname_prefix Name to attach to rows read in. Default is `NULL`
#'                       meaning names are not filled in.
#' @return Provides a function that can be customized.
#' @keywords internal
#' @importFrom utils read.table
#'
#' @noRd
read_psych = function(file, header = FALSE, sep = " ", skip = 0,
                      drop_columns = NULL, colname_prefix = NULL,
                      rowname_prefix = NULL) {

    # Read in data as a data.frame
    generic_df = read.table(file, header = header, sep = sep, skip = skip)

    # Convert to matrix
    generic_matrix = as.matrix(generic_df)

    if(!is.null(drop_columns) & is.integer(drop_columns)) {
        generic_matrix = generic_matrix[,-drop_columns]
    }

    # Number of Items
    j = ncol(generic_matrix)

    # Number of Responses
    n = nrow(generic_matrix)

    # Labeling
    if(!is.null(colname_prefix)) {
        colnames(generic_matrix) = sprintf(paste0(colname_prefix, "%0",
                                                  nchar(j), "i"), seq_len(j))
    } else {
        colnames(generic_matrix) = NULL
    }

    if(!is.null(rowname_prefix)) {
        rownames(generic_matrix) = sprintf(paste0(rowname_prefix,"%0",
                                                  nchar(n), "d"), seq_len(n))
    }

    generic_matrix
}


#' Import an Items Matrix
#'
#' Allows for a dichotmous items matrix to be imported with standard styling.
#'
#' @param file         name, url, or [connections] of the file to be read
#' @param header       a logical value to indicate if the data contains a
#'                     description field as the first row. Default is `FALSE`.
#' @param sep          the field separator character. Values on each line of the
#'                     file are separated by this character. Default is white
#'                     space given by `sep = " "`.
#' @param skip         the number of lines of the data file that should be
#'                     skipped before the data is read. Default is 0.
#' @param drop_columns columns of the item matrix that should be dropped. To
#'                     prevent the first and third columns specify
#'                     `drop_columns = c(1,3)`. By default, all columns are included.
#' @return A `matrix` labeled with row names as `SubjectXX` and column names as `ItemYY`
#' @export
#' @details
#' This function is designed to specifically read in dichotomous item matrices
#' into _R_. The matrix must be structured with the appropriate separator.
read_items = function(file, header = FALSE, sep = " ", skip = 0,
                      drop_columns = NULL) {
    read_psych(file, header = header, sep = sep, skip = skip,
               drop_columns = drop_columns,
               colname_prefix = "Item", rowname_prefix =  "Subject")
}

#' Import a Q Matrix
#'
#' Allows for a dichotmous Q Matrix to be imported with standard styling.
#'
#' @inheritParams read_items
#' @return A `matrix` labeled with row names as `ItemYY` and column names as `SkillZZ`
#' @export
#' @details
#' This function is designed to specifically read in dichotomous Q matrix
#' into _R_. The matrix must be structured with the appropriate separator.
read_qmatrix = function(file, header = FALSE, sep = " ", skip = 0,
                        drop_columns = NULL) {

    a = read_psych(file, header = header, sep = sep, skip = skip,
               drop_columns = drop_columns,
               colname_prefix = "Trait", rowname_prefix =  "Item")

    create_q_matrix(a)
}
