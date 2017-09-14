#' @importFrom ggplot2 autoplot
#' @export autoplot
autoplot <- autoplot


#' Convert the Number of Seconds to Time
#'
#' @param seconds A vector containing the amount of time in seconds.
#' @return A vector list with entries in `day`, `hour`, `minute`, and `seconds`
#' @export
#' @details
#' This function is based off of the lubridate function of a similar name.
convert_seconds_to_time = function(seconds) {
    span = as.double(seconds)
    remainder = abs(span)

    structure(list("day"     = remainder %/% (86400), # 60 s * 60 m * 24 h
                   "hour"    = remainder %/% (3600),  # 60 s * 60 m
                   "minute"  = remainder %/% (60),    # 60 s
                   "second" = remainder %%  (60)),   # Number of seconds
              class = "seconds_to_time")
}

format.seconds_to_time = function(x, ...){
    paste0("Days: ", x$day,", Hours: ", x$hour, ", Minutes: ", x$minute, ", Seconds: ", x$second, collapse="\n")
}

#' @export
print.seconds_to_time = function(x, ...) {
    cat(format(x))
}
