#' Where is the max()
#'
#' Determines the location, i.e., index of the (first) maximum of a numeric vector, break ties randomly; borrowed from nnet package
#'
#' @param x numeric vector
#' @return index of max value, break ties randomly

which.is.max <- function (x)
{
  y<- seq_along(x)[x == max(x)]
  if (length(y) > 1L)
    sample(y, 1L)
  else y
}
