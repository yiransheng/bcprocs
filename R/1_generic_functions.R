#' Vector-at-a-time sampling for procedures in this package
#'
#' @param proc S3 Class object for a specific selection procedure (\code{Rbc1}, \code{Rbc2}, \code{Rhbc1} or \code{Rhbc2})
#' @param vec an observation vector or 0/1 (alternative: TRUE/FALSE) of size k + 1
#' @return S3 Class object of the same procedure type
#' @export
vector.at.a.time <- function (proc, vec, ...) UseMethod("vector.at.a.time")
vector.at.a.time.default <- function () "Unknown class"



#' Generic function to determin whether or not a procedure has finished
#'
#' @param proc S3 Class object for a specific selection procedure (\code{Rbc1}, \code{Rbc2}, \code{Rhbc1} or \code{Rhbc2})
#' @return \code{TRUE} or \code{FALSE}
#' @export
is.done <- function (proc, ...) UseMethod("is.done")
is.done.default <- function () "Unknown class"



#' Generic function to run procedure with entire data set instead of sequentially (a.k.a. vector.at.a.time)
#'
#' @param proc S3 Class object for a specific selection procedure (\code{Rbc1}, \code{Rbc2}, \code{Rhbc1} or \code{Rhbc2})
#' @param data data.frame object
#' @return S3 Class object of the same procedure type
#' @export
run.experiment <- function(proc, data, ...) UseMethod("run.experiment")
run.experiment.default <- function (x, data) {
  if (class(data) != 'data.frame') {
    stop('Require data to be data.frame')
  }
  for (m in 1:x$nmax) {
    x <- vector.at.a.time(x, data[m, ])
    if(is.done(x)) break
  }
  x
}

#' Generic function to run stage one of two stage procedures
#'
#' @param proc S3 Class object for a specific selection procedure (\code{Rbc1}, \code{Rbc2}, \code{Rhbc1} or \code{Rhbc2})
#' @param \code{...} addtional arguments
#' @return S3 Class object of the same procedure type
#' @export
stage.one <- function(proc, ...) UseMethod("stage.one")
stage.one.default <- function(...) "Unknown class"

stage.two <- function(proc, ...) UseMethod("stage.two")
stage.two.default <- function(...) "Unknown class"
