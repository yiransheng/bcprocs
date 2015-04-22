#' S3 Class constructor for Rbc1 procedure
#'
#' @param k number of treatments
#' @param n number of maximum observations to take from each treatment
#' @param c design parameter c, cutoff point for choosing a treatment over control
#' @return S3 Class object of the same procedure type
#' @examples
#' # run on a full dataset
#' # data is a data.frame of dimention (128 ,4), with 0/1 values
#'
#' proc <- Rbc1(3, 13, 128, data)
#'
#' # alternative
#' proc <- Rbc1(3, 13, 128)
#'
#'
#' summary(run.experiment(proc, data))
#' # Simulation Example
#' # Choose: k=3, P1=0.95, P2=0.95, theta_0 = 0.2, delta_1 = 0.05, delta_2 = 0.2
#' # From Table 2: we have, n = 128, c=13
#' n <- 128
#' k <- 3
#' proc <- Rbc1(k=k, c=13, n=n)
#' repeat {
#'   p <- c(0.2, runif(k))         # population bernoulie parameters
#'   jmax <- which.max(p[-1]) + 1  # best treatment
#'   if (p[jmax] >= p[1] + 0.2) {  # ensure assumptions in Section (2) are satisfied AND there's a best treatment to select
#'     invalid <- sapply(p[-1], function(x) {
#'       x > p[1] + 0.05 & x < p[1] + 0.2
#'     })
#'     if(all(!invalid)) break
#'   }
#' }
#' step <- 1
#' repeat {
#'   sampleObs <- sapply(p, function(prob) { # generate B(1, p) binary r.v. as observations, from population probs.
#'     rbinom(1, 1, prob)
#'   })
#'   proc <- vector.at.a.time(proc, sampleObs)
#'   step <- step + 1
#'   if (step > n | ( is.done(proc) )) break
#' }
#' summary(proc)
#' @name Rbc1

Rbc1(k, n, c, data=NULL) %when% {
  k > 1L
  c >= 0L
  (length(n) == 1L || length(n) == k+1L)
  n > 0L
} %as% {
  if(length(n) == 1L) n <- rep(n, k+1L)
  me <- structure(list(
    c=c, # cutoff to select control
    k=k, # number of treatments
    M=0L, # M'th step
    n=as.integer(n), # number of observations to take from each treatment
    nmax=max(n),
    nobs = rep(0L, k+1L), # total number of observations from treatment i at step M
    y=rep(0L, k+1L), # number of successes from treatment i
    decision=NULL # final selection, 0 for control and k indicates kth treatment
  ), class = "Rbc1")
  if(!is.null(data)) {
    return( run.experiment(me, data) )
  } else {
    return(me)
  }
}

Rbc1(k, n, c, data=NULL) %when% {
  (!is.numeric(k) || k <= 1L)
} %as% {
  stop("k: number of treatments must at least be 2")
}

Rbc1(k, n, c, data=NULL) %when% {
  (!is.numeric(c) || c < 0L)
} %as% {
  stop("c must be greater or equal to 0")
}

Rbc1(k, n, c, data=NULL) %when% {
  (!is.numeric(n) || n <=0L || length(n) != k+1L)
} %as% {
  stop("n must either be a postive integer or a vector of positive integers length k+1")
}


#' Vector-at-a-time sampling for Rbc1 procedure
#'
#' @param x procedure object of class Rbc1
#' @param vec an observation vector or 0/1 (alternative: TRUE/FALSE) of size k + 1
#' @return S3 Class object of the same procedure type (Rbc1)
#' @examples
#' proc <- Rbc1(2, 56, 8)
#' # proc$k == 2
#' # is.done(proc) == FALSE
#' proc <- vector.at.a.time(proc, c(0,1,1))
#' summary(proc)
#' @name vector.at.a.time.Rbc1
vector.at.a.time.Rbc1 <- function (x, vec) {
  # update experiment observations
  x$M <- x$M + 1L
  for (i in 1L:(x$k+1L)) {
    veci <- vec[i] ; veci.bool <- !!veci
    if(
      x$nobs[i] < x$n[i] &&
      !is.na(veci) &&
      veci == veci.bool # make sure observation is either 0 or 1 or TRUE or FALSE
    ) {
      # update number of observations and successes
      x$nobs[i] <- x$nobs[i] + 1L
      x$y[i] <- x$y[i] + veci.bool  # integer + boolean -> integer
    }
  }
  # required number of observations has been reached
  if (is.done(x)) {
    theta <- x$y / x$nobs
    # decision rules:
    jmax <- which.is.max(theta[-1]) + 1
    if (theta[jmax] <= theta[1] + x$c/x$nmax) {
      x$decision = 0
    } else {
      x$decision = jmax - 1
    }
  }
  x
}


#' Check if Rbc1 procedure has completed
#'
#' @param x procedure object of class Rbc1
#' @return Boolean
#' @examples
#' proc <- Rbc1(2, 56, 8)
#' is.done(proc)
#' # FALSE
#' @name is.done.Rbc1
is.done.Rbc1 <- function(x) {
  all(x$nobs == x$n)
}

#' Summary of Rbc1 procedure
#'
#' @param x procedure object of class Rbc1
#' @return Table
#' @examples
#' # Sample output
#' Procedure Rbc1, see Buzaianu and Chen (2008)
#' Procedure completed: select Treatement 2
#' At Step:  128
#' Control Treatment 1 Treatment 2 Treatment 3
#' No. Obs       128.0000   128.00000 128.0000000   128.00000
#' No. Successes  24.0000    12.00000  93.0000000    76.00000
#' p.esti.         0.1875     0.09375   0.7265625     0.59375
summary.Rbc1 <- function(x) {
  cat('Procedure Rbc1, see Buzaianu and Chen (2008)\n')
  isComplete <- is.done(x)
  if (isComplete) {
    cat('Procedure completed: ')
    if (x$decision == 0) {
      cat('select Control.\n')
    } else {
      cat('select Treatement', x$decision, '\n')
    }
  } else {
    cat('Procedure not completed\n')
  }
  cat('At Step: ', x$M, '\n')
  z <- as.table( rbind(x$nobs, x$y) )
  z <- rbind(z, x$y/x$nobs)
  rownames(z) <- c('No. Obs', 'No. Successes', 'p.esti.')
  colnames(z) <- c('Control', paste('Treatment', seq(1, x$k)))
  z
}
