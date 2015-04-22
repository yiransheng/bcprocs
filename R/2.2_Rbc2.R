#' S3 Class constructor for Rbc2 procedure
#'
#' @param k number of treatments
#' @param n number of maximum observations to take from each treatment
#' @param c design parameter c, cutoff point for choosing a treatment over control
#' @return S3 Class object of the same procedure type
#' @examples
#' # Simulation Example
#' # Choose: k=3, P1=0.95, P2=0.95, theta_0 = 0.2, delta_1 = 0.05, delta_2 = 0.2
#' # From Table 2: we have, n = 128, c=13
#' n <- 128
#' k <- 3
#' proc <- Rbc2(k=k, c=13, n=n)
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
#' @name Rbc2
Rbc2(k, n, c) %when% {
  k > 1L
  c >= 0L
  n > 0L
} %as% {
  structure(list(
    c=c, # cutoff to select control
    k=k, # number of treatments
    M=0L, # M'th step
    nmax=n, # maximum allowed number of observations to take from each treatment
    terminated = rep(FALSE, k+1L), # vector to record which treatments we will no longer take observations from
    nobs = rep(0L, k+1L), # total number of observations from treatment i at step M
    y=rep(0L, k+1L), # number of successes from treatment i
    decision=NULL # final selection, 0 for control and k indicates kth treatment
  ), class = "Rbc2")
}

Rbc2(k, n, c) %when% {
  (!is.numeric(k) || k <= 1L)
} %as% {
  stop("k: number of treatments must at least be 2")
}

Rbc2(k, n, c) %when% {
  (!is.numeric(c) || c < 0L)
} %as% {
  stop("c must be greater or equal to 0")
}

Rbc2(k, n, c) %when% {
  (!is.numeric(n) || n <=0L)
} %as% {
  stop("n must a positive integer")
}


#' Vector-at-a-time sampling for Rbc2 procedure
#'
#' @param x procedure object of class Rbc2
#' @param vec an observation vector or 0/1 (alternative: TRUE/FALSE) of size k + 1
#' @return S3 Class object of the same procedure type (Rbc2)
#' @examples
#' proc <- Rbc2(2, 56, 8)
#' # proc$k == 2
#' # is.done(proc) == FALSE
#' proc <- vector.at.a.time(proc, c(0,1,1))
#' summary(proc)
#' @name vector.at.a.time.Rbc2
vector.at.a.time.Rbc2 <- function (x, vec) {
  if (x$M >= x$nmax) return(x)
  if (sum(x$terminated) == x$k + 1L) return(x)
  effectiveObs <- vec[!x$terminated]   # allow NA observations only for terminated treatments
  if (any(is.na(effectiveObs)) || any(effectiveObs != TRUE & effectiveObs != FALSE)) {
    return(x)
  }
  # update experiment observations
  x$M <- x$M + 1
  for (i in 1L:(x$k+1L)) {
    if(!x$terminated[i]) {
      # update number of observations and successes for non-terminated treatments
      x$nobs[i] <- x$nobs[i] + 1L
      x$y[i] <- x$y[i] + !!vec[i] # coerce vec[i] to bool
    }
  }
  jmax <- which.is.max(x$y[-1])+1 # most promising treatment so far (excluding control)
  # stopping rules
  shouldStop <- FALSE
  # condition (3.4) & (3.5)
  if (x$y[jmax] > x$y[1] + x$c + x$nmax - x$nobs[1]) {
    shouldStop <- TRUE
    for (i in 2L:(x$k+1L)) {
      if (i != jmax) {
        shouldStop <- shouldStop & (x$y[jmax] >= x$y[i] + x$nmax - x$nobs[i])
      }
    }
    x$decision <- jmax - 1
  }
  # condition (3.6)
  if (x$y[1] + x$c >= x$y[jmax] + x$nmax - x$M) {
    shouldStop <- TRUE
    x$decision <- 0
  }
  if (shouldStop) {
    x$terminated <- rep(TRUE, x$k+1L) # mark all treatments terminated
    return(x)
  } else {
    x$decision <- NULL
  }
  # sampling rules
  for (i in 1:(x$k+1L)) {
    # no more than nmax observations from treatment i should be taken
    if (x$nobs[i] >= x$nmax) {
      x$terminated[i] <- TRUE
    }
    # stop taking observations from treatment i when:
    if (x$y[jmax] > x$y[1] + x$c + x$nmax - x$nobs[1]) {
      x$terminated[1] <- TRUE # stop taking observations from control
    }
    for (j in (i+1L):(x$k+1L)) {
      if(j > i & i > 1L & i<x$k+1L) {
        d <- x$y[i] - x$y[j]
        if(d > 0) {
          x$terminated[j] <- x$terminated[j] || (d >= x$nmax - x$nobs[j])
        } else {
          x$terminated[i] <- x$terminated[i] || (-d >= x$nmax - x$nobs[i])
        }
      }
    }
  }
  x
}


#' Check if Rbc2 procedure has completed
#'
#' @param x procedure object of class Rbc2
#' @return Boolean
#' @examples
#' proc <- Rbc2(2, 56, 8)
#' is.done(proc)
#' # FALSE
#' @name is.done.Rbc2
is.done.Rbc2 <- function(x) {
  x$M >= x$nmax || all(x$terminated)
}

#' Summary of Rbc2 procedure
#'
#' @param x procedure object of class Rbc2
#' @return Table
#' @examples
#' # Sample output
#' Procedure Rbc2, see Buzaianu and Chen (2008)
#' Procedure completed: select Treatement 3
#' Used 365 observations in total, saving: 28.71%
#' At Step:  97
#' Control Treatment 1 Treatment 2 Treatment 3
#' No. Obs            91          97          80          97
#' No. Successes      21          45          15          76
#' Stopped             1           1           1           1
summary.Rbc2 <- function(x) {
  cat('Procedure Rbc2, see Buzaianu and Chen (2008)\n')
  isComplete <- sum(x$terminated) == x$k+1
  if (isComplete) {
    cat('Procedure completed: ')
    if (x$decision == 0) {
      cat('select Control.\n')
    } else {
      cat('select Treatement', x$decision, '\n')
    }
    cat('Used', sum(x$nobs), 'observations in total, saving:', sprintf("%1.2f%%", 100-100*sum(x$nobs)/(x$nmax*x$k+x$nmax)), '\n')
  } else {
    cat('Procedure not completed\n')
  }
  cat('At Step: ', x$M, '\n')
  z <- as.table( rbind(x$nobs, x$y) )
  z <- rbind(z, x$terminated)
  rownames(z) <- c('No. Obs', 'No. Successes', 'Stopped')
  colnames(z) <- c('Control', paste('Treatment', seq(1, x$k)))
  z
}

