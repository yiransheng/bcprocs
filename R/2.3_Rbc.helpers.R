#' Helper function to compute the size(\eqn{\alpha}) of a Rbc1 procedure
#'
#' @param theta0 control Bernoulii parameter
#' @param k number of treatments
#' @param c design constant c
#' @param n number of maximum observations from each treatment, either a integer of a vector of intergerst with length k+1
#' @param d0 design constant \eqn{\delta_0}, \eqn{\theta \leq \theta_0 + \delta_0} is considered ineffective
#' @examples
#' size(0.2, 3, 8, 60) # roughly 0.05

Size <- function (theta0, k, c, n, d0=0) {
  size <- 0
  theta1 <- theta0 + d0
  for (x1 in 0:n) {
    y <- x1 + c
    if (y <= n) {
      # max(x(i)) > x1 + c
      size <- size + dbinom(x1, n, theta0) * (1 - pbinom(y, n, theta1)^k)
    }
  }
  size
}

#' Helper function to compute the power(\eqn{1-\beta}) of a Rbc1 procedure
#'
#' @param theta0 control Bernoulii parameter
#' @param k number of treatments
#' @param c design constant c
#' @param n number of maximum observations from each treatment, either a integer of a vector of intergerst with length k+1
#' @param d0 design constant \eqn{\delta_0}, \eqn{\theta \leq \theta_0 + \delta_0} is considered ineffective
#' @param d1 design constant \eqn{\delta_1}, \eqn{\theta \leq \theta_0 + \delta_1} is considered at most marginally effective
#' @param d2 design constant \eqn{\delta_2}, \eqn{\theta \geq \theta_0 + \delta_1} is considered effective
#' @examples
#' power(0.2, 3, 8, 60) # roughly 0.75
Power <- function (theta0, k, c, n, d0=0, d1=0.05, d2=0.2) {
  power <- 0
  theta2 <- theta0 + d2
  theta1 <- theta0 + d1
  for (xk in 0:n) {
    # choosing k correcly given max(xi)=xk, breaking tie randomly
    pcsk <- 0
    for (tie in 1:k) {
      pcsk <- pcsk + (
        choose(k-1, tie-1) *
          pbinom(xk-1, n, theta1)^(k-tie) *
          dbinom(xk, n, theta1)^(tie-1) *
          dbinom(xk, n, theta2)
      ) / tie # we break the tie randomly, so prob of choosing correctly is inverse to number of ties
    }
    # xk >= x(k-1) >= ... >= x1, xk >= x0 + c
    power <- power + (
      pbinom(xk-c-1, n, theta0) * # prob of not selecting control
      pcsk                        # prob of correctly select pi_k
    )
  }
  power
}

solve.c <- function (n, k=2, a=0.05, b=0.80) {
  c <- 0
  repeat {
    s <- size(0.2, k, c, n)
    if(s <= a) {
      break
    } else {
      c <- c + 1
    }
  }
  return(c)
}

solve.n <- function (k=2, a=0.05, b=0.80) {
  n <- 1
  repeat {
    c <- solve.c(n, k, a, b)
    p <- power(0.2, k, c, n)
    if (p >= b) {
      break
    } else {
      n <- n+1
    }
  }
  return(c(n,c))
}

size.slow <- function (theta0, k, c, n, d0=0) {
  p <- c(theta0, rep(theta0+d0, k))
  mbinom.fn <- multivariate.binom.gen(rep(0, k+1), n)
  size <- 0
  repeat {
    x <- mbinom.fn()
    if(is.null(x)) break
    if(max(x[-1]) > x[1]+c) {
      probs <- Map(function(p, k) { dbinom(k, n, p) }, p, x) # prob vector of observing x, assuming independence
      size <- size + prod(unlist(probs))
    }
  }
  size
}

power.slow <- function (theta0, k, c, n, d0=0, d1=0.05, d2=0.20) {
  p <- c(theta0, rep(theta0+d1, k-1), theta0+d2)
  mbinom.fn <- multivariate.binom.gen(rep(0, k+1), n)
  power <- 0
  repeat {
    x <- mbinom.fn()
    if(is.null(x)) break
    xm <- max(x[-1])
    if(xm > x[1]+c && x[k+1]==xm) {
      ties <- sum(x[-1] == xm)
      probs <- Map(function(p, k) { dbinom(k, n, p) }, p, x) # prob vector of observing x, assuming independence
      power <- power + prod(unlist(probs))/ties
    }
  }
  power
}

# generator function to iterate over the sample space of multivariate binominal distribution
# Input: a vector of non-negative integers, eg. c(1,0, 1, 0) and binominal param n
# Output: a generator iterates over start until c(n, n, ..., n), after which returns NULL

multivariate.binom.gen(start, n) %as% {
  value <- as.integer(start)
  value[1] <- value[1] - 1
  first <- value
  k <- length(value)
  nxt <- function (vec) {
    vec[1] <- vec[1] + 1
    k <- length(vec)
    for (i in 1:k) {
      if(vec[i] > n) {
        if(i == k) stop('Out bound.')
        vec[i] <- vec[i] - n - 1
        vec[i+1] <- vec[i+1] + 1
      }
    }
    vec
  }
  function (reset=FALSE) {
    if(reset) {
      value <<- first
      return(value)
    } else if ( all(value >= rep(n, k)) ) {
      return(NULL)
    } else {
      value <<- nxt(value)
      return(value)
    }
  }
}
