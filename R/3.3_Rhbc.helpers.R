Sizeh <- function (theta0, k, c1, c2, n1, n2, d0=0) {
  theta1 <- theta0 + d0
  p.pass.I <- 0   # prob. of passing stage I with x01, xk1
  p.pass.II <- 0  # conditional prob. of not selecting control after stage II, given x01, xk1
  size <- 0
  for (x01 in 0:n1-c1-1) {
    for (xk1 in (x01+c1+1):n1) {
      p.pass.I <-
        dbinom(x01, n1, theta0) * k * dbinom(xk1, n1, theta1) * pbinom(xk1, n1, theta1)^(k-1)
      p.pass.II <- 0
      for (xk2 in 0:n2) {
        p.pass.II <- p.pass.II +
          pbinom(xk1+xk2-c2-x01-1, n2, theta0) * dbinom(xk2, n2, theta1)
      }
      size <- size + p.pass.I * p.pass.II
    }
  }
  size
}

Powerh <- function (theta0, k, c1, c2, n1, n2, d0=0, d1=0.05, d2=0.2) {
  power <- 0
  theta2 <- theta0 + d2
  theta1 <- theta0 + d1
  for (x0 in 0:(n1-c1-1)) {
    for (xk in (x0+c1+1):n1) {
      pcsk <- 0
      for (tie in 1:k) {
        pcsk <- pcsk + (
          choose(k-1, tie-1) *
            pbinom(xk-1, n1, theta1)^(k-tie) *
            dbinom(xk, n1, theta1)^(tie-1) *
            dbinom(xk, n1, theta2)
        ) / tie # we break the tie randomly, so prob of choosing correctly is inverse to number of ties
      }
      for (x02 in 0:n2) {
        power <- power + ( 1 - pbinom(x0+x02+c2-xk, n2, theta2) ) * dbinom(x0, n1, theta0) * dbinom(x02, n2, theta0) *pcsk
      }
    }
  }
  power
}
