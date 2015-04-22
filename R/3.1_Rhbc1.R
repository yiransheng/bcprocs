Rhbc1(k, n1, c1, n2, c2, data=NULL) %when% {
  k > 1L
  c1 >= 0L
  n1 > 0L
  c2 >= 0L
  n2 > 0L
} %as% {
  me <- structure(list(
    c=c1, # cutoff to select control
    c2=c2,
    k=k, # number of treatments
    M=0L, # M'th step
    nmax=n1,
    n=rep(as.integer(n1), k+1L),
    n2=n2, # maximum allowed number of observations for stage 2
    nobs = rep(0L, k+1L), # total number of observations from treatment i at step M
    y=rep(0L, k+1L), # number of successes from treatment i
    nobs2 = c(0L, 0L),
    y2 = c(0L, 0L),
    stage = 1L,
    decision=NULL, # final selection, 0 for control and k indicates kth treatment
    candidate=NULL # decision after stage I
  ), class = c("Rhbc1", "Rbc1"))
  if(!is.null(data)) {
    return( run.experiment(me, data) )
  } else {
    return(me)
  }
}

Rhbc1(k, n1, c1, n2, c2, data=NULL) %when% {
  (!is.numeric(k) || k <= 1L)
} %as% {
  stop("k: number of treatments must at least be 2")
}

Rhbc1(k, n1, c1, n2, c2, data=NULL) %when% {
  (!is.numeric(c1) || c1 < 0L  || !is.numeric(c2) || c2 < 0L)
} %as% {
  stop("c1 and c2 must be greater or equal to 0")
}

Rhbc1(k, n1, c1, n2, c2, data=NULL) %when% {
  (!is.numeric(n1) || n1 <=1L || !is.numeric(n2) || n2 < 1L)
} %as% {
  stop("n1, n2 must positive integers")
}

is.done.Rhbc1 <- function(x, stage=1L) {
  if (x$stage == 0L) return(TRUE)
  if (x$stage == 2L && stage == 2L) return(FALSE)
  if (x$stage == 1L && stage == 1L) {
    is.stage.done <- getS3method('is.done', 'Rbc1')
    return(is.stage.done(x))
  }
  if (x$stage == 2L && stage==1L) return(TRUE)
}

stage.one.Rhbc1(x, vec) %when% {
  x$stage == 1L
  "Rbc1" %in% class(x)
} %as% {
  x <- vector.at.a.time(x, vec)
  if (is.done(x, 1L)) {
    class(x) <- c("Rhbc1")   # when stage 1 is complete, remove "Rbc1" class
    x$candidate <- x$decision
    if(x$decision == 0L) {
      x$stage = 0L
    } else {
      x$stage = 2L
    }
  }
  x
}

stage.one.Rhbc1(x, vec) %as% {
  warning("Incorrect stage or corrupted Rhbc2 object.")
  x
}

stage.two.Rhbc1(proc, x12, xk2) %when% {
  proc$stage == 2L
  proc$decision != 0L
  !any(is.na(x12))
  !any(is.na(xk2))
  all(x12 == !!x12)  # x12 are all 0 or 1 | TRUE or FALSE
  all(xk2 == !!xk2)  # xk2 are all 0 or 1 | TRUE or FALSE
} %as% {
  N <- proc$nmax + proc$n2
  if ( (length(x12) + proc$nobs[1]) < N ) stop("Insufficient data")
  if ( (length(xk2) + proc$nobs[(proc$decision+1L)]) < N ) stop("Insufficient data")

  j <- proc$decision
  n2.1 <- N - proc$nobs[1]
  n2.j <- N - proc$nobs[j+1]
  proc$nobs2 <- c(n2.1, n2.j)
  x12 <- x12[1:n2.1]
  xk2 <- xk2[1:n2.j]
  y11 <- proc$y[1]
  yk1 <- proc$y[j+1]
  y12 <- sum(x12) # number of 2nd stage success
  yk2 <- sum(xk2) # number of 2nd stage success
  proc$y2 <- c(y12, yk2)
  if(yk2+yk1 - y11 - y12 > proc$c2) {
    proc$decision <- j
  } else {
    proc$decision <- 0L
  }
  proc$stage <- 0L   # mark proc done
  proc
}

stage.two.Rhbc1(proc, x12, xk2) %when% {
  proc$stage != 2L
} %as% {
  stop("Incorrect stage")
}

stage.two.Rhbc1(proc, ...)  %as% {
  stop("Invalid data.")
}

run.experiment.Rhbc1<- function(x, data) {
  if (x$stage == 0L) return(x)
  if (class(data) != 'data.frame') {
    stop('Require data to be data.frame')
  }
  if (x$stage == 2L) {
    return( stage.two(x, data[ ,1], data[ ,2]) )
  }
  is.stage.done <- getS3method('is.done', 'Rbc1')
  for (m in 1:x$nmax) {
    x <- stage.one(x, data[m, ])
    if(is.stage.done(x)) break
  }
  x
}

summary.Rhbc1 <- function(x) {
  cat('Procedure Rhbc1, see Buzaianu and Chen (2009)\n')
  isStageIComplete <- is.done(x, stage=1L)
  if (isStageIComplete) {
    cat('Stage I completed: ')
    if (x$decision == 0) {
      cat('select Control.\n')
    } else {
      cat('select Treatement', x$decision, '\n')
    }
  } else {
    cat('Stage I not completed\n')
  }
  nobs2 <- rep(NA, x$k+1L)
  y2 <- rep(NA, x$k+1L)
  if (x$stage == 2L) {
    cat('Ready to start Stage II.\n')
  } else if (x$stage == 0L) {
    nobs2[1] <- x$nobs2[1]
    nobs2[x$candidate+1] <- x$nobs2[2]
    y2[1] <- x$y2[1]
    y2[x$candidate+1] <- x$y2[2]
    cat('Stage II complted, select: ')
    if (x$decision == 0L) {
      cat('Control\n')
    } else {
      cat('Treatement', x$decision, '\n')
    }
  }

  z <- as.table( rbind(x$nobs, x$y, nobs2, y2) )
  rownames(z) <- c('No. Obs(I)', 'No. Successes(I)', 'No. Obs(II)', 'No. Successes(II)')
  colnames(z) <- c('Control', paste('Treatment', seq(1, x$k)))
  z
}
