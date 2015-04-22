# Experimental
# Multistage bc procedure exploration

RMbc1.state(k, ns, s, c, scale=1L) %when% {
  k > 1L
  ns >= 1L
  s > 1L
  c >= 0L
} %as% {
  structure(list(
    k=k,
    smax=s,
    ns=ns,
    decision=NULL,
    stage=0L,
    scale=scale
  ), class="RMbc1")
}

RMbc1.state(k, ns, s, c, scale=1L) %as% { stop('Invalid parameters in constructing RMbc1.') }

RMbc2.state(k, ns, s, c, scale=1L) %when% {
  k > 1L
  ns >= 1L
  s > 1L
  c >= 0L
} %as% {
  structure(list(
    k=k,
    smax=s,
    c=c,
    ns=as.integer(ns),
    terminated=rep(FALSE, k+1L),
    nobs=rep(0L, k+1L),
    success=rep(0L, k+1L),
    decision=NULL,
    scale=scale,
    stage=0L
  ), class="RMbc2")
}

RMbc2.state(k, ns, s, c, scale=1L) %as% { stop('Invalid parameters in constructing RMbc2.') }


RMbc2(k, ns, s, c, scale=1L) %as% {

  state <- RMbc2.state(k, ns, s, c, scale=1L)
  is.done <- function() {
    state$stage >= state$smax ||
    all(state$terminated)
  }
  get.state <- function() {
    state
  }
  next.state <- function(data) {
    if( is.done() ) return(state)
    # skip invalid data: missing values or abnormal number of success
    if( any(is.na(data[!state$terminated]) | data[!state$terminated] > state$ns) ) return(state)
    # update state
    state$stage <<- state$stage + 1L
    state$nobs <<- state$nobs + state$ns * (!state$terminated)
    state$success <<- state$success + data * (!state$terminated)
    # stoping rules
    jmax <- which.is.max(state$success[-1]) + 1
    max.success <- state$success[jmax]
    if(max.success <= state$success[1] + state$c - state$scale * (state$ns * state$smax - state$ns * state$stage) ) {
      state$terminated <<- rep(TRUE, state$k+1L)
      state$decision <<- 0L
    } else {
      state$terminated <<- (
        state$terminated |
        (max.success - state$success > c(state$c, rep(0, state$k)) + state$scale * (state$ns * state$smax - state$nobs))
      )
      if(all(state$success[!state$terminated] == max.success)) {
        state$terminated <<- rep(TRUE, state$k+1L)
        state$decision <<- jmax - 1
      }
    }
    state
  }
  list(next.state=next.state, is.done=is.done, get.state=get.state)
}
