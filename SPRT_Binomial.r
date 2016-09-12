#
# SPRT_Binomial.r
#
# Sequential probability ratio test.
# 3/11/15, WAH@QD
# 3/13/15: Draft
#
# REFERENCE
# A Conditional Sequential Test for the Equality of Two Binomial Proportions
# William Q. Meeker, Jr
# Journal of the Royal Statistical Society. Series C (Applied Statistics)
# Vol. 30, No. 2 (1981), pp. 109-115
#------------------------------------------------------------------------------#
binomial.SPRT <- function(x, y, t.1, alpha=0.05, beta=0.10, n.0) {
  #
  # Meeker's SPRT for matched `x` (treatment) and `y` (control), 
  # both indicator responses, likelihood ratio t.1, error rates alpha and beta,
  # and (optionally) truncation after trial n.0.
  #
  # The return variable contains these elements:
  # * outcome:   "continue," "reject null," or "accept null".
  # * index:     Index at which the outcome decision was made (or NA)
  # * limit.low: Lower test limit (determined by alpha and beta)
  # * limit.high:Upper test limit
  # * truncated: If the test was truncated, the value of `n.0`; NA otherwise
  # * x.1:       Original data `x`, cumulative
  # * r:         Cumulative sum of x+y
  # * stats:     Series of cumulative sums of log probability ratios
  # * limits:    Two rows giving lower and upper critical limits, respectively
  #
  # Input checks.
  #
  if (length(x) != length(y)) 
    warning("Treatment and control lengths differ; only earliest data are used.",
            immediate. = TRUE)
  if (t.1 <= 1)
    stop("The odds ratio should exceed 1.")
  if (alpha > 0.5 || beta > 0.5) 
    warning("Unrealistic values of alpha or beta were supplied.")
  #
  # Input cleaning.
  #
  if (!missing(n.0)) n.0 <- floor(n.0)
  #
  # Auxiliary functions.
  #
  g <- function(x, r, n, t.1, t.0=1) {
    #
    # Meeker's (1981) function `g`, the log probability ratio.
    #
    -log(h(x, r, n, t.1)) + log(h(x, r, n, t.0))
  }
  h <- function(x, r, n, t=1) {
    #
    # Reciprocal of Meeker's (1981) function `h`: the conditional probability of 
    # `x` given `r` and `n`, when the odds ratio is `t`.
    #
    # `x` is his "x1", the number of positives in `n` control trials.
    # `r` is the total number of positives.
    # `n` is the number of (control, treatment) pairs.
    # `t` is the odds ratio.
    #
    f(r, n, t, offset=f.term.log(x, r, n, t))
  }
  f <- function(r, n, t, offset=0) {
    #
    # Meeker's (1981) function exp(F(r,n,t)), proportional to the probability of 
    #  `r` (=x1+x2) in `n` paired trials with an odds ratio of `t`.
    #
    # This function does *not* vectorize over its arguments.
    #
    if (any(is.na(c(r,n,t)))) return (NA)
    sum(f.term(max(0, r-n):min(n,r), r, n, t, offset))
  }
  f.term <- function(j, r, n, t, offset=0) exp(f.term.log(j, r, n, t, offset))
  f.term.log <- function(j, r, n, t, offset=0) {
    #
    # Up to an additive constant, the log probability that (x1, x1+x2) = (j, r) 
    # in `n` paired trials with odds ratio of `t`.
    #
    # `offset` is used to adjust the result to avoid under/overflow.
    #
    lchoose(n,j) + lchoose(n,r-j) + j*log(t) - offset
  }
  log.f <- function(r, n, t, offset=0) {
    #
    # A protected vesion of log(f), Meeker's function `F`.
    #
    z <- f(r, n, t, offset); ifelse(z > 0, log(z), NA)
  }
  c.lower.upper <- function(r, n, t1, t0=1, alpha=0.05, beta=0.10) {
    #
    # Meeker's (1981) functions c_L(r,n) and c_U(r,n), the  critical values for x1.
    # 0 <= r <= 2n; t1 >= t0 > 0.
    #
    offset <- f.term.log(ceiling(r/2), r, n, t1)         # Estimate scale of terms
    z <- log.f(r, n, t1, log.f(r, n, t0, offset)+offset) # Common term
    a <- -log(alpha/(1-beta))
    b <- log(beta/(1-alpha))
    (c(lower=b, upper=1+a) + z) / log(t1 / t0)           # Eventually, take the floor
  }
  #
  # The calculations.
  #
  l <- log(beta / (1-alpha))
  u <- -log(alpha / (1-beta))
  n <- 1:min(length(x), length(y))
  if (!missing(n.0)) n <- n[n <= n.0]
  x.1 <- cumsum(x[n])
  r <- x.1 + cumsum(y[n])
  stats <- mapply(g, x=x.1, r=r, n=n, t.1=t.1)
  limits <- floor(mapply(c.lower.upper, r, n, 
                         MoreArgs=list(t1=t.1, alpha=alpha, beta=beta)))
  #
  # Perform the test by finding the first index, if any, at which `stats`
  # falls outside the open interval (l, u).
  #
  k <- which(stats >= u | stats <= l)
  if (length(k) < 1) {
    k <- NA
    outcome <- "continue"
  } else {
    k <- min(k)
    outcome <- ifelse(stats[k] >= u, "reject null", "accept null")
  }
  if (!missing(n.0) && is.na(k)) {
    #
    # Truncate at trial n.0, using Meeker's H0-conservative formula (2.2).
    # Leave k=NA to indicate the decision was made due to truncation.
    #
    c.l <- c.lower.upper(r, n.0, t.1, alpha, beta)
    c.l <- floor(mean(c.l) - 1/2)
    outcome <- ifelse(x.1[n.0] <= c.l, "accept null", "reject null")
    truncated <- n.0
  } else {truncated <- NA}
  return (list(outcome=outcome, index=k, limit.low=l, limit.high=u,
               truncated=truncated, x.1=x.1, r=r, stats=stats, limits=limits))
} # binomial.SPRT
#------------------------------------------------------------------------------#
#
# Plot the output of binomial.SPRT using the log probability ratios.
#
plot.SPRT <- function(stats, use.points=FALSE, col="#c02020", # removed add=FALSE
                      main, sub, ...) {
  #
  # Plot the binomial SPRT statistics, showing the cumulative log 
  # probability ratio and the uper and lower limits.
  #
  # To show individual values, set `use.points=TRUE`.  (This is automatically
  # done with fewer than 100 points.)
  #
  # Set `add=TRUE` to overwrite the existing plot (useful for simulations).
  #
  # `stats` is the object returned by `binomial.SPRT`.
  # The optional parameters are passed to the first `plot` call; the best
  # use is main=<text> to display a title.
  #
  y <- stats$stats
  n <- length(y)
  n.stop <- ifelse(is.na(stats$truncated), n+1/2, stats$truncated-1/2)
  l <- stats$limit.low; u <- stats$limit.high
  k <- stats$index
  suppress <- "#000000" # Color of irrelevant data (was #00000018)
  
  if (!add) {
    if (missing(main)) main="Sequential Binomial Proportion Test"
    if (missing(sub)) sub=paste("The decision is to", stats$outcome)
    plot(c(0, n), range(c(y, l, u)), type="n", bty="n",
       xlab="Trial", ylab="Log probability ratio", main=main, sub=sub, ...)
    rect(1/2, l, n.stop, u, col=gray(0.97), border=NA) # Indifference region
    abline(h = c(l, u, 0), lty=c(3, 3, 1), lwd=2, col=c("Black", "Black", "Gray"))
  }
  
  if (add) {
    lines(y, lwd=1, col=suppress)
  } else {
    lines(y, lwd=3, col=suppress)
  }
  
  if (n < 100 || use.points) points(y, pch=16, col=suppress)
  
  if(is.na(k)) k <- length(y)
  #abline(v=k, lwd=2, lty=3, col="Gray")
  if (add) {
    lines(y[1:k], lwd=1, col=col)
  } else {
    lines(y[1:k], lwd=3, col=col)
  }
  
  if (n < 100 || use.points) points(y[1:k], pch=16, col=col)
  
  if (!add) {
    text(1, l, "No difference", adj=c(0, -1/2))
    text(1, u, "Significant difference", adj=c(0, 1.25))
    #text(1, c(l, u), c("No difference", "Significant difference"), pos=c(3,1))
  }
}# plot.SPRT
#
# A simple plot in Meeker's style: cumulative sums of treatment successes,
# bounded by the critical limits.
#
# The first trial (if any) at which the cusum reaches either limit determines
# the decision, *regardless* of why might happen in successive trials.
#
plot.SPRT.I <- function(s, main, ...) {
  #
  # A "Method I" plot showing the cumulative treatment results and the 
  # critical region.
  #
  # `s`    is the output of binomial.SPRT.
  # `main` is the plot title; it has a reasonable default.
  # `...`  are optional arguments passed to `plot`.
  #
  n <- 1:length(s$x.1)
  c.L0 <- ifelse(s$limits[1, ] < 0, NA, s$limits[1, ])
  colors <- c("Black", "#c02020", "#2020c0")
  pch <- c(16, 1, 1)
  if (missing(main)) main="Sequential Binomial Proportion Test"
  plot(range(n), range(c(s$x.1,c.L0,s$limits[2, ]), na.rm=TRUE)+c(-1,1)/2, 
       type="n", bty="n", main=main, 
       xlab="Trial", ylab="Cumulative successes", ...)
  polygon(x=c(n, rev(n)), y=c(pmax(s$limits[1, ],0), rev(s$limits[2, ])), 
          border=NA, col=gray(0.97))
  lines(c.L0, col=colors[2])
  lines(s$limits[2, ], col=colors[3])
  lines(s$x.1, col=colors[1], lty=3)
  if (max(n) < 300) {
    points(c.L0, col=colors[2], pch=pch[2])
    points(s$limits[2, ], col=colors[3], pch=pch[3])
    points(s$x.1, col=colors[1], pch=pch[1])
  }
  legend("topleft", c("Test", "CL", "CU"), bty="n",
         col=colors, lty=c(3,1,1), lwd=2)
}# plot.SPRT.I
