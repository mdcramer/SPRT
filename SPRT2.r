#
# SPRT2.r
# Classical SPRT
# See https://en.wikipedia.org/wiki/Sequential_probability_ratio_test
# and Kendall, Stuart, & Ord Vol. II Chapter 24.
# 9 October 2015, WAH @QD
# Updated 14 October 2015: fixed the position of the decision point in the
#                          plots and included the threshold lines in them.
setwd("~/Rank Dynamics/Development/R Workspace/SPRT") # Set working directory
#==============================================================================#
SPRT <- function(x, p.0=1/2, p.A=1/2, alpha=0.05, beta=0.05, LR=0) {
  #
  # `x`     is an array of counts: positive counts for successes, negative
  #         counts for failures.
  # `p.0`   is the null hypothesis (chance of success).
  # `p.A`   is the alternative hypothesis.
  # `alpha` is an upper bound on the false positive rate. 
  # `beta`  is an upper bound on the false negative rate.
  # `LR`    is the previous log likelihood ratio.
  #
  likelihood.log <- function(y, p) y * ifelse(y > 0, log(p), -log(1-p))
  lr <- cumsum(likelihood.log(x, p.A) - likelihood.log(x, p.0)) + LR
  a <- log(beta/(1-alpha))
  b <- log((1-beta)/alpha)
  i.stop <- which(lr <= a | lr >= b)
  if (length(i.stop) >= 1) {
    i.stop <- min(i.stop)
    result <- ifelse(lr[i.stop] <= a, "Accept", "Reject")
  } else {
    i.stop <- Inf
    result <- "Continue"
  }
  return(list(Stop=i.stop, Result=result))
}
#==============================================================================#
#
# Actual data.
#
x.df <- read.csv("Ebay_interleave.csv")
#x.df <- read.csv("shopwiki_interleave.csv")
riffle <- function(x, y) as.vector(rbind(x, y))
x <- riffle(-x.df$Control, x.df$Test) # For each record, deduct all control events before processing test events
truex <- c(read.csv("Ebay_interleave_stream.csv"))
#truex <- c(read.csv("shopwiki_interleave_stream.csv"))
truex <- truex$Stream
#
# Run the test at an odds ratio of 'odds'.
#
odds <- 1.1
alpha <- 0.05
beta <- 0.05
p.0 <- 1/2
(p.A <- 1/ ( 1 + 1/odds ))
a <- log(1-beta) - log(alpha)
b <- log(beta) - log(1-alpha)
s <- log(p.0) - log(p.A) + log(1-p.A) - log(1-p.0)
u <- log(1-p.0) - log(1-p.A)
result <- SPRT(x, p.0=p.0, p.A = p.A, alpha=alpha, beta=beta)
#
# Display the estimate of p.
#
n.C <- sum(x.df$Control)
n.T <- sum(x.df$Test)
cat("Estimated proportion is", format(n.T / (n.C + n.T), digits=3), "\n")
#
# Compute a p-value as if this were a designed experiment.
#
cat("p-value is", format(pbinom(n.T, n.T+n.C, 1/2, lower.tail=FALSE), digits=3), "\n")
#
# Plot everything.
#
par(mfrow=c(1,1))
action <- ifelse(result$Result=="Continue", "Sampling", "The Null Hypothesis")
action <- paste(result$Result, action)
n <- cumsum(abs(x))
y <- cumsum(x)
usr <- par("usr")
pin <- par("pin")
asp <- diff(usr[1:2]) / diff(usr[3:4]) * pin[2] / pin[1] # Aspect ratio
mu <- -2*u/s - 1                                         # Slope of lines
srt <- atan(mu * asp) * 180 / pi
plot(cumsum(abs(truex)), cumsum(truex), type="l", col="Gray",
     xlab="Trials = Clicks (1 dot = 1 day)", ylab="Cumulative Treatment Minus Control Clicks",
     ylim=c(-25, 200),
     main="SPRT: Significant User Preference for Dynamically Ranked Results")
mtext(bquote(paste(H[0], ": p=", .(p.0), "; ", H[A], ": p=", .(round(p.A, 3)),
                   " - Clicks beyond top 10 results for queries with 200+ results")))
x.0 <- crossprod(range(n), c(1,2)/3)
text(x.0, 2*a/s + mu*x.0, "Accept", pos=1, col="Blue", srt=srt)
x.0 <- crossprod(range(n), c(2,1)/3)
text(x.0, 2*b/s + mu*x.0, "Reject", pos=3, col="Red", srt=srt)
abline(c(a, -u - s/2) * 2/s, col="Blue", lwd=2, lty=3)
abline(c(b, -u - s/2) * 2/s, col="Red", lwd=2, lty=3)
points(cumsum(x.df$Trials), cumsum(x.df$Test - x.df$Control), pch=16, cex=0.7)
col <- ifelse(result$Result=="Reject", "Red", "Blue")
if (result$Result != "Continue") {
  abline(v = n[result$Stop], col=col) #horizontal line at accept/reject decision
  i <- cumsum(x.df$Trials) >=  n[result$Stop]
  points(cumsum(x.df$Trials)[i], cumsum(x.df$Test - x.df$Control)[i], pch=16, cex=0.8, col=col)
}
#==============================================================================#
# # 
# # Code for testing SPRT.
# #
# p.0 <- 0.50
# p.A <- 0.55
# alpha <- 0.05
# beta <- 0.05
# a <- log(1-beta) - log(alpha)
# b <- log(beta) - log(1-alpha)
# s <- log(p.0) - log(p.A) + log(1-p.A) - log(1-p.0)
# u <- log(1-p.0) - log(1-p.A)
# 
# par(mfcol=c(2,5))
# set.seed(17)
# n <- 500
# for (i in 1:5) {
#   #
#   # Example 1: generate a series according to the null.
#   #
#   x <- 2*(runif(n) < 1/2) - 1
#   result <- SPRT(x, p.A=p.A, alpha=alpha, beta=beta)
#   plot(cumsum(x), cex=1/2, main="Null", sub=result$Result)
#   abline(c(a, -u - s/2)* 2/s, col="Blue", lwd=2, lty=3)
#   abline(c(b, -u - s/2)* 2/s, col="Red", lwd=2, lty=3)
#   col <- ifelse(result$Result=="Reject", "Red", "Blue")
#   if (result$Result != "Continue") abline(v = result$Stop, col=col)
#   #
#   # Example 2: generate a series according to the alternative.
#   #
#   x <- 2*(runif(n) < p.A) - 1
#   result <- SPRT(x, p.A=p.A)
#   plot(cumsum(x), cex=1/2, main="Alternate", sub=result$Result)
#   abline(c(a, -u - s/2)* 2/s, col="Blue", lwd=2, lty=3)
#   abline(c(b, -u - s/2)* 2/s, col="Red", lwd=2, lty=3)
#   col <- ifelse(result$Result=="Reject", "Red", "Blue")
#   if (result$Result != "Continue") abline(v = result$Stop, col=col)
# }
