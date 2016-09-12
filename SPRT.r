#
# Classical SPRT
# See https://en.wikipedia.org/wiki/Sequential_probability_ratio_test
# 9 October 2015, WAH @QD
#
setwd("~/Rank Dynamics/Development/R Workspace/SPRT") # Set working directory

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
x.df <- read.csv("ebay_interleave.csv")
riffle <- function(x, y) as.vector(rbind(x, y))
x <- riffle(-x.df$Control, x.df$Test) # For each record, deduct all control events before processing test events
#
# Run the test at an odds ratio of 1.1.
#
odds <- 1.1
(p.A <- 1/ ( 1 + 1/odds ))
result <- SPRT(x, p.A = p.A)
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
plot(cumsum(abs(x)), cumsum(x), type="l", col="Gray",
     xlab="Trials", ylab="Tests - Controls",
     main=paste("SPRT: Conclusion Is To", action),
     sub=paste("H0: p = 1/2; HA: p =", format(p.A, digits=3)))
points(x.df$Trials, cumsum(x.df$Test - x.df$Control), pch=16, cex=0.7)
col <- ifelse(result$Result=="Reject", "Red", "Blue")
if (result$Result != "Continue") abline(v = result$Stop, col=col)
#==============================================================================#
# # 
# # Code for testing SPRT.
# #
# par(mfcol=c(2,5))
# set.seed(17)
# n <- 500
# p.A <- 0.55
# for (i in 1:25) {
#   #
#   # Example 1: generate a series according to the null.
#   #
#   x <- 2*(runif(n) < 1/2) - 1
#   result <- SPRT(x, p.A=p.A)
#   plot(cumsum(x), cex=1/2, main="Null", sub=result$Result)
#   col <- ifelse(result$Result=="Reject", "Red", "Blue")
#   if (result$Result != "Continue") abline(v = result$Stop, col=col)
#   #
#   # Example 2: generate a series according to the alternative.
#   #
#   x <- 2*(runif(n) < p.A) - 1
#   result <- SPRT(x, p.A=p.A)
#   plot(cumsum(x), cex=1/2, main="Alternate", sub=result$Result)
#   col <- ifelse(result$Result=="Reject", "Red", "Blue")
#   if (result$Result != "Continue") abline(v = result$Stop, col=col)
# }
