#
# ABTest.r
#
# Read SurfCanyon data and apply the binomial SPRT to them.
# 3/11/2015, WAH@QD
# 3/13/2015: Draft
#
setwd("~/Rank Dynamics/Development/R Workspace/SPRT") # Set working directory
source("SPRT_Binomial.r") # exports binomial.SPRT(), plot.SPRT(), plot.SPRT.I()
#------------------------------------------------------------------------------#
#
# Impute binary streams associated with daily summaries.
# Example: 
#    trials <- c(5, 3, 3)
#    successes <- c(2, 0, 3)
#    stream(trials,successes, 1/2)
#
# Output is "[1] 1 0 1 0 0 0 0 0 1 1 1", showing a binary array of 5+3+3 values
# of which 2 of the first 5 are ones, 0 of the next 3 are ones, and all 3 of 
# the last 3 are ones.
#
# Omit `start` to position the ones randomly within each group.  Otherwise,
# the fractional value of `start` establishes the origin of an equally-spaced
# group of ones within each group.
#
stream <- function(trials, successes, start) {
  # Strip NA values
  successes <- successes[!is.na(successes)]
  trials <- trials[(length(trials) - length(successes) + 1):length(trials)]
  #
  # Inputs are summaries (total counts of trials and successes).
  # Output is a single binary array.
  #
  if (missing(start)) {
    #
    # Position the successes randomly.
    #
    f <- function(x, y, start) {
      sample(c(rep(0, x - y), rep(1, y)), x)
    }
  } else {
    #
    # Position the successes at regular intervals (Bresenham's algorithm).
    #
    f <- function(x, y, start) {
      #
      # x is a trial count and y is a success count. 
      # Presumably x >= y >= 0 are integers.
      #
      if (x == 0) return (integer(0))
      diff(floor((0:x + start*x) * (y / x)))
    }
  }
  as.vector(unlist(mapply(f, trials, successes, runif(length(trials)))))
}
#------------------------------------------------------------------------------#
#
# Read and pre-process the data:
# CSV form of [ab test] sheet, sorted by increasing time.
#
#fn.data <- "EBay_ABTest.csv" 
fn.data <- "ebay.csv"
#fn.data <- "ebay_interleave_sales.csv"
x.raw <- read.csv(fn.data)
seed <- round(runif(1)*1000)
set.seed(359) # Allows reproducible results
#control <- stream(x.raw$ControlUsers, x.raw$ControlSales)
#treatment <- stream(x.raw$TestUsers, x.raw$TestSales)
control <- stream(x.raw$Users, x.raw$Action) # for ebay.csv
treatment <- stream(x.raw$Users2, x.raw$Action2) # for ebay.csv
#control <- stream(x.raw$Clicks, x.raw$Static)
#treatment <- stream(x.raw$Clicks, x.raw$Dynamic)

#
# Display the data.
# Inputs are parallel binary vectors of treatment and control outcomes.
# The most likely modifications to make are to `xlab`, `ylab`, and `main`
# to annotate the axes and the plot, respectively.
#
plot.ABTest <- function(treatment, control,
                        col.control="#2020c0", col.treatment="#c02020") {
  suppress <- "#0000001a"
  s.control <- cumsum(control)
  s.treatment <- cumsum(treatment)
  d.control <- 1:length(control)
  d.treatment <- 1:length(treatment)
  n <- min(length(treatment), length(control))
  
  plot(range(c(d.control, d.treatment)), 
       range(c(s.control, s.treatment)), type="n", 
       xlab="Trial", 
       ylab="Cumulative sales events",
       main="Data Summary")
  lines(d.control, s.control, col=suppress, lwd=2)
  lines(d.treatment, s.treatment, col=suppress, lwd=2)
  lines(d.control[1:n], s.control[1:n], col=col.control, lwd=2)
  lines(d.treatment[1:n], s.treatment[1:n], col=col.treatment, lwd=2)
  if (n < 300) {
    points(d.control, s.control, col=suppress, pch=1)
    points(d.treatment, s.treatment, col=suppress, pch=16)
    points(d.control[1:n], s.control[1:n], col=col.control, pch=1)
    points(d.treatment[1:n], s.treatment[1:n], col=col.treatment, pch=16)
  }
  legend("topleft", c("Test", "Control"), bty="n",
         col=c(col.treatment, col.control), lwd=2)
} # plot.ABTest
#plot.ABTest(treatment, control)
#------------------------------------------------------------------------------#
#
# Compute the test statistics.
# The calculation will gradually slow as the input length increases.
# Vectors of length 20,000 will take just a few seconds.
#
alpha <- 0.05  # Acceptable Type I ("false positive") error rate
beta <- 0.05   # Acceptable Type II ("false negative") error rate
q <- 1.35      # Relative odds ratio to detect
stats <- binomial.SPRT(treatment, control, q, alpha, beta)
add <- FALSE

#
# Display the test results.
# Two different methods are available.
#
s <- paste0("Alpha=", round(alpha,2), "; Beta=", round(beta,2),
           "; Target Odds Ratio=", round(q, 2))
#for (f in c(plot.SPRT, plot.SPRT.I)) {
for (f in c(plot.SPRT)) { # Only make one plot
  f(stats, sub=s, cex.sub=0.8)
}

# These lines are for simulation
# sim.seed <- round(runif(1)*1000)
# set.seed(sim.seed) # Allows reproducible results
# add <- TRUE
# sim.continue <- 0
# sim.reject <- 0
# sim.accept <- 0
# for (sim.count in 1:100) {
#   control.successrate <- 3 / 1000
#   treatment.successrate <- control.successrate * 1 # Set to odds ratio
# #  control <- (runif(sum(x.raw$ControlUsers))<control.successrate)*1
# #  treatment <- (runif(sum(x.raw$TestUsers))<treatment.successrate)*1
#   control <- (runif(25000)<control.successrate)*1
#   treatment <- (runif(25000)<treatment.successrate)*1
#   stats <- binomial.SPRT(treatment, control, q, alpha, beta)
#   if (stats$outcome == "continue") { sim.continue <- sim.continue + 1 }
#   if (stats$outcome == "reject null") { sim.reject <- sim.reject + 1 }
#   if (stats$outcome == "accept null") { sim.accept <- sim.accept + 1 }
#   s <- paste0("Alpha=", round(alpha,2), "; Beta=", round(beta,2),
#               "; Target Odds Ratio=", round(q, 2))
#   for (f in c(plot.SPRT)) {
#     f(stats, sub=s, cex.sub=0.8)
#   }
# }

# Compute and plot rolling averages
#rolling.number <- 3000
#treatment.rolling <- character(0)
#for (i in rolling.number:length(treatment)) {
#  treatment.rolling <- c(treatment.rolling, sum(treatment[(i-rolling.number):i]))
#}
#x <- seq(rolling.number,length(treatment))
#plot(treatment.rolling, col="red", type="l", ylim=c(0,max(as.numeric(treatment.rolling))),
#     xlab="Users", ylab="Rolling Average Sales")
#control.rolling <- character(0)
#for (i in rolling.number:length(control)) {
#  control.rolling <- c(control.rolling, sum(control[(i-rolling.number):i]))
#}
#lines(control.rolling, col="blue", type="l")

