# Evaluate an R distribution function on the rows of a CSV file of hexadecimal floats.
#   Rscript eval.R <function> <in.csv> <out.txt>
# Upper tail on the log scale (lower.tail = FALSE, log.p = TRUE), as in the audit spec.
# Input and output are C99 hexadecimal floats ("%a"), so every double passes exactly.
args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
d <- read.csv(args[2], header = FALSE, colClasses = "character")
X <- unname(lapply(d, as.numeric))
f <- switch(fn,
  pt     = function(x, df)     pt(x, df, lower.tail = FALSE, log.p = TRUE),
  pf     = function(x, d1, d2) pf(x, d1, d2, lower.tail = FALSE, log.p = TRUE),
  pchisq = function(x, df)     pchisq(x, df, lower.tail = FALSE, log.p = TRUE),
  pbeta  = function(x, a, b)   pbeta(x, a, b, lower.tail = FALSE, log.p = TRUE),
  pgamma = function(x, a)      pgamma(x, a, lower.tail = FALSE, log.p = TRUE),
  pbinom = function(k, n, p)   pbinom(k, n, p, lower.tail = FALSE, log.p = TRUE),
  qt     = function(lp, df)    qt(lp, df, lower.tail = FALSE, log.p = TRUE),
  qbeta  = function(lp, a, b)  qbeta(lp, a, b, lower.tail = FALSE, log.p = TRUE),
  stop("unknown function"))
y <- suppressWarnings(do.call(f, X))
writeLines(sprintf("%a", y), args[3])
