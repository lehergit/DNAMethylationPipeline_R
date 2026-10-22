# legacy_analysis.R
# A deliberately messy, "grew over time" version of the same analysis.
# It stands in for the kind of script a lab accumulates: magic numbers,
# copy-pasted blocks, one-letter names, no checks, no comments.
# Written for this demo; it is NOT code from any real lab.
# The refactored pipeline must produce the same numbers as this.

legacy_analysis <- function(d, g) {
  x <- d
  x[x < 0.0001] <- 0.0001
  x[x > 0.9999] <- 0.9999
  m <- log2(x/(1-x))
  r <- c()
  for (i in 1:nrow(m)) {
    a <- m[i, g == "tumor"]
    b <- m[i, g == "normal"]
    t <- t.test(a, b)
    r <- rbind(r, c(mean(d[i, g == "tumor"]) - mean(d[i, g == "normal"]), t$p.value))
  }
  q <- p.adjust(r[, 2], method = "BH")
  out <- data.frame(probe = rownames(d), dB = r[, 1], p = r[, 2], q = q)
  out$hit <- out$q < 0.05 & abs(out$dB) >= 0.1
  out
}
