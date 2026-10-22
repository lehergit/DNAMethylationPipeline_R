# test_pipeline.R
# Run from the project folder with:  Rscript tests/test_pipeline.R
# Base R only, nothing to install.

source("R/pipeline.R")
source("R/simulate.R")
source("legacy/legacy_analysis.R")

passed <- 0
failed <- 0
check <- function(description, condition) {
  ok <- isTRUE(condition)
  if (ok) passed <<- passed + 1 else failed <<- failed + 1
  cat(if (ok) "  PASS  " else "  FAIL  ", description, "\n", sep = "")
}
section <- function(title) cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
fails_with_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

# ---------------------------------------------------------------------------
section("beta <-> M-value conversion")
check("beta 0.5 is M-value 0", beta_to_m(0.5) == 0)
check("beta 0.75 is M-value log2(3)", abs(beta_to_m(0.75) - log2(3)) < 1e-12)
check("M-value 0 is beta 0.5", m_to_beta(0) == 0.5)
check("round trip recovers the original betas",
      all(abs(m_to_beta(beta_to_m(c(0.05, 0.3, 0.5, 0.8, 0.95))) - c(0.05, 0.3, 0.5, 0.8, 0.95)) < 1e-9))
check("exact 0 and 1 stay finite (no -Inf / Inf)", all(is.finite(beta_to_m(c(0, 1)))))
check("matrix shape is preserved", all(dim(beta_to_m(matrix(0.5, 3, 4))) == c(3, 4)))
check("bigger beta gives bigger M-value", beta_to_m(0.9) > beta_to_m(0.2))

# ---------------------------------------------------------------------------
section("simulate_methylation(): the fake data is well-formed")
sim <- simulate_methylation(seed = 42)
check("2000 probes x 20 samples", all(dim(sim$beta) == c(2000, 20)))
check("all beta values are strictly between 0 and 1", all(sim$beta > 0 & sim$beta < 1))
check("no missing values", !anyNA(sim$beta))
check("10 normal and 10 tumor samples", all(table(sim$group) == c(10, 10)))
check("exactly 100 planted DMPs", length(sim$true_dmps) == 100)
check("planted DMPs are real probe names", all(sim$true_dmps %in% rownames(sim$beta)))
check("annotation has one row per probe", nrow(sim$annotation) == 2000)
check("same seed gives identical data", identical(sim$beta, simulate_methylation(seed = 42)$beta))
check("different seed gives different data", !identical(sim$beta, simulate_methylation(seed = 43)$beta))
check("n_dmp = 0 plants nothing", length(simulate_methylation(n_dmp = 0)$true_dmps) == 0)

# ---------------------------------------------------------------------------
section("filter_probes(): quality control")
toy <- rbind(
  good      = c(0.2, 0.3, 0.25, 0.7, 0.8, 0.75),
  too_many_na = c(NA, NA, NA, 0.7, 0.8, 0.75),
  flat      = rep(0.5, 6),
  one_na_ok = c(0.2, NA, 0.25, 0.7, 0.8, 0.75)
)
qc <- filter_probes(toy, max_na_fraction = 0.2)
check("keeps the good probe", "good" %in% rownames(qc$beta))
check("keeps a probe with one missing value (1/6 <= 0.2)", "one_na_ok" %in% rownames(qc$beta))
check("removes the probe with too many NAs", "too_many_na" %in% qc$removed)
check("removes the constant probe", "flat" %in% qc$removed)
check("reports exactly which were removed", setequal(qc$removed, c("too_many_na", "flat")))
check("clean simulated data loses nothing", length(filter_probes(sim$beta)$removed) == 0)

# ---------------------------------------------------------------------------
section("test_differential(): per-probe statistics on tiny, checkable data")
tiny <- rbind(
  clear_gain = c(0.10, 0.12, 0.11, 0.80, 0.78, 0.82),
  clear_loss = c(0.85, 0.80, 0.83, 0.20, 0.22, 0.18),
  no_change  = c(0.50, 0.55, 0.45, 0.52, 0.48, 0.50)
)
colnames(tiny) <- paste0("s", 1:6)
grp <- c("normal", "normal", "normal", "tumor", "tumor", "tumor")
tr <- test_differential(tiny, grp, control = "normal")
check("one result row per probe", nrow(tr) == 3)
check("gain has positive delta beta (~ +0.70)", abs(tr$delta_beta[1] - 0.70) < 0.01)
check("loss has negative delta beta (~ -0.64)", tr$delta_beta[2] < -0.6)
check("no-change probe has delta beta near 0", abs(tr$delta_beta[3]) < 0.05)
check("clear changes have tiny p-values", all(tr$p_value[1:2] < 0.001))
check("no-change probe has a large p-value", tr$p_value[3] > 0.3)
check("FDR is never smaller than the raw p-value", all(tr$fdr >= tr$p_value - 1e-15))
check("all p-values are valid probabilities", all(tr$p_value >= 0 & tr$p_value <= 1))
check("swapping which group is control flips the sign",
      abs(test_differential(tiny, grp, control = "tumor")$delta_beta[1] + tr$delta_beta[1]) < 1e-12)

check("rejects a group vector of the wrong length",
      fails_with_error(test_differential(tiny, c("a", "b", "a"))))
check("rejects three groups",
      fails_with_error(test_differential(tiny, c("a", "a", "b", "b", "c", "c"))))
check("rejects a control that isn't a group",
      fails_with_error(test_differential(tiny, grp, control = "healthy")))

# ---------------------------------------------------------------------------
section("call_dmps(): needs BOTH significance and effect size")
fake <- data.frame(probe = c("a", "b", "c", "d"),
                   delta_beta = c(0.30, 0.02, -0.25, 0.30),
                   fdr = c(0.001, 0.001, 0.001, 0.20))
called <- call_dmps(fake)
check("big change + low FDR is a hit", called$is_dmp[1])
check("tiny change is NOT a hit even if significant", !called$is_dmp[2])
check("big negative change + low FDR is a hit", called$is_dmp[3])
check("big change with high FDR is NOT a hit", !called$is_dmp[4])

# ---------------------------------------------------------------------------
section("full pipeline vs. the planted truth")
fit <- run_pipeline(sim$beta, sim$group, sim$annotation, control = "normal")
hits <- fit$results$probe[fit$results$is_dmp]
found <- intersect(hits, sim$true_dmps)
sensitivity <- length(found) / length(sim$true_dmps)
precision <- length(found) / length(hits)
check(sprintf("sensitivity at least 0.85 (got %.2f)", sensitivity), sensitivity >= 0.85)
check(sprintf("precision at least 0.90 (got %.2f)", precision), precision >= 0.90)
check("hits are mostly true: false discovery rate under 10%", (length(hits) - length(found)) / length(hits) < 0.10)
check("the 20 strongest probes are all truly different",
      all(head(fit$results$probe[order(fit$results$p_value)], 20) %in% sim$true_dmps))

# ---------------------------------------------------------------------------
section("null simulation: with NO real differences, we should find ~nothing")
null_calls <- integer(20); raw_rate <- numeric(20)
for (s in 1:20) {
  nul <- simulate_methylation(n_dmp = 0, seed = 1000 + s)
  r <- run_pipeline(nul$beta, nul$group, control = "normal")$results
  null_calls[s] <- sum(r$is_dmp)
  raw_rate[s] <- mean(r$p_value < 0.05)
}
check(sprintf("raw p < 0.05 for about 5%% of probes (got %.3f)", mean(raw_rate)),
      mean(raw_rate) > 0.04 && mean(raw_rate) < 0.06)
check(sprintf("average false hits per null run is under 2 (got %.2f)", mean(null_calls)), mean(null_calls) < 2)
check("no null run produces more than 10 false hits", max(null_calls) <= 10)

# ---------------------------------------------------------------------------
section("power simulation: bigger effects are easier to detect")
sens_at <- function(shift) {
  s <- simulate_methylation(m_shift = shift, seed = 7)
  r <- run_pipeline(s$beta, s$group, control = "normal")$results
  length(intersect(r$probe[r$is_dmp], s$true_dmps)) / length(s$true_dmps)
}
sens_small <- sens_at(0.5); sens_mid <- sens_at(1.0); sens_big <- sens_at(2.0)
check(sprintf("sensitivity rises with effect size (%.2f < %.2f <= %.2f)", sens_small, sens_mid, sens_big),
      sens_small < sens_mid && sens_mid <= sens_big)
check("a tiny effect (0.5) is mostly missed", sens_small < 0.3)

# ---------------------------------------------------------------------------
section("enrichment: finds a planted CpG-island signal, ignores noise")
check("planted island bias gives odds ratio above 2", fit$enrichment$odds_ratio > 2)
check("planted island bias is statistically significant", fit$enrichment$p_value < 0.001)
check("contingency table counts all probes", sum(fit$enrichment$table) == 2000)

no_bias_p <- sapply(1:30, function(s) {
  d <- simulate_methylation(island_bias = 1, seed = 500 + s)
  enrichment_test(run_pipeline(d$beta, d$group, control = "normal")$results, d$annotation)$p_value
})
check(sprintf("with no planted bias, p < 0.05 only rarely (%d of 30)", sum(no_bias_p < 0.05)),
      sum(no_bias_p < 0.05) <= 4)

# ---------------------------------------------------------------------------
section("refactor check: new pipeline matches the legacy script exactly")
old <- legacy_analysis(sim$beta, sim$group)
new <- fit$results
check("same number of probes", nrow(old) == nrow(new))
check("same probe order", identical(old$probe, new$probe))
check("delta beta identical (tolerance 1e-12)", max(abs(old$dB - new$delta_beta)) < 1e-12)
check("p-values identical (tolerance 1e-12)", max(abs(old$p - new$p_value)) < 1e-12)
check("adjusted p-values identical (tolerance 1e-12)", max(abs(old$q - new$fdr)) < 1e-12)
check("hit calls identical", identical(old$hit, new$is_dmp))

# ---------------------------------------------------------------------------
cat(sprintf("\nResult: %d passed, %d failed\n", passed, failed))
cat(sprintf("Planted-truth run: sensitivity %.2f, precision %.2f\n", sensitivity, precision))
if (failed > 0) quit(status = 1)
