# pipeline.R
# The analysis steps, one small function each.
#
# The big idea: a methylation "beta value" is a probe's methylation level,
# a number from 0 (never methylated) to 1 (always methylated). We compare
# beta values between two groups (say, tumor vs. normal tissue) and ask
# which probes differ by more than chance would explain.

# ---- Step 1: move between beta values and M-values --------------------------
# Beta values are easy to interpret but statistically awkward, because they
# are squeezed between 0 and 1. M-values (the log2 odds of beta) behave much
# better in t-tests, so we test on M-values and report on beta values.
# `eps` keeps exact 0s and 1s from turning into infinity.
beta_to_m <- function(beta, eps = 1e-4) {
  beta <- pmin(pmax(beta, eps), 1 - eps)
  log2(beta / (1 - beta))
}

m_to_beta <- function(m) {
  2^m / (1 + 2^m)
}

# ---- Step 2: quality control ------------------------------------------------
# Drop probes that cannot be tested sensibly: too many missing values, or
# no variation at all. We return what was removed so nothing disappears
# silently.
filter_probes <- function(beta, max_na_fraction = 0.2) {
  na_fraction <- rowMeans(is.na(beta))
  spread <- apply(beta, 1, sd, na.rm = TRUE)
  keep <- na_fraction <= max_na_fraction & !is.na(spread) & spread > 0
  list(
    beta = beta[keep, , drop = FALSE],
    removed = rownames(beta)[!keep]
  )
}

# ---- Step 3: test every probe -----------------------------------------------
# For each probe: Welch t-test on M-values between the two groups, plus the
# change in average beta value (delta beta) so results stay interpretable.
# Then correct for testing thousands of probes at once (Benjamini-Hochberg,
# which controls the false discovery rate, or FDR).
test_differential <- function(beta, group, control = levels(factor(group))[1]) {
  group <- factor(group)
  stopifnot(
    "group must have one entry per sample (column)" = ncol(beta) == length(group),
    "group must have exactly two levels" = nlevels(group) == 2,
    "control must be one of the group levels" = control %in% levels(group)
  )
  is_case <- group != control
  m <- beta_to_m(beta)

  stats <- t(vapply(seq_len(nrow(m)), function(i) {
    tt <- t.test(m[i, is_case], m[i, !is_case])
    c(t = unname(tt$statistic), p = tt$p.value)
  }, numeric(2)))

  out <- data.frame(
    probe      = rownames(beta),
    delta_beta = rowMeans(beta[, is_case, drop = FALSE], na.rm = TRUE) -
                 rowMeans(beta[, !is_case, drop = FALSE], na.rm = TRUE),
    t_stat     = stats[, 1],
    p_value    = stats[, 2],
    stringsAsFactors = FALSE
  )
  out$fdr <- p.adjust(out$p_value, method = "BH")
  rownames(out) <- NULL
  out
}

# ---- Step 4: call the differentially methylated probes (DMPs) ---------------
# A probe counts as a hit only if it is BOTH statistically convincing (low FDR)
# AND big enough to matter (|delta beta| above a floor). Requiring both is
# standard practice, because a tiny difference can be "significant" in a big
# study without meaning much.
call_dmps <- function(results, fdr_cutoff = 0.05, min_delta_beta = 0.1) {
  results$is_dmp <- results$fdr < fdr_cutoff & abs(results$delta_beta) >= min_delta_beta
  results
}

# ---- Step 5: are the hits clustered in one kind of region? -------------------
# Fisher's exact test on a 2x2 table: (hit / not hit) x (in region / not).
# An odds ratio above 1 means hits are over-represented in that region.
enrichment_test <- function(results, annotation,
                            column = "island_relation", category = "Island") {
  merged <- merge(results[, c("probe", "is_dmp")],
                  annotation[, c("probe", column)], by = "probe")
  in_category <- merged[[column]] == category
  tab <- table(
    dmp         = factor(merged$is_dmp, levels = c(TRUE, FALSE)),
    in_category = factor(in_category, levels = c(TRUE, FALSE))
  )
  ft <- fisher.test(tab)
  list(table = tab, odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}

# ---- Step 6: a volcano plot -------------------------------------------------
# x = how big the change is, y = how convincing it is. Hits are colored.
volcano_plot <- function(results, main = "Differential methylation") {
  colors <- ifelse(results$is_dmp, "firebrick", "grey70")
  plot(results$delta_beta, -log10(results$p_value),
       pch = 19, cex = 0.6, col = colors,
       xlab = "Change in methylation (delta beta)",
       ylab = "-log10(p-value)", main = main)
  abline(v = c(-0.1, 0.1), lty = 2, col = "grey40")
}

# ---- Everything together ----------------------------------------------------
run_pipeline <- function(beta, group, annotation = NULL,
                         fdr_cutoff = 0.05, min_delta_beta = 0.1,
                         control = levels(factor(group))[1]) {
  qc <- filter_probes(beta)
  results <- test_differential(qc$beta, group, control = control)
  results <- call_dmps(results, fdr_cutoff, min_delta_beta)
  enrichment <- if (!is.null(annotation) && any(results$is_dmp)) {
    enrichment_test(results, annotation)
  } else NULL
  list(results = results, removed_probes = qc$removed, enrichment = enrichment)
}
