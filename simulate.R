# simulate.R
# Builds a fake methylation dataset where WE decide which probes truly differ.
# That lets us check the pipeline against a known answer. Nothing here comes
# from real patients or any lab.
source("R/pipeline.R")

simulate_methylation <- function(n_probes = 2000, n_per_group = 10, n_dmp = 100,
                                 m_shift = 2, noise_sd = 0.6,
                                 island_bias = 4, seed = 42) {
  set.seed(seed)
  probes <- sprintf("cg%06d", seq_len(n_probes))

  # Where each probe sits relative to a CpG island, and which gene it maps to.
  island <- sample(c("Island", "Shore", "Shelf", "OpenSea"), n_probes,
                   replace = TRUE, prob = c(0.25, 0.25, 0.15, 0.35))
  gene <- sprintf("GENE%03d", sample(500, n_probes, replace = TRUE))

  # Each probe's baseline methylation level (in M-value space).
  baseline_m <- beta_to_m(rbeta(n_probes, 2, 2))

  # Choose the truly different probes. island_bias > 1 makes Island probes
  # more likely to be chosen; island_bias = 1 means no preference.
  weights <- ifelse(island == "Island", island_bias, 1)
  dmp_idx <- if (n_dmp > 0) sort(sample(n_probes, n_dmp, prob = weights)) else integer(0)
  direction <- sample(c(-1, 1), n_probes, replace = TRUE)
  shift <- numeric(n_probes)
  shift[dmp_idx] <- m_shift * direction[dmp_idx]

  group <- rep(c("normal", "tumor"), each = n_per_group)
  n_samples <- length(group)
  m <- matrix(rnorm(n_probes * n_samples, 0, noise_sd), nrow = n_probes) + baseline_m
  m[, group == "tumor"] <- m[, group == "tumor"] + shift

  beta <- m_to_beta(m)
  dimnames(beta) <- list(probes, sprintf("sample%02d", seq_len(n_samples)))

  list(
    beta = beta,
    group = group,
    annotation = data.frame(probe = probes, gene = gene, island_relation = island,
                            stringsAsFactors = FALSE),
    true_dmps = probes[dmp_idx]
  )
}
