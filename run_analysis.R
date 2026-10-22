# run_analysis.R
# Run the whole analysis from start to finish:  Rscript run_analysis.R
# It simulates a dataset, runs the pipeline, prints a summary, and saves a
# results table and a volcano plot into output/.

source("R/pipeline.R")
source("R/simulate.R")

cat("Step 1: simulate 2,000 probes x 20 samples (10 normal, 10 tumor)\n")
sim <- simulate_methylation(n_probes = 2000, n_per_group = 10, n_dmp = 100, seed = 42)

cat("Step 2: run the pipeline\n")
fit <- run_pipeline(sim$beta, sim$group, sim$annotation, control = "normal")
res <- fit$results

cat("Step 3: compare against the truth we planted\n")
called <- res$probe[res$is_dmp]
true_hits <- sim$true_dmps
found <- intersect(called, true_hits)
cat(sprintf("  probes tested:            %d (removed by QC: %d)\n", nrow(res), length(fit$removed_probes)))
cat(sprintf("  truly different probes:   %d\n", length(true_hits)))
cat(sprintf("  probes the pipeline hit:  %d\n", length(called)))
cat(sprintf("  correct hits:             %d\n", length(found)))
cat(sprintf("  sensitivity (found / true): %.2f\n", length(found) / length(true_hits)))
cat(sprintf("  precision (correct / hits): %.2f\n", length(found) / length(called)))

if (!is.null(fit$enrichment)) {
  cat(sprintf("\nStep 4: CpG-island enrichment: odds ratio %.2f, p = %.2g\n",
              fit$enrichment$odds_ratio, fit$enrichment$p_value))
}

dir.create("output", showWarnings = FALSE)
write.csv(res[order(res$p_value), ], "output/dmp_results.csv", row.names = FALSE)
png("output/volcano.png", width = 900, height = 700, res = 120)
volcano_plot(res, main = "Tumor vs. normal (simulated data)")
dev.off()
cat("\nSaved output/dmp_results.csv and output/volcano.png\n")
