# Differential Methylation Pipeline

A small, tested R pipeline for finding differentially methylated probes (DMPs) between two groups, such as tumor vs. normal tissue. It takes a matrix of methylation beta values, tests every probe, calls the ones that differ, checks whether the hits cluster in a particular genomic region, and draws a volcano plot.

The project also includes a simulator that builds a fake dataset where you decide which probes truly differ. That lets you check the pipeline against a known answer. **All data here is simulated. Nothing comes from real patients or any lab.**

Everything uses base R, so there is nothing to install.

![Volcano plot of simulated tumor vs. normal data](output/volcano.png)

## The idea

A methylation beta value is a probe's methylation level, from 0 (never methylated) to 1 (always methylated). Beta values are easy to interpret but awkward for statistics because they are squeezed between 0 and 1. So the pipeline tests on **M-values** (the log2 odds of beta), which behave much better in t-tests, and reports on **delta beta** (the change in average beta) so the results stay interpretable.

## Pipeline steps

Each step is one small function in `R/pipeline.R`.

| Step | Function | What it does |
|---|---|---|
| 1 | `beta_to_m()`, `m_to_beta()` | Convert between beta and M-values, clamping exact 0s and 1s so they don't become infinity |
| 2 | `filter_probes()` | Quality control: drops probes with too many missing values or no variation, and returns the list of what was removed so nothing disappears silently |
| 3 | `test_differential()` | Welch t-test per probe on M-values, delta beta, and Benjamini-Hochberg FDR correction across all probes |
| 4 | `call_dmps()` | A probe is a hit only if it passes **both** an FDR cutoff (default 0.05) and a minimum \|delta beta\| (default 0.1) |
| 5 | `enrichment_test()` | Fisher's exact test on whether hits are over-represented in a region category, such as CpG islands |
| 6 | `volcano_plot()` | Plots change size against significance, with hits colored |
| | `run_pipeline()` | Runs steps 2 to 5 together |

Requiring both significance and effect size is standard practice, because a tiny difference can be statistically significant in a large study without meaning much.

## Project structure

```
.
├── run_analysis.R             # Runs the whole analysis end to end
├── R/
│   ├── pipeline.R             # The analysis functions
│   └── simulate.R             # Builds fake data with a known set of true DMPs
├── legacy/
│   └── legacy_analysis.R      # A deliberately messy "before" version, for comparison
├── tests/
│   └── test_pipeline.R        # Automated checks
└── output/                    # Created when you run the analysis
    ├── dmp_results.csv
    └── volcano.png
```

## Getting started

You only need [R](https://www.r-project.org/). Run everything from the project folder.

### Run the analysis

```bash
Rscript run_analysis.R
```

This simulates 2,000 probes across 20 samples (10 normal, 10 tumor) with 100 planted DMPs, runs the pipeline, and prints how well it recovered the planted truth: probes tested, probes removed by QC, correct hits, sensitivity, precision, and the CpG-island enrichment result. It saves two files:

- `output/dmp_results.csv`: one row per probe, sorted by p-value, with columns `probe`, `delta_beta`, `t_stat`, `p_value`, `fdr`, and `is_dmp`
- `output/volcano.png`: the volcano plot above

### Run the tests

```bash
Rscript tests/test_pipeline.R
```

Prints a `PASS` or `FAIL` line for every check, then a summary, and exits with a non-zero status if anything fails.

### Use it on your own data

```r
source("R/pipeline.R")

# beta: probes x samples matrix of values between 0 and 1, with probe names as row names
# group: one label per sample (column), exactly two distinct values
fit <- run_pipeline(beta, group, control = "normal")

fit$results          # per-probe statistics and hit calls
fit$removed_probes   # probes dropped by QC
fit$enrichment       # NULL unless you pass an annotation table and there are hits
```

To test for region enrichment, pass an `annotation` data frame with a `probe` column and an `island_relation` column.

## How the tests work

The tests are built around checking the code against situations where the right answer is known:

- **Unit checks** on tiny, hand-checkable inputs for conversions, QC, per-probe statistics, and hit calling (including that a tiny but significant change is *not* called a hit)
- **Input validation:** wrong-length group vectors, three groups, and a control that isn't a group all raise errors
- **Planted-truth run:** the full pipeline must reach at least 0.85 sensitivity and 0.90 precision, and the 20 strongest probes must all be truly different
- **Null simulation:** with no real differences across 20 runs, about 5% of raw p-values fall below 0.05 and the pipeline calls almost nothing
- **Power simulation:** sensitivity rises with effect size, and a very small effect is mostly missed
- **Enrichment:** a planted CpG-island signal is detected, and with no planted bias significant results stay rare
- **Refactor check:** the new pipeline matches the legacy script to within 1e-12 on delta beta, p-values, adjusted p-values, and hit calls

## About the legacy script

`legacy/legacy_analysis.R` stands in for the kind of script that accumulates over time: magic numbers, copy-pasted blocks, one-letter names, and no checks. It was written for this demo and is not code from any real lab. The refactored pipeline produces the same numbers, and the test suite proves it.

## Notes

- Group sizes in the simulator, effect size, noise, and the strength of the CpG-island bias are all adjustable through `simulate_methylation()` arguments (`n_per_group`, `n_dmp`, `m_shift`, `noise_sd`, `island_bias`, `seed`).
- This is a teaching and demonstration project, not a validated tool for clinical or research use.

## Author

Leher Harlalka
