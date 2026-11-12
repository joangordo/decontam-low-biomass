# -----------------------------------------------------------------------------
# Simulation study: when do decontam's methods work?
#
# Full factorial design, `n_reps` simulated studies per scenario. Results are
# cached in analysis/results/ so the report renders in seconds.
# The run is resumable: finished scenarios are skipped.
#
# Run from the project root:   Rscript analysis/run_simulation.R
# (about 15 minutes on one core)
# -----------------------------------------------------------------------------

for (f in list.files("R", full.names = TRUE)) source(f)

n_reps <- 20
design <- expand.grid(
  n_neg = c(1, 2, 4, 8),
  biomass_ratio = c(100, 10, 1),   # sample DNA / contaminant DNA (median)
  conc_noise = c(0.1, 0.5),        # precise vs imprecise DNA quantification
  leak = c(0, 0.005),              # well-to-well leakage into the blanks
  KEEP.OUT.ATTRS = FALSE
)
design$scenario <- seq_len(nrow(design))

# Each scenario is saved as soon as it finishes, so an interrupted run can be
# resumed by simply running the script again.
cache_dir <- "analysis/results/scenarios"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
for (sc in seq_len(nrow(design))) {
  out_file <- file.path(cache_dir, sprintf("scenario_%02d.rds", sc))
  if (file.exists(out_file)) next
  p <- design[sc, ]
  res_sc <- lapply(seq_len(n_reps), function(rep) {
    sim <- simulate_study(n_samples = 40, n_neg = p$n_neg,
                          biomass_ratio = p$biomass_ratio,
                          conc_noise = p$conc_noise, leak = p$leak,
                          seed = 1e4 * sc + rep)
    ev <- evaluate_methods(sim, thresholds = c(0.1, 0.25, 0.5))
    cbind(p[rep(1, nrow(ev)), ], rep = rep, ev, row.names = NULL)
  })
  saveRDS(do.call(rbind, res_sc), out_file)
  message(sprintf("scenario %2d/%d done (%.1f min)", sc, nrow(design),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

files <- file.path(cache_dir, sprintf("scenario_%02d.rds", design$scenario))
if (all(file.exists(files))) {
  results <- do.call(rbind, lapply(files, readRDS))
  saveRDS(list(results = results, design = design, n_reps = n_reps,
               r_version = R.version.string,
               decontam_version = as.character(utils::packageVersion("decontam"))),
          "analysis/results/simulation_results.rds")
  message("Saved ", nrow(results), " rows from ", length(files), " scenarios.")
}
