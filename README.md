# Contaminant removal in low-biomass 16S data

**How well do common contaminant-removal methods work when we know the right
answer?** This project evaluates the R package
[decontam](https://benjjneb.github.io/decontam/) and several simple filters in
two settings where every read can be labelled as genuine or contaminant: a
public mock-community dilution series and a simulation study. It ends with a
practical guide and a reusable cleaning workflow.

📄 **[Read the reports](https://joangordo.github.io/decontam-low-biomass/)**

| Part | Question |
|---|---|
| [1. Benchmark](https://joangordo.github.io/decontam-low-biomass/analysis/01_mock_dilution_benchmark.html) | How do decontam and simple filters perform on a real mock-community dilution series with known composition? |
| [2. Simulation study](https://joangordo.github.io/decontam-low-biomass/analysis/02_simulation_study.html) | How do the number of negative controls, DNA quantification and cross-contamination change the answer? (48 scenarios × 20 studies) |
| [3. Decision guide](https://joangordo.github.io/decontam-low-biomass/analysis/03_decision_guide.html) | Which method and threshold to use, an end-to-end workflow per sample type, and auto-generated methods text |

## Key findings

**On real data (mock-community dilution series)**

- Removing every taxon seen in the negative control deleted 3 of the 9
  genuine taxa and up to 34% of the true reads in a sample.
- decontam's frequency method removed 14% of contaminant reads in the most
  diluted samples at the default threshold (0.1) and 71% at 0.5.
- Among the ASVs it can score, the frequency method separates contaminants
  from genuine taxa almost perfectly (AUC 0.99), but 675 of 946 ASVs appear
  in a single sample and cannot be scored.

**In simulation (realistic, slightly leaky negative controls)**

- **DNA quantification matters more than the method.** In very-low-biomass
  studies the frequency method removed 76% of the contamination with precise
  quantification and 8% with imprecise quantification.
- **High thresholds remove real biology.** With realistic communities, a
  threshold of 0.5 removed about 25 rare genuine taxa per study, something
  the mock community could not show.
- **Prevalence-based methods are vulnerable to cross-contamination.** With 8
  leaky controls, the prevalence method removed 32% of genuine reads; more
  controls made it worse, not better.

**Recommendation:** decontam's frequency method with a threshold of 0.1–0.25,
applied per sample type and reported explicitly. It was the best safe choice
in 75% of the simulated designs; a threshold of 0.5 never met the "don't
remove real biology" rule. When biomass is very low and quantification
imprecise, no method removed more than about a fifth of the contamination:
that has to be fixed in the lab, not in software.

## Background

In low-biomass samples (biopsies, swabs, blood, ascitic fluid) the DNA that
comes from reagents and the lab environment can rival or exceed the DNA from
the sample itself. Typical cleaning workflows combine read-depth cut-offs,
abundance pre-filters and decontam at a fixed threshold, but they are rarely
evaluated, because on real samples there is no ground truth. This repository
builds such a workflow from scratch and evaluates it on open data where the
truth is known.

The question comes from my own experience analysing low-biomass 16S data
during a research internship. No data, code or documents from that work are
used here: everything is built on public data and the published literature.

### Design choices

| Common shortcut | This repository |
|---|---|
| Assume metadata and count-table rows are in the same order | IDs are matched and validated (`load_dataset()`) |
| Abundance pre-filter computed over samples *and* controls | Computed over true samples only |
| "Control-dominant" filter described but not implemented | Implemented and tested (`flag_control_dominant()`) |
| Low-read rules applied by hand | `flag_low_read_samples()` with both rules |
| One fixed threshold, often unreported | Threshold treated as a parameter and swept |
| Extra fits drawn over decontam's models | Only decontam's two models are shown |
| No evaluation | Scored against a mock-community ground truth and in simulation |
| Prevalence or frequency, either way | Frequency first, threshold 0.1–0.25 (Part 3) |
| Manual reporting | `clean_low_biomass()` returns every decision and drafts the methods text |

## Repository structure

```
├── R/
│   ├── io.R                 # loading and validation
│   ├── preprocess.R         # simple filters and read-depth rules
│   ├── decontam_methods.R   # decontam wrapper
│   ├── evaluate.R           # metrics against the ground truth
│   ├── simulate.R           # simulation model and per-study scoring
│   └── workflow.R           # clean_low_biomass(): end-to-end, per sample type
├── analysis/
│   ├── 01_mock_dilution_benchmark.qmd
│   ├── 02_simulation_study.qmd
│   ├── 03_decision_guide.qmd
│   ├── run_simulation.R     # regenerates results/simulation_results.rds
│   └── results/
├── data/                    # processed CSVs + data README (source and licence)
├── data-raw/                # builds data/ from the original files
├── docs/                    # rendered site (GitHub Pages)
├── tests/testthat/          # unit tests (hand-checkable examples + real data)
└── .github/workflows/       # runs the tests on every push
```

## Reproduce

Requires R (≥ 4.3) and [Quarto](https://quarto.org/).

```r
source("install.R")                      # CRAN + Bioconductor packages
testthat::test_dir("tests/testthat")     # run the tests
```

```bash
Rscript analysis/run_simulation.R        # optional: ~6 min, results are cached
quarto render                            # writes the site to docs/
```

The processed data and simulation results are included, so the reports
render without re-running anything. The simulation is fully seeded and
resumable.

## Using the workflow on your own data

```r
for (f in list.files("R", full.names = TRUE)) source(f)

res <- clean_low_biomass(
  counts, samples,                       # samples x taxa matrix + metadata
  group = "sample_type",                 # each type cleaned against its own controls
  method = "frequency", threshold = 0.1,
  low_read_rule = c(blood = "neg_median")
)
res$summary                              # what was done, per sample type
res$counts                               # cleaned counts (true samples only)
```

## Data and citation

Data from Karstens L. et al. (2019), *mSystems* 4(4): e00290-19,
[doi:10.1128/mSystems.00290-19](https://doi.org/10.1128/mSystems.00290-19)
(open-access article, CC BY 4.0). See [`data/README.md`](data/README.md) for the
source files and their terms.

decontam: Davis N.M. et al. (2018), *Microbiome* 6: 226,
[doi:10.1186/s40168-018-0605-2](https://doi.org/10.1186/s40168-018-0605-2).

Code released under the MIT licence.
