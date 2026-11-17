# Data

Processed from the supplementary material of:

> Karstens L, Asquith M, Davin S, Fair D, Gregory WT, Wolfe AJ, Braun J,
> McWeeney S (2019). Controlling for contaminants in low-biomass 16S rRNA gene
> sequencing experiments. *mSystems* 4(4): e00290-19.
> https://doi.org/10.1128/mSystems.00290-19

Original files: <https://github.com/lakarstens/ControllingContaminants16S>.
Raw reads: SRA accession SRP155048. The article is open access under CC BY 4.0;
the source repository above has no explicit licence, so the data are used here
for research and teaching, with full citation. The files here are a format conversion (phyloseq `.RData` to CSV) made by
`data-raw/extract_karstens.R`, with no other changes.

| File | Content |
|---|---|
| `counts.csv` | ASV read counts, one row per sample (9 dilutions + 1 extraction blank) |
| `samples.csv` | `sample_id`, description, DNA concentration (NanoDrop, ng/µL), `is_neg`, dilution step, total reads |
| `taxonomy.csv` | Taxonomy per ASV and `expected` = TRUE for the 9 ASVs of the mock community |
