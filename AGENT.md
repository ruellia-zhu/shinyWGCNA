# shinyWGCNA Maintenance Notes

## Project Overview

- This repository is a single-file Shiny app for WGCNA workflows on a Linux HPC server.
- The main application file is `shinyWGCNA.R`.
- `README.md` is intentionally minimal. Put agent-facing maintenance guidance in this file.

## Common Maintenance Areas

- Dependencies are loaded near the top of `shinyWGCNA.R`. Follow the existing pattern: install missing CRAN packages with `install.packages()` and Bioconductor packages with `BiocManager::install()`.
- Uploaded tables are read through the shared upload helpers near the top of `shinyWGCNA.R`. Update those helpers first when adding support for new table formats.
- Expression matrix uploads must keep one gene ID column followed by numeric sample columns.
- Trait matrix uploads must keep one sample ID column followed by numeric trait columns.
- ShinyWGCNA compatibility wrappers such as `call_getpower()`, `call_powertest()`, `call_getMt()`, and `call_getKME()` protect this app from upstream function signature differences. Keep new compatibility logic close to those wrappers.

## Version And UI Updates

- Update the file header `Update date` and version history comments when changing user-visible behavior.
- Update `customLogo` `badgeText` when bumping the app version.
- Keep Data import labels aligned with the internal ShinyWGCNA values:
  - `Format` label `count` maps to internal value `count`.
  - `Format` label `FPKM/TPM/CPM` maps to internal value `normalized count`.
  - `Normalized method` label `FPKM/TPM/CPM` maps to internal value `raw`.
  - `Normalized method` label `log10(FPKM/TPM/CPM)` maps to internal value `logarithm`.

## Update Checklist

- Check expression matrix upload with CSV, tab-delimited text, and XLSX files.
- Check trait matrix upload with CSV, tab-delimited text, and XLSX files.
- Confirm Data import method switching:
  - `count` selects `VST` and uses expression cutoff `10`.
  - `FPKM/TPM/CPM` selects raw/logarithm choices and uses expression cutoff `1`.
- Confirm the app version badge matches the current release.
- Confirm the Before Use page does not reintroduce old note text unless requested.

## Test Suggestions

- Run the app in an R/Shiny environment and load a small expression matrix with at least one gene column and two sample columns.
- Run Module-trait with a small trait table whose first column matches expression sample IDs.
- If `Rscript` is available, run an R syntax check before handing off changes.
