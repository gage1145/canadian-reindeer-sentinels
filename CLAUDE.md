# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RT-QuIC prion surveillance project detecting chronic wasting disease (CWD) in Canadian reindeer using PrioScout sentinel devices. Samples are swabs from salt blocks placed in reindeer pens; they are tested using RT-QuIC assays and classified as Detected, Suspect, or Not Detected.

## Environment Setup

This project uses `renv` for reproducible package management with R 4.4.3.

```r
# Restore packages from lockfile
renv::restore()
```

The `quicR` package is a specialized package for RT-QuIC data processing — it handles raw MARS Excel output parsing and metric calculation.

## Running the Analysis Pipeline

Run scripts in this order:

```r
# 1. Curate raw MARS Excel files → data/data.csv
# 2. Calculate metrics → data/calculations.csv
source("analysis.R")

# 3. Generate plate view PNGs and summary figures → figures/
source("figures.R")

# 4. Salt block interference analysis (run independently)
source("salt_analysis.R")
```

Render reports:

```bash
quarto render report.qmd          # Round 1 report (HTML + PDF)
quarto render report_round2.qmd   # Round 2 report (PDF)
quarto render analysis.qmd        # Exploratory analysis (HTML)
```

## Data Flow

```
raw/*.xlsx (MARS instrument output)
    → analysis.R (quicR::get_quic + calculate_metrics)
    → data/data.csv (raw fluorescence curves)
    → data/calculations.csv (MPR, MS, AUC, TtT, RAF metrics + metadata)
    → data/results.csv (Detected/Suspect/Not Detected classifications)

raw/salt_block/*.xlsx
    → salt_analysis.R
    → data/salt_*.csv
```

Metadata is stored in `data/metadata.xlsx` and joined to calculations during `analysis.R`.

## Classification Criteria

- **Detected:** p ≤ 0.05 for all three metrics (MPR, MS, AUC) AND ≥50% threshold crossing
- **Suspect:** threshold positive but negative for one metric
- **Not Detected:** otherwise

Threshold for normalized fluorescence: 3× background (updated from 2× in Round 1).

## Key Metrics

| Metric | Description |
|--------|-------------|
| MPR | Maxpoint Ratio |
| MS | Max Slope |
| AUC | Area Under Curve |
| TtT | Time-to-Threshold |
| RAF | Rate of Amyloid Formation |

## Experimental Context

- 4 pens of CWD-inoculated reindeer
- PrioScout sentinel devices with salt blocks
- Weekly swabbing for 56 months
- RT-QuIC: 72 hours, 42°C, 4 replicates per sample
- Blinded sample identifiers are `random_id` (e.g., 117–228 for Round 2)

## Report Styling

Reports use Calibri font and custom kableExtra table styling with color-coded results (red = Detected, purple = Suspect). PDF outputs use 0.75-inch margins with `tbl-colwidths` set per table.
