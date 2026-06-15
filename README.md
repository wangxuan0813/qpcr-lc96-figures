# qPCR LC96 Figures

Codex skill for end-to-end qPCR analysis from Roche LightCycler/LC96 project files to publication-ready figures.

The skill parses `.lc96p`/`.lc96u` source files, extracts RDML Cq/Cp values, summarizes technical replicates, computes Delta Ct and relative expression, flags common QC issues, and exports WT/KO stage-wise qPCR bar plots with source data.

## Features

- Parse Roche LightCycler project containers without modifying raw files.
- Extract sample, target, well, Cq/Cp, call, and plate metadata.
- Summarize technical replicates before biological replicate statistics.
- Generate QC tables for missing Cq, Negative calls, high Cq, high technical SD, and unclear groups.
- Compute Delta Ct, DeltaDeltaCt, fold-change, and log2 fold-change tables.
- Draw WT/KO stage-wise grouped bars with mean + SEM, biological replicate points, source CSVs, Welch-test summaries, and SVG/PDF/PNG/TIFF exports.
- Support known label corrections such as swapped WT/KO directories through `--swap-label-dir`.

## Install As A Codex Skill

Clone the repository into your Codex skills directory:

```powershell
git clone https://github.com/wangxuan0813/qpcr-lc96-figures.git "$env:USERPROFILE\.codex\skills\qpcr-lc96-figures"
```

Restart Codex so the skill is discovered. Then ask Codex to use `$qpcr-lc96-figures`.

## Python Dependencies

Create or use a Python environment with:

```powershell
pip install -r requirements.txt
```

## Quick Start

Parse a qPCR directory:

```powershell
python scripts/parse_lc96_qpcr.py --input "path\to\qpcr" --out "path\to\outputs"
```

Apply a known WT/KO label swap while parsing:

```powershell
python scripts/parse_lc96_qpcr.py --input "path\to\qpcr" --out "path\to\outputs" --swap-label-dir "20251014 AD-B WT KO"
```

Draw stage-wise WT/KO bars:

```powershell
python scripts/plot_stage_wtko_bars.py --analysis-dir "path\to\outputs" --rel-dir "relative\dir\containing\lc96p" --out "path\to\figures"
```

Choose genes and include FEV:

```powershell
python scripts/plot_stage_wtko_bars.py --analysis-dir "path\to\outputs" --rel-dir "..." --genes ZNT8 PAX6 INS NKX6-1 PAX4 TPH1 --include-fev --out "..."
```

## Main Outputs

The parser writes:

- `all_raw_cq_long.csv`
- `file_inventory.csv`
- `sample_target_mean_cq_all.csv`
- `qc_flags_all.csv`
- `delta_ct_all.csv`
- `delta_ct_summary_all.csv`
- `relative_expression_auto_all.csv`

The plotter writes:

- editable `.svg`
- `.pdf`
- high-resolution `.png` and `.tiff`
- `source_stage_bar_sample_values.csv`
- `source_stage_bar_welch_tests.csv`
- `figure_qa_notes.txt`

## Data Safety

Raw `.lc96p` and `.lc96u` files are read-only inputs. The scripts write outputs only to the directory passed through `--out`.

This repository does not include experimental raw data.

## License

MIT License.
