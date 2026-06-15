---
name: qpcr-lc96-figures
description: End-to-end qPCR workflow for Roche LightCycler/LC96 raw files and publication figures. Use when Codex needs to parse .lc96p/.lc96u qPCR projects, extract RDML Cq/Cp values, audit sample/target/plate metadata, compute technical-replicate means, Delta Ct, DeltaDeltaCt, fold-change/log2FC, QC flags, or create Python/matplotlib publication-ready qPCR figures such as WT/KO stage-wise bars, dose-response plots, or gene-by-day heatmaps from qPCR source files.
---

# qPCR LC96 Figures

## Core Contract

1. Preserve raw `.lc96p/.lc96u` files. Never rename, move, or overwrite instrument files.
2. Start from raw LC96 projects when available; use old Excel/GraphPad files only as reference unless the user explicitly asks to use them.
3. Audit before interpreting:
   - samples and inferred groups,
   - targets and reference gene,
   - missing Cq,
   - Negative calls,
   - technical replicate SD,
   - known label corrections such as swapped WT/KO.
4. Average technical replicates before biological replicate statistics.
5. Keep flagged data visible. Do not silently remove high-Cq, Negative, or high technical-SD rows.
6. For manuscript figures, export SVG first with editable text, plus PDF and high-DPI PNG/TIFF when possible.

## Workflow

### 1. Parse and normalize LC96 projects

Run the parser against a qPCR directory:

```powershell
python scripts/parse_lc96_qpcr.py --input "path\to\qpcr" --out "path\to\outputs"
```

Optional controls:

```powershell
python scripts/parse_lc96_qpcr.py --input "path\to\qpcr" --out "path\to\outputs" --swap-label-dir "20251014 AD-B WT KO"
```

Outputs:

- `all_raw_cq_long.csv`: one row per reaction/well.
- `file_inventory.csv`: one row per `.lc96p`.
- `sample_target_mean_cq_all.csv`: technical replicate summary.
- `qc_flags_all.csv`: Negative/high-Cq/high technical SD/group issues.
- `delta_ct_all.csv`: sample-level Delta Ct.
- `delta_ct_summary_all.csv`: group summaries.
- `relative_expression_auto_all.csv`: auto-selected control contrasts; always review `control_group`.

### 2. Decide what can be interpreted

Use `file_inventory.csv` and `qc_flags_all.csv`.

Rules of thumb:

- **Use** files with clear groups, reference gene present, no missing Cq, and low/moderate QC flags.
- **Use cautiously** files with high-Cq/Negative markers in only a few targets, especially FEV or late EC markers.
- **Rerun** files with all/most missing Cq, wrong acquisition, wrong annealing temperature, bad target labels, or ambiguous group labels.

### 3. Plot WT/KO stage-wise bars

For stage-wise WT/KO qPCR:

```powershell
python scripts/plot_stage_wtko_bars.py --analysis-dir "path\to\outputs" --rel-dir "relative\dir\containing\lc96p" --out "path\to\figures"
```

Common options:

```powershell
python scripts/plot_stage_wtko_bars.py --analysis-dir "path\to\outputs" --rel-dir "..." --genes ZNT8 PAX6 INS NKX6-1 PAX4 TPH1 --out "..."
python scripts/plot_stage_wtko_bars.py --analysis-dir "path\to\outputs" --rel-dir "..." --include-fev --out "..."
```

The script saves source data and exports `.svg`, `.pdf`, `.png`, and `.tiff`.

## Figure Decisions

Before plotting, state:

- Core conclusion.
- Control group.
- Reference gene.
- Biological vs technical replicate definition.
- Whether FEV/high-Cq markers belong in the main panel or a supplement.
- Whether p-value stars are exploratory or should be omitted.

For `nature-figure` style output, use Python/matplotlib and follow that skill's backend rules when it is available.

## References

- Read `references/lc96-rdml-notes.md` when LC96 XML parsing or group inference is uncertain.
- Read `references/qpcr-figure-notes.md` before making manuscript-style qPCR figures.
