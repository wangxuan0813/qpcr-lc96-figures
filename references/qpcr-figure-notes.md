# qPCR Figure Notes

## Recommended qPCR Figure Types

### WT/KO stage-wise bars

Use when comparing genotype at each differentiation stage.

- Panel = stage/day.
- x = genes.
- bars = WT and KO.
- points = biological replicates.
- y = relative expression calibrated to same-day WT mean for each gene.
- error = SEM unless user requests SD.

### Gene-by-day heatmap

Use when summarizing direction across many days.

- rows = genes.
- columns = stages/days.
- value = log2 fold-change KO/WT or treatment/control.
- use diverging palette centered at 0.

### Dose-response qPCR

Use when comparing TPEN, ZnSO4, GSK126, or similar concentration series.

- x = dose.
- panel or color = gene.
- y = log2 fold-change or fold-change vs 0UM/DMSO.
- do not connect biologically unrelated conditions unless it is truly a dose series.

## Manuscript Defaults

- Use SVG as primary editable output.
- Also export PDF and 600 dpi PNG/TIFF.
- Show biological replicate points whenever feasible.
- Put high-Cq/QC-sensitive targets such as FEV in a supplement or separate panel if they distort the y-axis.
- Make p-value stars optional and label them exploratory unless multiple-testing correction and planned comparisons are specified.

## Caption Minimum

Include:

```text
Reference gene:
Calibration group:
Technical replicate handling:
Biological replicate n:
Center and error:
Statistical test:
QC exclusions or flags:
```

If no rows are excluded, say flagged rows were retained.
