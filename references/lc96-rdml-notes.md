# LC96 / RDML Notes

## File Structure

Roche LightCycler 96 `.lc96p` files are ZIP containers. Useful entries commonly include:

- `rdml_data.xml`: samples, targets, reactions, Cq, amplification points.
- `calculated_data.xml`: call status such as Positive/Negative and fitted curve details.
- `instrument_data.xml`, `app_data.xml`, `module_data.xml`: metadata.

`.lc96u` files may contain setup/user information and may not include calculated Cq data.

## Parsing Rules

- Read `.lc96p` as ZIP; parse `rdml_data.xml`.
- Extract samples from `<sample id=...><description>...`.
- Extract target IDs from `<react><data><tar id=...>`.
- Extract Cq from `<react><data><cq>`.
- Extract call status from `calculated_data.xml` when available.
- If row/column fields are absent, map reaction ID 1-96 to A1-H12.
- Strip target prefixes such as `SYBR Green I@`.
- Preserve both the raw sample name and inferred group.

## Reference Gene Selection

Preferred reference genes:

1. GAPDH
2. ACTB / BACT
3. B2M
4. HPRT1
5. RPLP0
6. TBP
7. 18S

Never assume a reference gene if none exists in the target list.

## Group Inference

Infer cautiously from sample names and directory names. Common patterns:

- `WT`, `KO`
- `+DOX`, `-DOX`, `DOX`, `CT`
- `DMSO`, `GSK126`, `OE`, `OE+GSK126`
- `0UM`, `0.5UM`, `1UM`, `2UM`, `5UM`, `10UM`, `50UM`
- stage labels such as `S5D1`, `S6D6`, `S7D10`

Always expose inferred groups in output tables for user review.

## QC Flags

Flag but do not delete:

- `Negative` call.
- mean Cq >= 35.
- technical replicate SD(Cq) > 0.5.
- missing Cq.
- group not inferred.

Rows with high-Cq EC markers such as FEV may be biologically interesting but should not be overinterpreted.
