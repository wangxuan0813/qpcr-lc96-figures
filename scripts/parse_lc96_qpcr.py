#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import zipfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import mean, stdev
from xml.etree import ElementTree as ET


REFERENCE_CANDIDATES = ("GAPDH", "ACTB", "BACT", "B2M", "HPRT1", "RPLP0", "TBP", "18S")
CONTROL_PRIORITY = ("WT", "0UM", "0", "DMSO", "-DOX", "CT", "SHNC", "SHCON", "CTRL", "CONTROL", "UNTREATED")


@dataclass(frozen=True)
class CqRecord:
    source_file: str
    rel_dir: str
    batch: str
    run_label: str
    well: str
    sample: str
    sample_raw_group: str
    group: str
    replicate: str
    target: str
    cq: float | None
    call: str | None


def strip_ns(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def children_by_name(node: ET.Element, name: str) -> list[ET.Element]:
    return [child for child in list(node) if strip_ns(child.tag) == name]


def child_text(node: ET.Element, name: str) -> str | None:
    for child in list(node):
        if strip_ns(child.tag) == name:
            return child.text
    return None


def attr_id(node: ET.Element | None, name: str) -> str | None:
    if node is None:
        return None
    for child in list(node):
        if strip_ns(child.tag) == name:
            return child.attrib.get("id")
    return None


def clean_target(target_id: str) -> str:
    return target_id.split("@", 1)[-1].replace("None", "").strip()


def well_name(row: int, col: int) -> str:
    return f"{chr(ord('A') + row - 1)}{col}"


def react_id_to_well(react_id: str) -> str:
    try:
        idx = int(react_id)
    except ValueError:
        return react_id
    row = ((idx - 1) // 12) + 1
    col = ((idx - 1) % 12) + 1
    return well_name(row, col)


def normalize_token(token: str) -> str:
    return re.sub(r"[^A-Za-z0-9.+-]+", "", token).upper()


def infer_sample_parts(sample: str, batch: str, swap_label_patterns: list[str]) -> tuple[str, str, str]:
    text = sample.strip()
    replicate = ""
    rep_match = re.search(r"(?:^|[-_\s])(?:R|REP)?(\d+)$", text, re.I)
    if rep_match:
        replicate = rep_match.group(1)

    base = re.sub(r"[-_\s]*(?:R|REP)?\d+$", "", text, flags=re.I).strip(" -_")
    batch_norm = normalize_token(batch)

    group = ""
    if re.search(r"^\+DOX\b|\+DOX", text, re.I):
        group = "+DOX"
    elif re.search(r"^-DOX\b|-DOX", text, re.I):
        group = "-DOX"
    elif re.search(r"OE.*GSK126|GSK126.*OE", text, re.I):
        group = "OE+GSK126"
    elif re.search(r"DMSO.*GSK126|GSK126.*DMSO", text, re.I):
        group = "DMSO+GSK126"
    elif re.search(r"\bDMSO\b", text, re.I):
        group = "DMSO"
    elif re.search(r"\bCT\b|^CT-", text, re.I):
        group = "CT"
    elif re.search(r"^DOX(?:-|$)", text, re.I):
        group = "DOX"
    elif re.search(r"^T-", text, re.I):
        group = "T"
    elif re.search(r"\bWT\b|WT-", text, re.I):
        group = "WT"
    elif re.search(r"\bKO\b|KO-", text, re.I):
        group = "KO"
    elif re.search(r"SH\s*PAX6|SHPAX6", text, re.I):
        group = "shPAX6"
    elif re.search(r"SH(?:NC|CON|CTRL)", text, re.I):
        group = "shCtrl"
    elif re.search(r"PAX6\s*OE|PAX6OE|\bOE\b", text, re.I):
        group = "OE"
    elif re.search(r"GSK126", text, re.I):
        group = "GSK126"
    else:
        conc_match = re.search(r"(?<![0-9.])(50|10|5|2|1|0\.5|0)\s*U?M\b", text, re.I)
        if conc_match:
            group = f"{conc_match.group(1)}UM"
        elif ("TPEN" in batch_norm or "ZNSO4" in batch_norm) and re.search(r"(?:^|[-_])(\d+(?:\.\d+)?)(?:[-_]|$)", base, re.I):
            conc_match = re.search(r"(?:^|[-_])(\d+(?:\.\d+)?)(?:[-_]|$)", base, re.I)
            if conc_match:
                group = f"{conc_match.group(1)}UM"

    if not group:
        local_dir = batch.split("\\")[-1]
        dir_conc = re.search(r"(?<![0-9.])(50|10|5|2|1|0\.5|0)\s*U?M\b", local_dir, re.I)
        if dir_conc and "ZNSO4" in batch_norm:
            group = f"{dir_conc.group(1)}UM_ZNSO4"
        elif "TPEN" in batch_norm and dir_conc:
            group = f"{dir_conc.group(1)}UM_TPEN"

    if not group:
        base_norm = normalize_token(base)
        stage_match = re.match(r"^(S\d+D?\d*)", base_norm, re.I)
        if stage_match:
            group = stage_match.group(1).upper()
        elif base_norm.startswith("ES"):
            group = "ES"
        else:
            cleaned = re.sub(r"(S\d+D\d+|D\d+|S\d+|HUSE8|AD-?B|MEL1)", "", text, flags=re.I)
            cleaned = re.sub(r"[-_\s]*(?:R|REP)?\d+$", "", cleaned, flags=re.I).strip(" -_")
            group = cleaned or "UNASSIGNED"

    raw_group = group
    if any(pattern.lower() in batch.lower() for pattern in swap_label_patterns):
        if re.search(r"\bWT-\d+\b|WT\d*$", text, re.I):
            group = "KO"
        elif re.search(r"\bKO-\d+\b|KO\d*$", text, re.I):
            group = "WT"
        raw_group = raw_group + " (source WT/KO swapped in analysis)"
    return raw_group, group, replicate


def safe_mean(values) -> float | None:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    return mean(vals) if vals else None


def safe_sd(values) -> float | None:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    return stdev(vals) if len(vals) > 1 else None


def parse_lc96p(path: Path, root: Path, swap_label_patterns: list[str]) -> list[CqRecord]:
    with zipfile.ZipFile(path) as zf:
        if "rdml_data.xml" not in zf.namelist():
            return []
        rdml_xml = zf.read("rdml_data.xml")
        calculated_xml = zf.read("calculated_data.xml") if "calculated_data.xml" in zf.namelist() else None

    rdml = ET.fromstring(rdml_xml)
    calculated = ET.fromstring(calculated_xml) if calculated_xml else None

    samples: dict[str, str] = {}
    for sample_node in children_by_name(rdml, "sample"):
        sid = sample_node.attrib.get("id")
        desc = child_text(sample_node, "description")
        if sid and desc:
            samples[sid] = desc.strip()

    calls: dict[str, str] = {}
    if calculated is not None:
        for react in calculated.iter():
            if strip_ns(react.tag) != "react":
                continue
            react_id = react.attrib.get("id")
            call = None
            for fact_graph in children_by_name(react, "factGraph"):
                call = child_text(fact_graph, "call")
                break
            if react_id:
                calls[react_id] = call or ""

    rel_dir = str(path.parent.relative_to(root))
    batch = rel_dir.split("\\")[0] if "\\" in rel_dir else rel_dir
    records: list[CqRecord] = []
    for react in rdml.iter():
        if strip_ns(react.tag) != "react":
            continue
        react_id = react.attrib.get("id", "")
        sample_id = attr_id(react, "sample")
        sample = samples.get(sample_id or "", "")
        data = children_by_name(react, "data")
        target_id = attr_id(data[0], "tar") if data else None
        cq_text = child_text(data[0], "cq") if data else None
        cq = None
        if cq_text:
            try:
                cq = float(cq_text)
            except ValueError:
                cq = None
        row = child_text(react, "pcrFormat_row")
        col = child_text(react, "pcrFormat_col")
        well = react_id_to_well(react_id)
        if row and col:
            well = well_name(int(row), int(col))
        raw_group, group, replicate = infer_sample_parts(sample, rel_dir, swap_label_patterns)
        run_match = re.search(r"S\d+D\d+|S\d+D?\d*|D\d+", sample, re.I)
        run_label = run_match.group(0).upper() if run_match else path.stem
        records.append(
            CqRecord(path.name, rel_dir, batch, run_label, well, sample, raw_group, group, replicate, clean_target(target_id or ""), cq, calls.get(react_id))
        )
    return records


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summarize(records: list[CqRecord], paths: list[Path], root: Path, out_dir: Path) -> None:
    raw_rows = [asdict(r) for r in records]
    write_csv(out_dir / "all_raw_cq_long.csv", raw_rows)

    file_rows = []
    for path in paths:
        rel_dir = str(path.parent.relative_to(root))
        file_records = [r for r in records if r.source_file == path.name and r.rel_dir == rel_dir]
        file_rows.append(
            {
                "rel_dir": rel_dir,
                "source_file": path.name,
                "reactions": len(file_records),
                "samples": len({r.sample for r in file_records}),
                "targets": ", ".join(sorted({r.target for r in file_records if r.target})),
                "groups_inferred": ", ".join(sorted({r.group for r in file_records if r.group})),
                "negative_calls": sum(1 for r in file_records if r.call == "Negative"),
                "missing_cq": sum(1 for r in file_records if r.cq is None),
            }
        )
    write_csv(out_dir / "file_inventory.csv", file_rows)

    sample_keys = sorted({(r.rel_dir, r.source_file, r.batch, r.run_label, r.sample, r.group, r.replicate, r.target) for r in records})
    sample_rows = []
    for rel_dir, source_file, batch, run_label, sample, group, replicate, target in sample_keys:
        rows = [
            r for r in records
            if (r.rel_dir, r.source_file, r.batch, r.run_label, r.sample, r.group, r.replicate, r.target)
            == (rel_dir, source_file, batch, run_label, sample, group, replicate, target)
        ]
        cqs = [r.cq for r in rows if r.cq is not None]
        calls = ";".join(sorted({r.call or "" for r in rows if r.call}))
        sample_rows.append(
            {
                "rel_dir": rel_dir,
                "source_file": source_file,
                "batch": batch,
                "run_label": run_label,
                "sample": sample,
                "group": group,
                "replicate": replicate,
                "target": target,
                "technical_n": len(cqs),
                "mean_cq": safe_mean(cqs),
                "sd_cq": safe_sd(cqs),
                "calls": calls,
            }
        )
    write_csv(out_dir / "sample_target_mean_cq_all.csv", sample_rows)

    qc_rows = []
    for row in sample_rows:
        flags = []
        if "Negative" in str(row["calls"]):
            flags.append("Negative call")
        if row["mean_cq"] is not None and row["mean_cq"] >= 35:
            flags.append("Mean Cq >= 35")
        if row["sd_cq"] is not None and row["sd_cq"] > 0.5:
            flags.append("Technical SD(Cq) > 0.5")
        if not row["group"] or row["group"] == "UNASSIGNED":
            flags.append("Group not inferred")
        if flags:
            qc_rows.append({**row, "qc_flag": "; ".join(flags)})
    write_csv(out_dir / "qc_flags_all.csv", qc_rows)

    delta_rows = []
    by_scope = defaultdict(list)
    for row in sample_rows:
        by_scope[(row["rel_dir"], row["source_file"])].append(row)
    for (rel_dir, source_file), rows in by_scope.items():
        targets = {str(r["target"]).upper(): r["target"] for r in rows}
        ref = next((targets[g] for g in REFERENCE_CANDIDATES if g in targets), None)
        if not ref:
            continue
        ref_by_sample = {(r["sample"], r["run_label"]): r["mean_cq"] for r in rows if r["target"] == ref}
        for row in rows:
            if row["target"] == ref:
                continue
            ref_cq = ref_by_sample.get((row["sample"], row["run_label"]))
            delta_cq = row["mean_cq"] - ref_cq if row["mean_cq"] is not None and ref_cq is not None else None
            delta_rows.append({**row, "ref_gene": ref, "ref_cq": ref_cq, "delta_cq": delta_cq})
    write_csv(out_dir / "delta_ct_all.csv", delta_rows)

    summary_rows = []
    summary_keys = sorted({(r["rel_dir"], r["source_file"], r["run_label"], r["target"], r["group"]) for r in delta_rows})
    for rel_dir, source_file, run_label, target, group in summary_keys:
        rows = [
            r for r in delta_rows
            if (r["rel_dir"], r["source_file"], r["run_label"], r["target"], r["group"]) == (rel_dir, source_file, run_label, target, group)
        ]
        dcts = [r["delta_cq"] for r in rows if r["delta_cq"] is not None]
        summary_rows.append(
            {
                "rel_dir": rel_dir,
                "source_file": source_file,
                "run_label": run_label,
                "target": target,
                "group": group,
                "n": len(dcts),
                "mean_delta_cq": safe_mean(dcts),
                "sd_delta_cq": safe_sd(dcts),
            }
        )
    write_csv(out_dir / "delta_ct_summary_all.csv", summary_rows)

    ddct_rows = []
    for rel_dir, source_file, run_label, target in sorted({(r["rel_dir"], r["source_file"], r["run_label"], r["target"]) for r in summary_rows}):
        rows = [
            r for r in summary_rows
            if (r["rel_dir"], r["source_file"], r["run_label"], r["target"]) == (rel_dir, source_file, run_label, target)
        ]
        controls = []
        for control_name in CONTROL_PRIORITY:
            controls = [r for r in rows if str(r["group"]).upper() == control_name]
            if controls:
                break
        if not controls:
            for control_name in CONTROL_PRIORITY:
                controls = [r for r in rows if control_name in str(r["group"]).upper()]
                if controls:
                    break
        if not controls:
            controls = rows[:1]
        control = controls[0]
        if control["mean_delta_cq"] is None:
            continue
        for row in rows:
            if row["mean_delta_cq"] is None:
                continue
            ddct = row["mean_delta_cq"] - control["mean_delta_cq"]
            ddct_rows.append(
                {
                    **row,
                    "control_group": control["group"],
                    "delta_delta_cq": ddct,
                    "relative_expression_vs_control": 2 ** (-ddct),
                    "control_selection_note": "auto",
                }
            )
    write_csv(out_dir / "relative_expression_auto_all.csv", ddct_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse Roche LightCycler/LC96 qPCR .lc96p files into analysis-ready CSVs.")
    parser.add_argument("--input", required=True, type=Path, help="Input qPCR directory")
    parser.add_argument("--out", required=True, type=Path, help="Output analysis directory")
    parser.add_argument("--swap-label-dir", action="append", default=[], help="Substring of relative directory where WT/KO source labels should be swapped")
    args = parser.parse_args()

    root = args.input.resolve()
    out_dir = args.out.resolve()
    paths = [p for p in sorted(root.rglob("*.lc96p")) if "node_modules" not in p.parts and "outputs" not in p.parts]
    records: list[CqRecord] = []
    for path in paths:
        records.extend(parse_lc96p(path, root, args.swap_label_dir))
    summarize(records, paths, root, out_dir)
    print(f"lc96p_files={len(paths)}")
    print(f"raw_records={len(records)}")
    print(f"outputs={out_dir}")


if __name__ == "__main__":
    main()
