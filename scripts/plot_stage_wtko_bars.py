from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


DEFAULT_GENES = ["ZNT8", "PAX6", "INS", "NKX6-1", "PAX4", "TPH1"]
COLORS = {"WT": "#7C8792", "KO": "#B64342"}
EDGE = "#272727"
GRID = "#D8DEE8"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw stage-wise WT/KO grouped qPCR bars from parse_lc96_qpcr.py outputs."
    )
    parser.add_argument("--analysis-dir", required=True, type=Path, help="Directory containing delta_ct_all.csv.")
    parser.add_argument("--rel-dir", required=True, help="Relative source directory to plot, as stored in rel_dir.")
    parser.add_argument("--out", required=True, type=Path, help="Output figure directory.")
    parser.add_argument("--genes", nargs="*", default=DEFAULT_GENES, help="Genes/targets to plot in order.")
    parser.add_argument("--include-fev", action="store_true", help="Append FEV to --genes when absent.")
    parser.add_argument("--days", nargs="*", default=None, help="Run labels/days to plot in order; default auto-sorted.")
    parser.add_argument("--reference", default="GAPDH", help="Reference gene name for the caption.")
    parser.add_argument("--title", default=None, help="Figure title.")
    parser.add_argument("--stem", default="wtko_stage_wise_gene_bars", help="Output filename stem.")
    parser.add_argument("--no-stars", action="store_true", help="Do not draw exploratory p-value stars.")
    return parser.parse_args()


def apply_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.size"] = 7
    plt.rcParams["axes.linewidth"] = 0.65
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["xtick.major.width"] = 0.6
    plt.rcParams["ytick.major.width"] = 0.6


def natural_key(value: str) -> tuple:
    parts = re.split(r"(\d+)", str(value))
    return tuple(int(part) if part.isdigit() else part.lower() for part in parts)


def sem(values: pd.Series) -> float:
    vals = values.dropna().astype(float)
    if len(vals) <= 1:
        return 0.0
    return float(vals.std(ddof=1) / np.sqrt(len(vals)))


def p_label(p_value: float) -> str:
    if not np.isfinite(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def load_and_calibrate(analysis_dir: Path, rel_dir: str, genes: list[str], days: list[str] | None, out_dir: Path) -> pd.DataFrame:
    delta_path = analysis_dir / "delta_ct_all.csv"
    if not delta_path.exists():
        raise FileNotFoundError(f"Missing {delta_path}; run parse_lc96_qpcr.py first.")

    delta = pd.read_csv(delta_path, encoding="utf-8-sig")
    required = {"rel_dir", "run_label", "target", "group", "sample", "delta_cq"}
    missing = required.difference(delta.columns)
    if missing:
        raise ValueError(f"{delta_path} is missing required columns: {sorted(missing)}")

    df = delta.loc[delta["rel_dir"].astype(str) == rel_dir].copy()
    if df.empty:
        candidates = sorted(delta["rel_dir"].dropna().astype(str).unique())
        preview = "\n".join(candidates[:25])
        raise ValueError(f"No rows matched --rel-dir {rel_dir!r}. Available rel_dir values include:\n{preview}")

    df = df.loc[df["target"].astype(str).isin(genes)].copy()
    if df.empty:
        raise ValueError(f"No rows matched requested genes: {genes}")

    if days is None:
        days = sorted(df["run_label"].dropna().astype(str).unique(), key=natural_key)
    df = df.loc[df["run_label"].astype(str).isin(days)].copy()
    if df.empty:
        raise ValueError(f"No rows matched requested days: {days}")

    df["delta_cq"] = pd.to_numeric(df["delta_cq"], errors="coerce")
    df["run_label"] = pd.Categorical(df["run_label"].astype(str), days, ordered=True)
    df["target"] = pd.Categorical(df["target"].astype(str), genes, ordered=True)
    df["group"] = df["group"].astype(str)

    if not {"WT", "KO"}.issubset(set(df["group"])):
        raise ValueError("Plot requires both WT and KO groups after parsing/group inference.")

    wt_mean = (
        df.loc[df["group"] == "WT"]
        .groupby(["run_label", "target"], observed=True)["delta_cq"]
        .mean()
        .rename("wt_mean_delta_cq")
        .reset_index()
    )
    df = df.merge(wt_mean, on=["run_label", "target"], how="left")
    df["delta_delta_cq"] = df["delta_cq"] - df["wt_mean_delta_cq"]
    df["relative_expression_vs_wt"] = np.power(2.0, -df["delta_delta_cq"])
    df["log2_expression_vs_wt"] = -df["delta_delta_cq"]
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "source_stage_bar_sample_values.csv", index=False, encoding="utf-8-sig")
    return df


def make_stats(df: pd.DataFrame, genes: list[str], days: list[str], out_dir: Path) -> pd.DataFrame:
    rows = []
    for day in days:
        for gene in genes:
            wt = df[(df["run_label"] == day) & (df["target"] == gene) & (df["group"] == "WT")][
                "log2_expression_vs_wt"
            ].astype(float)
            ko = df[(df["run_label"] == day) & (df["target"] == gene) & (df["group"] == "KO")][
                "log2_expression_vs_wt"
            ].astype(float)
            p_value = np.nan
            if len(wt) >= 2 and len(ko) >= 2:
                p_value = float(stats.ttest_ind(wt, ko, equal_var=False, nan_policy="omit").pvalue)
            rows.append({"day": day, "target": gene, "welch_p_log2": p_value, "p_label": p_label(p_value)})
    stat_df = pd.DataFrame(rows)
    stat_df.to_csv(out_dir / "source_stage_bar_welch_tests.csv", index=False, encoding="utf-8-sig")
    return stat_df


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )


def plot_day(
    ax: plt.Axes,
    df: pd.DataFrame,
    stat_df: pd.DataFrame,
    day: str,
    genes: list[str],
    show_stars: bool,
) -> None:
    sub = df[df["run_label"] == day].copy()
    x = np.arange(len(genes))
    width = 0.34
    rng = np.random.default_rng(sum(ord(char) for char in day))
    max_y = 0.0

    for group in ["WT", "KO"]:
        g = sub[sub["group"] == group]
        summary = (
            g.groupby("target", observed=True)["relative_expression_vs_wt"]
            .agg(mean="mean", sem=sem)
            .reindex(genes)
        )
        xpos = x + (-width / 2 if group == "WT" else width / 2)
        ax.bar(
            xpos,
            summary["mean"].astype(float),
            yerr=summary["sem"].astype(float),
            width=width,
            color=COLORS[group],
            edgecolor=EDGE,
            linewidth=0.55,
            error_kw={"elinewidth": 0.7, "capsize": 2.5, "capthick": 0.7},
            label=group,
            zorder=2,
        )
        for idx, gene in enumerate(genes):
            vals = g[g["target"] == gene]["relative_expression_vs_wt"].dropna().astype(float).values
            jitter = rng.normal(0, 0.025, len(vals))
            ax.scatter(
                np.full(len(vals), xpos[idx]) + jitter,
                vals,
                s=9,
                color="white",
                edgecolor=EDGE,
                linewidth=0.45,
                zorder=3,
            )
            if len(vals):
                mean_y = summary.loc[gene, "mean"]
                sem_y = summary.loc[gene, "sem"]
                max_y = max(max_y, float(np.nanmax(vals)), float(mean_y + sem_y))

    stats_day = stat_df[stat_df["day"] == day].set_index("target")
    if show_stars:
        for idx, gene in enumerate(genes):
            label = stats_day.loc[gene, "p_label"] if gene in stats_day.index else ""
            if label:
                gene_vals = sub[sub["target"] == gene]["relative_expression_vs_wt"].dropna().astype(float)
                if gene_vals.empty:
                    continue
                y = float(gene_vals.max()) * 1.12 + 0.05
                ax.plot([x[idx] - width / 2, x[idx] + width / 2], [y, y], color=EDGE, lw=0.55)
                ax.text(x[idx], y * 1.03, label, ha="center", va="bottom", fontsize=6)
                max_y = max(max_y, y * 1.18)

    ax.axhline(1.0, color="#8A8A8A", lw=0.65, ls="--", zorder=1)
    ax.set_title(day, loc="left", fontsize=8, fontweight="bold", pad=4)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_ylim(0, max(2.4, max_y * 1.22))
    ax.grid(axis="y", color=GRID, lw=0.45, zorder=0)
    ax.tick_params(axis="x", length=0)


def build_figure(
    df: pd.DataFrame,
    stat_df: pd.DataFrame,
    genes: list[str],
    days: list[str],
    out_dir: Path,
    stem: str,
    title: str,
    reference: str,
    show_stars: bool,
) -> None:
    apply_style()
    n_panels = len(days)
    n_cols = min(3, n_panels)
    n_rows = int(np.ceil(n_panels / n_cols))
    fig_w = max(3.1 * n_cols, 4.2)
    fig_h = max(2.55 * n_rows + 0.75, 3.2)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)
    axes_flat = axes.flatten()

    for idx, (ax, day) in enumerate(zip(axes_flat, days)):
        plot_day(ax, df, stat_df, day, genes, show_stars=show_stars)
        add_panel_label(ax, chr(ord("a") + idx))
    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    for row_idx in range(n_rows):
        axes[row_idx, 0].set_ylabel("relative expression vs same-day WT")

    handles, legend_labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles[:2], legend_labels[:2], loc="upper right", bbox_to_anchor=(0.975, 0.965), ncol=2)
    fig.suptitle(title, x=0.075, y=0.985, ha="left", fontsize=10, fontweight="bold")
    stat_note = " * Welch test on log2 expression, exploratory." if show_stars else ""
    fig.text(
        0.075,
        0.02,
        f"Bars show mean + SEM; points are biological replicates. Values are {reference}-normalized and calibrated to the same-day WT mean for each gene.{stat_note}",
        ha="left",
        va="bottom",
        fontsize=6,
        color="#4D4D4D",
    )
    fig.subplots_adjust(left=0.075, right=0.975, top=0.89, bottom=0.17, wspace=0.30, hspace=0.58)

    base = out_dir / stem
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_qa(df: pd.DataFrame, genes: list[str], days: list[str], args: argparse.Namespace) -> None:
    qc_lines = [
        "Figure: stage-wise WT/KO grouped qPCR bars.",
        "Backend: Python/matplotlib.",
        f"Input analysis directory: {args.analysis_dir}",
        f"Source rel_dir: {args.rel_dir}",
        f"Genes: {', '.join(genes)}",
        f"Days/run labels: {', '.join(days)}",
        f"Normalization: {args.reference}; each gene/day calibrated to same-day WT mean.",
        "Statistics: exploratory Welch t-test on log2 expression; no multiple-testing correction displayed.",
        f"Sample rows plotted: {len(df)}",
    ]
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "figure_qa_notes.txt").write_text("\n".join(qc_lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    genes = list(dict.fromkeys(args.genes))
    if args.include_fev and "FEV" not in genes:
        genes.append("FEV")

    df = load_and_calibrate(args.analysis_dir, args.rel_dir, genes, args.days, args.out)
    days = list(df["run_label"].cat.categories)
    stat_df = make_stats(df, genes, days, args.out)
    title = args.title or "Stage-wise WT and KO qPCR expression"
    build_figure(
        df=df,
        stat_df=stat_df,
        genes=genes,
        days=days,
        out_dir=args.out,
        stem=args.stem,
        title=title,
        reference=args.reference,
        show_stars=not args.no_stars,
    )
    write_qa(df, genes, days, args)


if __name__ == "__main__":
    main()
