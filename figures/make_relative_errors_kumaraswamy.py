"""Figura oficial do capitulo: erros relativos vs ordem de quadratura N.

Le `outputs/csv/chapter_campaign_multi_results.csv` e gera PDF + PNG em
`outputs/figures/`.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = PROJECT_ROOT / "outputs" / "csv" / "chapter_campaign_multi_results.csv"
OUT_DIR = PROJECT_ROOT / "outputs" / "figures"
OUT_PDF = OUT_DIR / "relative_errors_kumaraswamy.pdf"
OUT_PNG = OUT_DIR / "relative_errors_kumaraswamy.png"

CASE_ORDER: list[str] = ["Kum(2,8)", "Kum(2,2.41)", "Kum(4,1.5)"]
CASE_TITLES: dict[str, str] = {
    "Kum(2,8)": "Kumaraswamy (2,8)",
    "Kum(2,2.41)": "Kumaraswamy (2,2.41)",
    "Kum(4,1.5)": "Kumaraswamy (4,1.5)",
}
METRICS: list[tuple[str, str]] = [
    ("rel_mean", "Mean"),
    ("rel_var", "Variance"),
    ("rel_q50", r"$Q_{0.50}$"),
]
LINE_WIDTH = 2.5
MARKER_SIZE = 7
Y_FLOOR = 1e-16


def _load_campaign_csv(path: Path) -> dict[str, dict[str, list[float]]]:
    if not path.is_file():
        raise FileNotFoundError(f"CSV da campanha nao encontrado: {path}")

    needed = {"label", "n", "rel_mean", "rel_var", "rel_q50"}
    by_label: dict[str, dict[str, list[float]]] = {
        label: {"n": [], "rel_mean": [], "rel_var": [], "rel_q50": []}
        for label in CASE_ORDER
    }

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or not needed.issubset(reader.fieldnames):
            missing = needed - set(reader.fieldnames or [])
            raise ValueError(f"Colunas em falta no CSV: {sorted(missing)}")

        for row in reader:
            label = row["label"].strip()
            if label not in by_label:
                continue
            n = int(row["n"])
            by_label[label]["n"].append(n)
            for col, _ in METRICS:
                by_label[label][col].append(float(row[col]))

    for label in CASE_ORDER:
        order = np.argsort(by_label[label]["n"])
        for key in ("n", "rel_mean", "rel_var", "rel_q50"):
            arr = np.asarray(by_label[label][key], dtype=float)
            by_label[label][key] = arr[order].tolist()

        if len(by_label[label]["n"]) == 0:
            raise ValueError(f"Sem linhas para {label} em {path}")

    return by_label


def _apply_springer_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.linewidth": 0.8,
            "grid.alpha": 0.25,
            "lines.linewidth": LINE_WIDTH,
            "lines.markersize": MARKER_SIZE,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )


def make_figure(data: dict[str, dict[str, list[float]]]) -> plt.Figure:
    _apply_springer_style()

    colors = {"Mean": "#1f77b4", "Variance": "#d62728", r"$Q_{0.50}$": "#2ca02c"}
    markers = {"Mean": "o", "Variance": "s", r"$Q_{0.50}$": "^"}

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=True)

    for ax, label in zip(axes, CASE_ORDER):
        n_vals = np.asarray(data[label]["n"], dtype=int)
        for col, legend_name in METRICS:
            y = np.maximum(np.asarray(data[label][col], dtype=float), Y_FLOOR)
            ax.semilogy(
                n_vals,
                y,
                marker=markers[legend_name],
                color=colors[legend_name],
                label=legend_name,
                linewidth=LINE_WIDTH,
                markersize=MARKER_SIZE,
            )

        ax.set_title(CASE_TITLES[label])
        ax.set_xlabel("Quadrature order ($N$)")
        ax.set_xticks(n_vals)
        ax.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.25)
        ax.set_xlim(n_vals.min() - 0.4, n_vals.max() + 0.4)

    axes[0].set_ylabel("Relative error")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=3,
        frameon=True,
        framealpha=0.9,
    )

    fig.tight_layout(w_pad=1.2)
    fig.subplots_adjust(bottom=0.22)
    return fig


def main() -> int:
    data = _load_campaign_csv(CSV_PATH)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig = make_figure(data)
    fig.savefig(OUT_PDF, format="pdf")
    fig.savefig(OUT_PNG, format="png")
    plt.close(fig)

    print(OUT_PDF.resolve())
    print(OUT_PNG.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
