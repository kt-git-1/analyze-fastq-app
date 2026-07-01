#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def _read_scores(path):
    with path.open(errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError("pca_scores.tsv has no rows: %s" % path)
    return rows


def _read_variance(path):
    if not path or not path.exists():
        return {}
    values = {}
    with path.open(errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            component = row.get("component")
            explained = row.get("explained_variance")
            if not component or explained in (None, ""):
                continue
            try:
                values[component] = float(explained)
            except ValueError:
                continue
    return values


def _axis_label(component, variance):
    if component in variance:
        return "%s (%.1f%%)" % (component, variance[component] * 100.0)
    return component


def plot_pca(scores_path, variance_path, out_prefix, x_component, y_component, label_samples):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib が見つかりません。conda/pipで matplotlib を入れてから再実行してください。"
        ) from exc

    rows = _read_scores(scores_path)
    variance = _read_variance(variance_path)
    missing = [component for component in (x_component, y_component) if component not in rows[0]]
    if missing:
        raise ValueError("pca_scores.tsv に列がありません: %s" % ", ".join(missing))

    samples = [row["sample"] for row in rows]
    x_values = [float(row[x_component]) for row in rows]
    y_values = [float(row[y_component]) for row in rows]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x_values, y_values, s=42, alpha=0.85, edgecolors="black", linewidths=0.35)
    ax.axhline(0, color="#999999", linewidth=0.7, zorder=0)
    ax.axvline(0, color="#999999", linewidth=0.7, zorder=0)
    ax.set_xlabel(_axis_label(x_component, variance))
    ax.set_ylabel(_axis_label(y_component, variance))
    ax.set_title("PCA scatter")
    ax.grid(True, color="#dddddd", linewidth=0.6, alpha=0.7)

    if label_samples:
        for sample, x_value, y_value in zip(samples, x_values, y_values):
            ax.annotate(sample, (x_value, y_value), xytext=(3, 3), textcoords="offset points", fontsize=7)

    fig.tight_layout()
    png_path = Path(str(out_prefix) + ".png")
    pdf_path = Path(str(out_prefix) + ".pdf")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def main():
    parser = argparse.ArgumentParser(description="pca_scores.tsv から PCA 散布図 PNG/PDF を作成します。")
    parser.add_argument("scores", type=Path, help="pca_scores.tsv")
    parser.add_argument(
        "--variance",
        type=Path,
        default=None,
        help="pca_variance.tsv。指定すると軸ラベルに説明率を表示します。",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=None,
        help="出力prefix。省略時は pca_scores.tsv と同じ場所の pca_PC1_PC2",
    )
    parser.add_argument("--x", default="PC1", help="X軸のPC列名")
    parser.add_argument("--y", default="PC2", help="Y軸のPC列名")
    parser.add_argument("--label-samples", action="store_true", help="点にサンプル名を表示します。")
    args = parser.parse_args()

    out_prefix = args.out_prefix
    if out_prefix is None:
        out_prefix = args.scores.parent / ("pca_%s_%s" % (args.x, args.y))
    variance = args.variance or args.scores.parent / "pca_variance.tsv"
    png_path, pdf_path = plot_pca(args.scores, variance, out_prefix, args.x, args.y, args.label_samples)
    print("PNG: %s" % png_path)
    print("PDF: %s" % pdf_path)


if __name__ == "__main__":
    main()
