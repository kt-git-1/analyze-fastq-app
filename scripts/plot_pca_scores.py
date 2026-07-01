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


def _read_metadata(path, sample_column):
    if not path:
        return {}
    with path.open(errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if sample_column not in (reader.fieldnames or []):
            raise ValueError("metadata TSV に sample列がありません: %s" % sample_column)
        return {row[sample_column]: row for row in reader if row.get(sample_column)}


def _group_values(samples, metadata, color_by):
    if not metadata:
        return []
    groups = []
    for sample in samples:
        row = metadata.get(sample)
        value = row.get(color_by, "") if row else ""
        groups.append(value or "未指定")
    return groups


def _scatter_by_group(ax, x_values, y_values, groups, legend_title):
    if not groups:
        ax.scatter(x_values, y_values, s=42, alpha=0.85, edgecolors="black", linewidths=0.35)
        return

    import matplotlib.pyplot as plt

    group_order = []
    seen = set()
    for group in groups:
        if group not in seen:
            seen.add(group)
            group_order.append(group)
    cmap = plt.get_cmap("tab20")
    color_count = max(1, min(20, len(group_order)))
    for idx, group in enumerate(group_order):
        indices = [row_idx for row_idx, value in enumerate(groups) if value == group]
        color = cmap(idx % color_count)
        ax.scatter(
            [x_values[row_idx] for row_idx in indices],
            [y_values[row_idx] for row_idx in indices],
            s=46,
            alpha=0.88,
            edgecolors="black",
            linewidths=0.35,
            label=group,
            color=color,
        )
    ax.legend(title=legend_title, loc="best", fontsize=8, title_fontsize=9, frameon=True)


def plot_pca(
    scores_path,
    variance_path,
    out_prefix,
    x_component,
    y_component,
    label_samples,
    metadata_path=None,
    color_by=None,
    sample_column="sample",
):
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
    metadata = _read_metadata(metadata_path, sample_column)
    if color_by and metadata:
        first_row = next(iter(metadata.values()))
        if color_by not in first_row:
            raise ValueError("metadata TSV に色分け列がありません: %s" % color_by)
    missing = [component for component in (x_component, y_component) if component not in rows[0]]
    if missing:
        raise ValueError("pca_scores.tsv に列がありません: %s" % ", ".join(missing))

    samples = [row["sample"] for row in rows]
    x_values = [float(row[x_component]) for row in rows]
    y_values = [float(row[y_component]) for row in rows]
    groups = _group_values(samples, metadata, color_by) if color_by else []

    fig, ax = plt.subplots(figsize=(8, 6))
    _scatter_by_group(ax, x_values, y_values, groups, color_by or "group")
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
    parser.add_argument("--metadata", type=Path, default=None, help="sample情報を持つTSV。sample列でpca_scores.tsvと結合します。")
    parser.add_argument("--color-by", default=None, help="metadata TSV内の色分け列名。例: group, population, batch")
    parser.add_argument("--sample-column", default="sample", help="metadata TSV内のsample名列。デフォルト: sample")
    args = parser.parse_args()

    out_prefix = args.out_prefix
    if out_prefix is None:
        out_prefix = args.scores.parent / ("pca_%s_%s" % (args.x, args.y))
    variance = args.variance or args.scores.parent / "pca_variance.tsv"
    color_by = args.color_by or ("group" if args.metadata else None)
    png_path, pdf_path = plot_pca(
        args.scores,
        variance,
        out_prefix,
        args.x,
        args.y,
        args.label_samples,
        metadata_path=args.metadata,
        color_by=color_by,
        sample_column=args.sample_column,
    )
    print("PNG: %s" % png_path)
    print("PDF: %s" % pdf_path)


if __name__ == "__main__":
    main()
