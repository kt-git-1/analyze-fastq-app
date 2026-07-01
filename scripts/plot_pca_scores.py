#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path
from itertools import combinations


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


def _auto_group_value(sample, mode):
    if mode == "prefix":
        match = re.match(r"^(.+?)-\d+_", sample)
        if match:
            return match.group(1)
        return sample.split("_", 1)[0]
    raise ValueError("unknown auto group mode: %s" % mode)


def _auto_group_values(samples, mode):
    return [_auto_group_value(sample, mode) for sample in samples]


def _group_colors(group_count):
    palette = [
        "#0072B2",  # blue
        "#D55E00",  # vermillion
        "#009E73",  # green
        "#CC79A7",  # reddish purple
        "#E69F00",  # orange
        "#56B4E9",  # sky blue
        "#F0E442",  # yellow
        "#000000",  # black
    ]
    if group_count <= len(palette):
        return palette[:group_count]

    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("tab20")
    return [cmap(idx % 20) for idx in range(group_count)]


def _scatter_by_group(ax, x_values, y_values, groups, legend_title):
    if not groups:
        ax.scatter(x_values, y_values, s=54, alpha=0.9, edgecolors="white", linewidths=0.7)
        return

    group_order = []
    seen = set()
    for group in groups:
        if group not in seen:
            seen.add(group)
            group_order.append(group)
    colors = _group_colors(len(group_order))
    for idx, group in enumerate(group_order):
        indices = [row_idx for row_idx, value in enumerate(groups) if value == group]
        ax.scatter(
            [x_values[row_idx] for row_idx in indices],
            [y_values[row_idx] for row_idx in indices],
            s=58,
            alpha=0.92,
            edgecolors="white",
            linewidths=0.75,
            label=group,
            color=colors[idx],
        )
    ax.legend(title=legend_title, loc="best", fontsize=8, title_fontsize=9, frameon=True)


def _label_indices(x_values, y_values, label_top_n):
    if label_top_n is None or label_top_n <= 0:
        return []
    ranked = sorted(
        range(len(x_values)),
        key=lambda idx: (x_values[idx] * x_values[idx]) + (y_values[idx] * y_values[idx]),
        reverse=True,
    )
    return ranked[:label_top_n]


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
    auto_group=None,
    label_top_n=None,
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
    if auto_group:
        groups = _auto_group_values(samples, auto_group)
        legend_title = "group"
    else:
        groups = _group_values(samples, metadata, color_by) if color_by else []
        legend_title = color_by or "group"

    fig, ax = plt.subplots(figsize=(8, 6))
    _scatter_by_group(ax, x_values, y_values, groups, legend_title)
    ax.axhline(0, color="#999999", linewidth=0.7, zorder=0)
    ax.axvline(0, color="#999999", linewidth=0.7, zorder=0)
    ax.set_xlabel(_axis_label(x_component, variance))
    ax.set_ylabel(_axis_label(y_component, variance))
    ax.set_title("PCA scatter")
    ax.grid(True, color="#dddddd", linewidth=0.6, alpha=0.7)

    label_target_indices = range(len(samples)) if label_samples else _label_indices(x_values, y_values, label_top_n)
    for idx in label_target_indices:
        ax.annotate(
            samples[idx],
            (x_values[idx], y_values[idx]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=7,
        )

    fig.tight_layout()
    png_path = Path(str(out_prefix) + ".png")
    pdf_path = Path(str(out_prefix) + ".pdf")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def plot_pca_pairs(
    scores_path,
    variance_path,
    out_dir,
    pc_count,
    label_samples=False,
    metadata_path=None,
    color_by=None,
    sample_column="sample",
    auto_group=None,
    label_top_n=None,
):
    outputs = []
    for left, right in combinations(range(1, pc_count + 1), 2):
        x_component = "PC%d" % left
        y_component = "PC%d" % right
        out_prefix = out_dir / ("pca_%s_%s" % (x_component, y_component))
        outputs.append(
            plot_pca(
                scores_path,
                variance_path,
                out_prefix,
                x_component,
                y_component,
                label_samples,
                metadata_path=metadata_path,
                color_by=color_by,
                sample_column=sample_column,
                auto_group=auto_group,
                label_top_n=label_top_n,
            )
        )
    return outputs


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
        help="出力prefix。--all-pairs指定時は出力ディレクトリとして扱います。省略時はpca_scores.tsvと同じ場所に出力します。",
    )
    parser.add_argument("--x", default="PC1", help="X軸のPC列名")
    parser.add_argument("--y", default="PC2", help="Y軸のPC列名")
    parser.add_argument("--label-samples", action="store_true", help="点にサンプル名を表示します。")
    parser.add_argument("--label-top-n", type=int, default=None, help="原点から遠い上位Nサンプルだけにラベルを表示します。")
    parser.add_argument("--metadata", type=Path, default=None, help="sample情報を持つTSV。sample列でpca_scores.tsvと結合します。")
    parser.add_argument("--color-by", default=None, help="metadata TSV内の色分け列名。例: group, population, batch")
    parser.add_argument("--sample-column", default="sample", help="metadata TSV内のsample名列。デフォルト: sample")
    parser.add_argument(
        "--auto-group",
        choices=["prefix"],
        default=None,
        help="metadataなしでsample名から自動色分けします。prefixは PE-AncientHorses-01_S1 -> PE-AncientHorses、ZYJ2_S1 -> ZYJ2 のように解釈します。",
    )
    parser.add_argument(
        "--all-pairs",
        type=int,
        default=None,
        help="PC1から指定PCまでの全ペアをまとめて出力します。例: --all-pairs 4 はPC1-4の6ペアを出力します。",
    )
    args = parser.parse_args()

    variance = args.variance or args.scores.parent / "pca_variance.tsv"
    color_by = args.color_by or ("group" if args.metadata else None)
    if args.all_pairs:
        out_dir = args.out_prefix or args.scores.parent
        for png_path, pdf_path in plot_pca_pairs(
            args.scores,
            variance,
            out_dir,
            args.all_pairs,
            label_samples=args.label_samples,
            metadata_path=args.metadata,
            color_by=color_by,
            sample_column=args.sample_column,
            auto_group=args.auto_group,
            label_top_n=args.label_top_n,
        ):
            print("PNG: %s" % png_path)
            print("PDF: %s" % pdf_path)
        return

    out_prefix = args.out_prefix
    if out_prefix is None:
        out_prefix = args.scores.parent / ("pca_%s_%s" % (args.x, args.y))
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
        auto_group=args.auto_group,
        label_top_n=args.label_top_n,
    )
    print("PNG: %s" % png_path)
    print("PDF: %s" % pdf_path)


if __name__ == "__main__":
    main()
