import importlib.util
from pathlib import Path


def _load_plot_module():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "plot_pca_scores.py"
    spec = importlib.util.spec_from_file_location("plot_pca_scores", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_read_metadata_and_group_values(tmp_path):
    plot_pca_scores = _load_plot_module()
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sample\tgroup\tbatch\n"
        "S1\tancient\tA\n"
        "S2\tmodern\tA\n"
    )

    rows = plot_pca_scores._read_metadata(metadata, "sample")
    groups = plot_pca_scores._group_values(["S1", "S2", "S3"], rows, "group")

    assert groups == ["ancient", "modern", "未指定"]


def test_auto_group_prefix_uses_sample_name_prefix():
    plot_pca_scores = _load_plot_module()

    groups = plot_pca_scores._auto_group_values(
        ["PE-AncientHorses-01_S1", "PE-AncientHorses-16_S16", "ZYJ2_S1"],
        "prefix",
    )

    assert groups == ["PE-AncientHorses", "PE-AncientHorses", "ZYJ2"]


def test_group_colors_are_distinct_for_two_groups():
    plot_pca_scores = _load_plot_module()

    assert plot_pca_scores._group_colors(2) == ["#0072B2", "#D55E00"]


def test_label_indices_selects_farthest_points():
    plot_pca_scores = _load_plot_module()

    assert plot_pca_scores._label_indices([0.0, 2.0, 1.0], [0.0, 0.0, 2.0], 2) == [2, 1]


def test_plot_pca_pairs_generates_all_pairs(monkeypatch, tmp_path):
    plot_pca_scores = _load_plot_module()
    calls = []

    def fake_plot_pca(scores, variance, out_prefix, x_component, y_component, label_samples, **kwargs):
        calls.append((x_component, y_component, out_prefix.name, kwargs.get("auto_group"), kwargs.get("label_top_n")))
        return tmp_path / ("%s.png" % out_prefix.name), tmp_path / ("%s.pdf" % out_prefix.name)

    monkeypatch.setattr(plot_pca_scores, "plot_pca", fake_plot_pca)

    outputs = plot_pca_scores.plot_pca_pairs(
        tmp_path / "pca_scores.tsv",
        tmp_path / "pca_variance.tsv",
        tmp_path / "plots",
        4,
        auto_group="prefix",
        label_top_n=5,
    )

    assert [call[:2] for call in calls] == [
        ("PC1", "PC2"),
        ("PC1", "PC3"),
        ("PC1", "PC4"),
        ("PC2", "PC3"),
        ("PC2", "PC4"),
        ("PC3", "PC4"),
    ]
    assert calls[0][2:] == ("pca_PC1_PC2", "prefix", 5)
    assert len(outputs) == 6
