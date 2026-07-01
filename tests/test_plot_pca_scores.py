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
