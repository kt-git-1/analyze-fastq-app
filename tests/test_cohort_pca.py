from pathlib import Path

import modules.cohort_pca as cohort_pca
from tests.conftest import touch


def test_read_pca_sites_supports_vcf_and_saved_table(tmp_path):
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\trs1\tA\tG\t.\t.\t.\n"
        "chr1\t11\t.\tC\tT\t.\t.\t.\n"
        "chr1\t12\tmulti\tA\tG,T\t.\t.\t.\n"
    )

    sites = cohort_pca.read_pca_sites(vcf)

    assert [site.site_id for site in sites] == ["rs1", "chr1:11:C:T"]

    saved = cohort_pca.write_sites_table(sites, tmp_path / "pca_sites.tsv")
    assert cohort_pca.read_pca_sites(saved) == sites


def test_run_cohort_pca_generates_resumable_outputs(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    cfg.args.pca_sites = tmp_path / "sites.vcf"
    cfg.args.force = False
    cfg.args.pca_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\trs1\tA\tG\t.\t.\t.\n"
        "chr1\t20\trs2\tC\tT\t.\t.\t.\n"
        "chr1\t30\trs3\tG\tA\t.\t.\t.\n"
    )
    touch(cfg.results_dir / "S1" / "dedup" / "S1.dedup.sorted.bam")
    touch(cfg.results_dir / "S2" / "dedup" / "S2.dedup.sorted.bam")
    touch(cfg.results_dir / "S3" / "dedup" / "S3.dedup.sorted.bam")

    def fake_extract(sample_bams, sites, out_path, **kwargs):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            "sample\tsite_id\tchrom\tpos\tref\talt\tallele\tbaseq\tmapq\n"
            "S1\trs1\tchr1\t10\tA\tG\tA\t35\t60\n"
            "S1\trs2\tchr1\t20\tC\tT\tC\t35\t60\n"
            "S1\trs3\tchr1\t30\tG\tA\tG\t35\t60\n"
            "S2\trs1\tchr1\t10\tA\tG\tG\t35\t60\n"
            "S2\trs2\tchr1\t20\tC\tT\tT\t35\t60\n"
            "S2\trs3\tchr1\t30\tG\tA\tG\t35\t60\n"
            "S3\trs1\tchr1\t10\tA\tG\tG\t35\t60\n"
            "S3\trs2\tchr1\t20\tC\tT\tC\t35\t60\n"
            "S3\trs3\tchr1\t30\tG\tA\tA\t35\t60\n"
        )
        return out_path

    monkeypatch.setattr(cohort_pca, "extract_pseudohaploid_calls", fake_extract)

    outputs = cohort_pca.run_cohort_pca(cfg, ["S1", "S2", "S3"])

    assert outputs["raw_calls"].exists()
    assert outputs["filtered_matrix"].exists()
    assert outputs["eigenstrat_geno"].exists()
    assert outputs["plink_tped"].exists()
    assert outputs["pca_scores"].exists()
    assert outputs["mds"].exists()
    assert "PC1" in outputs["pca_scores"].read_text()
