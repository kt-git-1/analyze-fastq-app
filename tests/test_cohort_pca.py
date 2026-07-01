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
    cfg.args.pca_engine = "python"
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


def test_filter_matrix_applies_missing_maf_and_sex_chrom_qc(tmp_path):
    matrix = tmp_path / "matrix.tsv"
    matrix.write_text(
        "sample\trs1\trs2\trs3\trsX\n"
        "S1\t0\t0\t\t0\n"
        "S2\t1\t0\t\t1\n"
        "S3\t1\t0\t\t1\n"
        "S4\t1\t\t\t1\n"
    )
    sites = [
        cohort_pca.PCASite("chr1", 10, "rs1", "A", "G"),
        cohort_pca.PCASite("chr1", 20, "rs2", "C", "T"),
        cohort_pca.PCASite("chr1", 30, "rs3", "G", "A"),
        cohort_pca.PCASite("chrX", 40, "rsX", "A", "C"),
    ]

    out = tmp_path / "filtered.tsv"
    stats = cohort_pca._filter_matrix(
        matrix,
        out,
        max_site_missing=0.5,
        max_sample_missing=0.75,
        min_maf=0.25,
        exclude_sex_chr=True,
        sites=sites,
        ploidy=1,
    )

    assert out.read_text().splitlines()[0] == "sample\trs1"
    assert stats.kept_sites == 1
    assert stats.monomorphic_removed_sites == 1
    assert stats.missingness_removed_sites == 1
    assert stats.sex_chr_removed_sites == 1


def test_plink_chr_set_args_uses_nonhuman_autosome_count():
    assert cohort_pca._plink_chr_set_args(
        [
            cohort_pca.PCASite("chr1", 10, "rs1", "A", "G"),
            cohort_pca.PCASite("chr31", 20, "rs31", "C", "T"),
            cohort_pca.PCASite("chrX", 30, "rsX", "G", "A"),
        ]
    ) == ["--chr-set", "31"]
    assert cohort_pca._plink_chr_set_args(
        [cohort_pca.PCASite("chr22", 10, "rs22", "A", "G")]
    ) == []


def test_run_cohort_pca_passes_extraction_qc(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    cfg.args.pca_sites = tmp_path / "sites.vcf"
    cfg.args.pca_engine = "python"
    cfg.args.pca_min_mapq = 12
    cfg.args.pca_min_baseq = 13
    cfg.args.pca_trim_ends = 4
    cfg.args.pca_max_sample_missing = 0.9
    cfg.args.pca_max_site_missing = 0.9
    cfg.args.pca_min_maf = 0.0
    cfg.args.pca_exclude_sex_chr = False
    cfg.args.pca_ld_window = 50
    cfg.args.pca_ld_step = 5
    cfg.args.pca_ld_r2 = 0.2
    cfg.args.pca_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\trs1\tA\tG\t.\t.\t.\n"
        "chr1\t20\trs2\tC\tT\t.\t.\t.\n"
    )
    touch(cfg.results_dir / "S1" / "dedup" / "S1.dedup.sorted.bam")
    touch(cfg.results_dir / "S2" / "dedup" / "S2.dedup.sorted.bam")
    captured = {}

    def fake_extract(sample_bams, sites, out_path, **kwargs):
        captured.update(kwargs)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            "sample\tsite_id\tchrom\tpos\tref\talt\tallele\tbaseq\tmapq\n"
            "S1\trs1\tchr1\t10\tA\tG\tA\t35\t60\n"
            "S1\trs2\tchr1\t20\tC\tT\tC\t35\t60\n"
            "S2\trs1\tchr1\t10\tA\tG\tG\t35\t60\n"
            "S2\trs2\tchr1\t20\tC\tT\tT\t35\t60\n"
        )
        return out_path

    monkeypatch.setattr(cohort_pca, "extract_pseudohaploid_calls", fake_extract)

    cohort_pca.run_cohort_pca(cfg, ["S1", "S2"])

    assert captured == {"min_mapq": 12, "min_baseq": 13, "trim_ends": 4}


def test_eigensoft_pipeline_writes_commands_and_parfiles(monkeypatch, tmp_path):
    matrix = tmp_path / "filtered.tsv"
    matrix.write_text("sample\trs1\trs2\nS1\t0\t1\nS2\t1\t0\nS3\t1\t1\n")
    sites = [
        cohort_pca.PCASite("chr1", 10, "rs1", "A", "G"),
        cohort_pca.PCASite("chr1", 20, "rs2", "C", "T"),
    ]
    commands = []

    def fake_run(cmd, cwd=None):
        commands.append(cmd)
        if cmd[0].endswith("plink") and "--indep-pairwise" in cmd:
            prune_prefix = Path(cmd[cmd.index("--out") + 1])
            prune_prefix.parent.mkdir(parents=True, exist_ok=True)
            Path(str(prune_prefix) + ".prune.in").write_text("rs1\nrs2\n")
        if cmd[0].endswith("smartpca"):
            pca_dir = tmp_path / "cohort" / "pca"
            pca_dir.mkdir(parents=True, exist_ok=True)
            (pca_dir / "cohort.evec").write_text("#eigvals: 2 1\nS1 0.1 0.2 Unknown\nS2 -0.1 0.0 Unknown\n")
            (pca_dir / "cohort.eval").write_text("2\n1\n")

    monkeypatch.setattr(cohort_pca, "_run_command", fake_run)

    eigenstrat_files, plink_files, scores, variance, mds, pruned = cohort_pca.run_eigensoft_pca(
        matrix,
        sites,
        tmp_path / "cohort",
        cohort_pca.PCAQCConfig(max_sample_missing=0.8, max_site_missing=0.7, min_maf=0.1),
    )

    assert commands[0][:3] == [str(cohort_pca.PLINK_BIN), "--tfile", str(tmp_path / "cohort" / "plink" / "cohort")]
    assert commands[1][commands[1].index("--mind") + 1] == "0.8"
    assert commands[1][commands[1].index("--geno") + 1] == "0.7"
    assert commands[1][commands[1].index("--maf") + 1] == "0.1"
    assert commands[2][0] == str(cohort_pca.PLINK_BIN)
    assert commands[3][0] == str(cohort_pca.PLINK_BIN)
    assert commands[4][0] == str(cohort_pca.CONVERTF_BIN)
    assert commands[5][0] == str(cohort_pca.SMARTPCA_BIN)
    assert "genotypename:" in (tmp_path / "cohort" / "eigensoft" / "convertf.par").read_text()
    assert "evecoutname:" in (tmp_path / "cohort" / "eigensoft" / "smartpca.par").read_text()
    assert scores.read_text().splitlines()[0] == "sample\tPC1\tPC2"
    assert variance.read_text().splitlines()[1] == "PC1\t2\t0.66666667"
    assert mds.exists()
    assert pruned == 2
    assert eigenstrat_files[0].name == "cohort.geno"
    assert plink_files[0].name == "cohort.tped"


def test_eigensoft_pipeline_sets_numchrom_for_nonhuman_autosomes(monkeypatch, tmp_path):
    matrix = tmp_path / "filtered.tsv"
    matrix.write_text("sample\trs1\trs31\nS1\t0\t1\nS2\t1\t0\nS3\t1\t1\n")
    sites = [
        cohort_pca.PCASite("chr1", 10, "rs1", "A", "G"),
        cohort_pca.PCASite("chr31", 20, "rs31", "C", "T"),
    ]

    def fake_run(cmd, cwd=None):
        if cmd[0].endswith("plink") and "--indep-pairwise" in cmd:
            prune_prefix = Path(cmd[cmd.index("--out") + 1])
            prune_prefix.parent.mkdir(parents=True, exist_ok=True)
            Path(str(prune_prefix) + ".prune.in").write_text("rs1\nrs31\n")
        if cmd[0].endswith("smartpca"):
            pca_dir = tmp_path / "cohort" / "pca"
            pca_dir.mkdir(parents=True, exist_ok=True)
            (pca_dir / "cohort.evec").write_text("#eigvals: 2 1\nS1 0.1 0.2 Unknown\nS2 -0.1 0.0 Unknown\n")
            (pca_dir / "cohort.eval").write_text("2\n1\n")

    monkeypatch.setattr(cohort_pca, "_run_command", fake_run)

    cohort_pca.run_eigensoft_pca(
        matrix,
        sites,
        tmp_path / "cohort",
        cohort_pca.PCAQCConfig(),
    )

    assert "numchrom: 31" in (tmp_path / "cohort" / "eigensoft" / "convertf.par").read_text()
    assert "numchrom: 31" in (tmp_path / "cohort" / "eigensoft" / "smartpca.par").read_text()


def test_eigensoft_pipeline_resumes_existing_outputs(monkeypatch, tmp_path):
    cohort_dir = tmp_path / "cohort"
    matrix = tmp_path / "filtered.tsv"
    matrix.write_text("sample\trs1\trs2\nS1\t0\t1\nS2\t1\t0\n")
    sites = [
        cohort_pca.PCASite("chr1", 10, "rs1", "A", "G"),
        cohort_pca.PCASite("chr1", 20, "rs2", "C", "T"),
    ]
    for path in (
        cohort_dir / "plink" / "cohort.tped",
        cohort_dir / "plink" / "cohort.tfam",
        cohort_dir / "plink" / "cohort.bed",
        cohort_dir / "plink" / "cohort.bim",
        cohort_dir / "plink" / "cohort.fam",
        cohort_dir / "plink" / "cohort.qc.bed",
        cohort_dir / "plink" / "cohort.qc.bim",
        cohort_dir / "plink" / "cohort.qc.fam",
        cohort_dir / "plink" / "cohort.prune.prune.in",
        cohort_dir / "plink" / "cohort.pruned.bed",
        cohort_dir / "plink" / "cohort.pruned.bim",
        cohort_dir / "plink" / "cohort.pruned.fam",
        cohort_dir / "eigenstrat" / "cohort.geno",
        cohort_dir / "eigenstrat" / "cohort.snp",
        cohort_dir / "eigenstrat" / "cohort.ind",
        cohort_dir / "pca" / "cohort.evec",
        cohort_dir / "pca" / "cohort.eval",
        cohort_dir / "pca" / "pca_scores.tsv",
        cohort_dir / "pca" / "pca_variance.tsv",
        cohort_dir / "pca" / "mds.tsv",
    ):
        touch(path)
    (cohort_dir / "plink" / "cohort.prune.prune.in").write_text("rs1\nrs2\n")
    commands = []
    monkeypatch.setattr(cohort_pca, "_run_command", lambda cmd, cwd=None: commands.append(cmd))

    _eigenstrat, _plink, scores, variance, mds, pruned = cohort_pca.run_eigensoft_pca(
        matrix,
        sites,
        cohort_dir,
        cohort_pca.PCAQCConfig(),
    )

    assert commands == []
    assert scores == cohort_dir / "pca" / "pca_scores.tsv"
    assert variance == cohort_dir / "pca" / "pca_variance.tsv"
    assert mds == cohort_dir / "pca" / "mds.tsv"
    assert pruned == 2


def test_run_cohort_pca_modern_uses_diploid_bam_calls(monkeypatch, make_config, tmp_path):
    cfg = make_config(data_type="modern")
    cfg.args.pca_sites = tmp_path / "sites.vcf"
    cfg.args.pca_engine = "python"
    cfg.args.pca_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\trs1\tA\tG\t.\t.\t.\n"
        "chr1\t20\trs2\tC\tT\t.\t.\t.\n"
    )
    for sample in ("S1", "S2", "S3"):
        touch(cfg.results_dir / sample / "dedup" / ("%s.dedup.sorted.bam" % sample))

    def fake_extract(sample_bams, sites, out_path, **kwargs):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            "sample\tsite_id\tchrom\tpos\tref\talt\tref_count\talt_count\tdepth\tdosage\n"
            "S1\trs1\tchr1\t10\tA\tG\t10\t0\t10\t0\n"
            "S1\trs2\tchr1\t20\tC\tT\t5\t5\t10\t1\n"
            "S2\trs1\tchr1\t10\tA\tG\t5\t5\t10\t1\n"
            "S2\trs2\tchr1\t20\tC\tT\t0\t10\t10\t2\n"
            "S3\trs1\tchr1\t10\tA\tG\t0\t10\t10\t2\n"
            "S3\trs2\tchr1\t20\tC\tT\t10\t0\t10\t0\n"
        )
        return out_path

    monkeypatch.setattr(cohort_pca, "extract_modern_diploid_calls", fake_extract)

    outputs = cohort_pca.run_cohort_pca(cfg, ["S1", "S2", "S3"])

    assert outputs["raw_calls"].name == "modern_genotype_calls.tsv"
    assert outputs["matrix"].read_text().splitlines()[1:] == [
        "S1\t0\t1",
        "S2\t1\t2",
        "S3\t2\t0",
    ]
    assert "data_type\tmodern" in outputs["qc_summary"].read_text()
    assert outputs["pca_scores"].exists()


def test_run_cohort_pca_auto_infers_ancient_from_short_reads(monkeypatch, make_config, tmp_path):
    cfg = make_config(data_type="auto")
    cfg.args.pca_sites = tmp_path / "sites.vcf"
    cfg.args.pca_engine = "python"
    cfg.args.pca_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\trs1\tA\tG\t.\t.\t.\n"
        "chr1\t20\trs2\tC\tT\t.\t.\t.\n"
    )
    for sample in ("S1", "S2"):
        touch(cfg.results_dir / sample / "dedup" / ("%s.dedup.sorted.bam" % sample))

    monkeypatch.setattr(cohort_pca, "_median_mapped_read_length", lambda bam: 65.0)

    def fake_extract(sample_bams, sites, out_path, **kwargs):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            "sample\tsite_id\tchrom\tpos\tref\talt\tallele\tbaseq\tmapq\n"
            "S1\trs1\tchr1\t10\tA\tG\tA\t35\t60\n"
            "S1\trs2\tchr1\t20\tC\tT\tC\t35\t60\n"
            "S2\trs1\tchr1\t10\tA\tG\tG\t35\t60\n"
            "S2\trs2\tchr1\t20\tC\tT\tT\t35\t60\n"
        )
        return out_path

    monkeypatch.setattr(cohort_pca, "extract_pseudohaploid_calls", fake_extract)

    outputs = cohort_pca.run_cohort_pca(cfg, ["S1", "S2"])

    assert outputs["raw_calls"].name == "pseudohaploid_raw_calls.tsv"
    assert "requested_data_type\tauto" in outputs["qc_summary"].read_text()
    assert "data_type\tancient" in outputs["qc_summary"].read_text()
    assert (cfg.results_dir / "cohort" / "auto_data_type_summary.tsv").exists()


def test_run_cohort_pca_auto_infers_modern_from_long_reads(monkeypatch, make_config, tmp_path):
    cfg = make_config(data_type="auto")
    cfg.args.pca_sites = tmp_path / "sites.vcf"
    cfg.args.pca_engine = "python"
    cfg.args.pca_sites.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\trs1\tA\tG\t.\t.\t.\n"
        "chr1\t20\trs2\tC\tT\t.\t.\t.\n"
    )
    for sample in ("S1", "S2", "S3"):
        touch(cfg.results_dir / sample / "dedup" / ("%s.dedup.sorted.bam" % sample))

    monkeypatch.setattr(cohort_pca, "_median_mapped_read_length", lambda bam: 150.0)

    def fake_extract(sample_bams, sites, out_path, **kwargs):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            "sample\tsite_id\tchrom\tpos\tref\talt\tref_count\talt_count\tdepth\tdosage\n"
            "S1\trs1\tchr1\t10\tA\tG\t10\t0\t10\t0\n"
            "S1\trs2\tchr1\t20\tC\tT\t5\t5\t10\t1\n"
            "S2\trs1\tchr1\t10\tA\tG\t5\t5\t10\t1\n"
            "S2\trs2\tchr1\t20\tC\tT\t0\t10\t10\t2\n"
            "S3\trs1\tchr1\t10\tA\tG\t0\t10\t10\t2\n"
            "S3\trs2\tchr1\t20\tC\tT\t10\t0\t10\t0\n"
        )
        return out_path

    monkeypatch.setattr(cohort_pca, "extract_modern_diploid_calls", fake_extract)

    outputs = cohort_pca.run_cohort_pca(cfg, ["S1", "S2", "S3"])

    assert outputs["raw_calls"].name == "modern_genotype_calls.tsv"
    assert "requested_data_type\tauto" in outputs["qc_summary"].read_text()
    assert "data_type\tmodern" in outputs["qc_summary"].read_text()


def test_diploid_dosage_from_counts():
    assert cohort_pca._diploid_dosage_from_counts(10, 0, 10) == "0"
    assert cohort_pca._diploid_dosage_from_counts(5, 5, 10) == "1"
    assert cohort_pca._diploid_dosage_from_counts(0, 10, 10) == "2"
    assert cohort_pca._diploid_dosage_from_counts(0, 0, 0) == ""


def test_export_plink_text_handles_modern_heterozygotes(tmp_path):
    matrix = tmp_path / "matrix.tsv"
    matrix.write_text("sample\trs1\nS1\t0\nS2\t1\nS3\t2\n")
    sites = [cohort_pca.PCASite("chr1", 10, "rs1", "A", "G")]

    tped, _tfam, _sample_map = cohort_pca.export_plink_text(matrix, sites, tmp_path / "plink")

    assert tped.read_text().strip().split("\t") == [
        "1",
        "rs1",
        "0",
        "10",
        "A",
        "A",
        "A",
        "G",
        "G",
        "G",
    ]


def test_eigensoft_pipeline_uses_short_plink_ids_and_restores_sample_names(monkeypatch, tmp_path):
    matrix = tmp_path / "filtered.tsv"
    long_sample = "PE-AncientHorses-01_S1"
    matrix.write_text("sample\trs1\trs2\n%s\t0\t1\nZYJ2_S9\t1\t0\n" % long_sample)
    sites = [
        cohort_pca.PCASite("chr1", 10, "rs1", "A", "G"),
        cohort_pca.PCASite("chr1", 20, "rs2", "C", "T"),
    ]

    def fake_run(cmd, cwd=None):
        if cmd[0].endswith("plink") and "--indep-pairwise" in cmd:
            prune_prefix = Path(cmd[cmd.index("--out") + 1])
            prune_prefix.parent.mkdir(parents=True, exist_ok=True)
            Path(str(prune_prefix) + ".prune.in").write_text("rs1\nrs2\n")
        if cmd[0].endswith("smartpca"):
            pca_dir = tmp_path / "cohort" / "pca"
            pca_dir.mkdir(parents=True, exist_ok=True)
            (pca_dir / "cohort.evec").write_text(
                "#eigvals: 2 1\n"
                "I000001 0.1 0.2 Unknown\n"
                "I000002 -0.1 0.0 Unknown\n"
            )
            (pca_dir / "cohort.eval").write_text("2\n1\n")

    monkeypatch.setattr(cohort_pca, "_run_command", fake_run)

    _eigenstrat, _plink, scores, _variance, _mds, _pruned = cohort_pca.run_eigensoft_pca(
        matrix,
        sites,
        tmp_path / "cohort",
        cohort_pca.PCAQCConfig(),
    )

    tfam_lines = (tmp_path / "cohort" / "plink" / "cohort.tfam").read_text().splitlines()
    assert tfam_lines[0].split()[:2] == ["0", "I000001"]
    assert long_sample not in tfam_lines[0]
    assert "I000001\t%s" % long_sample in (
        tmp_path / "cohort" / "plink" / "sample_id_map.tsv"
    ).read_text()
    assert scores.read_text().splitlines()[1].startswith("%s\t" % long_sample)
