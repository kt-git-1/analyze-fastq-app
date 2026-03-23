import sys
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests
from tqdm import tqdm

from config import (
    PipelineConfig,
    parse_args,
    setup_logging,
    cleanup_intermediate_file,
)
from modules.fastq_parser import group_fastqs_by_run, parse_fastq_general, merge_lanes_by_cat
from modules.ena_downloader import ENADownloader, verify_file_md5
from modules.bwa_mapper import BWAMapper
from modules.softclipper import SoftClipper
from modules.bam_processor import BAMProcessor
from modules.analyzers import (
    MapDamageAnalyzer,
    QualimapAnalyzer,
    HaplotypeCaller,
)

logger = logging.getLogger(__name__)


# ============================================================
# サンプル単位の解析関数
# ============================================================

def process_sample(
    sample_acc: str,
    runs: List[Tuple[str, List[Path]]],
    config: PipelineConfig,
) -> Tuple[str, bool, str]:
    """
    1 サンプル分の解析パイプラインを実行する。

    各ランを独立にマッピング → ソフトクリップ → CleanSam し、
    その後サンプル単位でマージ → 重複除去 → QC → VCF を行う。

    Returns
    -------
    tuple[str, bool, str]
        (sample_acc, 成功したか, 失敗時のステップ名 or "")
    """
    done_flag = config.results_dir / sample_acc / ".done"
    force: bool = getattr(config.args, "force", False)

    if done_flag.exists() and not force:
        logger.info("スキップ (完了済み): %s", sample_acc)
        return sample_acc, True, ""

    bwa_mapper = BWAMapper(config)
    softclipper = SoftClipper(config)
    bam_processor = BAMProcessor(config)
    mapdamage_analyzer = MapDamageAnalyzer(config)
    qualimap_analyzer = QualimapAnalyzer(config)
    haplotypecaller = HaplotypeCaller(config)

    # ステップ数: ラン単位 3 ステップ × N ラン + サンプル単位 4 ステップ
    total_steps = 3 * len(runs) + 4
    current_step = 0

    def _log_progress(step_name: str) -> None:
        nonlocal current_step
        current_step += 1
        pct = int(current_step / total_steps * 100)
        logger.info(
            "[%s] %d%% (%d/%d) %s",
            sample_acc, pct, current_step, total_steps, step_name,
        )

    logger.info("サンプルの解析を開始します: %s (%d ラン, %d ステップ)", sample_acc, len(runs), total_steps)

    # ------------------------------------------------------------------
    # Phase 1: ラン単位の処理
    # ------------------------------------------------------------------
    run_clean_bams: List[Path] = []
    failed_runs: List[str] = []

    for run_idx, (run_id, fastq_files) in enumerate(runs, 1):
        # 1. BWA マッピング
        _log_progress(f"BWA mapping ({run_id}) [{run_idx}/{len(runs)}]")
        bam_file = bwa_mapper.run_mapping_pipeline(sample_acc, run_id, fastq_files)
        if not bam_file:
            logger.error("マッピング失敗 → ラン %s をスキップ", run_id)
            failed_runs.append(run_id)
            current_step += 2
            continue

        # 2. Soft clipping
        _log_progress(f"Soft clipping ({run_id})")
        softclipped_bam = softclipper.run_softclipping(sample_acc, run_id, bam_file)
        if not softclipped_bam:
            logger.error("Soft clipping 失敗 → ラン %s をスキップ", run_id)
            failed_runs.append(run_id)
            current_step += 1
            continue

        cleanup_intermediate_file(bam_file, logger)

        # 3. CleanSam (ラン単位)
        _log_progress(f"CleanSam ({run_id})")
        try:
            clean_bam = bam_processor.process_run_bam(sample_acc, run_id, softclipped_bam)
        except Exception:
            logger.exception("CleanSam 失敗 → ラン %s をスキップ", run_id)
            failed_runs.append(run_id)
            continue

        cleanup_intermediate_file(softclipped_bam, logger)

        if clean_bam:
            run_clean_bams.append(clean_bam)
        else:
            failed_runs.append(run_id)

    if failed_runs:
        logger.warning(
            "失敗したラン: %s (%d/%d)",
            ", ".join(failed_runs),
            len(failed_runs),
            len(runs),
        )

    if not run_clean_bams:
        logger.error("有効なランがありません: %s", sample_acc)
        return sample_acc, False, "all runs failed"

    # ------------------------------------------------------------------
    # Phase 2: サンプル単位の処理
    # ------------------------------------------------------------------

    # 4. マージ + 重複除去
    _log_progress("Merge + MarkDuplicates")
    try:
        dedup_bam = bam_processor.merge_and_dedup(sample_acc, run_clean_bams)
    except Exception:
        logger.exception("merge/dedup に失敗しました: %s", sample_acc)
        return sample_acc, False, "merge/dedup"

    if not dedup_bam:
        return sample_acc, False, "merge/dedup"

    for bam in run_clean_bams:
        cleanup_intermediate_file(bam, logger)

    # 5. mapDamage
    _log_progress("mapDamage")
    mapdamage_result = mapdamage_analyzer.run_mapdamage(sample_acc, dedup_bam)
    if not mapdamage_result:
        logger.error("mapDamage に失敗しました: %s", sample_acc)
        return sample_acc, False, "mapDamage"

    # 6. Qualimap
    _log_progress("Qualimap")
    qualimap_result = qualimap_analyzer.run_qualimap(sample_acc, dedup_bam)
    if not qualimap_result:
        logger.error("Qualimap に失敗しました: %s", sample_acc)
        return sample_acc, False, "Qualimap"

    # 7. HaplotypeCaller
    _log_progress("HaplotypeCaller")
    vcf_file = haplotypecaller.run_haplotypecaller(sample_acc, dedup_bam)
    if not vcf_file:
        logger.error("HaplotypeCaller に失敗しました: %s", sample_acc)
        return sample_acc, False, "HaplotypeCaller"

    # 8. 最終クリーンアップ
    dedup_bam_index = Path(str(dedup_bam) + ".bai")
    cleanup_intermediate_file(dedup_bam_index, logger)
    cleanup_intermediate_file(dedup_bam, logger)

    # 9. チェックポイント (.done) を記録
    done_flag.parent.mkdir(parents=True, exist_ok=True)
    done_flag.touch()

    logger.info("サンプルの解析を完了しました: %s", sample_acc)
    return sample_acc, True, ""


# ============================================================
# メインパイプライン
# ============================================================

def main() -> None:
    """NGS 解析パイプライン全体のエントリポイント。"""

    # 1. 設定とロギングを初期化
    args = parse_args()
    config = PipelineConfig(args)
    log_file = config.logs_dir / f"pipeline_{config.project_accession}.log"
    setup_logging(log_file=log_file)

    logger.info(
        "解析パイプラインを開始します: project_accession=%s",
        config.project_accession,
    )

    # 2. （オプション）HTTPS でダウンロードしてから解析
    if getattr(args, "download_via_https", False):
        logger.info("ena_download_https でダウンロードを実行します")
        out_dir = config.raw_data_dir.parent
        subprocess.run(
            [
                sys.executable,
                "-m",
                "modules.ena_download_https",
                "--project",
                config.project_accession,
                "--out",
                str(out_dir),
                "--workers",
                str(config.args.workers),
            ],
            check=True,
        )
        config.fastq_dir = config.raw_data_dir
        logger.info("ダウンロード完了。解析を開始します")

    # 3. リファレンスゲノムの存在を確認
    if not config.reference_genome.exists():
        logger.error("リファレンスゲノムが見つかりません: %s", config.reference_genome)
        sys.exit(1)

    # 4. リファレンスゲノムのインデックスを作成
    fai = config.reference_genome.with_suffix(".fai")
    if not fai.exists():
        logger.info("リファレンスゲノムのインデックスを作成します")
        subprocess.run(["samtools", "faidx", str(config.reference_genome)], check=True)

    # 5. 外部ツール・ファイルの事前検証
    config.validate_environment()

    session = requests.Session()

    # 6. データソースを決定
    if config.fastq_dir:
        if not config.fastq_dir.exists() or not config.fastq_dir.is_dir():
            logger.error(
                "指定された FASTQ ディレクトリが存在しないか、ディレクトリではありません: %s",
                config.fastq_dir,
            )
            sys.exit(1)

        logger.info("ローカルの FASTQ ファイルを使用します: %s", config.fastq_dir)
        sample_to_runs = group_fastqs_by_run(config.fastq_dir)

        if not sample_to_runs:
            logger.error("FASTQ ファイルが見つかりませんでした: %s", config.fastq_dir)
            sys.exit(1)

        for sa, runs in sample_to_runs.items():
            logger.info("  %s: %d ラン検出", sa, len(runs))

    else:
        logger.info("ENA からデータをダウンロードします")
        ena_downloader = ENADownloader(config)
        response_data = ena_downloader.get_api_response(
            config.project_accession, session
        )
        sample_to_files = ena_downloader.parse_response_with_checksums(response_data)

        downloaded_fastqs_root = config.results_dir / "ena_fastq"
        downloaded_fastqs_root.mkdir(parents=True, exist_ok=True)

        for sample_acc, url_md5_pairs in sample_to_files.items():
            logger.info("ENA sample %s をダウンロード中...", sample_acc)

            sample_dir = downloaded_fastqs_root / sample_acc
            sample_dir.mkdir(parents=True, exist_ok=True)

            ftp_urls = [u for u, _ in url_md5_pairs]
            fastq_files = ena_downloader.download_sample_data(sample_acc, ftp_urls)

            if not fastq_files:
                logger.warning("FASTQ がダウンロードできませんでした: %s", sample_acc)
                continue

            url_to_md5 = {
                Path(u.rstrip("/").split("/")[-1]): md5
                for u, md5 in url_md5_pairs
                if md5
            }
            for f in fastq_files:
                dest = sample_dir / f.name
                f.rename(dest)
                expected_md5 = url_to_md5.get(Path(f.name))
                if expected_md5 and not verify_file_md5(dest, expected_md5):
                    logger.warning("MD5 不一致のため削除します: %s", dest)
                    dest.unlink(missing_ok=True)

        sample_to_runs = group_fastqs_by_run(downloaded_fastqs_root)
        if not sample_to_runs:
            logger.error("FASTQ ファイルが見つかりませんでした: %s", downloaded_fastqs_root)
            sys.exit(1)

    # 7. 各サンプルの解析を実行
    parallel_samples: int = getattr(args, "parallel_samples", 1)
    total_samples = len(sample_to_runs)

    logger.info(
        "解析対象: %d サンプル (並列数: %d)",
        total_samples,
        parallel_samples,
    )

    succeeded: List[str] = []
    failed: List[Tuple[str, str]] = []

    if parallel_samples <= 1:
        for sample_acc, runs in tqdm(
            sample_to_runs.items(), desc="progress", unit="sample"
        ):
            acc, ok, step = process_sample(sample_acc, runs, config)
            if ok:
                succeeded.append(acc)
            else:
                failed.append((acc, step))
    else:
        with ThreadPoolExecutor(max_workers=parallel_samples) as executor:
            future_to_acc = {
                executor.submit(process_sample, acc, runs, config): acc
                for acc, runs in sample_to_runs.items()
            }
            with tqdm(total=total_samples, desc="progress", unit="sample") as pbar:
                for future in as_completed(future_to_acc):
                    acc = future_to_acc[future]
                    try:
                        _, ok, step = future.result()
                    except Exception:
                        logger.exception("サンプル %s で予期しない例外が発生しました", acc)
                        ok, step = False, "unexpected exception"
                    if ok:
                        succeeded.append(acc)
                    else:
                        failed.append((acc, step))
                    pbar.update(1)

    # 8. サマリーを出力
    logger.info("=" * 60)
    logger.info(
        "パイプライン完了: 成功 %d / %d サンプル",
        len(succeeded),
        total_samples,
    )
    if failed:
        logger.warning("失敗サンプル一覧:")
        for acc, step in failed:
            logger.warning("  %s  (失敗ステップ: %s)", acc, step)
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
