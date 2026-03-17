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
from modules.fastq_parser import parse_fastq_general, merge_lanes_by_cat
from modules.ena_downloader import ENADownloader
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
    fastq_files: List[Path],
    config: PipelineConfig,
) -> Tuple[str, bool, str]:
    """
    1 サンプル分の解析パイプライン（マッピング → QC → VCF）を実行する。

    各ワーカーでモジュールを独立にインスタンス化するため、
    スレッド間で状態を共有しない。

    完了したサンプルは <results_dir>/<sample_acc>/.done にフラグを残し、
    再実行時にスキップされる。--force で上書き可能。

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

    logger.info("サンプルの解析を開始します: %s", sample_acc)

    # 1. BWA マッピング
    bam_file = bwa_mapper.run_mapping_pipeline(sample_acc, fastq_files)
    if not bam_file:
        logger.error("マッピングに失敗しました: %s", sample_acc)
        return sample_acc, False, "BWA mapping"

    # 2. Soft clipping
    softclipped_bam = softclipper.run_softclipping(sample_acc, bam_file)
    if not softclipped_bam:
        logger.error("Soft clippingに失敗しました: %s", sample_acc)
        return sample_acc, False, "Soft clipping"

    cleanup_intermediate_file(bam_file, logger)

    # 3. BAM processing (dedup, index)
    dedup_bam = bam_processor.run_bam_processing(sample_acc, softclipped_bam)
    if not dedup_bam:
        logger.error("BAM processingに失敗しました: %s", sample_acc)
        return sample_acc, False, "BAM processing"

    # 4. mapDamage analysis (dedup 済み BAM を使用)
    mapdamage_result = mapdamage_analyzer.run_mapdamage(sample_acc, dedup_bam)
    if not mapdamage_result:
        logger.error("MapDamage analysisに失敗しました: %s", sample_acc)
        return sample_acc, False, "mapDamage"

    # 5. BWA の中間ファイルを削除
    bam_dir = dedup_bam.parent.parent / "bam_files"
    if bam_dir.exists():
        for pattern in ("*.bam", "*.bai", "*.truncated"):
            for intermediate in bam_dir.glob(pattern):
                cleanup_intermediate_file(intermediate, logger)

    # 6. Qualimap analysis
    qualimap_result = qualimap_analyzer.run_qualimap(sample_acc, dedup_bam)
    if not qualimap_result:
        logger.error("Qualimap analysisに失敗しました: %s", sample_acc)
        return sample_acc, False, "Qualimap"

    # 7. HaplotypeCaller
    vcf_file = haplotypecaller.run_haplotypecaller(sample_acc, dedup_bam)
    if not vcf_file:
        logger.error("HaplotypeCallerに失敗しました: %s", sample_acc)
        return sample_acc, False, "HaplotypeCaller"

    # 8. 中間ファイルを削除
    dedup_bam_index = Path(str(dedup_bam) + ".bai")
    cleanup_intermediate_file(dedup_bam_index, logger)
    cleanup_intermediate_file(dedup_bam, logger)
    cleanup_intermediate_file(softclipped_bam, logger)

    # 9. チェックポイント (.done) を記録
    done_flag.parent.mkdir(parents=True, exist_ok=True)
    done_flag.touch()

    logger.info("サンプルの解析を完了しました: %s", sample_acc)
    return sample_acc, True, ""


# ============================================================
# メインパイプライン
# ============================================================

def main() -> None:
    """
    NGS解析パイプライン全体のエントリポイント。
    ユーザーの入力に応じて、この関数はENAからリードデータをダウンロードするか、
    またはローカルディレクトリ内の事前ダウンロード済みFASTQファイルを処理します。
    """
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
        out_dir = config.raw_data_dir.parent  # base_dir/raw_data
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

    # 6. データソースを決定: ENAからのダウンロードか、ローカルのFASTQディレクトリ
    if config.fastq_dir:
        if not config.fastq_dir.exists() or not config.fastq_dir.is_dir():
            logger.error(
                "指定されたFASTQディレクトリが存在しないか、ディレクトリではありません: %s",
                config.fastq_dir,
            )
            sys.exit(1)

        logger.info("ローカルのFASTQファイルを使用します: %s", config.fastq_dir)

        sample_to_reads = parse_fastq_general(config.fastq_dir)
        if not sample_to_reads:
            logger.error("FASTQファイルが見つかりませんでした: %s", config.fastq_dir)
            sys.exit(1)

        merged_dir = config.results_dir / "merged_fastq"
        sample_to_fastqs = merge_lanes_by_cat(sample_to_reads, merged_dir, logger)

    else:
        logger.info("ENAからデータをダウンロードします")
        ena_downloader = ENADownloader(config)
        response_data = ena_downloader.get_api_response(
            config.project_accession, session
        )
        sample_to_ftp_urls = ena_downloader.parse_response_data(response_data)

        downloaded_fastqs_root = config.results_dir / "ena_fastq"
        downloaded_fastqs_root.mkdir(parents=True, exist_ok=True)

        for sample_acc, ftp_urls in sample_to_ftp_urls.items():
            logger.info(f"ENA sample {sample_acc} をダウンロード中...")
            fastq_files = ena_downloader.download_sample_data(sample_acc, ftp_urls)

            if not fastq_files:
                logger.warning(f"FASTQファイルがダウンロードできませんでした: {sample_acc}")
                continue

            sample_dir = downloaded_fastqs_root / sample_acc
            sample_dir.mkdir(parents=True, exist_ok=True)

            for f in fastq_files:
                f.rename(sample_dir / f.name)

        sample_to_reads = parse_fastq_general(downloaded_fastqs_root)
        if not sample_to_reads:
            logger.error("FASTQファイルが見つかりませんでした: %s", downloaded_fastqs_root)
            sys.exit(1)

        merged_dir = config.results_dir / "merged_fastq"
        sample_to_fastqs = merge_lanes_by_cat(sample_to_reads, merged_dir, logger)
        if not sample_to_fastqs:
            logger.error("レーンマージに失敗しました: %s", merged_dir)
            sys.exit(1)

    # 7. 各サンプルの解析を実行
    parallel_samples: int = getattr(args, "parallel_samples", 1)
    total_samples = len(sample_to_fastqs)

    logger.info(
        "解析対象: %d サンプル (並列数: %d)",
        total_samples,
        parallel_samples,
    )

    succeeded: List[str] = []
    failed: List[Tuple[str, str]] = []

    if parallel_samples <= 1:
        # --- 逐次実行 ---
        for sample_acc, fastq_files in tqdm(
            sample_to_fastqs.items(), desc="progress", unit="sample"
        ):
            acc, ok, step = process_sample(sample_acc, fastq_files, config)
            if ok:
                succeeded.append(acc)
            else:
                failed.append((acc, step))
    else:
        # --- 並列実行 ---
        with ThreadPoolExecutor(max_workers=parallel_samples) as executor:
            future_to_acc = {
                executor.submit(process_sample, acc, fqs, config): acc
                for acc, fqs in sample_to_fastqs.items()
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