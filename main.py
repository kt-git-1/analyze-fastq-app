import re
import sys
import shlex
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional

import requests
from tqdm import tqdm

from config import (
    PipelineConfig,
    parse_args,
    setup_logging,
    cleanup_intermediate_file,
    parse_fastq_general,
    merge_lanes_by_cat,
)
from modules.ena_downloader import ENADownloader
from modules.bwa_mapper import BWAMapper
from modules.softclipper import SoftClipper
from modules.bam_processor import BAMProcessor
from modules.analyzers import (
    MapDamageAnalyzer,
    QualimapAnalyzer,
    HaplotypeCaller,
)


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
    logger = setup_logging(log_file=log_file)

    logger.info(
        "解析パイプラインを開始します: project_accession=%s",
        config.project_accession,
    )

    # 2. リファレンスゲノムの存在を確認
    if not config.reference_genome.exists():
        logger.error("リファレンスゲノムが見つかりません: %s", config.reference_genome)
        sys.exit(1)

    # 3. リファレンスゲノムのインデックスを作成
    fai = config.reference_genome.with_suffix(".fai")
    if not fai.exists():
        logger.info("リファレンスゲノムのインデックスを作成します")
        subprocess.run(["samtools", "faidx", str(config.reference_genome)], check=True)

    # 4. モジュールをインスタンス化
    ena_downloader = ENADownloader(config)
    bwa_mapper = BWAMapper(config)
    softclipper = SoftClipper(config)
    bam_processor = BAMProcessor(config)
    mapdamage_analyzer = MapDamageAnalyzer(config)
    qualimap_analyzer = QualimapAnalyzer(config)
    haplotypecaller = HaplotypeCaller(config)

    session = requests.Session()

    # 5. データソースを決定: ENAからのダウンロードか、ローカルのFASTQディレクトリ
    if config.fastq_dir:
        # 6. ローカルのFASTQディレクトリを処理
        if not config.fastq_dir.exists() or not config.fastq_dir.is_dir():
            logger.error(
                "指定されたFASTQディレクトリが存在しないか、ディレクトリではありません: %s",
                config.fastq_dir,
            )
            sys.exit(1)

        logger.info("ローカルのFASTQファイルを使用します: %s", config.fastq_dir)

        # ★ 6-1. いろいろな命名パターンをまとめてパース
        sample_to_reads = parse_fastq_general(config.fastq_dir)

        if not sample_to_reads:
            logger.error("FASTQファイルが見つかりませんでした: %s", config.fastq_dir)
            sys.exit(1)

        # ★ 6-2. レーンをマージ（複数レーンがあれば cat、1ファイルならそのまま）
        merged_dir = config.results_dir / "merged_fastq"
        sample_to_fastqs = merge_lanes_by_cat(sample_to_reads, merged_dir, logger)

    else:
        # 7. ENAからデータをダウンロード
        logger.info("ENAからデータをダウンロードします")
        response_data = ena_downloader.get_api_response(
            config.project_accession, session
        )
        sample_to_ftp_urls = ena_downloader.parse_response_data(response_data)

        # ENAダウンロードしたFASTQを一時的に保存する辞書
        downloaded_fastqs_root = config.results_dir / "ena_fastq"
        downloaded_fastqs_root.mkdir(parents=True, exist_ok=True)

        # ENAダウンロードしたFASTQファイル全てを一時保存
        for sample_acc, ftp_urls in sample_to_ftp_urls.items():
            logger.info(f"ENA sample {sample_acc} をダウンロード中...")
            fastq_files = ena_downloader.download_sample_data(sample_acc, ftp_urls)

            if not fastq_files:
                logger.warning(f"FASTQファイルがダウンロードできませんでした: {sample_acc}")
                continue

            # サンプル毎の保存ディレクトリ
            sample_dir = downloaded_fastqs_root / sample_acc
            sample_dir.mkdir(parents=True, exist_ok=True)

            # ダウンロードしたFASTQをサンプルディレクトリへ移動
            for f in fastq_files:
                f.rename(sample_dir / f.name)

        # 8. ダウンロード済みFASTQをまとめてパース
        sample_to_reads = parse_fastq_general(downloaded_fastqs_root)
        if not sample_to_reads:
            logger.error("FASTQファイルが見つかりませんでした: %s", downloaded_fastqs_root)
            sys.exit(1)            

        # 9. レーンマージ
        merged_dir = config.results_dir / "merged_fastq"
        sample_to_fastqs = merge_lanes_by_cat(sample_to_reads, merged_dir, logger)
        if not sample_to_fastqs:
            logger.error("レーンマージに失敗しました: %s", merged_dir)
            sys.exit(1)

    # 10. 各サンプルの解析を実行
    for sample_acc, fastq_files in tqdm(
        sample_to_fastqs.items(), desc="progress", unit="sample"
    ):
        logger.info("サンプル/ラインの解析を実行します: %s", sample_acc)

        # 11. BWAマッピングを実行
        bam_file = bwa_mapper.run_mapping_pipeline(sample_acc, fastq_files)
        if not bam_file:
            logger.error("マッピングに失敗しました: %s", sample_acc)
            continue

        # 12. Soft clippingを実行
        softclipped_bam = softclipper.run_softclipping(sample_acc, bam_file)
        if not softclipped_bam:
            logger.error("Soft clippingに失敗しました: %s", sample_acc)
            continue

        cleanup_intermediate_file(bam_file, logger)

        # 13. BAM processing (dedup, index)
        dedup_bam = bam_processor.run_bam_processing(sample_acc, softclipped_bam)
        if not dedup_bam:
            logger.error("BAM processingに失敗しました: %s", sample_acc)
            continue

        # 14. mapDamage analysisを実行
        mapdamage_result = mapdamage_analyzer.run_mapdamage(sample_acc, softclipped_bam)
        if not mapdamage_result:
            logger.error("MapDamage analysisに失敗しました: %s", sample_acc)
            continue

        # 15. BWAの中間ファイルを削除
        bam_dir = dedup_bam.parent.parent / "bam_files"
        if bam_dir.exists():
            for pattern in ("*.bam", "*.bai", "*.truncated"):
                for intermediate in bam_dir.glob(pattern):
                    cleanup_intermediate_file(intermediate, logger)

        # 16. Qualimap analysisを実行
        qualimap_result = qualimap_analyzer.run_qualimap(sample_acc, dedup_bam)
        if not qualimap_result:
            logger.error("Qualimap analysisに失敗しました: %s", sample_acc)
            continue

        # 17. HaplotypeCallerを実行
        vcf_file = haplotypecaller.run_haplotypecaller(sample_acc, dedup_bam)
        if not vcf_file:
            logger.error("HaplotypeCallerに失敗しました: %s", sample_acc)
            continue

        # 18. 重複除去済みBAMとそのインデックスを削除
        dedup_bam_index = Path(str(dedup_bam) + ".bai")
        cleanup_intermediate_file(dedup_bam_index, logger)
        cleanup_intermediate_file(dedup_bam, logger)
        cleanup_intermediate_file(softclipped_bam, logger)

        logger.info("サンプル/ラインの解析を完了しました: %s", sample_acc)

    logger.info("パイプラインを完了しました!")


if __name__ == "__main__":
    main()