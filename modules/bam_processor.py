import subprocess
import logging
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


def _remove_intermediate(path: Path) -> None:
    """中間 BAM を削除する。存在しなければ何もしない。"""
    try:
        if path.exists():
            path.unlink()
            logger.info("中間ファイルを削除しました: %s", path)
    except OSError as exc:
        logger.warning("中間ファイルを削除できませんでした: %s (%s)", path, exc)


class BAMProcessor:
    def __init__(self, config):
        self.config = config

    # ------------------------------------------------------------------
    # ラン単位: CleanSam のみ（RG は BWA で設定済み）
    # ------------------------------------------------------------------

    def process_run_bam(
        self,
        sample_acc: str,
        run_id: str,
        softclipped_bam: Path,
    ) -> Optional[Path]:
        """
        1 ラン分の BAM に CleanSam を適用する。

        RG ヘッダーは BWA MEM 時点で設定済みなので
        AddOrReplaceReadGroups は不要。
        """
        logger.info("ラン BAM 処理を開始: %s / %s", sample_acc, run_id)
        run_bam_dir = self.config.results_dir / sample_acc / "runs" / run_id / "bam_files"
        run_bam_dir.mkdir(parents=True, exist_ok=True)

        clean_bam = run_bam_dir / f"{run_id}.clean.bam"
        picard_jar = str(self.config.args.picard_jar)

        self._run_cmd(
            [
                "java", "-Xmx" + self.config.args.java_mem, "-jar", picard_jar,
                "CleanSam",
                "I=" + str(softclipped_bam),
                "O=" + str(clean_bam),
                "VALIDATION_STRINGENCY=LENIENT",
            ],
            "Picard CleanSam",
        )

        logger.info("ラン BAM 処理完了: %s / %s", sample_acc, run_id)
        return clean_bam

    # ------------------------------------------------------------------
    # サンプル単位: merge → MarkDuplicates → sort → index
    # ------------------------------------------------------------------

    def merge_and_dedup(
        self,
        sample_acc: str,
        run_bams: List[Path],
    ) -> Optional[Path]:
        """
        複数ランの BAM をマージし、重複除去・sort・index を行う。

        MarkDuplicates はラン横断で重複を検出するため、
        ラン単位 RG が設定された状態でマージしてから実行する。
        """
        logger.info("サンプル BAM マージ+dedup を開始: %s (%d ラン)", sample_acc, len(run_bams))

        sample_dedup_dir = self.config.results_dir / sample_acc / "dedup"
        sample_dedup_dir.mkdir(parents=True, exist_ok=True)

        picard_jar = str(self.config.args.picard_jar)

        # --- merge ---
        merged_bam = sample_dedup_dir / f"{sample_acc}.merged.bam"
        valid_bams = [str(b) for b in run_bams if b is not None and b.exists()]
        if not valid_bams:
            logger.error("マージ対象の BAM がありません: %s", sample_acc)
            return None

        if len(valid_bams) == 1:
            import shutil
            shutil.move(valid_bams[0], str(merged_bam))
        else:
            self._run_cmd(
                ["samtools", "merge", "-f", str(merged_bam), *valid_bams],
                "samtools merge",
            )

        # --- MarkDuplicates ---
        marked_bam = sample_dedup_dir / f"{sample_acc}.marked.bam"
        metrics_file = sample_dedup_dir / f"{sample_acc}.marked_dup_metrics.txt"
        self._run_cmd(
            [
                "java", "-Xmx" + self.config.args.java_mem, "-jar", picard_jar,
                "MarkDuplicates",
                "I=" + str(merged_bam),
                "O=" + str(marked_bam),
                "M=" + str(metrics_file),
                "REMOVE_DUPLICATES=true",
                "VALIDATION_STRINGENCY=LENIENT",
            ],
            "Picard MarkDuplicates",
        )
        _remove_intermediate(merged_bam)

        # --- sort ---
        final_bam = sample_dedup_dir / f"{sample_acc}.dedup.sorted.bam"
        self._run_cmd(
            ["samtools", "sort", "-o", str(final_bam), str(marked_bam)],
            "samtools sort (dedup)",
        )
        _remove_intermediate(marked_bam)

        # --- index ---
        self._run_cmd(
            ["samtools", "index", str(final_bam)],
            "samtools index",
        )

        logger.info("サンプル BAM dedup 完了: %s", sample_acc)
        return final_bam

    # ------------------------------------------------------------------
    # ユーティリティ
    # ------------------------------------------------------------------

    def _run_cmd(self, cmd: list, step_name: str) -> None:
        try:
            logger.info("Running: %s", " ".join(str(c) for c in cmd))
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.error("Error in %s: %s", step_name, e)
            if e.stderr:
                logger.error("stderr: %s", e.stderr)
            raise
