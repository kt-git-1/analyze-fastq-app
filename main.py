import sys
import subprocess
import logging
import shutil
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

from config import (
    PipelineConfig,
    parse_args,
    setup_logging,
    cleanup_intermediate_file,
)
from modules.fastq_parser import FastqRun, group_fastqs_by_run, parse_fastq_general, merge_lanes_by_cat
from modules.bam_parser import group_bams_by_sample
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


FAILURE_HINTS = {
    "BWA mapping": "FASTQファイル、参照ゲノムのBWA index、bwa/samtoolsのインストールを確認してください。",
    "Soft clipping": "入力BAM、pysamで読み込めるBAM形式、CIGAR情報を確認してください。",
    "CleanSam": "Picard jar、Javaメモリ設定、入力BAMの整合性を確認してください。",
    "merge/dedup": "ランごとのBAM、samtools、Picard MarkDuplicates、Javaメモリ設定を確認してください。",
    "mapDamage": "mapDamage、参照ゲノム、BAM index、フィルタ後BAMを確認してください。",
    "Qualimap": "Qualimap、Java設定、BAM index、マッピング済みリード数を確認してください。",
    "HaplotypeCaller": "GATK、reference.fa/.fai/.dict、BAM index、スレッド数を確認してください。",
    "samtools index": "samtools、入力BAMの破損、BAMを書き込んだディレクトリの権限を確認してください。",
    "BAM missing": "指定したBAMディレクトリ、ファイル名パターン、サンプル名の抽出設定を確認してください。",
    "BAM empty": "入力BAMの作成元、前段のマッピング結果、ファイル転送の途中失敗を確認してください。",
    "all runs failed": "各ランのFASTQ、AdapterRemoval、BWA mapping、Soft clipping、CleanSamのログを確認してください。",
    "unexpected exception": "直前の例外ログと対象サンプルの入力ファイルを確認してください。",
}


def _progress_enabled(config: PipelineConfig) -> bool:
    return not getattr(config.args, "no_progress", False)


def _failure_hint(step_name: str) -> str:
    return FAILURE_HINTS.get(step_name, "直前の詳細ログと入力ファイルを確認してください。")


def _log_failure_hint(sample_acc: str, step_name: str) -> None:
    logger.error(
        "確認ポイント: サンプル %s / 失敗ステップ %s。%s",
        sample_acc,
        step_name,
        _failure_hint(step_name),
    )


def _is_nonempty_file(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _has_files(path: Path) -> bool:
    return path.exists() and path.is_dir() and any(p.is_file() for p in path.rglob("*"))


def _shorten_middle(text: str, width: int) -> str:
    if width <= 0:
        return ""
    if len(text) <= width:
        return text
    if width <= 3:
        return "." * width
    left = (width - 3) // 2
    right = width - 3 - left
    return text[:left] + "..." + text[-right:]


def _progress_bar(done: int, total: int, width: int = 30) -> str:
    if total <= 0:
        return "[" + "-" * width + "]"
    filled = int(width * min(done, total) / total)
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def _format_duration(seconds: float) -> str:
    seconds = max(0, int(seconds))
    minutes, sec = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return "約 %d時間%d分" % (hours, minutes)
    if minutes:
        return "約 %d分%d秒" % (minutes, sec)
    return "約 %d秒" % sec


def _set_stream_log_level(level: int) -> None:
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
            handler.setLevel(level)


class AnalysisDashboard:
    """解析パイプライン用の固定表示ダッシュボード。"""

    def __init__(
        self,
        project: str,
        input_mode: str,
        total_samples: int,
        parallel_samples: int,
        threads_per_sample: int,
        enabled: bool = True,
    ) -> None:
        self.project = project
        self.input_mode = input_mode
        self.total_samples = total_samples
        self.parallel_samples = parallel_samples
        self.threads_per_sample = threads_per_sample
        self.enabled = enabled
        self.completed_samples = 0
        self.succeeded_samples = 0
        self.failed_samples = 0
        self.sample_acc = "-"
        self.current_step = "-"
        self.current_step_index = 0
        self.total_steps = 0
        self.input_summary = "-"
        self.started_at = time.monotonic()
        self.events: List[str] = []
        self._rendered_lines = 0
        self._last_rendered_at = 0.0
        self._lock = threading.Lock()
        self._closed = threading.Event()
        self._heartbeat_thread: Optional[threading.Thread] = None
        if self.enabled:
            self._heartbeat_thread = threading.Thread(
                target=self._heartbeat,
                name="analysis-dashboard-heartbeat",
                daemon=True,
            )
            self._heartbeat_thread.start()

    def _heartbeat(self) -> None:
        while not self._closed.wait(1.0):
            self.render()

    def add_event(self, message: str) -> None:
        with self._lock:
            self._add_event_unlocked(message)
            self._render_unlocked(force=True)

    def _add_event_unlocked(self, message: str) -> None:
        self.events.append("%s  %s" % (time.strftime("%H:%M:%S"), message))
        self.events = self.events[-3:]

    def start_sample(self, sample_acc: str, total_steps: int, input_summary: str) -> None:
        with self._lock:
            self.sample_acc = sample_acc
            self.current_step = "-"
            self.current_step_index = 0
            self.total_steps = total_steps
            self.input_summary = input_summary
            self._add_event_unlocked("解析開始: %s" % sample_acc)
            self._render_unlocked(force=True)

    def start_step(self, sample_acc: str, step_name: str, step_index: int, total_steps: int) -> None:
        with self._lock:
            if self.sample_acc == sample_acc and self.current_step not in ("-", step_name):
                self._add_event_unlocked("完了: %s" % self.current_step)
            self.sample_acc = sample_acc
            self.current_step = step_name
            self.current_step_index = step_index
            self.total_steps = total_steps
            self._add_event_unlocked("開始: %s" % step_name)
            self._render_unlocked(force=True)

    def skip_steps(self, sample_acc: str, step_index: int, total_steps: int) -> None:
        with self._lock:
            self.sample_acc = sample_acc
            self.current_step_index = step_index
            self.total_steps = total_steps
            self._render_unlocked(force=True)

    def finish_sample(self, sample_acc: str, ok: bool, failed_step: str = "") -> None:
        with self._lock:
            if ok and self.sample_acc == sample_acc and self.current_step != "-":
                self._add_event_unlocked("完了: %s" % self.current_step)
            self.completed_samples += 1
            if ok:
                self.succeeded_samples += 1
                self._add_event_unlocked("解析完了: %s" % sample_acc)
            else:
                self.failed_samples += 1
                self._add_event_unlocked("解析失敗: %s / %s" % (sample_acc, failed_step))
            self._render_unlocked(force=True)

    def render_text(self) -> str:
        with self._lock:
            return self._render_text_unlocked()

    def _render_text_unlocked(self) -> str:
        overall_percent = int(round(self.completed_samples / self.total_samples * 100)) if self.total_samples else 100
        step_percent = int(round(self.current_step_index / self.total_steps * 100)) if self.total_steps else 0
        terminal_width = shutil.get_terminal_size((100, 20)).columns
        sample_width = max(18, min(36, terminal_width - 32))
        step_width = max(24, min(54, terminal_width - 12))
        event_lines = self.events or ["-"]

        return "\n".join(
            [
                "解析パイプライン: %s" % self.project,
                "",
                "全体進捗  %s  %3d%%  %d/%d samples" % (
                    _progress_bar(self.completed_samples, self.total_samples),
                    overall_percent,
                    self.completed_samples,
                    self.total_samples,
                ),
                "サンプル  %s" % _shorten_middle(self.sample_acc, sample_width),
                "現在      %s" % _shorten_middle(self.current_step, step_width),
                "ステップ  %s  %3d%%  %d/%d" % (
                    _progress_bar(self.current_step_index, self.total_steps),
                    step_percent,
                    self.current_step_index,
                    self.total_steps,
                ),
                "入力      %s" % self.input_summary,
                "並列数    %d samples" % self.parallel_samples,
                "スレッド  %d / sample" % self.threads_per_sample,
                "経過時間  %s" % _format_duration(time.monotonic() - self.started_at),
                "成功/失敗 %d / %d" % (self.succeeded_samples, self.failed_samples),
                "",
                "最近のイベント",
            ]
            + ["  %s" % line for line in event_lines]
        )

    def render(self, force: bool = False) -> None:
        with self._lock:
            self._render_unlocked(force=force)

    def _render_unlocked(self, force: bool = False) -> None:
        if not self.enabled:
            return
        now = time.monotonic()
        if not force and now - self._last_rendered_at < 0.2:
            return
        text = self._render_text_unlocked()
        if self._rendered_lines:
            sys.stderr.write("\033[%dF\033[J" % self._rendered_lines)
        sys.stderr.write(text + "\n")
        sys.stderr.flush()
        self._rendered_lines = text.count("\n") + 1
        self._last_rendered_at = now

    def close(self) -> None:
        self._closed.set()
        if self._heartbeat_thread is not None:
            self._heartbeat_thread.join(timeout=2)
        self.render(force=True)
        if self.enabled:
            sys.stderr.write("\n")
            sys.stderr.flush()


def _ensure_bam_index(bam_path: Path) -> bool:
    """Create a BAM index next to the BAM if it does not already exist."""
    bai_path = Path(str(bam_path) + ".bai")
    alternate_bai_path = bam_path.with_suffix(".bai")
    if bai_path.exists() or alternate_bai_path.exists():
        return True

    logger.info("BAM index が見つからないため作成します: %s", bai_path)
    try:
        subprocess.run(["samtools", "index", str(bam_path)], check=True)
        return True
    except subprocess.CalledProcessError:
        logger.exception("BAM index の作成に失敗しました: %s", bam_path)
        logger.error("確認ポイント: %s", _failure_hint("samtools index"))
        return False


def _existing_run_clean_bam(config: PipelineConfig, sample_acc: str, run_id: str) -> Optional[Path]:
    clean_bam = (
        config.results_dir
        / sample_acc
        / "runs"
        / run_id
        / "bam_files"
        / f"{run_id}.clean.bam"
    )
    if _is_nonempty_file(clean_bam):
        logger.info("再開: 既存CleanSam BAMを使用します: %s", clean_bam)
        return clean_bam
    return None


def _existing_sample_dedup_bam(config: PipelineConfig, sample_acc: str) -> Optional[Path]:
    dedup_bam = config.results_dir / sample_acc / "dedup" / f"{sample_acc}.dedup.sorted.bam"
    if _is_nonempty_file(dedup_bam):
        logger.info("再開: 既存dedup BAMを使用します: %s", dedup_bam)
        return dedup_bam
    return None


# ============================================================
# サンプル単位の解析関数
# ============================================================

def process_sample(
    sample_acc: str,
    runs: List[FastqRun],
    config: PipelineConfig,
    progress_position: Optional[int] = None,
    dashboard: Optional[AnalysisDashboard] = None,
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
    if dashboard is not None:
        dashboard.start_sample(sample_acc, total_steps, "FASTQ  %d runs" % len(runs))

    def _log_progress(step_name: str) -> None:
        nonlocal current_step
        current_step += 1
        pct = int(current_step / total_steps * 100)
        if dashboard is not None:
            dashboard.start_step(sample_acc, step_name, current_step, total_steps)
        logger.info(
            "[%s] %d%% (%d/%d) %s",
            sample_acc, pct, current_step, total_steps, step_name,
        )

    def _skip_progress(steps: int) -> None:
        nonlocal current_step
        current_step += steps
        if dashboard is not None:
            dashboard.skip_steps(sample_acc, current_step, total_steps)

    logger.info("サンプルの解析を開始します: %s (%d ラン, %d ステップ)", sample_acc, len(runs), total_steps)

    # ------------------------------------------------------------------
    # Phase 1: ラン単位の処理
    # ------------------------------------------------------------------
    run_clean_bams: List[Path] = []
    failed_runs: List[str] = []
    dedup_bam = _existing_sample_dedup_bam(config, sample_acc) if not force else None

    if dedup_bam is not None:
        if not _ensure_bam_index(dedup_bam):
            return sample_acc, False, "samtools index"
        _skip_progress(3 * len(runs) + 1)
    else:
        for run_idx, run in enumerate(runs, 1):
            run_id = run.run_id
            fastq_files = run.fastq_files
            clean_bam = _existing_run_clean_bam(config, sample_acc, run_id) if not force else None
            if clean_bam is not None:
                run_clean_bams.append(clean_bam)
                _skip_progress(3)
                continue

            # 1. BWA マッピング
            _log_progress(f"BWA mapping ({run_id}) [{run_idx}/{len(runs)}]")
            bam_file = bwa_mapper.run_mapping_pipeline(
                sample_acc,
                run_id,
                fastq_files,
                rg_library=run.rg_library,
            )
            if not bam_file:
                logger.error("マッピング失敗 → ラン %s をスキップ", run_id)
                _log_failure_hint(sample_acc, "BWA mapping")
                failed_runs.append(run_id)
                _skip_progress(2)
                continue

            # 2. Soft clipping
            _log_progress(f"Soft clipping ({run_id})")
            softclipped_bam = softclipper.run_softclipping(
                sample_acc,
                run_id,
                bam_file,
                progress_enabled=dashboard is None,
            )
            if not softclipped_bam:
                logger.error("Soft clipping 失敗 → ラン %s をスキップ", run_id)
                _log_failure_hint(sample_acc, "Soft clipping")
                failed_runs.append(run_id)
                _skip_progress(1)
                continue

            cleanup_intermediate_file(bam_file, logger)

            # 3. CleanSam (ラン単位)
            _log_progress(f"CleanSam ({run_id})")
            try:
                clean_bam = bam_processor.process_run_bam(sample_acc, run_id, softclipped_bam)
            except Exception:
                logger.exception("CleanSam 失敗 → ラン %s をスキップ", run_id)
                _log_failure_hint(sample_acc, "CleanSam")
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
            _log_failure_hint(sample_acc, "all runs failed")
            return sample_acc, False, "all runs failed"

    # ------------------------------------------------------------------
    # Phase 2: サンプル単位の処理
    # ------------------------------------------------------------------

    # 4. マージ + 重複除去
    if dedup_bam is None:
        _log_progress("Merge + MarkDuplicates")
        try:
            dedup_bam = bam_processor.merge_and_dedup(sample_acc, run_clean_bams)
        except Exception:
            logger.exception("merge/dedup に失敗しました: %s", sample_acc)
            _log_failure_hint(sample_acc, "merge/dedup")
            return sample_acc, False, "merge/dedup"

        if not dedup_bam:
            _log_failure_hint(sample_acc, "merge/dedup")
            return sample_acc, False, "merge/dedup"

        for bam in run_clean_bams:
            cleanup_intermediate_file(bam, logger)

    # 5. mapDamage
    mapdamage_dir = config.results_dir / sample_acc / "mapdamage"
    if not force and _has_files(mapdamage_dir):
        logger.info("再開: 既存mapDamage出力を使用します: %s", mapdamage_dir)
        _skip_progress(1)
    else:
        _log_progress("mapDamage")
        mapdamage_result = mapdamage_analyzer.run_mapdamage(sample_acc, dedup_bam)
        if not mapdamage_result:
            logger.error("mapDamage に失敗しました: %s", sample_acc)
            _log_failure_hint(sample_acc, "mapDamage")
            return sample_acc, False, "mapDamage"

    # 6. Qualimap
    qualimap_dir = config.results_dir / sample_acc / "qualimap"
    if not force and _has_files(qualimap_dir):
        logger.info("再開: 既存Qualimap出力を使用します: %s", qualimap_dir)
        _skip_progress(1)
    else:
        _log_progress("Qualimap")
        qualimap_result = qualimap_analyzer.run_qualimap(sample_acc, dedup_bam)
        if not qualimap_result:
            logger.error("Qualimap に失敗しました: %s", sample_acc)
            _log_failure_hint(sample_acc, "Qualimap")
            return sample_acc, False, "Qualimap"

    # 7. HaplotypeCaller
    _log_progress("HaplotypeCaller")
    vcf_file = haplotypecaller.run_haplotypecaller(sample_acc, dedup_bam)
    if not vcf_file:
        logger.error("HaplotypeCaller に失敗しました: %s", sample_acc)
        _log_failure_hint(sample_acc, "HaplotypeCaller")
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


def process_sample_from_dedup_bam(
    sample_acc: str,
    dedup_bam: Path,
    config: PipelineConfig,
    progress_position: Optional[int] = None,
    dashboard: Optional[AnalysisDashboard] = None,
) -> Tuple[str, bool, str]:
    """
    既存 dedup BAM から QC と VCF 出力のみを実行する。

    FASTQ 由来の前処理は行わず、入力 BAM は削除しない。
    """
    done_flag = config.results_dir / sample_acc / ".done"
    force: bool = getattr(config.args, "force", False)

    if done_flag.exists() and not force:
        logger.info("スキップ (完了済み): %s", sample_acc)
        return sample_acc, True, ""

    if not dedup_bam.exists() or not dedup_bam.is_file():
        logger.error("BAM ファイルが見つかりません: %s", dedup_bam)
        _log_failure_hint(sample_acc, "BAM missing")
        return sample_acc, False, "BAM missing"

    if dedup_bam.stat().st_size == 0:
        logger.error("BAM ファイルが空です: %s", dedup_bam)
        _log_failure_hint(sample_acc, "BAM empty")
        return sample_acc, False, "BAM empty"

    if not _ensure_bam_index(dedup_bam):
        return sample_acc, False, "samtools index"

    mapdamage_analyzer = MapDamageAnalyzer(config)
    qualimap_analyzer = QualimapAnalyzer(config)
    haplotypecaller = HaplotypeCaller(config)

    total_steps = 3
    current_step = 0
    if dashboard is not None:
        dashboard.start_sample(sample_acc, total_steps, "BAM  既存dedup BAM")

    def _log_progress(step_name: str) -> None:
        nonlocal current_step
        current_step += 1
        pct = int(current_step / total_steps * 100)
        if dashboard is not None:
            dashboard.start_step(sample_acc, step_name, current_step, total_steps)
        logger.info(
            "[%s] %d%% (%d/%d) %s",
            sample_acc, pct, current_step, total_steps, step_name,
        )

    logger.info("既存 dedup BAM から解析を開始します: %s (%s)", sample_acc, dedup_bam)

    mapdamage_dir = config.results_dir / sample_acc / "mapdamage"
    if not force and _has_files(mapdamage_dir):
        logger.info("再開: 既存mapDamage出力を使用します: %s", mapdamage_dir)
        current_step += 1
        if dashboard is not None:
            dashboard.skip_steps(sample_acc, current_step, total_steps)
    else:
        _log_progress("mapDamage")
        mapdamage_result = mapdamage_analyzer.run_mapdamage(sample_acc, dedup_bam)
        if not mapdamage_result:
            logger.error("mapDamage に失敗しました: %s", sample_acc)
            _log_failure_hint(sample_acc, "mapDamage")
            return sample_acc, False, "mapDamage"

    qualimap_dir = config.results_dir / sample_acc / "qualimap"
    if not force and _has_files(qualimap_dir):
        logger.info("再開: 既存Qualimap出力を使用します: %s", qualimap_dir)
        current_step += 1
        if dashboard is not None:
            dashboard.skip_steps(sample_acc, current_step, total_steps)
    else:
        _log_progress("Qualimap")
        qualimap_result = qualimap_analyzer.run_qualimap(sample_acc, dedup_bam)
        if not qualimap_result:
            logger.error("Qualimap に失敗しました: %s", sample_acc)
            _log_failure_hint(sample_acc, "Qualimap")
            return sample_acc, False, "Qualimap"

    _log_progress("HaplotypeCaller")
    vcf_file = haplotypecaller.run_haplotypecaller(sample_acc, dedup_bam)
    if not vcf_file:
        logger.error("HaplotypeCaller に失敗しました: %s", sample_acc)
        _log_failure_hint(sample_acc, "HaplotypeCaller")
        return sample_acc, False, "HaplotypeCaller"

    done_flag.parent.mkdir(parents=True, exist_ok=True)
    done_flag.touch()

    logger.info("既存 dedup BAM からの解析を完了しました: %s", sample_acc)
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
    progress_enabled = _progress_enabled(config)
    setup_logging(log_file=log_file, use_tqdm=progress_enabled)

    logger.info(
        "解析パイプラインを開始します: project_accession=%s",
        config.project_accession,
    )

    # 2. （オプション）HTTPS でダウンロードしてから解析
    if getattr(args, "download_via_https", False):
        logger.info("ena_download_https でダウンロードを実行します")
        out_dir = config.raw_data_dir.parent
        download_cmd = [
            sys.executable,
            "-m",
            "modules.ena_download_https",
            "--project",
            config.project_accession,
            "--out",
            str(out_dir),
            "--workers",
            str(config.args.workers),
        ]
        if not progress_enabled:
            download_cmd.append("--no-progress")
        subprocess.run(download_cmd, check=True)
        config.fastq_dir = config.raw_data_dir
        logger.info("ダウンロード完了。解析を開始します")

    # 3. リファレンスゲノムの存在を確認
    if not config.reference_genome.exists():
        logger.error("リファレンスゲノムが見つかりません: %s", config.reference_genome)
        sys.exit(1)

    # 4. リファレンスゲノムのインデックスを作成
    fai = Path(str(config.reference_genome) + ".fai")
    if not fai.exists():
        logger.info("リファレンスゲノムのインデックスを作成します")
        subprocess.run(["samtools", "faidx", str(config.reference_genome)], check=True)

    # 5. 外部ツール・ファイルの事前検証
    config.validate_environment()

    session = requests.Session()

    # 6. データソースを決定
    bam_mode = config.bam_dir is not None
    sample_to_bams: Dict[str, Path] = {}

    if config.bam_dir:
        if not config.bam_dir.exists() or not config.bam_dir.is_dir():
            logger.error(
                "指定された BAM ディレクトリが存在しないか、ディレクトリではありません: %s",
                config.bam_dir,
            )
            sys.exit(1)

        logger.info("既存 dedup BAM ファイルを使用します: %s", config.bam_dir)
        try:
            sample_to_bams = group_bams_by_sample(config.bam_dir, config.args.bam_pattern)
        except ValueError as exc:
            logger.error("BAM ファイルの検出に失敗しました: %s", exc)
            sys.exit(1)

        if not sample_to_bams:
            logger.error("BAM ファイルが見つかりませんでした: %s", config.bam_dir)
            sys.exit(1)

        for sa, bam in sample_to_bams.items():
            logger.info("  %s: %s", sa, bam)

    elif config.fastq_dir:
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
    total_samples = len(sample_to_bams) if bam_mode else len(sample_to_runs)

    logger.info(
        "解析対象: %d サンプル (並列数: %d)",
        total_samples,
        parallel_samples,
    )

    succeeded: List[str] = []
    failed: List[Tuple[str, str]] = []
    dashboard = AnalysisDashboard(
        config.project_accession,
        "BAM" if bam_mode else "FASTQ",
        total_samples,
        parallel_samples,
        config.args.threads,
        enabled=progress_enabled,
    )
    if progress_enabled:
        _set_stream_log_level(logging.WARNING)
    dashboard.render(force=True)

    try:
        if parallel_samples <= 1:
            sample_items = list(sample_to_bams.items() if bam_mode else sample_to_runs.items())
            for sample_acc, sample_input in sample_items:
                if bam_mode:
                    acc, ok, step = process_sample_from_dedup_bam(
                        sample_acc, sample_input, config, dashboard=dashboard,
                    )
                else:
                    acc, ok, step = process_sample(
                        sample_acc, sample_input, config, dashboard=dashboard,
                    )
                dashboard.finish_sample(acc, ok, step)
                if ok:
                    succeeded.append(acc)
                else:
                    failed.append((acc, step))
        else:
            with ThreadPoolExecutor(max_workers=parallel_samples) as executor:
                sample_items = list(sample_to_bams.items() if bam_mode else sample_to_runs.items())
                if bam_mode:
                    future_to_acc = {
                        executor.submit(
                            process_sample_from_dedup_bam, acc, bam, config, idx + 1, dashboard,
                        ): acc
                        for idx, (acc, bam) in enumerate(sample_items)
                    }
                else:
                    future_to_acc = {
                        executor.submit(process_sample, acc, runs, config, idx + 1, dashboard): acc
                        for idx, (acc, runs) in enumerate(sample_items)
                    }
                for future in as_completed(future_to_acc):
                    acc = future_to_acc[future]
                    try:
                        _, ok, step = future.result()
                    except Exception:
                        logger.exception("サンプル %s で予期しない例外が発生しました", acc)
                        _log_failure_hint(acc, "unexpected exception")
                        ok, step = False, "unexpected exception"
                    if ok:
                        succeeded.append(acc)
                    else:
                        failed.append((acc, step))
                    dashboard.finish_sample(acc, ok, step)
    finally:
        dashboard.close()
        if progress_enabled:
            _set_stream_log_level(logging.INFO)

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
        logger.warning("失敗ステップ別:")
        for step, count in Counter(step for _, step in failed).most_common():
            logger.warning("  %s: %d サンプル", step, count)
            logger.warning("    確認: %s", _failure_hint(step))
    logger.info("出力先: %s", config.results_dir)
    logger.info("ログファイル: %s", log_file)
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
