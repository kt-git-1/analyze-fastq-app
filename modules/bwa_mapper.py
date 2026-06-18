import subprocess
import logging
import shutil
import threading
from pathlib import Path
from typing import List, Optional

from modules.logging_utils import format_command, log_command_start, log_tool_output

logger = logging.getLogger(__name__)

class BWAMapper:
    def __init__(self, config):
        self.config = config

    @staticmethod
    def _adapter_output_path(temp_prefix: Path, suffix: str) -> Path:
        return Path(str(temp_prefix) + suffix)

    def run_mapping_pipeline(
        self,
        sample_acc: str,
        run_id: str,
        fastq_files: List[Path],
        rg_library: Optional[str] = None,
    ) -> Optional[Path]:
        """
        1 ラン分の AdapterRemoval + BWA MEM を実行し、sorted BAM を返す。

        RG ヘッダーは run_id を ID に、sample_acc を SM に設定する。
        出力先は results/<sample_acc>/runs/<run_id>/bam_files/ 。
        """
        logger.info("マッピングパイプラインを開始します: %s / %s", sample_acc, run_id)

        run_dir = self.config.results_dir / sample_acc / "runs" / run_id
        bam_dir = run_dir / "bam_files"
        bam_dir.mkdir(parents=True, exist_ok=True)

        effective_rg_library = rg_library or getattr(self.config.args, "rg_library", "unknown")
        rg_center = getattr(self.config.args, "rg_center", "unknown")
        rg_string = (
            f"@RG\\tID:{run_id}\\tSM:{sample_acc}"
            f"\\tLB:{effective_rg_library}\\tCN:{rg_center}\\tPL:ILLUMINA"
        )

        # --- paired-end ---
        if len(fastq_files) >= 2:
            fastq1, fastq2 = fastq_files[0], fastq_files[1]
            logger.info(
                "AdapterRemoval (ペアエンド, data_type=%s): %s / %s",
                self.config.data_type, sample_acc, run_id,
            )
            temp_prefix = self.config.temp_dir / sample_acc / run_id
            temp_prefix.parent.mkdir(parents=True, exist_ok=True)
            cmd = self._build_adapter_removal_cmd_paired(fastq1, fastq2, temp_prefix)
            try:
                self._run_streaming(cmd, f"AdapterRemoval {sample_acc}/{run_id}")
            except subprocess.CalledProcessError as e:
                logger.error("AdapterRemoval に失敗しました: %s / %s: %s", sample_acc, run_id, e)
                return None

            if self.config.data_type == "ancient":
                return self._map_ancient_pe(
                    temp_prefix, bam_dir, run_id, rg_string, sample_acc,
                )
            else:
                return self._map_modern_pe(
                    temp_prefix, bam_dir, run_id, rg_string, sample_acc,
                )

        # --- single-end ---
        elif len(fastq_files) == 1:
            fastq1 = fastq_files[0]
            logger.info(
                "AdapterRemoval (シングルエンド, data_type=%s): %s / %s",
                self.config.data_type, sample_acc, run_id,
            )
            temp_prefix = self.config.temp_dir / sample_acc / run_id
            temp_prefix.parent.mkdir(parents=True, exist_ok=True)
            cmd = self._build_adapter_removal_cmd_single(fastq1, temp_prefix)
            try:
                self._run_streaming(cmd, f"AdapterRemoval {sample_acc}/{run_id}")
            except subprocess.CalledProcessError as e:
                logger.error("AdapterRemoval に失敗しました: %s / %s: %s", sample_acc, run_id, e)
                return None

            temp_single = self._adapter_output_path(temp_prefix, ".truncated")
            if not temp_single.exists():
                logger.error("AdapterRemoval 出力が見つかりません: %s", temp_single)
                return None

            dest_single = bam_dir / f"{run_id}.truncated"
            shutil.move(str(temp_single), str(dest_single))

            bam_file = bam_dir / f"{run_id}.sorted.bam"
            bam = self._run_bwa_se_and_sort(dest_single, bam_file, rg_string, run_id)
            if bam is None:
                logger.error("BWA マッピングに失敗しました: %s / %s", sample_acc, run_id)
            return bam

        else:
            logger.error("FASTQ ファイルがありません: %s / %s", sample_acc, run_id)
            return None

    # ------------------------------------------------------------------
    # ancient paired-end: collapsed + non-collapsed → merge
    # ------------------------------------------------------------------

    def _map_ancient_pe(
        self,
        temp_prefix: Path,
        bam_dir: Path,
        run_id: str,
        rg_string: str,
        sample_acc: str,
    ) -> Optional[Path]:
        collapsed_file = self._adapter_output_path(temp_prefix, ".collapsed.truncated")
        r1_file = self._adapter_output_path(temp_prefix, ".pair1.truncated")
        r2_file = self._adapter_output_path(temp_prefix, ".pair2.truncated")
        partial_bams: List[Path] = []

        if collapsed_file.exists():
            collapsed_bam = bam_dir / f"{run_id}.collapsed.sorted.bam"
            bam = self._run_bwa_se_and_sort(collapsed_file, collapsed_bam, rg_string, run_id)
            if bam is not None:
                partial_bams.append(bam)
        else:
            logger.warning("collapsed ファイルが見つかりません: %s", collapsed_file)

        if r1_file.exists() and r2_file.exists():
            pe_bam = bam_dir / f"{run_id}.pe.sorted.bam"
            bam = self._run_bwa_pe_and_sort(r1_file, r2_file, pe_bam, rg_string, run_id)
            if bam is not None:
                partial_bams.append(bam)
        else:
            if not r1_file.exists():
                logger.warning("R1 ファイルが見つかりません: %s", r1_file)
            if not r2_file.exists():
                logger.warning("R2 ファイルが見つかりません: %s", r2_file)

        if not partial_bams:
            logger.error("BWA マッピングに失敗しました: %s / %s", sample_acc, run_id)
            return None

        merged_bam = bam_dir / f"{run_id}.merged.sorted.bam"
        if not self._merge_bams(partial_bams, merged_bam):
            logger.error("BAM マージに失敗しました: %s / %s", sample_acc, run_id)
            return None
        return merged_bam

    # ------------------------------------------------------------------
    # modern paired-end
    # ------------------------------------------------------------------

    def _map_modern_pe(
        self,
        temp_prefix: Path,
        bam_dir: Path,
        run_id: str,
        rg_string: str,
        sample_acc: str,
    ) -> Optional[Path]:
        temp_r1 = self._adapter_output_path(temp_prefix, ".pair1.truncated")
        temp_r2 = self._adapter_output_path(temp_prefix, ".pair2.truncated")
        if not (temp_r1.exists() and temp_r2.exists()):
            logger.error("AdapterRemoval 出力が見つかりません: %s / %s", sample_acc, run_id)
            return None

        dest_r1 = bam_dir / f"{run_id}.R1.truncated"
        dest_r2 = bam_dir / f"{run_id}.R2.truncated"
        shutil.move(str(temp_r1), str(dest_r1))
        shutil.move(str(temp_r2), str(dest_r2))

        bam_file = bam_dir / f"{run_id}.sorted.bam"
        bam = self._run_bwa_pe_and_sort(dest_r1, dest_r2, bam_file, rg_string, run_id)
        if bam is None:
            logger.error("BWA マッピングに失敗しました: %s / %s", sample_acc, run_id)
        return bam

    # ------------------------------------------------------------------
    # ユーティリティ: stderr ストリーミング実行
    # ------------------------------------------------------------------

    def _run_streaming(self, cmd: list, step_name: str) -> None:
        """コマンドを実行し stderr をリアルタイムでログに出力する。"""
        log_command_start(logger, cmd, step_name)
        proc = subprocess.Popen(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True,
        )
        for line in iter(proc.stderr.readline, ""):
            line = line.rstrip()
            if line:
                log_tool_output(logger, step_name, line)
        proc.stderr.close()
        proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)

    # ------------------------------------------------------------------
    # AdapterRemoval コマンド構築
    # ------------------------------------------------------------------

    def _build_adapter_removal_cmd_paired(self, fastq1, fastq2, temp_prefix):
        cmd = [
            "AdapterRemoval",
            "--file1", str(fastq1),
            "--file2", str(fastq2),
            "--trimns", "--trimqualities",
            "--basename", str(temp_prefix),
            "--threads", str(self.config.args.threads),
        ]
        if self.config.data_type == "ancient":
            cmd.extend(["--minquality", "20", "--minlength", "30", "--collapse"])
        elif self.config.data_type == "modern":
            cmd.extend(["--minquality", "25", "--minlength", "25"])
        return cmd

    def _build_adapter_removal_cmd_single(self, fastq1, temp_prefix):
        cmd = [
            "AdapterRemoval",
            "--file1", str(fastq1),
            "--trimns", "--trimqualities",
            "--basename", str(temp_prefix),
            "--threads", str(self.config.args.threads),
        ]
        if self.config.data_type == "ancient":
            cmd.extend(["--minquality", "20", "--minlength", "30"])
        elif self.config.data_type == "modern":
            cmd.extend(["--minquality", "25", "--minlength", "25"])
        return cmd

    # ------------------------------------------------------------------
    # BWA MEM + samtools sort
    # ------------------------------------------------------------------

    def _run_bwa_se_and_sort(
        self, fastq: Path, bam_out: Path, rg_string: str, label: str,
    ) -> Optional[Path]:
        logger.info("single-end BWA: %s → %s", fastq.name, bam_out.name)
        bwa_cmd = [
            "bwa", "mem",
            "-t", str(self.config.args.threads),
            "-K", "100000000", "-Y",
            "-R", rg_string,
            str(self.config.reference_genome),
            str(fastq),
        ]
        return self._pipe_bwa_sort(bwa_cmd, bam_out, label)

    def _run_bwa_pe_and_sort(
        self,
        fastq1: Path,
        fastq2: Path,
        bam_out: Path,
        rg_string: str,
        label: str,
    ) -> Optional[Path]:
        logger.info("pair-end BWA: %s, %s → %s", fastq1.name, fastq2.name, bam_out.name)
        bwa_cmd = [
            "bwa", "mem",
            "-t", str(self.config.args.threads),
            "-K", "100000000", "-Y",
            "-R", rg_string,
            str(self.config.reference_genome),
            str(fastq1),
            str(fastq2),
        ]
        return self._pipe_bwa_sort(bwa_cmd, bam_out, label)

    def _pipe_bwa_sort(
        self, bwa_cmd: list, bam_out: Path, label: str,
    ) -> Optional[Path]:
        """bwa mem | samtools view -b | samtools sort のパイプラインを実行する。"""
        view_cmd = ["samtools", "view", "-b", "-"]
        sort_cmd = ["samtools", "sort", "-o", str(bam_out), "-"]
        try:
            logger.info("外部ツールを実行します: BWA %s", label)
            logger.debug("実行コマンド (BWA %s): %s", label, format_command(bwa_cmd))
            logger.debug("実行コマンド (samtools view %s): %s", label, format_command(view_cmd))
            logger.debug("実行コマンド (samtools sort %s): %s", label, format_command(sort_cmd))
            bwa_proc = subprocess.Popen(
                bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )

            def _drain_bwa_stderr() -> None:
                for raw in bwa_proc.stderr:
                    line = raw.decode("utf-8", errors="replace").rstrip()
                    if line:
                        log_tool_output(logger, "BWA %s" % label, line)
                bwa_proc.stderr.close()

            stderr_thread = threading.Thread(target=_drain_bwa_stderr, daemon=True)
            stderr_thread.start()

            view_proc = subprocess.Popen(
                view_cmd,
                stdin=bwa_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            bwa_proc.stdout.close()
            sort_proc = subprocess.Popen(
                sort_cmd,
                stdin=view_proc.stdout,
                stderr=subprocess.PIPE,
                text=True,
            )
            view_proc.stdout.close()

            def _drain_text_stderr(proc, step_name: str) -> None:
                for line in iter(proc.stderr.readline, ""):
                    line = line.rstrip()
                    if line:
                        log_tool_output(logger, step_name, line)
                proc.stderr.close()

            view_stderr_thread = threading.Thread(
                target=_drain_text_stderr,
                args=(view_proc, "samtools view %s" % label),
                daemon=True,
            )
            sort_stderr_thread = threading.Thread(
                target=_drain_text_stderr,
                args=(sort_proc, "samtools sort %s" % label),
                daemon=True,
            )
            view_stderr_thread.start()
            sort_stderr_thread.start()

            sort_proc.wait()
            if sort_proc.returncode != 0:
                raise subprocess.CalledProcessError(sort_proc.returncode, sort_cmd)
            view_proc.wait()
            if view_proc.returncode != 0:
                raise subprocess.CalledProcessError(view_proc.returncode, view_cmd)
            bwa_proc.wait()
            if bwa_proc.returncode != 0:
                raise subprocess.CalledProcessError(bwa_proc.returncode, bwa_cmd)

            stderr_thread.join(timeout=30)
            view_stderr_thread.join(timeout=30)
            sort_stderr_thread.join(timeout=30)

            logger.info("BWA マッピング完了: %s → %s", label, bam_out.name)
            return bam_out
        except subprocess.CalledProcessError as e:
            logger.error("BWA マッピング失敗: %s → %s: %s", label, bam_out.name, e)
            return None

    # ------------------------------------------------------------------
    # samtools merge (ラン内 collapsed + PE マージ用)
    # ------------------------------------------------------------------

    def _merge_bams(self, bam_files: List[Path], merged_bam: Path) -> bool:
        bam_paths = [str(b) for b in bam_files if b is not None and b.exists()]
        if not bam_paths:
            logger.error("統合対象の BAM がありません")
            return False

        if len(bam_paths) == 1:
            shutil.move(bam_paths[0], str(merged_bam))
            return True

        merge_cmd = ["samtools", "merge", "-f", str(merged_bam), *bam_paths]
        try:
            log_command_start(logger, merge_cmd, "samtools merge")
            subprocess.run(merge_cmd, check=True)
            logger.info("BAM マージ完了: %s", merged_bam.name)
            return True
        except subprocess.CalledProcessError as e:
            logger.error("samtools merge 失敗: %s", e)
            return False
