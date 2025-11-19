import subprocess
import logging
import shutil
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

class BWAMapper:
    def __init__(self, config):
        self.config = config
    
    def run_mapping_pipeline(self, sample_acc, fastq_files):
        """
        BWAマッピングパイプラインを実行します。
        ペアエンドかシングルエンドかを判定して、AdapterRemovalとBWA MEMを実行します。
        """
        logger.info(f"マッピングパイプラインを開始します: {sample_acc}")
        
        sample_dir = self.config.results_dir / sample_acc
        sample_dir.mkdir(parents=True, exist_ok=True)

        bam_dir = sample_dir / "bam_files"
        bam_dir.mkdir(parents=True, exist_ok=True)
        
        # ペアエンド判定
        # pair-end
        if len(fastq_files) >= 2:
            fastq1, fastq2 = fastq_files[0], fastq_files[1]
            # AdapterRemoval (ペアエンド)
            logger.info(f"AdapterRemoval (ペアエンド, data_type={self.config.data_type})を実行します: {sample_acc}")
            temp_prefix = self.config.temp_dir / sample_acc
            adapter_removal_cmd = self.build_adapter_removal_cmd_paired(fastq1, fastq2, temp_prefix)
            try:
                subprocess.run(adapter_removal_cmd, check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"AdapterRemovalに失敗しました: {sample_acc}: {e}")
                return None
            # BWA MEM (ペアエンド)
            # ancient
            if self.config.data_type == "ancient":
                collapsed_file = temp_prefix.with_suffix(".collapsed.truncated")
                R1_file = temp_prefix.with_suffix(".pair1.truncated")
                R2_file = temp_prefix.with_suffix(".pair2.truncated")
                partial_bams = []
                # collapsed
                if collapsed_file.exists():
                    collapsed_bam = bam_dir / f"{sample_acc}.collapsed.sorted.bam"
                    bam = self.run_bwa_single_end_and_sort(collapsed_file, collapsed_bam, sample_acc)
                    if bam is not None:
                        partial_bams.append(bam)
                else:
                    logger.warning(f"collapsedファイルが見つかりませんでした: {collapsed_file}")
                # non-collapsed R1
                if R1_file.exists():
                    R1_bam = bam_dir / f"{sample_acc}.R1.sorted.bam"
                    bam = self.run_bwa_single_end_and_sort(R1_file, R1_bam, sample_acc)
                    if bam is not None:
                        partial_bams.append(bam)
                else:
                    logger.warning(f"R1ファイルが見つかりませんでした: {R1_file}")
                # non-collapsed R2
                if R2_file.exists():
                    R2_bam = bam_dir / f"{sample_acc}.R2.sorted.bam"
                    bam = self.run_bwa_single_end_and_sort(R2_file, R2_bam, sample_acc)
                    if bam is not None:
                        partial_bams.append(bam)
                else:
                    logger.warning(f"R2ファイルが見つかりませんでした: {R2_file}")

                if not partial_bams:
                    logger.error(f"BWAマッピングに失敗しました: {sample_acc}")
                    return None

                merged_bam = bam_dir / f"{sample_acc}.merged.sorted.bam"
                if not self.merge_bams(partial_bams, merged_bam):
                    logger.error(f"BAMマージに失敗しました: {sample_acc}")
                    return None

                return merged_bam

            # modern
            else:
                temp_R1 = temp_prefix.with_suffix(".R1.truncated")
                temp_R2 = temp_prefix.with_suffix(".R2.truncated")
                if not (temp_R1.exists() and temp_R2.exists()):
                    logger.error(f"AdapterRemoval出力ファイルが見つかりませんでした: {sample_acc}")  
                    return None

                dest_R1 = bam_dir / f"{sample_acc}.R1.truncated"
                dest_R2 = bam_dir / f"{sample_acc}.R2.truncated"
                shutil.move(str(temp_R1), str(dest_R1))
                shutil.move(str(temp_R2), str(dest_R2))

                bam_file = bam_dir / f"{sample_acc}.sorted.bam"
                bam = self.run_bwa_paired_and_sort(dest_R1, dest_R2, bam_file, sample_acc)
                if bam is None:
                    logger.error(f"BWAマッピングに失敗しました: {sample_acc}")
                return bam

        # single-end
        elif len(fastq_files) == 1:
            fastq1 = fastq_files[0]
            # AdapterRemoval (シングルエンド)
            logger.info(f"AdapterRemoval (シングルエンド, data_type={self.config.data_type})を実行します: {sample_acc}")
            temp_prefix = self.config.temp_dir / sample_acc
            adapter_removal_cmd = self.build_adapter_removal_cmd_single(fastq1, temp_prefix)
            try:
                subprocess.run(adapter_removal_cmd, check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"AdapterRemovalに失敗しました: {sample_acc}: {e}")
                return None
            # BWA MEM (シングルエンド)
            temp_single = temp_prefix.with_suffix(".truncated")
            if not temp_single.exists():
                logger.error(f"AdapterRemoval出力ファイルが見つかりませんでした: {sample_acc}")
                return None

            dest_single = bam_dir / f"{sample_acc}.truncated"
            shutil.move(str(temp_single), str(dest_single))

            bam_file = bam_dir / f"{sample_acc}.sorted.bam"
            bam = self.run_bwa_single_end_and_sort(dest_single, bam_file, sample_acc)
            if bam is None:
                logger.error(f"BWAマッピングに失敗しました: {sample_acc}")
            return bam

        # no FASTQ
        else:
            logger.error(f"FASTQファイルが見つかりませんでした: {sample_acc}")
            return None

    def validate_config(self):
        """
        設定の妥当性を検証
        """
        if not self.config.reference_genome.exists():
            raise ValueError(f"Reference genomeが見つかりません: {self.config.reference_genome}")
        
        # 必要なツールの存在確認
        required_tools = ['bwa', 'samtools', 'gatk', 'qualimap']
        for tool in required_tools:
            if not shutil.which(tool):
                raise ValueError(f"必要なツールが見つかりません: {tool}") 

    def build_adapter_removal_cmd_paired(self, fastq1, fastq2, temp_prefix):
        """
        ancient/modernに応じてAdapterRemoval (ペアエンド)のコマンドを構築します。
        """
        base_cmd = [
            "AdapterRemoval",
            "--file1", str(fastq1),
            "--file2", str(fastq2),
            "--trimns", "--trimqualities",
            "--basename", str(temp_prefix),
            "--threads", str(self.config.args.threads)
        ]

        if self.config.data_type == "ancient":
            base_cmd.extend([
                "--minquality", "20", 
                "--minlength", "30",
                "--collapse"
            ])
        elif self.config.data_type == "modern":
            base_cmd.extend([
                "--minquality", "25", 
                "--minlength", "25"
            ])
        return base_cmd
        
    def build_adapter_removal_cmd_single(self, fastq1, temp_prefix):
        """
        ancient/modernに応じてAdapterRemoval (シングルエンド)のコマンドを構築します。
        """
        base_cmd = [
            "AdapterRemoval",
            "--file1", str(fastq1),
            "--trimns", "--trimqualities",
            "--basename", str(temp_prefix),
            "--threads", str(self.config.args.threads)
        ]

        if self.config.data_type == "ancient":
            base_cmd.extend([
                "--minquality", "20", 
                "--minlength", "30",
            ])
        elif self.config.data_type == "modern":
            base_cmd.extend([
                "--minquality", "25", 
                "--minlength", "25"
            ])
        return base_cmd


    def run_bwa_single_end_and_sort(self, fastq: Path, bam_out: Path, sample_acc: str) -> Optional[Path]:
        """
        single-end FASTQ を bwa mem + samtools view + samtools sort で
        sorted BAM にする。成功したら bam_out を返し、失敗したら None。
        """
        logger.info(f"single-end BWAマッピングを実行します: {fastq.name} → {bam_out.name}")
        bwa_cmd = [
            "bwa", "mem",
            "-t", str(self.config.args.threads),
            "-K", "100000000",
            "-Y",
            "-R", f"@RG\\tID:{sample_acc}\\tSM:{sample_acc}\\tPL:ILLUMINA",
            str(self.config.reference_genome),
            str(fastq),
        ]
        # bwa mem → samtools view → samtools sort
        view_cmd = ["samtools", "view", "-b", "-"]
        sort_cmd = ["samtools", "sort", "-o", str(bam_out), "-"]

        try:
            # bwa → view → sort をパイプでつなぐ
            bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
            view_proc = subprocess.Popen(view_cmd, stdin=bwa_proc.stdout, stdout=subprocess.PIPE)
            bwa_proc.stdout.close()
            sort_proc = subprocess.Popen(sort_cmd, stdin=view_proc.stdout)
            view_proc.stdout.close()

            sort_proc.wait()
            if sort_proc.returncode != 0:
                raise subprocess.CalledProcessError(sort_proc.returncode, sort_cmd)

            logger.info(f"single-end BWAマッピングが完了しました: {fastq.name} → {bam_out.name}")
            return bam_out
        except subprocess.CalledProcessError as e:
            logger.error(f"single-end BWAマッピングに失敗しました: {fastq.name} → {bam_out.name}: {e}")
            return None

    def run_bwa_paired_and_sort(self, fastq1: Path, fastq2: Path, bam_out: Path, sample_acc: str) -> Optional[Path]:
        """
        pair-end FASTQ を bwa mem + samtools view + samtools sort で
        sorted BAM にする。
        """
        logger.info(f"pair-end BWAマッピングを実行します: {fastq1.name}, {fastq2.name} → {bam_out.name}")
        bwa_cmd = [
            "bwa", "mem",
            "-t", str(self.config.args.threads),
            "-K", "100000000",
            "-Y",
            "-R", f"@RG\\tID:{sample_acc}\\tSM:{sample_acc}\\tPL:ILLUMINA",
            str(self.config.reference_genome),
            str(fastq1),
            str(fastq2),
        ]
        view_cmd = ["samtools", "view", "-b", "-"]
        sort_cmd = ["samtools", "sort", "-o", str(bam_out), "-"]

        try:
            # bwa → view → sort をパイプでつなぐ
            bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
            view_proc = subprocess.Popen(view_cmd, stdin=bwa_proc.stdout, stdout=subprocess.PIPE)
            bwa_proc.stdout.close()
            sort_proc = subprocess.Popen(sort_cmd, stdin=view_proc.stdout)
            view_proc.stdout.close()

            sort_proc.wait()
            if sort_proc.returncode != 0:
                raise subprocess.CalledProcessError(sort_proc.returncode, sort_cmd)

            logger.info(f"pair-end BWAマッピングが完了しました: {fastq1.name}, {fastq2.name} → {bam_out.name}")
            return bam_out
        except subprocess.CalledProcessError as e:
            logger.error(f"pair-end BWAマッピングに失敗しました: {fastq1.name}, {fastq2.name} → {bam_out.name}: {e}")
            return None

    def merge_bams(self, bam_files, merged_bam: Path) -> bool:
        """
        samtools merge で複数の BAM を統合
        """
        bam_files = [str(b) for b in bam_files if b is not None and Path(b).exists()]
        if not bam_files:
            logger.error("統合対象のBAMファイルがありません")
            return False

        merge_cmd = [
            "samtools", "merge",
            "-f",
            str(merged_bam),
            *bam_files,
        ]

        logger.info(f"samtools merge を実行します: {' '.join(merge_cmd)}")

        try:
            subprocess.run(merge_cmd, check=True)
            logger.info(f"統合BAMを作成しました: {merged_bam}")
            
            logger.info(f"統合BAMのインデックスを作成します: {merged_bam}")
            index_cmd = ["samtools", "index", str(merged_bam)]
            subprocess.run(index_cmd, check=True)
            logger.info(f"統合BAMのインデックスを作成しました: {merged_bam}")
            
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools merge に失敗しました: {e}")
            return False