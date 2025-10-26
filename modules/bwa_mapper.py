import subprocess
import logging
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)

class BWAMapper:
    def __init__(self, config):
        self.config = config
    
    def run_mapping_pipeline(self, sample_acc, fastq_files):
        """BWA mapping pipeline (paired-end or single-end)"""
        logger.info(f"Starting mapping pipeline for {sample_acc}")
        
        sample_dir = self.config.results_dir / sample_acc
        sample_dir.mkdir(parents=True, exist_ok=True)

        bam_dir = sample_dir / "bam_files"
        bam_dir.mkdir(parents=True, exist_ok=True)
        
        # ペアエンド判定
        if len(fastq_files) >= 2:
            fastq1, fastq2 = fastq_files[0], fastq_files[1]
            # AdapterRemoval (ペアエンド)
            logger.info(f"Running AdapterRemoval (paired-end) for {sample_acc}")
            temp_prefix = self.config.temp_dir / sample_acc
            adapter_removal_cmd = [
                "AdapterRemoval",
                "--file1", str(fastq1),
                "--file2", str(fastq2),
                "--trimns", "--trimqualities",
                "--minquality", "25", "--minlength", "25",
                "--basename", str(temp_prefix),
                "--threads", str(self.config.args.threads)
            ]
            try:
                subprocess.run(adapter_removal_cmd, check=True)
                pair1_file = temp_prefix.with_suffix(".pair1.truncated")
                pair2_file = temp_prefix.with_suffix(".pair2.truncated")
                if pair1_file.exists() and pair2_file.exists():
                    shutil.move(str(pair1_file), str(bam_dir / f"{sample_acc}.pair1.truncated"))
                    shutil.move(str(pair2_file), str(bam_dir / f"{sample_acc}.pair2.truncated"))
                else:
                    logger.error(f"AdapterRemoval output files not found for {sample_acc}")
                    return None
            except subprocess.CalledProcessError as e:
                logger.error(f"AdapterRemoval failed for {sample_acc}: {e}")
                return None
            # BWA MEM (ペアエンド)
            bam_file = bam_dir / f"{sample_acc}.bam"
            bwa_cmd = [
                "bwa", "mem", "-t", str(self.config.args.threads), "-K", "100000000", "-Y",
                "-R", f"@RG\\tID:{sample_acc}\\tSM:{sample_acc}\\tPL:ILLUMINA",
                str(self.config.reference_genome),
                str(bam_dir / f"{sample_acc}.pair1.truncated"),
                str(bam_dir / f"{sample_acc}.pair2.truncated")
            ]
        elif len(fastq_files) == 1:
            fastq1 = fastq_files[0]
            # AdapterRemoval (シングルエンド)
            logger.info(f"Running AdapterRemoval (single-end) for {sample_acc}")
            temp_prefix = self.config.temp_dir / sample_acc
            adapter_removal_cmd = [
                "AdapterRemoval",
                "--file1", str(fastq1),
                "--trimns", "--trimqualities",
                "--minquality", "25", "--minlength", "25",
                "--basename", str(temp_prefix),
                "--threads", str(self.config.args.threads)
            ]
            try:
                subprocess.run(adapter_removal_cmd, check=True)
                single_file = temp_prefix.with_suffix(".truncated")
                if single_file.exists():
                    shutil.move(str(single_file), str(bam_dir / f"{sample_acc}.truncated"))
                else:
                    logger.error(f"AdapterRemoval output file not found for {sample_acc}")
                    return None
            except subprocess.CalledProcessError as e:
                logger.error(f"AdapterRemoval failed for {sample_acc}: {e}")
                return None
            # BWA MEM (シングルエンド)
            bam_file = bam_dir / f"{sample_acc}.bam"
            bwa_cmd = [
                "bwa", "mem", "-t", str(self.config.args.threads), "-K", "100000000", "-Y",
                "-R", f"@RG\\tID:{sample_acc}\\tSM:{sample_acc}\\tPL:ILLUMINA",
                str(self.config.reference_genome),
                str(bam_dir / f"{sample_acc}.truncated")
            ]
        else:
            logger.error(f"No FASTQ files found for {sample_acc}")
            return None

        # BWA MEM 実行（共通）
        try:
            with open(bam_file, 'wb') as bam_out:
                subprocess.run(bwa_cmd, stdout=bam_out, check=True)
            logger.info(f"BWA mapping completed for {sample_acc}")
            return bam_file
        except subprocess.CalledProcessError as e:
            logger.error(f"BWA mapping failed for {sample_acc}: {e}")
            return None 

    def validate_config(self):
        """設定の妥当性を検証"""
        if not self.config.reference_genome.exists():
            raise ValueError(f"Reference genome not found: {self.config.reference_genome}")
        
        # 必要なツールの存在確認
        required_tools = ['bwa', 'samtools', 'gatk', 'qualimap']
        for tool in required_tools:
            if not shutil.which(tool):
                raise ValueError(f"Required tool not found: {tool}") 