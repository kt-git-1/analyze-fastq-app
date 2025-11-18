import subprocess
import logging
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)

class BAMProcessor:
    def __init__(self, config):
        self.config = config
    
    def run_bam_processing(self, sample_acc, softclipped_bam):
        """Sort, deduplicate, and index BAM file"""
        logger.info(f"Running BAM processing for {sample_acc}")
        base_name = sample_acc
        sample_bam_dir = self.config.results_dir / sample_acc / "bam_files"
        sample_dedup_dir = self.config.results_dir / sample_acc / "dedup"
        sample_bam_dir.mkdir(parents=True, exist_ok=True)
        sample_dedup_dir.mkdir(parents=True, exist_ok=True)

        def run_cmd(cmd, step_name):
            try:
                logger.info(f"Running: {' '.join(str(c) for c in cmd)}")
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Error in {step_name}: {e}")
                if hasattr(e, 'output'):
                    logger.error(f"Command output: {e.output}")
                raise

        # Step 1: CleanSam
        clean_bam = sample_bam_dir / f"{base_name}.clean.bam"
        run_cmd([
            "java", "-Xmx" + self.config.args.java_mem, "-jar", "/usr/local/bin/picard.jar",
            "CleanSam", "I=" + str(softclipped_bam), "O=" + str(clean_bam),
            "VALIDATION_STRINGENCY=LENIENT"
        ], "Picard CleanSam")

        # Step 2: AddOrReplaceReadGroups
        grouped_bam = sample_bam_dir / f"{base_name}.grouped.bam"
        run_cmd([
            "java", "-Xmx" + self.config.args.java_mem, "-jar", "/usr/local/bin/picard.jar",
            "AddOrReplaceReadGroups", "I=" + str(clean_bam), "O=" + str(grouped_bam),
            "RGLB=Bayanbulag", "RGSM=" + base_name, "RGPU=tile", "RGPL=ILLUMINA",
            "RGID=" + base_name, "RGDS=" + base_name, "RGCN=KNZWUNIV",
            "VALIDATION_STRINGENCY=LENIENT"
        ], "Picard AddOrReplaceReadGroups")

        # Step 3: MarkDuplicates (REMOVE_DUPLICATESでdedupも同時に)
        marked_bam = sample_bam_dir / f"{base_name}.marked.bam"
        metrics_file = sample_bam_dir / f"{base_name}.marked_dup_metrics.txt"
        run_cmd([
            "java", "-Xmx" + self.config.args.java_mem, "-jar", "/usr/local/bin/picard.jar",
            "MarkDuplicates", "I=" + str(grouped_bam), "O=" + str(marked_bam),
            "M=" + str(metrics_file), "REMOVE_DUPLICATES=true", "VALIDATION_STRINGENCY=LENIENT"
        ], "Picard MarkDuplicates")

        # Step 4: Sort deduplicated BAM
        final_dedup_bam = sample_dedup_dir / f"{base_name}.marked.dedup.sorted.bam"
        run_cmd([
            "samtools", "sort", "-o", str(final_dedup_bam), str(marked_bam)
        ], "samtools sort (dedup)")

        # Step 5: Index
        run_cmd(["samtools", "index", str(final_dedup_bam)], "samtools index")

        logger.info(f"BAM processingが完了しました: {sample_acc}")
        return final_dedup_bam