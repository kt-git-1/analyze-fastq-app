import os
import shlex
import subprocess
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class MapDamageAnalyzer:
    def __init__(self, config):
        self.config = config
    
    def run_mapdamage(self, sample_acc, softclipped_bam):
        """Run mapDamage analysis"""
        logger.info(f"mapDamageを実行します: {sample_acc}")
        
        sample_outdir = self.config.results_dir / sample_acc / "mapdamage"
        sample_outdir.mkdir(parents=True, exist_ok=True)
        
        # Filter BAM (mapQ>=30, POS>=300)
        filtered_bam = self.config.temp_dir / f"{sample_acc}_filtered.bam"
        sorted_filtered_bam = self.config.temp_dir / f"{sample_acc}_filtered.sorted.bam"
        
        # Filter and sort using an awk pipeline
        try:
            threshold = 300
            filter_cmd = (
                f"samtools view -h -q 30 {shlex.quote(str(softclipped_bam))} | "
                f"awk 'BEGIN{{OFS=\"\\t\"}} /^@/ || $4>={threshold}' | "
                f"samtools view -b - | "
                f"samtools sort -o {shlex.quote(str(sorted_filtered_bam))}"
            )
            subprocess.run(filter_cmd, shell=True, executable='/bin/bash', check=True)
            # インデックス作成（BAI）
            subprocess.run(["samtools", "index", str(sorted_filtered_bam)], check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"フィルタリングに失敗しました: {sample_acc}: {e}")
            return None
        
        # Run mapDamage
        mapdamage_cmd = [
            "mapDamage", "-i", str(sorted_filtered_bam), 
            "-r", str(self.config.reference_genome), 
            "-d", str(sample_outdir), "--merge-libraries"
        ]
        
        try:
            subprocess.run(mapdamage_cmd, check=True)
            logger.info(f"mapDamageが完了しました: {sample_acc}")
            return sample_outdir
        except subprocess.CalledProcessError as e:
            logger.error(f"mapDamageに失敗しました: {sample_acc}: {e}")
            return None

class QualimapAnalyzer:
    def __init__(self, config):
        self.config = config
    
    def run_qualimap(self, sample_acc, dedup_bam):
        """
        Run Qualimap quality check
        
        - スキップ条件: BAM ファイルが存在しない／空サイズ、またはマッピングされたリード数がゼロ
        - 出力ディレクトリ: <results_dir>/<sample_acc>/qualimap
        - 成功時は出力ディレクトリの Path を、スキップまたは失敗時は None を返す
        """

        logger.info(f"Qualimapを実行します: {sample_acc}")
        
        # 1. ファイル存在とサイズのチェック
        try:
            if dedup_bam is None or not dedup_bam.exists() or dedup_bam.stat().st_size == 0:
                logger.warning(f"Qualimapをスキップします: {sample_acc}: BAMファイルが見つかりませんまたは空です ({dedup_bam})")
                return None
        except Exception as e:
            logger.error(f"BAMファイルのステータスを取得に失敗しました: {sample_acc}: {e}")
            return None

        # 2. マッピングされたリード数のチェック
        try:
            idxstats = subprocess.run(
                ["samtools", "idxstats", str(dedup_bam)],
                capture_output=True, text=True, check=True
            )
            mapped_reads = sum(
                int(line.split()[2]) for line in idxstats.stdout.strip().split("\n") if line
            )
            if mapped_reads == 0:
                logger.warning(f"Qualimapをスキップします: {sample_acc}: {dedup_bam} にマッピングされたリードがありません")
                return None
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools idxstatsに失敗しました: {sample_acc}: {e}")
            return None
        except FileNotFoundError:
            logger.error("samtoolsが見つかりません。samtoolsがインストールされているか、PATHに追加されているか確認してください。")
            return None

        # 3. 出力ディレクトリを作成
        sample_outdir = self.config.results_dir / sample_acc / "qualimap"
        sample_outdir.mkdir(parents=True, exist_ok=True)

        # 4. Qualimap 実行
        qualimap_cmd = [
            "qualimap", "bamqc",
            "-bam", str(dedup_bam),
            "-outdir", str(sample_outdir),
            "-outformat", "HTML",
            "--java-mem-size=8G",
        ]
        env = os.environ.copy()
        # Java 9+ no longer recognizes -XX:MaxPermSize, which Qualimap still adds
        # in its launch script. The following option tells the JVM to ignore
        # any unsupported -XX arguments so Qualimap can run under Java 17.
        env["JAVA_TOOL_OPTIONS"] = env.get("JAVA_TOOL_OPTIONS", "") + " -XX:+IgnoreUnrecognizedVMOptions"
        try:
            subprocess.run(qualimap_cmd, check=True, env=env)
            logger.info(f"Qualimapが完了しました: {sample_acc}")
            return sample_outdir
        except subprocess.CalledProcessError as e:
            logger.error(f"Qualimapに失敗しました: {sample_acc}: {e}")
            return None

class HaplotypeCaller:
    def __init__(self, config):
        self.config = config
    
    def run_haplotypecaller(self, sample_acc, dedup_bam):
        """Run GATK HaplotypeCaller"""
        logger.info(f"HaplotypeCallerを実行します: {sample_acc}")
        
        vcf_dir = self.config.results_dir / sample_acc / "vcf_files"
        vcf_dir.mkdir(parents=True, exist_ok=True)
        vcf_file = vcf_dir / f"{sample_acc}.vcf"
        
        haplotypecaller_cmd = [
            "gatk", "HaplotypeCaller",
            "-R", str(self.config.reference_genome),
            "-I", str(dedup_bam),
            "-O", str(vcf_file),
            "--output-mode", "EMIT_VARIANTS_ONLY",
            "-stand-call-conf", "30",
            "--native-pair-hmm-threads", str(self.config.args.threads)
        ]
        
        try:
            subprocess.run(haplotypecaller_cmd, check=True)
            logger.info(f"HaplotypeCallerが完了しました: {sample_acc}")
            return vcf_file
        except subprocess.CalledProcessError as e:
            logger.error(f"HaplotypeCallerに失敗しました: {sample_acc}: {e}")
            return None