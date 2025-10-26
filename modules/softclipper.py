import pysam
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

logger = logging.getLogger(__name__)

class SoftClipper:
    def __init__(self, config):
        self.config = config
        self.clip_length = 5
        self.batch_size = 10000
    
    def calculate_query_length(self, cigartuples):
        if cigartuples is None:
            return 0
        return sum(length for op, length in cigartuples if op in {0, 1, 4})
    
    def process_read(self, read):
        if read.is_unmapped or not read.cigartuples:
            return read
        
        expected_length = self.calculate_query_length(read.cigartuples)
        actual_length = len(read.query_sequence) if read.query_sequence else 0
        
        if actual_length != expected_length:
            logger.warning(f"CIGAR mismatch for {read.query_name}")
        
        new_cigar = []
        total_clip = self.clip_length
        remaining_clip = self.clip_length
        
        if read.cigartuples[0][0] == 4:
            total_clip += read.cigartuples[0][1]
            read.cigartuples = read.cigartuples[1:]
        
        new_cigar.append((4, total_clip))
        
        adjusted_cigar = []
        for op, length in read.cigartuples:
            if op == 0:
                if length > remaining_clip:
                    adjusted_cigar.append((0, length - remaining_clip))
                    remaining_clip = 0
                else:
                    remaining_clip -= length
            elif op == 2:
                adjusted_cigar.append((2, length))
            else:
                adjusted_cigar.append((op, length))
        
        final_expected_length = self.calculate_query_length(new_cigar + adjusted_cigar)
        if actual_length > final_expected_length:
            read.query_sequence = read.query_sequence[:final_expected_length]
            if read.query_qualities:
                read.query_qualities = read.query_qualities[:final_expected_length]
        elif actual_length < final_expected_length:
            logger.warning(f"Length mismatch after adjustment for {read.query_name}")
            return None
        
        read.cigartuples = new_cigar + adjusted_cigar
        read.is_paired = True
        
        return read
    
    def run_softclipping(self, sample_acc, bam_file):
        """Apply softclipping to BAM file"""
        logger.info(f"Running softclipping for {sample_acc}")
        
        input_bam = bam_file
        sample_dir = self.config.results_dir / sample_acc / "softclipped"
        sample_dir.mkdir(parents=True, exist_ok=True)
        output_bam = sample_dir / f"{sample_acc}_softclipped.bam"
        
        try:
            with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
                 pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
                
                with ThreadPoolExecutor(max_workers=self.config.args.threads) as executor:
                    while True:
                        batch = [read for _, read in zip(range(self.batch_size), in_bam)]
                        if not batch:
                            break
                        
                        processed_reads = list(executor.map(self.process_read, batch))
                        valid_reads = [r for r in processed_reads if r is not None]
                        
                        for r in valid_reads:
                            out_bam.write(r)
            
            logger.info(f"Softclipping completed for {sample_acc}")
            return output_bam
        except Exception as e:
            logger.error(f"Softclipping failed for {sample_acc}: {e}")
            return None 