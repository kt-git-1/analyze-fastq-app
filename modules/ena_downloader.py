import requests
import logging
from ftplib import FTP
from urllib.parse import urlparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logger = logging.getLogger(__name__)

class ENADownloader:
    def __init__(self, config):
        self.config = config
    
    def get_api_response(self, project_accession, session):
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project_accession}&result=read_run&fields=sample_accession,submitted_ftp&format=tsv"
        try:
            response = session.get(url, timeout=10)
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            logger.error(f"Error fetching data for project {project_accession}: {e}")
            sys.exit(1)
    
    def parse_response_data(self, response_data):
        sample_to_ftp_urls = {}
        lines = response_data.strip().split('\n')[1:]
        for line in lines:
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            sample_acc, ftp_urls = parts[0], parts[1]
            for ftp_url in ftp_urls.split(';'):
                ftp_url = ftp_url.strip()
                if not ftp_url:
                    continue
                sample_to_ftp_urls.setdefault(sample_acc, []).append(ftp_url)
        return sample_to_ftp_urls
    
    def download_from_ftp(self, ftp_url, destination):
        if not ftp_url.startswith('ftp://'):
            ftp_url = 'ftp://' + ftp_url
        try:
            parse = urlparse(ftp_url)
            ftp_server = parse.netloc
            ftp_path = parse.path
            filename = os.path.basename(ftp_path)
            logger.info(f"Downloading {filename} to {destination}")
            with FTP(ftp_server) as ftp:
                ftp.login()
                ftp.cwd(os.path.dirname(ftp_path))
                with open(destination, 'wb') as f:
                    ftp.retrbinary('RETR ' + filename, f.write)
            logger.info(f"Downloaded {filename} to {destination}")
            return destination
        except Exception as e:
            logger.error(f"Error downloading from {ftp_url}: {e}")
            raise
    
    def download_sample_data(self, sample_acc, ftp_urls):
        """Download all FASTQ files for a sample"""
        sample_dir = self.config.raw_data_dir / sample_acc
        sample_dir.mkdir(exist_ok=True)
        fastq_files = []
        
        with ThreadPoolExecutor(max_workers=self.config.args.workers) as executor:
            future_to_url = {
                executor.submit(
                    self.download_from_ftp,
                    url,
                    sample_dir / os.path.basename(urlparse(url if url.startswith('ftp://') else 'ftp://' + url).path)
                ): url
                for url in ftp_urls if url.endswith('.gz')
            }
            
            for future in as_completed(future_to_url):
                url = future_to_url[future]
                try:
                    dest = future.result()
                    fastq_files.append(dest)
                except Exception as e:
                    logger.error(f"Failed download {url}: {e}")
        
        return fastq_files 