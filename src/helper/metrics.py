import os
import subprocess
from helper.path_define import fastq_path_lane1
from helper.config import TOOLS, PATHS
from helper.logger import setup_logger


COVERAGE_FILE = os.path.join(PATHS["fastq_directory"], "coverage.txt")
logger = setup_logger(os.path.join(PATHS["logs"], "metrics.log"))

def get_fastq_coverage(name):
    if os.path.exists(COVERAGE_FILE):
        with open(COVERAGE_FILE, 'r') as f:
            for line in f:
                sample, coverage = line.strip().split('\t')
                if sample == name:
                    return float(coverage)


    input_fastq = fastq_path_lane1(name)
    result = subprocess.run(
        f"{TOOLS['seqkit']} stats {input_fastq} | tail -n 1",
        shell=True,
        capture_output=True,
        text=True
    )

    stats_values = result.stdout.split()
    
    if len(stats_values) < 5:
        raise ValueError("Segkit error!")

    total_bases = int(stats_values[4].replace(",", ""))  

    genome_size = 3200000000
    coverage = total_bases / genome_size

    with open(COVERAGE_FILE, 'a') as f:
        f.write(f"{name}\t{coverage}\n")

    return coverage


