import random
import subprocess
import os
from helper.config import PATHS, TOOLS, PARAMETERS
from helper.logger import setup_logger
from helper.path_define import get_vcf_ref, ground_truth_vcf

logger = setup_logger(os.path.join(PATHS["logs"], "file_utils.log"))

def extract_lane1_fq(fq_path, r1_path):
    logger.info(f"Extracting land 1 for {fq_path}...")
    if os.path.exists(r1_path):
        logger.info("Land 1 fastq already exists.")
        return
    if not os.path.exists(fq_path):
        logger.error(f"Sample {fq_path} cannot read.")
        raise RuntimeError(f"Failed to read sample: {fq_path}")
    
    cmd = f'{TOOLS["seqkit"]} grep -rp ":1:" {fq_path} -o {r1_path}'
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f"Extracted land 1 for {fq_path}")


def save_results_to_csv(file_path, df):
    os.makedirs(os.path.dirname(file_path), exist_ok=True) 
    df.to_csv(file_path, index=True)
    logger.info(f"Saved: {file_path}")


def create_vcf_list(from_dir, info=None):
    vcf_list = os.path.join(from_dir, f"{info}.vcf.list")
    with open(vcf_list, "w") as f:
        for root, _, files in os.walk(from_dir):
            for file in files:
                if file.endswith(".vcf.gz") and (info==None or f"{info}" in file) and (not "all" in file):
                    f.write(os.path.join(root, file) + "\n")
    logger.info(f"VCF list saved at {vcf_list}")
    return vcf_list


def merge_vcf_list(vcf_list, merged_vcf):
    command = [
        TOOLS["bcftools"], "concat",
        "--threads", f"{PARAMETERS['threads']}",
        "-a", "--rm-dups", "all",
        "-O", "z",
        "-o", merged_vcf
    ] + [line.strip() for line in open(vcf_list)]

    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"VCF merge failed: {process.stderr}")
        raise RuntimeError(f"VCF merge failed: {process.stderr}")
    logger.info(f"Merged VCF file created at {merged_vcf}")
    return merged_vcf


def extract_vcf(sample_name, chromosome):
    vcf_reference = get_vcf_ref(chromosome)
    output_vcf_path = ground_truth_vcf(sample_name, chromosome)

    # Xây dựng lệnh bcftools
    vcf_command = [
        TOOLS["bcftools"], "view", vcf_reference,
        "--samples", sample_name,
        "-m", "2", "-M", "2",
        "-Oz", "-o", output_vcf_path,
        f"--threads={PARAMETERS['threads']}"
    ]

    try:
        result = subprocess.run(vcf_command, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to extract VCF: {result.stderr}")
            raise RuntimeError(f"Failed to extract VCF: {result.stderr}")
        logger.info(f"VCF file extracted successfully to {output_vcf_path}.")
    except Exception as e:
        logger.error(f"Error extracting VCF: {e}")
        raise
