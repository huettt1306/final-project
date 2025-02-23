import subprocess
import os
import shutil
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import basevar_outdir, bamlist_dir, vcf_list_path, basevar_vcf
from helper.logger import setup_logger
from concurrent.futures import ThreadPoolExecutor

# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "basevar_pipeline.log"))


REF = PATHS["ref"]
REF_FAI = PATHS["ref_fai"]
BCFTOOLS = TOOLS["bcftools"]
TABIX = TOOLS["tabix"]
DELTA = PARAMETERS["basevar"]["delta"]

def load_reference_fai(in_fai, chroms=None):
    ref = []
    with open(in_fai) as fh:
        for r in fh:
            col = r.strip().split()
            if chroms is not None and len(chroms):
                if col[0] in chroms:
                    ref.append([col[0], 1, int(col[1])])
            else:
                ref.append([col[0], 1, int(col[1])])
    return ref


def run_basevar_step(fq, chromosome):
    ref_fai = load_reference_fai(REF_FAI, [chromosome])
    bamlist_path = bamlist_dir(fq)
    outdir = basevar_outdir(fq)

    print(f"Run basevar step for {fq} {chromosome} with ref_fai {ref_fai} and delta {DELTA}")

    for chr_id, reg_start, reg_end in ref_fai:
        for i in range(reg_start - 1, reg_end, DELTA):
            start = i + 1
            end = min(i + DELTA, reg_end)
            region = f"{chr_id}:{start}-{end}"
            outfile_prefix = f"{chr_id}_{start}_{end}"

            print(f"Starting BaseVar for {fq} region {region}")
            
            command = [
                TOOLS['basevar'], "basetype",
                "-t", f"{PARAMETERS['threads']}",
                "-R", REF,
                "-L", bamlist_path,
                "-r", region,
                "--min-af=0.001",
                "--output-vcf", f"{outdir}/{outfile_prefix}.vcf.gz",
                "--output-cvg", f"{outdir}/{outfile_prefix}.cvg.tsv.gz",
                "--smart-rerun"
            ]

            log_file = f"{outdir}/{outfile_prefix}.log"
            try:
                with open(log_file, "w") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
                    logger.info(f"Done BaseVar for region {region}")
            except subprocess.CalledProcessError as e:
                logger.error(f"BaseVar failed for region {region} with error: {e}")
                continue  # B? qua l?i và ti?p t?c pipeline
            except Exception as e:
                logger.error(f"Unexpected error occurred for region {region}: {e}")
                continue  # B? qua l?i và ti?p t?c pipeline


def create_vcf_list(fq, chromosome):
    vcf_list = vcf_list_path(fq, chromosome)
    print(f"Creating VCF list for chromosome {chromosome}")
    with open(vcf_list, "w") as f:
        for root, _, files in os.walk(basevar_outdir(fq)):
            for file in files:
                if file.endswith(".vcf.gz") and f"{chromosome}_" in file:
                    f.write(os.path.join(root, file) + "\n")
    print(f"VCF list saved at {vcf_list}")
    return vcf_list

def merge_vcf_files(fq, chromosome):
    vcf_list_path = create_vcf_list(fq, chromosome)
    merged_vcf = basevar_vcf(fq, chromosome)

    command = [
        BCFTOOLS, "concat",
        "--threads", f"{PARAMETERS['threads']}",
        "-a", "--rm-dups", "all",
        "-O", "z",
        "-o", merged_vcf
    ] + [line.strip() for line in open(vcf_list_path)]

    print(f"Merging VCF files for {fq} chromosome {chromosome}")
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"VCF merge failed: {process.stderr}")
        raise RuntimeError(f"VCF merge failed: {process.stderr}")
    print(f"Merged VCF file created at {merged_vcf}")
    return merged_vcf

def index_vcf_file(vcf_path):
    command = [TABIX, "-f", "-p", "vcf", vcf_path]

    print(f"Indexing VCF file {vcf_path}")
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"VCF indexing failed: {process.stderr}")
        raise RuntimeError(f"VCF indexing failed: {process.stderr}")
    print(f"VCF file indexed at {vcf_path}")


def run_basevar_chr(fq, chromosome): 
    print(f"Run basevar for {fq} {chromosome}")
    if os.path.exists(basevar_vcf(fq, chromosome)):
        print(f"Đã có kết quả basevar cho mẫu {fq} với {chromosome}.")
        return

    # Step 1: Run BaseVar for the chromosome
    run_basevar_step(fq, chromosome)

    # Step 2: Merge VCF files for the chromosome
    merged_vcf = merge_vcf_files(fq, chromosome)

    # Step 3: Index the merged VCF file
    index_vcf_file(merged_vcf)

    print(f"Completed processing for {fq} chromosome {chromosome}")


def run_basevar(fq):
    os.makedirs(f"{basevar_outdir(fq)}_final",exist_ok=True)
    os.makedirs(f"{basevar_outdir(fq)}",exist_ok=True)

    with ThreadPoolExecutor(max_workers=PARAMETERS["threads"]) as executor:
        executor.map(lambda chr: run_basevar_chr(fq, chr), PARAMETERS["chrs"])

    print(f"Completed BaseVar pipeline for {fq}")
    shutil.rmtree(basevar_outdir(fq))
    print(f"Temporary directory {basevar_outdir(fq)} deleted.")
