import subprocess
import os
import shutil
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import basevar_outdir, bamlist_dir, basevar_vcf, samid
from helper.logger import setup_logger
from concurrent.futures import ThreadPoolExecutor
from helper.file_utils import create_vcf_list, merge_vcf_list

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
            log_file = f"{outdir}/{outfile_prefix}.log"

            print(f"Starting BaseVar for {fq} region {region}")

            try:
                with open(log_file, "w") as log:
                    subprocess.run([TOOLS['basevar'], "basetype",
                        "-R", REF,
                        "-L", bamlist_path,
                        "-r", region,
                        "--min-af=0.001",
                        "--output-vcf", f"{outdir}/{outfile_prefix}.vcf.gz",
                        "--output-cvg", f"{outdir}/{outfile_prefix}.cvg.tsv.gz",
                        "--smart-rerun"
                    ], stdout=log, stderr=subprocess.STDOUT, check=True)
                    print(f"Done BaseVar for region {region}")

            except Exception as e:
                logger.error(f"Unexpected error occurred for region {region}: {e}")
                continue    


def merge_vcf_files(fq, chromosome):
    print(f"Creating vcf_list for {fq} {chromosome}")
    vcf_list = create_vcf_list(basevar_outdir(fq), f"{chromosome}")
    merged_vcf = basevar_vcf(fq, chromosome)

    print(f"Merging vcf_list for {fq} {chromosome}")
    merge_vcf_list(vcf_list, merged_vcf)

    print(f"Indexing VCF file {merged_vcf}")
    subprocess.run([TABIX, "-f", merged_vcf], check=True)


def run_basevar_chr(fq, chromosome): 
    finish_flag = os.path.join(f"{basevar_outdir(fq)}_final", f"basevar_{chromosome}.finish")

    # Step 0: Verify the flag
    if os.path.exists(finish_flag):
        print(f"Basevar result for {chromosome} {samid(fq)} already exist. Skip basevar for {chromosome}...")
        return

    print(f"Run basevar for {fq} {chromosome}")

    # Step 1: Run BaseVar for the chromosome
    run_basevar_step(fq, chromosome)

    # Step 2: Merge VCF files for the chromosome
    merge_vcf_files(fq, chromosome)

    # Create finish flag
    with open(finish_flag, "w") as flag:
        flag.write(f"Completed BaseVar for {fq} {chromosome}.")
    print(f"Completed BaseVar for {fq} chromosome {chromosome}")


def run_basevar(fq):
    os.makedirs(f"{basevar_outdir(fq)}_final",exist_ok=True)
    os.makedirs(f"{basevar_outdir(fq)}",exist_ok=True)

    with ThreadPoolExecutor(max_workers=PARAMETERS['threads']) as executor:
        executor.map(lambda chr: run_basevar_chr(fq, chr), PARAMETERS["chrs"])

    #shutil.rmtree(basevar_outdir(fq))
    print(f"Temporary directory {basevar_outdir(fq)} deleted.")

    print(f"Completed BaseVar pipeline for {fq}.")

