import subprocess
import os
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import vcf_prefix, get_vcf_path, filtered_tsv_path, filtered_vcf_path, chunks_path, norm_vcf_path, dbsnp_dir
from helper.logger import setup_logger
from concurrent.futures import ThreadPoolExecutor

logger = setup_logger(os.path.join(PATHS["logs"], "reference_panel_pipeline.log"))


BCFTOOLS = TOOLS["bcftools"]
BGZIP = TOOLS["bgzip"]
TABIX = TOOLS["tabix"]
GLIMPSE_CHUNK = TOOLS["GLIMPSE_chunk"]
reference_path = PATHS["reference_path"]

def check_reference_panel(chromosome):
    prefix = vcf_prefix(chromosome)
    required_files = [
        f"{reference_path}/{prefix}.biallelic.snp.maf0.001.vcf.gz",
        f"{reference_path}/{prefix}.biallelic.snp.maf0.001.sites.vcf.gz",
        f"{reference_path}/{prefix}.biallelic.snp.maf0.001.sites.tsv.gz",
        f"{reference_path}/{prefix}.chunks.txt"
    ]
    for file in required_files:
        if not os.path.exists(file):
            return False
    return True

def download_reference_panel(chromosome):
    url = f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/{vcf_prefix(chromosome)}.vcf.gz"
    vcf_path = get_vcf_path(chromosome)
    index_path = f"{vcf_path}.tbi"

    if os.path.exists(vcf_path) and os.path.exists(index_path):
        print(f"Reference panel for {chromosome} already exists. Skipping download.")
        return vcf_path

    print(f"Downloading reference panel for chromosome {chromosome}...")

    commands = [
        ["wget", "-c", url, "-P", reference_path],
        ["wget", "-c", f"{url}.tbi", "-P", reference_path]
    ]
    
    for command in commands:
        process = subprocess.run(command, capture_output=True, text=True)
        if process.returncode != 0:
            logger.error(f"Error downloading file: {process.stderr}")
            raise RuntimeError(f"Error downloading file: {process.stderr}")

    print(f"Downloaded reference panel for chromosome {chromosome}.")
    return vcf_path

def normalize_and_filter_reference(chromosome):
    vcf_path = get_vcf_path(chromosome)
    output_vcf = norm_vcf_path(chromosome)

    if os.path.exists(output_vcf) and os.path.exists(f"{output_vcf}.tbi"):
        print(f"Filtered VCF already exists. Skipping normalization and filtering.")
        return output_vcf

    print(f"Normalizing and filtering VCF for chromosome {chromosome}...")
    command = [
        BCFTOOLS, "norm", "-m", "-any", vcf_path, "-Ou",
        "--threads", f"{PARAMETERS['threads']}", "|",
        BCFTOOLS, "view", "-m", "2", "-M", "2", "-v", "snps", "-i", "'MAF>0.001'",
        "--threads", f"{PARAMETERS['threads']}", "-Oz", "-o", output_vcf
    ]

    process = subprocess.run(" ".join(command), shell=True, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error normalizing and filtering: {process.stderr}")
        raise RuntimeError(f"Error normalizing and filtering: {process.stderr}")

    index_command = [BCFTOOLS, "index", "-f", output_vcf]
    subprocess.run(index_command, check=True)

    print(f"Filtered VCF created at {output_vcf}.")
    return output_vcf

def process_snp_sites(chromosome):
    vcf_path = norm_vcf_path(chromosome)
    filtered_vcf = filtered_vcf_path(chromosome)
    tsv_output = filtered_tsv_path(chromosome)

    if os.path.exists(filtered_vcf) and os.path.exists(f"{filtered_vcf}.tbi") and os.path.exists(tsv_output):
        print(f"SNP site files already exist. Skipping SNP site processing.")
        return filtered_vcf, tsv_output

    print(f"Processing SNP sites for chromosome {chromosome}...")
    commands = [
        [BCFTOOLS, "view", "-G", "-m", "2", "-M", "2", "-v", "snps", vcf_path, "-Oz", "-o", filtered_vcf, "--threads", f"{PARAMETERS['threads']}",],
        [BCFTOOLS, "index", "-f", filtered_vcf],
    ]

    for command in commands:
        process = subprocess.run(command, check=True)
        if process.returncode != 0:
            logger.error(f"Error processing SNP sites: {process.stderr}")
            raise RuntimeError(f"Error processing SNP sites: {process.stderr}")

    print(f"Running bcftools query and bgzip for chromosome {chromosome}...")
    with open(tsv_output, "wb") as output_file:
        query_command = [BCFTOOLS, "query", "-f", "%CHROM\\t%POS\\t%REF,%ALT\\n", filtered_vcf]
        query_process = subprocess.Popen(query_command, stdout=subprocess.PIPE)
        subprocess.run([BGZIP, "-c"], stdin=query_process.stdout, stdout=output_file, check=True)
        query_process.stdout.close()
        query_process.wait()

    print(f"Running tabix for chromosome {chromosome}...")
    subprocess.run([TABIX, "-s1", "-b2", "-e2", tsv_output], check=True)

    print(f"Processed SNP sites for chromosome {chromosome}.")
    return filtered_vcf, tsv_output


def chunk_reference_genome(chromosome):
    vcf_path = get_vcf_path(chromosome)
    chunks_output = chunks_path(chromosome)

    print(chunks_output)

    if os.path.exists(chunks_output):
        print(f"Chunk file already exists: {chunks_output}. Skipping chunking.")
        return chunks_output

    print(f"Chunking reference genome for chromosome {chromosome}...")

    subprocess.run([GLIMPSE_CHUNK, 
        "--input", vcf_path, "--region", chromosome,
        "--window-mb", "2", "--buffer-mb", "0.2",
        "--output", chunks_output, "--sequential"
    ], capture_output=True, text=True, check=True)


    print(f"Chunk file created at {chunks_output}.")
    return chunks_output

def prepare_gatk_bundle():
    dbsnp = dbsnp_dir()

    if not os.path.exists(dbsnp):
        print(f"{dbsnp} not found. Compressing {dbsnp}...")
        
        subprocess.run([TOOLS['bgzip'], dbsnp], check=True)
        subprocess.run([TOOLS['tabix'], "-f", f"{dbsnp}.gz"], check=True)
        
        print(f"{dbsnp} has been compressed and indexed.")
    else:
        print(f"{dbsnp} already exists. No action needed.")


def prepare_reference_panel(chromosome):
    print(f"Preparing {chromosome}")
    os.makedirs(reference_path, exist_ok=True)

    if check_reference_panel(chromosome):
        print(f"Reference panel for {chromosome} already exists. Skipping.")
        return
    
    # Step 1: Download reference panel
    #download_reference_panel(chromosome)

    # Step 2: Normalize and filter reference panel
    normalize_and_filter_reference(chromosome)

    # Step 3: Process SNP sites
    process_snp_sites(chromosome)

    # Step 4: Chunk reference genome
    chunk_reference_genome(chromosome)

    print(f"Reference panel preparation completed for {chromosome}.")


def run_prepare_reference_panel():
    # Step 0: verify gatk bundle
    print("Preparing reference....")
    #prepare_gatk_bundle()

    prepare_reference_panel("chr22")
    return

    with ThreadPoolExecutor(max_workers=4) as executor:  
        executor.map(prepare_reference_panel, PARAMETERS["chrs"])
