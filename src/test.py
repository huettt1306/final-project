from pipeline.generate import generate_single_sample, generate_nipt_sample
from pipeline.alignment import run_alignment_pipeline
from pipeline.basevar import run_basevar
from pipeline.glimpse import run_glimpse
from statistic.statistic import run_statistic

from pipeline.reference_panel_prepare import run_prepare_reference_panel
from helper.config import PARAMETERS, TRIO_DATA, PATHS, TOOLS
from helper.metrics import get_fastq_coverage
from helper.logger import setup_logger
from helper.file_utils import extract_lane1_fq, extract_vcf
from helper.converter import convert_cram_to_fastq
from helper.path_define import ground_truth_vcf, fastq_path_lane1, fastq_path_lane2, cram_path, get_vcf_ref
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor


logger = setup_logger(os.path.join(PATHS["logs"], "main.log"))


def pipeline_for_sample(fastq_dir):
    logger.info(f"Run all pipeline for sample in {fastq_dir}")
    run_alignment_pipeline(fastq_dir)
    run_basevar(fastq_dir)
    run_glimpse(fastq_dir)

def prepare_data(name):
    print(f"Preparing data for {name}")
    if not os.path.exists(fastq_path_lane1(name)):
        convert_cram_to_fastq(cram_path(name), fastq_path_lane1(name), fastq_path_lane2(name))
    #return get_fastq_coverage(name)

def process_trio(trio_name, trio_info):
    """
    Xử lý một trio (bao gồm các bước pipeline cho từng mẫu).
    """
    logger.info(f"######## PROCESSING TRIO: {trio_name} ########")

    child_name = trio_info["child"]
    mother_name = trio_info["mother"]
    father_name = trio_info["father"]

    
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.map(prepare_data, [child_name, mother_name])

        #mother_avg_coverage = future_mother.result()
        #child_avg_coverage = future_child.result()


    for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
        logger.info(f"######## PROCESSING index {index} ########")

        for coverage in PARAMETERS["coverage"]:
            #pipeline_for_sample(generate_single_sample(mother_name, coverage, index))

            for ff in PARAMETERS["ff"]:
                pipeline_for_sample(generate_nipt_sample(child_name, mother_name, father_name, coverage, ff, index))


def main():    
    run_prepare_reference_panel()
    if len(sys.argv) < 2:
        logger.error("Please provide a sample name to process.")
        sys.exit(1)
    
    id = sys.argv[1]
    sample = f"/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/sample{id}/PGT-sample{id}.fastq.gz"
    
    
    pipeline_for_sample(sample)

def main1():    
    trios = ["VN047", "VN051", "VN066", "VN078", "VN91"]
    chromosomes = ["chr20"]
    
    for chromosome in chromosomes:
        vcf_reference = get_vcf_ref(chromosome)

        if not os.path.exists(vcf_reference):
            url = f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.recalibrated_variants.vcf.gz"
            subprocess.run(["wget", "-c", url, "-P", PATHS["vcf_directory"]], capture_output=True, text=True)
            

        subprocess.run([TOOLS["tabix"], "-f", vcf_reference], capture_output=True, text=True, check=True)
        for trio_name in trios:
            print(f"Processing {trio_name}...")
            trio_info = TRIO_DATA[trio_name]
            
            child_name = trio_info["child"]
            mother_name = trio_info["mother"]
            father_name = trio_info["father"]

            with ThreadPoolExecutor(max_workers=3) as executor:
                executor.map(extract_vcf, [child_name, mother_name, father_name], [chromosome, chromosome, chromosome])

        #os.remove(vcf_reference)

if __name__ == "__main__":
    main1()
