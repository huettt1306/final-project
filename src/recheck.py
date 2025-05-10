from helper.config import PARAMETERS, TRIO_DATA, PATHS, TOOLS
from helper.logger import setup_logger
from helper.path_define import fastq_single_path, bam_dir, bamlist_dir
import os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor
from helper.file_utils import create_vcf_list, merge_vcf_list
from steps.basevar import run_basevar




logger = setup_logger(os.path.join(PATHS["logs"], "recheck.log"))
REF = PATHS["ref"]
REF_FAI = PATHS["ref_fai"]
BCFTOOLS = TOOLS["bcftools"]
TABIX = TOOLS["tabix"]
DELTA = PARAMETERS["basevar"]["delta"]



def main():
    """
    Xử lý một trio (bao gồm các bước pipeline cho từng mẫu).
    """
    logger.info(f"######## PROCESSING RECHECK ########")

    
    for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
        for coverage in PARAMETERS["coverage"]:
            logger.info(f"######## PROCESSING index {index} and coverage {coverage} ########")

            alldir = os.path.join(PATHS["result_directory"], f"{coverage}x", "all", f"sample_{index}")
            os.makedirs(alldir, exist_ok=True)
            fqdir = os.path.join(alldir, "all.fastq.gz")
            bam_list = bamlist_dir(fqdir)
            os.makedirs(os.path.dirname(bam_list), exist_ok=True)
            
            with open(bam_list, "w") as f_out:
                for trio_name, trio_info in TRIO_DATA.items():
                    mother_name = trio_info["mother"]

                    fastq_path = os.path.join(fastq_single_path(mother_name, coverage, index), f"{mother_name}.fastq.gz")
                    
                    f_out.write(bam_dir(fastq_path) + "\n")
                
            run_basevar(fqdir)

            

            

if __name__ == "__main__":
    main()
