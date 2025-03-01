import os
from helper.config import PATHS
from helper.file_utils import extract_vcf

def cram_path(name):
    return os.path.join(PATHS["cram_directory"], f"{name}.final.cram")

def fastq_path(name):
    return os.path.join(PATHS["fastq_directory"], f"{name}.fastq.gz")

def fastq_path_lane1(name):
    return os.path.join(PATHS["fastq_directory"], f"{name}_1.fastq.gz")

def fastq_path_lane2(name):
    return os.path.join(PATHS["fastq_directory"], f"{name}_2.fastq.gz")

def fastq_single_path(name, coverage, index):
    return os.path.join(PATHS["result_directory"], f"{coverage}x", name, f"sample_{index}")

def fastq_nipt_path(child, mother, father, coverage, ff, index):
    return os.path.join(PATHS["result_directory"], f"{coverage}x", f"{child}_{mother}_{father}", f"{ff:.2f}", f"sample_{index}")

def base_dir(fq):
    return os.path.dirname(fq)

def samid(fq):
    return os.path.basename(fq).replace(".fastq.gz", "")

def tmp_outdir(fq):
    return os.path.join(base_dir(fq), "1tmp_files")

def batch1_final_outdir(fq):
    return os.path.join(base_dir(fq), "batch1_final_files")

def bamlist_dir(fq):
    return os.path.join(batch1_final_outdir(fq), "bam.list")

def basevar_outdir(fq):
    return os.path.join(base_dir(fq), "basevar_output")

def basevar_vcf(fq, chromosome):
    return os.path.join(f"{basevar_outdir(fq)}_final", f"{samid(fq)}_basevar_{chromosome}.vcf.gz")

def vcf_list_path(fq, chromosome):
    return os.path.join(basevar_outdir(fq), f"{samid(fq)}_basevar_{chromosome}.vcf.list")

def glimpse_outdir(fq):
    return os.path.join(base_dir(fq), "glimpse_output")

def vcf_prefix(chromosome):
    if(chromosome == "chrX"):
        return "CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2"
    return f"CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.filtered.shapeit2-duohmm-phased"

def get_vcf_path(chromosome):
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.vcf.gz")

def get_tsv_path(chromosome):
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.vcf.gz")

def norm_vcf_path(chromosome): 
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.biallelic.snp.maf0.001.vcf.gz")

def filtered_vcf_path(chromosome): 
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.biallelic.snp.maf0.001.sites.vcf.gz")

def filtered_tsv_path(chromosome):
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.biallelic.snp.maf0.001.sites.tsv.gz")

def chunks_path(chromosome):
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.chunks.txt")

def glimpse_vcf(fq, chromosome):
    return os.path.join(glimpse_outdir(fq), "imputed", f"{samid(fq)}_glimpse.{chromosome}_imputed.vcf.gz")

def glimpse_annot(fq, chromosome):
    return os.path.join(glimpse_outdir(fq), "annotated", f"{samid(fq)}_glimpse.{chromosome}_annotated.vcf.gz")

def get_vcf_ref(chromosome):
    return os.path.join(PATHS["vcf_directory"], f"20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.recalibrated_variants.vcf.gz")

def ground_truth_vcf(name, chromosome):
    path = os.path.join(PATHS["vcf_directory"], f"{name}_{chromosome}.vcf.gz")
    if os.path.exists(path):
        print(f"Ground truth VCF already exists: {path}")
        return path
    
    extract_vcf(name, get_vcf_ref(chromosome), path)
    return path

def statistic_outdir(fq, chromosome="all"):
    if chromosome == "all":
        return os.path.join(base_dir(fq), "statistic_output")
    return os.path.join(base_dir(fq), "statistic_output", f"{chromosome}")

def statistic_variants(fq, chromosome="all"):
    return os.path.join(statistic_outdir(fq, chromosome), f"{samid(fq)}_variants.csv")

def statistic_summary(fq, chromosome="all"):
    return os.path.join(statistic_outdir(fq, chromosome), f"{samid(fq)}_summary.csv")
    
def statistic_rare_summary(fq, chromosome="all"):
    return os.path.join(statistic_outdir(fq, chromosome), "rare_summary.csv")
    
def dbsnp_dir():
    return os.path.join(PATHS["gatk_bundle_dir"], "Homo_sapiens_assembly38.dbsnp138.vcf.gz")