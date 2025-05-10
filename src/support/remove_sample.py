import subprocess
from concurrent.futures import ThreadPoolExecutor

# Đường dẫn tới file chứa danh sách sample cần loại bỏ
REMOVE_SAMPLE_FILE = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/phased/remove_samples.txt"

# Hàm trả về đường dẫn tới file VCF tương ứng với nhiễm sắc thể
def get_path(chromosome: str) -> str:
    if(chromosome == "chrX"):
        return "/home/huettt/Documents/nipt/NIPT-human-genetics/working/phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.diploid_recal.vcf.gz"
    return f"/home/huettt/Documents/nipt/NIPT-human-genetics/working/phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz"

def get_new_path(chromosome: str) -> str:
    if(chromosome == "chrX"):
        return "/home/huettt/Documents/nipt/NIPT-human-genetics/working/phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.diploid_recal-rmsp.vcf.gz"
    return f"/home/huettt/Documents/nipt/NIPT-human-genetics/working/phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.filtered.shapeit2-duohmm-phased-rmsp.vcf.gz"

# Hàm xử lý từng nhiễm sắc thể: lọc sample và index lại
def process_chr(chr_name: str):
    input_vcf = get_path(chr_name)
    output_vcf = get_new_path(chr_name)

    print(f"[chr{chr_name}] Filtering...")

    try:
        # Lọc các sample
        subprocess.run([
            "bcftools", "view",
            "--force-samples",
            "-S", f"^{REMOVE_SAMPLE_FILE}",
            "-Oz",
            "-o", output_vcf,
            input_vcf
        ], check=True)

        print(f"[chr{chr_name}] Indexing...")

        # Tạo chỉ mục
        subprocess.run(["bcftools", "index", output_vcf], check=True)

        print(f"[chr{chr_name}] Done.")

    except subprocess.CalledProcessError as e:
        print(f"[chr{chr_name}] Error: {e}")

# Danh sách các nhiễm sắc thể
chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# Chạy song song với ThreadPoolExecutor
if __name__ == "__main__":
    process_chr("chr22")
