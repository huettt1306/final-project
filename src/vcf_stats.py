import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from helper.config import PARAMETERS, TRIO_DATA, PATHS

# ==== Cấu hình ====
vcf_dir = PATHS["vcf_directory"]
output_dir = os.path.join(vcf_dir, "stats")
os.makedirs(output_dir, exist_ok=True)

chromosomes = [str(i) for i in range(1, 22)] + ["X"]
summary_data = []

print("📊 Bắt đầu thống kê toàn bộ thông số biến thể...")

for chrom in chromosomes:
    vcf_file = os.path.join(vcf_dir, f"chr{chrom}_variants.vcf.gz")

    # === 1. Tổng số biến thể (SNPs + InDels) ===
    stats_file = os.path.join(output_dir, f"chr{chrom}.stats.txt")
    subprocess.run(f"bcftools stats {vcf_file} > {stats_file}", shell=True, check=True)

    total_variants = 0
    with open(stats_file) as f:
        for line in f:
            if line.startswith("SN") and "number of records:" in line:
                total_variants = int(line.strip().split()[-1])
                break

    # === 2. Số SNPs ===
    snps_file = os.path.join(output_dir, f"chr{chrom}.snps.vcf.gz")
    subprocess.run(f"bcftools view -v snps {vcf_file} -Oz -o {snps_file}", shell=True)
    subprocess.run(f"bcftools index {snps_file}", shell=True)
    snps_stats = subprocess.check_output(f"bcftools index -n {snps_file}", shell=True)
    num_snps = int(snps_stats.decode().strip())

    # === 3. Số InDels ===
    indels_file = os.path.join(output_dir, f"chr{chrom}.indels.vcf.gz")
    subprocess.run(f"bcftools view -v indels {vcf_file} -Oz -o {indels_file}", shell=True)
    subprocess.run(f"bcftools index {indels_file}", shell=True)
    indels_stats = subprocess.check_output(f"bcftools index -n {indels_file}", shell=True)
    num_indels = int(indels_stats.decode().strip())

    summary_data.append({
        "Chromosome": chrom,
        "TotalVariants": total_variants,
        "SNPs": num_snps,
        "InDels": num_indels
    })

# === Lưu bảng thống kê biến thể ===
df = pd.DataFrame(summary_data)
df["Chromosome"] = pd.Categorical(df["Chromosome"], categories=chromosomes, ordered=True)
df.sort_values("Chromosome", inplace=True)
df.to_csv(os.path.join(output_dir, "variant_summary.csv"), index=False)

# === Vẽ biểu đồ số biến thể ===
plt.figure(figsize=(12, 6))
plt.bar(df["Chromosome"], df["SNPs"], label="SNPs", alpha=0.7)
plt.bar(df["Chromosome"], df["InDels"], bottom=df["SNPs"], label="InDels", alpha=0.7)
plt.xlabel("Nhiễm sắc thể")
plt.ylabel("Số biến thể")
plt.title("Phân bố SNPs và InDels theo nhiễm sắc thể")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "variant_distribution.png"))
plt.show()

print("✅ Đã hoàn tất phần 2 - Thống kê biến thể!")

# === PHẦN 3: Tính đa dạng kiểu gen (toàn bộ dữ liệu) ===
print("📊 Bắt đầu phân tích đa dạng kiểu gen (Heterozygosity, MAF)...")

# Tính heterozygosity cho toàn bộ dữ liệu (không cần mẫu riêng biệt)
het_results = []
for chrom in chromosomes:
    vcf_file = os.path.join(vcf_dir, f"chr{chrom}_variants.vcf.gz")
    out_prefix = os.path.join(output_dir, f"chr{chrom}")
    het_file = out_prefix + ".het"

    subprocess.run(f"vcftools --gzvcf {vcf_file} --het --out {out_prefix}", shell=True)
    if os.path.exists(het_file):
        het_df = pd.read_csv(het_file, delim_whitespace=True)
        het_df["CHROM"] = chrom
        het_results.append(het_df)

# Gộp dữ liệu heterozygosity
if het_results:
    het_df_all = pd.concat(het_results, ignore_index=True)
    het_df_all.to_csv(os.path.join(output_dir, "heterozygosity_all_chromosomes.csv"), index=False)

# Tính tần suất alen (MAF) cho toàn bộ dữ liệu
maf_all = []
for chrom in chromosomes:
    vcf_file = os.path.join(vcf_dir, f"chr{chrom}_variants.vcf.gz")
    out_prefix = os.path.join(output_dir, f"chr{chrom}_maf")
    subprocess.run(f"vcftools --gzvcf {vcf_file} --freq2 --out {out_prefix}", shell=True)
    maf_file = out_prefix + ".frq"
    if os.path.exists(maf_file):
        maf_df = pd.read_csv(maf_file, delim_whitespace=True, comment="#")
        maf_df["CHROM"] = chrom
        maf_all.append(maf_df)

if maf_all:
    maf_df_all = pd.concat(maf_all, ignore_index=True)
    maf_df_all.to_csv(os.path.join(output_dir, "maf_all_chromosomes.csv"), index=False)

print("✅ Đã hoàn tất phần 3 - Đa dạng kiểu gen!")

# === Gợi ý bước tiếp theo ===
print("🎉 Tất cả thống kê đã được lưu trong thư mục:", output_dir)
print("📁 Bao gồm: Tổng biến thể, SNP/InDel, dị hợp tử, MAF và biểu đồ.")
