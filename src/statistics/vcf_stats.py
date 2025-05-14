import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

vcf_dir = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/vcf_"
output_dir = os.path.join(vcf_dir, "stats")
os.makedirs(output_dir, exist_ok=True)

chromosomes = [str(i) for i in range(1, 23)] + ["X"]
summary_data = []

print("Bắt đầu thống kê toàn bộ thông số biến thể...")

def stats_all():
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

    df = pd.DataFrame(summary_data)
    df["Chromosome"] = pd.Categorical(df["Chromosome"], categories=chromosomes, ordered=True)
    df.sort_values("Chromosome", inplace=True)
    df.to_csv(os.path.join(output_dir, "variant_summary.csv"), index=False)


def plot_all():
    df = pd.read_csv(os.path.join(output_dir, "variant_summary.csv"))

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

    print("Thống kê biến thể hoàn thành!")

plot_all()