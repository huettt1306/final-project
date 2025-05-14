import os
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import numpy as np

vcf_dir = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/vcf_"
output_dir = os.path.join(vcf_dir, "maf_analysis")
os.makedirs(output_dir, exist_ok=True)
chromosomes = [str(i) for i in range(1, 22)] + ["X"]

maf_output_file = os.path.join(output_dir, "all_maf.txt")

print("🔬 Đang tính MAF bằng bcftools...")
with open(maf_output_file, "w") as outfile:
    for chrom in chromosomes:
        vcf_path = os.path.join(vcf_dir, f"chr{chrom}_variants.vcf.gz")
        if not os.path.exists(vcf_path):
            continue
        cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF\\n' {vcf_path}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        outfile.write(result.stdout)

print("✅ Đã trích xuất tần suất alen (AF). Đang vẽ phân bố...")

maf_df = pd.read_csv(maf_output_file, sep="\t", header=None,
                     names=["CHROM", "POS", "REF", "ALT", "AF"])

maf_df = maf_df[pd.to_numeric(maf_df["AF"], errors='coerce').notnull()]
maf_df["AF"] = maf_df["AF"].astype(float)
maf_df = maf_df[(maf_df["AF"] >= 0.001) & (maf_df["AF"] <= 0.5)]  # Giữ MAF hợp lệ

counts, bin_edges = np.histogram(maf_df["AF"], bins=1000, range=(0.001, 0.5))
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

plt.figure(figsize=(10, 5))
plt.plot(bin_centers, counts, color="mediumseagreen", linewidth=1.2)
plt.xlabel("Tần suất alen nhỏ (MAF)")
plt.ylabel("Số biến thể")
plt.title("Phân bố MAF trong toàn bộ biến thể")

plt.xticks(ticks=[i/10 for i in range(0, 6)], labels=[f"{i/10:.1f}" for i in range(0, 6)])
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "maf_distribution_lineplot.png"))
plt.show()
