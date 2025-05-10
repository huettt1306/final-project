import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF
from helper.config import PARAMETERS, TRIO_DATA, PATHS


vcf_dir = PATHS["vcf_directory"]
chromosomes = [str(i) for i in range(1, 23)] + ["X"]

# === Hàm so sánh kiểu gen mẹ và con ===
def compare_genotypes(mom_gt, child_gt):
    if mom_gt is None and child_gt is None:
        return "missing"

    mom_gt = tuple(sorted(mom_gt))
    child_gt = tuple(sorted(child_gt))

    if mom_gt == child_gt:
        return "same"
    elif len(set(mom_gt) & set(child_gt)) == 1:
        return "partial"
    else:
        return "different"


# === Tính thống kê cho từng trio ===
trio_stats = {}

for trio_name, trio_info in TRIO_DATA.items():
    mom_id = trio_info["mother"]
    child_id = trio_info["child"]
    same_genotype = 0
    partial_shared = 0
    different = 0

    for chrom in chromosomes:
        vcf_path = os.path.join(vcf_dir, f"chr{chrom}_variants.vcf.gz")
        if not os.path.exists(vcf_path):
            continue

        vcf = VCF(vcf_path)
        if mom_id not in vcf.samples or child_id not in vcf.samples:
            continue

        mom_idx = vcf.samples.index(mom_id)
        child_idx = vcf.samples.index(child_id)

        for variant in vcf:
            mom_gt = variant.genotypes[mom_idx][:2]
            child_gt = variant.genotypes[child_idx][:2]

            result = compare_genotypes(mom_gt, child_gt)
            if result == "same":
                same_genotype += 1
            elif result == "partial":
                partial_shared += 1
            elif result == "different":
                different += 1

    total = same_genotype + partial_shared + different
    trio_stats[trio_name] = {
        "same_genotype": same_genotype,
        "shared_allele_only": partial_shared,
        "different": different,
        "total": total
    }

df = pd.DataFrame.from_dict(trio_stats, orient="index")
df.index.name = "Trio"
df.reset_index(inplace=True)
df.to_csv("trio_genotype_comparison.csv", index=False)

plt.figure(figsize=(10, 6))

x = range(len(df))
plt.bar(x, df["same_genotype"], label="Giống hoàn toàn", color="skyblue")
plt.bar(x, df["shared_allele_only"], bottom=df["same_genotype"], label="Giống 1 alen", color="orange")
plt.bar(x, df["different"], bottom=df["same_genotype"] + df["shared_allele_only"], label="Khác hoàn toàn", color="lightgray")

plt.xticks(x, df["Trio"])
plt.ylabel("Số biến thể")
plt.xlabel("Bộ trio")
plt.title("So sánh kiểu gen mẹ - con theo từng bộ trio")
plt.legend()
plt.tight_layout()
plt.savefig("trio_genotype_comparison_stacked.png")
plt.show()