from helper.config import PARAMETERS, TRIO_DATA, PATHS
from helper.path_define import statistic_outdir, fastq_single_path, statistic_summary

import os
import pandas as pd
import matplotlib.pyplot as plt

# Dictionary lưu dữ liệu số lượng basevar_GT (true + false) theo coverage và chromosome
coverage_chr_data = {cov: {} for cov in PARAMETERS["coverage"]}

def summary_df(output_dir):
    for coverage in PARAMETERS["coverage"]:
        print(f"\U0001F4CA Processing coverage = {coverage}")

        for chr_ in PARAMETERS["chrs"]:
            total = 0
            for trio_name, trio_info in TRIO_DATA.items():
                sample = trio_info["mother"] 
                fq_path = os.path.join(fastq_single_path(sample, coverage), f"{sample}.fastq.gz")
                df_path = statistic_summary(fq_path, chr_)

                if not os.path.exists(df_path):
                    continue

                try:
                    df = pd.read_csv(df_path)
                    total += (df["basevar_GT_true"].sum() + df["basevar_GT_false"].sum()) * 1.8
                except Exception as e:
                    print(f"Error reading {df_path}: {e}")

            coverage_chr_data[coverage][chr_] = total


def plot_all(output_dir):
    plt.figure(figsize=(12, 6))

    markers = ['o', 's', 'D', '^']
    for i, (coverage, chr_totals) in enumerate(coverage_chr_data.items()):
        chrs = list(chr_totals.keys())
        values = [chr_totals[chr_] for chr_ in chrs]
        plt.plot(chrs, values, label=f"{coverage}x", marker=markers[i % len(markers)], linewidth=2)

    plt.grid(True, which='both', axis='y', linestyle='--', linewidth=0.5)

    plt.xlabel("NST")
    plt.ylabel("Số biến thể")
    plt.title("Số biến thể theo NST và Coverage")
    plt.legend(title="Coverage")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "basevar_GT_by_chr_and_coverage.pdf"))
    plt.close()


def main():
    output_dir = os.path.join(PATHS["plot_directory"])
    os.makedirs(output_dir, exist_ok=True)

    summary_df(output_dir)
    plot_all(output_dir)


if __name__ == "__main__":
    main()
