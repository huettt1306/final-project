from helper.config import PARAMETERS, TRIO_DATA, PATHS
from helper.path_define import fastq_single_path, statistic_summary

import os
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

def sort_nst(nst):
    if nst == "chrX":
        return (23, "X")
    elif nst == "chrY":
        return (24, "Y")
    elif nst == "all":
        return (26, "all")
    else:
        num = int(nst.replace("chr", ""))
        return (num, "")

def load_ground_truth_stats(path):
    df = pd.read_csv(path, sep=",")
    gt_stats = defaultdict(lambda: defaultdict(dict))
    for _, row in df.iterrows():
        chr_, sample, maf = str(row["Chr"]), row["Sample"], int(row["MAF"])
        gt_stats[chr_][sample][maf] = {
            "total_gt": row["Total_GT"],
            "alt_gt": row["ALT_GT"],
        }
    return gt_stats


def recal_stats(df):
    df['call_rate_gt'] = (df['GT_true'] + df['GT_false']) / df['total_GT']
    df['accuracy_gt'] = df['GT_true'] / (df['GT_true'] + df['GT_false'])

    df['call_rate_alt'] = (df['ALT_true'] + df['ALT_false']) / df['total_ALT']
    df['accuracy_alt'] = df['ALT_true'] / (df['ALT_true'] + df['GT_false'])

    return df


def compute_bin_stats(df, gt_stats, chr_, sample, method):
    rows = []
    for maf in range(0, 51):
        gt_true = df.at[maf, f"{method}_GT_true"]
        gt_false = df.at[maf, f"{method}_GT_false"]
        alt_true = df.at[maf, f"{method}_ALT_true"]
        alt_false = df.at[maf, f"{method}_ALT_false"]

        gt_info = gt_stats.get(chr_, {}).get(sample, {}).get(maf, None)
        if gt_info is None:
            continue

        total_gt = gt_info["total_gt"]
        alt_gt = gt_info["alt_gt"]

        rows.append({
            "NST": chr_,
            "MAF": maf,
            "GT_true": gt_true,
            "GT_false": gt_false,
            "ALT_true": alt_true,
            "ALT_false": alt_false,
            "total_GT": total_gt,
            "total_ALT": alt_gt,
        })

    return pd.DataFrame(rows)


def compute_cumsum_stats(df, gt_stats, chr_, sample, method):
    rows = []
    for threshold in range(0, 51):
        gt_true = df.loc[threshold:50, f"{method}_GT_true"].sum()
        gt_false = df.loc[threshold:50, f"{method}_GT_false"].sum()
        alt_true = df.loc[threshold:50, f"{method}_ALT_true"].sum()
        alt_false = df.loc[threshold:50, f"{method}_ALT_false"].sum()

        total_gt = 0
        alt_gt = 0
        for maf in range(threshold, 51):
            stats = gt_stats.get(chr_, {}).get(sample, {}).get(maf, None)
            if stats:
                total_gt += stats["total_gt"]
                alt_gt += stats["alt_gt"]

        rows.append({
            "NST": chr_,
            "MAF_≥": threshold,
            "GT_true": gt_true,
            "GT_false": gt_false,
            "ALT_true": alt_true,
            "ALT_false": alt_false,
            "total_GT": total_gt,
            "total_ALT": alt_gt
        })

    return pd.DataFrame(rows)


def recheck_basevar_merge(output_dir=os.path.join(PATHS["plot_directory"]), coverage=0.1, index=1):
    gt_stats = load_ground_truth_stats(PATHS["ground_truth_stat"])
    fqdir = os.path.join(PATHS["result_directory"], f"{coverage}x", "all", f"sample_{index}", "all.fastq.gz")

    all_dfs = []

    for chr_id in PARAMETERS["chrs"]:
        file_path = statistic_summary(fqdir, chr_id)
        df = pd.read_csv(file_path)
        df["NST"] = chr_id

        def compute_total_alt(sample):
            maf_dict = gt_stats.get(chr_id, {}).get(sample, {})
            return sum(m["alt_gt"] for m in maf_dict.values())

        df["total_alt"] = df["Sample"].apply(compute_total_alt)
        df["call_rate"] = (df["ALT_correct"] + df["ALT_wrong"]) / df["total_alt"]
        df["accuracy"] = df["ALT_correct"] / (df["ALT_correct"] + df["ALT_wrong"])

        all_dfs.append(df)

    merged_df = pd.concat(all_dfs, ignore_index=True)

    mean_by_chr = merged_df.groupby("NST")[["call_rate", "accuracy"]].mean().reset_index()

    output_file = os.path.join(output_dir, f"recheck_basevar_{coverage}.csv")
    merged_df.to_csv(output_file, index=True)

    mean_file = os.path.join(output_dir, f"recheck_basevar_{coverage}_mean_by_chr.csv")
    mean_by_chr.to_csv(mean_file, index=False)

    return merged_df, mean_by_chr




def process_summary(df_path, gt_stats, chr_, sample):
    df = pd.read_csv(df_path)
    df["MAF"] = df.index 
    df["NST"] = chr_      

    results_bin = {}
    results_maf = {}

    for method in ["glimpse"]:
        bin_df = compute_bin_stats(df, gt_stats, chr_, sample, method)
        maf_df = compute_cumsum_stats(df, gt_stats, chr_, sample, method)
        results_bin[method] = bin_df
        results_maf[method] = maf_df

    return results_bin, results_maf


def plot_stats_multiple(df_by_coverage, method, type, Ox):
    stats = ['call_rate_gt', 'accuracy_gt', 'call_rate_alt', 'accuracy_alt']
    statnames = [
        "Tỷ lệ gọi kiểu gen",
        "Độ chính xác của kiểu gen",
        "Tỷ lệ gọi kiểu gen có allen thay thế ",
        "Độ chính xác của allen thay thế"
    ]
    xlabel = 'MAF tối thiểu' if type == "maf" else 'MAF'
    
    coverages = list(df_by_coverage.keys())
    
    for i, stat in enumerate(stats):
        fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharey=True)
        axes = axes.flatten()

        for idx, coverage in enumerate(coverages):
            df_bin = df_by_coverage[coverage]
            ax = axes[idx]

            for nst in df_bin['NST'].unique():
                nst_data = df_bin[df_bin['NST'] == nst]
                maf_values = nst_data[Ox]
                accuracy_values = nst_data[stat]
                
                if nst != 'all':
                    ax.plot(maf_values, accuracy_values, label=nst, linestyle='-', alpha=0.6, linewidth=1)
                else:
                    ax.plot(maf_values, accuracy_values, label=nst, color="blue", linestyle='-', alpha=0.9, linewidth=3)

            ax.set_title(f'Coverage {coverage}x', fontsize=14)
            ax.set_xlabel(xlabel, fontsize=12)
            if idx % 2 == 0:
                ax.set_ylabel(statnames[i], fontsize=12)
            ax.grid(True)

        handles, labels = ax.get_legend_handles_labels()

        plt.suptitle(f'{statnames[i]} theo {xlabel}', fontsize=18)

        plt.tight_layout(rect=[0, 0.1, 1, 0.95])

        fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=10)

        plt.savefig(os.path.join(PATHS["plot_directory"], "png", "single", f"{stat}_{method}_{type}_multiple.png"), dpi=180)
        plt.close()



def plot_gt_stacked(data, maf_value, method="glimpse"):
    coverages = []
    gt_trues = []
    gt_falses = []
    gt_not_calleds = []

    alt_trues = []
    alt_falses = []
    alt_not_calleds = []

    for item in data:
        coverage = item["coverage"]
        df = item["avg_maf"]

        row = df[df["MAF_≥"] == maf_value]
        if row.empty:
            continue

        row = row.iloc[0]
        gt_true = row["GT_true"]
        gt_false = row["GT_false"]
        gt_not_called = row["total_GT"] - gt_true - gt_false

        alt_true = row["ALT_true"]
        alt_false = row["ALT_false"]
        alt_not_called = row["total_ALT"] - alt_true - alt_false

        coverages.append(str(coverage))
        gt_trues.append(gt_true)
        gt_falses.append(gt_false)
        gt_not_calleds.append(gt_not_called)

        alt_trues.append(alt_true)
        alt_falses.append(alt_false)
        alt_not_calleds.append(alt_not_called)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3), sharey=True, constrained_layout=True)

    colors = ["#4daf4a", "#e41a1c", "#ff7f00"]
    labels = ["True", "False", "Not Called"]

    gt_bottom1 = gt_trues
    gt_bottom2 = [gt_trues[i] + gt_falses[i] for i in range(len(gt_trues))]

    ax1.bar(coverages, gt_trues, color=colors[0], width=0.4)
    ax1.bar(coverages, gt_falses, bottom=gt_trues, color=colors[1], width=0.4)
    ax1.bar(coverages, gt_not_calleds, bottom=gt_bottom2, color=colors[2], width=0.4)
    ax1.set_ylabel('Số lượng GT')
    ax1.set_title(f'Kết quả gọi kiểu gen (GT)\nMAF ≥ {maf_value}')

    alt_bottom1 = alt_trues
    alt_bottom2 = [alt_trues[i] + alt_falses[i] for i in range(len(alt_trues))]

    ax2.bar(coverages, alt_trues, color=colors[0], width=0.4)
    ax2.bar(coverages, alt_falses, bottom=alt_trues, color=colors[1], width=0.4)
    ax2.bar(coverages, alt_not_calleds, bottom=alt_bottom2, color=colors[2], width=0.4)
    ax2.set_ylabel('Số lượng ALT')
    ax2.set_title(f'Kết quả gọi alen thay thế (ALT)\nMAF ≥ {maf_value}')
    ax2.set_xlabel('Coverage')

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors[0]),
        plt.Rectangle((0, 0), 1, 1, color=colors[1]),
        plt.Rectangle((0, 0), 1, 1, color=colors[2])
    ]
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.01))

    plt.tight_layout(rect=[0, 0.05, 1, 1]) 
    plt.savefig(
        os.path.join(PATHS["plot_directory"], "png", f"gt_alt_{method}_{maf_value}.png"),
        dpi=180
    )
    plt.close()


def plot_basevar_recheck(one_sample_path, five_sample_path):
    five_sample = pd.read_csv(five_sample_path) 
    one_sample = pd.read_csv(one_sample_path)
    one_sample = one_sample[one_sample["MAF_≥"] == 0]
    one_sample = one_sample[one_sample["NST"] != "all"]

    five_sample = five_sample.iloc[five_sample["NST"].apply(sort_nst).argsort()]
    one_sample = one_sample.iloc[one_sample["NST"].apply(sort_nst).argsort()]

    plt.figure(figsize=(12, 5))
    plt.plot(one_sample["NST"], one_sample["call_rate_alt"], label="one_sample", marker='o')
    plt.plot(five_sample["NST"], five_sample["call_rate"], label="five_sample", marker='s')
    plt.xlabel("NST")
    plt.ylabel("Call Rate")
    plt.title("Tỉ lệ gọi biến thể trên các NST khi thực hiện basevar trên một mẫu và trên nhiều mẫu")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        os.path.join(PATHS["plot_directory"], "png", f"basevar_recheck_callrate.png"),
        dpi=180
    )

    plt.figure(figsize=(12, 5))
    plt.plot(one_sample["NST"], one_sample["accuracy_alt"], label="one_sample", marker='o')
    plt.plot(five_sample["NST"], five_sample["accuracy"], label="five_sample", marker='s')
    plt.xlabel("NST")
    plt.ylabel("Accuracy")
    plt.title("Độ chính xác của biến thể trên các NST khi thực hiện basevar trên một mẫu và trên nhiều mẫu")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        os.path.join(PATHS["plot_directory"], "png", f"basevar_recheck_accuracy.png"),
        dpi=180
    )


def summary_df(output_dir):
    gt_stats = load_ground_truth_stats(PATHS["ground_truth_stat"])  
    for coverage in PARAMETERS["coverage"]:
        print(f"Processing coverage = {coverage}")
        
        all_stats = {"bin": [], "maf": []}

        for trio_name, trio_info in TRIO_DATA.items():
            sample = trio_info["mother"] 
            fq_path = os.path.join(fastq_single_path(sample, coverage), f"{sample}.fastq.gz")

            df_bin = []
            df_maf = []

            for chr_ in PARAMETERS["chrs"]:
                df_path = statistic_summary(fq_path, chr_)
                if not os.path.exists(df_path):
                    continue

                stats_bin, stats_maf = process_summary(df_path, gt_stats, chr_, sample)
                
                df_bin.append(stats_bin["glimpse"])
                df_maf.append(stats_maf["glimpse"])

            df_bin = pd.concat(df_bin, ignore_index=True)
            numeric_cols = df_bin.select_dtypes(include='number').columns.drop("MAF", errors="ignore")
            bin_all = df_bin.groupby("MAF")[numeric_cols].sum().reset_index()
            bin_all["NST"] = "all"
            df_bin = pd.concat([df_bin, bin_all], ignore_index=True)
            df_bin = recal_stats(df_bin)

            df_maf = pd.concat(df_maf, ignore_index=True)
            numeric_cols = df_maf.select_dtypes(include='number').columns.drop("MAF_≥", errors="ignore")
            maf_all = df_maf.groupby("MAF_≥")[numeric_cols].sum().reset_index()
            maf_all["NST"] = "all"
            df_maf = pd.concat([df_maf, maf_all], ignore_index=True)
            df_maf = recal_stats(df_maf)

            all_stats["bin"].append(df_bin)
            all_stats["maf"].append(df_maf)

        df_bin_all = pd.concat(all_stats["bin"], ignore_index=True)
        numeric_cols = df_bin_all.select_dtypes(include='number').columns.drop("MAF", errors="ignore")
        bin_mean = df_bin_all.groupby(["NST", "MAF"])[numeric_cols].mean().reset_index()
        bin_mean.to_csv(os.path.join(output_dir, f"glimpse_bin_coverage_{coverage}.csv"))
        

        df_maf_all = pd.concat(all_stats["maf"], ignore_index=True)
        numeric_cols = df_maf_all.select_dtypes(include='number').columns.drop("MAF_≥", errors="ignore")
        maf_mean = df_maf_all.groupby(["NST", "MAF_≥"])[numeric_cols].mean().reset_index()
        maf_mean.to_csv(os.path.join(output_dir, f"glimpse_maf_coverage_{coverage}.csv"))



def plot_all(output_dir):
    pl_stats = []
    df_by_coverage_bin = {}
    df_by_coverage_maf = {}

    for coverage in PARAMETERS["coverage"]:
        print(f"Processing coverage = {coverage}")

        for method in ["glimpse"]:
            output_file_bin = os.path.join(output_dir, f"{method}_bin_coverage_{coverage}.csv")
            output_file_maf = os.path.join(output_dir, f"{method}_maf_coverage_{coverage}.csv")
            avg_bin = pd.read_csv(output_file_bin)
            avg_maf = pd.read_csv(output_file_maf)

            df_by_coverage_bin[coverage] = avg_bin
            df_by_coverage_maf[coverage] = avg_maf

            print(f"{method} → Loaded data for coverage {coverage}")

            if method == "glimpse":
                pl_stats.append({
                    "coverage": coverage,
                    "avg_maf": avg_maf
                })

    for method in ["glimpse"]:
        plot_stats_multiple(df_by_coverage_bin, method, "bin", "MAF")
        plot_stats_multiple(df_by_coverage_maf, method, "maf", "MAF_≥")

    print("All statistics saved successfully.")
    
    plot_gt_stacked(pl_stats, 0)
    plot_gt_stacked(pl_stats, 1)
    plot_gt_stacked(pl_stats, 5)



def main():
    output_dir = os.path.join(PATHS["plot_directory"])
    os.makedirs(output_dir, exist_ok=True)

    summary_df(output_dir)
    plot_all(output_dir)
    #recheck_basevar_merge()
    #plot_basevar_recheck(os.path.join(output_dir, "basevar_maf_coverage_0.1.csv"), os.path.join(output_dir, "recheck_basevar_0.1_mean_by_chr.csv"))



if __name__ == "__main__":
    main()
