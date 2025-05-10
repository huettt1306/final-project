from helper.config import PARAMETERS, TRIO_DATA, PATHS
from helper.path_define import statistic_outdir, fastq_nipt_path, statistic_summary

import os
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def load_ground_truth_stats(path):
    # Load ground truth file into dictionary with structure: gt_stats[chr][sample][maf] = stats dict
    df = pd.read_csv(path, sep=",")
    gt_stats = defaultdict(lambda: defaultdict(dict))
    for _, row in df.iterrows():
        chr_, sample, maf = str(row["Chr"]), row["Sample"], int(row["MAF"])
        gt_stats[chr_][sample][maf] = {
            "total_gt": row["Total_GT"],
            "alt_gt": row["ALT_GT"],
            "het_gt": row["HET_GT"],
            "hom_alt": row["HOM_ALT"]
        }
    return gt_stats


def compute_bin_stats(df, gt_stats, chr_, mother, child):
    rows = []
    for maf in range(0, 51):
        gt_both_correct = df.at[maf, "gt_both_correct"]
        gt_mom_correct = df.at[maf, "gt_mom_correct"]
        gt_child_correct = df.at[maf, "gt_child_correct"]
        gt_mom_correct_child_wrong = df.at[maf, "gt_mom_correct_child_wrong"]
        gt_child_correct_mom_wrong = df.at[maf, "gt_child_correct_mom_wrong"]
        gt_both_wrong = df.at[maf, "gt_both_wrong"]
        gt_mom_wrong = gt_both_wrong + gt_child_correct_mom_wrong
        gt_child_wrong = gt_both_wrong + gt_mom_correct_child_wrong

        alt_both_correct = df.at[maf, "alt_both_correct"]
        alt_mom_correct = df.at[maf, "alt_mom_correct"]
        alt_child_correct = df.at[maf, "alt_child_correct"]
        alt_mom_correct_child_wrong = df.at[maf, "alt_mom_correct_child_wrong"]
        alt_child_correct_mom_wrong = df.at[maf, "alt_child_correct_mom_wrong"]
        alt_both_wrong = df.at[maf, "alt_both_wrong"]
        alt_mom_wrong = alt_both_wrong + alt_child_correct_mom_wrong
        alt_child_wrong = alt_both_wrong + alt_mom_correct_child_wrong
        

        mother_info = gt_stats.get(chr_, {}).get(mother, {}).get(maf, None)
        if mother_info is None:
            continue
        total_gt_mother = mother_info["total_gt"]
        total_alt_mother = mother_info["alt_gt"]

        child_info = gt_stats.get(chr_, {}).get(child, {}).get(maf, None)
        if child_info is None:
            continue
        total_gt_child = child_info["total_gt"]
        total_alt_child = child_info["alt_gt"]

        rows.append({
            "NST": chr_,
            "MAF": maf,

            "gt_both_correct": gt_both_correct,
            "gt_mom_correct": gt_mom_correct,
            "gt_child_correct": gt_child_correct,
            "gt_mom_correct_child_wrong": gt_mom_correct_child_wrong,
            "gt_child_correct_mom_wrong": gt_child_correct_mom_wrong,
            "gt_both_wrong": gt_both_wrong,
            "gt_mom_wrong": gt_mom_wrong,
            "gt_child_wrong": gt_child_wrong,
            "total_gt_child": total_gt_child,
            "total_gt_mother": total_gt_mother,
            
            "alt_both_correct": alt_both_correct,
            "alt_mom_correct": alt_mom_correct,
            "alt_child_correct": alt_child_correct,
            "alt_mom_correct_child_wrong": alt_mom_correct_child_wrong,
            "alt_child_correct_mom_wrong": alt_child_correct_mom_wrong,
            "alt_both_wrong": alt_both_wrong,
            "alt_mom_wrong": alt_mom_wrong,
            "alt_child_wrong": alt_child_wrong,
            "total_alt_child": total_alt_child,
            "total_alt_mother": total_alt_mother,
        })
    return pd.DataFrame(rows)


def compute_cumsum_stats(df, gt_stats, chr_, mother, child):
    rows = []
    for threshold in range(0, 51):
        gt_both_correct = df.loc[threshold:50, "gt_both_correct"].sum()
        gt_mom_correct = df.loc[threshold:50, "gt_mom_correct"].sum()
        gt_child_correct = df.loc[threshold:50, "gt_child_correct"].sum()
        gt_mom_correct_child_wrong = df.loc[threshold:50, "gt_mom_correct_child_wrong"].sum()
        gt_child_correct_mom_wrong = df.loc[threshold:50, "gt_child_correct_mom_wrong"].sum()
        gt_both_wrong = df.loc[threshold:50, "gt_both_wrong"].sum()
        gt_mom_wrong = gt_both_wrong + gt_child_correct_mom_wrong
        gt_child_wrong = gt_both_wrong + gt_mom_correct_child_wrong

        alt_both_correct = df.loc[threshold:50, "alt_both_correct"].sum()
        alt_mom_correct = df.loc[threshold:50, "alt_mom_correct"].sum()
        alt_child_correct = df.loc[threshold:50, "alt_child_correct"].sum()
        alt_mom_correct_child_wrong = df.loc[threshold:50, "alt_mom_correct_child_wrong"].sum()
        alt_child_correct_mom_wrong = df.loc[threshold:50, "alt_child_correct_mom_wrong"].sum()
        alt_both_wrong = df.loc[threshold:50, "alt_both_wrong"].sum()
        alt_mom_wrong = alt_both_wrong + alt_child_correct_mom_wrong
        alt_child_wrong = alt_both_wrong + alt_mom_correct_child_wrong

        total_gt_mother = 0
        alt_gt_mother = 0
        total_gt_child = 0
        alt_gt_child = 0

        for maf in range(threshold, 51):
            stats_mother = gt_stats.get(chr_, {}).get(mother, {}).get(maf, None)
            if stats_mother:
                total_gt_mother += stats_mother["total_gt"]
                alt_gt_mother += stats_mother["alt_gt"]

            stats_child = gt_stats.get(chr_, {}).get(child, {}).get(maf, None)
            if stats_child:
                total_gt_child += stats_child["total_gt"]
                alt_gt_child += stats_child["alt_gt"]

        rows.append({
            "NST": chr_,
            "MAF_≥": threshold,

            "gt_both_correct": gt_both_correct,
            "gt_mom_correct": gt_mom_correct,
            "gt_child_correct": gt_child_correct,
            "gt_mom_correct_child_wrong": gt_mom_correct_child_wrong,
            "gt_child_correct_mom_wrong": gt_child_correct_mom_wrong,
            "gt_both_wrong": gt_both_wrong,
            "gt_mom_wrong": gt_mom_wrong,
            "gt_child_wrong": gt_child_wrong,
            "total_gt_child": total_gt_child,
            "total_gt_mother": total_gt_mother,
            
            "alt_both_correct": alt_both_correct,
            "alt_mom_correct": alt_mom_correct,
            "alt_child_correct": alt_child_correct,
            "alt_mom_correct_child_wrong": alt_mom_correct_child_wrong,
            "alt_child_correct_mom_wrong": alt_child_correct_mom_wrong,
            "alt_both_wrong": alt_both_wrong,
            "alt_mom_wrong": alt_mom_wrong,
            "alt_child_wrong": alt_child_wrong,
            "total_alt_child": alt_gt_child,
            "total_alt_mother": alt_gt_mother,
        })
    return pd.DataFrame(rows)

def process_summary(df_path, gt_stats, chr_, mother, child):
    df = pd.read_csv(df_path)
    df["MAF"] = df.index  
    df["NST"] = chr_    

    bin_df = compute_bin_stats(df, gt_stats, chr_, mother, child)
    maf_df = compute_cumsum_stats(df, gt_stats, chr_, mother, child)

    return bin_df, maf_df


def plot_stats(df_bin, coverage, ff, type, Ox):
    stats = ['call_rate_gt', 'accuracy_gt', 'call_rate_alt', 'accuracy_alt']  
    statnames = [
        "Tỷ lệ gọi kiểu gen",          
        "Độ chính xác kiểu gen",      
        "Tỷ lệ gọi kiểu gen có allen thay thế",       
        "Độ chính xác của allen thay thế"      
    ]

    for stat in stats:
        fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        handles = []
        labels = []

        for i, relation in enumerate(['mother', 'child']):
            ax = axs[i]
            for nst in df_bin['NST'].unique():
                nst_data = df_bin[df_bin['NST'] == nst]
                maf_values = nst_data[Ox]
                accuracy_values = nst_data[f"{stat}_{relation}"]

                if nst != 'all':
                    line, = ax.plot(maf_values, accuracy_values, label=nst, linestyle='-', alpha=0.6, linewidth=1)
                else:
                    line, = ax.plot(maf_values, accuracy_values, label=nst, color="blue", linestyle='-', alpha=0.9, linewidth=3)

                if i == 0:
                    handles.append(line)
                    labels.append(nst)

            ax.set_title(
                f"{statnames[stats.index(stat)]} so với hệ gen của {'mẹ' if relation == 'mother' else 'thai nhi'}\n"
                f"(Độ bao phủ {coverage}x, tỷ lệ DNA thai nhi {ff * 100:.0f}%)",
                fontsize=11
            )
            ax.set_xlabel('Tần số alen (MAF)', fontsize=10)
            if i == 0:
                ax.set_ylabel('Giá trị', fontsize=10)
            ax.grid(True)

        fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=9)
        fig.tight_layout(pad=3.0, rect=[0, 0.1, 1, 1]) 

        filename = f"{stat}_{coverage}_{ff}_{type}.png"
        plt.savefig(os.path.join(PATHS["plot_directory"], "png", "nipt", filename), dpi=180)
        plt.close()
        

def plot_stats_multiple(data_by_coverage, stat, relation, Ox, type):
    statnames = {
        'call_rate_gt': "Tỷ lệ gọi kiểu gen thành công",
        'accuracy_gt': "Độ chính xác của kiểu gen",
        'call_rate_alt': "Tỷ lệ gọi kiểu gen có allen thay thế thành công",
        'accuracy_alt': "Độ chính xác của allen thay thế"
    }

    relation_vn = "so với mẹ" if relation == "mother" else "so với con"
    type_vn = "theo MAF" if type == "bin" else "theo MAF tối thiểu"

    coverages = list(data_by_coverage.keys()) 
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharey=True)

    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = [axes]

    for idx, coverage in enumerate(coverages):
        ax = axes[idx]
        ff_datasets = data_by_coverage[coverage]

        for ff, df in ff_datasets.items():
            nst_data = df[df['NST'] == 'all']
            maf_values = nst_data[Ox]
            y_values = nst_data[f"{stat}_{relation}"]

            ax.plot(maf_values, y_values, label=f"{int(ff*100)}%", linestyle='-', alpha=0.8, linewidth=2)

        ax.set_title(f'Coverage {coverage}x', fontsize=14)
        ax.set_xlabel('MAF' if Ox == "MAF" else 'MAF tối thiểu', fontsize=12)
        if idx % 2 == 0:
            ax.set_ylabel(statnames[stat], fontsize=12)
        ax.grid(True)

    plt.suptitle(f'{statnames[stat]} ({relation_vn}) {type_vn}', fontsize=18)
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=10)

    filename = f"{stat}_{relation}_{type}.png"
    plt.savefig(os.path.join(PATHS["plot_directory"], "png", "nipt", filename), dpi=180)
    plt.close()




def plot_gt_stacked(data, maf_value):
    grouped_data = defaultdict(dict)
    ffs = set()
    coverages_set = set()

    for item in data:
        coverage = item["coverage"]
        ff = item["ff"]
        df = item["avg_maf"]

        row = df[df["MAF_≥"] == maf_value]
        if row.empty:
            continue

        row = row.iloc[0]
        grouped_data[coverage][ff] = {
            "gt_both": row["gt_both_correct"],
            "gt_mother": row["gt_mom_correct_child_wrong"],
            "gt_child": row["gt_child_correct_mom_wrong"],
            "gt_false": row["gt_both_wrong"],
            "gt_not_called": row["total_gt_mother"] - row["gt_both_correct"] - row["gt_mom_correct_child_wrong"] - row["gt_child_correct_mom_wrong"] - row["gt_both_wrong"],

            "alt_both": row["alt_both_correct"],
            "alt_mother": row["alt_mom_correct_child_wrong"],
            "alt_child": row["alt_child_correct_mom_wrong"],
            "alt_false": row["alt_both_wrong"],
            "alt_not_called": row["total_alt_mother"] - row["alt_both_correct"] - row["alt_mom_correct_child_wrong"] - row["alt_child_correct_mom_wrong"] - row["alt_both_wrong"]
        }

        ffs.add(ff)
        coverages_set.add(coverage)

    coverages = sorted(coverages_set)
    ffs = sorted(ffs)
    num_ffs = len(ffs)

    bar_width = 0.15
    spacing = 0.05
    x_positions = list(range(len(coverages)))

    colors = {
        "both": "#4daf4a",
        "mother": "#377eb8",
        "child": "#984ea3",
        "false": "#e41a1c",
        "not_called": "#ff7f00"
    }

    def draw_subplot(ax, mode):
        for i, ff in enumerate(ffs):
            values_list = {k: [] for k in ["both", "mother", "child", "false", "not_called"]}
            positions = []

            for j, cov in enumerate(coverages):
                values = grouped_data[cov].get(ff)
                if not values:
                    continue

                prefix = "gt_" if mode == "GT" else "alt_"
                values_list["both"].append(values[prefix + "both"])
                values_list["mother"].append(values[prefix + "mother"])
                values_list["child"].append(values[prefix + "child"])
                values_list["false"].append(values[prefix + "false"])
                values_list["not_called"].append(values[prefix + "not_called"])

                offset = (i - num_ffs / 2) * (bar_width + spacing) + bar_width / 2
                positions.append(x_positions[j] + offset)

            bottoms = values_list["both"]
            ax.bar(positions, values_list["both"], color=colors["both"], width=bar_width)
            for key in ["mother", "child", "false", "not_called"]:
                ax.bar(positions, values_list[key], bottom=bottoms, color=colors[key], width=bar_width)
                bottoms = [bottoms[k] + values_list[key][k] for k in range(len(bottoms))]

        ax.set_xticks(x_positions)
        ax.set_xticklabels([str(cov) for cov in coverages])
        ax.set_xlabel("Độ bao phủ (Coverage)")
        ax.set_ylabel(f"Số lượng {'kiểu gen' if mode == 'GT' else 'biến thể'}")
        ax.set_title(f"{'Kiểu gen' if mode == 'GT' else 'Kiểu gen có allen thay thế'} với MAF ≥ {maf_value} khi tỉ lệ DNA thai nhi là {ff}", fontsize=11)
        ax.grid(True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    draw_subplot(ax1, mode="GT")
    draw_subplot(ax2, mode="ALT")

    ax1.set_ylim(0, max(ax1.get_ylim()[1], ax2.get_ylim()[1]))
    ax2.set_ylim(ax1.get_ylim())


    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors["both"]),
        plt.Rectangle((0, 0), 1, 1, color=colors["mother"]),
        plt.Rectangle((0, 0), 1, 1, color=colors["child"]),
        plt.Rectangle((0, 0), 1, 1, color=colors["false"]),
        plt.Rectangle((0, 0), 1, 1, color=colors["not_called"])
    ]
    labels = ["Cả hai đúng", "Chỉ mẹ đúng", "Chỉ con đúng", "Cả hai sai", "Không gọi được"]
    fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=9)

    plt.tight_layout(rect=[0, 0.1, 1, 1])  
    plt.savefig(os.path.join(PATHS["plot_directory"], "png", "nipt", f"gt_alt_{maf_value}.png"), dpi=180)
    plt.close()


def stats():
    output_dir = os.path.join(PATHS["plot_directory"])
    os.makedirs(output_dir, exist_ok=True)

    ground_truth_stats_file = PATHS["ground_truth_stat"]
    gt_stats = load_ground_truth_stats(ground_truth_stats_file)

    for coverage in PARAMETERS["coverage"]:
        for ff in PARAMETERS["ff"]:
            print(f"📊 Processing coverage = {coverage}")
            
            all_stats_bin = []
            all_stats_maf = []

            for trio_name, trio_info in TRIO_DATA.items():
                mother = trio_info["mother"]  
                father = trio_info["father"]
                child = trio_info["child"]

                fq_path = os.path.join(fastq_nipt_path(child, mother, father, coverage, ff), f"{child}_{mother}_{father}.fastq.gz")

                df_bin = []
                df_maf = []

                for chr_ in PARAMETERS["chrs"]:
                    df_path = statistic_summary(fq_path, chr_)
                    if not os.path.exists(df_path):
                        continue

                    stats_bin_chr, stats_maf_chr = process_summary(df_path, gt_stats, chr_, mother, child)
                    df_bin.append(stats_bin_chr)
                    df_maf.append(stats_maf_chr)


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

                all_stats_bin.append(df_bin)
                all_stats_maf.append(df_maf)

            df_bin = pd.concat(all_stats_bin, ignore_index=True)
            df_maf = pd.concat(all_stats_maf, ignore_index=True)

            avg_bin = df_bin.groupby(["MAF", "NST"], as_index=False).mean()
            avg_maf = df_maf.groupby(["MAF_≥", "NST"], as_index=False).mean()

            avg_bin.to_csv(os.path.join(output_dir, f"bin_coverage_{coverage}x_{ff}.csv"), index=False)
            avg_maf.to_csv(os.path.join(output_dir, f"maf_coverage_{coverage}x_{ff}.csv"), index=False)

            print(f"✅ Done saving stats for coverage {coverage}")


def recal_stats(df):           
    df['gt_mom_wrong'] = df["gt_both_wrong"] + df["gt_child_correct_mom_wrong"]
    df['alt_mom_wrong'] = df["alt_both_wrong"] + df["alt_child_correct_mom_wrong"]
    df['gt_child_wrong'] = df["gt_both_wrong"] + df["gt_mom_correct_child_wrong"]
    df['alt_child_wrong'] = df["alt_both_wrong"] + df["alt_mom_correct_child_wrong"]
    

    df['call_rate_gt_mother'] = (df['gt_mom_correct'] + df['gt_mom_wrong']) / df['total_gt_mother']
    df['accuracy_gt_mother'] = df['gt_mom_correct'] / (df['gt_mom_correct'] + df['gt_mom_wrong'])
    df['call_rate_alt_mother'] = (df['alt_mom_correct'] + df['alt_mom_wrong']) / df['total_alt_mother']
    df['accuracy_alt_mother'] = df['alt_mom_correct'] / (df['alt_mom_correct'] + df['alt_mom_wrong'])


    df['call_rate_gt_child'] = (df['gt_child_correct'] + df['gt_child_wrong']) / df['total_gt_child']
    df['accuracy_gt_child'] = df['gt_child_correct'] / (df['gt_child_correct'] + df['gt_child_wrong'])
    df['call_rate_alt_child'] = (df['alt_child_correct'] + df['alt_child_wrong']) / df['total_alt_child']
    df['accuracy_alt_child'] = df['alt_child_correct'] / (df['alt_child_correct'] + df['alt_child_wrong'])

    return df


def plot_all():
    output_dir = os.path.join(PATHS["plot_directory"])

    data_bin_by_stat = {stat: {'mother': {}, 'child': {}} for stat in ['call_rate_gt', 'accuracy_gt', 'call_rate_alt', 'accuracy_alt']}
    data_maf_by_stat = {stat: {'mother': {}, 'child': {}} for stat in ['call_rate_gt', 'accuracy_gt', 'call_rate_alt', 'accuracy_alt']}

    pl_stats = []

    for coverage in PARAMETERS["coverage"]:
        for ff in PARAMETERS["ff"]:
            print(f"📊 Đang xử lý coverage = {coverage}, ff = {ff}")

            output_file_bin = os.path.join(output_dir, f"bin_coverage_{coverage}x_{ff}.csv")
            output_file_maf = os.path.join(output_dir, f"maf_coverage_{coverage}x_{ff}.csv")

            avg_bin = pd.read_csv(output_file_bin)
            avg_maf = pd.read_csv(output_file_maf)

            for stat in data_bin_by_stat:
                data_bin_by_stat[stat]['mother'].setdefault(coverage, {})[ff] = avg_bin
                data_bin_by_stat[stat]['child'].setdefault(coverage, {})[ff] = avg_bin
                data_maf_by_stat[stat]['mother'].setdefault(coverage, {})[ff] = avg_maf
                data_maf_by_stat[stat]['child'].setdefault(coverage, {})[ff] = avg_maf

            pl_stats.append({
                "coverage": coverage,
                "ff": ff,
                "avg_maf": avg_maf
            })

    for stat in data_bin_by_stat:
        for relation in ['mother', 'child']:
            plot_stats_multiple(data_bin_by_stat[stat][relation], stat, relation, "MAF", "bin")
            plot_stats_multiple(data_maf_by_stat[stat][relation], stat, relation, "MAF_≥", "maf")

    print("🎉 Tất cả biểu đồ đã được lưu.")
    plot_gt_stacked(pl_stats, 0)
    plot_gt_stacked(pl_stats, 1)
    plot_gt_stacked(pl_stats, 5)



if __name__ == "__main__":
    stats()
    plot_all()
