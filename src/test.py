from helper.converter import convert_cram_to_fastq
from helper.path_define import cram_path, fastq_path
from helper.metrics import compare_fastq_sequences, evaluate_vcf
import sys, os
import pandas as pd
from helper.file_utils import process_vcf
from statistic.GT import get_af_gt_different

def compare_single_variants(sample1, sample2):
    """
    So sánh biến thể giữa các phương pháp và lưu kết quả ra file CSV.
    """
    try:
        sample1_df = process_vcf(sample1, "sample1")
        sample2_df = process_vcf(sample2, "sample2")

        merged_df = pd.merge(sample1_df, sample2_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")

        return merged_df
    except Exception as e:
        raise


def update_stats(stats, af, field):
    """
    Kiểm tra xem af đã có trong stats chưa, nếu chưa thì thêm mới, nếu có thì cập nhật trường 'field'.
    """
    if af < 0:
        return stats
    
    if af not in stats:
        stats[af] = {
            "GT Dffference": 0,
        }

    # Cập nhật giá trị của trường 'field'
    stats[af][field] += 1
    return stats


def calculate_af_single_statistics(df):
    """
    Tính toán thống kê cho từng giá trị AF và trả về dưới dạng dictionary.
    Mỗi entry trong dictionary sẽ chứa các thống kê cho một giá trị AF.
    """
    stats = {}

    for _, row in df.iterrows():
        stats = update_stats(stats, get_af_gt_different(row, "sample1", "sample2"), "GT Dffference")
    print(stats)
    return stats


def main():
    file1 = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/sample1/glimpse_output/imputed_file_merged/glimpse.chr20_imputed.vcf.gz"
    file2 = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/sample2/glimpse_output/imputed_file_merged/glimpse.chr20_imputed.vcf.gz"
    
    b1 = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/sample1/basevar_output/NIPT_basevar_chr20.vcf.gz"
    b2 = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/sample2/basevar_output/NIPT_basevar_chr20.vcf.gz"
    df = compare_single_variants(file1, file2)
    calculate_af_single_statistics(df)


if __name__ == "__main__":
    main()


