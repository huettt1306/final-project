import pandas as pd
import os
from statistics.GT import get_af_gt, get_af_gt_true, get_af_gt_false, get_af_gt_not_given
from statistics.ALT import get_af_alt, get_af_alt_true, get_af_alt_false, get_af_alt_not_given
from helper.config import PATHS
from helper.logger import setup_logger


logger = setup_logger(os.path.join(PATHS["logs"], "single_statistic_pipeline.log"))


def compare_single_variants(ground_truth_df, basevar_df, glimpse_df):
    logger.info("Comparing single variants...")

    merged_df = pd.merge(glimpse_df, basevar_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
    merged_df = pd.merge(ground_truth_df, merged_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")

    merged_df.drop(columns=['AF_BaseVar'], inplace=True)
    merged_df.drop(columns=['AF_Glimpse'], inplace=True)
    merged_df.rename(columns={'AF_Truth': 'AF'}, inplace=True)

    merged_df["AF"].fillna(-1.0, inplace=True)
    merged_df.fillna(False, inplace=True)

    logger.info(f"Finished compare single variants.")
    return merged_df
    

def update_stats(stats, af, field):
    if af < 0:
        return stats
    
    if af not in stats:
        stats[af] = {
            "GT Truth": 0,
            "GT Glimpse": 0,
            "GT Glimpse True": 0,
            "GT Glimpse False": 0,
            "GT Truth not found": 0,

            "ALT Truth": 0,
            "ALT Basevar": 0,
            "ALT Glimpse": 0,
            "ALT Glimpse True": 0,
            "ALT Glimpse False": 0,
            "ALT Truth not found": 0,
        }

    stats[af][field] += 1
    return stats


def calculate_af_single_statistics(df):
    stats = {}

    for _, row in df.iterrows():
        stats = update_stats(stats, get_af_gt(row, "Truth"), "GT Truth")
        stats = update_stats(stats, get_af_gt(row, "Glimpse"), "GT Glimpse")
        stats = update_stats(stats, get_af_gt_true(row, "Glimpse", "Truth"), "GT Glimpse True")
        stats = update_stats(stats, get_af_gt_false(row, "Glimpse", "Truth"), "GT Glimpse False")
        stats = update_stats(stats, get_af_gt_not_given(row, "Truth", "Glimpse"), "GT Truth not found")

        stats = update_stats(stats, get_af_alt(row, "Truth"), "ALT Truth")
        stats = update_stats(stats, get_af_alt(row, "BaseVar"), "ALT Basevar")
        stats = update_stats(stats, get_af_alt(row, "Glimpse"), "ALT Glimpse")
        stats = update_stats(stats, get_af_alt_true(row, "Glimpse", "Truth"), "ALT Glimpse True")
        stats = update_stats(stats, get_af_alt_false(row, "Glimpse", "Truth"), "ALT Glimpse False")
        stats = update_stats(stats, get_af_alt_not_given(row, "Truth", "Glimpse"), "ALT Truth not found")

    return stats







