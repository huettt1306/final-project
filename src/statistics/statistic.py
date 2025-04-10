import pandas as pd
import os
from helper.file_utils import save_results_to_csv, process_vcf
from helper.path_define import ground_truth_vcf, statistic_summary, glimpse_vcf, basevar_vcf, samid
from helper.config import PATHS, PARAMETERS
from helper.logger import setup_logger
from statistics.single_stats import compare_single_variants, calculate_af_single_statistics
from statistics.nipt_stats import compare_nipt_variants, calculate_af_nipt_statistics

logger = setup_logger(os.path.join(PATHS["logs"], "statistic_pipeline.log"))

def process_dataframe(df):
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    df = df.sort_values(by="AF (%)", ascending=True)

    for col in df.columns:
        if col != "AF (%)":
            df[f'Total {col}'] = df[col][::-1].cumsum()[::-1]

    return df


def generate_summary_statistics(df, output_file, type="single"):
    try:
        logger.info("Generating summary statistics....")

        stats_data = {"AF (%)": []}

        if type == "nipt":
            stats = calculate_af_nipt_statistics(df)
        else:
            stats = calculate_af_single_statistics(df)

        for af_percent, stat_values in stats.items():
            stats_data["AF (%)"].append(int(af_percent))

            for key, value in stat_values.items():
                if key not in stats_data:
                    stats_data[key] = []
                stats_data[key].append(value)

        summary_df = process_dataframe(pd.DataFrame(stats_data))

        save_results_to_csv(output_file, summary_df)

        logger.info("Summary statistics generated successfully")
        return summary_df

    except Exception as e:
        logger.error(f"Error generating summary statistics: {e}")
        raise



def statistic(fq, chromosome):
    sample_name = samid(fq)

    try:
        logger.info(f"Starting statistics for sample {sample_name} on chromosome {chromosome}")
        
        basevar_df = process_vcf(basevar_vcf(fq, chromosome), "BaseVar")
        glimpse_df = process_vcf(glimpse_vcf(fq, chromosome), "Glimpse")

        if "_" not in sample_name:
            ground_truth_df = process_vcf(ground_truth_vcf(chromosome), "Truth", sample_name)
            df = compare_single_variants(ground_truth_df, basevar_df, glimpse_df)
            return generate_summary_statistics(df, statistic_summary(fq, chromosome))
            
        else:
            child, mom, dad = sample_name.split("_")
            child_df = process_vcf(ground_truth_vcf(chromosome), "Child", child)
            mom_df = process_vcf(ground_truth_vcf(chromosome), "Mother", mom)
            dad_df = process_vcf(ground_truth_vcf(chromosome), "Father", dad)
            df = compare_nipt_variants(child_df, mom_df, dad_df, basevar_df, glimpse_df)
            return generate_summary_statistics(df, statistic_summary(fq, chromosome), "nipt")
            
    except Exception as e:
        logger.error(f"Error in statistic function for sample {sample_name} and chromosome {chromosome}: {e}")
        raise


def run_statistic(fq):
    try:
        logger.info(f"Starting full statistics pipeline for sample {fq}")
        all_data = [] 

        for chromosome in PARAMETERS["chrs"]:
            df = statistic(fq, chromosome)

            if df is not None and not df.empty:
                all_data.append(df)  

        all_df = pd.concat(all_data).groupby("AF (%)", as_index=False).sum() 

        save_results_to_csv(statistic_summary(fq), all_df)
        logger.info(f"Completed full statistics pipeline for sample {fq}")

        return all_df
    
    except Exception as e:
        logger.error(f"Error in run_statistic pipeline for sample {fq}: {e}")
        raise