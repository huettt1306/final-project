import os, sys, pandas as pd
from helper.config import PARAMETERS, TRIO_DATA, PATHS
from helper.logger import setup_logger
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from helper.path_define import ground_truth_vcf, fastq_single_path, fastq_nipt_path, basevar_vcf, glimpse_vcf, samid, statistic_summary
from helper.record import get_key_from_record, compare_with_ground_truth
from cyvcf2 import VCF
from helper.file_utils import save_results_to_csv


logger = setup_logger(os.path.join(PATHS["logs"], "statistic.log"))

def read_vcf_by_variant(vcf_path):
    vcf = VCF(vcf_path)

    variant_data = {}

    for record in vcf:
        key = get_key_from_record(record)

        af = dict(record.INFO).get('AF', -1)
        if af is None or not (0 <= af <= 1):
            continue 

        maf = min(af, 1 - af)
        if maf < 0.001:
            continue  

        sample_gts = {}
        for i, sample_gt in enumerate(record.genotypes):
            sample_name = vcf.samples[i]
            # sample_gt: [allele1, allele2, phased] OR just [allele1, None, phased] OR [allele1]
            alleles = [a for a in sample_gt[:2] if a is not None and a >= 0]
            gt_set = set(alleles) if alleles else None
            sample_gts[sample_name] = gt_set


        variant_data[key] = {
            "MAF": int(maf * 100),
            "samples": sample_gts
        }

    return variant_data



def stats_single(ground_truth, chromosome, fq):
    print(f"Start stats single for {fq}")
    basevar = basevar_vcf(fq, chromosome)
    glimpse = glimpse_vcf(fq, chromosome)
    name = samid(fq)

    columns = [
        "basevar_GT_true", "basevar_GT_false",
        "basevar_ALT_true", "basevar_ALT_false",
        "glimpse_GT_true", "glimpse_GT_false",
        "glimpse_ALT_true", "glimpse_ALT_false"
    ]
    df = pd.DataFrame(0, index=range(0, 51), columns=columns)

    def update_stats(result, prefix):
        maf = result["MAF"]  
        df.loc[maf, f"{prefix}_GT_true"] += result["GT_correct"]
        df.loc[maf, f"{prefix}_GT_false"] += result["GT_wrong"]
        df.loc[maf, f"{prefix}_ALT_true"] += result["ALT_correct"]
        df.loc[maf, f"{prefix}_ALT_false"] += result["ALT_wrong"]

    for record in VCF(basevar):
        result = compare_with_ground_truth(record, ground_truth, name)
        if result:
            update_stats(result, "basevar")

    for record in VCF(glimpse):
        result = compare_with_ground_truth(record, ground_truth, name)
        if result:
            update_stats(result, "glimpse")

    df.index.name = "MAF"
    print(df)
    save_results_to_csv(statistic_summary(fq, chromosome), df)

    print(f"Stats single saved to csv.")
    return df  


def stats_nipt(ground_truth, chromosome, fq):
    print(f"Start stats nipt for {fq}")
    glimpse = glimpse_vcf(fq, chromosome)
    sample_name = samid(fq)
    child, mom, dad = sample_name.split("_")

    columns = [
        # GT
        "gt_both_correct",
        "gt_mom_correct",
        "gt_child_correct",
        "gt_mom_correct_child_wrong",
        "gt_child_correct_mom_wrong",
        "gt_both_wrong",

        # ALT
        "alt_both_correct",
        "alt_mom_correct",
        "alt_child_correct",
        "alt_mom_correct_child_wrong",
        "alt_child_correct_mom_wrong",
        "alt_both_wrong"
    ]
    df = pd.DataFrame(0, index=range(0, 51), columns=columns)

    def update_stats(result_child, result_mom):
        maf = result_child["MAF"]

        # --- GT ---
        gt_child_ok = result_child["GT_correct"] == 1
        gt_mom_ok = result_mom["GT_correct"] == 1

        if gt_child_ok and gt_mom_ok:
            df.loc[maf, "gt_both_correct"] += 1
        elif gt_mom_ok and not gt_child_ok:
            df.loc[maf, "gt_mom_correct_child_wrong"] += 1
        elif gt_child_ok and not gt_mom_ok:
            df.loc[maf, "gt_child_correct_mom_wrong"] += 1
        elif not gt_child_ok and not gt_mom_ok:
            df.loc[maf, "gt_both_wrong"] += 1

        df.loc[maf, "gt_mom_correct"] += gt_mom_ok
        df.loc[maf, "gt_child_correct"] += gt_child_ok

        # --- ALT ---
        alt_child_ok = result_child["ALT_correct"] == 1
        alt_mom_ok = result_mom["ALT_correct"] == 1

        if alt_child_ok and alt_mom_ok:
            df.loc[maf, "alt_both_correct"] += 1
        elif alt_mom_ok and not alt_child_ok:
            df.loc[maf, "alt_mom_correct_child_wrong"] += 1
        elif alt_child_ok and not alt_mom_ok:
            df.loc[maf, "alt_child_correct_mom_wrong"] += 1
        elif not alt_child_ok and not alt_mom_ok:
            df.loc[maf, "alt_both_wrong"] += 1

        df.loc[maf, "alt_mom_correct"] += alt_mom_ok
        df.loc[maf, "alt_child_correct"] += alt_child_ok

    for record in VCF(glimpse):
        result_mom = compare_with_ground_truth(record, ground_truth, mom)
        result_child = compare_with_ground_truth(record, ground_truth, child)

        if result_mom and result_child:
            update_stats(result_child, result_mom)

    df.index.name = "MAF"
    print(df)
    save_results_to_csv(statistic_summary(fq, chromosome), df)
    print(f"Stats nipt saved to csv")
    return df



def statistic(chromosome) :
    print(f"Loading ground truth for {chromosome}")
    ground_truth = read_vcf_by_variant(ground_truth_vcf(chromosome))
    print(f"Loaded ground truth for {chromosome}")
    
    for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
        logger.info(f"######## PROCESSING index {index} ########")

        for trio_name, trio_info in TRIO_DATA.items():
            father_name = trio_info["father"]
            mother_name = trio_info["mother"]
            child_name = trio_info["child"]

            with ThreadPoolExecutor(max_workers=4) as executor:
                executor.map(
                    lambda coverage: stats_single(
                        ground_truth, chromosome,
                        os.path.join(fastq_single_path(mother_name, coverage, index), f"{mother_name}.fastq.gz")
                    ),
                    PARAMETERS["coverage"]
                )

            for coverage in PARAMETERS["coverage"]:
                with ThreadPoolExecutor(max_workers=3) as executor:
                    executor.map(
                        lambda ff: stats_nipt(
                            ground_truth, chromosome,
                            os.path.join(fastq_nipt_path(child_name, mother_name, father_name, coverage, ff, index), f"{child_name}_{mother_name}_{father_name}.fastq.gz")
                        ),
                        PARAMETERS["ff"]
                    )


def main():
    if len(sys.argv) < 2:
        logger.error("Please provide a chr to process.")
        sys.exit(1)

    statistic(sys.argv[1])


if __name__ == "__main__":
    main()