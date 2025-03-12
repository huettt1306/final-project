import os, sys, argparse
from steps.alignment import run_alignment_pipeline
from steps.basevar import run_basevar
from steps.glimpse import run_glimpse
from steps.reference_panel_prepare import run_prepare_reference_panel
from helper.config import PARAMETERS

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input_path',
        type=str,
        help="path/to/fastq"
    )

    parser.add_argument(
        '--gender',
        type=int,
        choices=[1, 2],
        required=True,
        help="1 for male, 2 for female"
    )

    return parser.parse_args()

def main():
    #run_prepare_reference_panel()
    
    args = parse_arguments()
    fastq_dir = args.input_path
    PARAMETERS["gender"] =  args.gender
    print(f"Start pipeline for sample in {fastq_dir}...")

    # Step 1: Alignment and statistics
    print(f"Start alignment for {fastq_dir}...")
    run_alignment_pipeline(fastq_dir)

    # Step2: SNP detection and allele frequency estimation
    print(f"Start basevar for {fastq_dir}...")
    run_basevar(fastq_dir)

    # Step3: Genotype imputation
    print(f"Start glimpse for {fastq_dir}...")
    run_glimpse(fastq_dir)

    print(f"Done pipeline for sample in {fastq_dir}.")
if __name__ == "__main__":
    main()
