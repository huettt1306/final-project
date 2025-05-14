import os, sys, argparse
from steps.alignment import run_multiple_alignment
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

    return parser.parse_args()

def main():
    #run_prepare_reference_panel()
    
    args = parse_arguments()
    fqlist = args.input_path
    print(f"Start pipeline for samples in {fqlist}...")

    # Step 1: Alignment and statistics
    print(f"Start alignment for {fqlist}...")
    run_multiple_alignment(fqlist)

    # Step2: SNP detection and allele frequency estimation
    print(f"Start basevar for {fqlist}...")
    run_basevar(fqlist)

    # Step3: Genotype imputation
    print(f"Start glimpse for {fqlist}...")
    run_glimpse(fqlist)

    print(f"Done pipeline for sample in {fqlist}.")
if __name__ == "__main__":
    main()
