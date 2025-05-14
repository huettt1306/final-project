import os
import glob
import csv
from collections import defaultdict
from cyvcf2 import VCF

vcf_dir = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/vcf_"

stats = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {
    'Total_GT': 0,
    'ALT_GT': 0,
    'HET_GT': 0,
    'HOM_ALT': 0
})))

chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

for chr_name in chrs:
    file = os.path.join(vcf_dir, f"{chr_name}_variants.vcf.gz")
    vcf = VCF(file)

    for record in vcf:
        af = dict(record.INFO).get('AF', None)
        if af is None or not (0 <= af <= 1):
            continue

        maf = min(af, 1 - af)
        if maf < 0.001:
            continue

        maf_bin = int(maf * 100)  
        for i, sample_gt in enumerate(record.genotypes):
            sample_name = vcf.samples[i]

            # sample_gt: [allele1, allele2, phased] OR [allele1, None, phased] OR [allele1]
            alleles = [a for a in sample_gt[:2] if a is not None and a >= 0]
            gt_set = set(alleles) if alleles else None

            if gt_set is not None:
                stats[chr_name][sample_name][maf_bin]['Total_GT'] += 1

                if any(a > 0 for a in gt_set):
                    stats[chr_name][sample_name][maf_bin]['ALT_GT'] += 1

                if len(alleles) == 2:
                    if alleles[0] != alleles[1]:
                        stats[chr_name][sample_name][maf_bin]['HET_GT'] += 1
                    elif alleles[0] > 0:
                        stats[chr_name][sample_name][maf_bin]['HOM_ALT'] += 1

with open(os.path.join(vcf_dir, "genotype_stats_per_chr.csv"), 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Chr', 'Sample', 'MAF', 'Total_GT', 'ALT_GT', 'HET_GT', 'HOM_ALT'])
    for chr_name in chrs:
        for sample in sorted(stats[chr_name]):
            for maf_bin in sorted(stats[chr_name][sample].keys()):
                s = stats[chr_name][sample][maf_bin]
                writer.writerow([
                    chr_name,
                    sample,
                    maf_bin,
                    s['Total_GT'],
                    s['ALT_GT'],
                    s['HET_GT'],
                    s['HOM_ALT']
                ])
