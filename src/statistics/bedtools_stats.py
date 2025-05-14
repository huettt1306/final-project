import gzip
import matplotlib.pyplot as plt
from collections import defaultdict

def compute_avg_cov_per_chrom(filepath):
    total_cov = defaultdict(int)
    total_len = defaultdict(int)

    with gzip.open(filepath, 'rt') as f:
        for line in f:
            chrom, start, end, cov = line.strip().split('\t')
            start, end, cov = int(start), int(end), int(cov)
            length = end - start
            total_cov[chrom] += cov * length
            total_len[chrom] += length

    avg_cov = {
        chrom: total_cov[chrom] / total_len[chrom]
        for chrom in total_cov
    }
    return avg_cov

file_paths = {
    "0.1x": "0.1x.coverage.bed.gz",
    "0.2x": "0.2x.coverage.bed.gz",
    "0.5x": "0.5x.coverage.bed.gz",
    "1x":   "1x.coverage.bed.gz"
}

all_results = {
    label: compute_avg_cov_per_chrom(path)
    for label, path in file_paths.items()
}

all_chroms = sorted(set(chrom for result in all_results.values() for chrom in result))
plot_data = {
    chrom: [all_results[label].get(chrom, 0) for label in file_paths]
    for chrom in all_chroms
}

import numpy as np
labels = list(file_paths.keys())
x = np.arange(len(all_chroms))
width = 0.2

plt.figure(figsize=(16, 6))
for i, label in enumerate(labels):
    values = [plot_data[chrom][i] for chrom in all_chroms]
    plt.bar(x + i*width, values, width, label=label)

plt.xticks(x + width * (len(labels) - 1) / 2, all_chroms, rotation=45, ha='right')
plt.ylabel("Coverage trung bình")
plt.xlabel("Nhiễm sắc thể")
plt.title("Coverage trung bình theo từng NST")
plt.legend(title="Mẫu")
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(
        "bedtoold_cov.png",
        dpi=300
    )
plt.show()
