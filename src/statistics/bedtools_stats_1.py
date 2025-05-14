import gzip
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

def compute_avg_cov_per_bin(filepath, bin_size=5_000_000):
    bin_cov = defaultdict(int)
    bin_len = defaultdict(int)

    with gzip.open(filepath, 'rt') as f:
        for line in f:
            chrom, start, end, cov = line.strip().split('\t')
            start, end, cov = int(start), int(end), int(cov)
            length = end - start
            midpoint = (start + end) // 2
            bin_index = midpoint // bin_size
            bin_cov[bin_index] += cov * length
            bin_len[bin_index] += length

    avg_cov = {
        bin_index: bin_cov[bin_index] / bin_len[bin_index]
        for bin_index in bin_cov
    }
    return avg_cov

file_paths = {
    "0.1x": "0.1x.coverage.bed.gz",
    "0.2x": "0.2x.coverage.bed.gz",
    "0.5x": "0.5x.coverage.bed.gz",
    "1x":   "1x.coverage.bed.gz"
}

all_results = {
    label: compute_avg_cov_per_bin(path)
    for label, path in file_paths.items()
}

all_bins = sorted(set(bin for res in all_results.values() for bin in res))

plot_data = {
    label: [all_results[label].get(bin, 0) for bin in all_bins]
    for label in file_paths
}

plt.figure(figsize=(16, 5))
for label, values in plot_data.items():
    plt.plot([bin * 5 for bin in all_bins], values, label=label, linewidth=1)

plt.title("Coverage trung bình theo từng bin 5MB trên toàn genome")
plt.xlabel("Vị trí trên genome (Mb)")
plt.ylabel("Coverage trung bình")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(title="Mẫu")
plt.tight_layout()
plt.savefig(
        "bedtoold_cov_1.png",
        dpi=300
    )
plt.show()
