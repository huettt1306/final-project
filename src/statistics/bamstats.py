import pysam

def collect_stats(bam_path):
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    total_reads = 0
    mapped_reads = 0
    properly_paired = 0
    mapq_ge_20 = 0
    mismatches = 0
    total_bases = 0
    insert_sizes = []

    for read in bamfile.fetch(until_eof=True):
        total_reads += 1

        if not read.is_unmapped:
            mapped_reads += 1

            if read.is_proper_pair:
                properly_paired += 1

            if read.mapping_quality >= 20:
                mapq_ge_20 += 1

            try:
                nm = read.get_tag("NM")
                mismatches += nm
                total_bases += read.query_length
            except KeyError:
                pass

            if read.template_length > 0:
                insert_sizes.append(abs(read.template_length))

    bamfile.close()

    avg_mismatch_rate = mismatches / total_reads if total_reads > 0 else 0
    avg_insert_size = sum(insert_sizes) / len(insert_sizes) if insert_sizes else 0

    return {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "properly_paired": properly_paired,
        "mapq_ge_20": mapq_ge_20,
        "mismatch_rate_per_read": avg_mismatch_rate,
        "avg_insert_size": avg_insert_size
    }

def print_comparison(stats1, stats2, label1="Before", label2="After"):
    print(f"{'Metric':<30}{label1:<15}{label2:<15}")
    print("-" * 60)
    for key in stats1:
        val1 = stats1[key]
        val2 = stats2[key]
        print(f"{key:<30}{val1:<15.5f}{val2:<15.5f}")

bam_before = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/0.1x/HG02019/sample_1/1tmp_files/HG02019.sorted.rmdup.bam"
bam_after = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/0.1x/HG02019/sample_1/1tmp_files/HG02019.sorted.rmdup.realign.bam"

stats_before = collect_stats(bam_before)
stats_after = collect_stats(bam_after)

print_comparison(stats_before, stats_after)
