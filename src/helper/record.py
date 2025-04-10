def get_key_from_record(record):
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = record.ALT[0]
    return (chrom, pos, ref, alt)

def compare_with_ground_truth(record, ground_truth, sample_name):
    key = get_key_from_record(record)

    if key not in ground_truth:
        return None

    if sample_name not in ground_truth[key]["samples"]:
        return None

    gt = record.genotypes[0]
    if gt is None:
        return None

    alleles = [a for a in gt[:2] if a is not None and a >= 0]
    gt_set = set(alleles) if alleles else None
    gt_truth = ground_truth[key]["samples"][sample_name]

    result = {
        "GT_correct": 0,
        "GT_wrong": 0,
        "ALT_correct": 0,
        "ALT_wrong": 0,
        "MAF": ground_truth[key]["MAF"],
    }

    # So sánh chính xác GT
    if gt_set == gt_truth:
        result["GT_correct"] += 1
    else:
        result["GT_wrong"] += 1

    # Đánh giá ALT
    if gt_set != {0}:  # dự đoán có biến thể
        if gt_truth != {0}:  # ground truth cũng có biến thể
            result["ALT_correct"] += 1
        else:
            result["ALT_wrong"] += 1

    return result

