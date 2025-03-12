import os
import subprocess
import shutil
from helper.config import TOOLS, PATHS, PARAMETERS
from helper.path_define import samid, tmp_outdir, batch1_final_outdir, bamlist_dir
from helper.logger import setup_logger

REF = PATHS["ref"]
GATK_BUNDLE_DIR = PATHS["gatk_bundle_dir"]
BWA = TOOLS["bwa"]
SAMTOOLS = TOOLS["samtools"]
GATK = TOOLS["gatk"]
JAVA = TOOLS["java"]
BEDTOOLS = TOOLS["bedtools"]
BGZIP = TOOLS["bgzip"]
TABIX=TOOLS["tabix"]

logger = setup_logger(os.path.join(PATHS["logs"], "alignment_pipeline.log"))


def run_bwa_alignment(fq, outdir):
    logger.info(f"Calculating {fq}. We'll save it in {outdir}")

    sai_file = os.path.join(outdir, f"{samid(fq)}.sai")
    bam_file = os.path.join(outdir, f"{samid(fq)}.bam")
    sorted_bam = os.path.join(outdir, f"{samid(fq)}.sorted.bam")
    rmdup_bam = os.path.join(outdir, f"{samid(fq)}.sorted.rmdup.bam")
    finish_flag = os.path.join(outdir, "bwa_sort_rmdup.finish")

    # Step 0: Verify the flag
    if os.path.exists(finish_flag):
        logger.info(f"BWA alignment for {fq} already exist. Skip alignment...")
        return

    # Step 1: BWA alignment
    logger.info("\nRunning BWA alignment...")    
    with open(sai_file, "w") as sai_out:
        subprocess.run([BWA, "aln", "-e", "10", "-t", f"{PARAMETERS['threads']}", "-i", "5", "-q", "0", REF, fq
        ], stdout=sai_out, check=True)

    with open(bam_file, "wb") as bam_out:
        bwa_process = subprocess.Popen([BWA, "samse", "-r",
            f"@RG\\tID:default\\tPL:COMPLETE\\tSM:{samid(fq)}",
            REF, sai_file, fq
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
        subprocess.run([SAMTOOLS, "view", "-h", "-Sb", "-@", f"{PARAMETERS['threads']}", "-"
        ], stdin=bwa_process.stdout, stdout=bam_out, check=True)
            
        bwa_process.stdout.close()
        bwa_process.wait() 

    logger.info("** BWA done **")

    # Step 2: Sorting BAM
    logger.info("Sorting BAM...")
    subprocess.run([SAMTOOLS, "sort", "-@", f"{PARAMETERS['threads']}", "-O", "bam", "-o", sorted_bam, bam_file], check=True)
    logger.info("** BAM sorted done **")

    # Step 3: Removing duplicates
    logger.info("Removing duplicates...")
    subprocess.run([SAMTOOLS, "markdup", "-@", f"{PARAMETERS['threads']}", sorted_bam, rmdup_bam], check=True)
    logger.info("** rmdup done **")

    # Step 4: Indexing BAM
    logger.info("Indexing BAM...")
    subprocess.run([SAMTOOLS, "index", "-@", f"{PARAMETERS['threads']}", rmdup_bam], check=True)
    logger.info("** index done **")

    # Step 5: Create finish flag
    with open(finish_flag, "w") as finish_file:
        finish_file.write("BWA alignment completed successfully.")
    logger.info(f"BWA alignment for {samid} completed successfully.")


def run_bwa_realign(sample_id, outdir):
    logger.info(f"\nStarting realign {samid} pipeline...")

    bam_file = os.path.join(outdir, f"{sample_id}.sorted.rmdup.bam")
    intervals_file = os.path.join(outdir, f"{sample_id}.indel_target_intervals.list")
    realigned_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.bam")
    finish_flag = os.path.join(outdir, "realigner.finish")

    # Step 0: Verify the flag
    if os.path.exists(finish_flag):
        logger.info(f"Realigner for {sample_id} already exist. Skip realign...")
        return

    # Step 1: RealignerTargetCreator
    logger.info("Running RealignerTargetCreator...")
    subprocess.run([JAVA, "-Xmx15g", "-jar", GATK,
        "-T", "RealignerTargetCreator",
        "-R", REF, "-I", bam_file,
        "-known", os.path.join(GATK_BUNDLE_DIR, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
        "-known", os.path.join(GATK_BUNDLE_DIR, "Homo_sapiens_assembly38.known_indels.vcf.gz"),
        "-o", intervals_file
    ], check=True)
    logger.info("** RealignerTargetCreator done **")

    # Step 2: IndelRealigner
    logger.info("Running IndelRealigner...")
    subprocess.run([JAVA, "-Xmx15g", "-jar", GATK,
        "-T", "IndelRealigner",
        "-R", REF, "-I", bam_file,
        "-known", os.path.join(GATK_BUNDLE_DIR, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
        "-known", os.path.join(GATK_BUNDLE_DIR, "Homo_sapiens_assembly38.known_indels.vcf.gz"),
        "--targetIntervals", intervals_file,
        "-o", realigned_bam
    ], check=True)
    logger.info("** IndelRealigner done **")

    # Create finish flag for IndelRealigner
    with open(finish_flag, "w") as flag:
        flag.write("Realign completed successfully.")
    logger.info(f"Realign pipeline for {samid} completed successfully.")


def run_bqsr(sample_id, outdir):
    realigned_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.bam")
    recal_table = os.path.join(outdir, f"{sample_id}.recal_data.table")
    bqsr_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bam")
    bsqr_flag = os.path.join(outdir, "bsqr.finish")

    # Step 0: Verify the flag
    if os.path.exists(bsqr_flag):
        logger.info(f"BSQR result for {samid} already exist. Skip bsqr...")
        return

    # Step 1: Index the realigned BAM
    logger.info("Indexing realigned BAM...")
    subprocess.run([SAMTOOLS, "index", "-@", f"{PARAMETERS['threads']}", realigned_bam], check=True)
    logger.info("** Index done **")

    # Step 2: BaseRecalibrator
    logger.info("Running BaseRecalibrator...")
    subprocess.run([JAVA, "-jar", GATK,
        "-T", "BaseRecalibrator",
        "-nct", "8",
        "-R", REF,
        "-I", realigned_bam,
        "--knownSites", os.path.join(GATK_BUNDLE_DIR, "Homo_sapiens_assembly38.dbsnp138.vcf.gz"),
        "--knownSites", os.path.join(GATK_BUNDLE_DIR, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
        "--knownSites", os.path.join(GATK_BUNDLE_DIR, "Homo_sapiens_assembly38.known_indels.vcf.gz"),
        "-o", recal_table
    ], check=True)
    logger.info("** BaseRecalibrator done **")

    # Step 3: PrintReads
    logger.info("Running PrintReads...")
    subprocess.run([JAVA, "-jar", GATK,
        "-T", "PrintReads",
        "-nct", "8",
        "-R", REF,
        "--BQSR", recal_table,
        "-I", realigned_bam,
        "-o", bqsr_bam
    ], check=True)
    logger.info("** PrintReads done **")
    
    # Step 4: Index the BQSR BAM
    logger.info("Indexing BQSR BAM...")
    subprocess.run([SAMTOOLS, "index", "-@", f"{PARAMETERS['threads']}", bqsr_bam], check=True)
    logger.info("** BAM index done **")

    # Create finish flag for BSQR
    with open(bsqr_flag, "w") as flag:
        flag.write("BSQR completed successfully.")
    logger.info("BQSR pipeline completed successfully.")


def run_bam_stats(sample_id, outdir):
    bqsr_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bam")
    bam_stats_file = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bamstats")
    bam_stats_flag = os.path.join(outdir, "bamstats.finish")

    # Verify the flag
    if os.path.exists(bam_stats_flag):
        logger.info(f"Stats result for {samid} already exist. Skip stats...")
        return

    # Run Samtools stats
    logger.info(f"Running Samtools stats for {samid}...")
    with open(bam_stats_file, "w") as stats_out:
        subprocess.run([SAMTOOLS, "stats", "-@", f"{PARAMETERS['threads']}", bqsr_bam], stdout=stats_out, check=True)
    logger.info("** bamstats done **")

    # Create finish flag
    with open(bam_stats_flag, "w") as flag:
        flag.write("bamstats completed successfully.")
    logger.info(f"Done statistics for {sample_id}.")


def run_bedtools(sample_id, outdir):
    bqsr_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bam")
    cvg_bed_gz = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz")
    finish_flag = os.path.join(outdir, "sorted_rmdup_realign_BQSR_cvg_bed.finish")

    # Step 0: Verify the flag
    if os.path.exists(finish_flag):
        logger.info(f"Bedtools result for {samid} already exist. Skip bedtools...")
        return

    # Step 1: Bedtools genome coverage
    logger.info("Running Bedtools genome coverage...")
    with open(cvg_bed_gz, "wb") as cvg_out:
        subprocess.run([BEDTOOLS, "genomecov", "-ibam", bqsr_bam, "-bga", "-split"], stdout=subprocess.PIPE, check=True)
        subprocess.run([BGZIP], stdin=subprocess.PIPE, stdout=cvg_out, check=True)
    logger.info("** sorted.rmdup.realign.BQSR.cvg.bed.gz done **")

    # Step 2: Index the compressed BED file
    logger.info("Indexing the compressed BED file with Tabix...")
    subprocess.run([TABIX, "-p", "bed", cvg_bed_gz], check=True)

    # Create finish flag
    with open(finish_flag, "w") as flag:
        flag.write("Bedtools pipeline completed successfully.")
    logger.info("Bedtools pipeline completed successfully.")


def move_final_output(fq, tmp_dir, final_dir):
    bam_list_file = bamlist_dir(fq)

    # Step 1: Move files to the final output directory, generate bam.list 
    logger.info("Moving final files to the output directory...")
    for file_suffix in [".bam", ".bam.bai", ".cvg.bed.gz", ".cvg.bed.gz.tbi"]:
        src_file = os.path.join(tmp_dir, f"{samid(fq)}.sorted.rmdup.realign.BQSR{file_suffix}")
        dst_file = os.path.join(final_dir, os.path.basename(src_file))
        if os.path.exists(src_file):
            os.rename(src_file, dst_file)
            if file_suffix == ".bam":
                with open(bam_list_file, "a") as bam_list:
                    bam_list.write(f"{dst_file}\n")

    # Step 2: Remove the temporary output directory
    logger.info("Removing temporary output directory...")
    shutil.rmtree(tmp_dir)
    logger.info(f"Temporary directory {tmp_dir} deleted.")


def run_alignment_pipeline(fq):
    os.makedirs(tmp_outdir(fq), exist_ok=True)
    os.makedirs(batch1_final_outdir(fq), exist_ok=True)
    finish_flag = os.path.join(batch1_final_outdir(fq), "alignment.finish")

    # Step 0: Verify the flag
    if os.path.exists(finish_flag):
        logger.info(f"Alignment result for {samid(fq)} already exist. Skip alignment pipeline...")
        return

    # Step 1: Run BWA to align and remove duplicates
    run_bwa_alignment(fq, tmp_outdir(fq))

    # Step 2: Realignment
    run_bwa_realign(samid(fq), tmp_outdir(fq))

    # Step 3: Recalibrate Base Quality Scores (BQSR)
    run_bqsr(samid(fq), tmp_outdir(fq))

    # Step 4: Generate BAM and coverage statistics
    run_bam_stats(samid(fq), tmp_outdir(fq))

    # Step 5: Bedtools
    run_bedtools(samid(fq), tmp_outdir(fq))

    # Step 6: Move final result file to batch1_final_files
    move_final_output(fq, tmp_outdir(fq), batch1_final_outdir(fq))

    # Create finish flag
    with open(finish_flag, "w") as flag:
        flag.write("Alignment pipeline completed successfully.")
    logger.info(f"Done alignment for {samid(fq)}")
