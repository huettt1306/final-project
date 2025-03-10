import subprocess
import os, re, shutil
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import bamlist_dir, glimpse_outdir, dbsnp_dir, glimpse_annot
from helper.path_define import filtered_vcf_path, filtered_tsv_path, chunks_path, norm_vcf_path, glimpse_vcf
from helper.logger import setup_logger
from concurrent.futures import ThreadPoolExecutor
from helper.file_utils import create_vcf_list, merge_vcf_list


# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "glimpse_pipeline.log"))


BCFTOOLS = TOOLS["bcftools"]
BGZIP = TOOLS["bgzip"]
TABIX = TOOLS["tabix"]
GLIMPSE_PHASE = TOOLS["GLIMPSE_phase"]
GLIMPSE_LIGATE = TOOLS["GLIMPSE_ligate"]
REF = PATHS["ref"]
MAP_PATH = PATHS["map_path"]

def compute_gls(fq, chromosome):
    glpath = os.path.join(glimpse_outdir(fq), "GL_file")
    bamlist = bamlist_dir(fq)
    with open(bamlist, "r") as bam_file:
        bam_lines = bam_file.readlines()

    for line in bam_lines:
        line = line.strip()
        filename = os.path.basename(line)
        name = filename.split(".")[0]

        output_vcf = os.path.join(glpath, f"{name}.{chromosome}.vcf.gz")
        command = [
            BCFTOOLS, "mpileup",
            "-f", REF, "-I", "-E", "-a", "FORMAT/DP",
            "-T", filtered_vcf_path(chromosome), "-r", chromosome, line, "-Ou"
        ]
        call_command = [
            BCFTOOLS, "call", "-Aim", "-C", "alleles",
            "-T", filtered_tsv_path(chromosome), "-Oz", "-o", output_vcf
        ]

        logger.info(f"Computing GL for sample {name}, chromosome {chromosome}")
        mpileup_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        call_process = subprocess.run(call_command, stdin=mpileup_process.stdout, capture_output=True, text=True)
        mpileup_process.stdout.close()

        if call_process.returncode != 0:
            logger.error(f"Error in GL computation for sample {name}: {call_process.stderr}")
            raise RuntimeError(f"Error in GL computation for sample {name}: {call_process.stderr}")

        index_command = [TABIX, "-f", output_vcf]
        subprocess.run(index_command, check=True)


def merge_gls(fq, chromosome):
    glpath = os.path.join(glimpse_outdir(fq), "GL_file")
    glmergepath = os.path.join(glimpse_outdir(fq), "GL_file_merged")
    gl_list_path = os.path.join(glpath, f"glimpse.{chromosome}_GL_list.txt")
    merged_vcf = os.path.join(glmergepath, f"glimpse.{chromosome}.vcf.gz")

    with open(gl_list_path, "w") as gl_list:
        for file in os.listdir(glpath):
            if file.endswith(f".{chromosome}.vcf.gz"):
                gl_list.write(os.path.join(glpath, file) + "\n")

    command = [
        BCFTOOLS, "merge", "-m", "none", 
        "-r", chromosome, "-Oz", "-o", merged_vcf, "-l", gl_list_path
    ]

    logger.info(f"Merging GL files for chromosome {chromosome}")
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error merging GLs for chromosome {chromosome}: {process.stderr}")
        raise RuntimeError(f"Error merging GLs for chromosome {chromosome}: {process.stderr}")

    index_command = [TABIX, "-f", merged_vcf]
    subprocess.run(index_command, check=True)

    logger.info(f"Merged GL file created at {merged_vcf}")

def phase_genome(fq, chromosome):
    glmergepath = os.path.join(glimpse_outdir(fq), "GL_file_merged")
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")

    map_file = os.path.join(MAP_PATH, f"{chromosome}.b38.gmap.gz")
    reference_vcf = norm_vcf_path(chromosome)
    chunk_file = chunks_path(chromosome)
    merged_vcf = os.path.join(glmergepath, f"glimpse.{chromosome}.vcf.gz")

    with open(chunk_file, "r") as chunks:
        for line in chunks:
            fields = line.strip().split()
            chunk_id = f"{int(fields[0]):02d}"
            input_region = fields[2]
            output_region = fields[3]
            output_vcf = os.path.join(imputed_path, f"glimpse.{chromosome}.{chunk_id}.imputed.vcf")

            command = [
                GLIMPSE_PHASE,
                "--input-gl", merged_vcf,
                "--reference", reference_vcf,
                "--map", map_file,
                "--input-region", input_region,
                "--output-region", output_region,
                "--output", output_vcf
            ]

            logger.info(f"Phasing chromosome {chromosome}, chunk {chunk_id}")
            process = subprocess.run(command, capture_output=True, text=True)
            if process.returncode != 0:
                logger.warning(f"Error phasing chromosome {chromosome}, chunk {chunk_id}: {process.stderr}")
                continue

            bgzip_command = [BGZIP, output_vcf]
            tabix_command = [TABIX, "-f", f"{output_vcf}.gz"]

            subprocess.run(bgzip_command, check=True)
            subprocess.run(tabix_command, check=True)

def extract_chunk_id(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    imputed_list = os.path.join(imputed_path, f"glimpse.{chromosome}_imputed_list.txt")

    # Hàm trích xuất chunk_id từ tên file
    def get_chunk_id(filename):
        match = re.search(rf"glimpse\.{chromosome}\.(\d+)\.imputed\.vcf\.gz", filename)
        return int(match.group(1)) if match else float('inf')

    # Lấy danh sách file phù hợp
    files = [
        file for file in os.listdir(imputed_path)
        if file.startswith(f"glimpse.{chromosome}.") and file.endswith(".imputed.vcf.gz")
    ]

    # Sắp xếp file theo chunk_id tăng dần
    sorted_files = sorted(files, key=get_chunk_id)

    # Ghi danh sách file đã sắp xếp vào imputed_list
    with open(imputed_list, "w") as imp_list:
        for file in sorted_files:
            imp_list.write(os.path.join(imputed_path, file) + "\n")


def ligate_genome(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    merged_path = os.path.join(glimpse_outdir(fq))
    os.makedirs(merged_path, exist_ok=True)

    imputed_list = os.path.join(imputed_path, f"glimpse.{chromosome}_imputed_list.txt")
    output_vcf = glimpse_vcf(fq, chromosome).replace(".gz", "")

    command = [
        GLIMPSE_LIGATE, "--input", imputed_list, "--output", output_vcf
    ]

    logger.info(f"Ligating genome for chromosome {chromosome}")
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error ligating chromosome {chromosome}: {process.stderr}")
        raise RuntimeError(f"Error ligating chromosome {chromosome}: {process.stderr}")

    bgzip_command = [BGZIP, output_vcf]
    tabix_command = [TABIX, "-f", f"{output_vcf}.gz"]

    subprocess.run(bgzip_command, check=True)
    subprocess.run(tabix_command, check=True)


def annotate(fq, chromosome):
    logger.info(f"Annotating genome for chromosome {chromosome}")
    
    # Đường dẫn đến file VCF đầu vào và đầu ra
    input = glimpse_vcf(fq, chromosome)
    output = glimpse_annot(fq, chromosome)
    
    # Bước 1: Thêm rsID từ dbSNP vào file VCF
    command_annotate = [
        BCFTOOLS, "annotate", "-a", f"{dbsnp_dir()}", "-c", "ID", "-O", "z", "-o", f"{output}", f"{input}"
    ]
    subprocess.run(command_annotate, check=True)
    logger.info(f"Added rsID annotations for chromosome {chromosome}")
    
    # Bước 2: Lọc các biến thể có ID bắt đầu bằng "rs" và lưu ra file mới
    filtered_output = output.replace(".vcf.gz", "_filtered.vcf")  # File tạm thời (chưa nén)
    
    # Sử dụng zgrep để lọc các dòng có ID bắt đầu bằng "rs" và thêm header
    with open(filtered_output, "w") as f:
        # Lấy header từ file VCF
        subprocess.run(["zgrep", "^#", output], stdout=f, check=True)
        # Lọc các dòng có ID bắt đầu bằng "rs"
        subprocess.run(["zgrep", "^[^#].*rs", output], stdout=f, check=True)
    
    # Nén file đã lọc
    subprocess.run([BGZIP, filtered_output], check=True)
    subprocess.run([TABIX, "-f", f"{filtered_output}.gz"], check=True)
    
    # Ghi đè file đầu ra bằng file đã lọc
    subprocess.run(["mv", f"{filtered_output}.gz", output], check=True)
    subprocess.run(["mv", f"{filtered_output}.gz.tbi", f"{output}.tbi"], check=True)
    
    logger.info(f"Saved filtered annotations to {output}")


def run_glimpse_chr(fq, chromosome):
    logger.info(f"Run glimpse for {fq} {chromosome}")

    if os.path.exists(glimpse_annot(fq, chromosome)):
        logger.info(f"Đã có kết quả glimpse cho mẫu {fq} với {chromosome}")
        return

    logger.info(f"Starting pipeline for chromosome {chromosome}...")

    # Step 1: Compute GLs
    compute_gls(fq, chromosome)

    # Step 2: Merge GLs
    merge_gls(fq, chromosome)

    # Step 3: Phase genome
    phase_genome(fq, chromosome)

    # Step 4: Ligate genome
    extract_chunk_id(fq, chromosome)
    ligate_genome(fq, chromosome)

    #step5: Annotate
    annotate(fq, chromosome)

    logger.info(f"Pipeline completed for chromosome {chromosome}.")


def run_glimpse(fq):    
    os.makedirs(os.path.join(glimpse_outdir(fq), "GL_file"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "GL_file_merged"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "imputed_file"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "imputed"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "annotated"), exist_ok=True)


    with ThreadPoolExecutor(max_workers=PARAMETERS["threads"]) as executor:
        executor.map(lambda chr: run_glimpse_chr(fq, chr), PARAMETERS["chrs"])

    
    shutil.rmtree(os.path.join(glimpse_outdir(fq), "GL_file"))
    shutil.rmtree(os.path.join(glimpse_outdir(fq), "GL_file_merged"))
    shutil.rmtree(os.path.join(glimpse_outdir(fq), "imputed_file"))
    logger.info(f"Deleted tmp dir for {fq}")

    imputed_list = create_vcf_list(os.path.join(glimpse_outdir(fq), "imputed"), "imputed")
    merge_vcf_list(imputed_list, os.path.join(glimpse_vcf(fq, "all")))

    annotated_list = create_vcf_list(os.path.join(glimpse_outdir(fq), "annotated"), "annotated")
    merge_vcf_list(annotated_list, os.path.join(glimpse_annot(fq, "all")))
    logger.info(f"Created merge vcf for all snp in {fq}")