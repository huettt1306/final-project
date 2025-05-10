import subprocess
import os, re, shutil
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import bamlist_dir, glimpse_outdir, dbsnp_dir, glimpse_annot, samid
from helper.path_define import filtered_vcf_path, filtered_tsv_path, chunks_path, norm_vcf_path, glimpse_vcf
from helper.logger import setup_logger
from helper.converter import convert_haploid_to_diploid, convert_diploid_to_haploid
from concurrent.futures import ThreadPoolExecutor
from helper.file_utils import create_vcf_list, merge_vcf_list


logger = setup_logger(os.path.join(PATHS["logs"], "glimpse_pipeline.log"))


BCFTOOLS = TOOLS["bcftools"]
BGZIP = TOOLS["bgzip"]
TABIX = TOOLS["tabix"]
GLIMPSE_PHASE = TOOLS["GLIMPSE_phase"]
GLIMPSE_LIGATE = TOOLS["GLIMPSE_ligate"]
REF = PATHS["ref"]
MAP_PATH = PATHS["map_path"]

def create_samples_file_arg(fq):
    samples_file = os.path.join(os.path.join(glimpse_outdir(fq), "imputed_file"), "samples.txt")
    with open(samples_file, "w") as sf:
        sf.write(f"{samid(fq)} {PARAMETERS['gender']}\n")


def compute_gls(fq, chromosome):
    glpath = os.path.join(glimpse_outdir(fq), "GL_file")
    bamlist = bamlist_dir(fq)
    with open(bamlist, "r") as bam_file:
        bam_lines = bam_file.readlines()

    for line in bam_lines:
        line = line.strip()
        name = os.path.basename(line).split(".")[0]
        output_vcf = os.path.join(glpath, f"{name}.{chromosome}.vcf.gz")
    
        try:
            print(f"Computing GL for sample {name}, chromosome {chromosome}")
            
            mpileup_process = subprocess.Popen([BCFTOOLS, "mpileup",
                "-f", REF, "-I", "-E", "-a", "FORMAT/DP",
                "-T", filtered_vcf_path(chromosome), "-r", chromosome, line, "-Ou"
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            call_process = subprocess.run([BCFTOOLS, "call", "-Aim", "-C", "alleles",
                "-T", filtered_tsv_path(chromosome), "-Oz", "-o", output_vcf
            ], stdin=mpileup_process.stdout, capture_output=True, text=True, check=True)
            
            mpileup_process.stdout.close()
            
            # Check if mpileup had any errors
            _, mpileup_stderr = mpileup_process.communicate()
            if mpileup_process.returncode != 0:
                print(f"Error in mpileup for sample {name}, chromosome {chromosome}:")
                print(mpileup_stderr.decode())
                
        except subprocess.CalledProcessError as e:
            print(f"Error in bcftools call for sample {name}, chromosome {chromosome}:")
            print(f"Command: {e.cmd}")
            print(f"Return code: {e.returncode}")
            print(f"Output: {e.stdout}")
            print(f"Error: {e.stderr}")
        except Exception as e:
            print(f"Unexpected error processing sample {name}, chromosome {chromosome}:")
            print(str(e))

        subprocess.run([TABIX, "-f", output_vcf], check=True)


def merge_gls(fq, chromosome):
    gl_list_path = create_vcf_list(os.path.join(glimpse_outdir(fq), "GL_file"), f"{chromosome}") 
    merged_vcf = os.path.join(os.path.join(glimpse_outdir(fq), "GL_file_merged"), f"{chromosome}.vcf.gz")

    print(f"Merging GL files for chromosome {chromosome}")
    subprocess.run([BCFTOOLS, "merge", 
        "-m", "none", "-r", chromosome, "-Oz", "-o", merged_vcf, "-l", gl_list_path
    ], capture_output=True, text=True, check=True)
    
    subprocess.run([TABIX, "-f", merged_vcf], check=True)

    print(f"Merged GL file created at {merged_vcf}")


def phase_genome(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    map_file = os.path.join(MAP_PATH, f"{chromosome}.b38.gmap.gz")
    reference_vcf = norm_vcf_path(chromosome)
    chunk_file = chunks_path(chromosome)
    merged_vcf = os.path.join(os.path.join(glimpse_outdir(fq), "GL_file_merged"), f"{chromosome}.vcf.gz")
    samples_file = os.path.join(imputed_path, "samples.txt")

    with open(chunk_file, "r") as chunks:
        for line in chunks:
            fields = line.strip().split()
            chunk_id = f"{int(fields[0]):02d}"
            input_region = fields[2]
            output_region = fields[3]
            output_vcf = os.path.join(imputed_path, f"{chromosome}.{chunk_id}.imputed.vcf")

            print(f"Phasing chromosome {chromosome}, chunk {chunk_id}")

            command = [GLIMPSE_PHASE,
                "--input-gl", merged_vcf,
                "--reference", reference_vcf,
                "--map", map_file,
                "--input-region", input_region,
                "--output-region", output_region,
                "--output", output_vcf
            ]
            if chromosome == "chrX":
                command.extend(["--samples-file", samples_file])

            process = subprocess.run(command, capture_output=True, text=True)
            if process.returncode != 0:
                logger.warning(f"Error phasing chromosome {chromosome}, chunk {chunk_id}: {process.stderr}")
                continue

            if os.path.exists(f"{output_vcf}.gz"):
                os.remove(f"{output_vcf}.gz")

            subprocess.run([BGZIP, output_vcf], check=True)
            if chromosome == "chrX" and PARAMETERS['gender'] == 1:
                convert_haploid_to_diploid(f"{output_vcf}.gz")
            subprocess.run([TABIX, "-f", f"{output_vcf}.gz"], check=True)


def extract_chunk_id(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    imputed_list = os.path.join(imputed_path, f"{chromosome}_imputed_list.txt")

    def get_chunk_id(filename):
        match = re.search(rf"{chromosome}\.(\d+)\.imputed\.vcf\.gz", filename)
        return int(match.group(1)) if match else float('inf')

    files = [
        file for file in os.listdir(imputed_path)
        if file.startswith(f"{chromosome}.") and file.endswith(".imputed.vcf.gz")
    ]

    sorted_files = sorted(files, key=get_chunk_id)

    with open(imputed_list, "w") as imp_list:
        for file in sorted_files:
            imp_list.write(os.path.join(imputed_path, file) + "\n")


def ligate_genome(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    imputed_list = os.path.join(imputed_path, f"{chromosome}_imputed_list.txt")
    output_vcf = glimpse_vcf(fq, chromosome).replace(".gz", "")

    print(f"Ligating genome for chromosome {chromosome}")
    subprocess.run([GLIMPSE_LIGATE, 
        "--input", imputed_list, "--output", output_vcf
    ], capture_output=True, text=True, check=True)

    if os.path.exists(f"{output_vcf}.gz"):
        os.remove(f"{output_vcf}.gz")

    subprocess.run([BGZIP, output_vcf], check=True)
    if chromosome == "chrX" and PARAMETERS['gender'] == 1:
        convert_diploid_to_haploid(f"{output_vcf}.gz")
    subprocess.run([TABIX, "-f", f"{output_vcf}.gz"], check=True)


def annotate(fq, chromosome):
    print(f"Annotating genome for chromosome {chromosome}")
    
    input = glimpse_vcf(fq, chromosome)
    output = glimpse_annot(fq, chromosome)
    
    subprocess.run([BCFTOOLS, "annotate", 
        "-a", f"{dbsnp_dir()}", "-c", "ID", "-O", "z", "-o", f"{output}", f"{input}"
    ], check=True)
    print(f"Added rsID annotations for chromosome {chromosome}")
    
    filtered_output = output.replace(".vcf.gz", "_filtered.vcf") 
    
    with open(filtered_output, "w") as f:
        subprocess.run(["zgrep", "^#", output], stdout=f, check=True)
        subprocess.run(["zgrep", "^[^#].*rs", output], stdout=f, check=True)
    
    if os.path.exists(f"{filtered_output}.gz"):
        os.remove(f"{filtered_output}.gz")

    subprocess.run([BGZIP, filtered_output], check=True)
    subprocess.run([TABIX, "-f", f"{filtered_output}.gz"], check=True)
    
    subprocess.run(["mv", f"{filtered_output}.gz", output], check=True)
    subprocess.run(["mv", f"{filtered_output}.gz.tbi", f"{output}.tbi"], check=True)
    
    print(f"Saved filtered annotations to {output}")


def run_glimpse_chr(fq, chromosome):
    finish_flag = os.path.join(f"{glimpse_outdir(fq)}", "imputed", f"glimpse_{chromosome}.finish")

    # Step 0: Verify the flag
    if os.path.exists(finish_flag):
        print(f"Glimpse result for {chromosome} {fq} already exist. Skip glimpse for {chromosome}...")
        #return

    print(f"Starting glimpse for chromosome {chromosome}...")

    # Step 1: Compute GLs
    compute_gls(fq, chromosome)

    # Step 2: Merge GLs
    merge_gls(fq, chromosome)

    # Step 3: Phase genome
    phase_genome(fq, chromosome)

    # Step 4: Ligate genome
    extract_chunk_id(fq, chromosome)
    ligate_genome(fq, chromosome)

    # Step 5: Annotate
    annotate(fq, chromosome)

    # Create finish flag
    with open(finish_flag, "w") as flag:
        flag.write(f"Completed Glimpse for {fq} {chromosome}.")
    print(f"Completed Glimpse for {fq} chromosome {chromosome}")


def run_glimpse(fq):    
    os.makedirs(os.path.join(glimpse_outdir(fq), "GL_file"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "GL_file_merged"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "imputed_file"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "imputed"), exist_ok=True)
    os.makedirs(os.path.join(glimpse_outdir(fq), "annotated"), exist_ok=True)
    finish_flag = os.path.join(glimpse_outdir(fq), "glimpse.finish")

    # Verify the flag
    if os.path.exists(finish_flag):
        print(f"Glimpse result for {fq} already exist. Skip glimpse step...")
        #return
    
    create_samples_file_arg(fq)
    with ThreadPoolExecutor(max_workers=PARAMETERS['threads']) as executor:
        executor.map(lambda chr: run_glimpse_chr(fq, chr), PARAMETERS["chrs"])
    
    shutil.rmtree(os.path.join(glimpse_outdir(fq), "GL_file"))
    shutil.rmtree(os.path.join(glimpse_outdir(fq), "GL_file_merged"))
    shutil.rmtree(os.path.join(glimpse_outdir(fq), "imputed_file"))
    print(f"Deleted tmp dir for {fq}")

    return

    imputed_list = create_vcf_list(os.path.join(glimpse_outdir(fq), "imputed"), "imputed")
    merge_vcf_list(imputed_list, os.path.join(glimpse_vcf(fq, "all")))

    annotated_list = create_vcf_list(os.path.join(glimpse_outdir(fq), "annotated"), "annotated")
    merge_vcf_list(annotated_list, os.path.join(glimpse_annot(fq, "all")))
    print(f"Created merge vcf for all snp in {fq}")

    with open(finish_flag, "w") as flag:
        flag.write("Glimpse pipeline completed successfully.")
    print(f"Glimpse pipeline for {fq} completed successfully.")

