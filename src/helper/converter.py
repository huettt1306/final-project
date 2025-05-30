import subprocess, os
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.logger import setup_logger

logger = setup_logger(os.path.join(PATHS["logs"], "conversion.log"))

def convert_cram_to_fastq(cram_path, output_fastq_path_1, output_fastq_path_2):
    logger.info(f"Converting CRAM file {cram_path} to FASTQ.GZ  using 8 threads...")
    command = f"{TOOLS['samtools']} fastq -@ {PARAMETERS['threads']} -1 {output_fastq_path_1} -2 {output_fastq_path_2} {cram_path}"
    subprocess.run(command, shell=True, check=True)


def convert_haploid_to_diploid(vcf_path):
    tmp1 = vcf_path.replace(".vcf.gz", ".tmp1.vcf.gz")
    tmp2 = vcf_path.replace(".vcf.gz", ".tmp2.vcf.gz")
    
    subprocess.run([TOOLS["bcftools"], "+setGT", 
        "-Oz", "-o", tmp1, vcf_path,
        "--", "-t", "q", "-i", 'GT=="0"', "-n", 'c:0/0'
    ], check=True)
    
    subprocess.run([TOOLS["bcftools"], "+setGT", tmp1,
        "-Oz", "-o", tmp2, tmp1,
        "--", "-t", "q", "-i", 'GT=="1"', "-n", 'c:1/1'
    ], check=True)
    
    os.replace(tmp2, vcf_path)
    
    if os.path.exists(tmp1):
        os.remove(tmp1)


def convert_diploid_to_haploid(vcf_path):
    tmp_vcf = vcf_path.replace(".vcf.gz", ".tmp.vcf")

    bcftools_process = subprocess.Popen([TOOLS["bcftools"], "view", vcf_path
        ], stdout=subprocess.PIPE,  text=True  )

    sed_process = subprocess.Popen(
        ["sed", "-E", r's/(\t[01])\/[01]/\1/g'],
        stdin=bcftools_process.stdout,  
        stdout=subprocess.PIPE,  
        text=True 
    )

    with open(tmp_vcf, "w") as output_file:
        for line in sed_process.stdout:
            output_file.write(line)

    bcftools_process.wait()
    sed_process.wait()

    subprocess.run([TOOLS["bgzip"], tmp_vcf], check=True)

    os.replace(f"{tmp_vcf}.gz", vcf_path)


