import subprocess

def get_md5sum2(fastq_fp, sample):
    proc = subprocess.run(
        ["md5sum", fastq_fp], check=True, capture_output=True, encoding="utf-8"
    )
    res = sample, proc.stdout
    return res
