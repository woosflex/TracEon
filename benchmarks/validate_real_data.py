import os
import subprocess
import time
import gzip
import shutil
import requests # Pip install requests
from tqdm import tqdm # Pip install tqdm

# Configuration
TRACEON_BINARY = "../build/traceon_driver"
DATA_DIR = "real_data"
LOOKUPS = 1000000

# Real Datasets (Public URLs)
DATASETS = [
    {
        "name": "E. coli Genome (FASTA)",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
        "filename": "ecoli.fasta",
        "type": "fasta"
    },
    {
        "name": "C. elegans WGS (FASTQ)",
        "url": "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR077/SRR077487/SRR077487_1.fastq.gz",
        "filename": "celegans_100k.fastq",
        "type": "fastq",
    }
]

def setup():
    if not os.path.exists(TRACEON_BINARY):
        print(f"Error: Binary {TRACEON_BINARY} not found. Please compile first.")
        exit(1)
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

def download_with_progress(url, dest_path, desc):
    """Download file with tqdm progress bar."""
    try:
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        block_size = 1024 * 1024 # 1MB

        with open(dest_path, "wb") as file, tqdm(
            desc=desc,
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for data in response.iter_content(block_size):
                size = file.write(data)
                bar.update(size)
        return True
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False

def download_and_extract(ds):
    local_path = os.path.join(DATA_DIR, ds['filename'])
    if os.path.exists(local_path):
        print(f"  [Skip] {ds['name']} already exists.")
        return local_path

    print(f"  [Proc] Processing {ds['name']}...")
    gz_path = local_path + ".gz"
    
    # 1. Download
    if not download_with_progress(ds['url'], gz_path, "    Downloading"):
        return None

    # 2. Extract
    print(f"    Extracting...", end=" ", flush=True)
    try:
        with gzip.open(gz_path, 'rb') as f_in:
            with open(local_path, 'wb') as f_out:
                if 'limit_lines' in ds:
                    # Stream limited lines for massive FASTQ files
                    # Note: Line counting is slow in python, so we just iterate
                    lines_to_read = ds['limit_lines']
                    pbar = tqdm(total=lines_to_read, desc="    Extracting (Lines)", unit=" lines")
                    for i, line in enumerate(f_in):
                        if i >= lines_to_read: break
                        f_out.write(line)
                        if i % 1000 == 0: pbar.update(1000)
                    pbar.close()
                else:
                    shutil.copyfileobj(f_in, f_out)
        print("Done.")
    except Exception as e:
        print(f"\nError extracting: {e}")
        return None
    finally:
        if os.path.exists(gz_path):
            os.remove(gz_path)
            
    return local_path

def run_traceon(filepath):
    bin_path = filepath + ".traceon"
    
    # 1. LOAD
    print("    Running LOAD...", end="", flush=True)
    cmd = [TRACEON_BINARY, "load", filepath]
    t0 = time.time()
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        print(f" Done ({time.time()-t0:.2f}s)")
    except subprocess.CalledProcessError as e:
        print(f"\n    LOAD FAILED: {e.output.decode()}")
        return

    # 2. SAVE
    print("    Running SAVE...", end="", flush=True)
    cmd = [TRACEON_BINARY, "save", filepath, bin_path]
    t0 = time.time()
    try:
        subprocess.check_output(cmd)
        print(f" Done ({time.time()-t0:.2f}s)")
    except subprocess.CalledProcessError as e:
        print(f" FAIL {e}")
        return

    # 3. LOOKUP (Restore + Query)
    print("    Running LOOKUP...", end="", flush=True)
    # Note: Using dummy args for prefix/max_index as driver auto-samples keys now
    cmd = [TRACEON_BINARY, "lookup", bin_path, str(LOOKUPS), "dummy", "0"]
    
    try:
        out = subprocess.check_output(cmd, text=True)
        # Parse output
        throughput = "0"
        for line in out.splitlines():
            if "Throughput:" in line:
                throughput = line.split(":")[1].strip().split()[0]
        print(f" Done! -> {throughput} OPS/s")
    except subprocess.CalledProcessError as e:
        print(f"\n    ERROR: Driver crashed! {e}")

def main():
    setup()
    print("=== TracEon Real-Data Validation ===\n")
    
    for ds in DATASETS:
        print(f"Dataset: {ds['name']}")
        fpath = download_and_extract(ds)
        if fpath:
            run_traceon(fpath)
        print("-" * 60)

if __name__ == "__main__":
    main()