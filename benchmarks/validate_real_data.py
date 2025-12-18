#!/usr/bin/env python3
"""
TracEon Real-Data Validation Suite
Validates TracEon against actual genomic datasets from public repositories.

Features:
- Quick validation mode (small datasets)
- Medium validation mode (moderate datasets)
- Full validation mode (large datasets)
- Interactive CLI menu
- Sequential downloads (one file at a time)
"""

import os
import sys
import subprocess
import time
import gzip
import shutil
import requests
import argparse
import logging
from tqdm import tqdm
from typing import Optional, Dict, Any, List
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Safe unicode handling
def safe_unicode(char: str, fallback: str) -> str:
    """Return unicode char if supported, else fallback."""
    try:
        if sys.stdout.encoding and 'utf' in sys.stdout.encoding.lower():
            return char
    except:
        pass
    return fallback

CHECK = safe_unicode('✓', '[OK]')
CROSS = safe_unicode('✗', '[FAIL]')
INFO = safe_unicode('ℹ', '[INFO]')
WARN = safe_unicode('⚠', '[WARN]')
BULLET = safe_unicode('•', '*')

# Configuration
POSSIBLE_BINARIES = [
    Path("./build/traceon_driver"),
    Path("../build/traceon_driver"),
    Path("./traceon_driver"),
    Path("build/traceon_driver")
]
DATA_DIR = Path("real_data")
DEFAULT_LOOKUPS = 1000000
CHUNK_SIZE = 1024 * 1024 * 8  # 8 MB chunks
DOWNLOAD_TIMEOUT = 1200  # 20 minutes
MAX_RETRIES = 3

# Global session for connection reuse
_session = None

def get_session() -> requests.Session:
    """Get or create a global requests session."""
    global _session
    if _session is None:
        _session = requests.Session()
        adapter = requests.adapters.HTTPAdapter(
            pool_connections=10,
            pool_maxsize=10,
            max_retries=3,
            pool_block=False
        )
        _session.mount('http://', adapter)
        _session.mount('https://', adapter)
    return _session

# ============================================================================
# STANDARDIZED DATASETS (Organized by consistent size categories)
# ============================================================================

# QUICK: 1-10 MB each, fast download (<1 minute)
QUICK_DATASETS = [
    {
        "name": "Phi X 174 Phage",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz",
        "filename": "phix174.fasta.gz",
        "type": "fasta",
        "size": "~2 KB compressed",
        "records": "1 genome",
        "description": "Tiny phage genome (5.4 Kbp) - instant test",
        "full_size_mb": 0.002,
        "limit_lines": None,
    },
    {
        "name": "Lambda Phage",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/GCF_000840245.1_ViralProj14204_genomic.fna.gz",
        "filename": "lambda_phage.fasta.gz",
        "type": "fasta",
        "size": "~12 KB compressed",
        "records": "1 genome",
        "description": "Bacteriophage lambda (48.5 Kbp)",
        "full_size_mb": 0.012,
        "limit_lines": None,
    },
    {
        "name": "E. coli K-12 (Minimal)",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
        "filename": "ecoli_quick.fasta.gz",
        "type": "fasta",
        "size": "~1.4 MB compressed",
        "records": "1 chromosome",
        "description": "Complete E. coli genome (4.6 Mbp)",
        "full_size_mb": 1.4,
        "limit_lines": None,
    }
]

# MEDIUM: 10-100 MB each, moderate download (1-5 minutes)
MEDIUM_DATASETS = [
    {
        "name": "Yeast S. cerevisiae",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz",
        "filename": "yeast.fasta.gz",
        "type": "fasta",
        "size": "~12 MB compressed",
        "records": "16 chromosomes",
        "description": "Complete yeast genome (12 Mbp)",
        "full_size_mb": 12,
        "limit_lines": None,
    },
    {
        "name": "Drosophila melanogaster (Fruit Fly)",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz",
        "filename": "fruit_fly.fasta.gz",
        "type": "fasta",
        "size": "~38 MB compressed",
        "records": "6 chromosomes + MT",
        "description": "Fruit fly genome (180 Mbp)",
        "full_size_mb": 38,
        "limit_lines": None,
    },
    {
        "name": "Human Chromosome 21",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr21.fna.gz",
        "filename": "human_chr21.fasta.gz",
        "type": "fasta",
        "size": "~12 MB compressed",
        "records": "1 chromosome",
        "description": "Human chromosome 21 (46 Mbp)",
        "full_size_mb": 12,
        "limit_lines": None,
    },
    {
        "name": "Arabidopsis thaliana",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz",
        "filename": "arabidopsis.fasta.gz",
        "type": "fasta",
        "size": "~30 MB compressed",
        "records": "5 chromosomes",
        "description": "Model plant genome (135 Mbp)",
        "full_size_mb": 30,
        "limit_lines": None,
    }
]

# FULL: 100-1000 MB each, comprehensive download (5-20 minutes)
FULL_DATASETS = [
    {
        "name": "Human Chromosome 1",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fna.gz",
        "filename": "human_chr1.fasta.gz",
        "type": "fasta",
        "size": "~75 MB compressed",
        "records": "1 chromosome",
        "description": "Human chromosome 1 (248 Mbp)",
        "full_size_mb": 75,
        "limit_lines": None,
    },
    {
        "name": "Zebrafish Genome (Sample)",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz",
        "filename": "zebrafish_sample.fasta.gz",
        "type": "fasta",
        "size": "~500 MB (sampled)",
        "records": "25 chromosomes (sample)",
        "description": "Zebrafish genome sample (1.4 Gbp total)",
        "full_size_mb": 750,
        "limit_lines": 2000000,  # ~500 MB sample
    },
    {
        "name": "Mouse Genome (Sample)",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz",
        "filename": "mouse_sample.fasta.gz",
        "type": "fasta",
        "size": "~790 MB (sampled)",
        "records": "21 chromosomes (sample)",
        "description": "Mouse genome sample (2.7 Gbp total)",
        "full_size_mb": 850,
        "limit_lines": 1500000,  # ~400 MB sample
    },
]

# ADVANCED: 1GB+ datasets for benchmarking (optional)
ADVANCED_DATASETS = [
    {
        "name": "Human Chromosomes 1-3 (Combined)",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
        "filename": "human_chr1_3.fasta.gz",
        "type": "fasta",
        "size": "~300 MB compressed",
        "records": "3 chromosomes",
        "description": "Human chromosomes 1-3 (800 Mbp total)",
        "full_size_mb": 300,
        "chunks": [
            "chr1.fna.gz",
            "chr2.fna.gz",
            "chr3.fna.gz"
        ],
    }
]

# ANSI color codes
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def print_header(text: str):
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*70}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{text.center(70)}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*70}{Colors.ENDC}\n")

def print_success(text: str):
    print(f"{Colors.GREEN}{CHECK} {text}{Colors.ENDC}")

def print_warning(text: str):
    print(f"{Colors.YELLOW}{WARN} {text}{Colors.ENDC}")

def print_error(text: str):
    print(f"{Colors.RED}{CROSS} {text}{Colors.ENDC}")

def print_info(text: str):
    print(f"{Colors.CYAN}{INFO} {text}{Colors.ENDC}")

def calculate_total_size(datasets: List[Dict[str, Any]]) -> float:
    """Calculate total download size for a set of datasets."""
    total = 0
    for ds in datasets:
        if 'full_size_mb' in ds:
            if 'limit_lines' in ds and ds['limit_lines']:
                # Estimate actual download size for sampled datasets
                total += ds['limit_lines'] * 0.25 / 1000  # Conservative estimate
            else:
                total += ds['full_size_mb']
    return total

def find_binary() -> Optional[Path]:
    """Locate traceon_driver binary."""
    for p in POSSIBLE_BINARIES:
        if p.exists():
            return p.resolve()
    return None

def get_file_size_mb(filepath: Path) -> float:
    """Get file size in MB."""
    if filepath.exists():
        return filepath.stat().st_size / (1024 * 1024)
    return 0

def download_file(url: str, dest_path: Path, desc: str, 
                  timeout: int = DOWNLOAD_TIMEOUT,
                  max_retries: int = MAX_RETRIES) -> bool:
    """Download file with progress bar and retry logic."""
    session = get_session()
    
    for attempt in range(max_retries):
        try:
            response = session.get(url, stream=True, timeout=(30, timeout))
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(dest_path, "wb") as file, \
                 tqdm(desc=desc, total=total_size, unit='iB', 
                      unit_scale=True, unit_divisor=1024, miniters=1) as bar:
                for data in response.iter_content(chunk_size=CHUNK_SIZE):
                    size = file.write(data)
                    bar.update(size)
            return True
            
        except requests.exceptions.Timeout:
            if attempt == max_retries - 1:
                print_error(f"Download timeout after {max_retries} attempts!")
                return False
            print_warning(f"Timeout, retrying ({attempt + 1}/{max_retries})...")
            time.sleep(2 ** attempt)
            
        except requests.exceptions.RequestException as e:
            if attempt == max_retries - 1:
                print_error(f"Download failed: {e}")
                return False
            print_warning(f"Error: {e}, retrying ({attempt + 1}/{max_retries})...")
            time.sleep(2 ** attempt)
            
        except Exception as e:
            logger.error(f"Unexpected error during download: {e}")
            return False
    
    return False

def download_chunked_dataset(ds: Dict[str, Any]) -> Optional[Path]:
    """Download a dataset that consists of multiple chunks/parts (sequential)."""
    base_url = ds['url']
    chunks = ds.get('chunks', [])
    local_path = DATA_DIR / ds['filename']
    
    print_info(f"Downloading {len(chunks)} chunks sequentially for {ds['name']}...")
    
    temp_files = []
    
    # Download chunks one at a time
    for i, chunk_name in enumerate(chunks):
        chunk_url = base_url + chunk_name
        temp_file = DATA_DIR / f"temp_chunk_{i}.gz"
        temp_files.append(temp_file)
        
        print_info(f"Downloading chunk {i+1}/{len(chunks)}: {chunk_name}")
        if not download_file(chunk_url, temp_file, f"  Chunk {i+1}/{len(chunks)}"):
            print_error(f"Failed to download chunk {i+1}")
            # Cleanup partial downloads
            for tf in temp_files:
                if tf.exists():
                    tf.unlink()
            return None
    
    # Combine chunks
    print_info("Combining chunks into single file...")
    try:
        with open(local_path, 'wb') as outfile:
            for temp_file in temp_files:
                with open(temp_file, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)
    except Exception as e:
        print_error(f"Failed to combine chunks: {e}")
        return None
    finally:
        # Cleanup temp files
        for temp_file in temp_files:
            if temp_file.exists():
                try:
                    temp_file.unlink()
                except Exception as e:
                    logger.debug(f"Failed to remove temp file: {e}")
    
    return local_path

def prepare_dataset(ds: Dict[str, Any], skip_existing: bool = True) -> Optional[Path]:
    """Download and prepare a dataset."""
    local_path = DATA_DIR / ds['filename']
    
    # Skip if exists
    if skip_existing and local_path.exists():
        size_mb = get_file_size_mb(local_path)
        print_info(f"Found cached: {ds['filename']} ({size_mb:.2f} MB)")
        return local_path
    
    print(f"\n{Colors.BOLD}Processing: {ds['name']}{Colors.ENDC}")
    print(f"  Description: {ds['description']}")
    print(f"  Expected size: {ds['size']}")
    print(f"  Records: {ds['records']}")
    
    # Special handling for chunked datasets
    if 'chunks' in ds:
        return download_chunked_dataset(ds)
    
    # Handle limited downloads (for large files)
    if 'limit_lines' in ds and ds['limit_lines'] is not None:
        temp_gz = local_path.with_suffix('.temp.gz')
        print_info(f"Downloading (will extract {ds['limit_lines']:,} lines)...")
        
        if not download_file(ds['url'], temp_gz, "  Downloading"):
            if temp_gz.exists():
                temp_gz.unlink()
            return None
        
        print_info(f"Extracting first {ds['limit_lines']:,} lines...")
        try:
            with gzip.open(temp_gz, 'rt') as f_in, \
                 gzip.open(local_path, 'wt') as f_out:
                for i, line in enumerate(f_in):
                    if i >= ds['limit_lines']:
                        break
                    f_out.write(line)
                    if i > 0 and i % 10000 == 0:
                        print(f"    Processed {i:,} lines...", end='\r')
            print(f"\n    Processed {min(i+1, ds['limit_lines']):,} lines... Done!")
            print_success(f"Created limited dataset: {ds['filename']}")
        except Exception as e:
            print_error(f"Extraction failed: {e}")
            if local_path.exists():
                local_path.unlink()
            return None
        finally:
            if temp_gz.exists():
                try:
                    temp_gz.unlink()
                except Exception as e:
                    logger.debug(f"Failed to remove temp file: {e}")
    else:
        # Download complete file
        print_info("Downloading complete file...")
        if not download_file(ds['url'], local_path, "  Downloading"):
            return None
        print_success(f"Downloaded: {ds['filename']}")
    
    size_mb = get_file_size_mb(local_path)
    print_info(f"Final size: {size_mb:.2f} MB")
    return local_path

def run_traceon_test(binary: Path, filepath: Path, 
                     lookups: int = DEFAULT_LOOKUPS) -> Optional[Dict[str, Any]]:
    """Run TracEon load/save/lookup test on a dataset."""
    bin_path = filepath.with_suffix('.traceon')
    results = {}
    
    try:
        # Test 1: Load (with GZIP decompression)
        print(f"\n{Colors.CYAN}[1/3] Testing LOAD (GZIP decompression)...{Colors.ENDC}")
        cmd = [str(binary), "load", str(filepath)]
        t0 = time.time()
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, 
                                        text=True, timeout=300)
        load_time = time.time() - t0
        results['load_time'] = load_time
        print_success(f"Load completed in {load_time:.3f}s")
        
        # Test 2: Save (binary cache)
        print(f"\n{Colors.CYAN}[2/3] Testing SAVE (binary cache)...{Colors.ENDC}")
        cmd = [str(binary), "save", str(filepath), str(bin_path)]
        t0 = time.time()
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True, timeout=300)
        save_time = time.time() - t0
        results['save_time'] = save_time
        bin_size = get_file_size_mb(bin_path)
        results['binary_size_mb'] = bin_size
        print_success(f"Save completed in {save_time:.3f}s")
        print_info(f"Binary cache size: {bin_size:.2f} MB")
        
        # Test 3: Lookup (random access)
        print(f"\n{Colors.CYAN}[3/3] Testing LOOKUP ({lookups:,} random queries)...{Colors.ENDC}")
        cmd = [str(binary), "lookup", str(bin_path), str(lookups), "dummy", "0"]
        output = subprocess.check_output(cmd, text=True, timeout=300)
        throughput = 0
        for line in output.splitlines():
            if "Throughput:" in line:
                throughput_str = line.split(":")[1].strip().split()[0]
                throughput = float(throughput_str.replace(",", ""))
                break
        results['throughput'] = throughput
        print_success(f"Throughput: {throughput:,.0f} OPS/s")
        
        return results
        
    except subprocess.CalledProcessError as e:
        print_error(f"Command failed: {e.output if hasattr(e, 'output') else str(e)}")
        return None
    except subprocess.TimeoutExpired:
        print_error("Command timed out")
        return None
    except Exception as e:
        print_error(f"Test failed: {e}")
        return None
    finally:
        # Cleanup binary cache
        if bin_path.exists():
            try:
                bin_path.unlink()
            except Exception as e:
                logger.debug(f"Failed to remove binary: {e}")

def print_results_table(all_results: Dict[str, Optional[Dict[str, Any]]]):
    """Print formatted results table."""
    print_header("VALIDATION RESULTS SUMMARY")
    
    print(f"{'Dataset':<30} | {'Load (s)':<10} | {'Save (s)':<10} | "
          f"{'Throughput':<15} | {'Binary (MB)':<12}")
    print("-" * 95)
    
    for name, results in all_results.items():
        if results:
            print(f"{name:<30} | {results['load_time']:<10.3f} | "
                  f"{results['save_time']:<10.3f} | "
                  f"{results['throughput']:>10,.0f} OPS/s | "
                  f"{results['binary_size_mb']:>10.2f}")
        else:
            print(f"{name:<30} | {'FAILED':^10} | {'FAILED':^10} | "
                  f"{'FAILED':^15} | {'FAILED':^12}")
    
    print("-" * 95)

def interactive_menu() -> str:
    """Display interactive menu for dataset selection."""
    print_header("TracEon Real-Data Validation Suite v1.0.0")
    
    # Calculate sizes for each mode
    quick_size = calculate_total_size(QUICK_DATASETS)
    medium_size = calculate_total_size(MEDIUM_DATASETS)
    full_size = calculate_total_size(FULL_DATASETS)
    
    print(f"{Colors.BOLD}Select validation mode:{Colors.ENDC}\n")
    print(f"  {Colors.GREEN}[1] Quick Validation{Colors.ENDC}")
    print(f"      {BULLET} {len(QUICK_DATASETS)} tiny datasets ({quick_size:.1f} MB total)")
    print(f"      {BULLET} Perfect for initial testing")
    print(f"      {BULLET} Completes in ~1-3 minutes\n")
    
    print(f"  {Colors.CYAN}[2] Medium Validation{Colors.ENDC}")
    print(f"      {BULLET} {len(MEDIUM_DATASETS)} moderate datasets ({medium_size:.1f} MB total)")
    print(f"      {BULLET} Good balance of speed and coverage")
    print(f"      {BULLET} Completes in ~5-15 minutes\n")
    
    print(f"  {Colors.YELLOW}[3] Full Validation{Colors.ENDC}")
    print(f"      {BULLET} {len(FULL_DATASETS)} large datasets ({full_size:.1f} MB total)")
    print(f"      {BULLET} Comprehensive testing with varied data types")
    print(f"      {BULLET} Completes in ~20-40 minutes\n")
    
    print(f"  {Colors.BOLD}[4] Custom Selection{Colors.ENDC}")
    print(f"      {BULLET} Choose specific datasets from all categories\n")
    
    print(f"  {Colors.BOLD}[5] Advanced Benchmarking{Colors.ENDC}")
    print(f"      {BULLET} Very large datasets for performance testing\n")
    
    print(f"  {Colors.RED}[0] Exit{Colors.ENDC}\n")
    
    while True:
        try:
            choice = input(f"{Colors.BOLD}Enter choice [0-5]: {Colors.ENDC}").strip()
            if choice in ['0', '1', '2', '3', '4', '5']:
                return choice
            else:
                print_warning("Invalid choice. Please enter 0-5.")
        except KeyboardInterrupt:
            print("\n")
            return '0'

def custom_selection_menu() -> List[Dict[str, Any]]:
    """Allow user to select specific datasets."""
    print_header("Custom Dataset Selection")
    print(f"{Colors.BOLD}Available datasets (categorized by size):{Colors.ENDC}\n")
    
    categories = [
        ("QUICK (<10 MB)", QUICK_DATASETS),
        ("MEDIUM (10-100 MB)", MEDIUM_DATASETS),
        ("FULL (100+ MB)", FULL_DATASETS),
        ("ADVANCED (300+ MB)", ADVANCED_DATASETS)
    ]
    
    offset = 1
    dataset_map = {}
    
    for category_name, datasets in categories:
        if datasets:
            print(f"{Colors.BOLD}{category_name}:{Colors.ENDC}")
            for i, ds in enumerate(datasets, offset):
                print(f"  [{i}] {ds['name']}")
                print(f"      {ds['description']}")
                print(f"      Size: {ds['size']}, Records: {ds['records']}\n")
                dataset_map[i] = ds
            offset += len(datasets)
    
    print(f"  [0] Cancel\n")
    
    selected = []
    print("Enter dataset numbers separated by commas (e.g., 1,3,5):")
    try:
        choices = input(f"{Colors.BOLD}Selection: {Colors.ENDC}").strip()
        if choices == '0':
            return []
        
        indices = [int(x.strip()) for x in choices.split(',')]
        for idx in indices:
            if idx in dataset_map:
                selected.append(dataset_map[idx])
            else:
                print_warning(f"Skipping invalid index: {idx}")
    except (ValueError, KeyboardInterrupt):
        print_warning("\nInvalid input. Canceling selection.")
        return []
    
    return selected

def main():
    """Main execution function."""
    
    parser = argparse.ArgumentParser(
        description='TracEon Real-Data Validation Suite',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive menu
  python validate_real_data.py

  # Quick validation (non-interactive)
  python validate_real_data.py --quick

  # Medium validation
  python validate_real_data.py --medium

  # Full validation
  python validate_real_data.py --full

  # Skip download prompt
  python validate_real_data.py --quick --yes
        """
    )
    
    parser.add_argument('--quick', action='store_true',
                       help='Run quick validation (small datasets)')
    parser.add_argument('--medium', action='store_true',
                       help='Run medium validation')
    parser.add_argument('--full', action='store_true',
                       help='Run full validation (large datasets)')
    parser.add_argument('--advanced', action='store_true',
                       help='Run advanced validation (very large datasets)')
    parser.add_argument('--yes', '-y', action='store_true',
                       help='Skip confirmation prompts')
    parser.add_argument('--lookups', type=int, default=DEFAULT_LOOKUPS,
                       help=f'Number of random lookups (default: {DEFAULT_LOOKUPS:,})')
    parser.add_argument('--clean', action='store_true',
                       help='Remove cached datasets and exit')
    
    args = parser.parse_args()
    
    # Handle cleanup
    if args.clean:
        if DATA_DIR.exists():
            print_info(f"Removing cached datasets in {DATA_DIR}...")
            shutil.rmtree(DATA_DIR)
            print_success("Cleanup complete!")
        else:
            print_info("No cached data found.")
        return
    
    # Find binary
    binary = find_binary()
    if not binary:
        print_error("traceon_driver binary not found!")
        print_info("Please compile TracEon first:")
        print("  cmake -B build -DCMAKE_BUILD_TYPE=Release")
        print("  cmake --build build --target traceon_driver -j")
        sys.exit(1)
    
    print_success(f"Found binary: {binary}")
        # Create data directory
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # Determine which datasets to use
    selected_datasets = []
    mode_name = "Custom"

    if args.quick:
        selected_datasets = QUICK_DATASETS
        mode_name = "Quick"
    elif args.medium:
        selected_datasets = MEDIUM_DATASETS
        mode_name = "Medium"
    elif args.full:
        selected_datasets = FULL_DATASETS
        mode_name = "Full"
    elif args.advanced:
        selected_datasets = ADVANCED_DATASETS
        mode_name = "Advanced"
    else:
        # Interactive mode
        choice = interactive_menu()
        
        if choice == '0':
            print("Exiting...")
            return
        elif choice == '1':
            selected_datasets = QUICK_DATASETS
            mode_name = "Quick"
        elif choice == '2':
            selected_datasets = MEDIUM_DATASETS
            mode_name = "Medium"
        elif choice == '3':
            selected_datasets = FULL_DATASETS
            mode_name = "Full"
        elif choice == '4':
            selected_datasets = custom_selection_menu()
            if not selected_datasets:
                print_warning("No datasets selected. Exiting...")
                return
            mode_name = "Custom"
        elif choice == '5':
            selected_datasets = ADVANCED_DATASETS
            mode_name = "Advanced"

    # Calculate and display total size
    total_size_mb = calculate_total_size(selected_datasets)

    # Confirm execution
    if not args.yes:
        print_header(f"{mode_name} Validation Mode")
        print(f"Will test {len(selected_datasets)} datasets ({total_size_mb:.1f} MB total):")
        for ds in selected_datasets:
            category = ""
            if ds in QUICK_DATASETS:
                category = " (Quick)"
            elif ds in MEDIUM_DATASETS:
                category = " (Medium)"
            elif ds in FULL_DATASETS:
                category = " (Full)"
            elif ds in ADVANCED_DATASETS:
                category = " (Advanced)"
            
            print(f"  {BULLET} {ds['name']}{category}: {ds['size']}")
        print()
        
        confirm = input(f"{Colors.BOLD}Proceed? [y/N]: {Colors.ENDC}").strip().lower()
        if confirm not in ['y', 'yes']:
            print("Cancelled.")
            return

    # Run validation
    print_header(f"Starting {mode_name} Validation")
    print_info(f"Lookups per dataset: {args.lookups:,}")
    print_info(f"Total download size: {total_size_mb:.1f} MB")
    print_info(f"Data directory: {DATA_DIR.resolve()}")
    print_info(f"Download mode: Sequential (one file at a time)\n")

    all_results = {}

    for i, ds in enumerate(selected_datasets, 1):
        print(f"\n{Colors.BOLD}{'─'*70}{Colors.ENDC}")
        print(f"{Colors.BOLD}Dataset {i}/{len(selected_datasets)}: {ds['name']}{Colors.ENDC}")
        print(f"{Colors.BOLD}{'─'*70}{Colors.ENDC}")
        
        # Prepare dataset (download if needed)
        filepath = prepare_dataset(ds, skip_existing=True)
        if not filepath:
            print_error(f"Failed to prepare dataset: {ds['name']}")
            all_results[ds['name']] = None
            continue
        
        # Run tests
        results = run_traceon_test(binary, filepath, args.lookups)
        all_results[ds['name']] = results
        
        if results:
            print_success(f"Validation passed for {ds['name']}")
        else:
            print_error(f"Validation failed for {ds['name']}")

    # Print summary
    print_results_table(all_results)

    # Final verdict
    failed = sum(1 for r in all_results.values() if r is None)
    passed = len(all_results) - failed

    print()
    if failed == 0:
        print_success(f"{CHECK} ALL TESTS PASSED ({passed}/{len(all_results)})")
        print_info("TracEon is working correctly with real-world data!")
    else:
        print_warning(f"{WARN} {passed}/{len(all_results)} tests passed, {failed} failed")
        print_info("Review errors above for details.")

    print()
    print_info(f"Cached data location: {DATA_DIR.resolve()}")
    print_info(f"To clean up: python validate_real_data.py --clean")
    
    # Clean up session cache
    global _session_cache
    _session_cache.clear()

if __name__ == "__main__":
    main()