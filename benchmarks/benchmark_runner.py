import os
import time
import subprocess
import random
import shutil
import sys
import threading
from pathlib import Path
import psutil
import gzip

# === CONFIGURATION ===
REPO_ROOT = Path(__file__).resolve().parent.parent
TRACEON_BINARY = str(REPO_ROOT / "build" / "traceon_driver")
SEQTK_BINARY = "seqtk"
SEQKIT_BINARY = "seqkit"
MAX_LOOKUPS = 500000
TARGET_LOOKUP_VOLUME_MB = 500

# === MATRIX CONFIGURATION ===
FILE_SIZES_MB = [10, 100, 500, 1000]

DATA_TYPES = [
    ("WGS", "fastq", 150),
    ("WGS_GZ", "fastq.gz", 150),
    ("PacBio", "fastq", 15_000),
    ("PacBio_GZ", "fastq.gz", 30_000),
    ("Nanopore", "fastq", 30_000),
    ("Nanopore_GZ", "fastq.gz", 30_000),
    ("RefGenome", "fasta", 1_000_000),
    ("RefGenome_GZ", "fasta.gz", 1_000_000)
]

# === SHIM SCRIPTS ===
BP_PARSE_SCRIPT = "import sys, gzip; from Bio import SeqIO; open_func = gzip.open if sys.argv[1].endswith('.gz') else open; [x for x in SeqIO.parse(open_func(sys.argv[1], 'rt'), sys.argv[2])]"

PYFASTX_SCRIPT = """
import sys, random, pyfastx
fname, n, prefix, m = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])

# pyfastx handles gzip automatically
if 'fastq' in fname: db = pyfastx.Fastq(fname)
else: db = pyfastx.Fasta(fname)

random.seed(42)
hits = 0
for _ in range(n):
    try: 
        _ = str(db[f"{prefix}{random.randint(0, m-1)}"]) 
        hits += 1
    except: pass

if hits == 0:
    sys.stderr.write("Error: pyfastx found 0 records!\\n")
    sys.exit(1)
"""

# === UTILITIES ===
def monitor_process(proc, res):
    peak = 0
    try:
        p = psutil.Process(proc.pid)
        while proc.poll() is None:
            try: peak = max(peak, p.memory_info().rss)
            except: break
            time.sleep(0.01)
    except: pass
    res['peak_mb'] = peak / (1024*1024)

def run_cmd(cmd):
    stats = {'peak_mb': 0}
    t0 = time.time()
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        t = threading.Thread(target=monitor_process, args=(p, stats))
        t.start()
        _, stderr = p.communicate()
        t.join()
        if p.returncode != 0:
            print(f"Command failed ({p.returncode}): {' '.join(cmd)}", file=sys.stderr)
            if stderr:
                print(stderr.strip(), file=sys.stderr)
        return p.returncode, time.time()-t0, stats['peak_mb']
    except OSError as exc:
        print(f"Could not run {' '.join(cmd)}: {exc}", file=sys.stderr)
        return -1, 0, 0

def run_required(cmd):
    result = run_cmd(cmd)
    if result[0] != 0:
        raise RuntimeError(f"Benchmark command failed: {' '.join(cmd)}")
    return result

def parse_val(out, key):
    try: return float([l.split(":")[1] for l in out.splitlines() if key in l][0])
    except: return 0.0

def generate_dataset(fname, fmt, target_mb, avg_len):
    if "fastq" in fmt:
        bytes_per_rec = (avg_len * 2) + 50
        prefix = "r"
    else:
        bytes_per_rec = avg_len + 20
        prefix = "c"
        
    num_seqs = int((target_mb * 1024 * 1024) / bytes_per_rec)
    if num_seqs < 1: num_seqs = 1

    print(f"  Generating {fname} ({num_seqs:,} recs, ~{avg_len}bp)...")
    
    buf_size = 10_000_000
    buf = "".join(random.choices("ACGT", k=buf_size))
    
    opener = gzip.open if fname.endswith(".gz") else open
    mode = "wt" if fname.endswith(".gz") else "w"

    with opener(fname, mode) as f:
        for i in range(num_seqs):
            s_len = avg_len
            if s_len > buf_size: 
                seq = "".join(random.choices("ACGT", k=s_len))
            else:
                start = random.randint(0, buf_size - s_len)
                seq = buf[start:start+s_len]
                
            if "fastq" in fmt:
                f.write(f"@{prefix}{i}\n{seq}\n+\n{'I'*s_len}\n")
            else:
                f.write(f">{prefix}{i}\n{seq}\n")
    return num_seqs, prefix

def run_benchmark(label, fmt, size_mb, avg_len):
    fname = f"bench_{label}_{size_mb}MB.{fmt}"
    bin_file = f"bench_{label}_{size_mb}MB.traceon"
    
    n_seqs, prefix = generate_dataset(fname, fmt, size_mb, avg_len)
    actual_size_mb = os.path.getsize(fname) / (1024 * 1024)
    
    MIN_ITERS = 100_000
    
    calc_iters = int(TARGET_LOOKUP_VOLUME_MB * 1024 * 1024 / avg_len)
    iters = min(MAX_LOOKUPS, max(MIN_ITERS, calc_iters))
    
    res = {}
    res['FileSize'] = actual_size_mb
    res['Iterations'] = iters

    # 1. PARSING
    res['SeqTk'] = run_required([SEQTK_BINARY, "fqchk", fname]) if "fastq" in fmt else (0,0,0)
    res['SeqKit'] = run_required([SEQKIT_BINARY, "stats", fname])
    
    clean_fmt = "fastq" if "fastq" in fmt else "fasta"
    res['BioPy'] = run_required([sys.executable, "-c", BP_PARSE_SCRIPT, fname, clean_fmt])
    res['TracEon'] = run_required([TRACEON_BINARY, "load", fname])

    # 2. LIFECYCLE
    run_required([TRACEON_BINARY, "save", fname, bin_file])
    p = subprocess.run(
        [TRACEON_BINARY, "restore", bin_file],
        check=True,
        capture_output=True,
        text=True,
    )
    out = p.stdout
    res['Restore'] = parse_val(out, "Restore_Time_s")

    # 3. LOOKUPS
    # NOTE: BioPython lookup test removed as per instructions
    res['PyFastX_L'] = run_required([sys.executable, "-c", PYFASTX_SCRIPT, fname, str(iters), prefix, str(n_seqs)])

    p = subprocess.run(
        [TRACEON_BINARY, "lookup", bin_file, str(iters), prefix, str(n_seqs)],
        check=True,
        capture_output=True,
        text=True,
    )
    out = p.stdout
    res['TracEon_L'] = parse_val(out, "Throughput")
    res['Record_Count'] = n_seqs

    # Cleanup
    if os.path.exists(fname): os.remove(fname)
    if os.path.exists(bin_file): os.remove(bin_file)
    if os.path.exists(fname+".fxi"): os.remove(fname+".fxi")
    if os.path.exists(fname+".fai"): os.remove(fname+".fai")
    
    return res

def main():
    global FILE_SIZES_MB, DATA_TYPES
    if "--quick" in sys.argv:
        # Reduced matrix for fast sanity checks (replaces quick_benchmark.py).
        FILE_SIZES_MB = [10]
        DATA_TYPES = DATA_TYPES[:3]
    if not os.path.exists(TRACEON_BINARY):
        print("Compile driver first!")
        sys.exit(1)

    missing = [tool for tool in (SEQTK_BINARY, SEQKIT_BINARY) if not shutil.which(tool)]
    if missing:
        print(f"Missing benchmark tools on PATH: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    random.seed(42)

    print(f"Starting Matrix Benchmark: {len(FILE_SIZES_MB)} Sizes x {len(DATA_TYPES)} Scenarios")

    for size_mb in FILE_SIZES_MB:
        print(f"\n" + "="*200)
        print(f"  FILE SIZE: {size_mb} MB")
        print(f"{'SCENARIO':<12} | {'RECS':<9} | {'SEQTK':<6} | {'SEQKIT':<6} | {'BIOPY':<6} | {'TRACEON':<7} | {'RESTORE':<7} | {'TRACEON OPS':<12} | {'PYFASTX OPS':<12} | {'MEM TRA':<8} | {'MEM PFX':<8} | {'vs FILE':<8}")
        print("-" * 200)
        
        for (label, fmt, avg_len) in DATA_TYPES:
            r = run_benchmark(label, fmt, size_mb, avg_len)
            
            recs = f"{r['Record_Count']:,}"
            stk = f"{r['SeqTk'][1]:.2f}" if r['SeqTk'][1] > 0 else "-"
            skt = f"{r['SeqKit'][1]:.2f}" if r['SeqKit'][1] > 0 else "-"
            bio = f"{r['BioPy'][1]:.2f}"
            tra = f"{r['TracEon'][1]:.2f}"
            rst = f"{r['Restore']:.2f}"
            
            ops_tra = f"{int(r['TracEon_L']):,}"
            
            iters = r['Iterations']
            pfx_dur = r['PyFastX_L'][1]
            
            if pfx_dur > 0:
                ops_pfx_val = iters/pfx_dur
                ops_pfx = f"{int(ops_pfx_val):,}"
            else:
                ops_pfx = "-"
            
            # Memory Stats
            tra_mem = r['TracEon'][2]
            pfx_mem = r['PyFastX_L'][2]
            file_sz = r['FileSize']
            
            tra_mem_str = f"{tra_mem:.0f}"
            pfx_mem_str = f"{pfx_mem:.0f}"
            
            # Compare TracEon RAM vs File Size (overhead)
            # Shows how much larger/smaller the loaded structure is vs disk
            if file_sz > 0:
                overhead_pct = ((tra_mem / file_sz) - 1) * 100
                vs_file = f"{overhead_pct:+.0f}%"
            else:
                vs_file = "?"
            
            print(f"{label:<12} | {recs:<9} | {stk:<6} | {skt:<6} | {bio:<6} | {tra:<7} | {rst:<7} | {ops_tra:<12} | {ops_pfx:<12} | {tra_mem_str:<8} | {pfx_mem_str:<8} | {vs_file:<8}")
            
    print("="*200)

if __name__ == "__main__":
    main()
