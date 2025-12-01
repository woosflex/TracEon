import os
import time
import subprocess
import random
import sys
import threading
import psutil

# === CONFIGURATION ===
TRACEON_BINARY = "./build/traceon_driver"
SEQTK_BINARY = "seqtk"
SEQKIT_BINARY = "seqkit"
MAX_LOOKUPS = 500000
TARGET_LOOKUP_VOLUME_MB = 500

# === MATRIX CONFIGURATION ===
FILE_SIZES_MB = [10, 100]
# , 100, 500, 1000
DATA_TYPES = [
    ("WGS", "fastq", 150),
    ("WES", "fastq", 100),
    ("Panel", "fastq", 150),
    ("PacBio", "fastq", 15_000),
    ("Nanopore", "fastq", 30_000),
    ("RefGenome", "fasta", 1_000_000)
]

# === SHIM SCRIPTS ===
BP_PARSE_SCRIPT = "import sys; from Bio import SeqIO; [x for x in SeqIO.parse(sys.argv[1], sys.argv[2])]"

BP_LOOKUP_SCRIPT = """
import sys, random
from Bio import SeqIO
try:
    rd = SeqIO.to_dict(SeqIO.parse(sys.argv[1], sys.argv[2]))
    random.seed(42)
    for _ in range(int(sys.argv[3])):
        k = f"{sys.argv[4]}{random.randint(0, int(sys.argv[5])-1)}"
        _ = rd.get(k)
except MemoryError: sys.exit(137)
except Exception: sys.exit(1)
"""

# Updated shim: Crashes if 0 hits are found, ensuring the speed is real.
PYFASTX_SCRIPT = """
import sys, random, pyfastx
fname, n, prefix, m = sys.argv[1], int(sys.argv[3]), sys.argv[4], int(sys.argv[5])

if fname.endswith('.fastq'): db = pyfastx.Fastq(fname)
else: db = pyfastx.Fasta(fname)

random.seed(42)
hits = 0
for _ in range(n):
    try: 
        _ = str(db[f"{prefix}{random.randint(0, m-1)}"]) 
        hits += 1
    except: pass

# SANITY CHECK: If we found nothing, the benchmark is invalid.
if hits == 0:
    sys.stderr.write("Error: pyfastx found 0 records! Check ID format.\\n")
    sys.exit(1) # Return error code
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
        p.communicate()
        t.join()
        return p.returncode, time.time()-t0, stats['peak_mb']
    except: return -1, 0, 0

def parse_val(out, key):
    try: return float([l.split(":")[1] for l in out.splitlines() if key in l][0])
    except: return 0.0

def generate_dataset(fname, fmt, target_mb, avg_len):
    if fmt == "fastq":
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
    
    with open(fname, "w") as f:
        for i in range(num_seqs):
            s_len = avg_len
            if s_len > buf_size: 
                seq = "".join(random.choices("ACGT", k=s_len))
            else:
                start = random.randint(0, buf_size - s_len)
                seq = buf[start:start+s_len]
                
            if fmt == "fastq":
                f.write(f"@{prefix}{i}\n{seq}\n+\n{'I'*s_len}\n")
            else:
                f.write(f">{prefix}{i}\n{seq}\n")
    return num_seqs, prefix

def run_benchmark(label, fmt, size_mb, avg_len):
    fname, bin = f"bench_{label}_{size_mb}MB.{fmt}", f"bench_{label}_{size_mb}MB.bin"
    n_seqs, prefix = generate_dataset(fname, fmt, size_mb, avg_len)
    
    # Increase min lookups to ensure stable measurement
    iters = min(MAX_LOOKUPS, max(1000, int(TARGET_LOOKUP_VOLUME_MB*1024*1024/avg_len)))
    
    res = {}

    # 1. PARSING
    res['SeqTk'] = run_cmd([SEQTK_BINARY, "fqchk", fname]) if fmt == "fastq" else (0,0,0)
    res['SeqKit'] = run_cmd([SEQKIT_BINARY, "stats", fname])
    res['BioPy'] = run_cmd([sys.executable, "-c", BP_PARSE_SCRIPT, fname, fmt])
    res['TracEon'] = run_cmd([TRACEON_BINARY, "load", fname])

    # 2. LIFECYCLE
    run_cmd([TRACEON_BINARY, "save", fname, bin])
    p = subprocess.Popen([TRACEON_BINARY, "restore", bin], stdout=subprocess.PIPE, text=True)
    out, _ = p.communicate()
    res['Restore'] = parse_val(out, "Restore_Time_s")

    # 3. LOOKUPS
    res['BioPy_L'] = run_cmd([sys.executable, "-c", BP_LOOKUP_SCRIPT, fname, fmt, str(iters), prefix, str(n_seqs)])
    res['PyFastX_L'] = run_cmd([sys.executable, "-c", PYFASTX_SCRIPT, fname, str(iters), prefix, str(n_seqs)])

    p = subprocess.Popen([TRACEON_BINARY, "lookup", bin, str(iters), prefix, str(n_seqs)], stdout=subprocess.PIPE, text=True)
    out, _ = p.communicate()
    res['TracEon_L'] = parse_val(out, "Throughput")

    # Cleanup
    if os.path.exists(fname): os.remove(fname)
    if os.path.exists(bin): os.remove(bin)
    if os.path.exists(fname+".fxi"): os.remove(fname+".fxi")
    if os.path.exists(fname+".fai"): os.remove(fname+".fai")
    
    return res

def main():
    if not os.path.exists(TRACEON_BINARY):
        print("Compile driver first!")
        sys.exit(1)

    print(f"Starting Matrix Benchmark: {len(FILE_SIZES_MB)} Sizes x {len(DATA_TYPES)} Scenarios")

    for size_mb in FILE_SIZES_MB:
        print(f"\n" + "="*180)
        print(f"  FILE SIZE: {size_mb} MB")
        print(f"{'SCENARIO':<12} | {'SEQTK':<6} | {'SEQKIT':<6} | {'BIOPY':<6} | {'TRACEON':<7} | {'RESTORE':<7} | {'TRACEON OPS/s':<15} | {'PYFASTX OPS/s':<15} | {'RAM DIFF':<15}")
        print("-" * 180)
        
        for (label, fmt, avg_len) in DATA_TYPES:
            r = run_benchmark(label, fmt, size_mb, avg_len)
            
            stk = f"{r['SeqTk'][1]:.2f}" if r['SeqTk'][1] > 0 else "-"
            skt = f"{r['SeqKit'][1]:.2f}" if r['SeqKit'][1] > 0 else "-"
            bio = f"{r['BioPy'][1]:.2f}"
            tra = f"{r['TracEon'][1]:.2f}"
            rst = f"{r['Restore']:.2f}"
            
            ops_tra = f"{int(r['TracEon_L']):,}"
            
            iters = min(MAX_LOOKUPS, max(1000, int(TARGET_LOOKUP_VOLUME_MB*1024*1024/avg_len)))
            pfx_dur = r['PyFastX_L'][1]
            
            if pfx_dur > 0:
                ops_pfx_val = iters/pfx_dur
                ops_pfx = f"{int(ops_pfx_val):,}"
                if ops_pfx_val > 10_000_000: ops_pfx += " (?)"
            else:
                ops_pfx = "-"
            
            bp_rc, _, bp_ram_val = r['BioPy_L']
            tra_ram_val = r['TracEon'][2]
            
            if bp_rc != 0: ram_str = f"{tra_ram_val:.0f}MB (BioPy CRASH)"
            elif bp_ram_val > 0:
                diff_pct = (1 - tra_ram_val / bp_ram_val) * 100
                ram_str = f"{tra_ram_val:.0f}MB (-{diff_pct:.0f}%)"
            else: ram_str = f"{tra_ram_val:.0f}MB (?)"
            
            print(f"{label:<12} | {stk:<6} | {skt:<6} | {bio:<6} | {tra:<7} | {rst:<7} | {ops_tra:<15} | {ops_pfx:<15} | {ram_str:<15}")
            
    print("="*180)

if __name__ == "__main__":
    main()