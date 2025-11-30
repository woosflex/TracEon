import os
import time
import subprocess
import random
import sys

# === CONFIGURATION ===
TRACEON_BINARY = "./build/traceon_driver"
SEQTK_BINARY = "seqtk"         # Ensure installed
TARGET_DATA_VOLUME_MB = 500    # Target 500MB data movement for lookups (keeps runtime sane)
MAX_LOOKUPS_CAP = 1_000_000    # Cap lookups to prevent infinite loops

# Try importing BioPython
try:
    from Bio import SeqIO
except ImportError:
    print("Error: BioPython is not installed. Please run: pip install biopython")
    sys.exit(1)

# === BIOLOGICAL SCENARIOS ===
SCENARIOS = [
    {
        "name": "WGS (Whole Genome - Illumina)",
        "type": "fastq",
        "prefix": "cluster_",
        "num_seqs": 1_000_000, # 1 Million reads
        "avg_len": 150,
        "desc": "High throughput, massive file size"
    },
    {
        "name": "WES (Whole Exome - Illumina)",
        "type": "fastq",
        "prefix": "exome_",
        "num_seqs": 250_000,
        "avg_len": 100,
        "desc": "High depth, shorter reads"
    },
    {
        "name": "TGE (Targeted Gene Panel)",
        "type": "fastq",
        "prefix": "panel_",
        "num_seqs": 50_000,
        "avg_len": 150,
        "desc": "Low volume, high precision"
    },
    {
        "name": "Long Reads (PacBio HiFi/Nanopore)",
        "type": "fastq",
        "prefix": "read_",
        "num_seqs": 10_000,
        "avg_len": 15_000, # 15kb reads
        "desc": "Very long sequences, memory bandwidth heavy"
    },
    {
        "name": "Reference Genome (Alignment/Mapping)",
        "type": "fasta",
        "prefix": "chr",
        "num_seqs": 24,    # Human (1-22, X, Y)
        "avg_len": 10_000_000, # Simulating 10MB chunks for test speed
        "desc": "Massive contiguous strings, Random Access heavy"
    }
]

def run_command(cmd):
    try:
        start = time.time()
        # Using a large buffer size for pipe to avoid bottlenecks
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, bufsize=1024*1024)
        out, err = proc.communicate()
        end = time.time()
        return out, proc.returncode, (end - start)
    except FileNotFoundError:
        return None, -1, 0

def parse_val(output, key):
    if not output: return 0.0
    for line in output.splitlines():
        if key in line:
            return float(line.split(":")[1].strip())
    return 0.0

def generate_data(filename, scenario):
    print(f"  Generating {filename} ({scenario['num_seqs']} seqs)...")
    # Pre-generate a 1MB buffer to slice from (faster than randomizing every char)
    buffer_len = 1_000_000
    dna_buffer = "".join(random.choices("ACGT", k=buffer_len))
    
    with open(filename, "w") as f:
        for i in range(scenario['num_seqs']):
            # Jitter length by +/- 10%
            length = int(scenario['avg_len'] * random.uniform(0.9, 1.1))
            start = random.randint(0, buffer_len - length if length < buffer_len else 0)
            
            if length > buffer_len:
                 # Construct massive seqs on the fly
                 seq = "".join(random.choices("ACGT", k=length))
            else:
                 seq = dna_buffer[start : start + length]

            if scenario['type'] == 'fasta':
                f.write(f">{scenario['prefix']}{i} synthetic_data\n{seq}\n")
            else:
                qual = "I" * len(seq) # Phred 40
                f.write(f"@{scenario['prefix']}{i} synthetic_data\n{seq}\n+\n{qual}\n")

def calculate_adaptive_iterations(scenario):
    """
    Adjusts lookup count to process ~500MB of data volume.
    Ensures benchmarks don't freeze on Reference Genomes but run long enough for WES/WGS.
    """
    avg_size_bytes = scenario['avg_len']
    target_bytes = TARGET_DATA_VOLUME_MB * 1024 * 1024
    
    count = int(target_bytes / avg_size_bytes)
    
    # Safety Clamps
    if count < 50: count = 50 
    if count > MAX_LOOKUPS_CAP: count = MAX_LOOKUPS_CAP
    
    return count

# === BENCHMARK CORE ===

def bench_biopython_parse(filename, fmt):
    start = time.time()
    # Simply iterate to force parsing
    for _ in SeqIO.parse(filename, fmt): pass
    return time.time() - start

def bench_biopython_lookup(filename, fmt, num_lookups, prefix, max_idx):
    # 1. Load (Dict) - The expensive part
    try:
        # Use simple dict for WGS/WES, might OOM on Genome if not careful
        record_dict = SeqIO.to_dict(SeqIO.parse(filename, fmt))
    except MemoryError:
        return -1 # Fail Gracefully
        
    # 2. Lookup - The fast part
    t0 = time.time()
    random.seed(42)
    for _ in range(num_lookups):
        k = f"{prefix}{random.randint(0, max_idx-1)}"
        _ = record_dict.get(k)
    return time.time() - t0

def run_scenario(s):
    base_name = f"bench_{s['name'].split()[0]}"
    fname = f"{base_name}.{s['type']}"
    cname = f"{base_name}.bin"
    
    iter_count = calculate_adaptive_iterations(s)
    
    print(f"\n[SCENARIO: {s['name']}]")
    print(f"  Configuration: {s['num_seqs']:,} seqs, avg len {s['avg_len']}")
    print(f"  Adaptive Lookups: {iter_count:,} iterations (~{TARGET_DATA_VOLUME_MB}MB volume)")

    if not os.path.exists(fname): generate_data(fname, s)
    
    res = {}
    
    # --- 1. PARSING & IO ---
    print("  > Benchmarking Parsing Speed...")
    
    # SeqTk (C Stream) - The Baseline
    out, rc, duration = run_command([SEQTK_BINARY, "fqchk", fname])
    res['SeqTk'] = duration if rc == 0 else None
    
    # BioPython (Python Stream)
    res['BioPy_Parse'] = bench_biopython_parse(fname, s['type'])
    
    # TracEon (C++ Load)
    out, rc, _ = run_command([TRACEON_BINARY, "load", fname])
    res['TracEon_Load'] = parse_val(out, "Load_Time_s")
    
    # --- 2. LIFECYCLE (Save/Restore) ---
    run_command([TRACEON_BINARY, "save", fname, cname])
    
    out, rc, _ = run_command([TRACEON_BINARY, "restore", cname])
    res['TracEon_Restore'] = parse_val(out, "Restore_Time_s")
    
    # --- 3. RANDOM ACCESS LOOKUPS ---
    print(f"  > Benchmarking Random Access (Alignment I/O Sim)...")
    
    # BioPython Dict
    t_bio = bench_biopython_lookup(fname, s['type'], iter_count, s['prefix'], s['num_seqs'])
    res['BioPy_Lookup_Rate'] = (iter_count / t_bio) if t_bio > 0 else 0

    # TracEon Lookup
    out, rc, _ = run_command([TRACEON_BINARY, "lookup", cname, str(iter_count), s['prefix'], str(s['num_seqs'])])
    res['TracEon_Lookup_Rate'] = parse_val(out, "Throughput")

    # Cleanup
    if os.path.exists(fname): os.remove(fname)
    if os.path.exists(cname): os.remove(cname)
    
    return res

def print_summary(all_res):
    print("\n" + "="*125)
    print(f"{'SCENARIO':<35} | {'METRIC':<18} | {'SeqTk':<10} | {'BioPython':<12} | {'TracEon':<12} | {'IMPROVEMENT'}")
    print("-" * 125)
    
    for name, r in all_res.items():
        # 1. Parsing Row
        stk = f"{r['SeqTk']:.4f}s" if r['SeqTk'] else "N/A"
        bio = f"{r['BioPy_Parse']:.4f}s"
        tra = f"{r['TracEon_Load']:.4f}s"
        
        # Improvement: TracEon Load vs BioPy Parse
        imp = f"{r['BioPy_Parse']/r['TracEon_Load']:.1f}x (vs Py)" if r['TracEon_Load'] > 0 else "-"
        print(f"{name:<35} | {'Parse Time':<18} | {stk:<10} | {bio:<12} | {tra:<12} | {imp}")
        
        # 2. Restore Row
        tra_rest = f"{r['TracEon_Restore']:.4f}s"
        print(f"{'':<35} | {'Restore Time':<18} | {'-':<10} | {'-':<12} | {tra_rest:<12} | -")
        
        # 3. Lookup Row
        bp_ops = r['BioPy_Lookup_Rate']
        tr_ops = r['TracEon_Lookup_Rate']
        if bp_ops > 0:
            speedup = f"{tr_ops/bp_ops:.1f}x"
        else:
            speedup = "Inf"
            
        print(f"{'':<35} | {'Lookups/sec':<18} | {'-':<10} | {int(bp_ops):<12,} | {int(tr_ops):<12,} | {speedup}")
        print("-" * 125)

def main():
    if not os.path.exists(TRACEON_BINARY):
        print(f"Error: {TRACEON_BINARY} not found. Please compile traceon_driver first.")
        return
    
    results = {}
    for s in SCENARIOS:
        results[s['name']] = run_scenario(s)
        
    print_summary(results)

if __name__ == "__main__":
    main()