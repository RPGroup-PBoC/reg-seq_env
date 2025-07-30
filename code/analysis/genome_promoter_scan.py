from Bio import SeqIO
import regseq2
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# ——— 1. parameters ———
FASTA   = "../../data/metadata/mg1655_genome.fasta"
MAX_BP  = 10_000          # adjust as needed
WINDOW  = 78
UP      = WINDOW - 1     # 77 bp upstream
DOWN    = WINDOW - 1     # 77 bp downstream
BASES   = ['A','C','G','T']
OUTCSV  = "promoter_scan_results.csv"

# ——— 2. load genome ———
record = next(SeqIO.parse(FASTA, "fasta"))
genome = str(record.seq[:MAX_BP])

# ——— 3. worker function ———
def scan_mutation(args):
    pos, alt = args
    orig = genome[pos]

    # build mutated chunk around the position
    mut = genome[:pos] + alt + genome[pos+1:]
    start = max(0, pos - UP)
    end   = min(len(mut), pos + DOWN + 1)
    chunk = mut[start:end]

    # fresh model instance per process
    model = regseq2.promoter_calculator.Promoter_Calculator()
    model.run(chunk)
    out = model.output()

    rows = []
    for direction in ("Forward", "Reverse"):
        key = f"{direction}_Predictions_per_TSS"
        for tss, pred in out.get(key, {}).items():
            # compute genomic window start
            ws = start + (tss - WINDOW)
            # only keep windows covering the mutated base
            if ws <= pos < ws + WINDOW:
                # base row info
                row = {
                    "mut_pos":      pos,
                    "orig_base":    orig,
                    "alt_base":     alt,
                    "direction":    direction,
                    "window_start": ws,
                    "model_TSS":    int(pred['TSS'])
                }
                # include only numeric metrics (drop sequences)
                for k, v in pred.items():
                    if k == 'TSS':
                        continue
                    if hasattr(v, "dtype") and np.isscalar(v):
                        row[k] = float(v)
                    elif isinstance(v, (int, float)):
                        row[k] = v
                rows.append(row)
    return rows


def main():
    # ——— 4. build task list ———
    tasks = [(pos, alt)
             for pos in range(MAX_BP)
             for alt in BASES if alt != genome[pos]]

    # ——— 5. parallel execution & collect ———
    all_rows = []
    with Pool(processes=cpu_count()) as pool:
        for result in tqdm(pool.imap_unordered(scan_mutation, tasks, chunksize=100),
                            total=len(tasks), desc="Mutational scan (parallel)"):
            all_rows.extend(result)

    # ——— 6. save to CSV ———
    df = pd.DataFrame(all_rows)
    df.to_csv(OUTCSV, index=False)
    print(f"Wrote {len(df)} rows to {OUTCSV}")


if __name__ == '__main__':
    main()

