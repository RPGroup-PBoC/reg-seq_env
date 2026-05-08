"""
Scan the E. coli genome with the de novo PWM to find putative
additional targets of the unknown salt-response TF.

Run after run_analysis() — needs denovo_pwm and genome in namespace.

Usage in notebook:
    from scan_genome_motif import scan_genome, plot_hits
    hits = scan_genome(
        denovo_pwm, genome_path='../data/mg1655_genome.fasta',
        gff_path=None,  # optional: GFF for gene annotations
        threshold_pct=99.9,
    )
"""

import numpy as np
import pandas as pd

BASE_TO_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def load_genome(fasta_path):
    parts = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                parts.append(line.strip().upper())
    return ''.join(parts)


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def score_sequence(seq, pwm):
    """Score a single sequence against a PWM."""
    s = 0.0
    for i, base in enumerate(seq):
        if base in BASE_TO_IDX:
            s += pwm[i, BASE_TO_IDX[base]]
    return s


def scan_genome(denovo_pwm, genome_path, threshold_pct=99.9,
                window_around_gene=400):
    """
    Slide the de novo PWM across the entire genome (both strands)
    and return positions with scores above the threshold.

    Parameters
    ----------
    denovo_pwm : ndarray (L x 4)
    genome_path : str
    threshold_pct : float
        Percentile of all scores to use as threshold.
        99.9 means top 0.1% of all genomic positions.
    window_around_gene : int
        For annotation, report the nearest gene within this
        many bp upstream.

    Returns
    -------
    hits_df : DataFrame with columns:
        position, strand, score, sequence
    threshold : float
    all_scores : ndarray (for diagnostics)
    """
    genome = load_genome(genome_path)
    L = denovo_pwm.shape[0]
    genome_len = len(genome)

    print(f"Scanning genome ({genome_len:,} bp) with {L}-bp motif...")

    # Score every position on both strands
    n_positions = genome_len - L + 1
    fwd_scores = np.zeros(n_positions)
    rev_scores = np.zeros(n_positions)

    for i in range(n_positions):
        subseq = genome[i:i + L]
        if 'N' in subseq:
            fwd_scores[i] = -np.inf
            rev_scores[i] = -np.inf
            continue
        fwd_scores[i] = score_sequence(subseq, denovo_pwm)
        rev_scores[i] = score_sequence(reverse_complement(subseq), denovo_pwm)

    # Combine: for each position, take best strand
    all_scores = np.maximum(fwd_scores, rev_scores)

    # Threshold
    valid = all_scores[all_scores > -np.inf]
    threshold = np.percentile(valid, threshold_pct)
    print(f"  Threshold ({threshold_pct}th percentile): {threshold:.2f}")
    print(f"  Score range: [{valid.min():.2f}, {valid.max():.2f}]")
    print(f"  Median score: {np.median(valid):.2f}")

    # Collect hits
    hits = []
    for i in range(n_positions):
        best_score = max(fwd_scores[i], rev_scores[i])
        if best_score >= threshold:
            if fwd_scores[i] >= rev_scores[i]:
                strand = '+'
                seq = genome[i:i + L]
            else:
                strand = '-'
                seq = reverse_complement(genome[i:i + L])
            hits.append({
                'position': i + 1,  # 1-based
                'strand': strand,
                'score': best_score,
                'sequence': seq,
            })

    hits_df = pd.DataFrame(hits).sort_values('score', ascending=False)
    print(f"  Found {len(hits_df)} positions above threshold")

    return hits_df, threshold, all_scores


def annotate_hits_with_genes(hits_df, genome_path, gene_table_path=None):
    """
    Annotate hits with nearby genes using a simple gene table.

    If no gene table is provided, attempts to build one from
    common E. coli annotation sources.

    gene_table should have columns: gene, start, end, strand
    """
    if gene_table_path is not None:
        genes = pd.read_csv(gene_table_path)
    else:
        print("  No gene table provided. Returning hits without gene annotation.")
        print("  To annotate, provide a CSV with columns: gene, start, end, strand")
        return hits_df

    annotated = []
    for _, hit in hits_df.iterrows():
        pos = hit['position']
        # Find genes within 500 bp downstream (hit is in promoter region)
        nearby = []
        for _, g in genes.iterrows():
            if g['strand'] == '+' and 0 < g['start'] - pos < 500:
                nearby.append((g['gene'], g['start'] - pos, '+'))
            elif g['strand'] == '-' and 0 < pos - g['end'] < 500:
                nearby.append((g['gene'], pos - g['end'], '-'))

        if nearby:
            # Pick closest
            nearby.sort(key=lambda x: x[1])
            hit_dict = hit.to_dict()
            hit_dict['nearest_gene'] = nearby[0][0]
            hit_dict['distance'] = nearby[0][1]
            hit_dict['gene_strand'] = nearby[0][2]
        else:
            hit_dict = hit.to_dict()
            hit_dict['nearest_gene'] = ''
            hit_dict['distance'] = np.nan
            hit_dict['gene_strand'] = ''
        annotated.append(hit_dict)

    return pd.DataFrame(annotated)


def merge_nearby_hits(hits_df, merge_distance=30):
    """
    Merge hits within merge_distance bp, keeping the highest score.
    This collapses the multiple overlapping positions from a single
    binding site into one representative hit.
    """
    if len(hits_df) == 0:
        return hits_df

    hits = hits_df.sort_values('position').reset_index(drop=True)
    merged = []
    current = hits.iloc[0].to_dict()

    for i in range(1, len(hits)):
        row = hits.iloc[i]
        if row['position'] - current['position'] < merge_distance:
            # Same site — keep higher score
            if row['score'] > current['score']:
                current = row.to_dict()
        else:
            merged.append(current)
            current = row.to_dict()
    merged.append(current)

    result = pd.DataFrame(merged).sort_values('score', ascending=False)
    print(f"  Merged {len(hits_df)} positions into {len(result)} distinct sites")
    return result


def filter_known_sites(hits_df, known_positions, exclude_distance=50):
    """
    Remove hits near the three known binding sites (yjbJ, ybaY, yqjE).

    known_positions : list of int
        Genomic positions of the known binding sites.
    exclude_distance : int
        Remove hits within this many bp of a known site.
    """
    mask = pd.Series(True, index=hits_df.index)
    for kp in known_positions:
        mask &= (hits_df['position'] - kp).abs() > exclude_distance
    filtered = hits_df[mask]
    n_removed = len(hits_df) - len(filtered)
    if n_removed > 0:
        print(f"  Removed {n_removed} hits near known binding sites")
    return filtered