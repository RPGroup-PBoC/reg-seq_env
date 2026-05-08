"""
Fixed version of score_denovo_motif.py

Key fixes:
1. compare_motifs: when de novo PWM is wider than known site sequences,
   extract extended genomic context around each known binding site before scoring.
2. compute_null_distribution: sample random sequences long enough for the de novo PWM.
3. Handles -inf gracefully in normalization and plotting.
4. compare_motifs now returns self_score and norm_score columns.

All other logic unchanged from the original.
"""
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'

BASE_MAP = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}
BASE_TO_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def load_genome(fasta_path):
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                seq_parts.append(line.strip().upper())
    return ''.join(seq_parts)


def reverse_complement(seq):
    comp = str.maketrans('ACGT', 'TGCA')
    return seq.translate(comp)[::-1]


def extract_site(genome, center, width, strand):
    half = width // 2
    start = int(center) - half
    end = start + width
    if start < 0 or end > len(genome):
        return None
    seq = genome[start:end]
    if strand == 'reverse' or strand == '-':
        seq = reverse_complement(seq)
    return seq


def reshape_shifts(shift_df, condition, promoters):
    cond_str = str(condition)
    mask = shift_df['ind'].astype(str).str.startswith(cond_str + '-')
    df = shift_df[mask].copy()
    if len(df) == 0:
        raise ValueError(f"No data for condition {condition}")
    print(f"  Condition {condition}: {len(df)} rows, reps: {sorted(df['ind'].unique())}")
    df['base_letter'] = df['base'].map(BASE_MAP)
    results = {}
    for prom in promoters:
        prom_df = df[df['promoter'] == prom]
        if len(prom_df) == 0:
            print(f"  WARNING: no data for '{prom}'")
            print(f"    Available: {sorted(df['promoter'].unique())[:10]}...")
            continue
        avg = prom_df.groupby(['pos', 'base_letter', 'wt_base'])['expression_shift'].mean().reset_index()
        pivot = avg.pivot_table(index=['pos', 'wt_base'], columns='base_letter', values='expression_shift').reset_index()
        pivot.columns.name = None
        pivot = pivot.sort_values('pos').reset_index(drop=True)
        for b in ['A', 'C', 'G', 'T']:
            if b not in pivot.columns:
                pivot[b] = 0.0
        results[prom] = pivot
        print(f"  {prom}: {len(pivot)} positions")
    return results


def get_wt_sequence(promoter_name, wt_sequences_csv):
    df = pd.read_csv(wt_sequences_csv)
    row = df[df['name'].str.contains(promoter_name, na=False)]
    if len(row) == 0:
        raise ValueError(f"Promoter {promoter_name} not found")
    return row.iloc[0]['promoter_seq']


def build_denovo_pwm(shift_matrices, binding_sites, wt_sequences_csv,
                     tss_in_seq=120, beta=1.0):
    first_key = list(binding_sites.keys())[0]
    site_width = binding_sites[first_key][1] - binding_sites[first_key][0]
    energy_matrices = []
    wt_site_seqs = {}
    site_details = {}

    for prom, (bs_start, bs_end) in binding_sites.items():
        if prom not in shift_matrices:
            continue
        pivot = shift_matrices[prom]
        wt_seq = get_wt_sequence(prom, wt_sequences_csv)
        wt_site = ''
        for pos in range(bs_start, bs_end):
            seq_idx = tss_in_seq + pos
            wt_site += wt_seq[seq_idx] if 0 <= seq_idx < len(wt_seq) else 'N'
        wt_site_seqs[prom] = wt_site

        E = np.zeros((site_width, 4))
        n_found = 0
        for i, pos in enumerate(range(bs_start, bs_end)):
            row = pivot[pivot['pos'] == pos]
            if len(row) == 0:
                continue
            n_found += 1
            row = row.iloc[0]
            for base, idx in BASE_TO_IDX.items():
                val = row.get(base, 0)
                E[i, idx] = 0.0 if pd.isna(val) else val
            seq_idx = tss_in_seq + pos
            if 0 <= seq_idx < len(wt_seq):
                wt_base = wt_seq[seq_idx]
                if wt_base in BASE_TO_IDX:
                    E[i, BASE_TO_IDX[wt_base]] = 0.0
        energy_matrices.append(E)
        site_details[prom] = {'wt_sequence': wt_site, 'energy_matrix': E, 'n_found': n_found}
        print(f"  {prom}: WT={wt_site}, {n_found}/{site_width} pos with data")

    if not energy_matrices:
        raise ValueError("No valid energy matrices")

    weights = np.array([np.sum(np.abs(E)) for E in energy_matrices])
    weights = weights / weights.sum()
    combined_E = sum(w * E for w, E in zip(weights, energy_matrices))

    preferences = combined_E
    freqs = np.zeros_like(preferences)
    for i in range(site_width):
        exp_vals = np.exp(beta * preferences[i, :])
        freqs[i, :] = exp_vals / exp_vals.sum()

    freqs = np.clip(freqs, 1e-4, 1.0)
    freqs = freqs / freqs.sum(axis=1, keepdims=True)
    pwm = np.log2(freqs / 0.25)

    bases = ['A', 'C', 'G', 'T']
    consensus = ''.join(bases[j] for j in np.argmax(freqs, axis=1))
    print(f"\n  Consensus: {consensus}")
    for prom, seq in wt_site_seqs.items():
        matches = sum(1 for a, b in zip(consensus, seq) if a == b)
        print(f"    {prom}: {seq}  ({matches}/{site_width} match)")

    return pwm, freqs, consensus, site_details


def build_known_pwms(binding_sites_csv, genome, site_width=20, min_sites=8,
                     filter_evidence=None):
    df = pd.read_csv(binding_sites_csv)
    if filter_evidence:
        if type(filter_evidence) == str:
            df = df[df['confidence_level'] == filter_evidence]
        elif type(filter_evidence) == list:
            df = df[[x in filter_evidence for x in df['confidence_level'].values]]

    tf_groups = defaultdict(list)
    for _, row in df.iterrows():
        try:
            center = float(row['bs_position'])
        except (ValueError, TypeError):
            continue
        tf_groups[row['tf']].append((center, str(row['strand'])))

    pwms, site_seqs = {}, {}
    # Also store genomic positions for later context extraction
    site_positions = {}
    for tf, positions in tf_groups.items():
        seqs = []
        pos_list = []
        for center, strand in positions:
            seq = extract_site(genome, center, site_width, strand)
            if seq and len(seq) == site_width and 'N' not in seq:
                seqs.append(seq)
                pos_list.append((center, strand))
        if len(seqs) < min_sites:
            continue
        L = len(seqs[0])
        counts = np.zeros((L, 4))
        for seq in seqs:
            for i, base in enumerate(seq):
                if base in BASE_TO_IDX:
                    counts[i, BASE_TO_IDX[base]] += 1
        ps = 0.5
        f = (counts + ps) / (counts.sum(axis=1, keepdims=True) + 4 * ps)
        pwms[tf] = np.log2(f / 0.25)
        site_seqs[tf] = seqs
        site_positions[tf] = pos_list

    return pwms, site_seqs, site_positions


def _score_seq_with_pwm(seq, pwm):
    """Score a sequence with a PWM, scanning all offsets and both strands.
    Returns best score. seq must be >= pwm length."""
    L_pwm = pwm.shape[0]
    if len(seq) < L_pwm:
        return -np.inf
    best = -np.inf
    for s_seq in [seq, reverse_complement(seq)]:
        for start in range(len(s_seq) - L_pwm + 1):
            sub = s_seq[start:start + L_pwm]
            sc = sum(pwm[i, BASE_TO_IDX.get(sub[i], 0)]
                     for i in range(L_pwm))
            best = max(best, sc)
    return best


def compute_null_distribution(genome, known_pwms, known_site_seqs,
                              denovo_pwm, n_random=1000, seed=42):
    """
    Empirical null: score random genomic sequences with the de novo
    PWM, normalize by each TF's self-score.

    FIX: sample sequences long enough for the de novo PWM.
    """
    rng = np.random.RandomState(seed)
    genome_len = len(genome)
    L_d = denovo_pwm.shape[0]

    # Sample random sequences at LEAST as long as the de novo PWM
    # Use the same length as the extended site sequences we'll use for scoring
    sample_len = max(L_d + 10, 30)  # some padding for scanning
    margin = 200
    random_positions = rng.randint(margin, genome_len - margin - sample_len,
                                   size=n_random)

    random_seqs = []
    for pos in random_positions:
        seq = genome[pos:pos + sample_len]
        if 'N' not in seq and len(seq) == sample_len:
            random_seqs.append(seq)
    print(f"  Sampled {len(random_seqs)} random {sample_len}-mers "
          f"(de novo PWM width: {L_d})")

    all_null_norms = []

    for tf, known_pwm in known_pwms.items():
        L_k = known_pwm.shape[0]

        # Self-score: true sites scored by own PWM
        true_scores = []
        for seq in known_site_seqs[tf]:
            s = sum(known_pwm[i, BASE_TO_IDX.get(seq[i], 0)]
                    for i in range(L_k))
            true_scores.append(s)
        self_score = np.mean(true_scores)
        if self_score <= 0:
            continue

        # Score random seqs with de novo PWM
        for seq in random_seqs:
            best_s = _score_seq_with_pwm(seq, denovo_pwm)
            if np.isfinite(best_s):
                all_null_norms.append(best_s / self_score)

    all_null_norms = np.array(all_null_norms)
    if len(all_null_norms) == 0 or not np.any(np.isfinite(all_null_norms)):
        print("  WARNING: null distribution empty or all non-finite")
        return all_null_norms, 0.0, 0.0

    threshold_95 = np.percentile(all_null_norms, 95)
    threshold_99 = np.percentile(all_null_norms, 99)

    print(f"  Null norm_score: median={np.median(all_null_norms):.3f}, "
          f"95th={threshold_95:.3f}, 99th={threshold_99:.3f}")

    return all_null_norms, threshold_95, threshold_99


def compare_motifs(denovo_pwm, known_pwms, known_site_seqs,
                   genome=None, known_site_positions=None):
    """
    Compare de novo PWM against all known TF motifs.

    FIX: when the de novo PWM is wider than known site sequences,
    extract extended genomic context around each known binding site
    so the de novo PWM can be scanned properly.

    Returns DataFrame with columns:
      tf, pwm_correlation, mean_site_score, self_score, norm_score,
      n_sites, corr_rank, score_rank, combined_rank
    """
    L_d = denovo_pwm.shape[0]
    rows = []

    for tf, known_pwm in known_pwms.items():
        # --- PWM correlation (handles size mismatch already) ---
        denovo_rc = denovo_pwm[::-1, ::-1]
        best_corr = -1
        for qpwm in [denovo_pwm, denovo_rc]:
            if qpwm.shape[0] <= known_pwm.shape[0]:
                small, big = qpwm, known_pwm
            else:
                small, big = known_pwm, qpwm
            for off in range(big.shape[0] - small.shape[0] + 1):
                a = small.flatten()
                b = big[off:off + small.shape[0], :].flatten()
                if np.std(a) < 1e-10 or np.std(b) < 1e-10:
                    continue
                best_corr = max(best_corr, np.corrcoef(a, b)[0, 1])

        # --- Site scoring with de novo PWM ---
        site_seqs_to_score = known_site_seqs.get(tf, [])
        L_site = len(site_seqs_to_score[0]) if site_seqs_to_score else 0

        # If de novo PWM is wider than site sequences, try to extend
        # using genomic context
        if L_d > L_site and genome is not None and known_site_positions is not None:
            extended_seqs = []
            extend_by = L_d - L_site + 10  # extra padding
            positions = known_site_positions.get(tf, [])
            for center, strand in positions:
                ext_width = L_site + extend_by
                seq = extract_site(genome, center, ext_width, strand)
                if seq and len(seq) == ext_width and 'N' not in seq:
                    extended_seqs.append(seq)
            if extended_seqs:
                site_seqs_to_score = extended_seqs

        scores = []
        for seq in site_seqs_to_score:
            best_s = _score_seq_with_pwm(seq, denovo_pwm)
            scores.append(best_s)

        mean_score = np.mean(scores) if scores else np.nan
        if not np.isfinite(mean_score):
            mean_score = np.nan

        # --- Self-score: known sites scored by their own PWM ---
        L_k = known_pwm.shape[0]
        self_scores = []
        for seq in known_site_seqs.get(tf, []):
            if len(seq) >= L_k:
                s = sum(known_pwm[i, BASE_TO_IDX.get(seq[i], 0)]
                        for i in range(L_k))
                self_scores.append(s)
        self_score = np.mean(self_scores) if self_scores else np.nan

        # Normalized score: de novo site score / self score
        norm_score = np.nan
        if (np.isfinite(mean_score) and np.isfinite(self_score)
                and self_score > 0):
            norm_score = mean_score / self_score

        rows.append({
            'tf': tf,
            'pwm_correlation': best_corr,
            'mean_site_score': mean_score,
            'self_score': self_score,
            'norm_score': norm_score,
            'n_sites': len(scores)
        })

    results = pd.DataFrame(rows)
    results['corr_rank'] = results['pwm_correlation'].rank(ascending=False)
    results['score_rank'] = results['norm_score'].rank(ascending=False,
                                                       na_option='bottom')
    results['combined_rank'] = (results['corr_rank'] + results['score_rank']) / 2
    return results.sort_values('combined_rank')


def run_analysis(shift_df, binding_sites, condition, highlight_tfs,
                 genome_path, binding_sites_csv, wt_sequences_csv,
                 tss_in_seq=120, beta=1.0, site_width_known=20,
                 min_sites=8, n_random=1000, output_dir='.'):
    import os

    print("=" * 60)
    print("STEP 1: Reshape expression shifts")
    print("=" * 60)
    shift_matrices = reshape_shifts(shift_df, condition, list(binding_sites.keys()))

    print("\n" + "=" * 60)
    print("STEP 2: Build de novo PWM")
    print("=" * 60)
    denovo_pwm, freqs, consensus, details = build_denovo_pwm(
        shift_matrices, binding_sites, wt_sequences_csv,
        tss_in_seq=tss_in_seq, beta=beta)
    wt_seqs = {p: d['wt_sequence'] for p, d in details.items()}

    print("\n" + "=" * 60)
    print("STEP 3: Build known TF PWMs")
    print("=" * 60)
    genome = load_genome(genome_path)
    print(f"  Genome: {len(genome):,} bp")
    known_pwms, known_seqs, known_positions = build_known_pwms(
        binding_sites_csv, genome,
        site_width=site_width_known, min_sites=min_sites)
    print(f"  PWMs for {len(known_pwms)} TFs")

    print("\n" + "=" * 60)
    print("STEP 4: Empirical null distribution")
    print("=" * 60)
    null_norms, threshold_95, threshold_99 = compute_null_distribution(
        genome, known_pwms, known_seqs, denovo_pwm, n_random=n_random)

    print("\n" + "=" * 60)
    print("STEP 5: Compare motifs")
    print("=" * 60)
    results = compare_motifs(denovo_pwm, known_pwms, known_seqs,
                             genome=genome,
                             known_site_positions=known_positions)

    # De novo self-score: score the WT site sequences against the de novo PWM
    denovo_self_scores = []
    for prom, d in details.items():
        wt_seq = d['wt_sequence']
        if len(wt_seq) >= denovo_pwm.shape[0] and 'N' not in wt_seq:
            s = _score_seq_with_pwm(wt_seq, denovo_pwm)
            if np.isfinite(s):
                denovo_self_scores.append(s)
    denovo_self_score = np.mean(denovo_self_scores) if denovo_self_scores else np.nan
    results['denovo_self_score'] = denovo_self_score
    print(f"\n  De novo self-score: {denovo_self_score:.3f} "
          f"(from {len(denovo_self_scores)} sites)")

    print(f"\n  {'TF':<15} {'Corr':>6} {'SiteScr':>8} {'SelfScr':>8} {'Norm':>6} {'Rank':>6} {'n':>5}")
    print(f"  {'-'*15} {'-'*6} {'-'*8} {'-'*8} {'-'*6} {'-'*6} {'-'*5}")
    shown = set()
    for _, row in results.head(15).iterrows():
        flag = ' *' if row['tf'] in highlight_tfs else ''
        site_str = (f"{row['mean_site_score']:>8.2f}"
                    if pd.notna(row['mean_site_score']) else '     NaN')
        self_str = (f"{row['self_score']:>8.2f}"
                    if pd.notna(row['self_score']) else '     NaN')
        norm_str = (f"{row['norm_score']:>6.3f}"
                    if pd.notna(row['norm_score']) else '   NaN')
        print(f"  {row['tf']:<15} {row['pwm_correlation']:>6.3f} "
              f"{site_str} {self_str} {norm_str} "
              f"{row['combined_rank']:>6.1f} {int(row['n_sites']):>5}{flag}")
        shown.add(row['tf'])
    for _, row in results.iterrows():
        if row['tf'] in highlight_tfs and row['tf'] not in shown:
            site_str = (f"{row['mean_site_score']:>8.2f}"
                        if pd.notna(row['mean_site_score']) else '     NaN')
            self_str = (f"{row['self_score']:>8.2f}"
                        if pd.notna(row['self_score']) else '     NaN')
            norm_str = (f"{row['norm_score']:>6.3f}"
                        if pd.notna(row['norm_score']) else '   NaN')
            print(f"  {row['tf']:<15} {row['pwm_correlation']:>6.3f} "
                  f"{site_str} {self_str} {norm_str} "
                  f"{row['combined_rank']:>6.1f} {int(row['n_sites']):>5} *")

    results_path = os.path.join(output_dir, 'denovo_motif_comparison.csv')
    results.to_csv(results_path, index=False)
    print(f"\n  Saved: {results_path}")

    return (denovo_pwm, known_pwms, freqs, consensus, results, details,
            known_seqs, null_norms, threshold_95, threshold_99,
            denovo_self_scores)