"""Information-footprint analysis for Reg-Seq barcode count data.

Loads per-sample DNA/RNA barcode counts and computes, for each position in a
promoter, the mutual information (in bits) between "is this base mutated away
from WT" and "was this read observed in the DNA or the RNA library". This is
the canonical implementation of `get_counts`/`compute_MI_from_counts`/
`compute_MI` — notebooks should import these three from here (`from regseq2
import compute_MI, ...`) rather than redefining them, so there is only ever
one version of the core MI math.

Scrambling / null-distribution functions (`scramble_counts`,
`get_counts_scrambled`, `compute_MI_scrambled`, `compute_MI_null_percentile`)
are still under active development and intentionally NOT here — they're
defined directly in the notebooks that use them (they reuse
`_prepare_promoter_arrays`/`_mi_footprint` from this module so their MI math
stays identical to `compute_MI_from_counts`, without duplicating it).

`compute_footprints.ipynb` is the one place these functions are additionally
shown inline (for explanation/derivation); it should be kept in sync with this
module by hand.
"""
import os
from functools import lru_cache
from glob import glob

import numpy as np
import pandas as pd

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_PACKAGE_DIR))
DATA_DIR = os.path.join(_REPO_ROOT, "data")

_EXCLUDED_PROMOTERS = ("galEp", "ybeDp2")


@lru_cache(maxsize=1)
def _load_df_map():
    """Barcode → promoter mapping table, read once per process."""
    return pd.read_csv(os.path.join(DATA_DIR, "mapping", "mapped_barcodes_filtered.csv"))


@lru_cache(maxsize=1)
def _load_wt_dict():
    """promoter name -> WT (unmutated) 160 nt sequence, read once per process."""
    df_wt = pd.read_csv(os.path.join(DATA_DIR, "metadata", "wt_sequences.csv"))[["promoter_seq", "promoter"]]
    return dict(zip(df_wt["promoter"].values, df_wt["promoter_seq"].values))


def get_counts(ind, promoter=None):
    """Load DNA + RNA barcode counts for sample `ind` and join them with the barcode → promoter map.

    Returns one row per barcode with columns: ct_0 (DNA reads), ct_1 (RNA reads),
    barcode, name (promoter id), promoter (160 nt mutated sequence).

    `promoter` can be None (all promoters), a str, or a list of strs.
    """
    barcode_dir = os.path.join(DATA_DIR, "extracted_barcodes")

    # DNA file: ct_0 is the input library count
    file_DNA = glob(f"{barcode_dir}/{ind}*DNA*.txt")[0]
    df_DNA = pd.read_csv(file_DNA, names=["ct_0", "barcode"], sep="\\s+")
    # RNA file: ct_1 is the readout count
    file_RNA = glob(f"{barcode_dir}/{ind}*RNA*.txt")[0]
    df_RNA = pd.read_csv(file_RNA, names=["ct_1", "barcode"], sep="\\s+")

    # outer-join keeps barcodes seen in either DNA or RNA; inner-join with df_map drops unmapped barcodes
    df_counts = df_DNA.merge(df_RNA, on="barcode", how="outer").fillna(0)
    df_counts = df_counts.merge(_load_df_map(), on="barcode", how="inner")

    # filter to requested promoter(s)
    if promoter is None:
        return df_counts
    elif type(promoter) == str:
        return df_counts[df_counts["name"] == promoter]
    elif type(promoter) == list:
        return df_counts[[x in promoter for x in df_counts["name"].values]]


def _transform_counts(c, transform="raw"):
    """Reshape count weights to dampen highly-expressed barcodes.

    "sqrt" / "log1p" / "pow:a" compress dynamic range; "clip:k" caps extreme values.
    """
    c = c.astype(np.float64)
    if transform == "raw":
        return c
    if transform == "sqrt":
        return np.sqrt(c)
    if transform == "log1p":
        return np.log1p(c)
    if transform.startswith("pow:"):
        a = float(transform.split(":")[1])
        return np.power(c, a)
    if transform.startswith("clip:"):
        k = float(transform.split(":")[1])
        return np.minimum(c, k)
    raise ValueError(f"Unknown transform: {transform}")


def _mi_2x2_vectorized(joint):
    """MI in bits for a stack of 2×2 joint distributions, shape (L, 2, 2), each summing to 1.

    Vectorized replacement for calling a per-position MI routine in a Python loop:
    all L positions are computed together with array ops instead of L separate calls.
    """
    px = joint.sum(axis=2)  # (L, 2): P(x) marginal per position
    pz = joint.sum(axis=1)  # (L, 2): P(z) marginal per position
    outer = px[:, :, None] * pz[:, None, :]  # (L, 2, 2): P(x)P(z) under independence
    ratio = np.divide(joint, outer, out=np.ones_like(joint), where=outer > 0)
    terms = np.where((joint > 0) & (outer > 0), joint * np.log2(ratio), 0.0)
    return terms.sum(axis=(1, 2))


def _h2_vectorized(p):
    """Bernoulli entropy in bits, elementwise; `p` may be scalar or array-valued."""
    p = np.asarray(p, dtype=np.float64)
    h = np.zeros_like(p)
    valid = (p > 0) & (p < 1)
    h[valid] = -p[valid] * np.log2(p[valid]) - (1.0 - p[valid]) * np.log2(1.0 - p[valid])
    return h


def _prepare_promoter_arrays(gdf, promoter, *, transform="sqrt", barcode_equal=False, require_both_counts=False):
    """Precompute the per-barcode data for one promoter that is independent of how
    ct_0/ct_1 pairs get assigned to barcodes: the WT-match matrix and the
    (transformed, optionally barcode-equalized) DNA/RNA weight vectors.

    Scrambling only changes which weight-vector entry lines up with which
    `is_match` row (see `_mi_footprint`), so this expensive part — parsing
    sequences, comparing against WT, transforming counts — only needs to run once
    per promoter no matter how many scrambles are drawn from it.

    Parameters
    ----------
    gdf : pandas.DataFrame
        Counts already filtered to a single `promoter` (e.g. one group from
        `df_counts.groupby("name")`, or `df_counts[df_counts["name"] == promoter]`).
    promoter : str
        Name of the promoter `gdf` belongs to (used to look up its WT sequence).

    Returns
    -------
    is_match : (N, L) float64 array
        1.0 where a barcode's sequence matches WT at that position, else 0.0.
        float64 (not int64/bool) so `is_match.T @ x` in `_mi_footprint` takes
        numpy's BLAS-backed matmul path instead of a slow generic one.
    x, y : (N,) float64 arrays
        Transformed (and optionally barcode-equalized) DNA/RNA weights.
    n_seqs : int
        Number of distinct mutated sequences observed for this promoter (unfiltered).
    """
    wt_seq = _load_wt_dict()[promoter]
    L = len(wt_seq)

    seqs = gdf["promoter"].astype(str)
    if (seqs.str.len() != L).any():
        raise ValueError(f"all sequences for {promoter} must have length {L}")
    seq_arr = np.frombuffer(
        "".join(seqs.values).encode("ascii"), dtype=np.uint8
    ).reshape(len(seqs), L)
    wt_arr = np.frombuffer(wt_seq.encode("ascii"), dtype=np.uint8)
    is_match = (seq_arr == wt_arr[None, :]).astype(np.float64)

    if require_both_counts:
        idx = (gdf["ct_0"].values > 0) & (gdf["ct_1"].values > 0)
    else:
        idx = np.ones(len(gdf), dtype=bool)

    is_match = is_match[idx, :]
    x = _transform_counts(gdf["ct_0"].values[idx], transform=transform)
    y = _transform_counts(gdf["ct_1"].values[idx], transform=transform)

    if barcode_equal:
        denom = x + y
        denom = np.where(denom > 0, denom, 1.0)
        x = x / denom
        y = y / denom

    n_seqs = gdf["promoter"].nunique()
    return is_match, x, y, n_seqs


def _mi_footprint(is_match, x, y, *, mismatch=None, a=0.0, normalize_by_entropy=False):
    """Core per-position MI computation from a WT-match matrix and DNA/RNA weight
    vectors (see `_prepare_promoter_arrays`). This is the only part that needs to
    be redone for each scramble, since scrambling just reorders/resamples which
    entries of `x`/`y` line up with the (fixed) rows of `is_match`.

    `mismatch` (1 - is_match) can be passed in precomputed when this is called
    repeatedly on the same `is_match` (e.g. once per scramble in
    `compute_MI_null_percentile`) to avoid rebuilding that (N, L) array every call;
    if omitted it's computed fresh, which is fine for a single call.

    Returns
    -------
    (mi, mi_norm, total, N) or None if the total read weight is <= 0.
    """
    L = is_match.shape[1]
    total = np.sum(x + y)
    if total <= 0:
        return None

    if mismatch is None:
        mismatch = 1.0 - is_match

    # Build joint weights per position. Indices: joint[pos, x_state, z]
    #   x_state: 0 = match WT, 1 = mismatch (mutated); z: 0 = DNA, 1 = RNA
    joint = np.zeros((L, 2, 2), dtype=np.float64)
    joint[:, 0, 0] = is_match.T @ x
    joint[:, 0, 1] = is_match.T @ y
    joint[:, 1, 0] = mismatch.T @ x
    joint[:, 1, 1] = mismatch.T @ y
    joint /= total

    mi = _mi_2x2_vectorized(joint)

    N = is_match.shape[0]
    mi = np.maximum(0.0, mi - a / N)

    if normalize_by_entropy:
        pz_dna = np.sum(x) / total
        Hz = float(_h2_vectorized(pz_dna))
        px_mut = joint[:, 1, 0] + joint[:, 1, 1]  # P(X=mut) per position
        Hx = _h2_vectorized(px_mut)
        denom_h = np.minimum(Hx, Hz)
        mi_norm = np.where(denom_h > 0, mi / denom_h, 0.0)
    else:
        mi_norm = None

    return mi, mi_norm, total, N


def compute_MI_from_counts(
    df_counts,
    ind=None,
    *,
    transform="sqrt",            # "raw"|"sqrt"|"log1p"|"pow:0.5"|"clip:200"
    barcode_equal=False,        # each barcode contributes total weight 1
    require_both_counts=False,   # filter to ct_0>0 & ct_1>0
    positions=None,             # override length/pos axis if needed
    normalize_by_entropy=False, # MI / min(H(X_i), H(Z))
    a=0.0,                      # background MI (in units of a/N) subtracted per position
):
    """Per-position mutual information between (mutation state) and (DNA vs. RNA),
    computed directly from an already-loaded counts dataframe.

    This is the shared core behind `compute_MI` (real counts) and
    `compute_MI_scrambled` (counts with ct_0/ct_1 scrambled via `scramble_counts`),
    so both produce information footprints through the exact same MI computation.
    Per-promoter setup (`_prepare_promoter_arrays`) and the MI computation itself
    (`_mi_footprint`) are both vectorized over numpy arrays rather than looping in
    Python over barcodes/positions.

    Parameters
    ----------
    df_counts : pandas.DataFrame
        Counts dataframe as returned by `get_counts` (or `get_counts_scrambled`/
        `scramble_counts`), with columns 'name', 'promoter', 'ct_0', 'ct_1'.
    ind : str | None
        Growth condition index, stored in the output for bookkeeping only.
    transform : str
        How count weights are reshaped before forming the joint (see `_transform_counts`).
    barcode_equal : bool
        If True, normalise each barcode’s (DNA, RNA) weights to sum to 1, so every
        barcode contributes equally regardless of total reads.
    require_both_counts : bool
        If True, drop barcodes with zero DNA or zero RNA reads.
    positions : array-like | None
        Override the position axis (defaults to 0..L-1).
    normalize_by_entropy : bool
        If True, also return MI normalised by min(H(X_i), H(Z)).
    a : float
        Background MI level (finite-sample bias), applied as `a / N` where `N` is the
        per-promoter barcode count (the 'N' column of the output). Subtracted from MI
        at every position, then clipped at 0 (`max(0, mi - a/N)`) since MI cannot be
        negative.

    Returns
    -------
    pandas.DataFrame with one row per (promoter, position).
    """
    df_out = []

    for promoter, gdf in df_counts.groupby("name"):
        if promoter in _EXCLUDED_PROMOTERS:
            continue

        is_match, x, y, n_seqs = _prepare_promoter_arrays(
            gdf, promoter,
            transform=transform, barcode_equal=barcode_equal, require_both_counts=require_both_counts,
        )
        result = _mi_footprint(is_match, x, y, a=a, normalize_by_entropy=normalize_by_entropy)
        if result is None:
            continue
        mi, mi_norm, total, N = result
        L = is_match.shape[1]

        # position axis
        if positions is None:
            pos = np.arange(L)
        else:
            if len(positions) != L:
                raise ValueError("positions length must match sequence length")
            pos = np.asarray(positions)

        out = pd.DataFrame({
            "promoter": promoter,
            "pos": pos,
            "mut_info": mi,
            "N": np.full(L, N, dtype=float),
            "N_reads": np.full(L, total, dtype=float),
            "N_seqs": np.full(L, n_seqs, dtype=float),
            "ind": ind,
            "transform": transform,
            "barcode_equal": barcode_equal,
            "require_both_counts": require_both_counts,
            "a": a,
        })
        if mi_norm is not None:
            out["mut_info_norm"] = mi_norm

        df_out.append(out)

    return pd.concat(df_out, ignore_index=True) if df_out else pd.DataFrame()


def compute_MI(ind, promoter_list=None, *, a=0.0, **kwargs):
    """Per-position mutual information between (mutation state) and (DNA vs. RNA),
    loading counts for sample `ind` via `get_counts`.

    See `compute_MI_from_counts` for the MI computation, the `a` background-
    subtraction parameter, and remaining keyword arguments (transform, barcode_equal,
    require_both_counts, positions, normalize_by_entropy).

    Parameters
    ----------
    ind : str
        Growth condition index (matches DNA/RNA barcode files).
    promoter_list : str | list[str] | None
        Restrict to one or more promoters; None → all promoters in `df_map`.
    a : float
        Background MI subtracted as `a / N` per position, clipped at 0.

    Returns
    -------
    pandas.DataFrame with one row per (promoter, position).
    """
    df_counts = get_counts(ind, promoter_list)
    return compute_MI_from_counts(df_counts, ind=ind, a=a, **kwargs)
