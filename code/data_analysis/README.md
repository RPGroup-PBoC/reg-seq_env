# Data analysis

Notebooks and helper modules for analysing the RegSeq2 barcode-count data: mutual-information / expression-shift footprints, *de novo* PWM discovery, comparison to RegulonDB-known TFs, and emerging-promoter prediction with the Salis-lab calculator.

## Notebooks

Run roughly in this order — each one consumes outputs from the ones above.

| File | Purpose |
| --- | --- |
| `sum_mapping.ipynb` | Summarise and filter the barcode → promoter mapping. Drops barcodes that map to multiple variants, computes per-promoter mutation rates, and compares the recovered pool against the original Twist order. |
| `compute_footprints.ipynb` | For every (condition, promoter), compute per-position mutual information (MI) and per-base expression shift from the DNA + RNA barcode counts. Writes `mi_all.csv` and `ex_all.csv`, plus footprint plots for selected promoters. Also fits the `a/N` background-MI scaling used downstream. |
| `compare_denovo_motifs.ipynb` | Build *de novo* PWMs from expression-shift footprints (one block each for salt, anaerobic, and TyrR conditions) and rank known E. coli TFs by similarity. Includes a genome-wide salt-regulon scan against PRECISE-1K fold changes and a yadI vs. CRP sanity check. |
| `de-novo_promoter.ipynb` | Two-stage pipeline. (1) Scan every (promoter, condition) pair with `denovo_scanner.run_scan` to find MI peaks above the `a/N` background that fall *outside* the canonical −55..0 promoter window, then aggregate recurrent peaks across conditions. (2) Salis-lab promoter-calculator analysis: for selected (promoter, gc, rep) targets, find the dominant expression-shift mutation and ask whether it creates / shifts a putative TSS. Exposes `compute_emerging_promoter_rates(promoter_names, gcs, reps, ...)`. |

## Helper modules

Imported by the notebooks above; not meant to be run on their own.

| File | Used by | Provides |
| --- | --- | --- |
| `denovo_scanner.py` | `de-novo_promoter.ipynb` | `run_scan`, `save_outputs`, `build_wt_dict`, `load_df_map` — the per-promoter MI peak caller. |
| `score_denovo_motif.py` | `compare_denovo_motifs.ipynb` | `run_analysis` — builds the de novo PWM from a footprint, scores it against every RegulonDB TF, computes a random-sequence null distribution. Also helpers like `build_known_pwms`, `extract_site`, `get_wt_sequence`. |
| `scan_genome_motif.py` | (helper for `compare_denovo_motifs.ipynb`) | `scan_genome`, `plot_hits` — slides a PWM across the MG1655 genome on both strands and returns positions above a percentile threshold. |

## Cached outputs

These CSVs are written by the notebooks and can be reloaded to skip the heavy computation. They are large enough that re-generating them takes minutes-to-hours.

| File | Written by | Contents |
| --- | --- | --- |
| `mi_all.csv` (~136 MB) | `compute_footprints.ipynb` | One row per (condition, promoter, position) with the per-position mutual information. |
| `ex_all.csv` (~237 MB) | `compute_footprints.ipynb` | One row per (condition, promoter, position, base) with the signed expression shift. |
| `emerging_promoter_tx_rates.csv` (~6 MB) | `de-novo_promoter.ipynb` | Salis-calculator Tx-rate predictions per putative TSS for each (promoter, gc, rep) target, with the constant library flanks. Columns: `promoter, gc, rep, pos, rates_w_mutation, rates_wo_mutation`. |
| `emerging_promoter_tx_rates_w_gen_flanks.csv` (~6 MB) | `de-novo_promoter.ipynb` | Same as above but using the genomic flanks looked up by `find_loc`. |
| `denovo_motif_comparison.csv` | `compare_denovo_motifs.ipynb` | Per-TF site / correlation scores against the de novo PWM. |
| `dgoR_predicted_TSS.pdf` | `de-novo_promoter.ipynb` | Example figure: predicted TSS rates for the dgoR promoter. |
