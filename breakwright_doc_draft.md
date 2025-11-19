# Breakwright

Breakwright is a three–stage workflow for **detecting, annotating, and visually vetting candidate misassemblies** in long-read assemblies.

It ties together:

- contig→reference PAF alignments,
- HiFiASM GFA graphs,
- read alignments back to the assembly (BAM),
- optional HiFiASM `lowQ.bed` tracks,

to give you **ranked break candidates** plus IGV-style figures for human inspection.

---

## Overview

The workflow has three scripts:

1. **`breakwright_paf_breakfinder.py`**  
   From contigs→reference PAF, proposes candidate breakpoints and reasons.

2. **`breakwright_gfa_annotator.py`**  
   Annotates those candidates with HiFiASM GFA evidence and exports small graph subgraphs.

3. **`breakwright_viz_plus.py`**  
   Uses read alignments, optional `lowQ.bed`, and optional contigs→reference PAF to:
   - compute coverage/read-support metrics,
   - classify strict spanning reads in reference space,
   - generate IGV-style panels around each break,
   - and write an augmented TSV with all support scores.

An example panel (output of `breakwright_viz_plus.py`) is shown below:

![Example Breakwright panel](/mnt/data/55ced7d9-3d3e-4ac9-ad13-81bbcc498fa4.png)

---

## 0. Inputs & Dependencies

### Required inputs

- **Assembly contigs**: from HiFiASM (or similar).
- **Reads**: long reads mapped back to contigs (`.bam + .bai`).
- **Reference genome**: for contig→reference mapping.
- **HiFiASM GFA**: e.g. `*.p_ctg.gfa`.

### Python dependencies

- Python ≥ 3.8
- `pysam`
- `numpy`
- `matplotlib`

Install via:

```bash
pip install pysam numpy matplotlib
````

---

## 1. `breakwright_paf_breakfinder.py`

### Purpose

Call **candidate misassembly breakpoints** from contig→reference PAF alignments.

### Input

PAF from `minimap2`:

```bash
minimap2 -x asm5 -t 32 ref.fa contigs.fa > contigs_vs_ref.paf
```

### Basic usage

```bash
python breakwright_paf_breakfinder.py \
  --paf contigs_vs_ref.paf \
  --outprefix breaks/soy \
  --min_mapq 20 \
  --min_aln_len 5000 \
  --min_identity 0.90 \
  --min_qgap 10000 \
  --min_tjump 100000 \
  --allow_overlap 1000 \
  --max_micro_overlap 5000 \
  --identity_drop 0.10 \
  --low_mapq_edge 30 \
  --min_tail_unmapped 20000 \
  --max_merge_dist 10000
```

### What it does (conceptually)

1. **Parse PAF into blocks per contig**
   Retain only high-quality blocks (`min_mapq`, `min_aln_len`, `min_identity`).
   Each block stores query span, reference span, strand, MAPQ, identity.

2. **Mark repetitive regions on the query**
   For each pair of blocks on the same contig:

   * if they overlap on the query by ≥ 80% of the shorter block **and**

     * map to different reference sequences, or
     * map far apart on the same sequence (≥ `min_tjump/2` bp),
       then both are marked `repetitive`, and their partner hits are recorded.

3. **Identify unmapped tails**
   If the leading or trailing unmapped region on a contig is ≥ `min_tail_unmapped` bp:

   * record `unmapped_lead` or `unmapped_tail` break.

4. **Examine adjacent blocks in query order**
   For each adjacent pair `(a,b)` on a contig, accumulate **reasons** that a break might exist at `a.qend`:

   * Reference behaviour

     * `switch_chr` – `a.tname != b.tname`
     * `strand_flip` – mapping orientation flips
     * `large_tjump` – midpoints on same reference differ by ≥ `min_tjump`

   * Query behaviour

     * `large_qgap` – gap on the contig between blocks ≥ `min_qgap`
     * `micro_overlap` – overlap on the contig exceeds `allow_overlap` but ≤ `max_micro_overlap`

   * Alignment quality changes

     * `identity_drop` – identity drops by ≥ `identity_drop`
     * `low_mapq_edge` – either block has MAPQ ≤ `low_mapq_edge`

   * Repeats

     * `repetitive_region` – either block overlaps a repetitive region

5. **Merge nearby cuts**
   Raw candidate cuts are merged if they’re within `max_merge_dist` bp on the contig.
   Merging unions:

   * `reason` terms,
   * `repetitive_region` flags,
   * and `repetitive_matches` partner metadata.

### Outputs

1. **Candidate breaks**

`<outprefix>_breaks.tsv`

Columns (per break):

* `qname` – contig ID
* `cut` – proposed break position (bp, contig coordinate)
* `reason` – comma-separated list of reasons
* `qlen` – contig length
* `qbeg`, `qend` – local query span between the blocks
* `tname1`, `tpos1` – reference locus for block A (midpoint)
* `tname2`, `tpos2` – reference locus for block B (midpoint)
* `repetitive_region` – yes/no
* `repetitive_matches` – partner mappings in the form `tname:start-end@qQSTART-QEND`

2. **Reason summary**

`<outprefix>_summary.tsv`

Two columns:

* `reason`
* `count` (how many breaks include that reason)

3. **Optional block dump**

If `--emit_all_blocks` is set:

`<outprefix>_blocks.tsv` – QC view of all retained PAF blocks.

### How to read it

* Breaks with **multiple structural reasons** (e.g. `switch_chr,large_qgap,identity_drop`) are prime misassembly candidates.
* Breaks flagged only as `repetitive_region` may simply reflect repeats and will often require GFA/read support to decide.

---

## 2. `breakwright_gfa_annotator.py`

### Purpose

Add **HiFiASM graph context** to each break candidate and export small subgraphs for manual inspection.

### Basic usage

```bash
python breakwright_gfa_annotator.py \
  --gfa asm.p_ctg.gfa \
  --breaks breaks/soy_breaks.tsv \
  --outprefix breaks/soy \
  --end_proximity_bp 10000 \
  --min_overlap_warn 50 \
  --subgraph_hops 2 \
  --subgraph_min_overlap 0 \
  --emit_subgraph_gfa
```

### What it does

1. **Parse GFA**

Build per-segment information:

* `seg_len[seg]` – length from sequence or `LN:i:` tag.
* `deg[seg]['left'|'right']` – number of links at each end.
* `ovls[seg]['left'|'right']` – list of overlap sizes for links on each end.
* `adj[seg]` – undirected adjacency list (segment, overlap_bp) for BFS.
* `s_lines`, `l_lines` – raw S/L lines (for mini-GFA export).

2. **Annotate breaks**

For each break:

* If `qname` not in GFA, send row to `_unmatched.tsv`.
* Otherwise compute:

  * `qlen` – contig length from GFA.
  * `dist_to_start`, `dist_to_end` – distance from cut to each end.
  * `deg_left`, `deg_right` – degrees at ends.
  * `min_ovl_left`, `min_ovl_right` – smallest overlap at each end (or `NA`).
  * `nearest_end` – `"left"` or `"right"` (closest end).
  * `nearest_end_degree` – degree at nearest end.
  * `nearest_junction_bp` – distance from cut to the **nearest end with degree≠1** (or `NA`).

  Then classify `gfa_flag`:

  * `junction_near_end` – nearest end within `end_proximity_bp` and degree > 1.
  * `tip_near_end`      – nearest end within `end_proximity_bp` and degree = 0.
  * `weak_overlap_end`  – overlaps at an end < `min_overlap_warn` bp (and no other flag).
  * `simple`            – none of the above.

3. **Export subgraphs**

For each annotated break, run a BFS starting at `qname`:

* up to `subgraph_hops` link steps,
* traversing only edges with `overlap_bp ≥ subgraph_min_overlap`.

Write:

* `<subgraph_dir>/<qname>_<cut>.txt` – enriched text dump:

  * #nodes, #edges, degree stats, junction list (deg≥3), node and edge lists.
* (optionally) `<subgraph_dir>/<qname>_<cut>.gfa` – mini-GFA with only relevant `S`/`L` lines.

Also write an index:

`<outprefix>_subgraphs.tsv` with:

* `qname`, `cut`, `n_nodes`, `hops`, `min_overlap`,
* `node_file`, `subgraph_gfa`, `nodes_csv`.

### Outputs

1. **Annotated breaks**

`<outprefix>_breaks_gfa.tsv` – main output for the next stage.

Adds:

* `qlen`, `dist_to_start`, `dist_to_end`
* `deg_left`, `deg_right`
* `min_ovl_left`, `min_ovl_right`
* `nearest_end`, `nearest_end_degree`
* `nearest_junction_bp`
* `gfa_flag`

2. **Unmatched breaks**

`<outprefix>_unmatched.tsv` – breaks whose `qname` is missing from the GFA (unless `--keep_only_matched`).

3. **Subgraphs & index**

As described above.

### How to read it

* `junction_near_end` with a small `nearest_junction_bp` means the contig end sits near a graph junction and is structurally interesting.
* `tip_near_end` suggests a dangling contig end near the break.
* `weak_overlap_end` highlights ends that are held together only by very short overlaps (potentially fragile joins).
* Graph subgraphs help distinguish simple paths from complex repeat tangles around the break.

---

## 3. `breakwright_viz_plus.py`

### Purpose

For each GFA-annotated break, generate:

* quantitative read/coverage metrics,
* optional `lowQ.bed`/GFA support scores,
* optional contig→reference spanning-read classification,
* IGV-style panels,
* a final **augmented breaks table**.

### Basic usage

```bash
python breakwright_viz_plus.py \
  --bam reads_to_contigs.bam \
  --breaks breaks/soy_breaks_gfa.tsv \
  --outdir viz \
  --window 25000 \
  --min_mapq 10 \
  --hard_min_mapq_for_panel 20 \
  --max_reads 10000 \
  --min_aln_len 500 \
  --dpi 300 \
  --gfa_keep_flags junction_near_end,weak_overlap_end \
  --gfa_max_nearest_junction_bp 20000 \
  --lowq_bed asm.lowQ.bed \
  --paf contigs_vs_ref.paf \
  --span_flank_bp 2000 \
  --ref_span_local_bp 50000
```

### What it does

1. **Filter breaks**

* Keep only rows whose `gfa_flag` is in `--gfa_keep_flags`, if provided.
* Keep only rows with `nearest_junction_bp ≤ --gfa_max_nearest_junction_bp`, if provided.
* Optional `--limit` to cap number of breaks rendered.

2. **Coverage & reads**

For each break:

* Define `[start,end] = [cut-window, cut+window]`, clamped to contig length.
* Using `bam.fetch`:

  * keep primary, non-duplicate, non-QC-fail reads,
  * with `MAPQ ≥ max(min_mapq, hard_min_mapq_for_panel)`,
  * `aligned_len ≥ min_aln_len`,
  * random subsampling to `max_reads` if needed.
* Build coverage array `cov[start:end]`.

3. **Read support metrics**

From `cov` and `all_reads`, compute:

* `n_reads_window`
* `n_reads_cross_cut` – reads overlapping the cut
* `n_reads_span_strict` – reads that strictly span `cut±span_flank_bp`
* `cov_at_cut`, `cov_median`
* `cov_adj_cut`, `cov_adj_dist` – position and distance of local coverage trough (if present)
* `cov_flag` – one of:

  * `no_reads`
  * `no_support_at_cut`
  * `low_cov_trough`
  * `ok`
* `read_support_score` (0–2):

  * 0 if no reads / no support at cut
  * 1 for normal coverage across cut
  * 2 if strong coverage trough

4. **LowQ bed support (optional)**

If `--lowq_bed` is provided:

* `lowq_overlap` – `"yes"` if any lowQ interval overlaps the window.
* `lowq_at_cut` – `"yes"` if lowQ interval covers the cut.
* `lowq_support` – 1 if `lowq_at_cut=="yes"`, else 0.

5. **GFA support score**

Combine existing GFA annotations:

* 2 points if `gfa_flag == junction_near_end` and `nearest_junction_bp ≤ threshold`.
* 1 point if `gfa_flag ∈ {weak_overlap_end, tip_near_end}` or `deg_left`/`deg_right` ≠ 1.
* 0 otherwise.

6. **Total support score**

```text
support_score_total = gfa_support_score
                      + read_support_score
                      + lowq_support
```

7. **Reference-span classification (optional)**

If `--paf` is provided:

* Build contig→reference blocks (`parse_paf_blocks`).
* For each strict spanning read:

  * lift positions at `cut±span_flank_bp` into reference space (`lift_pos_to_ref`).
  * classify each read:

    * `diff_chr`, `same_chr_far`, `same_chr_local`, `unmapped`, or `no_paf`.
* Aggregate:

  * `n_span_ref_same_chr_local`
  * `n_span_ref_same_chr_far`
  * `n_span_ref_diff_chr`
  * `n_span_ref_unmapped`
  * `span_ref_class` summarizing the pattern (e.g. `reads_span_diff_chr`).

If no PAF is provided but you have strict spanners, `span_ref_class = "strict_spanners_no_paf"`.

8. **QC notes**

* `qc_note = "contig_not_in_bam"` – contig missing from BAM header.
* `qc_note = "fetch_error"` – BAM fetch failure.
* `qc_note = "no_reads_in_window"` – window is empty.
  Only the first two skip plotting; `no_reads_in_window` can still produce an empty panel if desired.

### Outputs

For each break (unless plotting is skipped):

1. **Figures**

* `<outdir>/<prefix>_<qname>_<cut>.pdf`
* `<outdir>/<prefix>_<qname>_<cut>.png`

2. **Read list**

* `<outdir>/<prefix>_<qname>_<cut>.reads.txt` with:

  ```text
  #lane  read_name  MAPQ  alen  start  end  span_ref_class
  ```

3. **Augmented break table**

* `<outdir>/<breaks_basename>_reads.tsv`

  (Note: this file is intentionally written to `--outdir`.)

Adds the following columns (if not already present):

```text
cov_flag
n_reads_window
n_reads_cross_cut
n_reads_span_strict
cov_at_cut
cov_median
cov_adj_cut
cov_adj_dist
read_support_score
gfa_support_score
lowq_overlap
lowq_at_cut
lowq_support
support_score_total
qc_note
n_span_ref_same_chr_local
n_span_ref_same_chr_far
n_span_ref_diff_chr
n_span_ref_unmapped
span_ref_class
```

You can sort/filter this table to prioritize which panels to inspect (e.g. high `support_score_total`, `span_ref_class == "reads_span_diff_chr"`).

### How to read the panel

Using the example above:

* **Top panel** – coverage vs contig coordinate.

  * Black dashed line: original `cut`.
  * Red dashed line: `cov_adj_cut` (coverage trough) if different.
* **Bottom panel** – top 50 reads, ranked by MAPQ and length.

  * Horizontal bars = individual reads.
  * Bar shade = MAPQ (darker is higher), color scale shown on the right.
  * Border colors for strict spanning reads:

    * Red: `diff_chr`
    * Orange: `same_chr_far`
    * Blue: `same_chr_local`
    * Black: `unmapped`/`no_paf`
  * Vertical colored ticks inside bands = SNPs vs a pseudo-reference (A/C/G/T), revealing haplotype patterns.
* **Annotation box** – textual summary:

  * `gfa_flag`, degrees (`L:deg_left R:deg_right`), `nearest_junction_bp`
  * original `reason` from PAF stage
  * `cov_flag`
  * `S:total (G:gfa,R:read,Q:lowQ)` support scores
  * `reads:cross`, `reads:strict`
  * `ref_span:span_ref_class (diff=…, far=…, local=…, unmap=…)`
  * `QC:` note if present.

A “classic” misassembly signature is:

* `cov_flag = low_cov_trough`,
* many strict spanning reads,
* `span_ref_class = reads_span_diff_chr` or `reads_span_same_chr_far`,
* `gfa_flag = junction_near_end` or `weak_overlap_end`,
* high `support_score_total`.

---

## 4. Suggested pipeline

A minimal end-to-end run might look like:

```bash
# 1) Call candidate breaks from PAF
python breakwright_paf_breakfinder.py \
  --paf contigs_vs_ref.paf \
  --outprefix breaks/soy

# 2) Annotate with HiFiASM GFA and export subgraphs
python breakwright_gfa_annotator.py \
  --gfa asm.p_ctg.gfa \
  --breaks breaks/soy_breaks.tsv \
  --outprefix breaks/soy \
  --emit_subgraph_gfa

# 3) Visualize breaks with read coverage and ref-span classification
python breakwright_viz_plus.py \
  --bam reads_to_contigs.bam \
  --breaks breaks/soy_breaks_gfa.tsv \
  --outdir viz \
  --paf contigs_vs_ref.paf \
  --lowq_bed asm.lowQ.bed \
  --gfa_keep_flags junction_near_end,weak_overlap_end \
  --gfa_max_nearest_junction_bp 20000
```

From here you can:

* review `viz/*pdf`/`*png` panels,
* sort `viz/soy_breaks_gfa_reads.tsv` by `support_score_total`,
* inspect GFA subgraphs in Bandage using the per-break mini-GFAs.

---

## 5. Citation

If you use Breakwright in a publication, please cite:

> Conrad & collaborators, Breakwright (in prep). A practical framework for identifying and validating misassemblies in long-read HiFiASM genomes.
