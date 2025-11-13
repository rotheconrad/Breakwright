![Breakwright Logo](images/breakwright_logo.png)

# BREAKWRIGHT: Precision-crafted genome correction.

This repository provides a reproducible workflow for detecting and correcting structural misassemblies in **genome assemblies** using a reference genome. It was developed using PacBio HiFi reads and assembled contigs from HiFiASM.

The pipeline filters low-value contigs, identifies potential misjoins (e.g., fused telomeres), visualizes breakpoints, and produces a corrected FASTA. It gives the user data to make informed curation decisions about which positions in assembled contigs to break.

---

## 🧬 Overview

This pipeline performs:

1. **Genome Assembly** - Assemble PacBio HiFi reads with HiFiASM. *note: this pipeline can be used for other genome assemblies as well but you may have to skip the gfa step*
    ```bash
    hifiasm -o sample hifi_reads.fa 
    ```
2. **Contig to reference alignment** - align assembled contigs to reference genome using minimap2 for PAF output file. Used to find proposed break points
    ```bash
    minimap2 -x asm5 ref.fa sample.bp.p_ctg.fa > contigs_vs_ref.paf
    ```
3. **Contig filtering** — remove short, redundant, or unmapped contigs based on reference alignment coverage (`01_breakwright_contig_paf_filter.py`).
    ```bash
    python 01_breakwright_contig_paf_filter.py \
        --assembly_fa sample.bp.p_ctg.fa \
        --contigs_vs_ref_paf contigs_vs_ref.paf \
        --outprefix sample_assembly_filtered.fa
    ```
4. **Unassembled reads to assembled contigs alignment** — sam to bam output from `minimap2`. Used to verify proposed break points.
    ```bash
    minimap2 -ax map-hifi sample_assembly_filtered.fa hifi_reads.fastq.gz > hifi_reads_to_contigs.sam
    samtools view -bS hifi_reads_to_contigs.sam | samtools sort -o hifi_reads_to_contigs.bam
    samtools index hifi_reads_to_contigs.bam
    ```
5. **Breakpoint detection** — identify possible chimeric joins using PAF-based inspection (`02_breakwright_paf_breakfinder.py`).
    ```bash
    python 02_breakwright_paf_breakfinder.py \
        --paf contigs_vs_ref.paf \
        --outprefix assembly_breaks.tsv
    ```
6. **Breakpoint verification** - Use the GFA file from the HiFiASM output to incorporate additional evidence for the proposed break points (`03_breakwright_gfa_annotator.py`).
    ```bash
    python 03_breakwright_gfa_annotator.py \
        --gfa sample.bp.p_ctg.gfa \
        --breaks assembly_breaks.tsv \
        --outprefix assembly_breaks_gfa.tsv
    ```
7. **Breakpoint visualization** — plot local alignment structure and supporting read coverage (`04_breakwright_dotplot.py` and `05_breakwright_viz_plus.py`).
    ```bash
    python 04_breakwright_dotplot.py \
        --paf contigs_vs_ref.paf \
        --breaks assembly_breaks_gfa.tsv \
        --outdir dotplots --mode full,per-chr --draw_chr_ticks --export_png
    
    python 05_breakwright_viz_plus.py \
        --bam hifi_reads_to_contigs.bam \
        --breaks assembly_breaks_gfa.tsv \
        --outdir viz --export_png
    ```
8. **Assembly correction** — split contigs at curated breakpoints (`06_breakwright_split_breaks.py`).
    ```bash
    python 06_breakwright_split_breaks.py \
      --assembly_fa sample_assembly_filtered.fa \
      --breaks_tsv curated_breaks.tsv \
      --out_fa sample_assembly_filtered_corrected.fa
    ```

---

## 🧩 Dependencies

All tools are open-source and installable via `conda`:

- `python>=3.9`
- `numpy`, `pandas`, `matplotlib`, `pysam`
- `hifiasm`
- `minimap2`
- `samtools`
- `seqkit` (for assembly stats)

---

## ⚙️ Create Conda Environment

```bash
conda create -p ./envs/binf python=3.10 -y
conda activate ./envs/binf

# Install dependencies
conda install -c bioconda hifiasm minimap2 seqkit samtools -y
pip install numpy pandas matplotlib
```

---

## 🧬 Step 1. Assemble genome with HiFiASM

Use `hifiasm` to assemble contigs.

```bash
hifiasm -o sample hifi_reads.fa
```

---

## 🧭 Step 2. Map contigs to the reference genome

Use `minimap2` to align assembled contigs to a high-quality reference genome (e.g., *Glycine max* Williams 82).

```bash
minimap2 -x asm5 -t 32 ref.fa sample.bp.p_ctg.fa > contigs_vs_ref.paf
```

---

## 🔍 Step 3. Filter contigs with `01_breakwright_contig_paf_filter.py`

Removes short contigs that are redundant or unmapped based on their novelty of reference coverage.

#### 🧠 Overview
This script filters contigs from an assembly FASTA based on their alignments to a reference genome (in PAF format).
It removes short or redundant contigs that either don't map or map to regions already covered by longer contigs.

```bash
01_breakwright_contig_paf_filter.py \
  --assembly_fa sample.bp.p_ctg.fa \
  --contigs_vs_ref_paf contigs_vs_ref.paf \
  --outprefix sample_assembly_filtered.fa \
  --min_len 10000 \
  --min_mapq 20 \
  --min_aln_len 5000 \
  --min_identity 0.9 \
  --novel_bp_thresh 10000 \
  --novel_frac_thresh 0.25
```

#### 💻 Example Command

```bash
01_breakwright_contig_paf_filter.py   --assembly_fa sample.bp.p_ctg.fa   --contigs_vs_ref_paf contigs_vs_ref.paf   --outprefix sample_assembly_filtered.fa   --min_len 10000   --min_mapq 20   --min_aln_len 5000   --min_identity 0.9   --novel_bp_thresh 10000   --novel_frac_thresh 0.25   --mode greedy
```

#### 📤 Outputs

| Output File | Description |
|--------------|-------------|
| `<outprefix>.kept.fa` | FASTA file of all retained contigs. |
| `<outprefix>.dropped.list` | List of dropped contig IDs (one per line). |
| `<outprefix>.decision.tsv` | Table with columns: `contig`, `length`, `status`, `novel_bp`, `total_aligned_bp`, and `decision`. |
| `<outprefix>.kept_coverage.bed` | BED file showing merged reference coverage contributed by kept contigs. Useful for IGV or `bedtools` visualization. |

#### ⚙️ Required Arguments

| Argument | Type | Description |
|-----------|------|-------------|
| `--assembly_fa` | *string (path)* | **Required.** Path to the input contig FASTA file (can be `.gz`). Used to read contigs and calculate lengths. |
| `--contigs_vs_ref_paf` | *string (path)* | **Required.** PAF file from `minimap2 -x asm5` mapping contigs to the reference genome. Used to assess coverage and redundancy. |
| `--outprefix` | *string (path prefix)* | **Required.** Prefix for output files. All generated files will use this as a base (e.g., `outprefix.kept.fa`, `outprefix.dropped.list`, etc.). |

#### 🧩 Optional Arguments (with Defaults)

| Argument | Default | Description |
|-----------|----------|-------------|
| `--min_len` | `10000` | Contigs shorter than this (in bp) are considered “short.” Short contigs must add novel reference coverage to be retained. |
| `--min_mapq` | `20` | Minimum mapping quality (MAPQ) required for an alignment block to be considered valid. |
| `--min_aln_len` | `5000` | Minimum alignment length (in bp). Shorter alignments are ignored. |
| `--min_identity` | `0.90` | Minimum identity (nmatch/alignment length). Alignments below this threshold are discarded. |
| `--novel_bp_thresh` | `10000` | Minimum number of **novel reference base pairs** (not already covered by longer contigs) needed for a short contig to be retained. |
| `--novel_frac_thresh` | `0.25` | Fraction of aligned bases that must be novel for a short contig to be kept. |
| `--mode` | `"greedy"` | Determines contig evaluation order. Options: `"greedy"` (by contig length) or `"score"` (by aligned bp). |

#### 🧮 Suggested Thresholds

- **10 kb** → conservative default for removing tiny fragments.  
- **5 kb** → lenient mode, keeps shorter but potentially unique contigs.  
- **1 kb** → only if you want to preserve very small unique segments (e.g., organellar inserts).  
  Consider raising `--novel_bp_thresh` in that case.

#### 🔍 Notes

- If PAF lacks reliable identities, rely on `--min_mapq` and `--min_aln_len` thresholds.  
- Use `--mode score` to rank contigs by total aligned bp instead of contig length.  
- The BED output helps visualize coverage vs. the reference genome.

---

## 🧬 Step 4. Map HiFi reads to the filtered contigs

This step is used to confirm breakpoints and support manual curation.

```bash
minimap2 -ax map-hifi -t 32 sample_assembly_filtered.fa hifi_reads.fastq.gz > hifi_reads_to_contigs.sam
samtools view -bS hifi_reads_to_contigs.sam | samtools sort -o hifi_reads_to_contigs.bam
samtools index hifi_reads_to_contigs.bam
```
---

## 🔧 Step 5. Identify candidate breaks with `paf_breakfinder.py`

This script parses the PAF of contigs aligned to the reference to detect structural inconsistencies (e.g., multi-chromosome mappings, inversions, large internal gaps).

#### 🧠 Overview
`paf_breakfinder.py` scans a **PAF** file of **contigs → reference** alignments to identify **candidate misassembly breakpoints**.
It detects hallmark patterns such as **reference chromosome switches**, **strand flips**, **large intrachromosomal jumps**, and **large internal gaps on the contig**.
It also flags **micro-overlaps**, **identity drops**, **edge low-MAPQ blocks**, and **unmapped leading/trailing contig tails**.

```bash
python 02_breakwright_paf_breakfinder.py \
  --paf contigs_vs_ref.paf \
  --outprefix assembly_breaks.tsv
```

#### 💻 Example

```bash
python 02_breakwright_paf_breakfinder.py   --paf contigs_vs_ref.paf   --outprefix assembly_breaks.tsv  --min_mapq 20 --min_aln_len 5000 --min_identity 0.9   --min_qgap 10000 --min_tjump 100000 --allow_overlap 1000   --max_micro_overlap 5000 --identity_drop 0.10 --low_mapq_edge 30   --min_tail_unmapped 20000 --max_merge_dist 10000
```

#### 📤 Outputs

| Output | Description |
|---|---|
| `<outprefix>_breaks.tsv` | Candidate breakpoints with columns: `qname, cut, reason, qlen, qbeg, qend, tname1, tpos1, tname2, tpos2`. |
| `<outprefix>_summary.tsv` | Per-contig summary of counts by reason (`switch_chr`, `strand_flip`, `large_tjump`, `large_qgap`, `micro_overlap`, `identity_drop`, `low_mapq_edge`, `unmapped_lead`, `unmapped_tail`). |
| `<outprefix>_blocks.tsv` | *(optional; if `--emit_all_blocks`)* Filtered alignment blocks retained after thresholds. |

#### ⚙️ Required Arguments

| Argument | Type | Description |
|---|---|---|
| `--paf` | *string (path)* | **Required.** PAF file from `minimap2 -x asm5 ref.fa contigs.fa` (contigs→reference). |
| `--outprefix` | *string (path prefix)* | **Required.** Prefix for outputs (e.g., `breaks/soy`). |

#### 🧩 Optional Arguments (with Defaults)

| Argument | Default | Description |
|---|---|---|
| `--min_mapq` | `20` | Minimum MAPQ to accept a PAF alignment block. |
| `--min_aln_len` | `5000` | Minimum alignment length (bp); shorter blocks are ignored. |
| `--min_identity` | `0.90` | Minimum identity (nmatch/alen) to retain an alignment block. |
| `--min_qgap` | `10000` | Minimum **contig gap** (bp) between adjacent blocks to flag a candidate break (`large_qgap`). |
| `--min_tjump` | `100000` | Minimum **reference jump** (bp) across adjacent blocks on same chromosome to flag (`large_tjump`). |
| `--allow_overlap` | `1000` | Allow up to this many bp of query overlap between adjacent blocks when evaluating gaps. |
| `--require_ordered` | `True` | If set, only consider adjacency in ascending query order for breakpoint inference. |
| `--max_merge_dist` | `10000` | Merge nearby breakpoints on the same contig if within this many bp into a single candidate. |
| `--emit_all_blocks` | `False` | If set, also emits a TSV with all filtered alignment blocks for QC. |
| **New:** `--max_micro_overlap` | `5000` | Flag `micro_overlap` if adjacent blocks overlap by **(allow_overlap, max_micro_overlap]** bp. |
| **New:** `--identity_drop` | `0.10` | Flag `identity_drop` if identity decreases by ≥ this fraction between adjacent blocks (e.g., 0.10 = 10%). |
| **New:** `--low_mapq_edge` | `30` | Flag `low_mapq_edge` if either adjacent block has MAPQ ≤ this value (even if ≥ `--min_mapq`). |
| **New:** `--min_tail_unmapped` | `20000` | Flag `unmapped_lead` / `unmapped_tail` if leading/trailing unmapped contig tails exceed this many bp. |

#### 🔍 Detection Heuristics (adjacent blocks on the same contig)
- **switch_chr** — consecutive blocks map to **different reference chromosomes**.  
- **strand_flip** — strand changes between adjacent blocks on the same chromosome.  
- **large_tjump** — same chromosome but reference midpoints jump by `≥ --min_tjump`.  
- **large_qgap** — contig gap (`next.qstart - prev.qend`) ≥ `--min_qgap` after allowing `--allow_overlap`.  
- **micro_overlap** — query overlap **in (allow_overlap, max_micro_overlap]** bp (suspicious micro-overlap).  
- **identity_drop** — drop in identity (prev.ident - next.ident) ≥ `--identity_drop`.  
- **low_mapq_edge** — either block has MAPQ ≤ `--low_mapq_edge` (despite passing `--min_mapq`).  
- **unmapped_lead / unmapped_tail** — leading (`first.qstart`) or trailing (`qlen - last.qend`) unmapped contig tails ≥ `--min_tail_unmapped`.

> Heuristics are conservative defaults; tune for organism and assembly specifics.

#### 🧠 Notes
- Blocks are sorted per contig by **query start** to infer adjacency.  
- Multiple reasons may be concatenated for a merged breakpoint window.  
- Use together with read support (HiFi mappings) and visualization for final curation.

#### Example output of breaks.tsv

| qname      | cut      | reason        | qlen     | qbeg    | qend     | tname1      | tpos1    | tname2      | tpos2    |
| ---------- | -------- | ------------- | -------- | ------- | -------- | ----------- | -------- | ----------- | -------- |
| ptg000004l | 41418232 | unmapped_lead | 44724916 | 0       | 41418233 | NC_038241.2 | 1651120  | NC_038241.2 | 1651120  |
| ptg000007l | 61937    | unmapped_lead | 5461000  | 0       | 61938    | NC_038251.2 | 19627494 | NC_038251.2 | 19627494 |
| ptg000007l | 353365   | unmapped_tail | 5461000  | 353365  | 5461000  | NC_038251.2 | 19637626 | NC_038251.2 | 19637626 |
| ptg000007l | 77458    | large_qgap    | 5461000  | 70277   | 106787   | NC_038251.2 | 19639015 | NC_038251.2 | 19634567 |
| ptg000007l | 127865   | large_qgap    | 5461000  | 112266  | 177907   | NC_038251.2 | 19640635 | NC_038251.2 | 19636409 |
| ptg000007l | 177907   | large_qgap    | 5461000  | 154954  | 218336   | NC_038251.2 | 19636409 | NC_038251.2 | 19640715 |
| ptg000007l | 242879   | large_qgap    | 5461000  | 234763  | 276006   | NC_038251.2 | 19633745 | NC_038251.2 | 19638869 |
| ptg000007l | 324101   | large_qgap    | 5461000  | 307828  | 353365   | NC_038251.2 | 19639711 | NC_038251.2 | 19637626 |
| ptg000011l | 8222744  | unmapped_tail | 31842626 | 8222744 | 31842626 | NC_038255.2 | 43978809 | NC_038255.2 | 43978809 |

---

## 📊 Step 6. Verify breakpoints with `03_breakwright_gfa_annotator.py`

#### Purpose
Annotate curated breakpoints with **assembly-graph context** from a HiFiASM GFA:
- Flag likely graph-supported breakpoints (e.g., nearby junctions, tips).
- Export **Bandage-friendly subgraphs** around each break:
  - Plain-text **subgraph summary** (`*.txt`)
  - **Mini-GFA** fragment (`*.gfa`) — contains only `S` and `L` records for the local neighborhood.

#### Example
```bash
python 03_breakwright_gfa_annotator.py   --gfa hifiasm.asm.p_ctg.gfa   --breaks breaks/soy_breaks.tsv   --outprefix breaks/soy_breaks_gfa   --subgraph_hops 2   --subgraph_min_overlap 100   --emit_subgraph_gfa
```

#### Outputs
- `<outprefix>_breaks_gfa.tsv` — breaks with minimal GFA columns appended.
- `<outprefix>_subgraphs.tsv` — index of per-break stats and file paths.
- `<outprefix>_subgraphs/<qname>_<cut>.txt` — **detailed subgraph summary** (nodes/edges/degree/junctions).
- `<outprefix>_subgraphs/<qname>_<cut>.gfa` — **mini-GFA** fragment (if `--emit_subgraph_gfa`).

#### Required

| Argument | Description |
|---|---|
| `--gfa` | HiFiASM `.gfa` (e.g., `hifiasm.asm.p_ctg.gfa`) |
| `--breaks` | TSV with `qname` and `cut` (may include other cols) |
| `--outprefix` | Prefix for output files |

#### Optional (defaults)

| Argument | Default | Description |
|---|---|---|
| `--subgraph_hops` | `2` | BFS depth in `L` links from the contig node |
| `--subgraph_min_overlap` | `0` | Traverse links with numeric overlap ≥ this many bp (if parseable) |
| `--subgraph_dir` | `<outprefix>_subgraphs/` | Directory for per-break artifacts |
| `--emit_subgraph_gfa` | *flag* | If set, write `<contig>_<cut>.gfa` with only the relevant `S`/`L` lines |
| `--sep` | `\t` | Breaks TSV delimiter |

#### Example Output of breaks_gfa.tsv

| qname | cut | reason | qlen | qbeg | qend | tname1 | tpos1 | tname2 | tpos2 | dist_to_start | dist_to_end | deg_left | deg_right | min_ovl_left | min_ovl_right | nearest_end | nearest_end_degree | nearest_junction_bp | gfa_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ptg000004l | 41418232 | unmapped_lead | 44724916 | 0 | 41418233 | NC_038241.2 | 1651120 | NC_038241.2 | 1651120 | 41418232 | 3306684 | 0 | 0 | NA | NA | right | 0 | 3306684 | simple |
| ptg000007l | 61937 | unmapped_lead | 5461000 | 0 | 61938 | NC_038251.2 | 19627494 | NC_038251.2 | 19627494 | 61937 | 5399063 | 0 | 0 | NA | NA | left | 0 | 61937 | simple |
| ptg000007l | 353365 | unmapped_tail | 5461000 | 353365 | 5461000 | NC_038251.2 | 19637626 | NC_038251.2 | 19637626 | 353365 | 5107635 | 0 | 0 | NA | NA | left | 0 | 353365 | simple |
| ptg000007l | 77458 | large_qgap | 5461000 | 70277 | 106787 | NC_038251.2 | 19639015 | NC_038251.2 | 19634567 | 77458 | 5383542 | 0 | 0 | NA | NA | left | 0 | 77458 | simple |
| ptg000007l | 127865 | large_qgap | 5461000 | 112266 | 177907 | NC_038251.2 | 19640635 | NC_038251.2 | 19636409 | 127865 | 5333135 | 0 | 0 | NA | NA | left | 0 | 127865 | simple |
| ptg000007l | 177907 | large_qgap | 5461000 | 154954 | 218336 | NC_038251.2 | 19636409 | NC_038251.2 | 19640715 | 177907 | 5283093 | 0 | 0 | NA | NA | left | 0 | 177907 | simple |
| ptg000007l | 242879 | large_qgap | 5461000 | 234763 | 276006 | NC_038251.2 | 19633745 | NC_038251.2 | 19638869 | 242879 | 5218121 | 0 | 0 | NA | NA | left | 0 | 242879 | simple |
| ptg000007l | 324101 | large_qgap | 5461000 | 307828 | 353365 | NC_038251.2 | 19639711 | NC_038251.2 | 19637626 | 324101 | 5136899 | 0 | 0 | NA | NA | left | 0 | 324101 | simple |
| ptg000011l | 8222744 | unmapped_tail | 31842626 | 8222744 | 31842626 | NC_038255.2 | 43978809 | NC_038255.2 | 43978809 | 8222744 | 23619882 | 0 | 0 | NA | NA | left | 0 | 8222744 | simple |

---

## 📊 Step 7a. Visualize breakpoints with `04_breakwright_dotplot.py`

This complements `break_viz_plus.py` by giving a macroscopic alignment view.

#### 🧠 Overview
`break_dotplot.py` renders **dotplot-style figures** from **PAF** (contigs→reference):
- **Full-genome dotplot** — concatenates reference chromosomes on the X-axis.
- **Per-chromosome dotplots** — one figure per reference chromosome.
- **Break markers** — plots **red “×”** markers at candidate break coordinates projected onto the dotplot.

```bash
python 04_breakwright_dotplot.py \
    --paf contigs_vs_ref.paf \
    --breaks assembly_breaks_gfa.tsv \
    --outdir dotplots --mode full,per-chr --draw_chr_ticks --export_png
```

![Example dot plot figure output](images/dotplot_example.png)

#### 💻 Examples

Full genome + per-chromosome with breaks:
```bash
python 04_breakwright_dotplot.py   --paf contigs_vs_ref.paf   --breaks breaks/soy_breaks.tsv   --outdir dotplots   --mode full,per-chr   --paf_min_mapq 20 --paf_min_aln 5000 --paf_min_id 0.9   --draw_chr_ticks --export_pdf --export_png
```

#### 📤 Outputs

- **Full genome:** `<outdir>/dotplot_full.pdf/.png`  
- **Per-chromosome:** `<outdir>/dotplot_<tname>.pdf/.png`

Each figure shows PAF block segments: **x = reference coordinate**, **y = contig coordinate**.  
If `--breaks` is provided, per-break **red “×”** markers are drawn at projected positions.

#### ⚙️ Required Arguments

| Argument | Type | Description |
|---|---|---|
| `--paf` | *string (path)* | **Required.** PAF from `minimap2 -x asm5 ref.fa contigs.fa` (contigs→reference). |

#### 🧩 Optional Arguments (with Defaults)

| Argument | Default | Description |
|---|---|---|
| `--breaks` | *none* | Optional TSV with columns `qname,cut` (and optionally `reason`). Used to plot **red ×** markers. |
| `--outdir` | `dotplots` | Output directory for figures. |
| `--mode` | `full,per-chr` | Which plots to produce: `full`, `per-chr`, or `full,per-chr`. |
| `--win_for_break_proj` | `200000` | When projecting breaks, search ±win on contig to find the adjacent PAF blocks. |
| `--paf_min_mapq` | `20` | Minimum PAF MAPQ to include a block. |
| `--paf_min_aln` | `5000` | Minimum PAF alignment length to include (bp). |
| `--paf_min_id` | `0.90` | Minimum identity (nmatch/alen) to include. |
| `--max_blocks` | `1000000` | Max blocks to render for the **full-genome** plot (subsample beyond). |
| `--max_blocks_chr` | `200000` | Max blocks to render for **per-chromosome** plots. |
| `--ref_order` | *none* | Optional text file listing reference chromosome names (one per line) to define X-axis order. Defaults to natural-sorted names found in PAF. |
| `--draw_chr_ticks` | Flag | If set, draws vertical ticks at chromosome boundaries on the full-genome plot. |
| `--dpi` | `300` | DPI for PNG export. |
| `--export_pdf` | Flag | Export PDF figure(s). |
| `--export_png` | Flag | Export PNG figure(s). |
| `--label_reason` | Flag | Include `reason` in the legend/title if present in breaks. |
| `--limit_chroms` | *none* | Limit number of reference chromosomes (for quick tests). |

#### 🔬 Break Projection Details
For each break (`qname,cut`):
1. Find adjacent PAF blocks on that contig around `cut` (within `±win_for_break_proj`).  
2. If two blocks straddle the cut: mark **two red ×** at the **end of the left block** and **start of the right block** in dotplot space.  
3. If only one block is found: mark a single × at the **nearest end** of that block to the `cut`.  
4. On the **full-genome** plot, reference **x** is offset by cumulative chromosome lengths (concatenation).

#### 💻 Examples

Full genome + per-chromosome with breaks:
```bash
python 04_breakwright_dotplot.py   --paf contigs_vs_ref.paf   --breaks breaks/soy_breaks.tsv   --outdir dotplots   --mode full,per-chr   --paf_min_mapq 20 --paf_min_aln 5000 --paf_min_id 0.9   --draw_chr_ticks --export_pdf --export_png
```

Per-chromosome only, natural order:
```bash
python 04_breakwright_dotplot.py   --paf contigs_vs_ref.paf   --outdir dotplots   --mode per-chr --export_png
```

Custom reference order:
```bash
python 04_breakwright_dotplot.py   --paf contigs_vs_ref.paf   --breaks breaks/soy_breaks.tsv   --ref_order ref_order.txt   --mode full --export_pdf --export_png
```

#### 🧠 Notes
- Large assemblies can contain millions of blocks; performance guards (`--max_blocks*`) subsample beyond thresholds.  
- The PAF field `tlen` (column 7) is used to estimate chromosome lengths for offsetting.  
- If a break contig has no nearby blocks, it cannot be projected and is skipped with a warning. 

---

## 📊 Step 7b. Visualize breakpoints with `04_breakwright_dotplot.py`

Generates clear, publication-ready images showing alignment structure and supporting read coverage across candidate breaks.

#### 🧠 Overview
`break_viz_plus.py` produces **publication-ready panels** to visually confirm candidate misassembly breakpoints.
Given a **PAF** of contigs→reference, a **BAM** of HiFi reads mapped to contigs, and a TSV of curated or predicted **breaks**,
it renders per-break figures showing **local contig alignment structure** and **read support** around the cut site.

```bash
python 05_breakwright_viz_plus.py \
    --bam hifi_reads_to_contigs.bam \
    --breaks assembly_breaks_gfa.tsv \
    --outdir viz --export_png
```

#### 💻 Example

```bash
python break_viz_plus.py   --paf contigs_vs_ref.paf   --bam reads_vs_contigs.bam   --breaks breaks/soy_breaks.tsv   --outdir break_viz   --win 50000   --paf_min_mapq 20 --paf_min_aln 5000 --paf_min_id 0.9   --min_mapq 0 --export_pdf --export_png --label_coords --label_reason   --telomere_motif TTTAGGG --motif_min_run 5
```

![Example viz plus figure output](images/viz_plus_example.png)

#### 📤 Outputs

For each break (row in `--breaks`), the script writes one or two files (depending on flags):
- `<outdir>/<qname>_<cut>.pdf` (if `--export_pdf`)  
- `<outdir>/<qname>_<cut>.png` (if `--export_png`)

Each figure shows:
- **Top overlay:** PAF block spans across the window, relative to query coordinates (contig).  
- **Bottom plot:** Smoothed **read coverage** across the same window.  
- Vertical line at the **cut** position; optional labels, motif marks.

#### ⚙️ Required Arguments

| Argument | Type | Description |
|---|---|---|
| `--paf` | *string (path)* | **Required.** PAF file from `minimap2 -x asm5 ref.fa contigs.fa` (contigs→reference). |
| `--bam` | *string (path)* | **Required.** BAM of HiFi reads mapped to the contigs (e.g., `minimap2 -x map-hifi`). Must be indexed (`.bai`). |
| `--breaks` | *string (path)* | **Required.** TSV of break candidates (e.g., from `paf_breakfinder.py`), with columns: `qname, cut` (plus optional metadata). |
| `--outdir` | *string (path)* | **Required.** Output directory for figures. |

#### 🧩 Optional Arguments (with Defaults)

| Argument | Default | Description |
|---|---|---|
| `--win` | `50000` | Window size (bp) on each side of the cut to display (total span = 2×win). |
| `--min_mapq` | `0` | Minimum MAPQ for reads to count toward coverage. |
| `--max_reads` | `1000000` | Upper bound on reads to consider when computing coverage (for performance safety). |
| `--paf_min_mapq` | `20` | Minimum MAPQ to include a PAF alignment block in the overlay. |
| `--paf_min_aln` | `5000` | Minimum PAF block length to display. |
| `--paf_min_id` | `0.90` | Minimum identity for a PAF block to display. |
| `--dpi` | `300` | Figure DPI for PNG export. |
| `--font_size` | `10` | Base font size (axis labels, tick labels). |
| `--export_pdf` | Flag | If set, export a PDF per breakpoint. |
| `--export_png` | Flag | If set, export a PNG per breakpoint. |
| `--label_coords` | Flag | If set, plot coordinate ticks and label the cut site. |
| `--label_reason` | Flag | If set, annotate the figure title with the break reason(s) if present in the TSV. |
| `--telomere_motif` | `TTTAGGG` | Repeat motif to annotate (optional; soybean default). |
| `--motif_min_run` | `5` | Minimum number of consecutive motif repeats to annotate. |
| `--region_list` | *none* | Optional TSV listing `qname,cut,win` to override per-break window size. |
| `--limit` | *none* | Limit number of breaks to render (for testing). |

#### 🧠 Notes
- Ensure `--bam` is coordinate-sorted and indexed (`.bai`).  
- For large BAMs, consider downsampling or providing a `--region_list` to focus on selected breaks.  
- Figures are rendered with **matplotlib** (single axis per figure; default color cycle). 

---

## ✂️ Step 8. Split contigs with `06_breakwright_split_breaks.py`

Finally, manually curated breakpoints can be applied to generate a corrected assembly.

#### 🧠 Overview
`split_breaks.py` applies **curated breakpoints** to an assembly FASTA and writes a **corrected FASTA**.
It preserves original contig names by default and only appends suffixes (`a`, `b`, `c`, …) to contigs that are **actually split**.

- Accepts break coordinates from **PAF context** (default: `--cut_coord paf`) where `cut` is an index between bases
  (e.g., using `qend` from PAF, 0-based, end-exclusive), or **1-based** cuts (`--cut_coord one-based`).  
- Drops or keeps very small fragments according to a length threshold.

```bash
python 06_breakwright_split_breaks.py \
  --assembly_fa sample_assembly_filtered.fa \
  --breaks_tsv curated_breaks.tsv \
  --out_fa sample_assembly_filtered_corrected.fa
```

#### 💻 Examples

##### PAF-style cuts (default)
```bash
python 06_breakwright_split_breaks.py   --assembly_fa sample_assembly_filtered.fa   --breaks_tsv curated_breaks.tsv   --out_fa sample_assembly_filtered_corrected.fa   --cut_coord paf   --min_fragment_len 1000
```

##### 1-based cuts, drop tiny fragments
```bash
python 06_breakwright_split_breaks.py   --assembly_fa sample_assembly_filtered.fa   --breaks_tsv curated_breaks.tsv   --out_fa sample_assembly_filtered_corrected.fa   --cut_coord one-based   --min_fragment_len 5000   --min_gap_between_cuts 100
```

**Behavior:**
- Default: preserves original contig names  
- Contigs with breaks get suffixes (`a`, `b`, `c`, …)  
- Writes mapping tables (`.map.tsv`) and optionally drops very short fragments.

#### ⚙️ Required Arguments

| Argument | Type | Description |
|---|---|---|
| `--assembly_fa` | *string (path)* | **Required.** Input assembly FASTA (supports `.gz`). |
| `--breaks_tsv` | *string (path)* | **Required.** TSV of curated breakpoints with columns including `qname` and `cut`. Other columns are ignored. |
| `--out_fa` | *string (path)* | **Required.** Output FASTA for the split assembly. |

#### 🧩 Optional Arguments (with Defaults)

| Argument | Default | Description |
|---|---|---|
| `--cut_coord` | `paf` | Coordinate convention for `cut`. Options: `paf` (0-based index like PAF `qend`) or `one-based`. |
| `--min_fragment_len` | `1000` | Drop fragments shorter than this many bp. Set `0` to keep everything. |
| `--min_gap_between_cuts` | `50` | Ignore multiple cuts closer than this many bp (de-duplicate noisy cuts). |
| `--suffix_style` | `letters` | Suffix style for split fragments: `letters` → `a,b,c,...` (default), or `num` → `.1,.2,.3`. |
| `--wrap` | `60` | FASTA line wrap width. |
| `--report_prefix` | `<out_fa basename>` | Prefix for sidecar reports. If not given, derived from `--out_fa`. |
| `--gzip_out` | `False` | If set, write `--out_fa` as gzip (`.gz`). |

#### 🔍 Notes
- Breaks outside `[0, len]` or violating `--min_gap_between_cuts` are recorded in `*.skipped.tsv` and ignored.  
- Fragment coordinates reported in `.map.tsv` are **0-based, end-exclusive** `[start, end)` for clarity and easy interval math.  
- To **retain all fragments**, set `--min_fragment_len 0`.

---

## 🧠 Notes and Best Practices

- For large genomes, allocate at least **128 GB RAM** and **32 threads** for `minimap2` and `ragtag`.
- Always verify breakpoints visually before splitting contigs.
- Keep intermediate `*.bam` and `*.paf` files — they are invaluable for QC.
- Adjust thresholds depending on expected heterozygosity and repeat content.

---

## 📘 Citation

If you use this workflow in your research, please cite this repository and the underlying tools:

- Li, H. (2021). *Minimap2: pairwise alignment for nucleotide sequences.*  
- Kolmogorov et al. (2020). *Assembly of long, error-free reads using HiFi data.*
- BBMap, MUMmer, Samtools, and associated utilities.

---

## 📂 Directory Layout

```
project/
├── ref/
│   └── Gmax_Wm82_ref.fa
├── assembly/
│   └── assembly.fa
├── filtered/
│   ├── soy.kept.fa
│   ├── soy.decision.tsv
│   └── soy.kept_coverage.bed
├── breaks/
│   ├── soy_breaks.tsv
│   └── soy_breaks_summary.tsv
├── break_viz/
│   ├── contigX_break1.pdf
│   └── contigX_break1.png
└── final/
    └── soy.corrected.fa
```
