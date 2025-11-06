# BREAKWRIGHT: Precision-crafted genome correction.

This repository provides a reproducible workflow for detecting and correcting structural misassemblies in **genome assemblies** using a reference genome. It was developed using PacBio HiFi reads and assembled contigs from HiFiASM.

The pipeline filters low-value contigs, identifies potential misjoins (e.g., fused telomeres), visualizes breakpoints, and produces a corrected FASTA. It gives the user data to make informed curation decisions about which positions in assembled contigs to break.

---

## 🧬 Overview

This pipeline performs:

1. **Contig filtering** — remove short, redundant, or unmapped contigs based on reference alignment coverage.  
2. **Reference-guided comparison** — align contigs to the reference genome with `minimap2` or `nucmer`.  
3. **Breakpoint detection** — identify possible chimeric joins using PAF-based inspection (`paf_breakfinder.py`).  
4. **Breakpoint visualization** — plot local alignment structure and supporting read coverage (`break_viz_plus.py`).  
5. **Assembly correction** — split contigs at curated breakpoints (`split_breaks.py`).  

Optionally, long reads (HiFi or ONT) can be remapped to filtered contigs to validate structural breaks.

---

## 🧩 Dependencies

All tools are open-source and installable via `conda`:

- `python>=3.9`
- `numpy`, `pandas`, `matplotlib`
- `minimap2`
- `samtools`
- `mummer` (for dotplots)
- `bbtools` (for assembly stats)

---

## ⚙️ Create Conda Environment

```bash
conda create -p ./envs/binf python=3.10 -y
conda activate ./envs/binf

# Install dependencies
conda install -c bioconda minimap2 mummer bbmap samtools -y
pip install numpy pandas matplotlib
```

---

## 🧭 Step 1. Map contigs to the reference genome

Use `minimap2` to align assembled contigs to a high-quality reference genome (e.g., *Glycine max* Williams 82).

```bash
minimap2 -x asm5 -t 32 -a ref.fa assembly.fa > contigs_vs_ref.sam
samtools view -bS contigs_vs_ref.sam | samtools sort -o contigs_vs_ref.bam
samtools index contigs_vs_ref.bam

# Also produce a PAF for downstream filtering
minimap2 -x asm5 -t 32 ref.fa assembly.fa > contigs_vs_ref.paf
```

---

## 🔍 Step 2. Filter contigs with `contig_filter_by_paf.py`

Removes short contigs that are redundant or unmapped based on their novelty of reference coverage.

#### 🧠 Overview
This script filters contigs from an assembly FASTA based on their alignments to a reference genome (in PAF format).
It removes short or redundant contigs that either don't map or map to regions already covered by longer contigs.

```bash
python contig_filter_by_paf.py \
  --assembly_fa assembly.fa \
  --contigs_vs_ref_paf contigs_vs_ref.paf \
  --outprefix filtered/soy \
  --min_len 10000 \
  --min_mapq 20 \
  --min_aln_len 5000 \
  --min_identity 0.9 \
  --novel_bp_thresh 10000 \
  --novel_frac_thresh 0.25
```

#### 💻 Example Command

```bash
python contig_filter_by_paf.py   --assembly_fa assembly.fa   --contigs_vs_ref_paf contigs_vs_ref.paf   --outprefix filtered/soy   --min_len 10000   --min_mapq 20   --min_aln_len 5000   --min_identity 0.9   --novel_bp_thresh 10000   --novel_frac_thresh 0.25   --mode greedy
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

## 🧬 Step 3. Map HiFi reads to the filtered contigs

This step is used to confirm breakpoints and support manual curation.

#### 🧠 Overview
`paf_breakfinder.py` scans a **PAF** file of **contigs → reference** alignments to identify **candidate misassembly breakpoints**.
It looks for hallmarks such as **reference chromosome switches**, **strand flips**, **large intrachromosomal jumps**, and **large internal gaps on the contig** between adjacent alignment blocks.

```bash
minimap2 -x map-hifi -t 32 -a soy.kept.fa hifi_reads.fastq.gz > reads_vs_contigs.sam
samtools view -bS reads_vs_contigs.sam | samtools sort -o reads_vs_contigs.bam
samtools index reads_vs_contigs.bam
```
#### 💻 Example

```bash
python paf_breakfinder.py   --paf contigs_vs_ref.paf   --outprefix breaks/soy   --min_mapq 20   --min_aln_len 5000   --min_identity 0.9   --min_qgap 10000   --min_tjump 100000   --allow_overlap 1000   --max_merge_dist 10000
```

#### 📤 Outputs

| Output | Description |
|---|---|
| `<outprefix>_breaks.tsv` | Candidate breakpoints with columns: `qname, cut, reason, qlen, qbeg, qend, tname1, tpos1, tname2, tpos2`. |
| `<outprefix>_summary.tsv` | Per-contig summary of counts by reason (e.g., `switch_chr`, `strand_flip`, `large_tjump`, `large_qgap`). |
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
| `--min_qgap` | `10000` | Minimum **contig gap** (bp) between adjacent blocks to flag a candidate break. |
| `--min_tjump` | `100000` | Minimum **reference jump** (bp) across adjacent blocks on same chromosome to flag a candidate break. |
| `--allow_overlap` | `1000` | Allow up to this many bp of query overlap between adjacent blocks when evaluating gaps (tolerates small overlaps). |
| `--require_ordered` | `True` | If set (default), only consider adjacent blocks in ascending query order for breakpoint inference. |
| `--max_merge_dist` | `10000` | Merge nearby breakpoints on the same contig if within this many bp into a single candidate. |
| `--emit_all_blocks` | `False` | If set, also emits a TSV with all filtered alignment blocks for QC. |

#### 🔍 Detection Heuristics (adjacent blocks on the same contig)
- **Reference switch (`switch_chr`)**: consecutive blocks map to **different reference chromosomes**.  
- **Strand flip (`strand_flip`)**: strand changes between adjacent blocks on the same chromosome.  
- **Large reference jump (`large_tjump`)**: same chromosome but reference positions jump by `≥ --min_tjump`.  
- **Large contig gap (`large_qgap`)**: gap between query end of one block and query start of the next exceeds `≥ --min_qgap` (after allowing small overlap).

> Heuristics are conservative defaults; tune thresholds for your organism and assembly characteristics.

#### 🧠 Notes
- The script sorts blocks per contig by **query start** to infer adjacency in the assembled sequence.  
- Break candidates within `--max_merge_dist` are merged to reduce redundancy.  
- Outputs are intended to feed **Breakwright** steps like `break_viz_plus.py` and manual curation before `split_breaks.py`.

---

## 🔧 Step 4. Identify candidate breaks with `paf_breakfinder.py`

This script parses the PAF of contigs aligned to the reference to detect structural inconsistencies (e.g., multi-chromosome mappings, inversions, large internal gaps).

```bash
python paf_breakfinder.py \
  --paf contigs_vs_ref.paf \
  --outprefix breaks/soy
```

**Outputs:**
- `soy_breaks.tsv` — candidate break positions
- `soy_breaks_summary.tsv` — statistics and confidence scoring

---

## 📊 Step 5. Visualize breakpoints with `break_viz_plus.py`

Generates clear, publication-ready images showing alignment structure and supporting read coverage across candidate breaks.

```bash
python break_viz_plus.py \
  --paf contigs_vs_ref.paf \
  --bam reads_vs_contigs.bam \
  --breaks soy_breaks.tsv \
  --outdir break_viz
```

**Outputs:**
- `break_viz/*.pdf` and `.png` — visual confirmation plots

---

## ✂️ Step 6. Split contigs with `split_breaks.py`

Finally, manually curated breakpoints can be applied to generate a corrected assembly.

```bash
python split_breaks.py \
  --assembly_fa soy.kept.fa \
  --breaks_tsv curated_breaks.tsv \
  --out_fa soy.corrected.fa
```

**Behavior:**
- Default: preserves original contig names  
- Contigs with breaks get suffixes (`a`, `b`, `c`, …)  
- Writes mapping tables (`.map.tsv`) and optionally drops very short fragments.

---

## 🧾 Example Workflow Summary

| Step | Script | Input | Output | Purpose |
|------|---------|--------|---------|----------|
| 1 | minimap2 | contigs, reference | PAF/SAM/BAM | Reference alignment |
| 2 | contig_filter_by_paf.py | PAF, FASTA | kept.fa | Filter redundant contigs |
| 3 | minimap2 | HiFi reads, contigs | BAM | Read support mapping |
| 4 | paf_breakfinder.py | PAF | breaks.tsv | Detect structural anomalies |
| 5 | break_viz_plus.py | PAF, BAM | plots | Visual confirmation |
| 6 | split_breaks.py | FASTA, curated breaks | corrected.fa | Generate final assembly |

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
