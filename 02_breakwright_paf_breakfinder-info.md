
## How `02_breakwright_paf_breakfinder.py` works

`02_breakwright_paf_breakfinder.py` scans high-quality contig→reference PAF alignments and looks for **abrupt changes between adjacent alignment blocks** along each contig. Each change-point is treated as a **candidate misassembly breakpoint** and annotated with one or more “reasons” (heuristics).

### 1. Filtering the PAF

We first read the PAF and keep only “good” alignment blocks:

* MAPQ ≥ `--min_mapq`
* Alignment length (column 11) ≥ `--min_aln_len`
* Identity `nmatch / aln_len` ≥ `--min_identity`

For each surviving PAF line we store:

* Query: `qname`, `qlen`, `qstart`, `qend`, `strand`
* Target: `tname`, `tstart`, `tend`
* Alignment stats: `nmatch`, `alen`, `mapq`, `ident`

For each contig (`qname`), blocks are **sorted by `qstart`** so we can walk left-to-right along the contig.

### 2. Leading / trailing unmapped tails

Before looking at internal structure, we ask whether the contig has large unmapped ends relative to the reference:

* `lead = first.qstart`
* `trail = qlen - last.qend`

If `lead ≥ --min_tail_unmapped`, we emit an **`unmapped_lead`** break at `cut = lead - 1`.

If `trail ≥ --min_tail_unmapped`, we emit an **`unmapped_tail`** break at `cut = last.qend`.

These are “edge” breakpoints where a long chunk of sequence does not align anywhere in the reference under the global filters.

### 3. Adjacent-block heuristics

For each contig we then walk **adjacent pairs** of alignment blocks `(a, b)`:

* `qgap = b.qstart - a.qend`
  (positive = gap on the contig; negative = overlap)
* `overlap = max(0, a.qend - b.qstart)`
* Reference midpoints:

  * `mid1 = (a.tstart + a.tend) / 2`
  * `mid2 = (b.tstart + b.tend) / 2`

We attach one or more **reasons** to the boundary between `a` and `b`.
If any reasons fire, we place a candidate break at:

* `cut = a.qend` (right edge of the left block)

Each candidate record stores the contig interval (`qbeg = a.qstart`, `qend = b.qend`) and the two reference midpoints (`tname1,tpos1` and `tname2,tpos2`).

### 4. Merging nearby cuts

Breaks called close together can reflect the same underlying issue. After scanning adjacent blocks, we:

* Sort all candidate cuts on a contig by `cut`
* Merge any whose `cut` positions are within `--max_merge_dist` bp
* Union their reason strings into a comma-separated list

The final output is a per-contig list of non-redundant break candidates with **one or more heuristic labels**.

---

## Heuristics: what is detected and why it matters

Below each bullet is both the **code logic** and **conceptual interpretation**.

---

### 1. `switch_chr` – block-to-block chromosome switches

**Detection (code):**

For adjacent blocks `a` and `b`:

```python
if a["tname"] != b["tname"]:
    reasons.append("switch_chr")
```

So this fires whenever the contig jumps from chromosome (or scaffold) `tname1` to a different `tname2` between consecutive high-quality blocks.

**Why it’s important for assemblies:**

* Strong signal of a **chimeric contig** that joins pieces from different chromosomes or unlinked scaffolds.
* Often arises from mis-joined long reads, graph tangles around repeats, or aggressive scaffolding with incorrect links.

**If it’s real biology instead of an error:**

* Could represent a true **inter-chromosomal translocation** or **transposition** if the contig spans a real breakpoint that the reference doesn’t fully match.
* In polyploid or highly rearranged genomes, might reflect **structural divergence** between sample and reference (especially across populations or subspecies).
* Also possible when the “reference” is itself a mosaic or fragmented; a single query contig can bridge pieces that are stored as separate chromosomes/scaffolds.

---

### 2. `strand_flip` – inversions between adjacent blocks

**Detection (code):**

If the target chromosome stays the same but orientation flips:

```python
elif a["strand"] != b["strand"]:
    reasons.append("strand_flip")
```

This labels a boundary where the contig goes from aligning `+` to `-` (or vice versa) on the **same** reference sequence.

**Why it’s important for assemblies:**

* Canonical signature of a **misoriented segment** or **inversion misassembly** inside a contig.
* Could also indicate a contig that was flipped around a repeat and rejoined incorrectly.

**If it’s real biology:**

* Simple explanation is a **bona fide inversion** in your sample relative to the reference, especially if:

  * Many reads support the arrangement.
  * Other contigs show the same pattern.
* Nested or tandem inversions, fold-back structures, or complex SVs can also show up as strand flips.

---

### 3. `large_tjump` – big jumps along the same chromosome

**Detection (code):**

When both chromosome and strand agree, we look at midpoints:

```python
mid1 = (a["tstart"] + a["tend"]) // 2
mid2 = (b["tstart"] + b["tend"]) // 2
if abs(mid2 - mid1) >= args.min_tjump:
    reasons.append("large_tjump")
```

So this flags a big spatial jump along a single chromosome (≥ `--min_tjump` bp) with no chromosome switch or strand flip.

**Why it’s important for assemblies:**

* Suggests a **long-range rearrangement** inside the contig:

  * Two far-apart regions on the reference welded together.
* Often points to mis-joins mediated by repeats, mis-threaded graph edges, or incorrect long-range scaffolding information.

**If it’s real biology:**

* Could be:

  * An **insertion** or **deletion** where a big chunk is present or absent in your sample.
  * A large **intra-chromosomal translocation**.
* In some crops / wild populations, long-distance rearrangements relative to the reference are common; multiple contigs with the same pattern strengthen the “real SV” interpretation.

---

### 4. `large_qgap` – big gaps on the contig between blocks

**Detection (code):**

We look at the gap on the query:

```python
qgap = b["qstart"] - a["qend"]

# Ignore big *overlaps* (qgap << 0) beyond allow_overlap
if not (qgap < -args.allow_overlap) and (qgap >= args.min_qgap):
    reasons.append("large_qgap")
```

So `large_qgap` fires when:

* The next block starts at least `min_qgap` bp after the previous one ends, **and**
* The blocks are not heavily overlapping (beyond the allowed overlap tolerance).

**Why it’s important for assemblies:**

* A large internal segment of the contig doesn’t map to the reference:

  * Could be an **unresolved gap**, **collapsed repeat region**, or **sequence absent from the reference**.
* If combined with other flags (e.g. `switch_chr` or `identity_drop`), it often marks a **junction between very different sequence contexts**.

**If it’s real biology:**

* May mark an **insertion** in your sample relative to the reference:

  * Novel transposon, introgressed segment, or hyper-divergent region that fails the mapping identity threshold.
* Could also reflect highly divergent **structural haplotypes** (e.g. large, divergent blocks in a pan-genome context).

---

### 5. `micro_overlap` – suspicious micro-overlaps between blocks

**Detection (code):**

We quantify query-space overlap:

```python
overlap = max(0, a["qend"] - b["qstart"])

if overlap > args.allow_overlap and overlap <= args.max_micro_overlap:
    reasons.append("micro_overlap")
```

So this fires when adjacent blocks **overlap** more than the “safe” tolerance, but not so much that they’re clearly redundant.

**Why it’s important for assemblies:**

* These micro-overlaps are a classic symptom of:

  * **Over-aggressive trimming / polishing** around joins.
  * Slightly duplicated segments when two pieces were joined with small redundant chunks.
* When seen consistently, can indicate a **systematic assembler/scaffolder artifact**.

**If it’s real biology:**

* Could represent **micro-duplications** or short tandem repeats that are slightly mis-aligned between sample and reference.
* In highly repetitive zones, small overlaps may simply reflect messy but real alignment patterns, especially if identity is high and MAPQ is moderate.

---

### 6. `identity_drop` – sudden drop in alignment identity

**Detection (code):**

We compare identity between successive blocks:

```python
if (a["ident"] - b["ident"]) >= args.identity_drop:
    reasons.append("identity_drop")
```

So if identity falls by ≥ `--identity_drop` (e.g. ≥ 0.10 = 10 percentage points) between block `a` and `b`, we flag the boundary.

**Why it’s important for assemblies:**

* A sharp identity drop often means the contig crosses a boundary between:

  * A **well-conserved region** and a **more divergent or mis-aligned region**.
* For misassemblies, this is a typical signature of:

  * **Collapsed paralogs** / gene family members.
  * Spurious extension of a contig into a **wrong copy** of a repeat or duplicated block.

**If it’s real biology:**

* Could be a genuine transition from:

  * **Conserved backbone** into a divergent haplotype block.
  * One **paralog** to another, where your assembly is “correct” but the reference is different.
* In pan-genome settings, these regions can correspond to **presence/absence variation** or a long tail of divergence that sits just below the global identity threshold.

---

### 7. `low_mapq_edge` – low-confidence blocks flanking a boundary

**Detection (code):**

```python
if (a["mapq"] <= args.low_mapq_edge) or (b["mapq"] <= args.low_mapq_edge):
    reasons.append("low_mapq_edge")
```

As long as blocks pass the global `--min_mapq` filter, we still treat **edges with “low but acceptable” MAPQ** as suspicious.

**Why it’s important for assemblies:**

* Low MAPQ typically means:

  * The alignment is **not uniquely placed**, usually due to repeats or high similarity to other loci.
* At a contig boundary, this suggests the contig is threading through **ambiguous terrain**, where mis-joins are more likely.

**If it’s real biology:**

* May simply reflect genuine **multi-copy regions** (segmental duplications, tandem arrays, transposon nests).
* In polyploids or hybrid genomes, mapping uncertainty is expected around **homeologous** or highly similar haplotypes.

---

### 8. `unmapped_lead` / `unmapped_tail` – long unmapped contig ends

**Detection (code):**

```python
lead = first["qstart"]
trail = first["qlen"] - last["qend"]

if lead >= args.min_tail_unmapped:
    # unmapped_lead at cut = lead - 1
if trail >= args.min_tail_unmapped:
    # unmapped_tail at cut = last["qend"]
```

So we flag contig ends that have no high-quality alignment for at least `--min_tail_unmapped` bp.

**Why it’s important for assemblies:**

* Large unmapped segments may be:

  * **Scaffolding artifacts** or leftover junk sequence.
  * Collapsed or highly erroneous regions that failed global filters.
* These tails are natural candidates for trimming in a misassembly-cleaning pipeline.

**If it’s real biology:**

* Could correspond to real **sequence not present in the reference**:

  * Novel telomeric/centromeric tracts, local introgressions, structural haplotypes, or population-specific regions.
* In pan-genome or highly divergent samples, these tails may be exactly the novel material you care about, so they’re flagged but shouldn’t be automatically discarded.

---

## Putting the heuristics together

A few additional points that help contextualize results:

* **Multiple reasons per break:**
  The script records a comma-separated list (e.g. `switch_chr,large_qgap,identity_drop`). Complex junctions often fire several heuristics; those are usually the most interesting.

* **Merging nearby cuts:**
  When multiple heuristics trigger within `--max_merge_dist`, they are merged into a **single breakpoint** with unioned reasons. This prevents over-fragmenting the contig at hyper-local messy regions.

* **Not all breaks are errors:**
  Every flagged cut is a **candidate** misassembly, but in highly rearranged or divergent genomes many will correspond to **real structural differences**. Downstream inspection (e.g. your new dot plots + read panels) is used to distinguish misassemblies from biological SV.
