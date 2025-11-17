# Breakwright – HiFiASM GFA annotator

`03_breakwright_gfa_annotator.py`

## 1. High-level role in the Breakwright pipeline

`02_breakwright_paf_breakfinder.py` flags **candidate breakpoints** based purely on how contigs align to a reference (PAF). This is powerful, but it only sees a *linear* view of the assembly.

`03_breakwright_gfa_annotator.py` takes those candidate breaks and asks:

> *“What does the assembly **graph** think is happening at this contig and near this break?”*

Using the HiFiASM `p_ctg.gfa`, it:

1. Reads the **segment–link graph** for all primary contigs.
2. For each Breakwright breakpoint (`qname`, `cut`), it:

   * Measures **distance to contig ends**.
   * Looks at **degree / branching** at each end.
   * Summarizes **overlap strengths** at each end.
   * Classifies the break into simple flags:

     * `simple`
     * `tip_near_end`
     * `junction_near_end`
     * `weak_overlap_end`
3. Optionally exports a **local subgraph** around each contig for inspection in Bandage or other GFA viewers.

So the GFA annotator adds **graph-level evidence** that helps distinguish:

* “Clean, linear” regions vs.
* **Tips, junctions, and weakly supported joins** near the candidate break.

---

## 2. Inputs and outputs

### Inputs

1. **HiFiASM GFA graph** (`--gfa`)

   Typically the primary contig graph:
   `*.p_ctg.gfa` from HiFiASM.

2. **Break table from PAF breakfinder** (`--breaks`)

   TSV produced by `02_breakwright_paf_breakfinder.py`, must at least contain:

   * `qname` – contig ID (matches segment IDs in the GFA)
   * `cut` – candidate break position on the contig

3. **Key parameters**

   * `--end_proximity_bp` (default: 10,000 bp)
     How close a break must be to a contig end to be considered “near end.”
   * `--min_overlap_warn` (default: 50 bp)
     Overlaps shorter than this at contig ends are considered “weak.”
   * `--subgraph_hops` (default: 2)
     BFS radius in graph edges for subgraph export.
   * `--subgraph_min_overlap` (default: 0)
     Only traverse links with at least this much overlap.
   * `--emit_subgraph_gfa`
     Also emit mini-GFAs per break.

### Outputs

1. **Annotated breaks table**
   `<outprefix>_breaks_gfa.tsv`

   Your original break table plus new GFA-based columns (see below).

2. **Unmatched breaks (optional)**
   `<outprefix>_unmatched.tsv`
   Breaks whose `qname` wasn’t found in the GFA (e.g. contigs missing from the graph).

3. **Subgraph node/edge dumps**
   Directory: `<outprefix>_subgraphs/` (or `--subgraph_dir`)

   * One `.txt` summary per break: `qname_cut.txt`
   * Optional mini-GFA per break: `qname_cut.gfa` (with `--emit_subgraph_gfa`)
   * Index: `<outprefix>_subgraphs.tsv` listing all subgraphs and node files.

---

## 3. How the GFA is parsed

### 3.1 Segments and lengths

`parse_gfa()` reads all `S` records:

* Extracts **segment ID** (contig name).
* Gets length from:

  * Sequence length (if sequence present), or
  * `LN:i:<len>` tag if sequence is `*`, or
  * Falls back to 0 if neither is present.
* Stores in `seg_len[seg]`.

These lengths are used later as **contig lengths** (`qlen`), so `cut`, `dist_to_start`, and `dist_to_end` match the same coordinate system as your contigs.

### 3.2 Links, degrees, and overlaps

For each `L` line:

```text
L  a  o1  b  o2  ovl
```

* `a`, `b` – segment IDs.
* `o1`, `o2` – `+` or `-`, tell you which end of each contig is joined.
* `ovl` – CIGAR-like overlap string (e.g. `123M`).

The code:

1. Determines **which end** of each segment participates:

   * For segment `a`:

     * `o1 == '+'` → **right end** of `a`
     * `o1 == '-'` → **left end** of `a`
   * For segment `b`:

     * `o2 == '+'` → **left end** of `b`
     * `o2 == '-'` → **right end** of `b`

2. Increments **end-specific degrees**:

   ```python
   deg[a]['right'] += 1   # or ['left']
   deg[b]['left']  += 1   # or ['right']
   ```

   Intuitively:

   * `deg[seg]['left']` = number of edges attached to the left end.
   * `deg[seg]['right']` = number of edges attached to the right end.

3. Parses the overlap length (`\d+M`) from the CIGAR-like field and stores:

   ```python
   ovls[seg]['left/right'].append(overlap_bp)
   ```

4. Builds an **undirected adjacency list** (`adj`) for subgraph BFS, ignoring orientation:

   ```python
   adj[a].append((b, overlap_bp))
   adj[b].append((a, overlap_bp))
   ```

It also saves raw `S` and `L` lines so the script can write **mini-GFA** subgraphs later.

---

## 4. What gets added to each breakpoint

`annotate_breaks()` takes the original break rows and adds GFA context.

For each break record:

1. **Contig length and distances**

   ```python
   qlen          # from seg_len[qname]
   dist_to_start = d_left  = min(cut, qlen)
   dist_to_end   = d_right = qlen - cut
   ```

   These are simple distances from the break to the contig’s left and right ends.

2. **Degrees and overlap strengths at ends**

   ```python
   deg_left  = deg[q]['left']
   deg_right = deg[q]['right']
   min_ovl_left  = min(ovls[q]['left'])  or "NA"
   min_ovl_right = min(ovls[q]['right']) or "NA"
   ```

   * `deg_*` = **how many graph edges** attach to that end.
   * `min_ovl_*` = **shortest** overlap at that end (small = weak join).

3. **Nearest end and distance**

   The script decides which end is closer to the break:

   ```python
   if d_left <= d_right:
       nearest_end = 'left'
       nearest_end_degree = deg_left
       nearest_dist = d_left
   else:
       nearest_end = 'right'
       nearest_end_degree = deg_right
       nearest_dist = d_right
   ```

4. **Distance to nearest “junction-like” end: `nearest_junction_bp`**

   ```python
   candidates = []
   if deg_left  != 1: candidates.append(d_left)
   if deg_right != 1: candidates.append(d_right)

   nearest_junction_bp = min(candidates) if candidates else "NA"
   ```

   * **Degree = 1** is treated as a “simple linear link.”
   * Degree **0 or >1** are “non-linear”:

     * `deg = 0` → **tip** (dead end).
     * `deg > 1` → **junction** / branch point.
   * `nearest_junction_bp` is therefore:

     > *distance from the break to the closest end that is either a tip or a junction.*

5. **High-level flag: `gfa_flag`**

   Default: `"simple"`.

   Then:

   * If the break is **near an end** (within `--end_proximity_bp`, default 10 kb):

     ```python
     if nearest_dist <= end_prox_bp:
         if nearest_deg == 0:
             gfa_flag = "tip_near_end"
         elif nearest_deg > 1:
             gfa_flag = "junction_near_end"
     ```

     * `tip_near_end` – The break lies close to a dead-end path of the graph.
     * `junction_near_end` – The break lies close to a branch point.

   * Otherwise, if no tip/junction call and **weak overlaps** on an end:

     ```python
     if min_ovl_left  < min_ovl_warn or
        min_ovl_right < min_ovl_warn:
         if gfa_flag == "simple":
             gfa_flag = "weak_overlap_end"
     ```

     * `weak_overlap_end` – The contig end is supported only by short overlaps.

   If none of the above, the flag stays `"simple"` – graph looks clean and linear from the GFA’s perspective.

### 4.1 New columns (summary)

For each break you now get:

* `qlen` – contig length (from GFA).
* `dist_to_start`, `dist_to_end` – distance from break to contig ends.
* `deg_left`, `deg_right` – number of graph edges on each end.
* `min_ovl_left`, `min_ovl_right` – minimum overlap length at each end.
* `nearest_end` – `left` or `right` (which end is closer to break).
* `nearest_end_degree` – degree of that nearest end.
* `nearest_junction_bp` – distance to nearest tip/junction end (deg≠1). `"NA"` if both ends are degree 1.
* `gfa_flag` – one of `simple`, `tip_near_end`, `junction_near_end`, `weak_overlap_end`.

---

## 5. Subgraph export: putting breaks in local graph context

For each annotated break, `write_subgraphs()`:

1. Runs a **BFS on the adjacency graph** starting at `qname`:

   * Up to `--subgraph_hops` steps (default 2 edges away).
   * Only following edges with `overlap_bp ≥ --subgraph_min_overlap`.

2. Collects:

   * `nodes` – set of segments in the local neighborhood.
   * `edges` – segment pairs + overlap size.
   * `deg` – per-node degree **within the subgraph**.

3. Writes a **human-readable `.txt`**:

   ```text
   # Breakwright subgraph dump
   contig: ptg000XXX
   position: 123456
   nodes: N
   edges: E
   degree_min/median/mean/max: ...
   junction_count (deg>=3): ...
   junction_list:
     - ptg000YYY
   nodes_list:
     - ...
   edges_list:
     - ptg000XXX -- ptg000YYY overlap_bp=1234
   ```

4. Optionally writes a **mini-GFA** (`q_cut.gfa`) with:

   * All `S` lines for nodes in the subgraph.
   * All `L` lines connecting only those nodes.

5. Updates `<outprefix>_subgraphs.tsv` with:

   * `qname`, `cut`, `n_nodes`, `hops`, `min_overlap`,
   * Paths to `.txt` and `.gfa`,
   * A `nodes_csv` list of all segments.

This makes it easy to:

* Open the `.gfa` in **Bandage**, center on `qname`, and visually inspect:

  * Branches.
  * Tips.
  * Repeats.
  * Local graph complexity around each Breakwright breakpoint.

---

## 6. Biological interpretation: what do these GFA flags *mean*?

The GFA annotator doesn’t call anything “real biology” vs “misassembly” by itself.
Instead, it tells you **what kind of graph structure** is near each break, which you can interpret together with PAF evidence and other data.

### 6.1 `simple`

* **Graph view**:
  Both ends of the contig have degree 1, with reasonably sized overlaps; break is not near a tip or junction.
* **Possible interpretations**:

  * A **clean, linear contig**; break might reflect:

    * Genuine **structural differences** between your sample and the reference (e.g. inversion, translocation).
    * A breakpoint being **mis-placed** due to noise in PAF heuristics.
  * High confidence that the assembly graph itself doesn’t show a problem.

### 6.2 `junction_near_end`

* **Graph view**:
  The break lies within `end_proximity_bp` of an end whose degree > 1.
* **Biological / assembly scenarios**:

  * **Collapsed repeat** or **branching repeat** region:

    * Several contigs share a repeat or duplicated region, creating a junction.
  * **Haplotype divergence**:

    * Two paths representing different haplotypes could join/branch here.
  * **Potential misjoin**:

    * Assembler chose one path through a complex graph; alternate paths in subgraph might indicate misassembled junction.

**Interpretation**:
If the PAF break reason is `switch_chr`, `large_tjump`, or `strand_flip`, and the GFA flag is `junction_near_end`, this is a **prime suspect** for a misassembly at a graph branch.

### 6.3 `tip_near_end`

* **Graph view**:
  The break is near a contig end with degree 0 in the GFA graph → a **tip** (dead end).
* **Biological / assembly scenarios**:

  * **Low-coverage ends** that never connect confidently.
  * Regions flanking **repeats** where the assembler stops rather than wrongly joining.
  * **Chimeric read pile-ups** pruning off into a small branch.

**Interpretation**:

* If PAF suggests an `unmapped_tail`, a `tip_near_end` GFA flag supports the idea that the contig ended in a problematic or under-supported region.
* Biologically, these are often **edge effects**, low coverage, or unresolved repeats rather than clean structural variants.

### 6.4 `weak_overlap_end`

* **Graph view**:
  At least one contig end is supported only by **short overlaps** (`min_overlap_warn`, default < 50 bp), but there’s no strong signal of tip/junction near the break.
* **Biological / assembly scenarios**:

  * Contig ends that rely on **minimal shared sequence**:

    * Repeats with only short unique “k-mer anchors.”
  * **Borderline joins** where assembler had enough support to connect, but with weak overlap lengths.
* **Interpretation**:

  * These breaks might indicate potential **fragile join points** or regions that could break under slightly different assembly parameters.
  * Combine with PAF information:

    * If the PAF reason is `large_qgap`, `identity_drop`, or `low_mapq_edge`, and GFA flag is `weak_overlap_end`, you may be looking at a **suspicious low-confidence join**.

---

## 7. Putting it together with PAF heuristics

The real power comes from **combining** Breakwright’s PAF-based reasons with GFA flags:

* **Strong misassembly suspicion**:

  * PAF: `switch_chr`, `large_tjump`, `strand_flip`, `large_qgap`
  * GFA: `junction_near_end` or `weak_overlap_end`
  * Interpretation: contig crosses a graph branch / weak join *and* disagrees with the reference → likely misjoin.

* **Graph-supported structural variation candidate**:

  * PAF: `switch_chr` or `large_tjump`
  * GFA: `simple`
  * Interpretation: contig is linear and well supported in the assembly graph, but diverges from reference → candidate true SV / complex locus.

* **Unresolved ends / low-coverage tails**:

  * PAF: `unmapped_lead` / `unmapped_tail`
  * GFA: `tip_near_end` or `weak_overlap_end`
  * Interpretation: contig ends in a poorly supported graph tip → likely technical limitation rather than biological rearrangement.

* **Complex graph neighborhoods**:

  * Subgraph BFS reveals high junction counts, many nodes, short overlaps:

    * Suggests **repeats, CNVs, or multi-copy regions**.
    * Breakpoints in these contexts may be **intrinsically ambiguous**, not uniquely resolvable with current data.
