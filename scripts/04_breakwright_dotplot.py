#!/usr/bin/env python3
import os, sys, argparse, gzip, random, re
from collections import defaultdict, namedtuple, OrderedDict, Counter

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

Block = namedtuple("Block", "qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq ident")

def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]

def parse_paf(path, min_mapq=0, min_aln=0, min_id=0.0):
    """
    Parse a PAF file into per-query Block lists and a dict of reference "lengths".
    Here length is the maximum aligned coordinate per reference (R's chromMax),
    not the full reference length, to match the R script behavior.
    """
    opener = gzip.open if path.endswith(".gz") else open
    per_q = defaultdict(list)
    tlen_by_t = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 12:
                continue
            try:
                qname = f[0]; qlen = int(f[1]); qs = int(f[2]); qe = int(f[3])
                strand = f[4]; tname = f[5]
                tlen = int(f[6]); ts = int(f[7]); te = int(f[8])
                nmatch = int(f[9]); alen = int(f[10]); mapq = int(f[11])
            except Exception:
                continue
            ident = (nmatch / alen) if alen > 0 else 0.0
            if mapq < min_mapq or alen < min_aln or ident < min_id:
                continue

            per_q[qname].append(Block(
                qname, qlen, qs, qe, strand,
                tname, tlen, ts, te, nmatch, alen, mapq, ident
            ))

            # R-style chromMax: max aligned coordinate on this reference
            tmax = max(ts, te)
            tlen_by_t[tname] = max(tlen_by_t.get(tname, 0), tmax)

    # sort each query's blocks by query start coordinate for determinism
    for q in per_q:
        per_q[q].sort(key=lambda b: (min(b.qstart, b.qend), max(b.qstart, b.qend)))
    return per_q, tlen_by_t

def read_breaks(path):
    if not path:
        return [], []
    header = None
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            f = line.rstrip("\n").split("\t")
            rec = dict(zip(header, f))
            if "qname" in rec and "cut" in rec:
                try:
                    rec["cut"] = int(rec["cut"])
                except Exception:
                    continue
                flag = rec.get("gfa_flag", "")
                rec["gfa_flag"] = flag
                nj = rec.get("nearest_junction_bp", "")
                try:
                    rec["nearest_junction_bp"] = int(nj) if nj not in ("", "NA") else None
                except Exception:
                    rec["nearest_junction_bp"] = None
                rows.append(rec)
    return rows, header or []

def build_ref_offsets(tlen_by_t, ref_order=None):
    """
    Build concatenated reference offsets, preserving the CLI semantics of the
    script (either use user-supplied order or natural sort).
    """
    if ref_order:
        names = [x.strip() for x in ref_order if x.strip()]
    else:
        names = sorted(tlen_by_t.keys(), key=natural_key)
    offsets = OrderedDict()
    off = 0
    for t in names:
        L = tlen_by_t.get(t, 0)
        offsets[t] = (off, L)
        off += L
    return offsets

def filter_refs(per_q, tlen_by_t, keep_ref_n=0, ref_ids_csv=None):
    """
    Mimic R behavior:
      - If ref_ids_csv is provided, keep only those reference IDs (in that order).
      - Else:
          * Sort references by chromMax (tlen_by_t) descending.
          * If keep_ref_n > 0, keep only top N.
          * If keep_ref_n == 0, keep all, but keep the length-descending order.
    """
    # Explicit list of reference IDs overrides everything else
    if ref_ids_csv:
        allowed = [x.strip() for x in ref_ids_csv.split(",") if x.strip()]
        allowed_set = set(allowed)

        # Filter reference lengths
        tlen_by_t2 = {t: L for t, L in tlen_by_t.items() if t in allowed_set}

        # Filter per-query blocks
        new_per_q = {}
        for q, blist in per_q.items():
            new_list = [b for b in blist if b.tname in allowed_set]
            if new_list:
                new_per_q[q] = new_list

        # Preserve order according to 'allowed'
        tlen_by_t3 = OrderedDict()
        for t in allowed:
            if t in tlen_by_t2:
                tlen_by_t3[t] = tlen_by_t2[t]
        return new_per_q, tlen_by_t3

    # No explicit ref list: sort by chromMax (tlen_by_t) descending
    sorted_refs = sorted(tlen_by_t.items(), key=lambda kv: kv[1], reverse=True)

    if keep_ref_n and keep_ref_n > 0:
        sorted_refs = sorted_refs[:keep_ref_n]

    allowed_order = [t for t, _ in sorted_refs]
    allowed_set = set(allowed_order)

    # Filter per-query blocks
    new_per_q = {}
    for q, blist in per_q.items():
        new_list = [b for b in blist if b.tname in allowed_set]
        if new_list:
            new_per_q[q] = new_list

    # Ordered dict of reference "lengths"
    tlen_by_t2 = OrderedDict((t, tlen_by_t[t]) for t in allowed_order)
    return new_per_q, tlen_by_t2

def filter_queries_by_agg_len(per_q, min_query_aln):
    """
    R-like query filtering: drop queries whose total alignment length (sum of
    b.alen) is below min_query_aln. If min_query_aln <= 0, do nothing.
    """
    if not min_query_aln or min_query_aln <= 0:
        return per_q
    new_per_q = {}
    for q, blist in per_q.items():
        total = sum(b.alen for b in blist)
        if total >= min_query_aln:
            new_per_q[q] = blist
    return new_per_q

def project_breaks(per_q, breaks, win_for_proj):
    """
    Project breakpoints onto nearest left/right alignments on each contig,
    in raw PAF query/target coordinates. We later map them to the dotplot
    coordinate system so that layout and markers are consistent.
    """
    pts = []
    for rec in breaks:
        q = rec["qname"]; cut = rec["cut"]
        blist = per_q.get(q, [])
        if not blist:
            continue
        left = None; right = None
        for b in blist:
            qs = min(b.qstart, b.qend)
            qe = max(b.qstart, b.qend)
            if qe <= cut and (cut - qe) <= win_for_proj:
                if (left is None) or (qe > min(left.qstart, left.qend)):
                    left = b
            if qs >= cut and (qs - cut) <= win_for_proj:
                if (right is None) or (qs < min(right.qstart, right.qend)):
                    right = b
        if left is None and right is None:
            continue
        gfa_flag = rec.get("gfa_flag", "")
        nj = rec.get("nearest_junction_bp", None)
        reason = rec.get("reason", "")
        if left is not None:
            pts.append({
                "tname": left.tname,
                "x": max(left.tstart, left.tend),
                "y": max(left.qstart, left.qend),
                "qname": q,
                "cut": cut,
                "gfa_flag": gfa_flag,
                "nearest_junction_bp": nj,
                "reason": reason
            })
        if right is not None:
            pts.append({
                "tname": right.tname,
                "x": min(right.tstart, right.tend),
                "y": min(right.qstart, right.qend),
                "qname": q,
                "cut": cut,
                "gfa_flag": gfa_flag,
                "nearest_junction_bp": nj,
                "reason": reason
            })
    return pts

def color_size_for_flag(flag, scale):
    ms = 5 * scale; color = "gray"
    if flag == "junction_near_end":
        ms = 9 * scale; color = "red"
    elif flag == "tip_near_end":
        ms = 8 * scale; color = "orange"
    elif flag == "weak_overlap_end":
        ms = 8 * scale; color = "purple"
    return ms, color

def filter_breaks(breaks_rows, keep_flags=None, max_nj=None):
    if not breaks_rows:
        return breaks_rows
    out = []
    keep_set = None
    if keep_flags:
        keep_set = {x.strip() for x in keep_flags.split(",") if x.strip()}
    for rec in breaks_rows:
        if keep_set and rec.get("gfa_flag","") not in keep_set:
            continue
        if max_nj is not None:
            nj = rec.get("nearest_junction_bp", None)
            if nj is None or nj > max_nj:
                continue
        out.append(rec)
    return out

def build_dotplot_geometry(per_q, t_offsets):
    """
    Core port of the R dotPlotly layout logic into Python.

    Steps (for each query):
      1. Swap queryStart/queryEnd for negative-strand alignments.
      2. Decide if query is reverse-complement (sum(qEnd - qStart) < 0).
      3. If reverse-complement, flip coordinates: qMax - coord + 1.
      4. Re-base so minimum coordinate is 1.
      5. Stack queries along Y.

    Returns
    -------
    segs : list of dict
        Each dict has keys: tname, qname, x1, x2, y1, y2.
    qinfo : dict
        qname -> {rev, qmax, qmin, offset, ytick, height}
    mean_ident : dict
        qname -> mean identity (float).
    main_ref_by_q : dict
        qname -> tname of reference with max aggregate alignment length.
    """
    # Flatten all blocks
    all_blocks = []
    for q, blist in per_q.items():
        all_blocks.extend(blist)
    if not all_blocks:
        return [], {}, {}, {}

    # Mean identity per query
    mean_ident = {}
    for q, blist in per_q.items():
        if blist:
            mean_ident[q] = sum(b.ident for b in blist) / float(len(blist))
        else:
            mean_ident[q] = 0.0

    # Aggregate alignment length per ref per query
    agg_len_per_ref = defaultdict(lambda: defaultdict(int))  # q -> ref -> total len
    longest_blk_start = {}  # q -> (alen, x_ref_start)
    for q, blist in per_q.items():
        for b in blist:
            agg_len_per_ref[q][b.tname] += b.alen
            off, _ = t_offsets.get(b.tname, (0, 0))
            x_ref = off + b.tstart
            if q not in longest_blk_start or b.alen > longest_blk_start[q][0]:
                longest_blk_start[q] = (b.alen, x_ref)

    # Rank references in concatenated genome (by current t_offsets order)
    ref_rank = {t: i for i, t in enumerate(t_offsets.keys())}

    # Main reference (by aggregate alignment length) per query
    main_ref_by_q = {}
    for q, d in agg_len_per_ref.items():
        best_t = None
        best_len = -1
        for t, L in d.items():
            if L > best_len:
                best_len = L
                best_t = t
        main_ref_by_q[q] = best_t

    # Order queries roughly like the R script:
    #   - by reference rank of main target
    #   - then by position of longest alignment on that target
    order_keys = {}
    for q in per_q.keys():
        main_t = main_ref_by_q.get(q, None)
        r0 = ref_rank.get(main_t, len(ref_rank))
        x0 = longest_blk_start.get(q, (0, 0))[1]
        order_keys[q] = (r0, x0)
    q_order = sorted(per_q.keys(), key=lambda q: order_keys[q])

    # Build "fixed" coordinates (swap for minus strand)
    fixed_by_q = {}
    for q, blist in per_q.items():
        recs = []
        for b in blist:
            if b.strand == "-":
                qs_fixed = b.qend
                qe_fixed = b.qstart
            else:
                qs_fixed = b.qstart
                qe_fixed = b.qend
            recs.append({"block": b, "qs": qs_fixed, "qe": qe_fixed})
        fixed_by_q[q] = recs

    # Determine reverse-complement per query: sum(qEnd - qStart) < 0
    revcomp = {}
    for q, recs in fixed_by_q.items():
        diffs = [r["qe"] - r["qs"] for r in recs]
        revcomp[q] = (sum(diffs) < 0)

    segs = []
    qinfo = {}
    y_offset = 0.0

    for q in q_order:
        recs = fixed_by_q.get(q, [])
        if not recs:
            continue

        # queryMax from fixed coordinates
        coords = [c for r in recs for c in (r["qs"], r["qe"])]
        qmax_val = max(coords)

        # Step 3: orient coordinates (flip if reverse-complement)
        for r in recs:
            if revcomp[q]:
                qs2 = qmax_val - r["qs"] + 1
                qe2 = qmax_val - r["qe"] + 1
            else:
                qs2 = r["qs"]
                qe2 = r["qe"]
            r["qs2"] = qs2
            r["qe2"] = qe2

        # Step 4: re-base by per-query minimum
        coords2 = [c for r in recs for c in (r["qs2"], r["qe2"])]
        qmin_val = min(coords2)

        qmax_final = 0.0
        for r in recs:
            qs3 = r["qs2"] - qmin_val + 1
            qe3 = r["qe2"] - qmin_val + 1
            r["qs3"] = qs3
            r["qe3"] = qe3
            qmax_final = max(qmax_final, qs3, qe3)

        # Step 5: stack queries along Y with cumulative offset
        for r in recs:
            b = r["block"]
            off_x, _ = t_offsets.get(b.tname, (0, 0))
            x1 = off_x + b.tstart
            x2 = off_x + b.tend
            y1 = r["qs3"] + y_offset
            y2 = r["qe3"] + y_offset
            segs.append({
                "tname": b.tname,
                "qname": q,
                "x1": x1,
                "x2": x2,
                "y1": y1,
                "y2": y2
            })

        ytick = y_offset + qmax_final
        qinfo[q] = {
            "rev": revcomp[q],
            "qmax": qmax_val,   # max coord before re-basing
            "qmin": qmin_val,   # min coord after orientation, before re-basing
            "offset": y_offset,
            "ytick": ytick,
            "height": qmax_final
        }
        y_offset += qmax_final

    return segs, qinfo, mean_ident, main_ref_by_q

def map_breaks_to_dotplot(breaks_pts, t_offsets, qinfo):
    """
    Convert raw PAF-space breakpoint projections into the stacked query /
    concatenated reference coordinate system used for plotting, using the
    same R-style orientation and re-basing as for alignments.
    """
    mapped = defaultdict(list)  # flag -> list of mapped points
    if not breaks_pts:
        return mapped

    for p in breaks_pts:
        q = p.get("qname")
        t = p.get("tname")
        if q not in qinfo or t not in t_offsets:
            continue
        info = qinfo[q]
        off_x, _ = t_offsets[t]
        x_base = off_x + p["x"]

        # Apply the same orientation pipeline as for alignments
        q_pos = p["y"]
        if info["rev"]:
            q_oriented = info["qmax"] - q_pos + 1
        else:
            q_oriented = q_pos

        q_final = q_oriented - info["qmin"] + 1
        y_plot = q_final + info["offset"]

        rec = dict(p)
        rec["x_plot"] = x_base
        rec["y_plot"] = y_plot
        mapped[p.get("gfa_flag", "")].append(rec)

    return mapped

def compute_colors_for_segments(segs, mean_ident, main_ref_by_q,
                                use_similarity=True, identity_on_target=False):
    """
    Build the color array for segments:
      - If not use_similarity: return None (monochrome plotting).
      - If identity_on_target: only segments whose tname is the main ref for
        that query get color; others become NaN (plotted as 'bad' color).
      - Else: every segment gets the query's mean identity.
    """
    if not use_similarity:
        return None
    cols = []
    for s in segs:
        q = s["qname"]
        if identity_on_target:
            if main_ref_by_q.get(q) == s["tname"]:
                cols.append(mean_ident.get(q, 0.0))
            else:
                cols.append(np.nan)
        else:
            cols.append(mean_ident.get(q, 0.0))
    return np.array(cols, dtype=float)

def plot_full_genome(segs, qinfo, t_offsets, breaks_pts,
                     mean_ident, main_ref_by_q,
                     max_blocks, draw_chr_ticks, h_lines,
                     legend_by_gfa, marker_scale,
                     use_similarity, identity_on_target,
                     out_pdf, out_png, dpi, plot_size):
    """
    Full-genome dotplot using R-style layout and coloring. Overlays breakpoints
    as × markers in the same coordinate system.
    """
    if not segs:
        return

    # Possibly subsample segments
    if max_blocks and len(segs) > max_blocks:
        segs_plot = random.sample(segs, max_blocks)
    else:
        segs_plot = segs

    seg_xy = [((s["x1"], s["y1"]), (s["x2"], s["y2"])) for s in segs_plot]
    colors = compute_colors_for_segments(
        segs_plot, mean_ident, main_ref_by_q,
        use_similarity=use_similarity,
        identity_on_target=identity_on_target
    )

    # X/Y limits
    total_x = 0.0
    for t, (off, L) in t_offsets.items():
        total_x = max(total_x, off + L)
    ymax = max(info["ytick"] for info in qinfo.values())

    fig, ax = plt.subplots(figsize=(plot_size, plot_size))
    ax.set_title("Full-genome dotplot (target concatenated, queries stacked)")
    ax.set_xlabel("Reference coordinate (concatenated bp)")
    ax.set_ylabel("Query scaffolds / contigs")

    if colors is None:
        # similarity off: plain monochrome
        lc = LineCollection(seg_xy, colors="black", linewidth=0.25)
        ax.add_collection(lc)
        cbar = None
    else:
        cmap = plt.get_cmap("Spectral").copy()
        cmap.set_bad("lightgrey")
        lc = LineCollection(seg_xy, cmap=cmap, linewidth=0.25)
        colors_ma = np.ma.masked_invalid(colors)
        lc.set_array(colors_ma)
        ax.add_collection(lc)
        cbar = fig.colorbar(lc, ax=ax)
        cbar.set_label("Mean Percent Identity (per query)")

    # Y ticks: one per query
    qnames = list(qinfo.keys())
    yticks = [qinfo[q]["ytick"] for q in qnames]
    ylabels = [q[:20] for q in qnames]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=4, rotation=15)

    # Optional horizontal lines between queries
    if h_lines:
        for q in qnames:
            ax.axhline(qinfo[q]["ytick"], linewidth=0.1, color="grey", alpha=0.7)

    # X ticks & (optional) vertical boundaries
    t_names = list(t_offsets.keys())
    centers = []
    for t in t_names:
        off, L = t_offsets[t]
        centers.append(off + L / 2.0)
        if draw_chr_ticks:
            ax.axvline(off, linewidth=0.3, color="grey", linestyle="--")
    if draw_chr_ticks and t_names:
        last_off, last_L = t_offsets[t_names[-1]]
        ax.axvline(last_off + last_L, linewidth=0.3, color="grey", linestyle="--")

    ax.set_xticks(centers)
    ax.set_xticklabels(t_names, fontsize=6, rotation=90)

    ax.set_xlim(0, total_x)
    ax.set_ylim(0, ymax)
    ax.margins(x=0.005, y=0.01)

    # Overlay breakpoints
    mapped_breaks = map_breaks_to_dotplot(breaks_pts, t_offsets, qinfo)
    by_flag = mapped_breaks

    for flag, plist in by_flag.items():
        ms, color = color_size_for_flag(flag, marker_scale)
        xs = [p["x_plot"] for p in plist]
        ys = [p["y_plot"] for p in plist]
        ax.plot(xs, ys, marker='x', markersize=ms, linestyle='None',
                color=color, label=flag if flag else "no_flag")

    if legend_by_gfa and breaks_pts:
        counts = Counter(p.get("gfa_flag","") or "no_flag" for p in breaks_pts)
        handles = []; labels = []
        for flag, cnt in counts.items():
            ms, color = color_size_for_flag(flag if flag != "no_flag" else "", marker_scale)
            h, = ax.plot([], [], marker='x', linestyle='None', markersize=ms, color=color)
            handles.append(h); labels.append(f"{flag} ({cnt})")
        ax.legend(handles, labels, title="gfa_flag (count)", loc="upper right", fontsize=8)

    if out_pdf:
        plt.savefig(out_pdf, dpi=dpi, bbox_inches="tight")
    if out_png:
        plt.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def plot_per_chr(segs, qinfo, t_offsets, tname, breaks_pts,
                 mean_ident, main_ref_by_q,
                 max_blocks_chr, legend_by_gfa, marker_scale,
                 use_similarity, identity_on_target,
                 out_pdf, out_png, dpi, plot_size):
    """
    Per-chromosome dotplot: stacked query coordinate system but only shows
    alignments and breakpoints for a single reference sequence. Y-axis is
    restricted to queries that have alignments to this chromosome.
    """
    chr_segs = [s for s in segs if s["tname"] == tname]
    if not chr_segs:
        return

    if max_blocks_chr and len(chr_segs) > max_blocks_chr:
        chr_segs = random.sample(chr_segs, max_blocks_chr)

    seg_xy = [((s["x1"], s["y1"]), (s["x2"], s["y2"])) for s in chr_segs]
    colors = compute_colors_for_segments(
        chr_segs, mean_ident, main_ref_by_q,
        use_similarity=use_similarity,
        identity_on_target=identity_on_target
    )

    xmin = min(min(s["x1"], s["x2"]) for s in chr_segs)
    xmax = max(max(s["x1"], s["x2"]) for s in chr_segs)

    # Restrict Y ticks and limits to queries that actually map to this chromosome
    qnames_chr_set = {s["qname"] for s in chr_segs}
    qnames_chr = [q for q in qinfo.keys() if q in qnames_chr_set]

    if qnames_chr:
        yticks = [qinfo[q]["ytick"] for q in qnames_chr]
        ylabels = [q[:20] for q in qnames_chr]
        ymin = min(qinfo[q]["offset"] for q in qnames_chr)
        ymax = max(qinfo[q]["offset"] + qinfo[q]["height"] for q in qnames_chr)
    else:
        # Fallback: use all queries (shouldn't really happen if chr_segs is non-empty)
        yticks = [info["ytick"] for info in qinfo.values()]
        ylabels = [q[:20] for q in qinfo.keys()]
        ymin = 0.0
        ymax = max(info["ytick"] for info in qinfo.values())

    fig, ax = plt.subplots(figsize=(plot_size * 0.75, plot_size * 0.6))
    ax.set_title(f"Dotplot: {tname} (target = {tname}, queries stacked)")
    ax.set_xlabel(f"{tname} coordinate (bp)")
    ax.set_ylabel("Query scaffolds / contigs")

    if colors is None:
        lc = LineCollection(seg_xy, colors="black", linewidth=0.35)
        ax.add_collection(lc)
    else:
        cmap = plt.get_cmap("Spectral").copy()
        cmap.set_bad("lightgrey")
        lc = LineCollection(seg_xy, cmap=cmap, linewidth=0.35)
        colors_ma = np.ma.masked_invalid(colors)
        lc.set_array(colors_ma)
        ax.add_collection(lc)
        cbar = fig.colorbar(lc, ax=ax)
        cbar.set_label("Mean Percent Identity (per query)")

    # Y ticks: only queries that map to this chromosome
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=4, rotation=15)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.margins(x=0.005, y=0.01)

    # Overlay breakpoints (only those on this target)
    mapped_breaks = map_breaks_to_dotplot(breaks_pts, t_offsets, qinfo)
    by_flag = defaultdict(list)
    for flag, plist in mapped_breaks.items():
        by_flag[flag] = [p for p in plist if p["tname"] == tname]

    for flag, plist in by_flag.items():
        if not plist:
            continue
        ms, color = color_size_for_flag(flag, marker_scale)
        xs = [p["x_plot"] for p in plist]
        ys = [p["y_plot"] for p in plist]
        ax.plot(xs, ys, marker='x', markersize=ms, linestyle='None',
                color=color, label=flag if flag else "no_flag")

    if legend_by_gfa and any(by_flag.values()):
        all_pts = [p for plist in by_flag.values() for p in plist]
        counts = Counter(p.get("gfa_flag", "") or "no_flag" for p in all_pts)
        handles = []
        labels = []
        for flag, cnt in counts.items():
            ms, color = color_size_for_flag(flag if flag != "no_flag" else "", marker_scale)
            h, = ax.plot([], [], marker='x', linestyle='None', markersize=ms, color=color)
            handles.append(h)
            labels.append(f"{flag} ({cnt})")
        ax.legend(handles, labels, title="gfa_flag (count)", loc="upper right", fontsize=8)

    if out_pdf:
        plt.savefig(out_pdf, dpi=dpi, bbox_inches="tight")
    if out_png:
        plt.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def build_arg_parser():
    epilog = (
        "-----------------------------\n"
        "REQUIRED ARGUMENTS\n"
        "  --paf <path>                 PAF contigs→reference file.\n"
        "\n"
        "OPTIONAL ARGUMENTS (with defaults)\n"
        "  --breaks <path>              TSV 'qname,cut[,gfa_flag,nearest_junction_bp,reason]' for markers.\n"
        "  --outdir dotplots            Output directory (default: dotplots).\n"
        "  --mode full,per-chr          Which plots to produce: full, per-chr, or both (comma-separated).\n"
        "  --win_for_break_proj 200000  Window on contig for break projection (default: 200000).\n"
        "  --paf_min_mapq 20            Minimum PAF MAPQ to include a block.\n"
        "  --paf_min_aln 10000          Minimum PAF alignment length (bp).\n"
        "  --paf_min_id 0.90            Minimum PAF identity.\n"
        "  --min_query_aln 400000       Min aggregate alignment bp per query.\n"
        "  --max_blocks 1000000         Max blocks for full-genome plot (subsample over).\n"
        "  --max_blocks_chr 200000      Max blocks per chromosome plot (subsample over).\n"
        "  --ref_order <path>           Optional file listing reference chromosome order.\n"
        "  --keep_ref 0                 Number of reference chromosomes to keep (top by length; 0=all).\n"
        "  --ref_ids <csv>              Comma-separated list of reference IDs to keep (overrides --keep_ref).\n"
        "  --draw_chr_ticks             Draw vertical lines at chromosome boundaries on full-genome plot.\n"
        "  --h_lines                    Draw horizontal lines between queries on full-genome plot.\n"
        "  --plot_size 15               Plot size in inches (full plot is plot_size x plot_size).\n"
        "  --dpi 300                    DPI for PNG export.\n"
        "  --label_reason               (reserved; not used yet).\n"
        "  -v, --verbose                Print parameter settings (default: on).\n"
        "\n"
        "SIMILARITY COLORING (R-like)\n"
        "  --similarity                 Color alignments by mean identity per query (default: on).\n"
        "  --no_similarity              Disable identity coloring (monochrome lines).\n"
        "  --identity_on_target         Use identity only for on-target alignments; others grey.\n"
        "\n"
        "GFA-AWARE FILTERING & STYLING\n"
        "  --gfa_keep_flags <csv>       Keep only breaks whose gfa_flag is in this list.\n"
        "  --gfa_max_nearest_junction_bp N   Keep only breaks with nearest_junction_bp ≤ N.\n"
        "  --marker_scale 1.0           Scale marker sizes.\n"
        "  --legend_by_gfa              Add legend summarizing counts by gfa_flag.\n"
        "\n"
        "OUTPUTS (always written)\n"
        "  <outdir>/dotplot_full.pdf and .png\n"
        "  <outdir>/dotplot_<tname>.pdf and .png\n"
    )
    ap = argparse.ArgumentParser(
        description="Full-genome and per-chromosome dotplots from PAF; mark breaks with × (GFA-aware filtering/styling).",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog
    )
    ap.add_argument("--paf", required=True, help="PAF contigs→reference mapping (minimap2).")
    ap.add_argument("--breaks", help="TSV with columns qname,cut[,gfa_flag,nearest_junction_bp,reason] for × markers.")
    ap.add_argument("--outdir", default="dotplots", help="Output directory (default: dotplots).")
    ap.add_argument("--mode", default="full,per-chr", help="Which plots to produce: full, per-chr, or both (comma-separated).")
    ap.add_argument("--win_for_break_proj", type=int, default=200000, help="Window on contig for break projection (default: 200000).")
    ap.add_argument("--paf_min_mapq", type=int, default=20, help="Minimum PAF MAPQ (default: 20).")
    ap.add_argument("--paf_min_aln", type=int, default=10000, help="Minimum PAF alignment length (bp) (default: 10000).")
    ap.add_argument("--paf_min_id", type=float, default=0.90, help="Minimum PAF identity (default: 0.90).")
    ap.add_argument("--min_query_aln", type=int, default=400000,
                    help="Minimum aggregate alignment length per query (bp) (default: 400000).")
    ap.add_argument("--max_blocks", type=int, default=1000000, help="Max blocks for full-genome plot (default: 1e6).")
    ap.add_argument("--max_blocks_chr", type=int, default=200000, help="Max blocks per chromosome plot (default: 2e5).")
    ap.add_argument("--ref_order", help="Optional file listing reference chromosome names in desired order (one per line).")
    ap.add_argument("--keep_ref", type=int, default=0,
                    help="Number of reference chromosomes to keep (top by length). 0 = keep all.")
    ap.add_argument("--ref_ids", help="Comma-separated list of reference IDs to keep (overrides --keep_ref).")
    ap.add_argument("--draw_chr_ticks", action="store_true",
                    help="Draw vertical lines at chromosome boundaries on full plot.")
    ap.add_argument("--h_lines", action="store_true",
                    help="Draw horizontal lines between queries on full plot.")
    ap.add_argument("--plot_size", type=float, default=25.0,
                    help="Plot size (inches) for full-genome plot; per-chr plots scale from this.")
    ap.add_argument("--dpi", type=int, default=300, help="DPI for PNG export (default: 300).")
    ap.add_argument("--label_reason", action="store_true",
                    help="(Reserved) Include 'reason' in titles if present.")
    # Verbose
    ap.add_argument("-v", "--verbose", action="store_true", default=True,
                    help="Print out parameter settings (default: on).")
    # Similarity coloring
    ap.add_argument("--similarity", dest="similarity", action="store_true",
                    help="Color alignments by mean identity per query.")
    ap.add_argument("--no_similarity", dest="similarity", action="store_false",
                    help="Disable identity coloring (monochrome lines).")
    ap.set_defaults(similarity=True)  # similarity on by default
    ap.add_argument("--identity_on_target", action="store_true",
                    help="Use identity only for on-target alignments; others grey.")
    # GFA-aware filtering
    ap.add_argument("--gfa_keep_flags", help="Comma list of gfa_flag values to keep (others drop).")
    ap.add_argument("--gfa_max_nearest_junction_bp", type=int,
                    help="Keep only breaks with nearest_junction_bp ≤ this (if present).")
    ap.add_argument("--marker_scale", type=float, default=1.0,
                    help="Scale marker sizes (default: 1.0).")
    ap.add_argument("--legend_by_gfa", action="store_true",
                    help="Add legend summarizing counts by gfa_flag.")
    return ap

def main():
    ap = build_arg_parser()
    args = ap.parse_args()

    # Verbose parameter printout (R-style)
    if args.verbose:
        print("PARAMETERS:")
        print(f"  paf (--paf): {args.paf}")
        print(f"  breaks (--breaks): {args.breaks}")
        print(f"  outdir (--outdir): {args.outdir}")
        print(f"  mode (--mode): {args.mode}")
        print(f"  win_for_break_proj (--win_for_break_proj): {args.win_for_break_proj}")
        print(f"  paf_min_mapq (--paf_min_mapq): {args.paf_min_mapq}")
        print(f"  paf_min_aln (--paf_min_aln): {args.paf_min_aln}")
        print(f"  paf_min_id (--paf_min_id): {args.paf_min_id}")
        print(f"  min_query_aln (--min_query_aln): {args.min_query_aln}")
        print(f"  keep_ref (--keep_ref): {args.keep_ref}")
        print(f"  ref_ids (--ref_ids): {args.ref_ids}")
        print(f"  ref_order (--ref_order): {args.ref_order}")
        print(f"  plot_size (--plot_size): {args.plot_size}")
        print(f"  draw_chr_ticks (--draw_chr_ticks): {args.draw_chr_ticks}")
        print(f"  h_lines (--h_lines): {args.h_lines}")
        print(f"  similarity (--similarity): {args.similarity}")
        print(f"  identity_on_target (--identity_on_target): {args.identity_on_target}")
        print(f"  max_blocks (--max_blocks): {args.max_blocks}")
        print(f"  max_blocks_chr (--max_blocks_chr): {args.max_blocks_chr}")
        print(f"  gfa_keep_flags (--gfa_keep_flags): {args.gfa_keep_flags}")
        print(f"  gfa_max_nearest_junction_bp (--gfa_max_nearest_junction_bp): {args.gfa_max_nearest_junction_bp}")
        print(f"  marker_scale (--marker_scale): {args.marker_scale}")
        print(f"  legend_by_gfa (--legend_by_gfa): {args.legend_by_gfa}")
        print(f"  dpi (--dpi): {args.dpi}")
        print()

    os.makedirs(args.outdir, exist_ok=True)

    # Parse PAF and filter by per-block thresholds (unchanged)
    per_q, tlen_by_t = parse_paf(
        args.paf,
        min_mapq=args.paf_min_mapq,
        min_aln=args.paf_min_aln,
        min_id=args.paf_min_id
    )

    # R-like reference filtering
    per_q, tlen_by_t = filter_refs(per_q, tlen_by_t,
                                   keep_ref_n=args.keep_ref,
                                   ref_ids_csv=args.ref_ids)

    # R-like query aggregate length filter (min_query_aln)
    per_q = filter_queries_by_agg_len(per_q, args.min_query_aln)

    # Build reference offsets (after ref filtering)
    ref_order = None
    if args.ref_order:
        with open(args.ref_order) as fh:
            ref_order = [ln.strip() for ln in fh if ln.strip()]
    t_offsets = build_ref_offsets(tlen_by_t, ref_order)

    # Build R-style dotplot geometry + identity stats
    segs, qinfo, mean_ident, main_ref_by_q = build_dotplot_geometry(per_q, t_offsets)

    # Read and filter breakpoints (unchanged interface)
    breaks_rows, _ = read_breaks(args.breaks) if args.breaks else ([], [])
    if breaks_rows and (args.gfa_keep_flags or args.gfa_max_nearest_junction_bp is not None):
        breaks_rows = filter_breaks(breaks_rows, args.gfa_keep_flags, args.gfa_max_nearest_junction_bp)

    breaks_pts = project_breaks(per_q, breaks_rows, args.win_for_break_proj) if breaks_rows else []

    modes = [m.strip() for m in args.mode.split(",") if m.strip()]
    if "full" in modes:
        out_pdf = os.path.join(args.outdir, "dotplot_full.pdf")
        out_png = os.path.join(args.outdir, "dotplot_full.png")
        plot_full_genome(
            segs, qinfo, t_offsets, breaks_pts,
            mean_ident, main_ref_by_q,
            args.max_blocks,
            args.draw_chr_ticks, args.h_lines,
            args.legend_by_gfa, args.marker_scale,
            use_similarity=args.similarity,
            identity_on_target=args.identity_on_target,
            out_pdf=out_pdf, out_png=out_png,
            dpi=args.dpi, plot_size=args.plot_size
        )

    if "per-chr" in modes:
        for t in list(t_offsets.keys()):
            out_pdf = os.path.join(args.outdir, f"dotplot_{t}.pdf")
            out_png = os.path.join(args.outdir, f"dotplot_{t}.png")
            plot_per_chr(
                segs, qinfo, t_offsets, t, breaks_pts,
                mean_ident, main_ref_by_q,
                args.max_blocks_chr,
                args.legend_by_gfa, args.marker_scale,
                use_similarity=args.similarity,
                identity_on_target=args.identity_on_target,
                out_pdf=out_pdf, out_png=out_png,
                dpi=args.dpi, plot_size=args.plot_size
            )

    print(f"[OK] Dotplots written to: {args.outdir}")

if __name__ == "__main__":
    main()
