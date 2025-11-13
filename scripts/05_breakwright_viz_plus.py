#!/usr/bin/env python3
import os, sys, argparse, random, math, re
from collections import defaultdict, Counter
import pysam
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    epilog = (
        "Example:\n"
        "  python break_viz_plus.py \n"
        "    --bam reads_to_contigs.bam \n"
        "    --breaks breaks/soy_breaks_gfa.tsv \n"
        "    --outdir viz --window 20000 --min_mapq 10 \n"
        "    --export_png --dpi 300 --annotate_gfa \n"
        "    --gfa_keep_flags junction_near_end,weak_overlap_end --gfa_max_nearest_junction_bp 20000\n"
    )
    ap = argparse.ArgumentParser(
        description="IGV-style panels around curated breaks with coverage and read segments (GFA-aware filters/annotations).",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog
    )
    ap.add_argument("--bam", required=True, help="HiFi read alignments to contigs (indexed .bam + .bai).")
    ap.add_argument("--breaks", required=True, help="Breaks TSV with at least qname,cut; optional gfa columns.")
    ap.add_argument("--outdir", required=True, help="Output directory for figures.")
    ap.add_argument("--window", type=int, default=20000, help="Half-window around cut (bp). Default: 20000.")
    ap.add_argument("--min_mapq", type=int, default=0, help="Minimum read MAPQ to include. Default: 0.")
    ap.add_argument("--max_reads", type=int, default=200000, help="Max reads per panel (subsample if exceeded). Default: 200000.")
    ap.add_argument("--dpi", type=int, default=300, help="DPI for PNG export. Default: 300.")
    ap.add_argument("--export_pdf", action="store_true", help="Export PDF figure(s).")
    ap.add_argument("--export_png", action="store_true", help="Export PNG figure(s).")
    ap.add_argument("--limit", type=int, help="Limit number of breaks to render.")
    ap.add_argument("--prefix", default="", help="Optional filename prefix.")
    ap.add_argument("--rng_seed", type=int, default=13, help="Random seed for subsampling reads. Default: 13.")
    # GFA-aware
    ap.add_argument("--gfa_keep_flags", help="Comma list; keep only rows with gfa_flag in this set.")
    ap.add_argument("--gfa_max_nearest_junction_bp", type=int, help="Keep only rows with nearest_junction_bp ≤ value (if present).")
    ap.add_argument("--annotate_gfa", action="store_true", help="Add a small GFA info box on the panel if columns present.")
    return ap.parse_args()

def read_breaks(path):
    header = None; rows = []
    with open(path, "rt") as fh:
        for line in fh:
            if not line.strip(): continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            f = line.rstrip("\n").split("\t")
            r = dict(zip(header, f))
            if "qname" in r and "cut" in r:
                try: r["cut"] = int(r["cut"])
                except: continue
                # parse optional numeric cols
                for k in ("nearest_junction_bp","deg_left","deg_right"):
                    if k in r:
                        try: r[k] = int(r[k])
                        except: pass
                rows.append(r)
    return rows, header or []

def filter_breaks(rows, keep_flags=None, max_nj=None):
    keep = None
    if keep_flags:
        keep = {x.strip() for x in keep_flags.split(",") if x.strip()}
    out = []
    for r in rows:
        if keep and r.get("gfa_flag","") not in keep: continue
        if max_nj is not None:
            nj = r.get("nearest_junction_bp", None)
            if nj is None or nj > max_nj: continue
        out.append(r)
    return out

def safe_fetch_limits(contig_len, start, end):
    s = max(0, start); e = min(contig_len, end) if contig_len is not None else max(start+1, end)
    if e <= s: e = s + 1
    return s, e

def get_contig_length(bam, qname):
    for sq in bam.header.get("SQ", []):
        if sq.get("SN") == qname:
            return sq.get("LN")
    return None

def compute_coverage_and_reads(bam, qname, start, end, min_mapq, max_reads, rng):
    length = end - start
    cov = np.zeros(length, dtype=np.int32)
    reads = []
    for col in bam.pileup(qname, start, end, truncate=True, min_base_quality=0, stepper="all"):
        pos = col.reference_pos
        if pos is None or pos < start or pos >= end: continue
        cov[pos - start] = col.nsegments
    it = bam.fetch(qname, start, end)
    for aln in it:
        if aln.is_unmapped or aln.mapping_quality < min_mapq: continue
        rstart = max(start, aln.reference_start if aln.reference_start is not None else start)
        rend = min(end, aln.reference_end if aln.reference_end is not None else start)
        if rend <= rstart: continue
        reads.append((rstart, rend))
    if len(reads) > max_reads:
        rng.shuffle(reads); reads = reads[:max_reads]
    return cov, reads

def plot_panel(qname, cut, start, end, cov, reads, annotate_text=None, dpi=300, out_png=None, out_pdf=None):
    width = 12; height = 6
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(width, height), gridspec_kw={"height_ratios":[1,1]}, sharex=True)
    xs = np.arange(start, end)
    ax1.plot(xs, cov, linewidth=1.0)
    ax1.set_ylabel("Coverage")
    ax1.axvline(cut, linestyle="--")
    ax1.set_title(f"{qname}:{cut}  [{start}-{end}]")
    for (rs, re) in reads:
        ax2.plot([rs, re], [0, 0], linewidth=2.0)
    ax2.set_ylim(-1, 1)
    ax2.set_yticks([])
    ax2.set_xlabel(f"Coordinate on {qname} (bp)")
    ax2.axvline(cut, linestyle="--")
    if annotate_text:
        ax2.text(cut, 0.6, annotate_text, fontsize=8, ha="left", va="bottom",
                 bbox=dict(boxstyle="round", alpha=0.2))
    fig.tight_layout()
    if out_pdf: fig.savefig(out_pdf, dpi=dpi, bbox_inches="tight")
    if out_png: fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.rng_seed)

    rows, header = read_breaks(args.breaks)
    rows = filter_breaks(rows, args.gfa_keep_flags, args.gfa_max_nearest_junction_bp)
    if args.limit is not None:
        rows = rows[:args.limit]

    bam = pysam.AlignmentFile(args.bam, "rb")
    contig_len_cache = {}

    for r in rows:
        qname = r["qname"]; cut = r["cut"]
        contig_len = contig_len_cache.get(qname)
        if contig_len is None:
            contig_len = get_contig_length(bam, qname)
            contig_len_cache[qname] = contig_len

        start, end = safe_fetch_limits(contig_len, cut - args.window, cut + args.window)
        cov, reads = compute_coverage_and_reads(bam, qname, start, end, args.min_mapq, args.max_reads, rng)

        annotate_text = None
        if args.annotate_gfa:
            parts = []
            if "gfa_flag" in r and r["gfa_flag"]:
                parts.append(r["gfa_flag"])
            if "deg_left" in r or "deg_right" in r:
                parts.append(f"L:{r.get('deg_left','')} R:{r.get('deg_right','')}")
            if "nearest_junction_bp" in r and r["nearest_junction_bp"] not in (None, "NA", ""):
                parts.append(f"NJ:{r['nearest_junction_bp']}bp")
            if "reason" in r and r["reason"]:
                parts.append(f"reason:{r['reason']}")
            annotate_text = "\n".join(parts) if parts else None

        base = f"{args.prefix}{qname}_{cut}"
        out_pdf = os.path.join(args.outdir, base + ".pdf") if args.export_pdf else None
        out_png = os.path.join(args.outdir, base + ".png") if args.export_png else None
        plot_panel(qname, cut, start, end, cov, reads, annotate_text, dpi=args.dpi, out_png=out_png, out_pdf=out_pdf)

    bam.close()
    print(f"[OK] Wrote figures to {args.outdir}")

if __name__ == "__main__":
    main()
