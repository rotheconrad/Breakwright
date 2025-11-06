#!/usr/bin/env python3
import sys, os, argparse, gzip
from collections import defaultdict, namedtuple

PAF = namedtuple("PAF", "qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq")

def parse_paf(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                yield PAF(
                    qname=f[0], qlen=int(f[1]), qstart=int(f[2]), qend=int(f[3]), strand=f[4],
                    tname=f[5], tlen=int(f[6]), tstart=int(f[7]), tend=int(f[8]),
                    nmatch=int(f[9]), alen=int(f[10]), mapq=int(f[11])
                )
            except Exception:
                continue

def normalize_block(p):
    qmin = min(p.qstart, p.qend)
    qmax = max(p.qstart, p.qend)
    tmin = min(p.tstart, p.tend)
    tmax = max(p.tstart, p.tend)
    strand = "+" if p.qend >= p.qstart else "-"
    return {
        "qname": p.qname, "qlen": p.qlen, "qstart": p.qstart, "qend": p.qend,
        "qmin": qmin, "qmax": qmax, "strand": strand,
        "tname": p.tname, "tlen": p.tlen, "tstart": p.tstart, "tend": p.tend,
        "tmin": tmin, "tmax": tmax, "nmatch": p.nmatch, "alen": p.alen, "mapq": p.mapq
    }

def choose_primary_blocks(blocks, min_mapq=20, min_aln_len=10000):
    return [b for b in blocks if b["mapq"] >= min_mapq and b["alen"] >= min_aln_len]

def sort_blocks(blocks):
    return sorted(blocks, key=lambda b: (b["qname"], b["qmin"], b["qmax"]))

def detect_breaks(blocks, gap_thresh=50000, ref_jump_thresh=1000000):
    out = []
    for i in range(len(blocks)-1):
        a = blocks[i]; b = blocks[i+1]
        if a["qname"] != b["qname"]:
            continue
        qgap = b["qmin"] - a["qmax"]
        ref_jump = (abs(b["tmin"] - a["tmax"]) if a["tname"] == b["tname"] else None)
        chg_chr = (a["tname"] != b["tname"])
        chg_ori = (a["strand"] != b["strand"])
        big_gap = (qgap >= gap_thresh)
        big_ref = (ref_jump is not None and ref_jump >= ref_jump_thresh)
        if chg_chr or chg_ori or big_gap or big_ref:
            cut = (a["qmax"] + b["qmin"]) // 2
            reason = []
            if chg_chr: reason.append("CHR_SWITCH")
            if chg_ori: reason.append("ORI_SWITCH")
            if big_gap: reason.append(f"QGAP>={gap_thresh}")
            if big_ref: reason.append(f"REF_JUMP>={ref_jump_thresh}")
            out.append({
                "qname": a["qname"], "cut": cut, "qlen": a["qlen"], "reason": ";".join(reason),
                "prev_ref": a["tname"], "prev_strand": a["strand"], "prev_qmax": a["qmax"], "prev_tmax": a["tmax"],
                "next_ref": b["tname"], "next_strand": b["strand"], "next_qmin": b["qmin"], "next_tmin": b["tmin"],
                "prev_mapq": a["mapq"], "next_mapq": b["mapq"],
                "prev_alen": a["alen"], "next_alen": b["alen"],
            })
    return out

def index_reads_vs_contigs(reads_paf_path, min_mapq=0, min_aln_len=0):
    # Build an index: contig -> { read_name -> (tmin, tmax, alen, mapq) } choosing the best (max alen) per (read, contig).
    idx = defaultdict(dict)
    if not reads_paf_path:
        return idx
    for p in parse_paf(reads_paf_path):
        if p.mapq < min_mapq or p.alen < min_aln_len:
            continue
        tmin = min(p.tstart, p.tend); tmax = max(p.tstart, p.tend)
        d = idx[p.tname]
        cur = d.get(p.qname)
        if (cur is None) or (p.alen > cur[2]):
            d[p.qname] = (tmin, tmax, p.alen, p.mapq)
    return idx

def count_spanning_reads(read_index, contig, cut, flank=1000):
    # Count unique reads around a cut using the read index.
    reads = read_index.get(contig, {})
    span = set(); left = set(); right = set(); endnear = set()
    L0 = cut - flank; R0 = cut + flank
    for rname, (tmin, tmax, alen, mapq) in reads.items():
        if tmin <= L0 and tmax >= R0:
            span.add(rname); continue
        left_cov = (tmin <= L0 <= tmax)
        right_cov = (tmin <= R0 <= tmax)
        if left_cov and not right_cov:
            left.add(rname)
        elif right_cov and not left_cov:
            right.add(rname)
        if abs(tmin - cut) <= flank or abs(tmax - cut) <= flank:
            endnear.add(rname)
    return {
        "reads_span": len(span),
        "reads_left_only": len(left - span),
        "reads_right_only": len(right - span),
        "reads_end_near": len(endnear - span),
    }

def write_blocks_table(blocks_by_q, out_path):
    with open(out_path, "w") as fh:
        fh.write("\t".join(["qname","qlen","qstart","qend","qmin","qmax","strand","tname","tlen","tstart","tend","tmin","tmax","nmatch","alen","mapq"])+"\n")
        for q, blocks in blocks_by_q.items():
            for b in blocks:
                fh.write("\t".join(map(str,[b["qname"],b["qlen"],b["qstart"],b["qend"],b["qmin"],b["qmax"],b["strand"],b["tname"],b["tlen"],b["tstart"],b["tend"],b["tmin"],b["tmax"],b["nmatch"],b["alen"],b["mapq"]]))+"\n")

def write_breaks_table(candidates, out_path, read_index=None, flank=1000):
    header = [
        "qname","cut","qlen","reason",
        "prev_ref","prev_strand","prev_qmax","prev_tmax","prev_mapq","prev_alen",
        "next_ref","next_strand","next_qmin","next_tmin","next_mapq","next_alen"
    ]
    if read_index is not None:
        header += ["reads_span","reads_left_only","reads_right_only","reads_end_near","flank_bp"]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for c in candidates:
            row = [
                c["qname"], str(c["cut"]), str(c["qlen"]), c["reason"],
                c["prev_ref"], c["prev_strand"], str(c["prev_qmax"]), str(c["prev_tmax"]), str(c["prev_mapq"]), str(c["prev_alen"]),
                c["next_ref"], c["next_strand"], str(c["next_qmin"]), str(c["next_tmin"]), str(c["next_mapq"]), str(c["next_alen"])
            ]
            if read_index is not None:
                counts = count_spanning_reads(read_index, c["qname"], c["cut"], flank=flank)
                row += [str(counts["reads_span"]), str(counts["reads_left_only"]), str(counts["reads_right_only"]), str(counts["reads_end_near"]), str(flank)]
            fh.write("\t".join(row) + "\n")

def main():
    ap = argparse.ArgumentParser(description="Suggest breakpoints from contigs→reference PAF, optionally add HiFi read support across cuts.")
    ap.add_argument("--contigs_paf", required=True, help="PAF from minimap2 -x asm5 (contigs→reference). Supports .gz")
    ap.add_argument("--reads_vs_contigs_paf", default=None, help="Optional: PAF from minimap2 -x map-hifi (reads→contigs). Supports .gz")
    ap.add_argument("--flank", type=int, default=1000, help="Flank (bp) to define spanning coverage around cut (default: 1000 bp)")
    ap.add_argument("--min_mapq", type=int, default=20, help="Keep alignment blocks with MAPQ ≥ this")
    ap.add_argument("--min_aln_len", type=int, default=10000, help="Keep alignment blocks with alignment length ≥ this (bp)")
    ap.add_argument("--gap_thresh", type=int, default=50000, help="Suggest break if query gap ≥ this (bp)")
    ap.add_argument("--ref_jump_thresh", type=int, default=1000000, help="Suggest break if jump on same ref ≥ this (bp)")
    ap.add_argument("-o", "--outprefix", required=True, help="Output prefix, e.g., results/benning")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.outprefix) or ".", exist_ok=True)

    # Load contigs→reference blocks
    by_q = defaultdict(list)
    for rec in parse_paf(args.contigs_paf):
        by_q[rec.qname].append(normalize_block(rec))

    filt_sorted = {}
    for q, blks in by_q.items():
        kept = choose_primary_blocks(blks, args.min_mapq, args.min_aln_len)
        if not kept: continue
        filt_sorted[q] = sort_blocks(kept)

    # Detect candidate breaks
    candidates = []
    for q, blks in filt_sorted.items():
        candidates.extend(detect_breaks(blks, args.gap_thresh, args.ref_jump_thresh))

    # Optional: index reads→contigs best alignments
    read_index = None
    if args.reads_vs_contigs_paf:
        read_index = index_reads_vs_contigs(args.reads_vs_contigs_paf, min_mapq=10, min_aln_len=max(1000, args.flank*2))

    # Write outputs
    blocks_path = f"{args.outprefix}.blocks.tsv"
    breaks_path = f"{args.outprefix}.breaks.tsv"
    write_blocks_table(filt_sorted, blocks_path)
    write_breaks_table(candidates, breaks_path, read_index=read_index, flank=args.flank)

    print(f"[OK] Blocks: {blocks_path}")
    print(f"[OK] Candidate breaks: {breaks_path}")
    if read_index is None:
        print("[HINT] Provide --reads_vs_contigs_paf to annotate each break with spanning-read support.")
    else:
        print(f"[OK] Annotated breaks with read support using flank={args.flank} bp.")

if __name__ == "__main__":
    main()
