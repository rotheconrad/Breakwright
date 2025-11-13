#!/usr/bin/env python3
import sys, os, argparse, gzip
from collections import defaultdict, Counter

def parse_paf(path, min_mapq, min_aln_len, min_identity):
    opener = gzip.open if path.endswith(".gz") else open
    blocks = defaultdict(list)
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname = f[0]
            try:
                qlen = int(f[1]); qstart = int(f[2]); qend = int(f[3])
                strand = f[4]; tname = f[5]
                tlen = int(f[6]); tstart = int(f[7]); tend = int(f[8])
                nmatch = int(f[9]); alen = int(f[10]); mapq = int(f[11])
            except Exception:
                continue
            ident = (nmatch/alen) if alen > 0 else 0.0
            if mapq < min_mapq or alen < min_aln_len or ident < min_identity:
                continue
            tmin = min(tstart, tend); tmax = max(tstart, tend)
            qmin = min(qstart, qend); qmax = max(qstart, qend)
            blocks[qname].append({
                "qname": qname, "qlen": qlen, "qstart": qmin, "qend": qmax,
                "strand": strand, "tname": tname, "tstart": tmin, "tend": tmax,
                "nmatch": nmatch, "alen": alen, "mapq": mapq, "ident": ident
            })
    for q in blocks:
        blocks[q].sort(key=lambda b: (b["qstart"], b["qend"]))
    return blocks

def build_arg_parser():
    epilog = (
        "-----------------------------\n"
        "REQUIRED ARGUMENTS\n"
        "  --paf <path>                 PAF from `minimap2 -x asm5 ref.fa contigs.fa` (contigs→reference).\n"
        "  --outprefix <prefix>         Prefix for outputs (writes <prefix>_breaks.tsv and <prefix>_summary.tsv).\n"
        "\n"
        "OPTIONAL ARGUMENTS (with defaults)\n"
        "  --min_mapq 20                Minimum MAPQ to accept a PAF block.\n"
        "  --min_aln_len 5000           Minimum alignment length (bp) to accept a PAF block.\n"
        "  --min_identity 0.90          Minimum identity (nmatch/aln length) to accept a PAF block.\n"
        "  --min_qgap 10000             Min contig gap (bp) between adjacent blocks to flag 'large_qgap'.\n"
        "  --min_tjump 100000           Min reference jump (bp) on same chromosome to flag 'large_tjump'.\n"
        "  --allow_overlap 1000         Allow up to this many bp of query overlap between adjacent blocks.\n"
        "  --require_ordered            Only consider adjacency in ascending query order (default: set/True).\n"
        "  --max_merge_dist 10000       Merge candidate breaks within this distance (bp).\n"
        "  --emit_all_blocks            Also emit filtered alignment blocks TSV for QC.\n"
        "  --max_micro_overlap 5000     Flag 'micro_overlap' if adjacent blocks overlap by (allow_overlap, max_micro_overlap] bp.\n"
        "  --identity_drop 0.10         Flag 'identity_drop' if identity decreases by ≥ this fraction (e.g., 0.10 = 10%).\n"
        "  --low_mapq_edge 30           Flag 'low_mapq_edge' if either adjacent block has MAPQ ≤ this (still ≥ min_mapq).\n"
        "  --min_tail_unmapped 20000    Flag 'unmapped_lead'/'unmapped_tail' for large unmapped contig tails.\n"
        "\n"
        "OUTPUTS\n"
        "  <prefix>_breaks.tsv          Candidate breakpoints: qname,cut,reason,qlen,qbeg,qend,tname1,tpos1,tname2,tpos2\n"
        "  <prefix>_summary.tsv         Per-contig counts by reason.\n"
        "  <prefix>_blocks.tsv          (optional) All filtered PAF blocks.\n"
        "\n"
        "EXAMPLE\n"
        "  python paf_breakfinder.py \\\n"
        "    --paf contigs_vs_ref.paf \\\n"
        "    --outprefix breaks/soy \\\n"
        "    --min_mapq 20 --min_aln_len 5000 --min_identity 0.9 \\\n"
        "    --min_qgap 10000 --min_tjump 100000 --allow_overlap 1000 \\\n"
        "    --max_micro_overlap 5000 --identity_drop 0.10 --low_mapq_edge 30 \\\n"
        "    --min_tail_unmapped 20000 --max_merge_dist 10000\n"
    )
    ap = argparse.ArgumentParser(
        description="Detect candidate misassembly breakpoints from contigs→reference PAF alignments.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog
    )
    ap.add_argument("--paf", required=True, help="PAF file from minimap2 mapping contigs to reference.")
    ap.add_argument("--outprefix", required=True, help="Prefix for output files.")
    ap.add_argument("--min_mapq", type=int, default=20, help="Minimum MAPQ to keep a PAF block (default: 20).")
    ap.add_argument("--min_aln_len", type=int, default=5000, help="Minimum alignment length (bp) to keep a block (default: 5000).")
    ap.add_argument("--min_identity", type=float, default=0.90, help="Minimum identity (nmatch/aln length) (default: 0.90).")
    ap.add_argument("--min_qgap", type=int, default=10000, help="Min contig gap (bp) to flag a breakpoint (default: 10000).")
    ap.add_argument("--min_tjump", type=int, default=100000, help="Min reference jump (bp) on same chrom to flag breakpoint (default: 100000).")
    ap.add_argument("--allow_overlap", type=int, default=1000, help="Allow this many bp of query overlap (default: 1000).")
    ap.add_argument("--require_ordered", action="store_true", default=True, help="Consider blocks in ascending query order only (default: True).")
    ap.add_argument("--max_merge_dist", type=int, default=10000, help="Merge breakpoints within this distance (default: 10000).")
    ap.add_argument("--emit_all_blocks", action="store_true", help="Emit all filtered PAF blocks TSV.")
    ap.add_argument("--max_micro_overlap", type=int, default=5000, help="Flag 'micro_overlap' if overlap in (allow_overlap, max_micro_overlap] bp (default: 5000).")
    ap.add_argument("--identity_drop", type=float, default=0.10, help="Flag 'identity_drop' if drop ≥ this fraction (default: 0.10).")
    ap.add_argument("--low_mapq_edge", type=int, default=30, help="Flag 'low_mapq_edge' if either block has MAPQ ≤ this (default: 30).")
    ap.add_argument("--min_tail_unmapped", type=int, default=20000, help="Flag 'unmapped_lead/tail' if tail ≥ this many bp (default: 20000).")
    return ap

def main():
    ap = build_arg_parser()
    args = ap.parse_args()

    blocks = parse_paf(args.paf, args.min_mapq, args.min_aln_len, args.min_identity)

    # Optional QC dump
    if args.emit_all_blocks:
        with open(f"{args.outprefix}_blocks.tsv", "w") as fh:
            fh.write("qname\tqlen\tqstart\tqend\tstrand\ttname\ttstart\ttend\tnmatch\talen\tmapq\tident\n")
            for q, blist in blocks.items():
                for b in blist:
                    fh.write("\t".join(map(str, [
                        b["qname"], b["qlen"], b["qstart"], b["qend"], b["strand"],
                        b["tname"], b["tstart"], b["tend"], b["nmatch"], b["alen"], b["mapq"], f'{b["ident"]:.4f}'
                    ])) + "\\n")

    breaks = []
    summary = Counter()

    for qname, blist in blocks.items():
        if len(blist) == 0:
            # No mappings; consider whole contig as unmapped — no break suggested here.
            continue

        # Leading/trailing unmapped tails
        first = blist[0]; last = blist[-1]
        lead = first["qstart"]
        trail = first["qlen"] - last["qend"]
        if lead >= args.min_tail_unmapped:
            # cut just before first mapped base
            breaks.append({
                "qname": qname, "cut": lead - 1, "reason": "unmapped_lead",
                "qlen": first["qlen"], "qbeg": 0, "qend": first["qstart"],
                "tname1": first["tname"], "tpos1": (first["tstart"]+first["tend"])//2,
                "tname2": first["tname"], "tpos2": (first["tstart"]+first["tend"])//2
            })
            summary["unmapped_lead"] += 1
        if trail >= args.min_tail_unmapped:
            # cut at last mapped base
            breaks.append({
                "qname": qname, "cut": last["qend"], "reason": "unmapped_tail",
                "qlen": last["qlen"], "qbeg": last["qend"], "qend": last["qlen"],
                "tname1": last["tname"], "tpos1": (last["tstart"]+last["tend"])//2,
                "tname2": last["tname"], "tpos2": (last["tstart"]+last["tend"])//2
            })
            summary["unmapped_tail"] += 1

        # Adjacent block heuristics
        local = []
        for i in range(len(blist) - 1):
            a = blist[i]; b = blist[i+1]

            # Query gap / overlap
            qgap = b["qstart"] - a["qend"]
            overlap = max(0, a["qend"] - b["qstart"])  # positive if they overlap

            reasons = []

            # Core reasons
            if a["tname"] != b["tname"]:
                reasons.append("switch_chr")
            elif a["strand"] != b["strand"]:
                reasons.append("strand_flip")
            else:
                mid1 = (a["tstart"] + a["tend"]) // 2
                mid2 = (b["tstart"] + b["tend"]) // 2
                if abs(mid2 - mid1) >= args.min_tjump:
                    reasons.append("large_tjump")

            # large_qgap (after allowing small overlaps)
            if not (qgap < -args.allow_overlap) and (qgap >= args.min_qgap):
                reasons.append("large_qgap")

            # Additional reasons
            if overlap > args.allow_overlap and overlap <= args.max_micro_overlap:
                reasons.append("micro_overlap")
            if (a["ident"] - b["ident"]) >= args.identity_drop:
                reasons.append("identity_drop")
            if (a["mapq"] <= args.low_mapq_edge) or (b["mapq"] <= args.low_mapq_edge):
                reasons.append("low_mapq_edge")

            if reasons:
                cut = a["qend"]
                local.append({
                    "qname": qname,
                    "cut": cut,
                    "reason": ",".join(reasons),
                    "qlen": a["qlen"],
                    "qbeg": a["qstart"],
                    "qend": b["qend"],
                    "tname1": a["tname"], "tpos1": (a["tstart"]+a["tend"])//2,
                    "tname2": b["tname"], "tpos2": (b["tstart"]+b["tend"])//2
                })
                for r in reasons:
                    summary[r] += 1

        # Merge nearby cuts (by distance on query)
        if local:
            local.sort(key=lambda x: x["cut"])
            merged = [local[0]]
            for c in local[1:]:
                if c["cut"] - merged[-1]["cut"] <= args.max_merge_dist:
                    # merge reasons
                    mr = set(merged[-1]["reason"].split(","))
                    cr = set(c["reason"].split(","))
                    merged[-1]["reason"] = ",".join(sorted(mr | cr))
                else:
                    merged.append(c)
            breaks.extend(merged)

    # Write outputs
    if os.path.dirname(args.outprefix):
        os.makedirs(os.path.dirname(args.outprefix), exist_ok=True)

    with open(f"{args.outprefix}_breaks.tsv", "w") as fh:
        fh.write("qname\tcut\treason\tqlen\tqbeg\tqend\ttname1\ttpos1\ttname2\ttpos2\n")
        for c in breaks:
            fh.write("\t".join(map(str, [c["qname"], c["cut"], c["reason"], c["qlen"], c["qbeg"], c["qend"], c["tname1"], c["tpos1"], c["tname2"], c["tpos2"]])) + "\n")

    with open(f"{args.outprefix}_summary.tsv", "w") as fh:
        fh.write("reason\tcount\n")
        for k, v in Counter({r:0 for r in ['switch_chr','strand_flip','large_tjump','large_qgap','micro_overlap','identity_drop','low_mapq_edge','unmapped_lead','unmapped_tail']}).items():
            pass  # ensure header order (optional; counts come from summary below)
        for k, v in summary.most_common():
            fh.write(f"{k}\t{v}\n")

    print(f"[OK] Break candidates: {args.outprefix}_breaks.tsv")
    print(f"[OK] Summary:         {args.outprefix}_summary.tsv")
    if os.path.dirname(args.outprefix):
        print(f"[OK] (Outputs in directory) {os.path.dirname(args.outprefix)}")
    
if __name__ == "__main__":
    main()
