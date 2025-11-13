#!/usr/bin/env python3
import sys, os, argparse, gzip
from collections import defaultdict

def fasta_iter(path):
    opener = gzip.open if path.endswith(".gz") else open
    name = None
    seq_chunks = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield (name, "".join(seq_chunks))
                name = line[1:].split()[0].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield (name, "".join(seq_chunks))

def load_fasta_lengths(path):
    d = {}
    for name, seq in fasta_iter(path):
        d[name] = len(seq)
    return d

def parse_paf(path, min_mapq, min_aln_len, min_identity):
    opener = gzip.open if path.endswith(".gz") else open
    blocks = defaultdict(list)  # qname -> list of dicts with tname, tmin, tmax, alen, mapq, ident
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
            ident = (nmatch / alen) if alen > 0 else 0.0
            if mapq < min_mapq or alen < min_aln_len or ident < min_identity:
                continue
            tmin = min(tstart, tend); tmax = max(tstart, tend)
            blocks[qname].append({"tname": tname, "tmin": tmin, "tmax": tmax,
                                  "alen": alen, "mapq": mapq, "ident": ident})
    return blocks

def merge_intervals(intervals, min_gap=0):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0][:]]
    for s,e in intervals[1:]:
        if s <= merged[-1][1] + min_gap:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s,e])
    return merged

def intervals_novel_len(new_intervals, covered_by_chrom):
    # Return bp of novel coverage that is not already in covered_by_chrom (list of merged intervals).
    novel = 0
    i = 0
    for (s,e) in new_intervals:
        pos = s
        while pos <= e:
            while i < len(covered_by_chrom) and covered_by_chrom[i][1] < pos:
                i += 1
            if i >= len(covered_by_chrom) or covered_by_chrom[i][0] > e:
                novel += (e - pos + 1)
                break
            cov_s, cov_e = covered_by_chrom[i]
            if cov_s > pos:
                novel += (min(e, cov_s - 1) - pos + 1)
                pos = cov_s
            pos = min(e+1, cov_e + 1)
    return novel

def add_coverage(covered_by_chrom, new_intervals):
    # Union new_intervals into covered_by_chrom in place.
    all_ints = covered_by_chrom + [i[:] for i in new_intervals]
    all_ints.sort(key=lambda x: (x[0], x[1]))
    merged = []
    for s,e in all_ints:
        if not merged or s > merged[-1][1]:
            merged.append([s,e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    covered_by_chrom[:] = merged

def contig_alignment_summary(blocks):
    # Summarize per contig: total aligned bp and intervals grouped by reference chrom.
    summary = {}
    for q, blist in blocks.items():
        alen_sum = sum(b["alen"] for b in blist)
        by_chrom = defaultdict(list)
        for b in blist:
            by_chrom[b["tname"]].append([b["tmin"], b["tmax"]])
        for t in by_chrom:
            by_chrom[t] = merge_intervals(by_chrom[t])
        summary[q] = {"aligned_bp": alen_sum, "by_chrom": by_chrom}
    return summary

def write_bed(path, cov):
    with open(path, "w") as fh:
        for chrom, ivals in cov.items():
            for s,e in ivals:
                fh.write(f"{chrom}\t{s-1}\t{e}\n")  # BED 0-based, half-open

def build_arg_parser():
    import argparse
    epilog = (
        "-----------------------------\n"
        "REQUIRED ARGUMENTS\n"
        "  --assembly_fa <path>         Input contig FASTA (can be .gz). Used to read contigs and lengths.\n"
        "  --contigs_vs_ref_paf <path>  PAF from `minimap2 -x asm5 ref.fa contigs.fa` (contigs→reference).\n"
        "  --outprefix <prefix>         Prefix for outputs (writes <prefix>.kept.fa, .dropped.list, .decision.tsv, .kept_coverage.bed).\n"
        "\n"
        "OPTIONAL ARGUMENTS (with defaults)\n"
        "  --min_len 10000              Contigs shorter than this (bp) are 'short'. Short contigs must add novel coverage to be kept.\n"
        "  --min_mapq 20                Minimum MAPQ to accept an alignment block from the PAF.\n"
        "  --min_aln_len 5000           Minimum alignment length (bp) to accept a block from the PAF.\n"
        "  --min_identity 0.90          Minimum identity (nmatch/alignment length) to accept a block.\n"
        "  --novel_bp_thresh 10000      Keep a short contig if it contributes at least this many novel reference bp.\n"
        "  --novel_frac_thresh 0.25     OR keep a short contig if at least this fraction of its aligned bp are novel.\n"
        "  --mode greedy                Selection order: 'greedy' (by contig length) or 'score' (by aligned bp).\n"
        "\n"
        "OUTPUTS\n"
        "  <prefix>.kept.fa             FASTA of contigs retained after filtering.\n"
        "  <prefix>.dropped.list        One dropped contig ID per line.\n"
        "  <prefix>.decision.tsv        Per-contig summary: contig, length, status, novel_bp, total_aligned_bp, decision.\n"
        "  <prefix>.kept_coverage.bed   Merged reference coverage contributed by kept contigs (BED, 0-based, half-open).\n"
        "\n"
        "EXAMPLE\n"
        "  python contig_filter_by_paf.py \\\n"
        "    --assembly_fa assembly.fa \\\n"
        "    --contigs_vs_ref_paf contigs_vs_ref.paf \\\n"
        "    --outprefix filtered/soy \\\n"
        "    --min_len 10000 \\\n"
        "    --min_mapq 20 \\\n"
        "    --min_aln_len 5000 \\\n"
        "    --min_identity 0.9 \\\n"
        "    --novel_bp_thresh 10000 \\\n"
        "    --novel_frac_thresh 0.25 \\\n"
        "    --mode greedy\n"
    )
    ap = argparse.ArgumentParser(
        description="Filter contigs using contig→reference PAF by removing short unmapped/redundant contigs.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog
    )
    # Required
    ap.add_argument("--assembly_fa", required=True, help="Input contig FASTA (supports .gz).")
    ap.add_argument("--contigs_vs_ref_paf", required=True, help="PAF from minimap2 -x asm5 ref.fa contigs.fa (contigs→reference).")
    ap.add_argument("--outprefix", required=True, help="Prefix for outputs (<prefix>.kept.fa, .dropped.list, .decision.tsv, .kept_coverage.bed).")
    # Optional
    ap.add_argument("--min_len", type=int, default=10000, help="Short contig threshold (bp). Short contigs must add novelty to be kept (default: 10000).")
    ap.add_argument("--min_mapq", type=int, default=20, help="Minimum MAPQ to accept a PAF block (default: 20).")
    ap.add_argument("--min_aln_len", type=int, default=5000, help="Minimum alignment length (bp) to accept a PAF block (default: 5000).")
    ap.add_argument("--min_identity", type=float, default=0.90, help="Minimum identity (nmatch/aln length) to accept a PAF block (default: 0.90).")
    ap.add_argument("--novel_bp_thresh", type=int, default=10000, help="Novel reference bp threshold to keep a short contig (default: 10000).")
    ap.add_argument("--novel_frac_thresh", type=float, default=0.25, help="Novelty fraction threshold to keep a short contig (default: 0.25).")
    ap.add_argument("--mode", choices=["greedy","score"], default="greedy", help="Selection order: 'greedy' (by contig length) or 'score' (by aligned bp). Default: greedy.")
    return ap

def main():
    ap = build_arg_parser()
    args = ap.parse_args()

    lengths = load_fasta_lengths(args.assembly_fa)
    blocks = parse_paf(args.contigs_vs_ref_paf, args.min_mapq, args.min_aln_len, args.min_identity)
    summary = contig_alignment_summary(blocks)

    contigs = list(lengths.keys())
    if args.mode == "greedy":
        contigs.sort(key=lambda q: lengths[q], reverse=True)
    else:
        contigs.sort(key=lambda q: summary.get(q, {}).get("aligned_bp", 0), reverse=True)

    covered = defaultdict(list)  # chrom -> merged intervals

    kept = []
    dropped = []
    decisions = []

    for q in contigs:
        qlen = lengths[q]
        s = summary.get(q)
        if s is None or s["aligned_bp"] == 0:
            if qlen >= args.min_len:
                kept.append(q)
                decisions.append([q, qlen, "unmapped", 0, 0, "KEPT_large_unmapped"])
            else:
                dropped.append(q)
                decisions.append([q, qlen, "unmapped", 0, 0, "DROPPED_short_unmapped"])
            continue

        total_aligned = 0
        total_novel = 0
        for chrom, ivals in s["by_chrom"].items():
            merged = ivals
            al_bp = sum(e-s+1 for s,e in merged)
            total_aligned += al_bp
            cov = covered[chrom]
            novel_bp = intervals_novel_len(merged, cov)
            total_novel += novel_bp

        novel_frac = (total_novel / total_aligned) if total_aligned > 0 else 0.0

        if qlen >= args.min_len:
            kept.append(q)
            for chrom, ivals in s["by_chrom"].items():
                add_coverage(covered[chrom], [iv[:] for iv in ivals])
            decisions.append([q, qlen, "mapped", total_novel, total_aligned, "KEPT_large"])
        else:
            if (total_novel >= args.novel_bp_thresh) or (novel_frac >= args.novel_frac_thresh):
                kept.append(q)
                for chrom, ivals in s["by_chrom"].items():
                    add_coverage(covered[chrom], [iv[:] for iv in ivals])
                decisions.append([q, qlen, "mapped", total_novel, total_aligned, "KEPT_small_novel"])
            else:
                dropped.append(q)
                decisions.append([q, qlen, "mapped", total_novel, total_aligned, "DROPPED_small_redundant"])

    kept_set = set(kept)
    with open(f"{args.outprefix}.kept.fa", "w") as outfa:
        for name, seq in fasta_iter(args.assembly_fa):
            if name in kept_set:
                outfa.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    outfa.write(seq[i:i+60] + "\n")

    with open(f"{args.outprefix}.dropped.list", "w") as fh:
        for q in dropped:
            fh.write(q + "\n")

    with open(f"{args.outprefix}.decision.tsv", "w") as fh:
        fh.write("contig\tlength\tstatus\tnovel_bp\ttotal_aligned_bp\tdecision\n")
        for row in decisions:
            fh.write("\t".join(map(str, row)) + "\n")

    write_bed(f"{args.outprefix}.kept_coverage.bed", covered)

    print(f"[OK] Kept contigs FASTA: {args.outprefix}.kept.fa")
    print(f"[OK] Dropped contigs list: {args.outprefix}.dropped.list")
    print(f"[OK] Decisions: {args.outprefix}.decision.tsv")
    print(f"[OK] Kept coverage BED: {args.outprefix}.kept_coverage.bed")

if __name__ == "__main__":
    main()
