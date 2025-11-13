#!/usr/bin/env python3
import sys, os, argparse, gzip
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser(
        description="Split contigs at curated breakpoints with optional GFA guardrails."
    )
    ap.add_argument("--fasta", required=True, help="Input contig FASTA.")
    ap.add_argument("--breaks", required=True, help="Breaks TSV (qname,cut[,gfa_flag,nearest_junction_bp,...]).")
    ap.add_argument("--outfasta", required=True, help="Output FASTA with applied cuts.")
    ap.add_argument("--sep", default="\t", help="Breaks TSV delimiter (default: TAB).")
    ap.add_argument("--require_gfa_flag", help="Comma list; only apply cuts whose gfa_flag ∈ this set.")
    ap.add_argument("--max_nearest_junction_bp", type=int, help="Apply only if nearest_junction_bp ≤ value (if present).")
    ap.add_argument("--report_only", action="store_true", help="Do not write FASTA; only a report of intended changes.")
    ap.add_argument("--min_segment_len", type=int, default=0, help="Discard segments < this many bp after cutting (default: 0).")
    ap.add_argument("--suffix_style", choices=["letters","numeric"], default="letters", help="Suffix style for split parts (default: letters).")
    ap.add_argument("--preserve_names", action="store_true", help="Keep original names for contigs without cuts (default behavior).")
    ap.add_argument("--outreport", help="Report TSV path (default: <outfasta>.report.tsv).")
    return ap.parse_args()

def read_breaks(path, sep="\t"):
    header = None; rows = []
    with open(path, "rt") as fh:
        for line in fh:
            if not line.strip(): continue
            if header is None:
                header = line.rstrip("\n").split(sep)
                continue
            f = line.rstrip("\n").split(sep)
            rec = dict(zip(header, f))
            if "qname" not in rec or "cut" not in rec: 
                continue
            try:
                rec["cut"] = int(rec["cut"])
            except:
                continue
            if "nearest_junction_bp" in rec:
                try: rec["nearest_junction_bp"] = int(rec["nearest_junction_bp"]) if rec["nearest_junction_bp"] not in ("", "NA") else None
                except: pass
            rows.append(rec)
    return rows

def pass_gfa_filters(rec, require_flags=None, max_nj=None):
    if require_flags:
        flag = rec.get("gfa_flag","")
        if flag not in require_flags:
            return False
    if max_nj is not None:
        nj = rec.get("nearest_junction_bp", None)
        if nj is None or nj > max_nj:
            return False
    return True

def load_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    seqs = {}
    order = []
    name = None
    buf = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                    order.append(name)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
            order.append(name)
    return order, seqs

def write_fasta(path, records):
    with open(path, "wt") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

def suffix_generator(style="letters"):
    if style == "numeric":
        i = 1
        while True:
            yield f"_{i}"; i += 1
    else:
        import string
        alphabet = list(string.ascii_lowercase)
        n = 1
        while True:
            for comb in _letter_combinations(alphabet, n):
                yield comb
            n += 1

def _letter_combinations(alphabet, length):
    if length == 1:
        for ch in alphabet:
            yield ch
    else:
        for prefix in _letter_combinations(alphabet, length-1):
            for ch in alphabet:
                yield prefix + ch

def split_by_cuts(seq, cuts, min_len=0):
    L = len(seq)
    coords = [0] + [c for c in cuts if 0 < c < L] + [L]
    kept = []; discarded = []
    for i in range(len(coords)-1):
        s, e = coords[i], coords[i+1]
        if e - s >= min_len:
            kept.append((s, e))
        else:
            discarded.append((s, e))
    return kept, discarded

def main():
    args = parse_args()
    require_set = None
    if args.require_gfa_flag:
        require_set = {x.strip() for x in args.require_gfa_flag.split(",") if x.strip()}

    breaks = read_breaks(args.breaks, sep=args.sep)
    cuts_by_q = defaultdict(list)
    for r in breaks:
        q = r["qname"]; cut = r["cut"]
        if not pass_gfa_filters(r, require_set, args.max_nearest_junction_bp):
            continue
        cuts_by_q[q].append(cut)
    for q in list(cuts_by_q.keys()):
        cuts_by_q[q] = sorted(set(cuts_by_q[q]))

    order, seqs = load_fasta(args.fasta)

    out_records = []
    report_rows = []
    for q in order:
        seq = seqs[q]
        qlen = len(seq)
        cuts = [c for c in cuts_by_q.get(q, []) if 0 < c < qlen]
        if not cuts:
            out_name = q
            if not args.report_only:
                out_records.append((out_name, seq))
            report_rows.append({
                "qname": q, "qlen": qlen, "n_cuts": 0, "applied_cuts_csv": "",
                "kept_segments_csv": f"0-{qlen}",
                "discarded_segments_csv": "",
                "final_names_csv": out_name
            })
            continue

        kept, discarded = split_by_cuts(seq, cuts, min_len=args.min_segment_len)
        suffix = suffix_generator(args.suffix_style)
        final_names = []
        if not args.report_only:
            for (s, e) in kept:
                suf = next(suffix)
                name = f"{q}{suf}"
                out_records.append((name, seq[s:e]))
                final_names.append(name)
        else:
            for (s, e) in kept:
                suf = next(suffix)
                final_names.append(f"{q}{suf}")

        kept_csv = ",".join([f"{s}-{e}" for (s,e) in kept])
        disc_csv = ",".join([f"{s}-{e}" for (s,e) in discarded]) if discarded else ""
        report_rows.append({
            "qname": q, "qlen": qlen, "n_cuts": len(cuts),
            "applied_cuts_csv": ",".join(map(str, cuts)),
            "kept_segments_csv": kept_csv,
            "discarded_segments_csv": disc_csv,
            "final_names_csv": ",".join(final_names)
        })

    outreport = args.outreport or (args.outfasta + ".report.tsv")
    with open(outreport, "wt") as fh:
        header = ["qname","qlen","n_cuts","applied_cuts_csv","kept_segments_csv","discarded_segments_csv","final_names_csv"]
        fh.write("\t".join(header) + "\n")
        for r in report_rows:
            fh.write("\t".join(str(r[h]) for h in header) + "\n")

    if not args.report_only:
        write_fasta(args.outfasta, out_records)

    print(f"[OK] cuts_applied={sum(1 for v in cuts_by_q.values() if v)} kept_records={len(out_records)}")
    print(f"[OK] report: {outreport}")
    if not args.report_only:
        print(f"[OK] outfasta: {args.outfasta}")

if __name__ == "__main__":
    main()
