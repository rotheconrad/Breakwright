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

def write_fa(path, entries, width=60):
    with open(path, "w") as fh:
        for name, seq in entries:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i+width] + "\n")

def parse_breaks_tsv(path):
    # Return dict contig -> sorted unique list of cut positions (1-based, between bases).
    cuts = defaultdict(set)
    header = None
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                if "qname" not in header or "cut" not in header:
                    raise SystemExit("ERROR: breaks TSV must have columns 'qname' and 'cut'")
                continue
            rec = dict(zip(header, f))
            try:
                q = rec["qname"]
                c = int(rec["cut"])
            except Exception:
                continue
            cuts[q].add(c)
    return {k: sorted(v) for k,v in cuts.items()}

def suffix_letters(i):
    # Zero-based index to letters: 0->a, 1->b, ... 25->z, 26->aa, ...
    letters = []
    while True:
        i, rem = divmod(i, 26)
        letters.append(chr(ord('a') + rem))
        if i == 0:
            break
        i -= 1  # Excel-style base-26
    return "".join(reversed(letters))

def zero_pad(n, width):
    s = str(n)
    if len(s) < width:
        s = "0"*(width-len(s)) + s
    return s

def main():
    ap = argparse.ArgumentParser(description="Split assembly contigs at curated breakpoints.")
    ap.add_argument("--assembly_fa", required=True, help="Input assembly FASTA (supports .gz)")
    ap.add_argument("--breaks_tsv", required=True, help="Curated breaks TSV with columns at least: qname, cut (1-based)")
    ap.add_argument("--out_fa", required=True, help="Output FASTA with split contigs (plain FASTA)")
    ap.add_argument("--map_tsv", default=None, help="Write mapping TSV (old->new coordinates). Default: out_fa+'.map.tsv'")
    ap.add_argument("--dropped_tsv", default=None, help="Write TSV of dropped segments (if --min_seg_len>0). Default: out_fa+'.dropped.tsv'")
    # Default behavior: preserve original names; only split contigs get suffixes a,b,c...
    # If you want global renaming like contig0001a/b, enable --prefix_mode.
    ap.add_argument("--prefix_mode", action="store_true", help="Use global prefix numbering (contig0001a/b/...). Default off.")
    ap.add_argument("--prefix", default="contig", help="Base name prefix when --prefix_mode is on (default: contig)")
    ap.add_argument("--pad", type=int, default=4, help="Zero-pad width for prefix numbering (default: 4 -> contig0001)")
    ap.add_argument("--min_seg_len", type=int, default=0, help="If >0, drop segments shorter than this length")
    args = ap.parse_args()

    if args.map_tsv is None:
        args.map_tsv = args.out_fa + ".map.tsv"
    if args.dropped_tsv is None:
        args.dropped_tsv = args.out_fa + ".dropped.tsv"

    cuts_by_ctg = parse_breaks_tsv(args.breaks_tsv)

    out_entries = []
    map_rows = []
    dropped_rows = []
    global_index = 0  # used only in prefix_mode

    for old_name, seq in fasta_iter(args.assembly_fa):
        L = len(seq)
        cut_list = [c for c in cuts_by_ctg.get(old_name, []) if 1 <= c < L]
        cut_list = sorted(set(cut_list))

        if not cut_list:
            # No breaks for this contig
            if args.prefix_mode:
                global_index += 1
                base = f"{args.prefix}{zero_pad(global_index, args.pad)}"
                new_name = f"{base}a"
                out_entries.append((new_name, seq))
                map_rows.append("\t".join([old_name, new_name, "1", str(L), str(L), "1"]))
            else:
                # Preserve original name by default
                out_entries.append((old_name, seq))
                map_rows.append("\t".join([old_name, old_name, "1", str(L), str(L), "1"]))
            continue

        # There are cuts: split into segments
        bounds = []
        prev = 1
        for c in cut_list:
            bounds.append((prev, c))
            prev = c + 1
        bounds.append((prev, L))

        if args.prefix_mode:
            global_index += 1
            base = f"{args.prefix}{zero_pad(global_index, args.pad)}"
        # Build kept segments with per-contig suffixes a,b,c,... (only for kept ones to avoid gaps)
        kept_idx = 0
        for (s, e) in bounds:
            seg_len = e - s + 1
            if args.min_seg_len > 0 and seg_len < args.min_seg_len:
                # record drop; do not advance kept_idx (so suffixes are contiguous among kept)
                prop_name = (f"{base}{suffix_letters(kept_idx)}" if args.prefix_mode
                             else f"{old_name}{suffix_letters(kept_idx)}")
                dropped_rows.append("\t".join([old_name, str(s), str(e), str(seg_len), prop_name]))
                continue
            # kept segment
            if args.prefix_mode:
                new_name = f"{base}{suffix_letters(kept_idx)}"
            else:
                new_name = f"{old_name}{suffix_letters(kept_idx)}"
            out_entries.append((new_name, seq[s-1:e]))
            map_rows.append("\t".join([old_name, new_name, str(s), str(e), str(seg_len), str(kept_idx+1)]))
            kept_idx += 1

    # Write outputs
    write_fa(args.out_fa, out_entries)
    with open(args.map_tsv, "w") as fh:
        fh.write("old_name\tnew_name\told_start\told_end\tnew_len\tsegment_index\n")
        for row in map_rows:
            fh.write(row + "\n")

    if args.min_seg_len > 0:
        with open(args.dropped_tsv, "w") as fh:
            fh.write("old_name\told_start\told_end\tseg_len\tproposed_new_name\n")
            for row in dropped_rows:
                fh.write(row + "\n")

    print(f"[OK] Wrote FASTA: {args.out_fa}")
    print(f"[OK] Mapping TSV: {args.map_tsv}")
    if args.min_seg_len > 0:
        print(f"[OK] Dropped segments TSV: {args.dropped_tsv}  (min_seg_len={args.min_seg_len})")

if __name__ == "__main__":
    main()
