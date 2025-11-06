#!/usr/bin/env python3
import sys, os, argparse, gzip, json
from collections import defaultdict, namedtuple

# Optional dependency for depth from BAM
try:
    import pysam
except Exception:
    pysam = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAF = namedtuple("PAF", "qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq")

def parse_paf(path):
    if not path:
        return
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

def load_fasta_dict(path):
    d = {}
    if not path:
        return d
    for name, seq in fasta_iter(path):
        d[name] = seq
    return d

def clamp(a, lo, hi):
    return max(lo, min(hi, a))

def rc(s):
    comp = str.maketrans("ACGTacgtNn","TGCAtgcaNn")
    return s.translate(comp)[::-1]

def write_fa(path, entries, width=60):
    with open(path, "w") as fh:
        for name, seq in entries:
            fh.write(">"+name+"\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i+width] + "\n")

def find_adjacent_blocks(contig, cut, paf_blocks):
    blks = [b for b in paf_blocks if b.qname == contig]
    if not blks:
        return None, None
    blks = sorted(blks, key=lambda p: min(p.qstart, p.qend))
    prev_blk, next_blk = None, None
    for b in blks:
        qmin = min(b.qstart, b.qend)
        qmax = max(b.qstart, b.qend)
        if qmax <= cut:
            prev_blk = b
        elif qmin >= cut and next_blk is None:
            next_blk = b
    return prev_blk, next_blk

def reads_overlapping_window(reads_paf, contig, beg, end, max_reads=400):
    out = []
    for p in reads_paf:
        if p.tname != contig:
            continue
        tmin = min(p.tstart, p.tend); tmax = max(p.tstart, p.tend)
        if tmax < beg or tmin > end:
            continue
        out.append((p.qname, tmin, tmax, p.mapq))
    best = {}
    for r, tmin, tmax, mq in out:
        cur = best.get(r)
        span = tmax - tmin
        if (cur is None) or (span > cur[2]-cur[1]):
            best[r] = (r, tmin, tmax, mq)
    rows = list(best.values())
    rows.sort(key=lambda x: (x[1], x[2]-x[1]), reverse=False)
    if len(rows) > max_reads:
        rows = rows[:max_reads]
    return rows

def motif_counts_multi(seq, motifs):
    s = seq.upper()
    counts = {}
    for m in motifs:
        mm = m.upper()
        rr = rc(mm)
        c = 0; i = 0
        while True:
            i = s.find(mm, i)
            if i == -1: break
            c += 1; i += 1
        c_rc = 0; j = 0
        while True:
            j = s.find(rr, j)
            if j == -1: break
            c_rc += 1; j += 1
        counts[mm] = {"forward": c, "reverse": c_rc}
    return counts

def plot_reads(contig, cut, beg, end, rows, out_png):
    fig = plt.figure(figsize=(10, 6))
    ax = plt.gca()
    y = 0
    for r, tmin, tmax, mq in rows:
        ax.plot([tmin, tmax], [y, y])
        y += 1
    # Vertical line at cut
    ax.axvline(cut, linestyle="--")
    # Coordinate labels: beg, cut, end
    ax.annotate(f"{beg:,}", xy=(beg, y), xycoords=("data","data"),
                xytext=(0, 8), textcoords="offset points", ha="left", va="bottom", fontsize=9)
    ax.annotate(f"{cut:,}", xy=(cut, y), xycoords=("data","data"),
                xytext=(0, 8), textcoords="offset points", ha="center", va="bottom", fontsize=9)
    ax.annotate(f"{end:,}", xy=(end, y), xycoords=("data","data"),
                xytext=(0, 8), textcoords="offset points", ha="right", va="bottom", fontsize=9)
    ax.set_xlim(beg, end)
    ax.set_ylim(-1, max(5, y+2))
    ax.set_xlabel(contig+" position (bp)")
    ax.set_ylabel("reads (subset)")
    ax.set_title("Reads vs "+contig+" around cut "+str(cut))
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    # Also export PDF
    out_pdf = out_png[:-4] + ".pdf" if out_png.lower().endswith(".png") else out_png + ".pdf"
    fig.savefig(out_pdf)
    plt.close(fig)

def compute_depth_track(bam_path, contig, beg, end, bin_size=200):
    if not pysam:
        return None, "pysam_not_available"
    if not os.path.exists(bam_path):
        return None, "bam_missing"
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception:
        return None, "bam_open_failed"
    # Simple binned depth (avoid per-base to keep memory small)
    bins = []
    x = []
    for s in range(beg, end+1, bin_size):
        e = min(end, s + bin_size - 1)
        cov = 0
        try:
            for pileupcol in bam.pileup(contig, s-1, e, truncate=True):
                pos = pileupcol.reference_pos + 1
                if pos < s or pos > e:
                    continue
                cov += pileupcol.nsegments
        except Exception:
            bam.close()
            return None, "pileup_failed"
        bins.append(cov / max(1, (e - s + 1)))
        x.append((s + e) // 2)
    bam.close()
    return (x, bins), None

def plot_depth(contig, beg, end, x, y, out_png, cut=None):
    fig = plt.figure(figsize=(10, 3))
    ax = plt.gca()
    ax.plot(x, y)
    if cut is not None:
        ax.axvline(cut, linestyle="--")
        ax.annotate(f"{cut:,}", xy=(cut, max(y) if len(y)>0 else 0), xycoords=("data","data"),
                    xytext=(0, 8), textcoords="offset points", ha="center", va="bottom", fontsize=9)
    ax.annotate(f"{beg:,}", xy=(beg, max(y)/10 if len(y)>0 else 0), xycoords=("data","data"),
                xytext=(0, 8), textcoords="offset points", ha="left", va="bottom", fontsize=9)
    ax.annotate(f"{end:,}", xy=(end, max(y)/10 if len(y)>0 else 0), xycoords=("data","data"),
                xytext=(0, 8), textcoords="offset points", ha="right", va="bottom", fontsize=9)
    ax.set_xlim(beg, end)
    ax.set_xlabel(contig+" position (bp)")
    ax.set_ylabel("depth")
    ax.set_title("Read depth")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    out_pdf = out_png[:-4] + ".pdf" if out_png.lower().endswith(".png") else out_png + ".pdf"
    fig.savefig(out_pdf)
    plt.close(fig)

def write_multifasta(path, entries):
    # entries: list of (name, seq)
    write_fa(path, entries, width=60)

def make_report_html(path, info):
    parts = []
    parts.append("<!DOCTYPE html><html><head><meta charset='utf-8'><title>Break visual report</title></head><body>")
    parts.append("<h2>{} @ {}</h2>".format(info['contig'], info['cut']))
    parts.append("<ul>")
    parts.append("<li>Window: {} - {}</li>".format(info['beg'], info['end']))
    parts.append("<li>Reads (subset) plotted: {}</li>".format(info.get('reads_plotted', 0)))
    parts.append("<li>Motif counts: {}</li>".format(info.get('motifs_text',"NA")))
    if info.get('prev_ref'):
        parts.append("<li>Prev ref block: {}</li>".format(info['prev_ref']))
    if info.get('next_ref'):
        parts.append("<li>Next ref block: {}</li>".format(info['next_ref']))
    parts.append("</ul>")
    if os.path.exists(os.path.join(os.path.dirname(path), "readplot.png")):
        parts.append("<p><img src='readplot.png' alt='reads around cut'></p>")
    if os.path.exists(os.path.join(os.path.dirname(path), "depthplot.png")):
        parts.append("<p><img src='depthplot.png' alt='depth around cut'></p>")
    parts.append("<p>Files in this folder:</p><pre>")
    parts.append("\n".join(info['files']))
    parts.append("</pre></body></html>")
    with open(path, "w") as fh:
        fh.write("".join(parts))

def write_index_html(outdir, entries):
    # entries: list of dicts with keys contig, cut, relpath
    parts = []
    parts.append("<!DOCTYPE html><html><head><meta charset='utf-8'><title>Break index</title></head><body>")
    parts.append("<h2>Break reports</h2><ul>")
    for e in entries:
        parts.append("<li><a href='{}'>{} @ {}</a></li>".format(e["relpath"], e["contig"], e["cut"]))
    parts.append("</ul></body></html>")
    with open(os.path.join(outdir, "index.html"), "w") as fh:
        fh.write("".join(parts))

def main():
    ap = argparse.ArgumentParser(description="Per-break visual confirmations and exports for publication.")
    ap.add_argument("--assembly_fa", required=True)
    ap.add_argument("--breaks_tsv", required=True, help="Curated breaks TSV from paf_breakfinder.py (keep only desired rows).")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank", type=int, default=20000, help="bp around cut to export/plot (default: 20 kb)")
    ap.add_argument("--reads_vs_contigs_paf", default=None, help="PAF from minimap2 -x map-hifi contigs.fa reads.fastq.gz")
    ap.add_argument("--contigs_vs_ref_paf", default=None, help="PAF from minimap2 -x asm5 ref.fa assembly.fa")
    ap.add_argument("--ref_fa", default=None, help="Reference FASTA to export windows at mapped positions (optional)")
    ap.add_argument("--max_reads", type=int, default=300, help="Cap reads plotted to avoid huge images")
    ap.add_argument("--bam_reads_vs_contigs", default=None, help="Optional BAM (sorted + indexed) for depth track (requires pysam)")
    ap.add_argument("--bin_size", type=int, default=200, help="Bin size for depth (bp)")
    ap.add_argument("--motifs", default="TTTAGGG", help="Comma-separated motifs to count in window (reverse complements counted separately)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    asm = load_fasta_dict(args.assembly_fa)
    reads_paf = list(parse_paf(args.reads_vs_contigs_paf)) if args.reads_vs_contigs_paf else []
    contig_ref_paf = list(parse_paf(args.contigs_vs_ref_paf)) if args.contigs_vs_ref_paf else []
    ref = load_fasta_dict(args.ref_fa) if args.ref_fa else {}

    motif_list = [m.strip() for m in args.motifs.split(",") if m.strip()]

    # Read curated breaks TSV
    hdr = None
    breaks = []
    with open(args.breaks_tsv) as fh:
        for line in fh:
            if not line.strip(): continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            rec = dict(zip(hdr, f))
            rec["cut"] = int(rec["cut"])
            rec["qlen"] = int(rec.get("qlen","0") or 0)
            breaks.append(rec)

    index_entries = []

    for rec in breaks:
        q = rec["qname"]; cut = rec["cut"]
        if q not in asm:
            print("[WARN] {} not in assembly; skipping".format(q))
            continue
        seq = asm[q]
        L = len(seq)
        beg = clamp(cut - args.flank, 1, L)
        end = clamp(cut + args.flank, 1, L)

        outdir = os.path.join(args.outdir, "{}_{}".format(q, cut))
        os.makedirs(outdir, exist_ok=True)

        # Contig windows
        left = seq[beg-1:cut]
        right = seq[cut:end]
        write_fa(os.path.join(outdir, "contig_left.fa"), [("{}:{}-{}".format(q, beg, cut), left)])
        write_fa(os.path.join(outdir, "contig_right.fa"), [("{}:{}-{}".format(q, cut+1, end), right)])
        write_fa(os.path.join(outdir, "contig_window.fa"), [("{}:{}-{}".format(q, beg, end), seq[beg-1:end])])

        files_list = ["contig_left.fa", "contig_right.fa", "contig_window.fa"]

        # Multi-FASTA bundle for MSA
        multi_entries = [("contig_left|{}:{}-{}".format(q,beg,cut), left),
                         ("contig_right|{}:{}-{}".format(q,cut+1,end), right)]

        # Reference windows (if provided)
        prev_ref_txt = None; next_ref_txt = None
        if contig_ref_paf:
            prev_blk, next_blk = find_adjacent_blocks(q, cut, contig_ref_paf)
            if prev_blk:
                prev_ref_txt = "{}:{}-{} ({})".format(prev_blk.tname, min(prev_blk.tstart,prev_blk.tend), max(prev_blk.tstart,prev_blk.tend), prev_blk.strand)
                if ref.get(prev_blk.tname):
                    tmin = min(prev_blk.tstart, prev_blk.tend); tmax = max(prev_blk.tstart, prev_blk.tend)
                    anchor = tmax if prev_blk.tend>prev_blk.tstart else tmin
                    tbeg = clamp(anchor - args.flank, 1, prev_blk.tlen)
                    tend = clamp(anchor + args.flank, 1, prev_blk.tlen)
                    rseq = ref[prev_blk.tname][tbeg-1:tend]
                    write_fa(os.path.join(outdir, "ref_prev_window.fa"), [("{}:{}-{}".format(prev_blk.tname, tbeg, tend), rseq)])
                    files_list.append("ref_prev_window.fa")
                    multi_entries.append(("ref_prev|{}:{}-{}".format(prev_blk.tname, tbeg, tend), rseq))
            if next_blk:
                next_ref_txt = "{}:{}-{} ({})".format(next_blk.tname, min(next_blk.tstart,next_blk.tend), max(next_blk.tstart,next_blk.tend), next_blk.strand)
                if ref.get(next_blk.tname):
                    tmin = min(next_blk.tstart, next_blk.tend); tmax = max(next_blk.tstart, next_blk.tend)
                    anchor = tmin if next_blk.tend>next_blk.tstart else tmax
                    tbeg = clamp(anchor - args.flank, 1, next_blk.tlen)
                    tend = clamp(anchor + args.flank, 1, next_blk.tlen)
                    rseq = ref[next_blk.tname][tbeg-1:tend]
                    write_fa(os.path.join(outdir, "ref_next_window.fa"), [("{}:{}-{}".format(next_blk.tname, tbeg, tend), rseq)])
                    files_list.append("ref_next_window.fa")
                    multi_entries.append(("ref_next|{}:{}-{}".format(next_blk.tname, tbeg, tend), rseq))

        # Write the multi-FASTA for MSA tools
        write_multifasta(os.path.join(outdir, "msa_bundle.fa"), multi_entries)
        files_list.append("msa_bundle.fa")

        # Motif scan (arbitrary list)
        motif_counts = {}
        window_seq = seq[beg-1:end]
        counts = {}
        for m in motif_list:
            s = window_seq.upper()
            mm = m.upper()
            rr = rc(mm)
            cf = 0; i = 0
            while True:
                i = s.find(mm, i)
                if i == -1: break
                cf += 1; i += 1
            cr = 0; j = 0
            while True:
                j = s.find(rr, j)
                if j == -1: break
                cr += 1; j += 1
            counts[mm] = {"forward": cf, "reverse": cr}
        with open(os.path.join(outdir, "motifs.tsv"), "w") as fh:
            fh.write("motif\tforward\treverse\n")
            for m, dct in counts.items():
                fh.write("{}\t{}\t{}\n".format(m, dct["forward"], dct["reverse"]))
        files_list.append("motifs.tsv")
        motifs_text = ", ".join(["{}:{}/{}".format(m, d["forward"], d["reverse"]) for m, d in counts.items()]) if counts else "NA"

        # Reads near window plot + table (PAF)
        reads_rows = []
        if reads_paf:
            reads_rows = reads_overlapping_window(reads_paf, q, beg, end, max_reads=args.max_reads)
            with open(os.path.join(outdir, "reads_near.tsv"), "w") as fh:
                fh.write("read\tcontig_start\tcontig_end\tmapq\n")
                for rname, s0, e0, mq in reads_rows:
                    fh.write("{}\t{}\t{}\t{}\n".format(rname, s0, e0, mq))
            files_list.append("reads_near.tsv")
            plot_reads(q, cut, beg, end, reads_rows, os.path.join(outdir, "readplot.png"))
            files_list.append("readplot.png")

        # Depth track (BAM)
        if args.bam_reads_vs_contigs:
            depth, err = compute_depth_track(args.bam_reads_vs_contigs, q, beg, end, bin_size=args.bin_size)
            if depth:
                x, y = depth
                plot_depth(q, beg, end, x, y, os.path.join(outdir, "depthplot.png"), cut=cut)
                files_list.append("depthplot.png")
            else:
                with open(os.path.join(outdir, "depthplot.txt"), "w") as fh:
                    fh.write("Depth not available: {}\n".format(err))
                files_list.append("depthplot.txt")

        summary = {
            "contig": q, "cut": cut, "beg": beg, "end": end,
            "reads_plotted": len(reads_rows),
            "prev_ref": prev_ref_txt, "next_ref": next_ref_txt,
            "motifs_text": motifs_text,
            "files": files_list
        }
        with open(os.path.join(outdir, "summary.json"), "w") as fh:
            json.dump(summary, fh, indent=2)
        make_report_html(os.path.join(outdir, "report.html"), summary)

        index_entries.append({"contig": q, "cut": cut, "relpath": "{}/report.html".format(os.path.basename(outdir))})

    # Top-level index
    write_index_html(args.outdir, index_entries)

    print("[OK] Wrote per-break folders under: {}".format(args.outdir))
    print("[OK] Index: {}".format(os.path.join(args.outdir, "index.html")))

if __name__ == "__main__":
    main()
