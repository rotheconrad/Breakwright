#!/usr/bin/env python3
import sys, os, argparse, re
from collections import deque, defaultdict

def parse_args():
    ap = argparse.ArgumentParser(
        description="Annotate Breakwright break candidates with HiFiASM GFA graph evidence (simple text parsing), and export Bandage subgraphs."
    )
    ap.add_argument("--gfa", required=True, help="HiFiASM GFA file (e.g., *.p_ctg.gfa).")
    ap.add_argument("--breaks", required=True, help="Break table TSV from paf_breakfinder.py (must include qname, cut).")
    ap.add_argument("--outprefix", required=True, help="Output prefix for annotated TSV(s).")
    ap.add_argument("--end_proximity_bp", type=int, default=10000, help="Distance (bp) considered 'near end' (default: 10000).")
    ap.add_argument("--min_overlap_warn", type=int, default=50, help="Overlap (bp) below which an end is flagged as weak (default: 50).")
    ap.add_argument("--keep_only_matched", action="store_true", help="If set, do not write *_unmatched.tsv.")
    ap.add_argument("--sep", default="\t", help="Column delimiter (default: TAB).")

    # Subgraph export options (same names; richer .txt now)
    ap.add_argument("--subgraph_hops", type=int, default=2, help="BFS radius in link steps from qname (default: 2).")
    ap.add_argument("--subgraph_min_overlap", type=int, default=0, help="Only traverse links with >= this overlap (bp) (default: 0).")
    ap.add_argument("--subgraph_dir", help="Directory to write per-break node lists (default: <outprefix>_subgraphs).")

    # NEW: optional mini-GFA per break
    ap.add_argument("--emit_subgraph_gfa", action="store_true",
                    help="If set, also emit a mini-GFA (*.gfa) per break with only relevant S/L lines.")
    return ap.parse_args()

def parse_gfa(path):
    """
    Return:
        seg_len: dict seg -> length
        deg:     dict seg -> {'left':degree, 'right':degree}
        ovls:    dict seg -> {'left':[overlaps], 'right':[overlaps]}
        adj:     dict seg -> list of (neighbor, overlap_bp)
        s_lines: dict seg -> raw S line (for mini-GFA emit)
        l_lines: list of raw L lines (for mini-GFA emit)
    End-side logic:
      L a o1 b o2 ovl
        a end = right if o1 == '+' else left
        b end = left  if o2 == '+' else right
    Overlap parsed from CIGAR-like field (e.g., '123M').
    """
    seg_len = {}
    deg = {}
    ovls = {}
    adj = defaultdict(list)
    s_lines = {}
    l_lines = []

    ovl_re = re.compile(r"(\d+)M")

    def ensure(seg):
        if seg not in deg:
            deg[seg] = {'left':0, 'right':0}
            ovls[seg] = {'left':[], 'right':[]}
        if seg not in seg_len:
            seg_len[seg] = 0

    with open(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip("\n").split("\t")
            rec = f[0]
            if rec == 'S':
                seg = f[1]
                s_lines[seg] = line.rstrip("\n")
                # sequence or LN:i:length
                seq = f[2] if len(f) > 2 else '*'
                L = None
                if seq != '*':
                    L = len(seq)
                else:
                    for tok in f[3:]:
                        if tok.startswith("LN:i:"):
                            try:
                                L = int(tok.split(":",2)[2])
                            except Exception:
                                pass
                            break
                if L is None:
                    L = 0
                seg_len[seg] = L
                ensure(seg)
            elif rec == 'L':
                if len(f) < 6:
                    continue
                a, o1, b, o2, ovl = f[1], f[2], f[3], f[4], f[5]
                ensure(a); ensure(b)
                end_a = 'right' if o1 == '+' else 'left'
                end_b = 'left'  if o2 == '+' else 'right'
                deg[a][end_a] += 1
                deg[b][end_b] += 1
                m = ovl_re.search(ovl)
                ovllen = int(m.group(1)) if m else 0
                ovls[a][end_a].append(ovllen)
                ovls[b][end_b].append(ovllen)
                # undirected adjacency for subgraphs
                adj[a].append((b, ovllen))
                adj[b].append((a, ovllen))
                l_lines.append(line.rstrip("\n"))
    return seg_len, deg, ovls, adj, s_lines, l_lines

def read_tsv(path, sep="\t"):
    with open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(sep)
        rows = []
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split(sep)
            if len(f) < len(header):
                f += [''] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return header, rows

def write_tsv(path, header, rows, sep="\t"):
    with open(path, "wt") as fh:
        fh.write(sep.join(header) + "\n")
        for r in rows:
            fh.write(sep.join(str(r.get(h, "")) for h in header) + "\n")

def annotate_breaks(seg_len, deg, ovls, breaks_rows, end_prox_bp, min_ovl_warn, sep="\t"):
    out_rows = []
    unmatched_rows = []

    new_cols = [
        "qlen","dist_to_start","dist_to_end",
        "deg_left","deg_right","min_ovl_left","min_ovl_right",
        "nearest_end","nearest_end_degree","nearest_junction_bp","gfa_flag"
    ]

    for rec in breaks_rows:
        q = rec.get("qname")
        cut_str = rec.get("cut")
        try:
            cut = int(cut_str)
        except Exception:
            unmatched_rows.append(rec); continue

        if q not in seg_len:
            unmatched_rows.append(rec); continue

        qlen = seg_len[q]
        d_left = max(0, min(cut, qlen))
        d_right = max(0, max(0, qlen - cut))

        degL = deg.get(q, {}).get('left', 0)
        degR = deg.get(q, {}).get('right', 0)
        minOvL = min(ovls.get(q, {}).get('left', ["NA"])) if ovls.get(q, {}).get('left') else "NA"
        minOvR = min(ovls.get(q, {}).get('right', ["NA"])) if ovls.get(q, {}).get('right') else "NA"

        if d_left <= d_right:
            nearest_end = 'left'; nearest_deg = degL; nearest_dist = d_left
        else:
            nearest_end = 'right'; nearest_deg = degR; nearest_dist = d_right

        candidates = []
        if degL != 1: candidates.append(d_left)
        if degR != 1: candidates.append(d_right)
        nearest_junction_bp = min(candidates) if candidates else "NA"

        flag = "simple"
        if nearest_dist <= end_prox_bp:
            if nearest_deg == 0:
                flag = "tip_near_end"
            elif nearest_deg > 1:
                flag = "junction_near_end"
        weak_end = False
        if isinstance(minOvL, int) and minOvL != "NA" and minOvL < min_ovl_warn: weak_end = True
        if isinstance(minOvR, int) and minOvR != "NA" and minOvR < min_ovl_warn: weak_end = True
        if weak_end and flag == "simple":
            flag = "weak_overlap_end"

        rec2 = dict(rec)
        rec2.update({
            "qlen": qlen,
            "dist_to_start": d_left,
            "dist_to_end": d_right,
            "deg_left": degL,
            "deg_right": degR,
            "min_ovl_left": minOvL,
            "min_ovl_right": minOvR,
            "nearest_end": nearest_end,
            "nearest_end_degree": nearest_deg,
            "nearest_junction_bp": nearest_junction_bp,
            "gfa_flag": flag
        })
        out_rows.append(rec2)

    return new_cols, out_rows, unmatched_rows

def bfs_nodes_and_edges(adj, start, hops, min_ovl):
    """
    Undirected BFS by hops with optional min-overlap threshold.
    Returns:
        nodes: set of node IDs
        edges: set of (min(u,v), max(u,v), ovl_bp) for uniqueness
        deg:   dict node -> degree within subgraph
    """
    nodes = set([start]) if start in adj or start else set([start])
    edges = set()
    deg = defaultdict(int)
    if start not in adj:
        return nodes, edges, deg

    q = deque([(start, 0)])
    while q:
        node, d = q.popleft()
        if d == hops:
            continue
        for nb, ovl in adj.get(node, []):
            if ovl < min_ovl:
                continue
            a, b = (node, nb) if node <= nb else (nb, node)
            edges.add((a, b, ovl))
            deg[node] += 1
            deg[nb] += 1
            if nb not in nodes:
                nodes.add(nb)
                q.append((nb, d+1))
    return nodes, edges, deg

def write_mini_gfa(gfa_path, nodes, s_lines, l_lines):
    with open(gfa_path, "wt") as out:
        out.write("H\tVN:Z:1.0\n")
        sub = set(nodes)
        for n in sorted(sub):
            if n in s_lines:
                out.write(s_lines[n] + "\n")
        for ln in l_lines:
            f = ln.split("\t")
            if len(f) >= 6 and f[0] == "L":
                a, o1, b, o2 = f[1], f[2], f[3], f[4]
                if a in sub and b in sub:
                    out.write(ln + "\n")

def write_subgraphs(adj, breaks_rows, outprefix, subgraph_dir, hops, min_ovl, s_lines, l_lines, emit_gfa=False, sep="\t"):
    os.makedirs(subgraph_dir, exist_ok=True)
    index_rows = []
    for rec in breaks_rows:
        q = rec.get("qname")
        cut = rec.get("cut")
        if not q:
            continue
        nodes, edges, deg = bfs_nodes_and_edges(adj, q, hops, min_ovl)

        # Enriched .txt: include stats, junctions (deg>=3), node/edge lists
        txt_path = os.path.join(subgraph_dir, f"{q}_{cut}.txt")
        with open(txt_path, "wt") as fh:
            fh.write("# Breakwright subgraph dump\n")
            fh.write(f"contig: {q}\n")
            fh.write(f"position: {cut}\n")
            fh.write(f"nodes: {len(nodes)}\n")
            fh.write(f"edges: {len(edges)}\n")
            # degree stats
            if deg:
                vals = list(deg.values()); vals.sort()
                dmin = min(vals); dmed = vals[len(vals)//2]; dmean = sum(vals)/len(vals); dmax = max(vals)
            else:
                dmin = dmed = dmean = dmax = 0
            fh.write(f"degree_min: {dmin}\n")
            fh.write(f"degree_median: {dmed}\n")
            fh.write(f"degree_mean: {dmean:.2f}\n")
            fh.write(f"degree_max: {dmax}\n")
            junc = sorted([n for n,v in deg.items() if v >= 3])
            fh.write(f"junction_count (deg>=3): {len(junc)}\n")
            if junc:
                fh.write("junction_list:\n")
                for j in junc:
                    fh.write(f"  - {j}\n")
            fh.write("nodes_list:\n")
            for n in sorted(nodes):
                fh.write(f"  - {n}\n")
            fh.write("edges_list:\n")
            for a,b,ov in sorted(edges):
                fh.write(f"  - {a} -- {b} overlap_bp={ov}\n")

        gfa_path = ""
        if emit_gfa:
            gfa_path = os.path.join(subgraph_dir, f"{q}_{cut}.gfa")
            write_mini_gfa(gfa_path, nodes, s_lines, l_lines)

        index_rows.append({
            "qname": q,
            "cut": cut,
            "n_nodes": len(nodes),
            "hops": hops,
            "min_overlap": min_ovl,
            "node_file": txt_path,
            "subgraph_gfa": gfa_path,
            "nodes_csv": ",".join(sorted(nodes))
        })
    header = ["qname","cut","n_nodes","hops","min_overlap","node_file","subgraph_gfa","nodes_csv"]
    write_tsv(outprefix + "_subgraphs.tsv", header, index_rows, sep=sep)
    return outprefix + "_subgraphs.tsv"

def main():
    args = parse_args()

    seg_len, deg, ovls, adj, s_lines, l_lines = parse_gfa(args.gfa)
    header, rows = read_tsv(args.breaks, sep=args.sep)

    new_cols, annotated, unmatched = annotate_breaks(
        seg_len, deg, ovls, rows,
        end_prox_bp=args.end_proximity_bp,
        min_ovl_warn=args.min_overlap_warn,
        sep=args.sep
    )

    outpath = args.outprefix + "_breaks_gfa.tsv"
    final_header = list(header)
    for c in new_cols:
        if c not in final_header:
            final_header.append(c)
    write_tsv(outpath, final_header, annotated, sep=args.sep)

    if not args.keep_only_matched and unmatched:
        unmp = args.outprefix + "_unmatched.tsv"
        write_tsv(unmp, header, unmatched, sep=args.sep)

    subdir = args.subgraph_dir or (args.outprefix + "_subgraphs")
    sub_index = write_subgraphs(
        adj, annotated, args.outprefix, subdir,
        args.subgraph_hops, args.subgraph_min_overlap,
        s_lines, l_lines, emit_gfa=args.emit_subgraph_gfa, sep=args.sep
    )

    print(f"[OK] Annotated: {outpath}")
    if not args.keep_only_matched:
        print(f"[OK] Unmatched: {args.outprefix}_unmatched.tsv (n={len(unmatched)})")
    print(f"[OK] Subgraph index: {sub_index}")
    print(f"[OK] Node lists dir: {subdir}")

if __name__ == "__main__":
    main()
