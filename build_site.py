#!/usr/bin/env python3
"""
build_site.py — build the static Klebsiella regulatory-network site.

Usage:
    python3 build_site.py PATH_TO_Merge_plus_TU [-o OUT] [--gff GFF_DIR_or_FILE]

Reads each *_net_predictions_TU file and writes:
    OUT/data/<accession>.json   (one compact network per genome)
    OUT/manifest.json           (genome list for the selector)
    OUT/index.html              (the viewer; copied from next to this script)

Optional gene names (--gff): give the folder (or a single file) with the
GenBank (GCA) GFF3 annotations. For each genome its GFF is matched by the
assembly number; CDS features map protein accession / locus_tag -> gene
symbol / locus / product, and each node gets a readable label + function.
Use the GenBank (GCA_...) GFF, not RefSeq (GCF_...): only GCA accessions
match the network node IDs.

MergePlusTU format (TSV, 7 columns, no header):
    regulator  target  ref_orgs  ref_ids  n_orgs(1-4|NULL)  strand(F|R)  tu_id
"""
import sys, os, re, csv, json, argparse, glob, shutil
from urllib.parse import unquote

SUFFIX = "_protein.faa_net_predictions_TU"

def accession(fname):
    b = os.path.basename(fname)
    return b[:-len(SUFFIX)] if b.endswith(SUFFIX) else b

def genome_num(acc):
    m = re.search(r'GC[AF]_?(\d+)', acc)
    return m.group(1) if m else acc

# ---------- GFF gene-name mapping ----------
def parse_gff(path):
    """CDS features -> {accession|locus_tag: (gene, locus, product)}."""
    idx = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            a = {}
            for kv in c[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1); a[k] = unquote(v)
            acc = a.get("Name") or a.get("protein_id")
            locus = a.get("locus_tag") or (a.get("ID", "").replace("cds-", "") or None)
            rec = (a.get("gene"), locus, a.get("product"))
            if acc:   idx[acc] = rec
            if locus: idx.setdefault(locus, rec)
    return idx

def label_for(nid, idx):
    """Return (display_label, product). Handles 'A/B' and 'A-B' complexes."""
    if "/" in nid or "-" in nid:
        labs, prods = [], []
        for tk in re.split(r'([/-])', nid):
            if tk in ("/", "-"):
                labs.append(tk); continue
            l, p = label_for(tk, idx)
            labs.append(l)
            if p: prods.append(p)
        seen = list(dict.fromkeys(prods))
        return "".join(labs), "; ".join(seen)
    rec = idx.get(nid)
    if not rec:
        return nid, ""
    gene, locus, prod = rec
    return (gene or locus or nid), (prod or "")

def find_gff(acc, gff_arg):
    if not gff_arg:
        return None
    num = genome_num(acc)
    if os.path.isdir(gff_arg):
        for f in os.listdir(gff_arg):
            if num in f and (f.endswith(".gff") or f.endswith(".gff3")):
                return os.path.join(gff_arg, f)
        return None
    return gff_arg if num in os.path.basename(gff_arg) else None

# ---------- network ----------
# ---------- Literature-based TF family table ----------
# Curated gene-symbol -> family from standard bacterial TF classifications
# (RegulonDB / EcoCyc DNA-binding-domain families). Used as a fallback for hubs
# that have a gene symbol (from GFF) but no Pfam/InterProScan annotation.
LIT_FAMILIES = {
    "crp":"CRP/FNR", "fnr":"CRP/FNR",
    "fis":"Fis",
    "fur":"Fur", "zur":"Fur",
    "cytr":"LacI", "cra":"LacI", "galr":"LacI", "gals":"LacI", "purr":"LacI", "laci":"LacI",
    "arca":"OmpR/RR", "phop":"OmpR/RR", "phob":"OmpR/RR", "yycf":"OmpR/RR", "walr":"OmpR/RR",
    "ompr":"OmpR/RR", "cpxr":"OmpR/RR", "kdpe":"OmpR/RR", "baer":"OmpR/RR", "creb":"OmpR/RR",
    "narl":"NarL/FixJ", "narp":"NarL/FixJ", "uhpa":"NarL/FixJ", "fixj":"NarL/FixJ", "rcsb":"NarL/FixJ",
    "lrp":"Lrp/AsnC", "asnc":"Lrp/AsnC",
    "lexa":"LexA",
    "soxs":"AraC/XylS", "mara":"AraC/XylS", "rob":"AraC/XylS", "arac":"AraC/XylS", "xyls":"AraC/XylS",
    "nsrr":"Rrf2", "iscr":"Rrf2",
    "ihfa":"IHF/HU", "ihfb":"IHF/HU", "hupa":"IHF/HU", "hupb":"IHF/HU",
    "hns":"H-NS", "stpa":"H-NS",
    "lysr":"LysR", "oxyr":"LysR", "cysb":"LysR", "ilvy":"LysR", "metr":"LysR", "nac":"LysR",
    "tetr":"TetR", "acrr":"TetR", "ramr":"TetR", "envr":"TetR",
    "gntr":"GntR", "fadr":"GntR",
    "marr":"MarR", "mprа":"MarR", "slya":"MarR", "hpar":"MarR",
    "dnaa":"DnaA", "fur_2":"Fur",
    "rpoD":"Sigma", "rpoS":"Sigma", "rpoN":"Sigma", "rpoH":"Sigma", "fliA":"Sigma",
    "csgd":"LuxR", "malt":"LuxR",
    "argr":"ArgR", "trpr":"TrpR", "rob_1":"AraC/XylS",
}
def norm_symbol(s):
    import re as _re
    return _re.sub(r'[_\-.]\d+$', '', s).lower()

# ---------- Pfam / InterProScan family mapping ----------
# Shorten common Pfam TF-family descriptions to compact family labels.
def short_family(desc):
    d = desc.lower()
    table = [
        ("crp family","CRP/FNR"), ("fis family","Fis"), ("laci family","LacI"),
        ("ferric uptake","Fur"), ("cyclic nucleotide","CRP/FNR"),
        ("response regulator","OmpR/RR"), ("marr","MarR"), ("tetr","TetR"),
        ("lysr","LysR"), ("arac","AraC"), ("gntr","GntR"), ("asnc","Lrp/AsnC"),
        ("lrp","Lrp/AsnC"), ("ihf","IHF/HU"), ("histone-like","IHF/HU"),
        ("sigma","Sigma"), ("two-component","OmpR/RR"), ("rok","ROK"),
        ("deor","DeoR"), ("iclr","IclR"), ("padr","PadR"), ("merr","MerR"),
        ("cold-shock","CSD"), ("ros/mucr","Ros/MucR"),
    ]
    for key, lab in table:
        if key in d:
            return lab
    return desc.split(",")[0].strip()[:24]

def parse_interproscan(path):
    """InterProScan TSV -> {protein_acc: (pfam_acc, family_label)} (first Pfam hit)."""
    fam = {}
    with open(path) as fh:
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            acc, analysis, sig, desc = c[0], c[3], c[4], c[5]
            if acc in fam:
                continue
            if analysis.lower().startswith("pfam") and sig.startswith("PF"):
                fam[acc] = (sig, short_family(desc), desc)
    return fam

def find_ipr(acc, ipr_arg):
    if not ipr_arg:
        return None
    num = genome_num(acc)
    if os.path.isdir(ipr_arg):
        for f in os.listdir(ipr_arg):
            if num in f and f.endswith(".tsv"):
                return os.path.join(ipr_arg, f)
        return None
    return ipr_arg if num in os.path.basename(ipr_arg) else None

def build_genome(path):
    outdeg, indeg = {}, {}
    tfs, order = set(), []
    idx = {}
    edges = []
    def node(x):
        if x not in idx:
            idx[x] = len(order); order.append(x)
        return idx[x]
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) != 7:
                continue
            tf, tg, _o, _i, n, _s, _tu = row
            si, ti = node(tf), node(tg)
            tfs.add(tf)
            edges.append((si, ti, -1 if n == "NULL" else int(n)))
            if tf != tg:
                outdeg[tf] = outdeg.get(tf, 0) + 1
                indeg[tg] = indeg.get(tg, 0) + 1
    if not order:
        return None, None
    def tier(x):
        if x in tfs:
            return 0 if indeg.get(x, 0) == 0 else 1
        return 2
    nodes = [[x, 1 if x in tfs else 0, outdeg.get(x, 0), indeg.get(x, 0), tier(x)]
             for x in order]
    g = {"n": nodes, "e": [[s, t, c] for (s, t, c) in edges],
         "stats": {"nodes": len(order), "tfs": len(tfs),
                   "targets": len(order) - len(tfs), "edges": len(edges)}}
    return g, order

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("-o", "--out", default="klebsiella-red")
    ap.add_argument("--gff", default="", help="GFF dir or file (GenBank/GCA) for gene names")
    ap.add_argument("--pfam", default="", help="InterProScan TSV dir or file for TF-family (Pfam) of hubs")
    ap.add_argument("--hubtop", type=int, default=15, help="number of top hubs for the family heatmap")
    ap.add_argument("--lit", action="store_true", help="fill hub families from the built-in literature table when Pfam is missing (needs --gff gene symbols)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.src, "*" + SUFFIX)))
    if args.only:
        want = set(args.only.split(",")); files = [f for f in files if accession(f) in want]
    if args.limit:
        files = files[:args.limit]

    datadir = os.path.join(args.out, "data"); os.makedirs(datadir, exist_ok=True)
    manifest = []
    for i, f in enumerate(files, 1):
        acc = accession(f)
        g, order = build_genome(f)
        if g is None:
            continue
        g["genome"] = acc
        named = 0
        gff = find_gff(acc, args.gff)
        if gff:
            gidx = parse_gff(gff)
            names = {}
            for nid in order:
                lab, prod = label_for(nid, gidx)
                if lab != nid or prod:
                    names[nid] = [lab, prod]; named += 1
            if names:
                g["names"] = names
        # --- TF family of the top hubs (for the family heatmap) ---
        # Primary source: Pfam (InterProScan). Fallback: literature table by gene symbol.
        hubfam = None
        _ipr = find_ipr(acc, args.pfam) if args.pfam else None
        fam = parse_interproscan(_ipr) if _ipr else {}
        names_map = g.get("names", {})
        if fam or (args.lit and names_map):
            rows = sorted(g["n"], key=lambda r: r[2], reverse=True)[:args.hubtop]  # by out-degree
            hf = []
            for r in rows:
                nid = r[0]; disp = names_map.get(nid, [nid])[0]
                pf = fam.get(nid)
                if pf:
                    hf.append([disp, nid, r[2], pf[1], pf[0], "pfam"])
                elif args.lit and norm_symbol(disp) in LIT_FAMILIES:
                    hf.append([disp, nid, r[2], LIT_FAMILIES[norm_symbol(disp)], None, "lit"])
                else:
                    hf.append([disp, nid, r[2], None, None, None])
            if any(h[3] for h in hf):
                hubfam = hf
                g["hubfam"] = hf

        with open(os.path.join(datadir, acc + ".json"), "w") as out:
            json.dump(g, out, separators=(",", ":"))
        manifest.append({"id": acc, "tfs": g["stats"]["tfs"],
                         "targets": g["stats"]["targets"], "edges": g["stats"]["edges"],
                         "named": named, "hubfam": hubfam})
        if i % 200 == 0:
            print(f"  ...{i}/{len(files)}", file=sys.stderr)

    manifest.sort(key=lambda m: m["id"])
    with open(os.path.join(args.out, "manifest.json"), "w") as out:
        json.dump({"genomes": manifest}, out, separators=(",", ":"))

    here = os.path.dirname(os.path.abspath(__file__))
    src_html = os.path.join(here, "index.html")
    dst_html = os.path.join(args.out, "index.html")
    if os.path.exists(src_html) and os.path.abspath(src_html) != os.path.abspath(dst_html):
        shutil.copy(src_html, dst_html)
    print(f"Done: {len(manifest)} genomes -> {args.out}/", file=sys.stderr)

if __name__ == "__main__":
    main()
