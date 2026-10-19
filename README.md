# Klebsiella regulatory network — static site

Viewer for the predicted regulatory networks (TF to regulated genes), one per
genome. Static site: just files, no backend.

## Structure
    klebsiella-red/
      index.html        the viewer (Cytoscape.js + JSZip, from CDN)
      manifest.json     genome list for the selector
      data/             one network per genome: <accession>.json
      build_site.py     generator (rebuilds data/ + manifest.json)

## Generate data
    python3 build_site.py  PATH/Merge_plus_TU  -o klebsiella-red

### Gene names (optional)
Add gene symbols / locus / function with the GenBank (GCA) GFF annotations:

    python3 build_site.py  PATH/Merge_plus_TU  -o klebsiella-red  --gff PATH/gff_folder

- Use the GCA_... (GenBank) GFF, NOT GCF_... (RefSeq): only GCA accessions match
  the network node IDs. On a test genome this gives ~100% coverage (~86% a real
  gene symbol, the rest a locus tag) plus a function for every node.
- Each genome GFF is matched by its assembly number.

Other options: --limit N (first N genomes) ; --only A,B (only those accessions)

## View locally (loads by fetch, so it needs a server)
    cd klebsiella-red && python3 -m http.server 8000   # http://localhost:8000

## Publish
    scp -r klebsiella-red  user@SERVER:/var/www/html/klebsiella   # -> /klebsiella/

## Viewer features
- Genome selector (searchable) + "Load local data" (.zip / files, in-browser,
  session only).
- Color by: Hierarchy / Hub / Confidence. Views: TF-TF backbone, Top 10 hubs
  (blue-to-orange rank ramp).
- Trace downstream cascade from any TF (button in the detail panel): shows the
  full TF to target chain, laid out in tiers and colored by depth level.
- Confidence filter: hide weaker edges (All / >=1 / >=2 / >=3 reference organisms;
  ">=1" drops TU-only predictions).
- Betweenness ("Bottleneck" color mode): marks relay/bottleneck regulators by the
  share of shortest paths through them (Brandes, computed in the browser); shown
  in the detail panel and the Nodes CSV.
- Compare strains: pick a second genome, open a TF with a gene symbol, and see
  its regulon conserved vs strain-specific (matched by gene symbol).
- Export: high-res PNG of the current view, and Edges / Nodes CSV (with gene
  names), respecting the current view and confidence filter.
- Progressive labels: gene name / locus appears as you zoom (biggest nodes first).
- Network topology panel (nodes, edges, components, giant component, degree
  distribution, etc.).
- Hub family heatmap (tab): rows = organisms, columns = hub rank, cell colored
  by TF family (Pfam). Needs InterProScan: regenerate with --pfam <dir>.
  Only genomes with InterProScan get family data.
- --lit: fill hub families from a built-in literature table (known E. coli TF
  classification) by gene symbol, for genomes that have GFF names but no Pfam.
  Heatmap marks literature cells with a dashed border vs solid for Pfam.
- Hierarchy colors: TFs are green (master = dark, intermediate = light), targets grey.
- With --gff: node label = gene name (accession fallback); function + accession
  in the detail panel; search by gene name, locus or accession.
