#!/usr/bin/env python3
import argparse
import pandas as pd
from collections import defaultdict
import re
import os
import glob
from Bio import SeqIO

# --- Weights ---
PLASMID_WEIGHTS = {'conjugative': 1.0, 'mobilizable': 0.75, 'non-mobilizable': 0.5}
PHAGE_WEIGHTS = {'virulent': 1.0, 'temperate': 0.5}
TRANSPOSON_WEIGHTS = {'conjugative': 1.0, 'composite': 0.75, 'unit': 0.5}
DEFAULT_WEIGHT = 0.1

# WHO/priority pathogens list (lowercased)
PATHOGENS = {
    'escherichia coli',
    'staphylococcus aureus',
    'klebsiella pneumoniae',
    'streptococcus pneumoniae',
    'acinetobacter baumannii',
    'pseudomonas aeruginosa',
    'mycobacterium tuberculosis',
    'enterococcus faecium',
    'enterobacter asburiae',
    'enterobacter cancerogenus',
    'enterobacter chengduensis',
    'enterobacter cloacae',
    'enterobacter hormaechei',
    'enterobacter kobei',
    'enterobacter roggenkampii',
    'streptococcus agalactiae',
    'salmonella enterica',
    'enterococcus faecalis',
    'proteus columbae',
    'proteus mirabilis',
    'proteus penneri',
    'proteus vulgaris',
    'enterococcus avium',
    'enterococcus gilvus',
    'enterococcus hirae',
    'enterococcus pallens',
    'serratia liquefaciens',
    'serratia marcescens',
    'serratia odorifera',
    'serratia rubidaea',
    'streptococcus pyogenes',
    'citrobacter amalonaticus',
    'citrobacter freundii',
    'citrobacter koseri',
    'citrobacter portucalensis',
    'citrobacter werkmanii',
    'citrobacter youngae',
    'haemophilus influenzae',
    'shigella boydii',
    'shigella dysenteriae',
    'shigella flexneri',
    'shigella sonnei',
    'morganella morganii'
}


# ---------------- Normalization helpers ----------------
_NODE_RE = re.compile(r'(NODE_\d+)')

def normalize_id(x: str) -> str:
    """
    Normalize any contig-like id to 'NODE_<num>' (strip sample prefix and trailing length/cov).
    Works with:
      - 'ERRxxxx:NODE_123_length_9999_cov_3.0'
      - 'NODE_123_length_9999_cov_3.0'
      - 'NODE_123'
    """
    s = str(x).strip()
    if ':' in s:
        s = s.split(':', 1)[1]
    m = _NODE_RE.search(s)
    return m.group(1).lower() if m else s.lower()

def normalize_node_id(node_id: str) -> str:
    return normalize_id(node_id)

# ---------------- Robust readers ----------------
def _try_read_any(file):
    """
    Robust reader for BLAST/ARG/microbe result tables:
    - Try pandas sniffed header
    - If missing qseqid, retry with tab + no header
    - If tab parsing fails (irregular cols), fall back to python engine or manual split
    """
    import csv

    # 1) try sniffed header
    try:
        df = pd.read_csv(file, sep=None, engine="python")
        if "qseqid" in df.columns or "Accession" in df.columns:
            return df
    except Exception:
        pass

    # 2) try tab with header=None
    try:
        df = pd.read_csv(file, sep="\t", header=None, engine="c")
        base_cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
                     "qstart","qend","sstart","send","evalue","bitscore"]
        if df.shape[1] >= 12:
            extras = [f"extra{i}" for i in range(df.shape[1] - 12)]
            df.columns = base_cols + extras
        return df
    except pd.errors.ParserError:
        pass

    # 3) last resort: manual csv.reader with tab delimiter
    rows = []
    with open(file) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    max_cols = max(len(r) for r in rows)
    df = pd.DataFrame([r + [None]*(max_cols-len(r)) for r in rows])
    return df

# ---------------- ARG + Data Loaders ----------------
def load_arg(file):
    """
    Load ARG BLAST hits for contigs/plasmids/phages.
    Accepts .txt/.tsv/.csv with or without headers; also supports FASTA (ids only).
    Keeps 'Subtype' and 'C1 Score' if present; creates C1_score if missing (0.0).
    Adds normalized 'node' column (NODE_<num>).
    """
    if file.lower().endswith((".fasta", ".fa")):
        records = list(SeqIO.parse(file, "fasta"))
        df = pd.DataFrame({"qseqid": [r.id for r in records]})
    else:
        df = _try_read_any(file)
    # standardize headers
    df.columns = [c.strip() for c in df.columns]
    # rename if needed
    if "C1 Score" in df.columns and "C1_score" not in df.columns:
        df = df.rename(columns={"C1 Score": "C1_score"})
    if "Subtype" not in df.columns and "Subtype" in [c.strip() for c in df.columns]:
        pass  # already covered
    if "qseqid" not in df.columns:
        raise ValueError(f"'qseqid' column missing in {file}; found {df.columns.tolist()}")
    # normalize node
    df["node"] = df["qseqid"].apply(normalize_id)
    # ensure C1_score
    if "C1_score" not in df.columns:
        df["C1_score"] = 0.0
    return df

def compute_ARG_score(arg_df):
    """
    Equation (1): ARG_Score_contig = sum over unique subtypes (max C1 per subtype).
    If 'Subtype' missing, fall back to sseqid as grouping key.
    Returns:
      arg_score_per_node (dict NODE->float),
      arg_count_per_node (dict NODE->int)
    """
    if arg_df.empty:
        return {}, {}
    key = "Subtype" if "Subtype" in arg_df.columns else "sseqid"
    # max C1 per (node, subtype)
    grp = arg_df.groupby(["node", key], dropna=False)["C1_score"].max().reset_index()
    # sum across subtypes -> ARG score per node
    arg_score = grp.groupby("node")["C1_score"].sum()
    # raw ARG counts per node (rows/hits on that node)
    arg_counts = arg_df.groupby("node").size()
    return arg_score.to_dict(), arg_counts.to_dict()

def load_microbes(file):
    """
    Parse CAT output and extract NODE -> species mapping.
    Only species-level info (s__) is extracted, case-insensitive.
    """
    try:
        df = pd.read_csv(file, sep="\t", header=0)
    except Exception:
        df = _try_read_any(file)

    # Clean column names
    df.columns = df.columns.map(str)
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_')

    # Identify contig/node column robustly
    node_col = None
    for col in df.columns:
        # Remove leading '#' and '_' to normalize
        col_clean = col.lstrip('#').lstrip('_')
        if col_clean in ['node', 'contig', 'contig_id']:
            node_col = col
            break
    if node_col is None:
        raise ValueError(f"No contig/node column found in {file}, got {df.columns.tolist()}")

    if 'lineage' not in df.columns:
        raise ValueError(f"No lineage column found in {file}, got {df.columns.tolist()}")

    mapping = {}
    for _, row in df.iterrows():
        contig_id = normalize_node_id(str(row[node_col]))
        lineage = str(row['lineage']).strip()

        species = None
        for part in lineage.split(';'):
            part = part.strip()
            if part.startswith('s__'):
                species = part.replace('s__','').strip()
                break

        # fallback: last level
        if not species:
            parts = [p.strip() for p in lineage.split(';') if p.strip()]
            species = parts[-1] if parts else lineage

        mapping[contig_id] = species.lower()  # store lowercase for WHO matching

    return mapping

def compute_host_pathogenicity(nodes, species_map):
    """
    Assign HP score:
    - 1.0 if exact match (case-insensitive) to WHO priority pathogens
    - 0.1 otherwise
    """
    hp = {}
    for node in nodes:
        sp = species_map.get(node, "")
        hp[node] = 1.0 if sp in {p.lower() for p in PATHOGENS} else 0.1
    return hp

def load_transposons(file):
    """
    contig_transposon_blast.txt: qseqid (contig), sseqid like 'TnXXXX_CJ|...' where suffix encodes type:
      _CJ -> conjugative, _CO -> composite, _U -> unit
    Weight per contig = mean of subtype weights for its hits (if none -> DEFAULT)
    """
    df = pd.read_csv(file, sep="\t", header=None)
    # name columns flexibly (>=12)
    base_cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
                 'qstart','qend','sstart','send','evalue','bitscore']
    if df.shape[1] >= 12:
        extras = [f"extra{i}" for i in range(df.shape[1] - 12)]
        df.columns = base_cols + extras
    weight_map = defaultdict(list)
    for _, row in df.iterrows():
        contig = normalize_id(row['qseqid'])
        sseq = str(row['sseqid'])
        if '_CJ' in sseq:
            ttype = 'conjugative'
        elif '_CO' in sseq:
            ttype = 'composite'
        elif '_U' in sseq:
            ttype = 'unit'
        else:
            ttype = None
        weight_map[contig].append(TRANSPOSON_WEIGHTS.get(ttype, DEFAULT_WEIGHT))
    return {contig: (sum(w)/len(w) if w else DEFAULT_WEIGHT) for contig, w in weight_map.items()}

def load_ARGs_in_transposon_nodes(file):
    """
    For completeness: read transposon_blast (ARGs located in transposon nodes).
    We only need nodes; used for potential future filters.
    """
    if file.endswith(".txt"):
        df = pd.read_csv(file, sep="\t", header=None)
        base_cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
                     "qstart","qend","sstart","send","evalue","bitscore"]
        if df.shape[1] >= 12:
            extras = [f"extra{i}" for i in range(df.shape[1] - 12)]
            df.columns = base_cols + extras
    elif file.endswith((".fasta", ".fa")):
        records = list(SeqIO.parse(file, "fasta"))
        df = pd.DataFrame({"qseqid": [r.id for r in records]})
    else:
        raise ValueError(f"Unsupported ARG-in-transposon file format: {file}")
    df["node"] = df["qseqid"].apply(normalize_id)
    return df

def load_plasmids(file):
    # Let pandas automatically detect the separator
    df = pd.read_csv(file, sep=None, engine='python')  # auto-detect separator
    
    # Strip spaces from column names
    df.columns = [str(c).strip() for c in df.columns]
    
    # find columns with tolerant names
    contig_col = None
    for cand in ["contig_id", "contig", "contigs"]:
        if cand in df.columns:
            contig_col = cand
            break
    if contig_col is None:
        # case-sensitive variants seen in your example
        for cand in ["contig_id", "contig_id_y", "contig_id_x"]:
            if cand in df.columns:
                contig_col = cand
                break
    if contig_col is None:
        raise ValueError(f"No contig_id column in {file}; got {df.columns.tolist()}")

    mob_col = None
    for cand in ["predicted_mobility", "predicted_mobility_y", "predicted_mobility_x"]:
        if cand in df.columns:
            mob_col = cand
            break
    if mob_col is None:
        raise ValueError(f"No predicted_mobility column in {file}; got {df.columns.tolist()}")

    weights = defaultdict(lambda: DEFAULT_WEIGHT)
    for _, row in df.iterrows():
        ids = str(row[contig_col]).split(";")
        mobility = str(row[mob_col]).strip().lower()
        w = PLASMID_WEIGHTS.get(mobility, DEFAULT_WEIGHT)
        for cid in ids:
            sid = normalize_id(cid)
            if sid:
                weights[sid] = max(weights[sid], w)
    return weights

def load_phages(file):
    """
    PT.tsv: columns like: Accession, Length, TYPE, PhaTYPScore
    Use TYPE to map to virulent/temperate weights.
    """
    df = _try_read_any(file)
    df.columns = [c.strip() for c in df.columns]
    if "Accession" not in df.columns:
        raise ValueError(f"No 'Accession' in {file}; got {df.columns.tolist()}")
    # TYPE field for mapping to weights (fallback to 'PhaTYP' if present)
    type_col = "TYPE" if "TYPE" in df.columns else ("PhaTYP" if "PhaTYP" in df.columns else None)
    weights = defaultdict(lambda: DEFAULT_WEIGHT)
    for _, row in df.iterrows():
        acc = normalize_id(row["Accession"])
        typ = str(row.get(type_col, "")).strip().lower() if type_col else ""
        weights[acc] = PHAGE_WEIGHTS.get(typ, DEFAULT_WEIGHT)
    return weights

# ---------------- RANK + ERR ----------------
def assign_rank(row):
    """
    Rank rules (unchanged):
      1: ARGs + ≥2 MGE types + pathogen
      2: ARGs + ≥2 MGE type 
      3: ARGs + ≥1 MGE type
      4: ARGs only
      5: else
    """
    mge_count = int(row['Plasmid_weight'] > DEFAULT_WEIGHT) + \
                int(row['Phage_weight'] > DEFAULT_WEIGHT) + \
                int(row['Transposon_weight'] > DEFAULT_WEIGHT)
    is_pathogen = row['HP'] > DEFAULT_WEIGHT
    if row['ARGs'] > 0 and mge_count >= 2 and is_pathogen:
        return 1
    elif row['ARGs'] > 0 and mge_count >= 2:
        return 2
    elif row['ARGs'] > 0 and mge_count >= 1:
        return 3
    elif row['ARGs'] > 0:
        return 4
    return 5

def compute_ERR(args):
    """
    Core computation following Methodology_MobiRes:

      Eq(1): ARG_Score_node = sum_over_subtypes max(C1)   (no normalization)
      Eq(2): HP_node ∈ {1.0, 0.1}
      Eq(3): MI_node = (W_plasmid + W_phage + W_transposon) / 3
      Eq(4): ERR_node = ARG_Score_node * MI_node * HP_node
      Eq(5): RR_sample = sum_over_nodes (ARGs_node * ERR_node)

    args: dict with file paths and outdir
    """
    # ARG hits on contigs + for counting plasmid/phage ARGs on same nodes
    df_arg = load_arg(args["arg_contig"])
    plasmid_df = load_arg(args["arg_plasmid"])
    phage_df = load_arg(args["arg_phage"])

    # Species -> HP
    species_map = load_microbes(args["microbes"])
    node_set = set(df_arg["node"])
    hp_dict = compute_host_pathogenicity(node_set, species_map)

    # ARG score (Eq 1) and simple ARG counts
    arg_score_dict, arg_count_dict = compute_ARG_score(df_arg)

    # Mobility weights
    plasmid_weights = load_plasmids(args["mob"])
    phage_weights = load_phages(args["phage"])
    transposon_weights = load_transposons(args["transposon"])
    _ = load_ARGs_in_transposon_nodes(args["arg_transposon"])  # kept for completeness

    # union of nodes seen anywhere
    all_nodes = set(arg_count_dict) | set(plasmid_weights) | set(phage_weights) | set(transposon_weights) | node_set

    rows = []
    for node in sorted(all_nodes):
        arg_count = int(arg_count_dict.get(node, 0))
        arg_score = float(arg_score_dict.get(node, 0.0))

        # HP per node (Eq 2) – default 0.1 if unknown
        hp = float(hp_dict.get(node, 0.1))

        # Raw mobility weights per node (no count-scaling; Eq 3)
        plasmid_weight = float(plasmid_weights.get(node, DEFAULT_WEIGHT))
        phage_weight   = float(phage_weights.get(node,   DEFAULT_WEIGHT))
        trans_weight   = float(transposon_weights.get(node, DEFAULT_WEIGHT))

        mobility_index = (plasmid_weight + phage_weight + trans_weight) / 3.0  # Eq 3

        # ERR per node (Eq 4)
        err_score = arg_score * mobility_index * hp

        rows.append({
            'Node': node,
            'ARGs': arg_count,
            'ARG_Score': arg_score,
            'Plasmid_weight': plasmid_weight,
            'Phage_weight': phage_weight,
            'Transposon_weight': trans_weight,
            'Mobility_Index': mobility_index,
            'HP': hp,
            'ERR_Score': err_score
        })

    df = pd.DataFrame(rows).sort_values("ERR_Score", ascending=False).reset_index(drop=True)
    df['ResCon_Rank'] = df.apply(assign_rank, axis=1)

    # flags for subset files
    df['has_ARG'] = df['ARGs'] > 0
    df['has_plasmid'] = df['Plasmid_weight'] > DEFAULT_WEIGHT
    df['has_phage'] = df['Phage_weight'] > DEFAULT_WEIGHT
    df['has_transposon'] = df['Transposon_weight'] > DEFAULT_WEIGHT
    df['has_pathogen'] = df['HP'] > DEFAULT_WEIGHT

    # Ensure output dir exists
    os.makedirs(args["outdir"], exist_ok=True)

    # --- Save subsets ---
    df[df['has_ARG'] & (df['has_plasmid'] | df['has_phage'] | df['has_transposon'] | df['has_pathogen'])] \
        .to_csv(os.path.join(args["outdir"], "nodes_ARG_plus_any_MGE_or_pathogen12.csv"), index=False)

    df[df['has_ARG'] & df['has_plasmid'] & df['has_phage'] & df['has_transposon'] & df['has_pathogen']] \
        .to_csv(os.path.join(args["outdir"], "nodes_ARG_plus_all_MGEs_and_pathogen12.csv"), index=False)

    df.to_csv(os.path.join(args["outdir"], 'node_ERR_with_HP.csv'), index=False)

    df.head(20).to_csv(os.path.join(args["outdir"], 'top_20_ERR_nodes_with_HP.csv'), index=False)

    # --- Sample summary (Eq 5) ---
    # RR_sample = sum over nodes (ARGs_node * ERR_node)
    sample_rr = float((df['ARGs'] * df['ERR_Score']).sum())
    total_ARGs = int(df['ARGs'].sum())
    total_ERR = float(df['ERR_Score'].sum())
    total_ARG_score = float(df['ARG_Score'].sum())

    sample_rank_any = 0
    subset_any = df[df['has_ARG'] & (df['has_plasmid'] | df['has_phage'] | df['has_transposon'] | df['has_pathogen'])]
    if not subset_any.empty:
        sample_rank_any = int(round(subset_any['ResCon_Rank'].mean()))

    summary_df = pd.DataFrame([{
        'Total_ARGs': total_ARGs,
        'Total_ARG_Score': total_ARG_score,
        'Total_ERR_Score': total_ERR,
        'Sample_ERR_Score': sample_rr,
        'Sample_ResCon_Rank_ARG_plus_any_MGE_or_pathogen': sample_rank_any
    }])
    summary_df.to_csv(os.path.join(args["outdir"], 'sample_ERR_summary.csv'), index=False)

    # --- Print Summary ---
    print(f"\n--- Sample Level Summary ---")
    print(f"Total ARGs: {total_ARGs}")
    print(f"Total ARG Score: {total_ARG_score:.6f}")
    print(f"Total ERR Score (sum of node ERRs): {total_ERR:.6f}")
    print(f"Sample ERR Score (Eq 5): {sample_rr:.6f}")
    print(f"Sample ResCon Rank (ARG + any MGE or pathogen): {sample_rank_any}")

# ---------------- MAIN ----------------
if __name__ == "__main__":
    # You can optionally allow CLI args for single-sample runs (kept for completeness)
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--base_dir', help="Base directory containing *_results folders (default: /home/psg/Downloads/Dairy)")
    cli = parser.parse_args()

    base_dir = cli.base_dir if cli.base_dir else "/home/psg/Downloads/Wastewater"

    for sample_dir in glob.glob(os.path.join(base_dir, "*_results")):
        sample_name = os.path.basename(sample_dir).replace("_results", "")
        files_needed = {
            "arg_contig":    os.path.join(sample_dir, "contig_blast.csv"),
            "arg_plasmid":   os.path.join(sample_dir, "plasmid_blast.txt"),
            "arg_phage":     os.path.join(sample_dir, "phage_blast.txt"),
            "transposon":    os.path.join(sample_dir, "contig_transposon_blast.txt"),
            "arg_transposon":os.path.join(sample_dir, "transposon_blast.txt"),
            "mob":           os.path.join(base_dir, f"{sample_name}_mob.csv"),
            "phage":         os.path.join(base_dir, f"{sample_name}_PT.tsv"),
            "microbes":      os.path.join(base_dir, f"{sample_name}_CAT.txt"),
            "outdir":        os.path.join(sample_dir, "ERR_output")
        }
        missing = [f for f, path in files_needed.items() if f != "outdir" and not os.path.exists(path)]
        if missing:
            print(f"[MISSING FILES] {sample_name}: {', '.join(missing)}")
            continue
        print(f"[PROCESSING] {sample_name}")
        compute_ERR(files_needed)
