#!/usr/bin/env python3
import argparse
import pandas as pd
from collections import defaultdict
import re
import os
import glob
from Bio import SeqIO

# --- Weights ---
PLASMID_WEIGHTS = {'conjugative': 1.0, 'mobilizable': 0.5, 'non-mobilizable': 0.1}
PHAGE_WEIGHTS = {'temperate': 1.0, 'virulent': 0.1}
TRANSPOSON_WEIGHTS = {'conjugative': 1.0, 'composite': 0.5, 'unit': 0.1}
IE_WEIGHTS = {'ICE': 1.0, 'IME': 0.5, 'CIME': 0.1}
DEFAULT_WEIGHT = 0.05

# WHO priority pathogens (lowercased)
PATHOGENS = {
    'escherichia coli','staphylococcus aureus','klebsiella pneumoniae','streptococcus pneumoniae',
    'acinetobacter baumannii','pseudomonas aeruginosa','mycobacterium tuberculosis',
    'enterococcus faecium','enterobacter asburiae','enterobacter cancerogenus','enterobacter chengduensis',
    'enterobacter cloacae','enterobacter hormaechei','enterobacter kobei','enterobacter roggenkampii',
    'streptococcus agalactiae','salmonella enterica','enterococcus faecalis','proteus columbae',
    'proteus mirabilis','proteus penneri','proteus vulgaris','enterococcus avium','enterococcus gilvus',
    'enterococcus hirae','enterococcus pallens','serratia liquefaciens','serratia marcescens',
    'serratia odorifera','serratia rubidaea','streptococcus pyogenes','citrobacter amalonaticus',
    'citrobacter freundii','citrobacter koseri','citrobacter portucalensis','citrobacter werkmanii',
    'citrobacter youngae','haemophilus influenzae','shigella boydii','shigella dysenteriae',
    'shigella flexneri','shigella sonnei','neisseria gonorrhoeae','morganella morganii'
}

# ---------------- Helpers ----------------
_NODE_RE = re.compile(r'(NODE_\d+)')

def normalize_node(x):
    s = str(x).strip()
    s = s.split(":")[0]
    s = s.split("_length_")[0]
    m = re.search(r'(NODE_\d+)', s, re.I)
    return m.group(1).upper() if m else s.upper()

def _try_read_any(file):
    import csv
    try:
        df = pd.read_csv(file, sep=None, engine="python")
        if "qseqid" in df.columns or "Accession" in df.columns:
            return df
    except Exception: pass
    try:
        df = pd.read_csv(file, sep="\t", header=None, engine="c")
        base_cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
                     "qstart","qend","sstart","send","evalue","bitscore"]
        if df.shape[1] >= 12:
            extras = [f"extra{i}" for i in range(df.shape[1] - 12)]
            df.columns = base_cols + extras
        return df
    except pd.errors.ParserError: pass
    rows=[]
    with open(file) as fh:
        reader=csv.reader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    max_cols=max(len(r) for r in rows)
    df=pd.DataFrame([r+[None]*(max_cols-len(r)) for r in rows])
    return df

# ---------------- ARG loaders ----------------
def load_arg(file):
    if file.lower().endswith((".fasta",".fa")):
        records=list(SeqIO.parse(file,"fasta"))
        df=pd.DataFrame({"qseqid":[r.id for r in records]})
    else:
        df=_try_read_any(file)
    df.columns=[c.strip() for c in df.columns]
    if "C1 Score" in df.columns and "C1_score" not in df.columns:
        df.rename(columns={"C1 Score":"C1_score"}, inplace=True)
    if "qseqid" not in df.columns:
        raise ValueError(f"'qseqid' column missing in {file}")
    df["node"]=df["qseqid"].apply(normalize_node)
    if "C1_score" not in df.columns:
        df["C1_score"]=0.0
    return df

def compute_ARG_score(arg_df):
    if arg_df.empty: return {}, {}
    key="Subtype" if "Subtype" in arg_df.columns else "sseqid"
    grp=arg_df.groupby(["node",key],dropna=False)["C1_score"].max().reset_index()
    arg_score=grp.groupby("node")["C1_score"].sum()
    arg_counts=arg_df.groupby("node").size()
    return arg_score.to_dict(), arg_counts.to_dict()

# ---------------- Microbes ----------------
def load_microbes(file):
    try: df=pd.read_csv(file, sep="\t", header=0)
    except Exception: df=_try_read_any(file)
    df.columns=df.columns.map(str).str.strip().str.lower().str.replace(' ','_')
    node_col=None
    for col in df.columns:
        col_clean=col.lstrip('#').lstrip('_')
        if col_clean in ['node','contig','contig_id']:
            node_col=col
            break
    if node_col is None:
        raise ValueError(f"No contig/node column found in {file}")
    if 'lineage' not in df.columns:
        raise ValueError(f"No lineage column found in {file}")
    mapping={}
    for _, row in df.iterrows():
        contig_id=normalize_node(str(row[node_col]))
        lineage=str(row['lineage']).strip()
        species=None
        for part in lineage.split(';'):
            part=part.strip()
            if part.startswith('s__'):
                species=part.replace('s__','').strip()
                break
        if not species:
            parts=[p.strip() for p in lineage.split(';') if p.strip()]
            species=parts[-1] if parts else lineage
        mapping[contig_id]=species.lower()
    return mapping

def compute_host_pathogenicity(nodes, species_map):
    hp={}
    for node in nodes:
        sp=species_map.get(node,"")
        hp[node]=1.0 if sp in PATHOGENS else 0.1
    return hp

# ---------------- Transposons ----------------
def load_transposons(file,arg_nodes=None):
    df=pd.read_csv(file, sep="\t", header=None)
    base_cols=['qseqid','sseqid','pident','length','mismatch','gapopen',
               'qstart','qend','sstart','send','evalue','bitscore']
    if df.shape[1]>=12:
        extras=[f"extra{i}" for i in range(df.shape[1]-12)]
        df.columns=base_cols+extras
    weight_map=defaultdict(list)
    for _, row in df.iterrows():
        contig=normalize_node(row['qseqid'])
        sseq=str(row['sseqid'])
        if '_CJ' in sseq: ttype='conjugative'
        elif '_CO' in sseq: ttype='composite'
        elif '_U' in sseq: ttype='unit'
        else: ttype=None
        weight_map[contig].append(TRANSPOSON_WEIGHTS.get(ttype,DEFAULT_WEIGHT))
    return {node:(sum(w)/len(w) if w else DEFAULT_WEIGHT) for node,w in weight_map.items()}

# ---------------- Plasmid / Phage ----------------
def load_plasmids(mob_file, plasmid_arg_file):
    df_arg=_try_read_any(plasmid_arg_file)
    df_arg["node"]=df_arg.iloc[:,0].apply(normalize_node)
    arg_count=df_arg["node"].value_counts().to_dict()
    total_args=sum(arg_count.values())
    df=pd.read_csv(mob_file, sep=None, engine="python")
    df.columns=[c.strip() for c in df.columns]
    contig_col=next(c for c in df.columns if "contig" in c.lower())
    mob_col=next(c for c in df.columns if "mob" in c.lower())
    weights={}
    for _,row in df.iterrows():
        node=normalize_node(row[contig_col])
        if node not in arg_count or total_args==0:
            weights[node]=DEFAULT_WEIGHT
            continue
        mobility=str(row[mob_col]).lower().strip()
        mob_score=PLASMID_WEIGHTS.get(mobility, DEFAULT_WEIGHT)
        weights[node]=(arg_count[node]/total_args)*mob_score
    return weights

def load_phages(phage_type_file, phage_arg_file):
    df_arg=_try_read_any(phage_arg_file)
    df_arg["node"]=df_arg.iloc[:,0].apply(normalize_node)
    arg_count=df_arg["node"].value_counts().to_dict()
    total_args=sum(arg_count.values())
    df=pd.read_csv(phage_type_file, sep=None, engine="python")
    df.columns=[c.strip() for c in df.columns]
    acc_col="Accession"
    type_col="TYPE"
    weights={}
    for _,row in df.iterrows():
        node=normalize_node(row[acc_col])
        if node not in arg_count or total_args==0:
            weights[node]=DEFAULT_WEIGHT
            continue
        ptype=str(row[type_col]).lower().strip()
        life_score=PHAGE_WEIGHTS.get(ptype,DEFAULT_WEIGHT)
        weights[node]=(arg_count[node]/total_args)*life_score
    return weights

# ----------------------------
# Load IE hits
# ----------------------------
def load_IE_hits(file):
    """
    Load IE hits per node.
    Columns expected: NODE_ID, IE_name, start, end, identity
    Adds:
      - 'node' normalized NODE_ID
      - 'IE_type' (ICE, IME, CIME)
    """
    df = pd.read_csv(file, sep="\t", header=None,
                     names=["NODE_ID", "IE_name", "start", "end", "identity"])
    df["node"] = df["NODE_ID"].apply(normalize_node)

    def classify_ie(ie_name):
        ie_name = str(ie_name)
        if "CIME" in ie_name:
            return "CIME"
        elif "_IME" in ie_name or "|IME|" in ie_name or "IME" in ie_name:
            return "IME"
        else:
            return "ICE"

    df["IE_type"] = df["IE_name"].apply(classify_ie)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df

# ----------------------------
# Load ARG hits in IE (robust)
# ----------------------------
def load_IE_ARG(file):
    """
    Load ARG hits with coordinates for IE nodes.
    Columns: qseqid with format NODE:START-END, ARG_info
    Adds:
      - 'node', 'ARG_start', 'ARG_end'
    """
    df = pd.read_csv(file, sep="\t", header=None, dtype=str, engine="python")
    df['node'] = df[0].apply(normalize_node)

    # Safe coordinate parsing
    df['ARG_start'] = 0
    df['ARG_end'] = 0
    for i, row in df.iterrows():
        try:
            # Only parse if ':' exists
            if ':' in row[0] and '-' in row[0]:
                coords = row[0].split(":")[1].split("-")
                df.at[i,'ARG_start'] = int(coords[0].strip())
                df.at[i,'ARG_end'] = int(coords[1].strip())
            else:
                df.at[i,'ARG_start'] = 0
                df.at[i,'ARG_end'] = 0
        except Exception:
            df.at[i,'ARG_start'] = 0
            df.at[i,'ARG_end'] = 0

    df['ARG_info'] = df[1]
    return df

# ----------------------------
# Compute IE weights
# ----------------------------
def compute_IE_weights(ie_df, arg_df, total_arg_count):
    """
    Correct IE weight:
    (Σ ARGs_in_IE_type × IE_weight + non_IE_ARGs × DEFAULT_WEIGHT)
    / total_ARGs_in_contig
    """

    weights = {}

    for node in arg_df["node"].unique():

        node_args = arg_df[arg_df["node"] == node]
        node_ies = ie_df[ie_df["node"] == node]

        if node_args.empty or total_arg_count == 0:
            weights[node] = DEFAULT_WEIGHT
            continue

        # Initialize ARG counters
        arg_counts = {"ICE": 0, "IME": 0, "CIME": 0}
        args_in_ie = set()

        # Count ARGs overlapping each IE
        for _, ie in node_ies.iterrows():
            overlapping_args = node_args[
                (node_args["ARG_start"] <= ie["end"]) &
                (node_args["ARG_end"] >= ie["start"])
            ]

            for idx in overlapping_args.index:
                arg_counts[ie["ie_type"]] += 1
                args_in_ie.add(idx)

        total_node_args = len(node_args)
        non_ie_args = total_node_args - len(args_in_ie)

        weighted_sum = (
            arg_counts["ICE"]  * 1.0 +
            arg_counts["IME"]  * 0.5 +
            arg_counts["CIME"] * 0.1 +
            non_ie_args       * DEFAULT_WEIGHT
        )

        weights[node] = weighted_sum / total_arg_count

    return weights

# ---------------- Rank + RR ----------------
def assign_rank(row):
    mge_count=int(row['Plasmid_weight']>DEFAULT_WEIGHT)+\
              int(row['Phage_weight']>DEFAULT_WEIGHT)+\
              int(row['Transposon_weight']>DEFAULT_WEIGHT)+\
              int(row['IE_weight']>DEFAULT_WEIGHT)
    is_pathogen=row['HP']>DEFAULT_WEIGHT
    if row['ARGs']>0 and mge_count>=2 and is_pathogen: return 1
    elif row['ARGs']>0 and mge_count>=2: return 2
    elif row['ARGs']>0 and mge_count>=1: return 3
    elif row['ARGs']>0: return 4
    return 5

def compute_RR(files_needed):
    
    df_arg = load_arg(files_needed["arg_contig"])
    arg_nodes_set = set(df_arg["node"])
    arg_score_dict, arg_count_dict = compute_ARG_score(df_arg)

    
    species_map = load_microbes(files_needed["microbes"])
    node_set = set(df_arg["node"])
    hp_dict = compute_host_pathogenicity(node_set, species_map)

    
    plasmid_weights = load_plasmids(files_needed["mob"], files_needed["arg_plasmid"])
    phage_weights = load_phages(files_needed["phage"], files_needed["arg_phage"])
    transposon_weights = load_transposons(files_needed["transposon"], arg_nodes=arg_nodes_set)

   
    ie_hits = load_IE_hits(files_needed["ie"])       # load IE hits per node
    ie_arg_df = load_IE_ARG(files_needed["ie_arg"])  # ARGs inside IEs
    total_ARGs_in_contig = int(df_arg.shape[0])
    ie_weights = compute_IE_weights(
        ie_hits,
        ie_arg_df,
        total_ARGs_in_contig
    )


    
    all_nodes = set(arg_count_dict) | set(plasmid_weights) | set(phage_weights) | \
                set(transposon_weights) | set(ie_weights) | node_set

    rows = []
    for node in sorted(all_nodes):
        arg_count = int(arg_count_dict.get(node, 0))
        arg_score = float(arg_score_dict.get(node, 0.0))
        hp = float(hp_dict.get(node, 0.1))
        plasmid_weight = float(plasmid_weights.get(node, DEFAULT_WEIGHT))
        phage_weight = float(phage_weights.get(node, DEFAULT_WEIGHT))
        trans_weight = float(transposon_weights.get(node, DEFAULT_WEIGHT))
        ie_weight = float(ie_weights.get(node, DEFAULT_WEIGHT))

        mobility_index = (plasmid_weight + phage_weight + trans_weight + ie_weight) / 4.0
        rr_score = arg_score * mobility_index * hp

        rows.append({
            'Node': node,
            'ARGs': arg_count,
            'ARG_Score': arg_score,
            'Plasmid_weight': plasmid_weight,
            'Phage_weight': phage_weight,
            'Transposon_weight': trans_weight,
            'IE_weight': ie_weight,
            'Mobility_Index': mobility_index,
            'HP': hp,
            'RR_Score': rr_score
        })

   
    df = pd.DataFrame(rows).sort_values("RR_Score", ascending=False).reset_index(drop=True)
    df['ResCon_Rank'] = df.apply(assign_rank, axis=1)

    os.makedirs(files_needed["outdir"], exist_ok=True)

    
    # Nodes with ARG + any MGE or pathogen
    df_any = df[
        (df['ARGs'] > 0) &
        ((df['Plasmid_weight'] > DEFAULT_WEIGHT) |
        (df['Phage_weight'] > DEFAULT_WEIGHT) |
        (df['Transposon_weight'] > DEFAULT_WEIGHT) |
        (df['HP'] > DEFAULT_WEIGHT))
    ]
    df_any.to_csv(os.path.join(files_needed["outdir"], "nodes_ARG_plus_any_MGE_or_pathogen.csv"), index=False)

    # Nodes with ARG + all MGEs + pathogen
    df_all = df[
        (df['ARGs'] > 0) &
        (df['Plasmid_weight'] > DEFAULT_WEIGHT) &
        (df['Phage_weight'] > DEFAULT_WEIGHT) &
        (df['Transposon_weight'] > DEFAULT_WEIGHT) &
        (df['HP'] > DEFAULT_WEIGHT)
    ]
    df_all.to_csv(os.path.join(files_needed["outdir"], "nodes_ARG_plus_all_MGEs_and_pathogen.csv"), index=False)

    df.to_csv(os.path.join(files_needed["outdir"], 'node_RR_with_HP.csv'), index=False)
    df.head(20).to_csv(os.path.join(files_needed["outdir"], 'top_20_RR_nodes_with_HP.csv'), index=False)
    
    
    sample_rr = float((df['ARGs'] * df['RR_Score']).sum())
    total_ARGs = int(df['ARGs'].sum())
    total_RR = float(df['RR_Score'].sum())
    total_ARG_score = float(df['ARG_Score'].sum())
    
    sample_rank_any = int(round(df_any['ResCon_Rank'].mean())) if not df_any.empty else 0
    summary_df = pd.DataFrame([{
        'Total_ARGs': total_ARGs,
        'Total_ARG_Score': total_ARG_score,
        'Total_RR_Score': total_RR,
        'Sample_RR_Score': sample_rr,
        'Sample_ResCon_Rank_ARG_plus_any_MGE_or_pathogen': sample_rank_any
    }])
    summary_df.to_csv(os.path.join(files_needed["outdir"], 'sample_RR_summary.csv'), index=False)

    
    print(f"\n--- Sample Level Summary ---")
    print(f"Total ARGs: {total_ARGs}")
    print(f"Total ARG Score: {total_ARG_score:.6f}")
    print(f"Total RR Score (sum of node RRs): {total_RR:.6f}")
    print(f"Sample RR Score (Eq 5): {sample_rr:.6f}")
    print(f"Sample ResCon Rank (ARG + any MGE or pathogen): {sample_rank_any}")


# ---------------- MAIN ----------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--base_dir',
        required=True,
        help="Base directory containing sample folders and metadata files"
    )
    parser.add_argument(
        '--sample',
        required=True,
        help="Sample name (e.g., SRR2329607)"
    )

    args = parser.parse_args()

    base_dir = args.base_dir
    sample_name = args.sample

    sample_dir = os.path.join(base_dir, f"{sample_name}_results")

    if not os.path.isdir(sample_dir):
        raise FileNotFoundError(f"Sample directory not found: {sample_dir}")

    files_needed = {
        "arg_contig": os.path.join(sample_dir, "contig_blast.csv"),
        "arg_plasmid": os.path.join(sample_dir, "plasmid_blast.txt"),
        "arg_phage": os.path.join(sample_dir, "phage_blast.txt"),
        "transposon": os.path.join(sample_dir, "contig_transposon_blast.txt"),
        "arg_transposon": os.path.join(sample_dir, "transposon_blast.txt"),
        "mob": os.path.join(base_dir, f"{sample_name}_mob.csv"),
        "phage": os.path.join(base_dir, f"{sample_name}_PT.tsv"),
        "microbes": os.path.join(base_dir, f"{sample_name}_CAT.txt"),
        "ie": os.path.join(base_dir, f"{sample_name}_IE_hits.tsv"),
        "ie_arg": os.path.join(base_dir, f"{sample_name}_IE_CARD.tsv"),
        "outdir": os.path.join(sample_dir, "RR_output")
    }

    missing = [k for k, v in files_needed.items()
               if k != "outdir" and not os.path.exists(v)]

    if missing:
        raise FileNotFoundError(
            f"[MISSING FILES] {sample_name}: {', '.join(missing)}"
        )

    print(f"[PROCESSING SINGLE SAMPLE] {sample_name}")
    compute_RR(files_needed)
