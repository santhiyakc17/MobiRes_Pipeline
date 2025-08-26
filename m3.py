
import argparse
import pandas as pd
from collections import defaultdict
import re

import os
import pandas as pd
from Bio import SeqIO


# --- Weights ---
PLASMID_WEIGHTS = {'conjugative': 1.0, 'mobilizable': 0.75, 'non-mobilizable': 0.5}
PHAGE_WEIGHTS = {'virulent': 1.0, 'temperate': 0.5}
TRANSPOSON_WEIGHTS = {'conjugative': 1.0, 'composite': 0.75, 'unit': 0.5}
DEFAULT_WEIGHT = 0.1

PATHOGENS = {
    'escherichia coli', 'salmonella enterica', 'klebsiella pneumoniae',
    'staphylococcus aureus', 'acinetobacter baumannii',
    'pseudomonas aeruginosa', 'enterococcus faecium'
}

def normalize_id(x):
    return str(x).split(':')[-1].split('_length_')[0].lower().strip()

def normalize_node_id(node_id):
    match = re.match(r'(NODE_\d+)', node_id, re.IGNORECASE)
    return match.group(1).lower() if match else node_id.lower().strip()

from Bio import SeqIO

def load_arg(file):
    if file.endswith((".fasta", ".fa")):
        # Handle FASTA files
        records = list(SeqIO.parse(file, "fasta"))
        df = pd.DataFrame({"qseqid": [r.id for r in records]})
    elif file.endswith(".txt"):
        # BLAST tabular output (no header)
        cols = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        df = pd.read_csv(file, sep="\t", names=cols)
    else:
        # For CSV or TSV with headers
        df = pd.read_csv(file, sep=None, engine="python")
        df.columns = df.columns.str.strip()
        if "qseqid" not in df.columns:
            raise ValueError(f"qseqid column missing in {file}, found {df.columns.tolist()}")

    # Normalize and standardize
    if "C1 Score" in df.columns:
        df.rename(columns={"C1 Score": "C1_score"}, inplace=True)

    df["node"] = df["qseqid"].apply(normalize_id)

    if "C1_score" not in df.columns:
        df["C1_score"] = 0.0

    return df


def compute_ARG_score(arg_df):
    total_args = len(arg_df)
    subtype_counts = arg_df.groupby(['node', 'sseqid']).size().reset_index(name='count')
    subtype_c1 = arg_df.groupby('sseqid')['C1_score'].max()
    subtype_counts['C1_score'] = subtype_counts['sseqid'].map(subtype_c1)
    subtype_counts['score'] = subtype_counts['count'] * subtype_counts['C1_score']
    per_node_scores = subtype_counts.groupby('node')['score'].sum()
    norm_scores = {node: (score / total_args) if total_args else 0.0 for node, score in per_node_scores.items()}
    return per_node_scores.to_dict(), norm_scores, total_args

def load_microbes(file):
    df = pd.read_csv(file, sep='\t')
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_')
    return {normalize_node_id(str(row['node'])): str(row['species']).strip().lower() for _, row in df.iterrows()}

def compute_host_pathogenicity(arg_df, species_map, total_species):
    hp_scores = defaultdict(float)
    for _, row in arg_df.iterrows():
        node = normalize_node_id(str(row['node']))
        species = species_map.get(node, '')
        arg_weight = 1.0 if species in PATHOGENS else 0.1
        hp_scores[node] += arg_weight
    for node in hp_scores:
        hp_scores[node] /= total_species if total_species else 1
    return hp_scores

def load_transposons(file):
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
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
    return {contig: sum(weights) / len(weights) for contig, weights in weight_map.items()}

from Bio import SeqIO

def load_ARGs_in_transposon_nodes(file):
    # Case 1: BLAST tabular result
    if file.endswith(".txt"):
        df = pd.read_csv(file, sep="\t", header=None)
        # assign column names if needed
        df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore"][:len(df.columns)]
    # Case 2: FASTA sequences
    elif file.endswith(".fasta") or file.endswith(".fa"):
        records = list(SeqIO.parse(file, "fasta"))
        df = pd.DataFrame({"qseqid": [r.id for r in records]})
    else:
        raise ValueError(f"Unsupported transposon file format: {file}")

    # Normalize IDs
    if "qseqid" not in df.columns:
        raise ValueError(f"'qseqid' column missing in {file}, found {df.columns.tolist()}")
    df["node"] = df["qseqid"].apply(normalize_id)
    return df

def load_plasmids(file):
    df = pd.read_csv(file)
    weights = defaultdict(lambda: DEFAULT_WEIGHT)
    for _, row in df.iterrows():
        sid = normalize_id(row.get('sample_id', row.get('contig_id', 'unknown')))
        mobility = str(row.get('predicted_mobility', '')).lower()
        weights[sid] = PLASMID_WEIGHTS.get(mobility, DEFAULT_WEIGHT)
    return weights

def load_phages(file):
    df = pd.read_csv(file, sep='\t')
    weights = defaultdict(lambda: DEFAULT_WEIGHT)
    for _, row in df.iterrows():
        acc = normalize_id(row['Accession'])
        typ = str(row.get('PhaTYP', '')).lower()
        weights[acc] = PHAGE_WEIGHTS.get(typ, DEFAULT_WEIGHT)
    return weights

def assign_rank(row):
    mge_count = int(row['Plasmid_weight'] > DEFAULT_WEIGHT) +                 int(row['Phage_weight'] > DEFAULT_WEIGHT) +                 int(row['Transposon_weight'] > DEFAULT_WEIGHT)
    is_pathogen = row['HP'] > DEFAULT_WEIGHT
    if row['ARGs'] > 0 and mge_count >= 2 and is_pathogen:
        return 1
    elif row['ARGs'] > 0 and mge_count >= 1 and is_pathogen:
        return 2
    elif row['ARGs'] > 0 and mge_count >= 1:
        return 3
    elif row['ARGs'] > 0:
        return 4
    return 5

def compute_ERR(args):
    df_arg = load_arg(args.arg_contig)
    plasmid_df = load_arg(args.arg_plasmid)
    phage_df = load_arg(args.arg_phage)

    species_map = load_microbes(args.microbes)
    total_species = len(set(species_map.values()))
    hp_dict = compute_host_pathogenicity(df_arg, species_map, total_species)

    raw_arg_scores, norm_arg_scores, total_args = compute_ARG_score(df_arg)
    plasmid_weights = load_plasmids(args.mob)
    phage_weights = load_phages(args.phage)
    transposon_weights = load_transposons(args.transposon)
    args_in_transposon = load_ARGs_in_transposon_nodes(args.arg_transposon)

    all_nodes = set(df_arg['node']) | set(plasmid_weights) | set(phage_weights) | set(transposon_weights)

    rows = []
    for node in sorted(all_nodes):
        arg_count = df_arg[df_arg['node'] == node].shape[0]
        arg_score = norm_arg_scores.get(node, 0.0)
        hp = hp_dict.get(node, 0.0)
        plasmid_arg_count = plasmid_df[plasmid_df['node'] == node].shape[0]
        phage_arg_count = phage_df[phage_df['node'] == node].shape[0]
        plasmid_type_score = plasmid_weights.get(node, DEFAULT_WEIGHT)
        phage_type_score = phage_weights.get(node, DEFAULT_WEIGHT)
        trans_type_score = transposon_weights.get(node, DEFAULT_WEIGHT)

        total_plasmid_args = len(plasmid_df)
        total_phage_args = len(phage_df)
        plasmid_weight = (plasmid_arg_count * plasmid_type_score) / total_plasmid_args if total_plasmid_args else DEFAULT_WEIGHT
        phage_weight = (phage_arg_count * phage_type_score) / total_phage_args if total_phage_args else DEFAULT_WEIGHT
        trans_weight = trans_type_score

        mobility_index = (plasmid_weight + phage_weight + trans_weight) / 3
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
            'ARGs_in_transposon': args_in_transposon.get(node, 0),
            'ERR_Score': err_score
        })

    df = pd.DataFrame(rows)
    df['ResCon_Rank'] = df.apply(assign_rank, axis=1)

    df['has_ARG'] = df['ARGs'] > 0
    df['has_plasmid'] = df['Plasmid_weight'] > DEFAULT_WEIGHT
    df['has_phage'] = df['Phage_weight'] > DEFAULT_WEIGHT
    df['has_transposon'] = df['Transposon_weight'] > DEFAULT_WEIGHT
    df['has_pathogen'] = df['HP'] > DEFAULT_WEIGHT

   # Ensure output dir exists
    os.makedirs(args.outdir, exist_ok=True)

   # --- Save subsets (5 files) ---
    df[df['has_ARG'] & (df['has_plasmid'] | df['has_phage'] | df['has_transposon'] | df['has_pathogen'])] \
        .to_csv(os.path.join(args.outdir, "nodes_ARG_plus_any_MGE_or_pathogen12.csv"), index=False)

    df[df['has_ARG'] & df['has_plasmid'] & df['has_phage'] & df['has_transposon'] & df['has_pathogen']] \
        .to_csv(os.path.join(args.outdir, "nodes_ARG_plus_all_MGEs_and_pathogen12.csv"), index=False)

    df.to_csv(os.path.join(args.outdir, 'node_ERR_with_HP.csv'), index=False)

    df.sort_values('ERR_Score', ascending=False).head(20) \
        .to_csv(os.path.join(args.outdir, 'top_20_ERR_nodes_with_HP.csv'), index=False)

    summary_df = pd.DataFrame([{
        'Total_ARGs': df['ARGs'].sum(),
        'Total_ERR_Score': df['ERR_Score'].sum(),
        'Sample_ERR_Score': df['ARGs'].sum() * df['ERR_Score'].sum(),
        'Sample_ResCon_Rank_ARG_plus_any_MGE_or_pathogen': int(round(df[df['has_ARG'] & 
            (df['has_plasmid'] | df['has_phage'] | df['has_transposon'] | df['has_pathogen'])]['ResCon_Rank'].mean())) 
            if not df[df['has_ARG'] & (df['has_plasmid'] | df['has_phage'] | df['has_transposon'] | df['has_pathogen'])].empty else 0
    }])
    summary_df.to_csv(os.path.join(args.outdir, 'sample_ERR_summary.csv'), index=False)

   # --- Print Summary ---
    print(f"\n--- Sample Level Summary ---")
    print(f"Total ARGs: {df['ARGs'].sum()}")
    print(f"Total ERR Score: {df['ERR_Score'].sum():.6f}")
    print(f"Sample ERR Score: {(df['ARGs'].sum() * df['ERR_Score'].sum()):.4f}")

    if not df[df['has_ARG'] & (df['has_plasmid'] | df['has_phage'] | df['has_transposon'] | df['has_pathogen'])].empty:
        sample_rank_any = int(round(df[df['has_ARG'] & 
                                       (df['has_plasmid'] | df['has_phage'] | 
                                        df['has_transposon'] | df['has_pathogen'])]['ResCon_Rank'].mean()))
    else:
        sample_rank_any = 0

    print(f"Sample ResCon Rank (ARG + any MGE or pathogen): {sample_rank_any}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--arg_contig', required=True)
    parser.add_argument('--arg_plasmid', required=True)
    parser.add_argument('--arg_phage', required=True)
    parser.add_argument('--transposon', required=True)
    parser.add_argument('--arg_transposon', required=True)
    parser.add_argument('--mob', required=True)
    parser.add_argument('--phage', required=True)
    parser.add_argument('--microbes', required=True)
    parser.add_argument('--outdir', required=True, help="Output directory for results")
    args = parser.parse_args()
    compute_ERR(args)

