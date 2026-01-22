#!/usr/bin/env python3
import pandas as pd


# Input/output paths 

blast_file = "blast_results/ARG_BLAST_contig.txt"
score_file = "input/CI_score.csv"
output_file = "blast_results/ARG_BLAST.csv"

# BLAST column names
cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


# Load BLAST results

blast_df = pd.read_csv(blast_file, sep="\t", names=cols)

# Extract subtype from sseqid
def extract_arg_subtype(value):
    if pd.isna(value):
        return None
    if "~~~" in value:
        return value.split("~~~")[1].strip()
    return value.strip()

blast_df["Subtype"] = blast_df["sseqid"].apply(extract_arg_subtype)


# Load CI score mapping

score_df = pd.read_csv(score_file, sep=None, engine="python")
score_df.columns = ["Class", "Subtype", "CI Score"]

# Normalize subtype names for safe merging
blast_df["Subtype_norm"] = blast_df["Subtype"].str.strip().str.lower()
score_df["Subtype_norm"] = score_df["Subtype"].str.strip().str.lower()


# Merge BLAST + CI scores

merged = pd.merge(
    blast_df,
    score_df[["Subtype_norm", "CI Score"]],
    on="Subtype_norm",
    how="left"
)

# Drop helper column
merged = merged.drop(columns=["Subtype_norm"])

# Replace NaN with 0 for unmatched ARGs
merged["CI Score"] = merged["CI Score"].fillna(0.0)

# Drop duplicates
merged = merged.drop_duplicates(subset=["qseqid", "Subtype"])


# Save result

merged.to_csv(output_file, index=False)
print(f" Saved: {output_file} with {merged['CI Score'].notna().sum()} ARGs matched to CI scores")

