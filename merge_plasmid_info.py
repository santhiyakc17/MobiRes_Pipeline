#!/usr/bin/env python3
import pandas as pd

# ================================
# Input files
# ================================
file1 = "mob_output/mobrecon_results/contig_report.txt"
file2 = "mob_output/mobrecon_results/mobtyper_results.txt"
outfile = "blast_results/mob_out.csv"

# ================================
# Load MOB-suite outputs
# ================================
df1 = pd.read_csv(file1, sep="\t")   # contig report
df2 = pd.read_csv(file2, sep="\t")   # mobtyper results

# Merge on primary_cluster_id
merged_df = pd.merge(df1, df2, on="primary_cluster_id", how="inner")

# Keep only plasmid rows
plasmid_df = merged_df[merged_df["molecule_type"] == "plasmid"]

# Group by plasmid cluster & collapse contigs
grouped = plasmid_df.groupby("primary_cluster_id").agg({
    "contig_id": lambda x: ";".join(sorted(set(x))),
    "predicted_mobility_y": "first",
    "rep_type(s)_y": "first",
    "relaxase_type(s)_y": "first",
    "mpf_type_y": "first",
    "mash_neighbor_identification_y": "first",
    "predicted_host_range_overall_name": "first"
}).reset_index()

# Save output
grouped.to_csv(outfile, index=False)
print(f"✅ Saved plasmid-cluster summary → {outfile} with {len(grouped)} plasmids")

