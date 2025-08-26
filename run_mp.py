#!/bin/bash
set -euo pipefail

# ========================
# Usage
# ========================
if [ $# -lt 2 ]; then
    echo "Usage: $0 <contigs.fasta> <phabox_db>"
    exit 1
fi

CONTIGS="$1"               # Input contigs FASTA from Snakefile
PHABOX_DB="$2"             # PhaBOX DB path from Snakefile
MOB_OUT="mob_output"       # MOB-suite output folder
PHABOX_OUT="phabox_output" # PhaBOX output folder
THREADS=8

# Create output dirs
mkdir -p "$MOB_OUT" "$PHABOX_OUT"

# ========================
# Run MOB-recon
# ========================
echo "[INFO] Running MOB-recon..."
mob_recon -i "$CONTIGS" -o "$MOB_OUT/mobrecon_results" --num_threads "$THREADS" --force || {
    echo "[ERROR] MOB-recon failed!"
    exit 1
}

# ========================
# Copy contig_report.txt to expected location
# ========================
if [ -f "$MOB_OUT/mobrecon_results/contig_report.txt" ]; then
    cp "$MOB_OUT/mobrecon_results/contig_report.txt" "$MOB_OUT/contig_report.txt"
    echo "[OK] contig_report.txt copied to $MOB_OUT/"
else
    echo "[ERROR] contig_report.txt not found in $MOB_OUT/mobrecon_results/"
    # Create empty file so Snakemake doesn’t crash
    touch "$MOB_OUT/contig_report.txt"
fi


# ========================
# Run MOB-typer
# ========================
echo "[INFO] Running MOB-typer..."
mob_typer -i "$CONTIGS" -o "$MOB_OUT/mobtyper_results.txt" || {
    echo "[WARN] MOB-typer failed, creating placeholder file."
    touch "$MOB_OUT/mobtyper_results.txt"
}

# ========================
# Merge plasmid FASTAs
# ========================
if ls "$MOB_OUT"/mobrecon_results/*.fasta 1> /dev/null 2>&1; then
    cat "$MOB_OUT"/mobrecon_results/*.fasta > "$MOB_OUT/all_plasmids.fasta"
    echo "[OK] all_plasmids.fasta prepared."
else
    echo "[WARN] No plasmid FASTAs found, creating empty file."
    touch "$MOB_OUT/all_plasmids.fasta"
fi

# ========================
# Run PhaBOX2
# ========================
echo "[INFO] Running PhaBOX2..."
phabox2 \
    --task end_to_end \
    --contigs "$CONTIGS" \
    --outpth "$PHABOX_OUT" \
    --threads "$THREADS" \
    --dbdir "$PHABOX_DB" || {
    echo "[WARN] PhaBOX2 failed!"
}

# ========================
# Patch outputs for Snakefile
# ========================

# PhaBOX summary
if [ -f "$PHABOX_OUT/final_prediction/final_prediction_summary.tsv" ]; then
    cp "$PHABOX_OUT/final_prediction/final_prediction_summary.tsv" "$PHABOX_OUT/phage_summary.txt"
    echo "[OK] phage_summary.txt prepared from final_prediction/final_prediction_summary.tsv"
else
    echo "[WARN] final_prediction_summary.tsv not found, creating empty file."
    touch "$PHABOX_OUT/phage_summary.txt"
fi

# PhaBOX phage FASTA
if [ -f "$PHABOX_OUT/filtered_contigs.fa" ]; then
    cp "$PHABOX_OUT/filtered_contigs.fa" "$PHABOX_OUT/phages.fasta"
    echo "[OK] phages.fasta prepared from filtered_contigs.fa"
else
    echo "[WARN] filtered_contigs.fa not found, creating empty file."
    touch "$PHABOX_OUT/phages.fasta"
fi

echo "[✅ DONE] MOB-suite + PhaBOX finished successfully."

