#!/usr/bin/env python3
import os
import subprocess
import argparse
from Bio import SeqIO

# Utility Function

def run_blast(query, db, output, evalue="1e-5"):
    """Run BLAST search and save tabular output."""
    os.makedirs(os.path.dirname(output), exist_ok=True)
    cmd = [
        "blastn",
        "-query", query,
        "-db", db,
        "-out", output,
        "-evalue", evalue,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    ]
    print("[INFO] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def parse_blast_hits(blast_file, identity_threshold=70.0):
    """Parse BLAST tabular results and filter hits by identity %."""
    hits = set()
    if not os.path.exists(blast_file):
        return hits
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                qseqid, pident = parts[0], float(parts[2])
                if pident >= identity_threshold:
                    hits.add(qseqid)
    return hits

def extract_sequences(hits, fasta, output, seq_type=""):
    """Extract sequences corresponding to BLAST hits into FASTA."""
    seqs = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    count = 0
    with open(output, "w") as out:
        for h in hits:
            if h in seqs:
                SeqIO.write(seqs[h], out, "fasta")
                count += 1
    print(f"✅ Extracted {count} {seq_type} sequences → {output}")


# Main

def main():
    parser = argparse.ArgumentParser(description="BLAST ARGs and Transposons")
    parser.add_argument("--contigs", required=True, help="Input contigs FASTA")
    parser.add_argument("--card_db", required=True, help="CARD DB prefix (without .nin/.nsq)")
    parser.add_argument("--tn_db", required=True, help="Transposon DB prefix")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Known inputs
    contig_fasta   = args.contigs
    plasmid_fasta  = "mob_output/all_plasmids.fasta"
    phage_fasta    = "phabox_output/phages.fasta"
    arg_db         = args.card_db
    tn_db          = args.tn_db
    identity_thr   = 90.0

    # Outputs
    contig_blast = os.path.join(args.outdir, "ARG_BLAST_contig.txt")
    plasmid_blast = os.path.join(args.outdir, "ARG_BLAST_plasmid.txt")
    phage_blast = os.path.join(args.outdir, "ARG_BLAST_phage.txt")
    tn_blast = os.path.join(args.outdir, "TN_BLAST_raw.txt")

    contig_out = os.path.join(args.outdir, "ARG_sequences_contig.fasta")
    plasmid_out = os.path.join(args.outdir, "ARG_sequences_plasmids.fasta")
    phage_out = os.path.join(args.outdir, "ARG_sequences_phage.fasta")
    tn_out = os.path.join(args.outdir, "ARG_sequences_transposon.fasta")

    # ---- ARG detection ----
    run_blast(contig_fasta, arg_db, contig_blast)
    extract_sequences(parse_blast_hits(contig_blast, identity_thr),
                      contig_fasta, contig_out, seq_type="ARG in contigs")

    if os.path.exists(plasmid_fasta):
        run_blast(plasmid_fasta, arg_db, plasmid_blast)
        extract_sequences(parse_blast_hits(plasmid_blast, identity_thr),
                          plasmid_fasta, plasmid_out, seq_type="ARG in plasmids")

    if os.path.exists(phage_fasta):
        run_blast(phage_fasta, arg_db, phage_blast)
        extract_sequences(parse_blast_hits(phage_blast, identity_thr),
                          phage_fasta, phage_out, seq_type="ARG in phages")

    # ---- Transposon detection ----
    run_blast(contig_fasta, tn_db, tn_blast)
    extract_sequences(parse_blast_hits(tn_blast, identity_thr),
                      contig_fasta, tn_out, seq_type="ARG in transposons")
    

    print("🎯 BLAST analysis complete.")

if __name__ == "__main__":
    main()

