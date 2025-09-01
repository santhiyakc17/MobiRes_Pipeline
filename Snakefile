
# MobiRes Snakefile


# --- Load config.yaml ---
configfile: "config.yaml"

CONTIGS   = config["contigs"]
MICROBES  = config["microbes"]
C1_SCORE  = config["c1_score"]
CARD_DB   = config["card_db"]
TN_DB     = config["transposon_db"]
PHABOX_DB = config["phabox_db"]
THREADS   = config["threads"]

rule all:
    input:
        "output/sample_ERR_summary.csv",
        "output/node_ERR_with_HP.csv",
        "output/top_20_ERR_nodes_with_HP.csv"

# --------------------------------
# Step 1: MOB-suite + PhaBOX
# --------------------------------
rule mob_phabox:
    input:
        contigs = CONTIGS
    output:
        mob = "mob_output/mobtyper_results.txt",
        contig_report = "mob_output/contig_report.txt",
        plasmids = "mob_output/all_plasmids.fasta",
        phage = "phabox_output/phage_summary.txt"
    log:
        "logs/mob_phabox.log"
    shell:
        """
        bash run_mp.py {input.contigs} {PHABOX_DB}
        """

# --------------------------------
# Step 2: BLAST (ARGs + Transposons)
# --------------------------------
rule blast_all:
    input:
        contigs = CONTIGS,
        card_db = [
            CARD_DB + ".nin",
            CARD_DB + ".nsq",
            CARD_DB + ".nhr"
        ],
        tn_db = [
            TN_DB + ".nin",
            TN_DB + ".nsq",
            TN_DB + ".nhr"
        ]
    output:
        contig_blast = "blast_results/ARG_BLAST_contig.txt",
        plasmid_blast = "blast_results/ARG_BLAST_plasmid.txt",
        phage_blast = "blast_results/ARG_BLAST_phage.txt",
        tn_blast = "blast_results/TN_BLAST_raw.txt",
        contig_fa = "blast_results/ARG_sequences_contig.fasta",
        plasmid_fa = "blast_results/ARG_sequences_plasmids.fasta",
        phage_fa = "blast_results/ARG_sequences_phage.fasta",
        tn_fa = "blast_results/ARG_sequences_transposon.fasta"
    params:
        card_prefix = CARD_DB,
        tn_prefix = TN_DB
    log:
        "logs/blast_all.log"
    shell:
        """
        python3 blast_all.py \
            --contigs {input.contigs} \
            --card_db {params.card_prefix} \
            --tn_db {params.tn_prefix} \
            --outdir blast_results
        """

# --------------------------------
# Step 3: Merge plasmid info
# --------------------------------
rule merge_plasmid_info:
    input:
        contig_report = "mob_output/contig_report.txt",
        mobtyper = "mob_output/mobtyper_results.txt"
    output:
        "blast_results/mob_out.csv"
    log:
        "logs/merge_plasmid_info.log"
    shell:
        """
        python3 merge_plasmid_info.py \
            --contig_report {input.contig_report} \
            --mobtyper {input.mobtyper} \
            --out {output}
        """

# --------------------------------
# Step 4: Merge ARG BLAST with C1 Score
# --------------------------------
rule merge_ARG_C1:
    input:
        arg_blast = "blast_results/ARG_BLAST_contig.txt",
        c1 = C1_SCORE
    output:
        "blast_results/ARG_BLAST.csv"
    log:
        "logs/merge_ARG_C1.log"
    shell:
        """
        python3 score.py \
            --arg_blast {input.arg_blast} \
            --c1_score {input.c1} \
            --out {output}
        """
        
# --------------------------------
# Step 5: Compute ERR + ResCon
# --------------------------------
rule compute_ERR:
    input:
        arg_contig = "blast_results/ARG_BLAST.csv",
        arg_plasmid = "blast_results/ARG_sequences_plasmids.fasta",
        arg_phage = "blast_results/ARG_sequences_phage.fasta",
        transposon = "blast_results/TN_BLAST_raw.txt",    # fixed
        arg_transposon = "blast_results/ARG_sequences_transposon.fasta",
        mob = "blast_results/mob_out.csv",
        phage = "phabox_output/phage_summary.txt",
        microbes = MICROBES
    output:
        "output/sample_ERR_summary.csv",
        "output/node_ERR_with_HP.csv",
        "output/top_20_ERR_nodes_with_HP.csv"
    log:
        "logs/compute_ERR.log"
    shell:
        """
        python3 m3.py \
            --arg_contig {input.arg_contig} \
            --arg_plasmid {input.arg_plasmid} \
            --arg_phage {input.arg_phage} \
            --transposon {input.transposon} \
            --arg_transposon {input.arg_transposon} \
            --mob {input.mob} \
            --phage {input.phage} \
            --microbes {input.microbes} \
            --outdir output &> {log}
        """

