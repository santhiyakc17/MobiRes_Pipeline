🧬
MobiRes is a computational pipeline designed to evaluate the resistome risk of microbial communities by integrating antibiotic resistance gene (ARG) profiles with the mobilome (mobile genetic elements, MGEs).

This pipeline integrates **MOB-suite**, **PhaBOX2**, **BLAST**, **DIAMOND**, and custom scripts to analyze ARG mobility and compute resistome risk. It is implemented with **Snakemake** for reproducibility.

---

⚙️ Requirements

* Linux / WSL
* Python >= 3.8 (with Biopython, Pandas, Scikit-learn, etc.)
* Snakemake >= 7
* BLAST+ (`makeblastdb`, `blastn`)
* DIAMOND (Buchfink et al., 2015, *Nature Methods*)
* MOB-suite >= 3.0.3
* PhaBOX2 (requires `prodigal-gv`, `diamond`, `mcl`)
* Conda/venv for reproducibility

---

📥 Downloading the PhaBOX database

Download using `wget` and place in the database folder for smooth activation:

```bash
wget https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2_1.zip
unzip phabox_db_v2_1.zip > /dev/null
```

---

🚀 Installation & Run

```bash
# Clone the repo
git clone https://github.com/santhiyakc17/MobiRes_Pipeline.git
cd MobiRes_Pipeline

# Create and activate virtual environment
conda create -n venv python=3.9 -y
conda activate venv

# Install dependencies
./setup_env.py

# Activate environment and run pipeline
source venv/bin/activate
snakemake -s Snakefile --cores 4
```

---


📄 References

MOB-suite: https://github.com/phac-nml/mob-suite

PhaBOX: https://github.com/KennthShang/PhaBOX

BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi

Diamond: https://github.com/bbuchfink/diamond


