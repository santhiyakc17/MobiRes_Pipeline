#!/bin/bash
set -e
cd "$(dirname "$0")"   # ensure inside MobiRes

# ----------------------------
# 1. Create venv (Python 3.9 required for pandas==1.0.5 + mob-suite)
# ----------------------------
echo "[*] Creating virtual environment at ./venv with Python 3.9"
python3.9 -m venv --without-pip venv
source venv/bin/activate

# ----------------------------
# 2. Install pip manually
# ----------------------------
echo "[*] Installing pip manually"
curl -sS https://bootstrap.pypa.io/pip/3.8/get-pip.py | python

# ----------------------------
# 3. Upgrade pip/setuptools/wheel
# ----------------------------
pip install --upgrade pip setuptools wheel packaging

# ----------------------------
# 4. Core scientific Python packages
# ----------------------------
# First install a working numpy (must be >=1.23 for numexpr, and also prebuilt)
pip install --only-binary=:all: "numpy==1.23.5"

# Now install scipy (compatible with mob-suite: <2.0, >=1.1.0)
pip install --only-binary=:all: "scipy==1.9.3"

# Rest of core deps
pip install "numexpr==2.8.4" "biopython==1.77" "pulp<=2.6" snakemake

pip install "pandas==1.4.4"


# ----------------------------
# 5. MOB-suite + dependencies
# ----------------------------
pip install ete3 pycurl PyQt5 tables
pip install --no-deps mob-suite==3.0.3
MOB_SETUP=$(python -c "import mob_suite, os; print(os.path.join(os.path.dirname(mob_suite.__file__), '..', 'setup.py'))")
if [ -f "$MOB_SETUP" ]; then
    sed -i 's/pandas<=1.0.5/pandas>=1.0.5/' "$MOB_SETUP"
    echo "[✓] Patched mob-suite setup.py to allow modern pandas"
fi

# ----------------------------
# 6. BLAST+
# ----------------------------
echo "[*] Installing BLAST+"
mkdir -p blast
cd blast
wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*-x64-linux.tar.gz
tar -xzf ncbi-blast-*-x64-linux.tar.gz --strip-components=1
rm ncbi-blast-*-x64-linux.tar.gz
cd ..
BLAST_PATH="export PATH=$(pwd)/blast/bin:\$PATH"
if ! grep -Fxq "$BLAST_PATH" venv/bin/activate; then
    echo "$BLAST_PATH" >> venv/bin/activate
fi
export PATH=$(pwd)/blast/bin:$PATH

# ----------------------------
# 7. PhaBOX2 + deps
# ----------------------------
pip uninstall -y phabox2 || true
if [ ! -d "PhaBOX" ]; then
    git clone https://github.com/KennthShang/PhaBOX.git
fi
cd PhaBOX
pip install .
cd ..

# Required deps (version-pinned for stability)
pip install "kcounter==0.3.6" "pyrodigal==0.5.4" "scikit-learn==1.0.2" \
            "networkx==2.5.1" "tqdm==4.66.1" "matplotlib==3.5.3" \
            "seaborn==0.11.2" "datasets==2.14.5" "joblib==1.2.0" \
            "torch==1.13.1" "transformers==4.36.2" "accelerate==0.26.1"

# ----------------------------
# 8. Install Diamond (latest version, into venv/bin)
# ----------------------------
echo "[*] Installing Diamond v0.9.14"
mkdir -p diamond
cd diamond
wget -q http://github.com/bbuchfink/diamond/releases/download/v0.9.14/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz
rm diamond-linux64.tar.gz
cd ..
ln -sf $(pwd)/diamond/diamond venv/bin/diamond
echo "[✓] Diamond installed at $(pwd)/venv/bin/diamond"

# ----------------------------
# 9. Link mcl into venv
# ----------------------------
if [ -x /usr/bin/mcl ]; then
    ln -sf /usr/bin/mcl venv/bin/mcl
    echo "[✓] mcl linked into venv"
else
    echo "[!] mcl not found globally, please install with: sudo apt install mcl"
fi

# ----------------------------
# 10. Patch MOB-suite for Biopython >=1.83
# ----------------------------
MOB_UTILS=$(python -c "import mob_suite, os; print(os.path.join(mob_suite.__path__[0], 'utils.py'))" 2>/dev/null || true)
if [ -f "$MOB_UTILS" ]; then
    sed -i 's/from Bio.SeqUtils import GC/from Bio.SeqUtils import gc_fraction as GC/' "$MOB_UTILS"
    sed -i 's/GC(/gc_fraction(/g' "$MOB_UTILS"
    echo "[✓] MOB-suite patched successfully"
else
    echo "[!] MOB-suite utils.py not found, patch skipped"
fi

# ----------------------------
# 11. Done
# ----------------------------
deactivate
echo ""
echo "[✓] Setup complete!"
echo "To activate, run:"
echo "    source venv/bin/activate"
echo ""
echo "Verify installs with:"
echo "    mob_typer --help"
echo "    blastn -version"
echo "    snakemake --version"
echo "    phabox2 --help"
echo "    mcl --version"
echo "    diamond --version"
