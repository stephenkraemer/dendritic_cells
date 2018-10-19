""" Reproduce DC analysis

Preparation:

1. Install the base environment
conda env create -f dendritic_cells.yml

"""
import subprocess

subprocess.run("""
source activate dendritic_cells
python3 -m dendritic_cells.wgbs.dmr_calling
""")