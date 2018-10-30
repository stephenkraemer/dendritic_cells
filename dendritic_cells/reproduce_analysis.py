""" Reproduce DC analysis

Preparation:
1. Make sure that conda is available
2. Create an environment for running the reproduce_analysis script:
    conda env create -f /home/stephen/projects/dendritic_cells/env.yaml
    pip install -e dendritic_cells
3. Create the environments required for the individual analysis steps
   (yaml files in dendritic_cells/envs)
"""

# DMR calling
# ======================================================================
import subprocess

subprocess.run("""
    source ~/.bashrc
    conda activate dc_dmr_calling
    python3 -m dendritic_cells.wgbs.dmr_calling
""", shell=True, check=True, executable='/bin/bash')

# DMR characterization
# ======================================================================

# Clustering and heatmaps
# ----------------------------------------------------------------------
subprocess.run("""
    source ~/.bashrc
    conda activate dc_dmr_characterization
    python3 -m dendritic_cells.wgbs.clustering_and_heatmaps
""")