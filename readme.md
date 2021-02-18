# Introduction

This repository documents the methylome analysis for the paper *Constitutive DNA methylation governs dendritic cell development from myeloid-restricted hematopoietic stem cells in homeostasis and autoimmune disease*.

# Used software packages

The versions of all used software packages are documented as conda environment files, in the envs directory. 

# Usage

This code is meant as a documentation of the performed analysis steps. As such, it contains the complete analysis code, hopefully in a relatively useful organization. However, it was not developed to be run on an external system without modifications. Currently, such an endeavor would be further complicated because we cannot publish the used raw data before the publication is accepted.

Steps to reproduce the analysis on an external system would be:

1. Install the conda environments specified in the envs folder
2. The analysis code is intended to be used as a python package of its own.
   1. Clone the analysis code 
   2. Install it in development mode with `pip -e install /path/to/dendritic_cells`
   3. Install the analysis code from another project, which we use as a dependency here: `pip install git+https://github.com/stephenkraemer/mouse_hema_meth@9003ea`
      - As this project is not yet published, access to the private repository would need to be requested from me.
3. Obtain the methylation calls from GEO (GSE164124)
   - again, access would need to be requested prior to publication of the manuscript.
4. As you go through the notebooks and scripts, adapt hardcoded paths (e.g. to the methylation calls) as necessary.

