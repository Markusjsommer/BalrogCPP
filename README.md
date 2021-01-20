# Bacterial Annotation by Learned Representation Of Genes (C++ version)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/installer/conda.svg)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/balrog.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/platforms.svg)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/license.svg)](https://anaconda.org/bioconda/balrog)

## Overview
Balrog is a prokaryotic gene finder based on a Temporal Convolutional Network. We took a data-driven approach to prokaryotic gene finding, relying on the large and diverse collection of already-sequenced genomes. By training a single, universal model of bacterial genes on protein sequences from many different species, we were able to match the sensitivity of current gene finders while reducing the overall number of gene predictions. Balrog does not need to be refit on any new genome.

Preprint available on bioRxiv [here](https://www.biorxiv.org/content/10.1101/2020.09.06.285304v1).

## Getting started
    # install via conda
    conda create -n balrog_env python=3.7
    conda activate balrog_env
    
    conda install pytorch=1.7.1 -c conda-forge
    conda install balrog -c conda-forge -c bioconda
    
    balrog --help

Conda can be very slow and are working on providing precompiled binaries for unix systems in the near future

