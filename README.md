# Bacterial Annotation by Learned Representation Of Genes (C++ version)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/installer/conda.svg)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/downloads.svg)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/platforms.svg)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/license.svg)](https://anaconda.org/bioconda/balrog)

## Overview
This repo is a work in progress...

Balrog is a prokaryotic gene finder based on a Temporal Convolutional Network. We took a data-driven approach to prokaryotic gene finding, relying on the large and diverse collection of already-sequenced genomes. By training a single, universal model of bacterial genes on protein sequences from many different species, we were able to match the sensitivity of current gene finders while reducing the overall number of gene predictions. Balrog does not need to be refit on any new genome.

Preprint available on bioRxiv [here](https://www.biorxiv.org/content/10.1101/2020.09.06.285304v1).

## Getting started
    # install via conda
    conda install balrog --channel bioconda
    
NOTE: the pretrained gene and TIS models can be found in the models directory in this repository (we will add these paths by default in the next release)

Conda can be very slow and are working on providing precompiled binaries for unix systems in the near future

We are also working on integrating MMseqs2 as described in the preprint to reduce false positive predictions
