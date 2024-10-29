# PseudoknotVisualizer
PseudoknotVisualizer is a Pymol Extension for visualization that assigns different colors to each Pseudoknot layer.

# Overview
This tool enables us to visually understand the RNA tertiary structures with pseudoknots.
This is essential for prediction of tertiary structures and selecting the best structure from the structure ensemble.

PseudoknotVisualizer is available at pymol, meaning that it is very easy to install.

# Demo
![pymol_demo_6T3R](https://github.com/user-attachments/assets/444dadcf-72f9-46cf-b442-484c611adcfd)


# How to Install

## Prerequisite: Installtion of RNAVIEW
You need to pre-install [rnaview](https://github.com/rcsb/RNAView).

After the installation of rnaview, your `~/.bashrc` should contain these two lines.
```~/.bashrc
# rnaview setting
export PATH=$PATH:/path/to/rnaview/bin
export RNAVIEW=/path/to/rnaview/
```

## Installation of PseudoknotVisualizer
### overview
1. clone this repository.
2. Rewrite config.py


# to do list
- Now intermediate dir should end with "/", it could cause an error.
- We are going to support "non-canonical bp" mode, extension for chimera.


