# PseudoknotVisualizer
PseudoknotVisualizer is a Pymol Extension for visualization that assigns different colors to each Pseudoknot layer.

# Overview
This tool enables us to visually understand the RNA tertiary structures with pseudoknots.
This is essential for prediction of tertiary structures and selecting the best structure from the structure ensemble.

PseudoknotVisualizer is available at pymol, meaning that it is easy to install.
However, the installation of dependent module "RNAVIEW" is a bit complicated.
RNAVIEW cannot work in Macbook.

This tool has two modes of use: CLI and GUI (using PyMOL).

# Demo
![pymol_demo_6T3R](https://github.com/user-attachments/assets/444dadcf-72f9-46cf-b442-484c611adcfd)


# How to Install
The dependent module "RNAVIEW"  is only available in Linux OS.
## Implementation on WSL & conda Recommended
```
wsl$ conda create -n pymol python=3.11.0
wsl$ conda activate pymol
wsl$ conda install -c conda-forge pymol-open-source
wsl$ pymol
```
Type `pymol` in conda pymol env, then open source pymol app starts.

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
1. Clone this repository.
2. Rewrite config.py; path and other enviromental variables.
3. Rewrite or create `~/.pymolrc.py` in order to load the extension at startup.

-----

1. Cloning
```sh
# 1. 
$ git clone git@github.com:TakumiOtagaki/PseudoknotVisualizer.git
```

2. Rewriting config.py
Add the two variables related to the RNAVIEW you installed earlier, RNAVIEW and PATH of RNAVIEW, to config.py.
```config.py
# RNAVIEW related variables
RNAVIEW="/path/to/RNAVIEW"
RNAVIEW_PATH="/path/to/RNAVIEW/bin"
# The path to this PseudoknotVisualizer repository.
PseudoKnotVisualizer_DIR = "/path/to/PseudoknotVisualizer"
...
```

3. Rewrite or create `~/.pymolrc.py`
This is optional, but if you want to automatically load the extension at startup, please followth einstructions below.

```sh
$ vim ~/.pymolrc.py
```
And write a few lines as follows.
```~/.pymolrc.py
# ~/.pymolrc.py
import sys
sys.path.insert(0, "/path/to/PseudoknotVisualizer")

from pymol import cmd
cmd.run("/path/to/PseudoknotVisualizer/PseudoknotVisualizer.py")
```


Now, you can use our extension easily.

# How to use
## Basic Usage
After the installation, Pymol automatically loads our extensions.
Before calling our extention, you should load the model. 
After loading models, it can be called and used as follows:
```
# pymol command line after loading model
pkv $pdb_object $chainID
```
 - pdb_object = a model, it can be multimer.
 - chainID = A, B, C, ...

For example, if you want to visualize the pseudoknots in 1kpd in PDB, run the followings:
```
# loading
fetch 1kpd
# at pymol command line.
pkv 1kpd, A
# OR
pseudoknotvisualizer 1kpd, A
```


## Changing Colors (Optinal)
If you want to change the color of each layer, modify PseudoknotVisualizer/colors.json. You can also add new lines.

Make sure to update colors.json before launching PyMOL.

```colors.json
{
    "1": "red",
    "2": "blue",
    "3": "green",
    "4": "yellow",
    "5": "purple",
    "6": "orange",
    "default": "gray"
}
```
If the number of layers is greater than 6, PseudoknotVisualizer will color the 7th and subsequent layers with a default color.

To increase this limit beyond 6, simply add entries like "7": "another color".

# For CLI user
After the installation (except for step 4), you can use our CLI.

## CLI Usage
```sh
$ python '/large/otgk/PseudoknotVisualizer/CLI_PseudoknotVisualizer.py' --help

usage: CLI_PseudoknotVisualizer.py [-h] -i INPUT -o OUTPUT -f {chimera,pymol} [-m MODEL] [-c CHAIN]

Visualize pseudoknots in RNA structure

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file containing RNA structure
  -o OUTPUT, --output OUTPUT
                        Output script file for visualization
  -f {chimera,pymol}, --format {chimera,pymol}
                        Format of RNA structure (chimera or pymol)
  -c CHAIN, --chain CHAIN
                        Chain ID for RNA structure, default is A

chimera options:
  Options specific to Chimera format

  -m MODEL, --model MODEL
                        Model ID (required if Chimera format is selected)
```
## Installation for CLI user
 1.	Run
    ```sh
   	conda create -n pymol python=3.11
    conda activate pymol
    ```
2.	Run
  ```
 pip install -r requirements.txt
 ```
3.	Complete steps 1 through 3 from the installation instructions above.

4.	Execute the command `python CLI_PseudoknotVisualizer.py -i input.pdb -o ...`

## Example of CLI usage
```sh
python PseudoknotVisualizer/CLI_PseudoknotVisualizer.py \
  -i test/1KPD.pdb \  # input pdb file.
  -o test/coloring_1KPD.0.A.pymol.txt \ # path of output script txtfile
  -c A \ # chain ID
  -f pymol \ # format: chimera or pymol
  -m 0 # model ID in your viewer.

# to do list 
- Now intermediate dir should end with "/", it could cause an error.
- We are going to support "non-canonical bp" mode, extension for chimera.


