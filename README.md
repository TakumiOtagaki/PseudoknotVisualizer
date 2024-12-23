# PseudoknotVisualizer
PseudoknotVisualizer is a **PyMOL Extension** for visualization that assigns **different colors to each Pseudoknot layer**.
<img src="https://github.com/TakumiOtagaki/PseudoknotVisualizer/blob/main/demo.png" alt="pymol_demo_6T3R" width="100%">

# Overview
This tool enables us to visually understand the RNA tertiary structures with pseudoknots.
This is essential for prediction of tertiary structures and selecting the best structure from the structure ensemble.

PseudoknotVisualizer is available at PyMOL, meaning that it is easy to install.
This tool has two modes of use: CLI and GUI (using PyMOL).

# Demo

<img src="https://github.com/TakumiOtagaki/PseudoknotVisualizer/blob/main/uncolored_6t3r.png" alt="pymol_demo_6T3R" width="30%"><img src="https://github.com/TakumiOtagaki/PseudoknotVisualizer/blob/main/colored_6t3r.png" alt="pymol_demo_6T3R" width="30%"><img src="https://github.com/TakumiOtagaki/PseudoknotVisualizer/blob/main/colored_6t3r.gif" alt="demo_gif" width="35.5%">

- Left: Before coloring pseudoknots.
- Right: After coloring
 - red: pseudoknot layer 1
 - blue: pseudoknot layer 2
 - green: pseudoknot layer 3
 - (gray: Main Layer)








# How to Install

## Prerequisite: Installtion of RNAView
You need to pre-install [rnaview](https://github.com/rcsb/RNAView).
The installation steps are like,
```
git clone git@github.com:rcsb/RNAView.git
cd RNAView
make
ls bin 
```

After the installation of RNAView, your `~/.bashrc` should contain these two lines.
```~/.bashrc
# ------------ RNAView setting ---------------
export PATH=$PATH:/path/to/RNAView/bin
export RNAVIEW=/path/to/RNAView/
```


## Prepairing "pymol" Conda Environment. (Recommended)
```
conda create -n pymol python=3.11.0
conda activate pymol
conda install pandas numpy
conda install -c conda-forge pymol-open-source
pymol
```
Type `pymol` in conda pymol env, then open source PyMOL app will start.



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
Add the two variables related to the RNAView you installed earlier, RNAVIEW and PATH of RNAVIEW, to config.py.
Please rewrite the three lines:
```config.py
# RNAVIEW related variables
RNAVIEW="/path/to/RNAVIEW"
RNAVIEW_PATH="/path/to/RNAVIEW/bin"
# The path to this PseudoknotVisualizer repository.
PseudoKnotVisualizer_DIR = "/path/to/PseudoknotVisualizer"
...
```

3. Rewrite or create `~/.pymolrc.py`
This is optional, but if you want to automatically load the extension at startup, please follow the instructions below.

```sh
$ vim ~/.pymolrc.py
```
And write a few lines as follows.
```~/.pymolrc.py
# ~/.pymolrc.py
import sys
pathtoPseudoknotVisualizer = "/path/to/PseudoknotVisualizer" # <-- Please modify this line!
sys.path.insert(0, pathtoPseudoknotVisualizer)

from pymol import cmd
cmd.run(pathtoPseudoknotVisualizer + "PseudoknotVisualizer.py")
```


Now, you can use our extension easily.
After this step, the PseudoknotVisualizer extention will be automatically loaded when PyMOL starts.

If you skipped this step, you have to load manually (not recommended).

# How to use
## Basic Usage
After loading models, it can be called and used as follows:
```
# PyMOL command line after loading model
pkv $pdb_object, $chainID
```
 - pdb_object = a model, it can be multimer.
 - chainID = A, B, C, ...

For example, if you want to visualize the pseudoknots in 1kpd in PDB, run the followings:
```
# loading
fetch 1kpd
# at PyMOL command line.
pkv 1kpd, A
# OR
pkv sele, A # if 1kpd is selected.
```
As you can see from this example, you can use "sele" to identify the model.

<img src="https://github.com/TakumiOtagaki/PseudoknotVisualizer/blob/main/casp15_examples.png" alt="pymol_demo_6T3R" width="50%">


## Changing Colors (Optinal)
If you want to change the color of each layer, modify PseudoknotVisualizer/colors.json. You can also add new lines.

Make sure to update colors.json before launching PyMOL.

```colors.json
{
    "1": "gray",
    "2": "red",
    "3": "blue",
    "4": "green",
    "5": "yellow",
    "6": "purple",
    "default": "gray"
}
```
If the number of layers is greater than 6, PseudoknotVisualizer will color the 7th and subsequent layers with a default color.

To increase this limit beyond 6, simply add entries like "7": "another color".

# For CLI user
After the installation (except for step 4), you can use our CLI.

## CLI Usage
```sh
$ python 'PseudoknotVisualizer/CLI_PseudoknotVisualizer.py' --help

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

Also you can fetch PDB file from Protein Data Bank.
```sh
$ python fetch_pdb.py
Enter PDB ID (q to quit): 1kpd
Enter output filename(if not provided, pdb_id.pdb will be created in current directory): 
PDB file for 1kpd downloaded as ./1kpd.pdb
```
Then, 1kpd.pdb is downloaded in current directory.

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
```


# PyMOL errors
## Case 1
PseudoknotVisualizer can not color accurately the specified molecule when the sequence index in PyMOL viewer does not start with 1.
If so, please check the sequence index (`your_start_index`) pushing "S" button and do as followings:
```sh
select rna_chain, your_pdb_id and chain your_chain_id
alter rna_chain, resi = int(resi) - (your_start_index - 1)
pkv your_pdb_id, A
```
Here, please rewrite 
 - your_pdb_id
 - your_chain_id
 - your_start_index.

Then, it will work.



# to do list
- We are going to support "non-canonical bp" mode, extension for chimera.


