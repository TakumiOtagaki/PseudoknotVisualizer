"""
PyMOL RC template for PseudoknotVisualizer

Usage:
  Copy this file to your home as ~/.pymolrc_pkv.py, then start PyMOL.

On macOS/Linux:
  cp pymolrc_template.py ~/.pymolrc_pkv.py

If your repository path differs from the default below, edit 'pathtoPKV'.
"""

# ~/.pymolrc_pkv.py
import sys
import pathlib
from pymol import cmd

# --------------------- please modify this line: the path of PseudoknotVisualizer repository -------------------------
pathtoPKV = pathlib.Path.home() / "PseudoknotVisualizer"  # Example: if the repo is under your home directory
# If you installed PseudoknotVisualizer elsewhere, set the path accordingly, e.g.:
# pathtoPKV = pathlib.Path("/Users/ootagakitakumi/PseudoknotVisualizer")
# --------------------------------------------------------------------------------------------------------------------

sys.path.append(str(pathtoPKV))
cmd.run( str(pathtoPKV /  "PseudoknotVisualizer.py"))


