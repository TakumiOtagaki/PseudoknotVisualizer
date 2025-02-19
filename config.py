from pathlib import Path

"""
The variable RNAVIEW is the path to the RNAView directory and RNAVIEW/bin/rnaview is the path to the RNAView executable.
"""

# RNAVIEW related paths
RNAVIEW = Path("/path/to/RNAView") # <-- please edit this path to your RNAVIEW path
RNAVIEW_PATH = RNAVIEW / "bin"

# PseudoknotVisualizer repository path
PseudoKnotVisualizer_DIR = Path("/path/to/PseudoknotVisualizer") # <-- please edit this path to your PseudoKnotVisualizer repository path

# ------------- Do not edit below this line. ------------------
INTEREMEDIATE_DIR = PseudoKnotVisualizer_DIR / "intermediate"