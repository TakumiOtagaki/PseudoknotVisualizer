from pathlib import Path

# The variable RNAVIEW is the path to the RNAView directory and RNAVIEW/bin/rnaview is the path to the RNAView executable.
RNAVIEW = Path("/dir/to/RNAView") # <-- please edit this path to your RNAVIEW directory.


# ------------- Do not edit below this line. ------------------
RNAVIEW_PATH = RNAVIEW / "bin"
PseudoKnotVisualizer_DIR = Path(__file__).parent
INTERMEDIATE_DIR = PseudoKnotVisualizer_DIR / "intermediate"