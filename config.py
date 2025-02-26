from pathlib import Path



# ------------------------ plsease edit this path to your RNAVIEW directory. ------------------------

# The variable RNAVIEW is the path to the RNAView directory and RNAVIEW/bin/rnaview is the path to the RNAView executable.
RNAVIEW = Path("/Users/ootagakitakumi/Applications/RNAView") # <-- please edit this path to your "RNAView" repo directory.

# ---------------------------------------------------------------------------------------------------



# ------------- Do not edit below this line. ------------------
RNAVIEW_PATH = RNAVIEW / "bin"
PseudoKnotVisualizer_DIR = Path(__file__).parent
INTERMEDIATE_DIR = PseudoKnotVisualizer_DIR / "intermediate"
# -------------------------------------------------------------