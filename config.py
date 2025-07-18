from pathlib import Path


# ------------- Do not edit below this line. ------------------
PseudoKnotVisualizer_DIR = Path(__file__).parent
INTERMEDIATE_DIR = PseudoKnotVisualizer_DIR / "intermediate"
# -------------------------------------------------------------

# RNAView configuration
RNAVIEW_DIR = PseudoKnotVisualizer_DIR / "RNAView"
RNAVIEW_EXEC = RNAVIEW_DIR / "bin/rnaview"

# DSSR configuration
DSSR_DIR = PseudoKnotVisualizer_DIR / "DSSR"
DSSR_EXEC = DSSR_DIR / "x3dna-dssr"