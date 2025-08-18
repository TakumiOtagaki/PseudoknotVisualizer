from pathlib import Path


# ------------- Do not edit below this line. ------------------
PseudoKnotVisualizer_DIR = Path(__file__).parent
INTERMEDIATE_DIR = PseudoKnotVisualizer_DIR / "intermediate"
# -------------------------------------------------------------

# RNAView configuration
# - RNAVIEW_DIR: Directory that contains the RNAView installation.
#   If you installed RNAView elsewhere, set it like: RNAVIEW_DIR = Path("/opt/RNAView")
RNAVIEW_DIR = PseudoKnotVisualizer_DIR / "RNAView"

# - RNAVIEW_EXEC: Path to the RNAView binary.
#   Example (custom path): RNAVIEW_EXEC = Path("/opt/RNAView/bin/rnaview")
RNAVIEW_EXEC = RNAVIEW_DIR / "bin/rnaview"

# DSSR configuration
# - DSSR_EXEC: Path to the x3dna-dssr binary. By default we expect it under this repo's DSSR/ folder.
#   Example (custom path): DSSR_EXEC = Path("/usr/local/bin/x3dna-dssr")
DSSR_EXEC = PseudoKnotVisualizer_DIR / "DSSR" / "x3dna-dssr"