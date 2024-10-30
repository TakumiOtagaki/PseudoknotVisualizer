from PseudoknotVisualizer import clear_intermediate_files, rnaview, colors
from coloring import CLI_coloring_canonical, load_colors_from_json
from argparser import argparser, args_validation
from config import RNAVIEW_PATH, RNAVIEW, PseudoKnotVisualizer_DIR, INTEREMEDIATE_DIR
from rna import PKextractor
from addressRNAviewOutput import extract_base_pairs
from Bio.PDB import PDBParser, PDBIO
import tempfile
import subprocess

def CLI_rnaview(pdb_file, chain_id):
    # chains = cmd.get_chains(pdb_object)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    chain_structure = structure[0]
    selected_chain = None

    for chain in chain_structure:
        if chain.id == chain_id:
            selected_chain = chain
            break
    
    # チェーンが見つからなかった場合
    if selected_chain is None:
        raise ValueError(f"Chain ID {chain_id} not found in {pdb_file}")
    
    io = PDBIO()
    io.set_structure(selected_chain)

    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=INTEREMEDIATE_DIR) as tmp_pdb:
            pdb_path = tmp_pdb.name # tmp.pdb is created and deleted automatically after the block.
            # cmd.save(pdb_path, pdb_object)
            io.save(pdb_path)

            result = subprocess.run(
                [rnaview, pdb_path],
                env={"RNAVIEW": RNAVIEW},
                cwd=INTEREMEDIATE_DIR,
                check=True
            )
            if result.returncode != 0:
                raise Exception("RNAVIEW failed")
    except Exception as e:
        raise Exception("RNAVIEW failed or Exporting PDB failed: " + str(e))

    result_file = INTEREMEDIATE_DIR + pdb_path.split("/")[-1] + ".out"
    valid_bps_df = extract_base_pairs(result_file) # pandas
    print(valid_bps_df)
    BPL = [(row["left_idx"], row["right_idx"]) for _, row in valid_bps_df.iterrows()]

    return BPL

def CLI_PseudoKnotVisualizer(pdb_file, chain_id, format, output_file, model_id):
    BPL = CLI_rnaview(pdb_file, chain_id)
    pdb_id = pdb_file.split("/")[-1].split(".")[0]
    PKlayers = PKextractor(BPL)

    with open(output_file, "w") as f:
        for depth, PKlayer in enumerate(PKlayers):
            # color = str(depth + 1)
            color = colors[str(depth + 1)]
            script = CLI_coloring_canonical(pdb_id, model_id, chain_id, PKlayer, color, format)
            f.write(script)        

    print("Coloring done.")
    print(f"Depth is {len(PKlayers)}")
    print(f"Output script is saved as {output_file}")
    
    # clear_intermediate_files()
    return


def main():
    args = argparser()
    args_validation(args)

    print("PseudoKnotVisualizer started.")
    CLI_PseudoKnotVisualizer(args.input, args.chain, args.format, args.output, args.model)


if __name__ == "__main__":
    main()