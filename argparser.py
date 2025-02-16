import argparse

def argparser():
    parser = argparse.ArgumentParser(description='Visualize pseudoknots in RNA structure')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file containing RNA structure')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output script file for visualization')
    parser.add_argument(
        '-f', '--format', choices=['chimera', 'pymol'], 
        required=True, help='Format of RNA structure (chimera or pymol)'
    )
    chimera_group = parser.add_argument_group('chimera options', 'Options specific to Chimera format')
    chimera_group.add_argument(
        '-m', '--model', type=int, default=None, 
        help='Model ID (required if Chimera format is selected)'
    )

    parser.add_argument('-c', '--chain', type=str, default='A', help='Chain ID for RNA structure, default is A')

    return parser.parse_args()

def args_validation(args):
    if args.format.lower() not in ['pymol', 'chimera']:
        raise ValueError("Output format must be either 'pymol' or 'chimera'")
    # if not args.input.endswith('.pdb'):
    #     raise ValueError("Input file must be a PDB file")
    if args.format.lower() == 'chimera' and args.model is None:
        raise ValueError("Model ID is required for Chimera format")
    return

