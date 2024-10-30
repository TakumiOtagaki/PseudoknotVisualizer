import argparse

def argparser():
    parser = argparse.ArgumentParser(description='Visualize pseudoknots in RNA sequences')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file containing RNA sequences')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output script file for visualization')
    parser.add_argument(
        '-f', '--format', choices=['chimera', 'pymol'], 
        required=True, help='Format of RNA sequences (chimera or pymol)'
    )
    chimera_group = parser.add_argument_group('chimera options', 'Options specific to Chimera format')
    chimera_group.add_argument(
        '-m', '--model', type=int, default=None, 
        help='Model ID (required if Chimera format is selected)'
    )


    parser.add_argument('-c', '--chain', type=str, default='A', help='Chain ID for RNA sequences')

    return parser.parse_args()

def args_validation(args):
    if args.format.lower() not in ['pymol', 'chimera']:
        raise ValueError("Output format must be either 'pymol' or 'chimera'")
    if not args.input.endswith('.pdb'):
        raise ValueError("Input file must be a PDB file")
    return
