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
    parser.add_argument(
        '-a', '--annotator', choices=['DSSR', 'RNAView'],
        default='RNAView', help='Base-pair annotator to use (DSSR or RNAView), default is RNAView'
    )
    parser.add_argument(
        '--include-all', action='store_true', default=False,
        help='Include all base pairs (canonical + non-canonical). Default: canonical only'
    )
    # Hidden legacy options for backward compatibility (do not show in --help)
    parser.add_argument('-p', dest='annotator', choices=['DSSR', 'RNAView'], help=argparse.SUPPRESS)
    parser.add_argument('--parser', dest='annotator', choices=['DSSR', 'RNAView'], help=argparse.SUPPRESS)

    return parser.parse_args()

def args_validation(args):
    if args.format.lower() not in ['pymol', 'chimera']:
        raise ValueError("Output format must be either 'pymol' or 'chimera'")
    # if not args.input.endswith('.pdb'):
    #     raise ValueError("Input file must be a PDB file")
    if args.format.lower() == 'chimera' and args.model is None:
        raise ValueError("Model ID is required for Chimera format")
    
    chosen = getattr(args, 'annotator', None)
    if chosen is None or chosen.upper() not in ['DSSR', 'RNAVIEW']:
        raise ValueError("Annotator must be either 'DSSR' or 'RNAView'")
    return
