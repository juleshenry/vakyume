import argparse
from equation_permuter import main

def parse_arguments():
    parser = argparse.ArgumentParser(description="Custom Constants Parser")
    parser.add_argument("--tab", default="    ", help="Custom indentation string")
    parser.add_argument("--typehint", default=": float", help="Default type annotation")
    parser.add_argument("--infile", default="in_file.py", help="Input file name")
    args = parser.parse_args()
    parser.add_argument("--outfile", default=f"permuted_{str(args.infile).split('/')[-1][:-3]}", help="Output file name")
    args = parser.parse_args()
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    main(args.tab, args.typehint, args.infile, args.outfile)