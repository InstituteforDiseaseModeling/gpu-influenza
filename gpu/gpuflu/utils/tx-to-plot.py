#! /usr/bin/env python3

from argparse import ArgumentParser
from graphviz import Digraph
import numpy as np

INFECTOR = 0
TEXPOSED = 1
TINFECTIOUS = 2
TRECOVERED = 3


def main(input_filename, output_filename, format):

    dot = Digraph( engine="twopi" )

    print(f"Loading data from '{input_filename}...")
    data = np.fromfile(input_filename, dtype=np.int32)
    rows = data.shape[0] // 4
    data = data.reshape((rows, 4))

    print("Creating graph (tree)...")
    for infected in range(rows):
        infector = data[infected, INFECTOR]
        if infector != -1:
            dot.node( str( infected ) )
            dot.node( str( infector ) )
            dot.edge( str(infector), str(infected) )

    print(f"Rendering graph/tree to '{output_filename}'...")
    dot.render( filename=output_filename, format=format )

    return

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", default="transmissions.bin")
    parser.add_argument("-o", "--output", default="txtree.sv")
    parser.add_argument("-f", "--format", default="png", help="jpeg, svg, tif, gif, tiff, png, jpg, pdf")

    args = parser.parse_args()

    main(args.input, args.output, args.format)
