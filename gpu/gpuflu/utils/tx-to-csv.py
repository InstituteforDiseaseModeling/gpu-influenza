#! /usr/bin/env python3

from argparse import ArgumentParser
import numpy as np


def main(input_filename, output_filename):
    data = np.fromfile(input_filename, dtype=np.int32)
    rows = data.shape[0] // 4
    data = data.reshape((rows,4))
    np.savetxt(output_filename, np.hstack((np.arange(0,rows).reshape((rows,1)),data)), fmt='%u', delimiter=',', header='id,infectedById,timeExposed,timeInfectious,timeRecovered', comments='')

    return


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default="transmissions.bin")
    parser.add_argument("-o", "--output", type=str, default="transmissions.csv")

    args = parser.parse_args()

    main(args.input, args.output)
