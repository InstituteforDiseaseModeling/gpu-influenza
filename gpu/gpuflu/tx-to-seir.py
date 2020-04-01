#! /usr/bin/env python3

from argparse import ArgumentParser
import numpy as np

INFECTOR = 0
TEXPOSED = 1
TINFECTIOUS = 2
TRECOVERED = 3

SUSCEPTIBLE = 0
EXPOSED = 1
INFECTIOUS = 2
RECOVERED = 3

def main(input_filename, output_filename):

    data = np.fromfile(input_filename, dtype=np.int32)
    rows = data.shape[0] // 4
    data = data.reshape((rows, 4))

    max_time_step = max([data[:,TEXPOSED].max(), data[:,TINFECTIOUS].max(), data[:,TRECOVERED].max()]) + 1
    seir = np.zeros((max_time_step,4))
    seir[:,SUSCEPTIBLE] = rows
    for row in range(rows):
        texp = data[row,TEXPOSED]
        if texp >= 0:
            tinf = data[row,TINFECTIOUS]
            trec = data[row,TRECOVERED]
            seir[texp,SUSCEPTIBLE] -= 1
            seir[texp,EXPOSED] += 1
            seir[tinf,SUSCEPTIBLE] -= 1
            seir[tinf,INFECTIOUS] += 1
            seir[trec,SUSCEPTIBLE] -= 1
            seir[trec,RECOVERED] += 1

    np.savetxt(output_filename, np.hstack((np.arange(0,max_time_step).reshape((max_time_step,1)),seir)), fmt='%u', delimiter=',', header='TIME,SUSCEPTIBLE,EXPOSED,INFECTIOUS,RECOVERED', comments='')

    return


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", default="transmissions.bin")
    parser.add_argument("-o", "--output", default="seir.csv")

    args = parser.parse_args()

    main(args.input, args.output)
