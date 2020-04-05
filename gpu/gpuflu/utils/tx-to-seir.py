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

    num_time_steps = max([data[:,TEXPOSED].max(), data[:,TINFECTIOUS].max(), data[:,TRECOVERED].max()]) + 1
    seir = np.zeros((num_time_steps,4))
    seir[:,SUSCEPTIBLE] = rows
    for row in range(rows):
        texp = max(data[row,TEXPOSED],    0)
        tinf = max(data[row,TINFECTIOUS], 0)
        trec = max(data[row,TRECOVERED],  0)
        if (texp == 0) and (tinf == 0) and (trec == 0):
            continue
        for t in range(texp, tinf):
            seir[t,EXPOSED] += 1
        for t in range(tinf, trec):
            seir[t,INFECTIOUS] += 1
        seir[trec,RECOVERED] += 1
    
    seir[:,RECOVERED] = np.cumsum(seir[:,RECOVERED])
    seir[:,SUSCEPTIBLE] = rows - seir[:,EXPOSED] - seir[:,INFECTIOUS] - seir[:,RECOVERED]

    np.savetxt(output_filename, np.hstack((np.arange(0,num_time_steps).reshape((num_time_steps,1)),seir)), fmt='%u', delimiter=',', header='TIME,SUSCEPTIBLE,EXPOSED,INFECTIOUS,RECOVERED', comments='')

    return


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", default="transmissions.bin")
    parser.add_argument("-o", "--output", default="seir.csv")

    args = parser.parse_args()

    main(args.input, args.output)
