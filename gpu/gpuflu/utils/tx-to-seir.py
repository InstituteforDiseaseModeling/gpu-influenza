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
        infector = data[row, INFECTOR]
        texp = data[row,TEXPOSED]
        tinf = data[row,TINFECTIOUS]
        trec = data[row,TRECOVERED]
        if infector >= 0:
            seir[texp,EXPOSED] += 1
            seir[tinf,INFECTIOUS] += 1
            seir[trec,RECOVERED] += 1
        else:
            if tinf >= 0:
                seir[0,EXPOSED] += 1
            elif trec >= 0:
                seir[0,INFECTIOUS] += 1
    
    seir[0,SUSCEPTIBLE] = rows - seir[0,EXPOSED] - seir[0,INFECTIOUS]

    for row in range(1, num_time_steps):
        seir[row,SUSCEPTIBLE] = seir[row-1,SUSCEPTIBLE]                        - seir[row,EXPOSED]      # S(t) = S(t-1)         - dE(t)
        seir[row,EXPOSED]     = seir[row-1,EXPOSED]     + seir[row,EXPOSED]    - seir[row,INFECTIOUS]   # E(t) = E(t-1) + dE(t) - dI(t)
        seir[row,INFECTIOUS]  = seir[row-1,INFECTIOUS]  + seir[row,INFECTIOUS] - seir[row,RECOVERED]    # I(t) = I(t-1) + dI(t) - dR(t)
        seir[row,RECOVERED]   = seir[row-1,RECOVERED]   + seir[row,RECOVERED]                           # R(t) = R(t-1) + dR(t)

    np.savetxt(output_filename, np.hstack((np.arange(0,num_time_steps).reshape((num_time_steps,1)),seir)), fmt='%u', delimiter=',', header='TIME,SUSCEPTIBLE,EXPOSED,INFECTIOUS,RECOVERED', comments='')

    return


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", default="transmissions.bin")
    parser.add_argument("-o", "--output", default="seir.csv")

    args = parser.parse_args()

    main(args.input, args.output)
