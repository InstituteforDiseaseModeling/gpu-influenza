#! /usr/bin/env python3

from argparse import ArgumentParser
import numpy as np


def main(population_size, output_filename):

    attributes = {}

    # put people in households of size 4
    households = np.arange(population_size, dtype=np.uint16)
    households = households // 4
    attributes["households"] = households

    # assign people to primary daytime locations of size ~512
    # e.g. 8192 >> 9 = 16, so people are assigned to locations 0..15 (each approximately 512 people)
    primaries = np.random.randint(0, population_size >> 9, population_size, dtype=np.uint16)
    attributes["primaries"] = primaries

    with open(output_filename, "wt") as fattributes:
        for key, data in attributes.items():
            filename = f"{key}.bin"
            fattributes.write(f"{key}:{filename}\n")
            with open(filename, "wb") as fdata:
                data.tofile(fdata)

    return

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--population", type=int, default=1024)
    parser.add_argument("-o", "--output", default="attributes.txt")

    args = parser.parse_args()

    main(args.population, args.output)
