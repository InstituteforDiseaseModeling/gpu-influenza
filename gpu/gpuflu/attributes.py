#! /usr/bin/env python3

from argparse import ArgumentParser
import numpy as np


def main(population_size, household_size, output_filename):

    attributes = {}

    # put people in households of size household_size (default 4)
    households = np.arange(population_size, dtype=np.uint16)    # 0, 1, 2, 3, 4, 5, ..., population_size-1
    households = households // household_size                   # 0, 0, 0, 0, 1, 1, ..., (population_size-1)//household_size
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
    parser.add_argument("-h", "--household", type=int, default=4)
    parser.add_argument("-o", "--output", default="attributes.txt")

    args = parser.parse_args()

    main(args.population, args.household, args.output)
