#! /usr/bin/env python3

from argparse import ArgumentParser
import numpy as np


def main(population_size=65536, output_filename="ages.bin"):

    print(f"Generating {population_size} samples...")
    ages = np.random.choice(np.arange(3, dtype=np.uint16), population_size, p=[0.25, 0.5, 0.25])
    print(f"Writing samples to '{output_filename}'...")
    with open(output_filename, "wb") as handle:
        ages.tofile(handle)

    return

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--population", type=int, default=65536)
    parser.add_argument("-o", "--output", type=str, default="ages.bin")

    args = parser.parse_args()

    main(args.population, args.output)
