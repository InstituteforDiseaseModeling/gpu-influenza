#! /usr/bin/env python3

from argparse import ArgumentParser
import numpy as np


def main(population_size=65536, ages_filename="ages.bin", daytime_filename="primary.bin"):

    print(f"Generating {population_size} samples...")
    ages = np.random.choice(np.arange(3, dtype=np.uint16), population_size, p=[0.25, 0.5, 0.25])
    print(f"Writing samples to '{ages_filename}'...")
    with open(ages_filename, "wb") as handle:
        ages.tofile(handle)

    primary = np.zeros_like(ages)

    num_kids = len(np.where(ages == 0)[0])
    num_schools = num_kids >> 10    # 1024 kids per school (too large for preschool, too small for high school)
    print(f"Assigning {num_kids} to {num_schools} schools...")
    primary[np.where(ages == 0)] = np.random.choice(num_schools, num_kids)  # 0...num_schools-1

    num_working = len(np.where(ages == 1)[0])
    num_workplaces = num_working >> 7   # 128 employees per workplace
    print(f"Assigning {num_working} adults to {num_workplaces} workplaces...")
    primary[np.where(ages == 1)] = (np.random.choice(num_workplaces, num_working) + num_schools)    # num_schools...(num_schools+num_workplaces-1)

    num_older = len(np.where(ages == 2)[0])
    num_hangouts = num_older >> 5   # 32 adults per hangout
    print(f"Assigning {num_older} retirees to {num_hangouts} hangouts...")
    primary[np.where(ages == 2)] = (np.random.choice(num_hangouts, num_older) + num_schools + num_workplaces)   # (num_schools+num_worplaces)...(num_schools+num_workplaces+num_hangouts)

    print(f"Writing samples to '{daytime_filename}'...")
    with open(daytime_filename, "wb") as handle:
        primary.tofile(handle)

    return

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--population", type=int, default=65536)
    parser.add_argument("-a", "--ages", type=str, default="ages.bin")
    parser.add_argument("-d", "--daytime", type=str, default="primary.bin")

    args = parser.parse_args()

    main(args.population, args.ages, args.daytime)
