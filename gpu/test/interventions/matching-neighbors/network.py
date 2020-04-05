#! /usr/bin/env python3

from argparse import ArgumentParser
import numpy as np

# TODO - parameterize R0, R0_std, and infectious_duration
R0 = 0.9162907                  # #define DEFAULT_R0                          (0.9162907f) // log(2.5)
R0_STD = 0.4054651              # #define DEFAULT_R0_STD                      (0.4054651f) // log(1.5)
INFECTIOUS_DURATION = 1.609438  # #define DEFAULT_INFECTIOUS_DURATION_MEAN    (1.609438f)  // log(5)


def main(population_size, connection_count, ages_filename, daytime_filename, network_filename):

    # ages = np.fromfile(ages_filename, dtype=np.uint16)
    # assert( len(ages) == population_size )
    print(f"Reading daytime locations from {daytime_filename}...")
    daytime = np.fromfile(daytime_filename, dtype=np.uint16)
    assert( len(daytime) == population_size )

    # initially build a random network, note, there is a small possibility of being connected to oneself
    print(f"Creating random network with {population_size*connection_count} ({population_size} x {connection_count}) links...")
    ids = np.arange(population_size, dtype=np.uint32)   # important, individual IDs are uint32
    network = np.random.choice(ids, population_size*connection_count)
    network = network.reshape((population_size, connection_count))
    weights = np.zeros((population_size, connection_count), dtype=np.float16)   # important, connection weights are float16/half

    temp = 3 * connection_count // 4
    for individual in range(population_size):
        if individual % 1024 == 0:
            print(f"Generating {temp} in-group links for individual {individual}...")
        # age = ages[individual]
        locale = daytime[individual]
        options = ids[daytime == locale]
        count = temp
        if len(options) <= temp:
            # print(f"Individual {individual} only has {len(options)-1} 'friends'.")
            count = len(options)-1
        network[individual,0:count] = np.random.choice(options, count, replace=False)

        # calculate connection weights
        connection = np.random.lognormal(R0, R0_STD) / (INFECTIOUS_DURATION * connection_count)
        weights[individual,:] = connection

    print(f"Writing network to '{network_filename}'")
    with open(network_filename, "wb") as handle:
        network.tofile(handle)
        weights.tofile(handle)

    return

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--population", type=int, default=65536)
    parser.add_argument("-c", "--connections", type=int, default=32)
    parser.add_argument("-a", "--ages", type=str, default="ages.bin")
    parser.add_argument("-d", "--daytime", type=str, default="primary.bin")
    parser.add_argument("-n", "--network", type=str, default="network.bin")

    args = parser.parse_args()

    main(args.population, args.connections, args.ages, args.daytime, args.network)
