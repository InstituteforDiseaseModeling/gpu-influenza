#!/usr/bin/python3

from argparse import ArgumentParser
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def main(filename, population_size, num_connections):

    # Get the edges from the given file
    network = np.fromfile(filename, np.uint32, population_size * num_connections)
    network = np.reshape(network, (population_size, num_connections))

    graph = nx.DiGraph()
    for i in range(population_size):
        for c in range(num_connections):
            graph.add_edge(i, network[i,c])

    layout = nx.layout.spring_layout(graph, k=1)
    # layout = nx.layout.circular_layout(graph)

    nx.draw(graph, layout, node_color='blue', arrowstyle='->', arrowsize=10, width=2, edge_cmap=plt.cm.Blues, alpha=0.3)
    nx.draw_networkx_labels(graph, layout, font_size=10, font_family='sans-serif')
    plt.show()
    plt.savefig('network')

    return


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('filename', default='network.bin')
    parser.add_argument('-p', '--population', required=True, type=int, help='Number of individuals.')
    parser.add_argument('-c', '--connection', required=True, type=int, help='Number of connections per individual.')
    args = parser.parse_args()

    main(args.filename, args.population, args.connection)
