#!/usr/bin/python3

from argparse import ArgumentParser
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def main(filename):

    # Get the edges from the given file
    tx = np.fromfile(filename, np.uint32)
    tx = np.reshape(tx, (tx.shape[0]//3, 3))
    print('Found {0} transmission events in "{1}"'.format(tx.shape[0], filename))

    graph = nx.DiGraph()
    for i in range(tx.shape[0]):
        graph.add_edge(tx[i,0], tx[i,1], time=tx[i,2])

    layout = nx.layout.spring_layout(graph, k=1)
    # layout = nx.layout.circular_layout(graph)

    nx.draw(graph, layout, node_color='blue', arrowstyle='->', arrowsize=10, width=2, edge_cmap=plt.cm.Blues, alpha=0.3)
    nx.draw_networkx_labels(graph, layout, font_size=10, font_family='sans-serif')
    plt.show()
    plt.savefig('network')

    return


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('filename', default='transmissions.bin')
    args = parser.parse_args()

    main(args.filename)
