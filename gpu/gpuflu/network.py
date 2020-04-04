import numpy
import networkx




def nxGraphToNetworkFile( G, connections, filename ):
    """
     Write a networkX graph into a compatible network file.
     Args:
      G           : Graph
      connections : Maximum number of connections per node to be stored in
                    the output file.
      filename    : Name of output file
    """

    # Get graph info and initialize temporary arrays
    invalid = 0
    N = G.number_of_nodes()
    edges = numpy.full( (N, connections), invalid, dtype = numpy.uint32 )
    weights = numpy.full( (N, connections), invalid, dtype = numpy.float16 )
    lastEdge = numpy.zeros( N, dtype=numpy.int )

    # Get edges and weights as numpy arrays
    for edge in list( G.edges() ):
        src  = edge[0]
        dst  = edge[1]
        last = lastEdge[src] 
        if last < connections:     
            edges  [ src, last ] = dst
            weights[ src, last ] = numpy.float16( G[src][dst]["weight"] )
        else:
            print("... [nxGraphToNetworkFile]: edge (",src,",",dst,") dropped")

        lastEdge[src] = last + 1

    # Save edges and weights in output file
    with open( filename, "bw" ) as out:
        edges  .tofile( out )
        weights.tofile( out )

    # Finalize and return
    return




def networkFileToNxGraph( filename, N, connections ):
    """
     Transform a gpu-influenza network file into a networkX graph.
     Args:
      filename    : Name of the network file.
      N           : Number of nodes
      connections : Number of connections per node.
     Returns:
      NetworkX graph.
    """

    # Read file
    nEdges  = N*connections
    edges   = numpy.fromfile( filename, dtype=numpy.uint32, count=nEdges )
    edges   = numpy.reshape( edges, newshape = (population, connections) )
    weights = numpy.fromfile( filename, dtype=numpy.float16, count=nEdges, 
                              offset=nEdges*4)
    weights = numpy.reshape( weights, newshape = (population, connections) )

    # Create graph
    G = networkx.DiGraph()
    G.add_nodes_from( range(0, N) )
    for i in range(0, N):
        for j in range(0, connections):
            G.add_edge( i, edges[i,j], weight=weights[i,j] )

    # Finalize and return
    return G

