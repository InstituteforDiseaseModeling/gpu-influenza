import pandas
import networkx
import random
import gpuflu.network




# Parameters
filename_in  = "./53033/people.txt"
filename_out = "king_county.bin"
w_household  = 1.0  # Weight for connections in same household
w_school     = 0.9  # ditto school
w_work       = 0.8  # ditto work
max_school_edges = 100
max_work_edges   = 50
seed = 10100


# Initialize random number generator 
random.seed(seed)


# Read input file
print( "... reading input file: ", filename_in )
people = pandas.read_table( filename_in, sep="\t", header=0, 
                            skipinitialspace=True )


# Create graph and add nodes
print( "... adding nodes to graph" )
G = networkx.DiGraph()
G.add_nodes_from( range(0, len(people) ) )


# Add household edges
print( "... adding household edges" )
households = people.groupby( ["sp_hh_id"] )
nHouseholds = len(households.groups)
k = 0
for householdId in households.groups:
    members = households.get_group( householdId ).index
    n = len(members)
    for ii in range(0,n):
        i = members[ii]
        for jj in range(i+1,n):
            j = members[jj]
            G.add_edge( i, j, weight = w_household )
            G.add_edge( j, i, weight = w_household )
    k = k + 1
    if (k % 100000 == 0) or (k == nHouseholds):
        print( "    ", k, "/", nHouseholds )


# Add school edges
print( "... adding school edges" )
schools = people.groupby( ["school_id"] )
nSchools = len(schools.groups)
k = 0
for schoolId in schools.groups:
    members = schools.get_group( schoolId ).index
    n = len(members)
    if n < max_school_edges:
        for ii in range(0,n):
            i = members[ii]
            for j in range(i+1,n):
                j = members[jj]
                G.add_edge( i, j, weight = w_school )
                G.add_edge( j, i, weight = w_school )
    else:
        for i in members:
            selected = random.sample( range(n), max_school_edges )
            for j in selected:
                if i != members[j]:
                    G.add_edge( i, members[j], weight = w_school )
    k = k + 1
    if (k % 100 == 0) or (k == nSchools):
        print( "    ", k, "/", nSchools )
    if (k == nSchools-1):
        break    # Last group is no-school; this could be made robustly
                 # by not processing groups with school_id="X"


# Add workplace edges
print( "... adding workplace edges" )
workplaces = people.groupby( ["work_id"] )
nWorkplaces = len(workplaces.groups)
k = 0
for workplaceId in workplaces.groups:
    members = workplaces.get_group( workplaceId ).index
    n = len(members)
    if n < max_work_edges:
        for ii in range(0,n):
            i = members[ii]
            for jj in range(ii+1,n):
                j = members[jj]
                G.add_edge( i, j, weight = w_work )
                G.add_edge( j, i, weight = w_work )
    else:
        for i in members:
            selected = random.sample( range(n), max_work_edges )
            for j in selected:
                if i != members[j]:
                    G.add_edge( i, members[j], weight = w_work )
    k = k + 1
    if (k % 10000 == 0) or (k == nWorkplaces):
        print( "    ", k, "/", nWorkplaces )
    if (k == nWorkplaces-1):
        break    # Last group is no-workplace; this could be made robustly
                 # by not processing groups with work_id="X"


# Generating gpu-influenza network file
#degree = sorted( [d for n, d in G.degree()] )
degree = [d for n, d in G.degree()] 
population  = len(people)
connections = max( degree ) 
print( "... generating gpu-influenza network file: ", filename_out )
print( "     population  = ", population  )
print( "     connections = ", connections )
gpuflu.network.nxGraphToNetworkFile( G, connections, filename_out )


# Finalize
print( "... graph information: \n", networkx.info(G) ) 
print( "... done" )


