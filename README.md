Thuy Do 7/2017
The code is based on the paper: https://people.eecs.berkeley.edu/~vazirani/pubs/arvcacm.pdf
"Geometry, Flows and Graph-Partitioning Algorithms"
This code follow the geometry approach, for expander flow approach I coded
in java.
The motivation of the code is to find a balance cut of a given graph (unweighted undirected)
given G = (V,E). V: set of vertices; E: set of edges.

The main file to run is ARV_main_entry.m.
To run in matlab: in the matlab command window type: ARV_main_entry
Note: The code uses cvx library from http://cvxr.com/cvx/download/ and plug in matlab for being ready before running the code

input: the adjacent matrix of a graph in csv file and a parameter c = 0.2. c should be <1/2; (e.g. graph_10_vertices.csv)
output: The information of the CUT obtained

for example:

There are 5/10 vertices in the cut 
The cut: 

cut =

    10     9     8     5     4

Number of edges in the graph: 14 
Number of edges in the part 1 of the cut: 6 
Number of edges in the part 2 of the cut: 6 
Number of edges crossing the cut: 2 