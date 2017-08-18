function A = read_csv_file_graph(fname)
% The function reads the file specified by the parameter fname (file name)
% for example fname = 'graph_5_vertices.csv'
% The file is a unweight and undirect graph stored in a csv file
% The output is an adjacent matrix of the graph
% for example of the input file
% 
%1,1,0,0,1
%1,1,1,0,0
%0,1,1,0,1
%0,0,0,1,1
%1,0,1,1,1
%=> the output matrix should be symmetric , and all entries in the main
%diagonal line are 1
%
A = csvread(fname)