%The code is based on the paper: https://people.eecs.berkeley.edu/~vazirani/pubs/arvcacm.pdf
%"Geometry, Flows and Graph-Partitioning Algorithms"
%This code is for the geometry approach, for expander flow approach I coded
%in java.
%find a balance cut of a graph (unweighted undirected)
% given G = (V,E)
% let the embbedding be v1, v2, ... vn in R^n where n = |V|
% minimize sum {|vi-vj|^2: (i,j) in E } 
% here |vi-vj| : the distance between vi and vj (Euclidian distance)
% |vi - vj| = sqrt [(vi1 - vj1)^2 + (vi2-vj2)^2 ... + (vin-vjn)^2]
% vi(vi1, vi2, ...vin) is a point
% subject to: 
%          + for every i |vi|^2 = 1 . We map vertices -> points in the unit sphere
%          + for every i, j, k |vi-vk|^2 + |vk-vj|^2 >= |vi-vk|^2
%          + sum{|vi-vj|^2: i<j} >= 4c(1-c)n^2
% here c: The cut(S, S_bar) is c-balanced if both |S|, |S_bar| >= c|V|
% so c in the range (0, 1/2]
%-------------------------------------------------------------------
%we use the matrix A for represent the graph G with:
%|V| columns and |E| row
%each row represents 1 edge in E.
%for example (i,j) in E => In A, at some row (assume row k) we set Aki = 1, Akj =
%-1, Akt =0 (for all t# i,j)
%for each edge (i,j) (assume i<j) we count 1 time


clear;
%A =[ 1 -1 0 0 0 0
 %    1 0 0 -1 0 0
%     0 1 -1 0 0 0
 %    0 1 0 -1 0 0
%     0 0 1 0 -1 0
 %    0 0 1 0 0 -1
%     0 0 0 1 -1 0
 %    0 0 0 0 1 -1 ];

%B =[ 0 -1 1
 %    -1 2 -1
  %   1 -1 0 ];

%C =[1 1 0 1 0 0
%   1 1 1 1 0 0
%    0 1 1 0 1 1
%    1 1 0 1 1 0
%    0 0 1 1 1 1
%    0 0 1 0 1 1];
%
%fname = 'graph_5_vertices.csv';
fname = 'graph_10_vertices.csv';
C = read_csv_file_graph(fname);
V = size(C,2); %E = size(A,1);

%display(V); 
c = 0.2;
%---------------------------------------------------
% compute M matrix M(i,j) = <vi,vj>
%---------------------------------------------------
M = ARV_Z_matrix(C,c);
display(M);
%---------------------------------------------------
% compute eigenvalues and eigenvectors of M
%---------------------------------------------------
[evec, eval] = eig(M);
%eval = 
%---------------------------------------------------
% compute u1, u2, .... un are corresponding orthonormal vectors of these
% eigenvectors 
%---------------------------------------------------
%display(evec);
%u = zeros(V);
%for i = 1:V
%    u(i,:) = orthonormal_vector(evec(:,i)');
%end
%display(u);
% compute the embedding v1, v2, .... vn 
% vi(k) = sqrt(eval(k)) * Uk(i)
%---------------------------------------------------
%arr_eval = [];
%d = 5;
arr_eval =[];
for i = 1:V
    arr_eval(i) = eval(i,i);
end
new_evec = [];
[sorted_arr_eval, index_] = sort(arr_eval, 'descend');
display(sorted_arr_eval);
positive_eval_num = 0;
for i=1:V
    if (sorted_arr_eval(i) > 0)
        positive_eval_num = positive_eval_num + 1;
    else
        break;
    end
end;
d = positive_eval_num;
for i=1:d
    %new_evec(:,i) = evec(:,index_(i));
    new_evec(:,i) = evec(:,index_(i))/norm(evec(:,index_(i)),2);
end
for i=1:V
    for k=1:d
        v(i,k) = sqrt(sorted_arr_eval(k)) * new_evec(i,k);
    end
end
%display(v);
% compute a good cut from the embedding
cut = ARV_find_good_cut(v,d);
display(cut);
