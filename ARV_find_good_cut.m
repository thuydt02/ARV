function [cut, num_edges, num_edges_cut1, num_edges_cut2, num_edges_crossing_cut] = ARV_find_good_cut(C, v, d)
%Thuy Do, 7/2017
%The function is an implementation of ARV algorithm
%To find a good cut from the embedding v
% 
% input: the embedding v1, v2, ... vn in R^d => v is a matrix
%        d is the number of dimension (=n in this situation)
%output: the bipartition of vertices of a graph which has the embedding v
% Assume the cut = (S1, S2) st S1 U S2 = V = set of vertices in the graph
% cut = S1;
% num_edges: the number of edges in the graph = num_edges_in_cut1 +
% num_edges_in_cut2 + num_edges_crossing_cut
% num_edges_in_cut1: the number of edges in S1
% num_edges_in_cut2: the number of edges in S2
% num_edges_crossing_cut: the number of edges crossing S1 and S2


%Pick a random line direct_vec through the origin u (O, Po) where O(0,0,...0); P0(randomly choose))
%O = zeros(1,d);
%direct_vec = ones(1,d);
direct_vec = [];
%sign = 1;
for i=1:d
    threshold = rand;
    new_rnd = rand;
    if (new_rnd>threshold)
        sign = 1;
    else
        sign = -1;
    end
    direct_vec(i) = sign * new_rnd;
end
%direct_vec = rand;
direct_vec = (1/norm(direct_vec))*direct_vec;
%display(direct_vec);
nn = size(v,1);
%Compute P and N
P = []; 
N = [];
splitor = 1/sqrt(d);
ind_P = 1;ind_N = 1;
%display(O);
%display(direct_vec);
%display(v);
%display(nn);
%for i=1:nn
%    if (distance_point_vector_n_dim(O, v(i,:),direct_vec)>=splitor)
%        P(ind_P,:) = v(i,:);ind_P = ind_P + 1;
%    %else
%    end
%    if (distance_point_vector_n_dim(O,v(i,:),direct_vec)<=splitor)
%        N(ind_N,:) = v(i,:);ind_N = ind_N + 1;
%    end
%end

for i=1:nn
    prj = dot(v(i,:),direct_vec);
    if (prj>=splitor)
        P(ind_P,:) = v(i,:);ind_P = ind_P + 1;
    %else
    end
    if (prj<=splitor)
        N(ind_N,:) = v(i,:);ind_N = ind_N + 1;
    end
end
display(P);
display(N);
%Compute S and T

S = P; T = N; delta = 1/sqrt(log(nn));
n1 = size(P,1); n2=size(N,1);
x = 1; y = 1; removed = 0;
%display(n1);display(n2); 
%display(S);display(T); 

while (x <= n1) 
    while (y <=n2) 
        if (dot(S(x,:)-T(y,:),S(x,:)-T(y,:))<= delta)
            % eliminate S(x,:), T(y,:)
            %S(ind_S,:) = P(x,:);T(ind_T,:) = P(x,:);
            [rest_S, ps] = removerows(S,'ind',[x]);
            [rest_T, pt] = removerows(T,'ind',[y]);
            S = rest_S; T = rest_T; removed = 1;
        else
            removed = 0;
        end
        if (removed == 1)
            y = 1;
        else
            y = y + 1;
        end
        n2= size(T,1);        
        if (x > size(S,1))
            break;
        end
    end
    n1 = size(S,1);x = x + 1;
end
%Choose random 0<=r<=delta
display(S);
display(T);

r = rand;
if (r>delta) 
    r = r * delta;
end
display(r);display(delta);
%Specify the cut

tmp_cut = [];
ind_cut = 0;
size_S = size(S,1);
for i=1:nn
    for j=1:size_S
        xx = S(j,:);
        if (dot(v(i,:)-xx,v(i,:)-xx) <= r)
            ind_cut = ind_cut + 1;
            tmp_cut(ind_cut) = i;
            break;
        end
    end
end

[sorted_cut, ind1] = sort(tmp_cut,'descend');
num_edges1 = 0;num_edges2 = 0;num_edges_crossing = 0;
num_edges = 0;
sorted_cut_bar = []; j=0;
for i = 1:nn
    if (ismember(i,sorted_cut) == 0)
        j = j + 1;
        sorted_cut_bar(j) = i;
    end
end
%display(sorted_cut);
%display(sorted_cut_bar);
size_sorted_cut = size(sorted_cut,2);
size_sorted_cut_bar = nn - size_sorted_cut;
for i=1:size_sorted_cut-1
    for j=i+1:size_sorted_cut
        if (C(sorted_cut(i),sorted_cut(j)) == 1)
            num_edges1 = num_edges1 + 1;
        end
    end
end

for i=1:size_sorted_cut_bar-1
    for j=i+1:size_sorted_cut_bar
        if (C(sorted_cut_bar(i),sorted_cut_bar(j)) == 1)
            num_edges2 = num_edges2 + 1;
        end
    end
end

for i=1:size_sorted_cut
    for j=1:size_sorted_cut_bar
        if (C(sorted_cut(i),sorted_cut_bar(j)) == 1)
            num_edges_crossing = num_edges_crossing + 1;
        end
    end
end
cut = sorted_cut;
num_edges_cut1 = num_edges1;
num_edges_cut2 = num_edges2;
num_edges_crossing_cut = num_edges_crossing;
num_edges = num_edges1 + num_edges2 + num_edges_crossing;
%cut = sorted_cut;

