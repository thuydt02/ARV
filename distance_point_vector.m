function D = distance_point_vector(Q, P, v)
% return the distance between a point Q and a vector v
% assume P is a point in the line with a direction vector v, 
% distance = |vec(PQ) x vec(v)|/|vec(v)|
% where vec(v) : vector v
% vec(PQ) x vec(v) = the cross product of PQ and v (the result is a vector, the direction product of 2 vectors)
%|vec(v)| = the length of vector v
PQ = [];
n = size(Q,2);
for i=1:n
    PQ(i) = Q(i) - P(i);    
end
D = norm(cross(PQ,v),2)/norm(v);
