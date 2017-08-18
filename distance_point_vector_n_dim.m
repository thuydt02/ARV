function D = distance_point_n_dim_vector(a,p,n)
% return the distance between a point p and a vector n
% assume a is a point in the line with a direction vector n, 
% distance = norm((a-p) - <(a-p),n>*n) wiki: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
%display(a);
%display(p);
ap = a-p; %vector a->p
%dim = size(p,2);
%for i=1:dim
%    ap(i) = a(i) - p(i);
%end
n = (1/norm(n))*n;
dot_apn = dot(ap, n);
D = norm(ap - dot_apn*n);
