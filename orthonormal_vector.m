function out_val = orthonormal_vector(A)
% for example A = [1 2 3 4];
% This function will randomly return a vector which is orthonormal to the input
% vector A.
% given a vector v in an n-dimension space, a vector u is orthonormal to v
% if v is orthogonal to u and |v| = 1 (|v| = (magnitude of v) = norm(v))
% v is orthogonal to u iif dot(v,u) = <v,u> = 0
% so we generate u as follow:
% - normalize v (means that v = v/norm(v)) to make sure that |v| == 1
% - randomly generate a vector in an n-dimension space, called u_prime
% - u = u_prime - <v, u_prime>*v.
% - normalize u
% can see here https://www.youtube.com/watch?v=TRktLuAktBQ&t=141s
n = size(A,2);
v = (1/norm(A))*A;
u_prime = rand(1,n);
u = u_prime - dot(v,u_prime)*v;
u = (1/norm(u))*u;
display(dot(u,v));
out_val = u;




%n = size(A,2);
%ret_val = A;
%tmp = 0;
%for i=1:n-1
%    tmp = tmp + A(i) * A(i);
%end
%ret_val(n) = -tmp/A(n);
%display(ret_val);
%out_val = (1/norm(ret_val,2))*ret_val;


