function  M= ARV_Z_matrix(C, c)
%Thuy Do 7/2017
%C is a adjacent matrix of an undirect and unweighted graph
%This function is to find the matriX M such that Mij = <vi,vj> where v1,
%v2, ... vn are points (embedding) on the unit sphere and vi is
%corresponding to vertex i in the graph.
% we want to find v1, v2 ..., vn. But have to find M firstly
%C =[1 1 0 1 0 0
%    1 1 1 1 0 0
%    0 1 1 0 1 1
%    1 0 0 1 1 0
%    0 0 1 1 1 1
%    0 0 1 0 1 1];
V = size(C,1); 
unit_matrix = ones(V);
%c = 0.2;
cvx_begin sdp
    cvx_solver sdpt3
    %variable X(V,V) symmetric  % X is the distance matrix; X(i,j) = |vi-vj|^2 the square of the distance between 2 points vi and vj
    variable Z(V,V) symmetric    
    maximize sum(sum(C.*Z,2),1)   
    
    subject to
    %---------------------------------------------------
    %
    for i = 1:V 
        Z(i,i) == 1;
    end
    
    %---------------------------------------------------
    sum(sum(unit_matrix.* Z,2),1) <= V*V*(1-2*c)*(1-2*c);
    %---------------------------------------------------
    %triangle inequality
    for i=1:V
       for j=1:V
           for k=1:V                
               if ((i ~=j)&&(i~=k)&&(j~=k))
                Z(i,j) + Z(j,k)- Z(i,k)<=1;
               end
           end           
       end
    end           
    %---------------------------------------------------
    %Z(i,j) = <vi,vj> dot product between 2 points in the sphere
    %--------------------------------------------------
    Z == semidefinite(V);  
cvx_end    
M = Z;
    %-------------------------------------------------
    %Points are well spread
    %-------------------------------------------------
    %s = 0;
    %for i = 1:V-1   
%        for j=i+1:V
%            s = s + Z(i,j);
 %       end
%    end
%    s <= V*V - 2*c*(1-c)*V*V;
    %------------------------------------------------
    % X(i,j) = |vi-vj|^2 = <vi,vi>-2*<vi,vj>+<vj,vj>
    %------------------------------------------------
    
%
    %for i=1:V-1
%        for j=i+1:V
 %           X(i,j) > 0; %.1e-05;
            %X(i,j) <= 3.14;
 %       end    
  %  end
    %--------------------------------------------------
%----------------------------------------------------
% compute M, actually M = Z
%----------------------------------------------------
%display(X);
%display(Z);


