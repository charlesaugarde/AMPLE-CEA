function [etpl,coord] = formCoord1D(nelsx,lx)

%One dimensional finite element grid generation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   06/05/2015
% Description:
% Function to generate a 2D finite element grid of linear quadrilateral 
% elements.
%
%--------------------------------------------------------------------------
% [etpl,coord] = FORMCOORD2D(nelsx,nelsy,lx,ly)
%--------------------------------------------------------------------------
% Input(s):
% nelsx - number of elements in the x direction
% lx    - length in the x direction
%--------------------------------------------------------------------------
% Ouput(s);
% etpl  - element topology
% coord - nodal coordinates
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

nels  = nelsx;                                                              % number of elements
nodes = (nelsx+1);                                                          % number of nodes

%% node generation
coord = zeros(nodes,1);                                                     % zero coordinates 
node  = 0;                                                                  % zero node counter
for i=0:nels
    node=node+1;
    coord(node,1)=lx*i/nelsx;
end
%% element generation
etpl = zeros(nels,2);                                                       % zero element topology
nel  = 0;                                                                   % zero element counter
for nelx=1:nelsx
    nel=nel+1;
    etpl(nel,1)= nel;
    etpl(nel,2)= nel+1;   
end

end