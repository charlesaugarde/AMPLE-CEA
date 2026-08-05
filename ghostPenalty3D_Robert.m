


function [Kj] = ghostPenalty3D(mesh)

%Ghost penalty
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/06/2021
% Description:
% Function to determine the penalty matrix associated with a Ghost penalty
% approach.  This function is based on the paper of:
%
% S. Sticko, G. Ludvigsson, and G. Kreiss. High-order cut finite elements
% for the elastic wave equation. Advances in Computational Mathematics,
% 46:45:1-28, 2020.
%
% Equation (29) is the key equation that is implemented within this
% function.
%
% Note that this function determines the stabilisation term but it needs to
% be factored before it can be applied to the stiffness of mass matrices.
%
%--------------------------------------------------------------------------
% [Kj] = GHOSTPENALTY(mesh)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh structured array. Function requires:
%           - coord : coordinates of the grid nodes (nodes,nD)
%           - etpl  : element topology (nels,nen)
%           - ftpl  : elements connected to boundary faces
%           - fntpl : nodes connected to boundary faces
%--------------------------------------------------------------------------
% Ouput(s);
% Kj     - Ghost penalty matrix
%--------------------------------------------------------------------------
% See also:
% DETBOUNDELEM - boundary element information
% GAUSSQUAD    - Gauss quadrature points
% SHAPEFUNC    - background grid basis functions
% DETELEMDER   - background grid basis function derivatives
%--------------------------------------------------------------------------

ngp  = 2;                                                                   % number of Gaiss points
mesh = detBoundElem(mesh);                                                  % determine boundary element information

coord = mesh.coord;                                                         % nodal coordinates
etpl  = mesh.etpl;                                                          % element topology
bF    = mesh.bF;                                                            % boundary face flag
nf    = sum(bF);                                                            % number of faces
fetpl = mesh.ftpl(bF,:);                                                    % elements connected to boundary faces
fntpl = mesh.fntpl(bF,:);                                                   % nodes connected to boundary faces

if sum(bF) > 0
    fprintf('');
end

nen        = size(etpl,2);                                                  % number of nodes/element
[nodes,nD] = size(coord);                                                   % number of nodes and dimensions
[GpLoc,w]  = GaussQuad(ngp,(nD-1));                                         % Gauss point positions & weights

sJ    = nen*nD*2;                                                           % size of [J]
neJ   = sJ^2;                                                               % no. entries in J
Jrow  = zeros(neJ*nf,1); Jcol = Jrow; Jval = Jrow;                          % zero stiffness storage
npCnt = 0;                                                                  % zero counter

for i = 1:nf
    fn  = fntpl(i,:);                                                       % face node numbers
    fc  = coord(fn.',:);                                                    % face node coordinates
    
    if nD == 2
        len = norm(fc(1,:)-fc(2,:));                                            % face length
        h = len;
        detJ = len / 2;
    elseif nD == 3
        area = norm(fc(2,:)-fc(1,:)) * norm(fc(4,:) - fc(1,:));                                            % face length
        h = max([norm(fc(2,:)-fc(1,:)),norm(fc(4,:) - fc(1,:))]);
        detJ = area / 4;
    end
    
    
    if nD == 2
        n   = [0 -1; 1 0]*(fc(1,:).'-fc(2,:).')/len;                            % normal to face
        nn  = [n(1) 0    0    n(2) ;                                            % normal matrix
            0    n(2) n(1) 0   ].';
    elseif nD == 3
        
        n = cross((fc(2,:)-fc(1,:)),(fc(4,:) - fc(1,:)));
        n = n / norm(n);
        
        
        nn = zeros(9,3);
        nn([1 5 8],1) = n;
        nn([4 2 7],2) = n;
        nn([9 6 3],3) = n;
    end
    m   = nn*nn.';                                                          % normal product
    J = zeros(nen*nD*2);                                                    % face penalty matrix
    for gp = 1:(ngp*(nD-1))
        if nD == 2
            N = shapefunc(2,GpLoc(gp),1);                                       % face shape function
        else
            N = shapefunc(4,GpLoc(gp,:),2);                                       % face shape function
        end
        x = N*coord(fntpl(i,:),:);                                          % global position of Gauss point
        eCp = coord(etpl(fetpl(i,1),:),:);                                  % positive element coordinates
        
        Gp  = detElemDer(eCp,x);                                            % positive element [G+]
        eCn = coord(etpl(fetpl(i,2),:),:);                                  % negative element coordinates
        Gn  = detElemDer(eCn,x);                                            % negative element [G-]
        G = [Gp -Gn];                                                       % combined [G] matrix
        J = J + ((h^3)/3)*(G.'*m*G)*w(gp)*detJ;                              % face penalty
        %J = J + (G.'*m*G)*w(gp);
    end
    if i == 1
        Jref = reduceJ(J);
    end
    edp = repmat((etpl(fetpl(i,1),:)-1)*nD,nD,1)+repmat((1:nD).',1,nen);    % degrees of freedom of nodes (matrix form)
    edp = reshape(edp,1,nen*nD);                                            % degrees of freedom of nodes (vector form)
    edn = repmat((etpl(fetpl(i,2),:)-1)*nD,nD,1)+repmat((1:nD).',1,nen);    % degrees of freedom of nodes (matrix form)
    edn = reshape(edn,1,nen*nD);                                            % degrees of freedom of nodes (vector form)
    ed  = [edp edn];                                                        % combined degrees of freedom
    Jrow(npCnt+1:npCnt+neJ) = repmat(ed.',sJ,1);                            % row position storage
    Jcol(npCnt+1:npCnt+neJ) = repmat(ed  ,sJ,1);                            % column position storage
    Jval(npCnt+1:npCnt+neJ) = J;                                            % stiffness storage
    fprintf("\n%.5f",norm(J-Jref))
    npCnt = npCnt+neJ;                                                      % entries counter
end
nDoF = nodes*nD;                                                            % number of degrees of freedom
Kj   = sparse(Jrow,Jcol,Jval,nDoF,nDoF);                                    % form the global stiffness matrix

end


%Shape function derivatives
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/06/2021
% Description:
% Function to determine the derivatives of a finite element's shape
% functions based on the element coordinates and a global point within the
% element.
%
%--------------------------------------------------------------------------
% [G] = DETELEMDER(eC,x)
%--------------------------------------------------------------------------
% Input(s):
% eC     - element coordinates (nen,nD)
% x      - global position for derivative evaluation (nD,1)
%--------------------------------------------------------------------------
% Ouput(s);
% G      - strain-displacement matrix
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------


function [G] = detElemDer(eC,x)

[nen,nD] = size(eC);                                                        % number of element nodes and dimensions
GpLoc  = findLocalCoord(x,eC);                                              % local position within the element
[~,dN] = shapefunc(nen,GpLoc,nD);                                           % local shape function derivatives
JT     = dN*eC;                                                             % Jacobian (volume)
dN     = JT\dN;                                                             % global shape function derivatives
if nD==2
    G = zeros(4,nen*nD);                                                    % stain displacement matrix
    G([1 4],1:2:end) = dN;
    G([3 2],2:2:end) = dN;
elseif nD == 3
    G = zeros(9,nen*nD);
    G([1 5 8],1:nD:end) = dN;
    G([4 2 7],2:nD:end) = dN;
    G([9 6 3],3:nD:end) = dN;
end

end


%Gauss point positions and weights
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/06/2021
% Description:
% Function to provide the positions and weights for and Gauss quadrature
% for 1D integration.
%
%--------------------------------------------------------------------------
% [GpLoc,w] = GAUSSQUAD(ngp)
%--------------------------------------------------------------------------
% Input(s):
% ngp    - number of Gauss points
%--------------------------------------------------------------------------
% Ouput(s);
% GpLoc  - Gauss point locations
% w      - Gauss point weights
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

function [GpLoc,w] = GaussQuad(ngp,nD)

if nD == 1
if ngp ==1                                                                  % one Gauss point
    GpLoc = 0;                                                              % location
    w     = 2;                                                              % weight
elseif ngp==2                                                               % two Gauss points
    GpLoc = [-1 1]/sqrt(3);
    w     = [1 1];
elseif ngp==3                                                               % three Gauss points
    GpLoc = [-1 0 1]*sqrt(3/5);
    w     = [5 8 5]/9;
end
elseif nD == 2
GpLoc = zeros(ngp,nD);
w = zeros(ngp,nD);
if ngp ==1                                                                  % one Gauss point
    gpLoc = 0;                                                              % location
    W     = 2;                                                              % weight
elseif ngp==2                                                               % two Gauss points
    gpLoc = [-1 1]/sqrt(3);
    W     = [1 1];
elseif ngp==3                                                               % three Gauss points
    gpLoc = [-1 0 1]*sqrt(3/5);
    W     = [5 8 5]/9;
end
for i = 1:ngp
    for j = 1:ngp
        GpLoc((i-1)*ngp + j,1) = gpLoc(i);
        GpLoc((i-1)*ngp + j,2) = gpLoc(j);
        
        w((i-1)*ngp + j,1) = W(i);
        w((i-1)*ngp + j,2) = W(j);
    end
end
end


end

