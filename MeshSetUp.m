function FEmesh = MeshSetUp(x,y)
%FUNCTION DICmesh = MeshSetUp(x,y,DICpara)
% Objective: To set up a Q4 uniform FE-mesh  
% ----------------------------------------------
%
%   INPUT: x,y       
%           
%   OUTPUT: Generated FE-mesh {coordinatesFEM, elementsFEM, ...}
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jin.yang@austin.utexas.edu
% Last time updated: 02/2020; 11/2022
% ==============================================


 
%% Build FE-Mesh
M = size(x,2);  N = size(x,1);   % N is vertically in image; M is horizontally in image;
coordinatesFEM = zeros(M*N ,2);

x = x'; y = y';
% I have transpose x and y because Matlab is read matrix in column direction
for i = 1:size(coordinatesFEM,1)
    coordinatesFEM(i,:)  = [x(i),y(i)];
    % x is horizontal position in the image
    % y is vertical position in the image
end

elementsFEM = zeros((M-1)*(N-1),4);
for j = 1:N-1
    for i = 1:M-1
        elementsFEM((j-1)*(M-1)+i ,:) = [(j-1)*(M)+i (j-1)*(M)+i+1 j*(M)+i+1 j*(M)+i];
    end
end
 

% ======== Assign BC values ==========
% -------- dirichlet BC --------
% dirichlet = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))' ;
%             linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))' ;
%             linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))' ;
%             linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))' ];
dirichlet = [];
% -------- neumann BC --------       
% neumann = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))', zeros(M-1,1), -ones(M-1,1) ;
%              linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))', -ones(N-1,1), zeros(N-1,1) ;
%              linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))', ones(N-1,1), zeros(N-1,1) ;
%              linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))', zeros(M-1,1), ones(M-1,1) ];
neumann = [];

 

%% Assign variables
FEmesh.coordinatesFEM = coordinatesFEM;
FEmesh.elementsFEM = elementsFEM;
FEmesh.dirichlet = dirichlet;
FEmesh.neumann = neumann;
FEmesh.x0 = x; FEmesh.y0 = y; FEmesh.M = M; FEmesh.N = N; 
 
        


