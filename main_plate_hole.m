% To mesh a plate with hole at center using Transfinite Interpolation (TFI)
%{
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Warning : On running this the workspace memory will be deleted. Save if
 any data present before running the code !!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--------------------------------------------------------------------------
 Code originally written by : Siva Srinivas Kolukula, PhD                |
                   Structural Mechanics Laboratory                       |
                   Indira Gandhi Center for Atomic Research              |
                   India                                                 |
 E-mail : allwayzitzme@gmail.com                                         |
 web-link: https://sites.google.com/site/kolukulasivasrinivas/           |
                                                                         |
 Code is modified by: Jin Yang, PhD (2019@Caltech)                       |
 Contact: Jin Yang, jyang526@wisc.edu   -or-   aldicdvc@gmail.com        |
--------------------------------------------------------------------------
%}
% Version 1 : 15 November 2013
% Modified version: June 02, 2021

%%
clc ; clear; close all;

%% Dr. Siva Srinivas Kolukula's code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dimensions of the plate
L = 350-20 - (40+20) ;                % Horizontal length of the plate
B = 350-20 - (0+20) ;                 % Vertical breadth of the plate

% Number of discretizations along xi and eta axis
m = 9 ;
n = 9 ;
%
% Model plate as two regions which lie in first quadrant
global R theta;
R = 55 ;               % Radius of the hole at center
%%%%%%%%%%%%%%%%%%%%%%%Dont change from here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = pi/2 ;          % Quarter angle of the hole
global O P1 P2 P3 P4 P5 CMP ;
O = [0,0] ;                 % Centre of plate and hole
P1 = O + [R 0.] ;           % Edge of the hole and plate
P2 = O + [L/2 0.] ;         % Edge of the plate
P3 = O + [L/2 B/2] ;        % Edge of the plate
P4 = O + [0. B/2] ;         % Edge of the plate
P5 = O + [0. R] ;           % Edge of the hole and plate
CMP = [R*cos(theta/2.) R*sin(theta/2.)] ;
% discretize along xi and eta axis
xi = linspace(0.,1,m) ;
eta = linspace(0.,1.,n) ;
% Number of Domains 
Domain = 2 ;
DX = cell(1,Domain) ;   
DY = cell(1,Domain) ;
for d = 1:Domain        % Loop for two domains lying in first coordinate
    % Initialize matrices in x and y axis
    X = zeros(m,n) ;
    Y = zeros(m,n) ;

    for i = 1:m
        Xi = xi(i) ;
        for j = 1:n
            Eta = eta(j) ;
        
            % Transfinite Interpolation 
            XY = (1-Eta)*Xb(Xi,d)+Eta*Xt(Xi,d)+(1-Xi)*Xl(Eta,d)+Xi*Xr(Eta,d)......
                -(Xi*Eta*Xt(1,d)+Xi*(1-Eta)*Xb(1,d)+Eta*(1-Xi)*Xt(0,d)+(1-Xi)*(1-Eta)*Xb(0,d)) ;
    
            X(i,j) = XY(1) ;
            Y(i,j) = XY(2) ;
        
        end 
    end
    DX{d} = X ;
    DY{d} = Y ;
end
% Arrange the coordinates for each domain
X1 = DX{1} ; Y1 = DY{1} ;       % Grid for first domain
X2 = DX{2} ; Y2 = DY{2} ;       % Grid for second domain
X = [X1 ;X2(m-1:-1:1,:)] ;      % Merge both the domains
Y = [Y1 ;Y2(m-1:-1:1,:)] ;
% Plot 1/4th of the plate
figure(1)
plotgrid(X,Y) ;
% break
% Plot other domains of plate by imaging coordinates
vec = [1 1 ; -1 1 ; -1 -1 ; 1 -1] ;
figure(2)
for quadrant = 1:4
    plotgrid(vec(quadrant,1)*X,vec(quadrant,2)*Y) ;
    hold on
end
    
axis on;


%% Generate finite element mesh nodal coordinates and Q4 elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O2 = [193, 543]; % hole center coordinate

%%%%% top right %%%%%
coord_x_1 = X + O2(1); 
coord_y_1 = Y + O2(2); 

M = size(coord_x_1,1); N = size(coord_x_1,2);
elementsFEM1 = zeros((M-1)*(N-1),4);
for j = 1:N-1
    for i = 1:M-1
        elementsFEM1((j-1)*(M-1)+i ,:) = [(j-1)*(M)+i (j-1)*(M)+i+1 j*(M)+i+1 j*(M)+i];
    end
end

%%%%% top left %%%%%
coord_x_2 = -X + O2(1); 
coord_y_2 = Y + O2(2);
elementsFEM2 = zeros((M-1)*(N-1),4);
for j = 1:N-1
    for i = 1:M-1
        elementsFEM2((j-1)*(M-1)+i ,:) = M*N + [(j-1)*(M)+i (j-1)*(M)+i+1 j*(M)+i+1 j*(M)+i];
    end
end

%%%%% bottom left %%%%%
coord_x_3 = -X + O2(1);  
coord_y_3 = -Y + O2(2);  
elementsFEM3 = zeros((M-1)*(N-1),4);
for j = 1:N-1
    for i = 1:M-1
        elementsFEM3((j-1)*(M-1)+i ,:) = 2*M*N + [(j-1)*(M)+i (j-1)*(M)+i+1 j*(M)+i+1 j*(M)+i];
    end
end

%%%%% bottom right %%%%%
coord_x_4 = X + O2(1);
coord_y_4 = -Y + O2(2); 
elementsFEM4 = zeros((M-1)*(N-1),4);
for j = 1:N-1
    for i = 1:M-1
        elementsFEM4((j-1)*(M-1)+i ,:) = 3*M*N + [(j-1)*(M)+i (j-1)*(M)+i+1 j*(M)+i+1 j*(M)+i];
    end
end

%%%%% Combine coordinates and elements %%%%%
coordinatesFEM = [coord_x_1(:) coord_y_1(:);  coord_x_2(:) coord_y_2(:);  
                  coord_x_3(:) coord_y_3(:);  coord_x_4(:) coord_y_4(:) ];
              
elementsFEM = [elementsFEM1; elementsFEM2; elementsFEM3; elementsFEM4];
 
MeshSkeleton = ones(size(coordinatesFEM,1),1);
PlotMesh_show(MeshSkeleton, coordinatesFEM, elementsFEM);


%% %%%%% Add more elements at the top and the bottom %%%%%
[row1,col1] = find( coordinatesFEM(:,1)>0 );
x_min = min(coordinatesFEM(row1,1));
x_max = max(coordinatesFEM(row1,1));
y_min = min(coordinatesFEM(row1,2));
y_max = max(coordinatesFEM(row1,2));

[row2,col2] = find(coordinatesFEM(:,2) == y_min);
[row3,col3] = find(coordinatesFEM(:,2) == y_max);

xList = unique(coordinatesFEM(row2,1));
yList1 = flipud([y_min:-16:40]');
yList2 = [y_max:16:1040-40]';

[xGrid1,yGrid1] = ndgrid(xList,yList1);
[xGrid2,yGrid2] = ndgrid(xList,yList2);
DICpara.winstepsize = 16;
DICpara.ImgSize = [400,1040];

DICmesh1 = MeshSetUp(xGrid1,yGrid1,DICpara);
DICmesh2 = MeshSetUp(xGrid2,yGrid2,DICpara);

sizeCoordFEM1 = size(coordinatesFEM,1);
sizeCoordFEM2 = size(DICmesh1.coordinatesFEM,1);


%%%%% Combine all the elements %%%%%
coordinatesFEM = [coordinatesFEM; 
                 DICmesh1.coordinatesFEM; 
                 DICmesh2.coordinatesFEM];
elementsFEM = [elementsFEM; 
               sizeCoordFEM1 + DICmesh1.elementsFEM; 
               sizeCoordFEM1 + sizeCoordFEM2 + DICmesh2.elementsFEM];
 

%%%%% Delete repeated nodes %%%%% 
for tempi = 2:size(coordinatesFEM,1)
    if coordinatesFEM(tempi,1)>0
        % Check whether it's repeated or not
        [min1,ind1] = min ( sum( abs( ones(tempi-1,1)*coordinatesFEM(tempi,:) - coordinatesFEM(1:tempi-1,:) ) , 2) );
        
        if min1 < 1e-9
         
            coordinatesFEM(tempi,:) = 0*coordinatesFEM(tempi,:); % Clear this coordinate point
            
            [row1,col1] = find(elementsFEM == tempi);
            for tempj = 1:length(row1)
                elementsFEM(row1(tempj),col1(tempj)) = ind1;
            end
             
        end
    
    end
end


%%%%% Reorder elements %%%%%
for tempi = 1:size(elementsFEM,1)
    tempxy_center = mean(coordinatesFEM(elementsFEM(tempi,:),:));
    tempxy = coordinatesFEM(elementsFEM(tempi,:),:) - ones(4,1)*tempxy_center;
    temptheta = atan2( tempxy(:,2), tempxy(:,1) );
    temptheta(temptheta<0) = temptheta(temptheta<0) + 2*pi;
    
    [~,temptheta_ind] = sort(temptheta);
    elementsFEM(tempi,:) = elementsFEM(tempi,temptheta_ind);
end

      
%%%%% Visualization %%%%%
MeshSkeleton = ones(size( coordinatesFEM,1),1);
PlotMesh_show(MeshSkeleton,  coordinatesFEM,  elementsFEM);
 

%%%%% Save mat file %%%%%
save('plate_hole.mat','coordinatesFEM','elementsFEM');



