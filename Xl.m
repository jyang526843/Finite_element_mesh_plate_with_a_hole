function xyl = Xl(s,Domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code originally written by : Siva Srinivas Kolukula, PhD                |
%                   Structural Mechanics Laboratory                       |
%                   Indira Gandhi Center for Atomic Research              |                                              |
% E-mail : allwayzitzme@gmail.com                                         |
% web-link: https://sites.google.com/site/kolukulasivasrinivas/           |
%                                                                         |
% Code is modified by: Jin Yang, PhD (2019@Caltech)                       |
% Contact: Jin Yang, jyang526@wisc.edu   -or-   aldicdvc@gmail.com        |
%--------------------------------------------------------------------------
% 
% Version 1 : 15 November 2013
% Modified version: June 02, 2021;  May 18, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global O P1 P2 P3 P4 P5 CMP ;

switch Domain
    case 1
     
        x = P1(1)+(P2(1)-P1(1))*s ;
        y = P1(2)+(P2(2)-P1(2))*s ;
        
    case 2
        
        x = P5(1)+(P4(1)-P5(1))*s ;
        y = P5(2)+(P4(2)-P5(2))*s ;    
        
end

xyl = [x ; y] ;