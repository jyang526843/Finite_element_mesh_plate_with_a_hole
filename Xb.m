function xyb = Xb(s,Domain)
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

global R ;
r = R ;
global O P1 P2 P3 P4 P5 CMP ;

switch Domain
    case 1
       
        x = O(1)+r*cos(pi/4*s) ;
        y = O(2)+r*sin(pi/4*s) ;
           
    case 2

        x = O(1)+r*cos(pi/4+pi/4*s) ;
        y = O(2)+r*sin(pi/4+pi/4*s) ;
       
end

xyb = [x ; y] ;

