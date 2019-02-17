
%RPR dynamic statics

%input
%Li,Lj---length of i,j
%xb,yb---position of B
%xd,yd-----position of D
%xsi,ysi---i centroid si
%xsj,ysj---j centroid sj
%ddxsi,ddysi---acceleration of si
%ddxsj,ddysj---acceleration of sj
%ddthetai ---i angular acceleration 

%phij ---angle between j and horizontal line
%mi,mj---mass of i,j
%Ji ---moment of inertia of i
%Fpxi,Fpyi----external force on si
%Ti---external moment on i
%Fpxj,Fpyj----external force on sj
%Tj---external moment on j
%SS---length of CG

%output
%Frxb,Fryb----B's reaction force
%Frxd,Fryd----D's reaction force
%Frc Mrc----C's reaction force(perpendicular constraint force, constraint moment)

function [Frxb,Fryb,Frxd,Fryd,Frc,Mrc] = fRPR2(xb,yb,xd,yd,xsi,ysi,xsj,ysj,...
    ddxsi,ddysi,ddxsj,ddysj,ddthetai,ddthetaj,phij,mi,mj,Ji,Jj,Fpxi,Fpyi,Ti,Fpxj,Fpyj,Tj,SS)

    g = 9.8;

    %total external force on i, external force,inertia force and gravity
    Fxi = Fpxi - mi*ddxsi;
    Fyi = Fpyi - mi*ddysi - mi*g;
    Mfi = Ti - Ji*ddthetai;

    %total external force on i, external force,inertia force and gravity
    Fxj = Fpxj - mj*ddxsj;
    Fyj = Fpyj - mj*ddysj - mj*g;
    Mfj = Tj - Jj*ddthetaj;   %without inertia moment

    %define a variable
    MT = Fxi*(yb - ysi) - Fyi*(xb - xsi) + Mfi;
    
    %reaction moment on D
    Mrc = MT;

    %reaction force on D
    Frc = (Fxj*(ysj - yd) - Fyj*(xsj - xd) - Mfj - MT)/SS;

    %B's reaction force
    Frxb = - Frc*sin(phij) - Fxi;
    Fryb = Frc*cos(phij) - Fyi;

    %D's reaction force
    Frxd = Frc*sin(phij) - Fxj;
    Fryd = - Frc*cos(phij) - Fyj;

end