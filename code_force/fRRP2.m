
%RRP dynamic statics

%input
%Li,Lj---length of i,j
%xb,yb---position of B
%xc,yc-----position of C
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

%output
%Frxb,Fryb----B's reaction force
%Frxc,Fryc----C's reaction force
%Frd Mrd----D's reaction force(perpendicular constraint force, constraint moment)

function [Frxb,Fryb,Frxc,Fryc,Frd, Mrd] = fRRP2(xb,yb,xc,yc,xsi,ysi,xsj,ysj,...
    ddxsi,ddysi,ddxsj,ddysj,ddthetai,phij,mi,mj,Ji,Fpxi,Fpyi,Ti,Fpxj,Fpyj,Tj)

    g = 9.8;

    %total external force on i, external force,inertia force and gravity
    Fxi = Fpxi - mi*ddxsi;
    Fyi = Fpyi - mi*ddysi - mi*g;
    Mfi = Ti - Ji*ddthetai;

    %total external force on i, external force,inertia force and gravity
    Fxj = Fpxj - mj*ddxsj;
    Fyj = Fpyj - mj*ddysj - mj*g;
    Mfj = Tj;   %without inertia moment

    %define a variable
    FT = Fxi*(yc - ysi) - Fyi*(xc - xsi) + Mfi;

    %reaction force on D
    Frd = ((Fxi + Fxj)*(yc - yb) - (Fyi + Fyj)*(xc - xb) - FT)/((xc - xb)*cos(phij) ...
    +(yc - yb)*sin(phij));

    %reaction moment on D
    Mrd = Fyj*(xc - xsj) - Fxj*(yc - ysj) - Mfj;

    %C's reaction force
    Frxc = Frd*sin(phij) - Fxj;
    Fryc = -Frd*cos(phij) - Fyj;

    %B's reaction force
    Frxb = Frd*sin(phij) - Fxj -Fxi;
    Fryb = -Frd*cos(phij) - Fyj -Fyi;

end