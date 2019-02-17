
% force analysis for single link

%input
%xa,ya-----position of A
%xb,yb-----position of B
%xsi,ysi---position of centroid si
%ddxsi,ddysi--acceleration of si
%ddthetai--angular acceleration of link i

%mi --- mass of i
%Ji --- moment of inertia of i
%Frxb,Fryb----B's reaction force
%Fpxi,Fpyi----si's external force
%Ti---si's external moment

%output
%Frxa,Frya----A's reaction force
%Ty---balance moment on this link

function [Ty,Frxa,Frya] = ...
    fcrank(xa,ya,xb,yb,xsi,ysi,ddxsi,ddysi,ddthetai,mi,Ji,Frxb,Fryb,Fpxi,Fpyi,Ti)

    g = 9.8;
    Fxi = Fpxi - mi*ddxsi;
    Fyi = Fpyi - mi*ddysi - mi*g;
    Mfi = Ti - Ji*ddthetai;
    
    Frxa = Frxb - Fxi;
    Frya = Fryb - Fyi;
    Ty = Fryb*(xb - xa) - Frxb*(yb - ya) + Fxi*(ysi - ya) - Fyi*(xsi - xa) - Mfi;

end