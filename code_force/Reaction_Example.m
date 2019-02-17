
% Dynamic statics

clc
clear
L1 = 0.08; %0.08m
L2 = 0.26;
L3DC = 0.3;
L3 = 0.4;
L5 = 0.46;
xd1 = 0.17;
ya1 = 0.09;

Ls1 = 0; %?1????
Ls2 = L2/2; %?2????
Ls3 = L3DC; %?3????
Ls5 = L5/2;

%????
m1 = 3.6;
m2 = 6;
m3 = 7.2;
m5 = 8.5;
m6 = 8.5;

g = 10;

Js1 = 0.03;
Js2 = 0.08;
Js3 = 0.1;
Js5 = 0.12;

%??????
Fpx6 = -4000;

%???A D???
xa = 0; ya = ya1; dxa = 0; dya = 0; ddxa = 0; ddya = 0;
xd = xd1; yd = 0; dxd = 0; dyd = 0; ddxd = 0; ddyd = 0;

omega1 = 40;
alpha1 = 0;
dr = pi/180;

%rad. to deg.
rd = 180/pi;
deg = 0:1:360;
m = length(deg);

%initialize matrices
Frxe = ones(m,1);Frye = ones(m,1);Frxf = ones(m,1);Fryf = ones(m,1);
Frf = ones(m,1);Mrf = ones(m,1);Mb = ones(m,1);
Frxa = ones(m,1);Frya = ones(m,1);Frxb = ones(m,1);Fryb = ones(m,1);
Frxc = ones(m,1);Fryc = ones(m,1);Frxd = ones(m,1);Fryd = ones(m,1);

for n = 1:m
    
    theta1 = deg(n)*dr; %????
    
    %????????????????
    
    %??RR???B???????????
    [xb,yb,dxb,dyb,ddxb,ddyb] = RR(xa,ya,dxa,dya,ddxa,ddya,theta1,omega1,alpha1,L1);
    
    %??RRR???C???2??3?????
    [xc,yc,dxc,dyc,ddxc,ddyc,theta2,theta3,omega2,omega3,alpha2,alpha3] = ...
    RRR(xb,yb,dxb,dyb,ddxb,ddyb,xd,yd,dxd,dyd,ddxd,ddyd,L2,L3DC,0);

    %??RR???E??????????
    [xe,ye,dxe,dye,ddxe,ddye] = RR(xd,yd,dxd,dyd,ddxd,ddyd,theta3,omega3,alpha3,L3);
    
    %??RRP????5??6?????,???????
    [xf,yf,dxf,dyf,ddxf,ddyf,xf1,yf1,dxf1,dyf1,ddxf1,ddyf1,theta5,omega5,alpha5,s6,v6,a6]=...
    RRP(xe,ye,dxe,dye,ddxe,ddye,xa,ya,dxa,dya,ddxa,ddya,0,0,0,L5,0);

    %????????
%     [theta1*rd theta2*rd theta3*rd theta5*rd];
%     [theta1*rd alpha2 alpha3 alpha5];
    
    %---------------------------------------------------------------------
    %?1???????
    [xs1,ys1,dxs1,dys1,ddxs1,ddys1] = RR(xa,ya,dxa,dya,ddxa,ddya,theta1,omega1,alpha1,Ls1);
    %?2?????????RR??,
    [xs2,ys2,dxs2,dys2,ddxs2,ddys2] = RR(xb,yb,dxb,dyb,ddxb,ddyb,theta2,omega2,alpha2,Ls2);
    %?3?????????C?
    [xs3,ys3,dxs3,dys3,ddxs3,ddys3] = RR(xd,yd,dxd,dyd,ddxd,ddyd,theta3,omega3,alpha3,Ls3);
    %?5??????,???
    [xs5,ys5,dxs5,dys5,ddxs5,ddys5] = RR(xe,ye,dxe,dye,ddxe,ddye,theta5,omega5,alpha5,Ls5);

%     [xs1,ys1,dxs1,dys1,ddxs1,ddys1];
% 
%     %????????
%     [deg(n) ddxs2,ddys2 ddxs3,ddys3 ddxs5,ddys5];
%     [deg(n) a6 ddxf ddyf];
%     [deg(n) xb yb xe ye xf yf];
%     [deg(n) s6 xf yf];
    
    %-----------------------------------------------------------------
    % [Frxb,Fryb,Frxc,Fryc,Frd,Mrd] = fRRP2(xb,yb,xc,yc,xsi,ysi,xsj,ysj,...
    % ddxsi,ddysi,ddxsj,ddysj,ddthetai,phij,mi,mj,Ji,Fpxi,Fpyi,Ti,Fpxj,Fpyj,Tj)
    
    %???????????????????
    %??fRRP??,?E??F??????
    [Frxe(n),Frye(n),Frxf(n),Fryf(n),Frf(n),Mrf(n)] = fRRP2(xe,ye,xf,yf,xs5,ys5,xf,...
    yf,ddxs5,ddys5,ddxf,ddyf,alpha5,0,m5,m6,Js5,0,0,0,Fpx6,0,0);

    %??f3PairLinkExternalForce???????????????????????????
    [Fcvtx3,Fcvty3,Mcvtf3] = f3PairLinkExternalForce(xe,ye,xs3,ys3,Frxe(n),Frye(n));

    %??fRRR???????????????????????
    [Frxb(n),Fryb(n),Frxc(n),Fryc(n),Frxd(n),Fryd(n)] = fRRR2(xb,yb,xc,yc,xd,yd,xs2,...
    ys2,xs3,ys3,ddxs2,ddys2,ddxs3,ddys3,alpha2,alpha3,m2,m3,Js2,Js3,0,0,0,Fcvtx3,Fcvty3,Mcvtf3);

    %??fCrank,????????
    [Mb(n),Frxa(n),Frya(n)] = fcrank( xa,ya,xb,yb,xs1,ys1,ddxs1,ddys1,0,m1,Js1,Frxb(n),Fryb(n),0,0,0);
    
end

fprintf('Results \n'); 
fprintf('Frxa\t\tFrya\tFrxb\t\tFryb\t\tFrxc\tFryc\t\tMb\n');

for n = 1:m
    
    fprintf( '%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',Frxa(n),Frya(n),Frxb(n),Fryb(n),Frxc(n),Fryc(n),Mb(n))
    
end

fprintf('Frxe\t\tFrye\tFrxf\t\tFryf\t\tFrf\t\tMrf\n');

for n = 1:m
    fprintf( '%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',Frxe(n),Frye(n),Frxf(n),Fryf(n),Frf(n),Mrf(n))
end

figure(1)
subplot(2,2,1);
plot(deg,Frxf, 'b',deg,Fryf, 'r');
legend('Frxf','Fryf');
title('Reaction of F');
xlabel('\theta_{AB}/\circ')
ylabel('F/N')
grid on; hold on;

subplot(2,2,2);
plot(deg,Frxb, 'b',deg,Fryb, 'r');
legend('Frxb','Fryb');
title('Reaction of B');
xlabel('\theta_{AB}/\circ')
ylabel('F/N')
grid on; hold on;

subplot(2,2,3);
plot(deg,Frxc, 'b',deg,Fryc, 'r');
legend('Frxc','Fryc');
title('Reaction of C');
xlabel('\theta_{AB}/\circ')
ylabel('F/N')
grid on; hold on;

subplot(2,2,4);
plot(deg,Frf, 'b');
legend('Driving moment');
title('Driving Moment');
xlabel('\theta_{AB}/\circ')
ylabel('M/Nm')
grid on; hold on;