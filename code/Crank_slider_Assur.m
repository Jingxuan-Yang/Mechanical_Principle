clc
clear
L1 = 100;
L2 = 300;
omega1 = 10;
alpha1 = 0;
dr = pi/180; %ratio of deg. to rad.
%A???
Ax = 0; Ay = 0; Adx = 0; Ady = 0; Addx = 0; Addy = 0;
%rad. to deg.
rd = 180/pi;
deg = 0:1:720;
m = length(deg);

for n = 1:m
theta1 = deg(n)*dr;
[Bx,By,Bdx,Bdy,Bddx,Bddy] = RR(Ax,Ay,Adx,Ady,Addx,Addy,theta1,omega1,alpha1,L1);

[Cx,Cy,Cdx,Cdy,Cddx,Cddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy,theta_BC,omega_BC,alpha_BC,s,v,a]=...
RRP(Bx,By,Bdx,Bdy,Bddx,Bddy,Ax,Ay,Adx,Ady,Addx,Addy,0,0,0,L2,0);

%????
theta2(n) = theta_BC;
s3(n) = s;

%????
omega2(n) = omega_BC;
v3(n) = v;

%?????
alpha2(n) = alpha_BC;
a3(n) = a;

end

%????
figure(1)
subplot(2,2,1);
[AX,H1,H2] = plotyy(deg,theta2*rd, deg, s3);
legend('theta2','s3');
set(get(AX(1),'ylabel'),'string','????/\circ');
set(get(AX(2),'ylabel'),'string','????/mm');
title('????');
xlabel('????\theta_1/\circ')
grid on; hold on;

subplot(2,2,2);
[AX,H1,H2] = plotyy(deg,omega2, deg, v3);
legend('omega2','v3');
title('????');
xlabel('????\theta_1/\circ')
ylabel('?????/rad\cdots^{-1}')
set(get(AX(2),'ylabel'),'string','????/mm\cdots^{-1}');
grid on; hold on;

subplot(2,2,3);
[AX,H1,H2] = plotyy(deg,alpha2, deg, a3);
legend('alpha2','a3');
title('?????');
xlabel('????\theta_1/\circ')
ylabel('?????/rad\cdots^{-2}')
set(get(AX(2),'ylabel'),'string','????/mm\cdots^{-2}');

grid on; hold on;