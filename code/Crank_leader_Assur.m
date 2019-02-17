clc
clear
L1 = 120.0;
L4 = 380.0;

omega1 = 1;
alpha1 = 0;
dr = pi/180; %ratio of deg. to rad.

Ax = 0; Ay = L4; Adx = 0; Ady = 0; Addx = 0; Addy = 0;
Cx = 0; Cy =0; Cdx = 0; Cdy = 0; Cddx = 0; Cddy = 0;

%rad. to deg.
rd = 180/pi;
deg = 0:1:360;
m = length(deg);

for n = 1:m
   theta1 = deg(n)*dr;
   [Bx,By,Bdx,Bdy,Bddx,Bddy] = RR(Ax,Ay,Adx,Ady,Addx,Addy,theta1,omega1,alpha1,L1);
   [s_B,theta_3] = RPR(Bx,By,Bdx,Bdy,Bddx,Bddy,Cx,Cy,Cdx,Cdy,Cddx,Cddy,0,0,600);
    %????
    theta3(n) = theta_3;
    s(n) = s_B;

end

plot(deg,theta3*rd, 'b', deg, s,'k');
legend('theta3','s');
grid on; hold on;