
% Profile design for disc cam with oscillating follower
% date: 2018/6/25
% designer: XuanYuan_huan

clc
clear

L_AB = 100.62;
L_AC = 360.0;
L_CD = 590.35;
L_DE = 159.40;
L = 578.59;             %vertical distance between C and E

n = 52;                 %r/min
omega_AB = 2*pi*n/60;   %rad/s
alpha_AB = 0;
dr = pi/180;            %ratio of deg. to rad.

Ax = 0; Ay = L_AC; Adx = 0; Ady = 0; Addx = 0; Addy = 0;
Cx = 0; Cy = 0; Cdx = 0; Cdy = 0; Cddx = 0; Cddy = 0;
Kx = 0; Ky = L; Kdx = 0; Kdy = 0; Kddx = 0; Kddy = 0;

%angle of DE
phi = 0; dphi =0; ddphi =0;

%rad. to deg.
rd = 180/pi;
deg = 0:1:360;
m = length(deg);

%initialize matrices
theta_AB = ones(m,1);
theta_BC = ones(m,1);s_BC = ones(m,1);
Bx = ones(m,1);By = ones(m,1);Bdx = ones(m,1);Bdy = ones(m,1);Bddx = ones(m,1); Bddy = ones(m,1);
Dx = ones(m,1);Dy = ones(m,1);Ddx = ones(m,1);Ddy = ones(m,1);Dddx = ones(m,1); Dddy = ones(m,1);
Ex = ones(m,1);Ey = ones(m,1);Edx = ones(m,1);Edy = ones(m,1);Eddx = ones(m,1); Eddy = ones(m,1);

for n = 1:m
    
   theta_AB(n) = deg(n)*dr;
   
   % A->B
   [Bx(n),By(n),Bdx(n),Bdy(n),Bddx(n),Bddy(n)] =...
       RR(Ax,Ay,Adx,Ady,Addx,Addy,theta_AB(n),omega_AB,alpha_AB,L_AB);
   
   %B,C->D
   [~,~,~,~,~,~,Dx(n),Dy(n),Ddx(n),Ddy(n),Dddx(n),Dddy(n),theta_CD,omega_CD,alpha_CD,s_B,v_B,a_B] =...
       RPR(Bx(n),By(n),Bdx(n),Bdy(n),Bddx(n),Bddy(n),Cx,Cy,Cdx,Cdy,Cddx,Cddy,0,0,L_CD);
   
   %D->E
   [~,~,~,~,~,~,Ex(n),Ey(n),Edx(n),Edy(n),Eddx(n),Eddy(n),theta_DE,omega_DE,alpha_DE,s_E,v_E,a_E] =...
        RRP(Dx(n),Dy(n),Ddx(n),Ddy(n),Dddx(n),Dddy(n),Kx,Ky,Kdx,Kdy,Kddx,Kddy,phi,dphi,ddphi,L_DE,0);

   %angle of BC and length of L_BC
   theta_BC(n) = theta_CD;
   s_BC(n) = s_B;

end

% length and angle of BC
figure(1)
plot(deg,theta_BC*rd, 'r', deg, s_BC,'k');
legend('theta3','s');
title('Length and angle of BC');
xlabel('\theta_{AB}/\circ');

% position of E
figure(2)
plot(Ex,Ey,'b');
title('Position of E');
xlabel('Ex/(mm)');
ylabel('Ey/(mm)');

% velocity of E
figure(3)
plot(theta_AB*rd,Edx,'b');
title('Velocity of E');
xlabel('\theta_{AB}/\circ');
ylabel('Ev/(mm/s)');

% acceleration of E
figure(4)
plot(theta_AB*rd,Eddx,'b');
title('Acceleration of E');
xlabel('\theta_{AB}/\circ');
ylabel('Ea/(mm/s^2)');

%for correction
figure(5)
plot(Dx,Dy,'b');
title('Position of D');
xlabel('Dx/(mm)');
ylabel('Dy/(mm)');

% push travel and oscillation angle of BC
figure(6)
plot(theta_BC*rd,Ex);
title('Push travel vs oscillation angle of BC');
xlabel('\theta_{BC}/\circ');
ylabel('Ex/(mm)');

%Position of E vs driving link's turning angle
figure(7)
plot(theta_AB*rd,Ex);
title('Position of E vs driving link''s turning angle');
xlabel('\theta_{AB}/\circ');
ylabel('Ex/(mm)');
