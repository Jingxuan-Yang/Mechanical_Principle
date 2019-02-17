clc
clear

%initialize the length of links
L1 = 26.5; L2 = 105.6; LCE = 65;
L3 = 67.5; L4 = 99.4;
L5 = 37.4; L6 = 28;

%angle BCE
alpha_BCE = 60;

%driving link
omega_AB = 1; %rad/s
alpha_AB = 0;

dr = pi/180; %ratio of deg. to rad.

%position, velocity and acceleration of known points
Ax = 0; Ay = 0; Adx = 0; Ady = 0; Addx = 0; Addy = 0;
Dx = L4; Dy =0; Ddx = 0; Ddy = 0; Dddx = 0; Dddy = 0;
Gx = 153.5; Gy = 41.7; Gdx = 0; Gdy = 0; Gddx = 0; Gddy = 0;

%rad. to deg.
rd = 180/pi;
deg = 0:1:360;
m = length(deg);

%initialize matrices
Bx = ones(m,1);By = ones(m,1);Bdx = ones(m,1);Bdy = ones(m,1);Bddx = ones(m,1); Bddy = ones(m,1);
Cx = ones(m,1);Cy = ones(m,1);Cdx = ones(m,1);Cdy = ones(m,1);Cddx = ones(m,1); Cddy = ones(m,1);
Ex = ones(m,1);Ey = ones(m,1);Edx = ones(m,1);Edy = ones(m,1);Eddx = ones(m,1); Eddy = ones(m,1);
Fx = ones(m,1);Fy = ones(m,1);Fdx = ones(m,1);Fdy = ones(m,1);Fddx = ones(m,1); Fddy = ones(m,1);

theta_AB = ones(m,1);
theta_BC = ones(m,1); omega_BC = ones(m,1); alpha_BC = ones(m,1);
theta_CD = ones(m,1); omega_CD = ones(m,1); alpha_CD = ones(m,1);
theta_EF = ones(m,1); omega_EF = ones(m,1); alpha_EF = ones(m,1);
theta_FG = ones(m,1); omega_FG = ones(m,1); alpha_FG = ones(m,1);

for n = 1:m

%convert deg. to rad.
theta_AB(n) = deg(n)*dr;

%A->B
[Bx(n),By(n),Bdx(n),Bdy(n),Bddx(n),Bddy(n)] =...
    RR(Ax,Ay,Adx,Ady,Addx,Addy,theta_AB(n),omega_AB,alpha_AB,L1);

%B,D->C
[Cx(n),Cy(n),Cdx(n),Cdy(n),Cddx(n),Cddy(n),...
    theta_BC(n), theta_CD(n),omega_BC(n), omega_CD(n),alpha_BC(n), alpha_CD(n)] = ...
        RRR(Bx(n),By(n),Bdx(n),Bdy(n),Bddx(n),Bddy(n),Dx,Dy,Ddx,Ddy,Dddx,Dddy,L2,L3,0);

%C->E
[Ex(n),Ey(n),Edx(n),Edy(n),Eddx(n),Eddy(n)] =...
    RR(Cx(n),Cy(n),Cdx(n),Cdy(n),Cddx(n),Cddy(n),...
        - alpha_BCE*dr + theta_BC(n),omega_BC(n),alpha_BC(n),LCE);

%G,E->F, GFE clockwise
% [Fx(n),Fy(n),Fdx(n),Fdy(n),Fddx(n),Fddy(n),...
%     theta_FG(n), theta_EF(n),omega_FG(n), omega_EF(n),alpha_FG(n), alpha_EF(n)] = ...
%         RRR(Gx,Gy,Gdx,Gdy,Gddx,Gddy,Ex(n),Ey(n),Edx(n),Edy(n),Eddx(n),Eddy(n),L6,L5,0);

%E,G->F, EFG counterclockwise
[Fx(n),Fy(n),Fdx(n),Fdy(n),Fddx(n),Fddy(n),...
    theta_EF(n), theta_FG(n),omega_EF(n), omega_FG(n),alpha_EF(n), alpha_FG(n)] = ...
        RRR(Ex(n),Ey(n),Edx(n),Edy(n),Eddx(n),Eddy(n),Gx,Gy,Gdx,Gdy,Gddx,Gddy,L5,L6,1);

end

%omit phase transition
theta_EF = unwrap(theta_EF);
theta_FG = unwrap(theta_FG);

%output results
fprintf('Results \n'); fprintf('\n');

%angle, angular velocity and angular acceleration of BC
fprintf('theta_AB|theta_BC|omega_BC|alpha_BC \n');

for n = 1:m
    
fprintf( '%6.2f %6.2f %6.2f %6.2f\n',theta_AB(n)*rd,theta_BC(n)*rd,omega_BC(n),alpha_BC(n));

end

%angle, angular velocity and angular acceleration of CD
fprintf('\n');
fprintf('theta_CD|omega_CD|alpha_CD \n');
fprintf('\n');

for n = 1:m
    
fprintf( '%6.2f %6.2f %6.2f\n',theta_CD(n)*rd,omega_CD(n),alpha_CD(n));

end

%angle, angular velocity and angular acceleration of CE
fprintf('\n');
fprintf('theta_CE|omega_CE|alpha_CE \n');
fprintf('\n');

for n = 1:m
    
fprintf( '%6.2f %6.2f %6.2f\n',theta_BC(n)*rd + 120,omega_CD(n),alpha_CD(n));

end

%angle, angular velocity and angular acceleration of EF
fprintf('\n');
fprintf('theta_EF|omega_EF|alpha_EF \n');
fprintf('\n');

for n = 1:m
    
fprintf( '%6.2f %6.2f %6.2f\n',theta_EF(n)*rd,omega_EF(n),alpha_EF(n));

end

%angle, angular velocity and angular acceleration of FG
fprintf('\n');
fprintf('theta_FG|omega_FG|alpha_FG \n');
fprintf('\n');

for n = 1:m
    
fprintf( '%6.2f %6.2f %6.2f\n',theta_FG(n)*rd,omega_FG(n),alpha_FG(n));

end

%position, velocity and acceleration of point E
fprintf('\n');
fprintf('   Ex  |  Ey  | Edx  | Edy | Eddx | Eddy \n');
fprintf('\n');

for n = 1:m
    
fprintf( '%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',Ex(n),Ey(n),Edx(n),Edy(n),Eddx(n),Eddy(n));

end

% position of E
figure(1)
plot(Ex,Ey,'b');
title('Position of E');
xlabel('Ex/(mm)')
ylabel('Ey/(mm)')

% velocity of E
Ev = ones(m,1);

for n=1:m
    
Ev(n) = sqrt(Edx(n)*Edx(n) + Edy(n)*Edy(n));

end

figure(2)
plot(theta_AB*rd,Ev,'b');
title('Velocity of E');
xlabel('\theta_1/\circ')
ylabel('Ev/(mm/s)')

% acceleration of E
Ea = ones(m,1);

for n=1:m
    
Ea(n) = sqrt(Eddx(n)*Eddx(n) + Eddy(n)*Eddy(n));

end

figure(3)
plot(theta_AB*rd,Ea,'b');
title('Acceleration of E');
xlabel('\theta_1/\circ')
ylabel('Ea/(mm/s^2)')

%Position of F, for correction
figure(4)
plot(Fx,Fy,'b');
title('Position of F');
xlabel('Fx/(mm)')
ylabel('Fy/(mm)')

%Position of C, for correction
figure(5)
plot(Cx,Cy,'b');
title('Position of C');
xlabel('Cx/(mm)')
ylabel('Cy/(mm)')
