%function: RPR

%input
%Bx,By,Bdx,Bdy,Bddx,Bddy----B: position, velocity and acceleration
%Dx,Dy,Ddx,Ddy,Dddx,Dddy----D: position, velocity and acceleration
%LBC,LDG,LEG----length of links BC, DG and EG

%output
%Cx,Cy,Cdx,Cdy,Cddx,Cddy----C: position, velocity and acceleration
%Ex,Ey,Edx,Edy,Eddx,Eddy----E: position, velocity and acceleration
%theta_EG,omega_EG,alpha_EG----EG: angle, angular velocity and angular acceleration
%s----displacement of the slider on the guide line

function [Cx,Cy,Cdx,Cdy,Cddx,Cddy,Ex,Ey,Edx,Edy,Eddx,Eddy,theta_EG,omega_EG,alpha_EG,s,v,a] = ...
    RPR(Bx,By,Bdx,Bdy,Bddx,Bddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy,LBC,LDG,LEG)

%angle of EG
A = Bx - Dx;
B = By - Dy;
C = LBC + LDG;
s = sqrt(A*A + B*B - C*C);
%theta_EG = atan((B*s + A*C)/(A*s - B*C));

%[-pi,pi],overflow
if atan2(A,-B) < 0
theta_EG = atan2(A,-B) + 2*pi - atan2(s,C);
else
theta_EG = atan2(A,-B) - atan2(s,C);
end

%position of point C
Cx = Bx - LBC*sin(theta_EG);
Cy = By - LBC*cos(theta_EG);

%position of point E
Ex = Cx + (LEG - s)*cos(theta_EG);
Ey = Cy + (LEG - s)*sin(theta_EG);

%angular velocity and velocity of slider on the guide line
G4 = (Bx - Dx)*cos(theta_EG) + (By - Dy)*sin(theta_EG);

omega_EG = ((Bdy - Ddy)*cos(theta_EG) - (Bdx - Ddx)*sin(theta_EG))/G4;
v = ((Bdx - Ddx)*(Bx - Dx) + (Bdy - Ddy)*(By - Dy))/G4;

%velocity of points C and E
Cdx = Bdx - omega_EG*LBC*cos(theta_EG);
Cdy = Bdy - omega_EG*LBC*sin(theta_EG);

Edx = Ddx - omega_EG*(LEG*sin(theta_EG) - LDG*cos(theta_EG));
Edy = Ddy + omega_EG*(LEG*cos(theta_EG) + LDG*sin(theta_EG));

%angular acceleration and accerelation of slider on the guide line
G5 = Bddx - Dddx + omega_EG^2*(Bx - Dx) + 2*v*omega_EG*sin(theta_EG);
G6 = Bddy - Dddy + omega_EG^2*(By - Dy) - 2*v*omega_EG*cos(theta_EG);

alpha_EG = (G6*cos(theta_EG)-G5*sin(theta_EG))/G4;
a = (G5*(Bx - Dx) + G6*(By - Dy))/G4;

%accerelation of points C and E
Cddx = Bddx - alpha_EG*LBC*cos(theta_EG) + omega_EG^2*LBC*sin(theta_EG);
Cddy = Bddy - alpha_EG*LBC*sin(theta_EG) - omega_EG^2*LBC*cos(theta_EG);

Eddx = Dddx - alpha_EG*(LEG*sin(theta_EG) - LDG*cos(theta_EG))...
            - omega_EG^2*(LEG*cos(theta_EG) + LDG*sin(theta_EG));
Eddy = Dddy + alpha_EG*(LEG*cos(theta_EG) + LDG*sin(theta_EG))...
            - omega_EG^2*(LEG*sin(theta_EG) - LDG*cos(theta_EG));

end