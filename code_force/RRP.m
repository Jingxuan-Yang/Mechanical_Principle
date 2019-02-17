% function: RRP

%input
%Bx,By,Bdx,Bdy,Bddx,Bddy----B: position, velocity and acceleration
%Kx,Ky,Kdx,Kdy,Kddx,Kddy----K (reference point): position, velocity and acceleration
%phi,dphi,ddphi----guide line: angle, angular velocity and angular acceleration
%LBC,LCD----length of links BC and CD

%output
%Cx,Cy,Cdx,Cdy,Cddx,Cddy----C: position, velocity and acceleration
%Dx,Dy,Ddx,Ddy,Dddx,Dddy----D: position, velocity and acceleration
%theta_BC,omega_BC,alpha_BC----BC: angle, angular velocity and angular acceleration
%s,v,a------slider's position, velocity and acceleration on the guide line

function [Cx,Cy,Cdx,Cdy,Cddx,Cddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy,theta_BC,omega_BC,alpha_BC,s,v,a] =...
    RRP(Bx,By,Bdx,Bdy,Bddx,Bddy,Kx,Ky,Kdx,Kdy,Kddx,Kddy,phi,dphi,ddphi,LBC,LCD)

%angle of BC
A0 = (Bx - Kx)*sin(phi) - (By - Ky)*cos(phi);
theta_BC = asin((A0+LCD)/LBC) + phi;

%position of point C
Cx = Bx + LBC*cos(theta_BC);
Cy = By + LBC*sin(theta_BC);

%displacement of the slider
s = (Cx - Kx + LCD*sin(phi))/cos(phi);

%position of the slider
Dx = Kx + s*cos(phi);
Dy = Ky + s*sin(phi);

%angular velocity and velocity of the slider on the guide line
Q1 = Kdx - Bdx - dphi*(s*sin(phi) + LCD*cos(phi));
Q2 = Kdy - Bdy + dphi*(s*cos(phi) - LCD*sin(phi));
Q3 = LBC*sin(theta_BC)*sin(phi) + LBC*cos(theta_BC)*cos(phi);

omega_BC = (-Q1*sin(phi) + Q2*cos(phi))/Q3;

v = -(Q1*LBC*cos(theta_BC) + Q2*LBC*sin(theta_BC))/Q3;

%velocity of point C
Cdx = Bdx - omega_BC*LBC*sin(theta_BC);
Cdy = Bdy + omega_BC*LBC*cos(theta_BC);

%velocity of point D
Ddx = Kdx + v*cos(phi) - s*dphi*sin(phi);
Ddy = Kdy + v*sin(phi) + s*dphi*cos(phi);

%angular velocity of BC and acceleration of the slider on the guide line
Q4 = Kddx - Bddx + omega_BC^2*LBC*cos(theta_BC) - ddphi*(s*sin(phi) + LCD*cos(phi))...
    -dphi^2*(s*cos(phi) - LCD*sin(phi)) - 2*v*dphi*sin(phi);

Q5 = Kddy - Bddy + omega_BC^2*LBC*sin(theta_BC) + ddphi*(s*cos(phi) - LCD*sin(phi))...
    -dphi^2*(s*sin(phi) + LCD*sin(phi)) + 2*v*dphi*cos(phi);

alpha_BC = (-Q4*sin(phi) + Q5*cos(phi))/Q3;
a = (-Q4*LBC*cos(theta_BC) - Q5*LBC*sin(theta_BC))/Q3;

%acceleration of point C
Cddx = Bddx - alpha_BC*LBC*sin(theta_BC) - omega_BC^2*LBC*cos(theta_BC);
Cddy = Bddy + alpha_BC*LBC*cos(theta_BC) - omega_BC^2*LBC*sin(theta_BC);

%acceleration of point D
Dddx = Kddx + a*cos(phi) - s*ddphi*sin(phi) - s*dphi^2*cos(phi) - 2*v*dphi*sin(phi);
Dddy = Kddy + a*sin(phi) + s*ddphi*cos(phi) - s*dphi^2*sin(phi) + 2*v*dphi*cos(phi);

end