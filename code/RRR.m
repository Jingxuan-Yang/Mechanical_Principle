%function: RRR

%input
%Bx,By,Bdx,Bdy,Bddx,Bddy----B: position, velocity and acceleration
%Dx,Dy,Ddx,Ddy,Dddx,Dddy----D: position, velocity and acceleration
%LBC,LCD----length of link BC and CD
%Flag----0: BCD is clockwise; 1: BCD is counterclockwise

%output
%Cx,Cy,Cdx,Cdy,Cddx,Cddy----C: position, velocity and acceleration
%theta_BC, theta_DC,omega_BC----BC: angle, angular velocity and angular acceleration
%omega_DC,alpha_BC,alpha_DC----DC: angle, angular velocity and angular acceleration

function [Cx,Cy,Cdx,Cdy,Cddx,Cddy,theta_BC,theta_DC,omega_BC,omega_DC,alpha_BC,alpha_DC] = ...
    RRR(Bx,By,Bdx,Bdy,Bddx,Bddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy,LBC,LCD,Flag)

LBD = sqrt((Dx-Bx)*(Dx-Bx) + (Dy-By)*(Dy-By));
A = 2*LBC*(Dx-Bx);
B = 2*LBC*(Dy-By);
C = LBC*LBC + LBD*LBD - LCD*LCD;

%clockwise
if Flag == 0
theta_BC = atan2(B,A) - atan2(-sqrt(A*A + B*B - C*C),C);
%theta_BC = 2*atan((B + sqrt(A*A + B*B -C*C))/(A + C));
end

%counterclockwise
if Flag == 1
% theta_BC = atan2(B,A) - atan2(sqrt(A*A + B*B - C*C),C);
theta_BC = 2*atan((B - sqrt(A*A + B*B -C*C))/(A + C));
end

% angle of DC
Cx = Bx + LBC*cos(theta_BC);
Cy = By + LBC*sin(theta_BC);
theta_DC = atan2((Cy - Dy),(Cx - Dx));

%angular velocity of BC and DC, velocity of point C
CBC = LBC*cos(theta_BC);
SBC = LBC*sin(theta_BC);
CDC = LCD*cos(theta_DC);
SDC = LCD*sin(theta_DC);

G1 = CBC*SDC - CDC*SBC;

omega_BC = (CDC*(Ddx - Bdx) + SDC*(Ddy - Bdy))/G1;
omega_DC = (CBC*(Ddx - Bdx) + SBC*(Ddy - Bdy))/G1;

Cdx = Bdx - omega_BC*LBC*sin(theta_BC);
Cdy = Bdy + omega_BC*LBC*cos(theta_BC);

%angular acceleration of BC and DC, acceleration of point C
G2 = Dddx - Bddx + omega_BC*omega_BC*CBC - omega_DC*omega_DC*CDC;
G3 = Dddy - Bddy + omega_BC*omega_BC*SBC - omega_DC*omega_DC*SDC;

alpha_BC = (G2*CDC + G3*SDC)/G1;
alpha_DC = (G2*CBC + G3*SBC)/G1;

Cddx = Bddx - omega_BC*omega_BC*LBC*cos(theta_BC) - alpha_BC*LBC*sin(theta_BC);
Cddy = Bddy - omega_BC*omega_BC*LBC*sin(theta_BC) + alpha_BC*LBC*cos(theta_BC);

end