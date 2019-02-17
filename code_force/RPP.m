%function: RPP

%input
%Bx,By,Bdx,Bdy,Bddx,Bddy----B: position, velocity and acceleration
%Kx,Ky,Kdx,Kdy,Kddx,Kddy----K (reference point): position, velocity and acceleration
%LBC----length of link BC
%phi_DK,omega_DK,alpha_DK----guide line DK: angle, angular velocity and angular acceleration
%delta----angle between two guide lines DK and CD

%output
%Cx,Cy,Cdx,Cdy,Cddx,Cddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy----C,D: position, velocity and acceleration
%Cs,Cv,Ca----slider C's position, velocity and acceleration on the guide line
%Ds,Dv,Da----slider D's position, velocity and acceleration on the guide line

function [Cs,Cv,Ca,Cx,Cy,Cdx,Cdy,Cddx,Cddy,Ds,Dv,Da,Dx,Dy,Ddx,Ddy,Dddx,Dddy] =...
    RPP(Bx,By,Bdx,Bdy,Bddx,Bddy,Kx,Ky,Kdx,Kdy,Kddx,Kddy,LBC,phi_DK,omega_DK,alpha_DK,delta)

%position of C
B1 = Bx - Kx + LBC*sin(phi_DK + delta);
B2 = By - Ky - LBC*cos(phi_DK + delta);
B3 = sin(delta);

Cx = Bx + LBC*sin(phi_DK + delta);
Cy = By - LBC*cos(phi_DK + delta);

Cs = (B2*cos(phi_DK) - B1*sin(phi_DK))/B3;
Ds = (B1*sin(phi_DK + delta) - B2*cos(phi_DK + delta));

%position of D
Dx = Kx + Ds*cos(phi_DK);
Dy = Ky + Ds*sin(phi_DK);

%velocity of C and D on guide lines DK and CD
B4 = Bdx - Kdx + omega_DK*(LBC*cos(phi_DK + delta)...
                + Cs*sin(phi_DK + delta) + Ds*sin(phi_DK));
            
B5 = Bdy - Kdy + omega_DK*(LBC*sin(phi_DK + delta)...
                - Cs*cos(phi_DK + delta) + Ds*cos(phi_DK));
            
Cv = (B5*cos(phi_DK) - B4*sin(phi_DK))/B3;
Dv = (B4*sin(phi_DK + delta) - B5*cos(phi_DK + delta))/B3;

%velocity of point C
Cdx = Bdx + omega_DK*LBC*cos(phi_DK + delta);
Cdy = Bdy + omega_DK*LBC*sin(phi_DK + delta);

%velocity of point D
Ddx = Kdx + Dv*cos(phi_DK) - Ds*omega_DK*sin(phi_DK);
Ddy = Kdy + Dv*sin(phi_DK) + Ds*omega_DK*cos(phi_DK);

%accerelation of C and D on guide lines DK and CD
B6 = Bddx - Kddx...
     +alpha_DK*(LBC*cos(phi_DK + delta)+ Cs*sin(phi_DK + delta) + Ds*sin(phi_DK))...
     -omega_DK^2*(LBC*sin(phi_DK + delta)- Cs*cos(phi_DK + delta) + Ds*cos(phi_DK))...
     +2*omega_DK*(Cv*sin(phi_DK + delta) + Dv*cos(phi_DK + delta));
 
B7 = Bddy - Kddy...
     +alpha_DK*(LBC*sin(phi_DK + delta)- Cs*cos(phi_DK + delta) + Ds*cos(phi_DK))...
     +omega_DK^2*(LBC*cos(phi_DK + delta)+ Cs*sin(phi_DK + delta) + Ds*sin(phi_DK))...
     -2*omega_DK*(Cv*cos(phi_DK + delta) + Dv*sin(phi_DK + delta));

Ca = (B7*cos(phi_DK) - B6*sin(phi_DK))/B3;
Da = (B7*sin(phi_DK + delta) - B6*cos(phi_DK + delta))/B3;

%accerelation of point C
Cddx = Bddx + alpha_DK*LBC*cos(phi_DK + delta) - omega_DK^2*LBC*sin(phi_DK + delta);
Cddy = Bddy + alpha_DK*LBC*sin(phi_DK + delta) + omega_DK^2*LBC*cos(phi_DK + delta);

%accerelation of point D
Dddx = Kddx + Da*cos(phi_DK) - Ds*alpha_DK*sin(phi_DK)...
        - 2*Dv*omega_DK*sin(phi_DK) - Ds*omega_DK^2*cos(phi_DK);
    
Dddy = Kddy + Da*sin(phi_DK) + Ds*alpha_DK*cos(phi_DK)...
        + 2*Dv*omega_DK*cos(phi_DK) - Ds*omega_DK^2*sin(phi_DK);

end