%function: PRP

%input
%LBC,LCD----length of links BC and CD
%phi_BK,omega_BK,alpha_BK----guide line BK: angle, angular velocity and angular acceleration
%phi_DM,omega_DM,alpha_DM----guide line DM: angle, angular velocity and angular acceleration
%Kx,Ky,Kdx,Kdy,Kddx,Kddy----K (reference point): position, velocity and acceleration
%Mx,My,Mdx,Mdy,Mddx,Mddy----M (reference point): position, velocity and acceleration

%output
%Bs,Bv,Ba----slider B's position, velocity and acceleration on the guide line
%Ds,Dv,Da----slider D's position, velocity and acceleration on the guide line
%Cx,Cy,Cdx,Cdy,Cddx,Cddy----C: position, velocity and acceleration
%Bx,By,Bdx,Bdy,Bddx,Bddy----B: position, velocity and acceleration
%Dx,Dy,Ddx,Ddy,Dddx,Dddy----D: position, velocity and acceleration

function [Bs,Bv,Ba,Ds,Dv,Da,Cx,Cy,Cdx,Cdy,Cddx,Cddy,...
            Bx,By,Bdx,Bdy,Bddx,Bddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy]=...
    PRP(LBC,LCD,phi_BK,omega_BK,alpha_BK,phi_DM,omega_DM,alpha_DM,...
            Kx,Ky,Kdx,Kdy,Kddx,Kddy,Mx,My,Mdx,Mdy,Mddx,Mddy)

%position of B and D on the guide lines BK and DM
C1 = Mx - Kx - LBC*sin(phi_BK) - LCD*sin(phi_DM);
C2 = My - Ky - LBC*cos(phi_BK) - LCD*cos(phi_DM);
C3 = sin(phi_DM - phi_BK);

Bs = (C1*sin(phi_DM) - C2*cos(phi_DM))/C3;
Ds = (C1*sin(phi_BK) - C2*cos(phi_BK))/C3;

%position of points B, C and D
Cx = Kx + Bs*cos(phi_BK) - LBC*sin(phi_BK);
Cy = Ky + Bs*sin(phi_BK) + LBC*cos(phi_BK);

Bx = Kx + Bs*cos(phi_BK);
By = Ky + Bs*sin(phi_BK);

Dx = Mx + Ds*cos(phi_DM);
Dy = My + Ds*sin(phi_DM);

%velocity of B and D on the guide lines BK and DM
C4 = Mdx - Kdx + omega_BK*(LBC*cos(phi_BK) + Bs*sin(phi_BK))...
     - omega_DM*(LCD*cos(phi_DM) + Ds*sin(phi_DM));
 
C5 = Mdy - Kdy + omega_BK*(LBC*sin(phi_BK) - Bs*cos(phi_BK))...
     - omega_DM*(LCD*sin(phi_DM) - Ds*cos(phi_DM));

Bv = (C4*sin(phi_DM) - C5*cos(phi_DM))/C3;
Dv = (C4*sin(phi_BK) - C5*cos(phi_BK))/C3;

%velocity of points B, C and D

Cdx = Kdx + Bv*cos(phi_BK) - omega_BK*(LBC*cos(phi_BK) + Bs*sin(phi_BK));
Cdy = Kdy + Bv*sin(phi_BK) - omega_BK*(LBC*sin(phi_BK) - Bs*cos(phi_BK));

Bdx = Kdx + Bv*cos(phi_BK) - omega_BK*Bs*sin(phi_BK);
Bdy = Kdy + Bv*sin(phi_BK) - omega_BK*Bs*cos(phi_BK);

Ddx = Mdx + Dv*cos(phi_DM) - omega_DM*Ds*sin(phi_DM);
Ddy = Mdy + Dv*sin(phi_DM) - omega_DM*Ds*cos(phi_DM);

%accerelation of B and D on the guide lines BK and DM
C6 = Kddx - Mddx + alpha_BK*(LBC*cos(phi_BK) + Bs*sin(phi_BK))...
     - alpha_DM*(LCD*cos(phi_DM) + Ds*sin(phi_DM))...
     - omega_BK^2*(LBC*sin(phi_BK) - Bs*cos(phi_BK))...
     + omega_DM^2*(LCD*sin(phi_DM) - Ds*cos(phi_DM))...
     + 2*(Bv*omega_BK*sin(phi_BK) - Dv*omega_DM*sin(phi_DM));
 
C7 = Kddy - Mddy + alpha_BK*(LBC*sin(phi_BK) + Bs*cos(phi_BK))...
     - alpha_DM*(LCD*sin(phi_DM) - Ds*cos(phi_DM))...
     + omega_BK^2*(LBC*cos(phi_BK) + Bs*sin(phi_BK))...
     - omega_DM^2*(LCD*cos(phi_DM) + Ds*sin(phi_DM))...
     - 2*(Bv*omega_BK*cos(phi_BK) - Dv*omega_DM*cos(phi_DM));

Ba = (C6*cos(phi_DM) - C7*sin(phi_DM))/C3;
Da = (C6*sin(phi_BK) - C7*cos(phi_BK))/C3;

%accerelation of points B, C and D
Cddx = Kddx + Ba*cos(phi_BK) - alpha_BK*(LBC*cos(phi_BK) + Bs*sin(phi_BK))...
       + omega_BK^2*(LBC*sin(phi_BK) - Bs*cos(phi_BK))...
       - 2*Bv*omega_BK*sin(phi_BK);
Cddy = Kddy + Ba*sin(phi_BK) - alpha_BK*(LBC*sin(phi_BK) - Bs*cos(phi_BK))...
       - omega_BK^2*(LBC*cos(phi_BK) + Bs*sin(phi_BK))...
       + 2*Bv*omega_BK*cos(phi_BK);

Bddx = Kddx + Ba*cos(phi_BK) - alpha_BK*Bs*sin(phi_BK)...
       - omega_BK^2*Bs*cos(phi_BK) - 2*Bv*omega_BK*sin(phi_BK);
Bddy = Kddy + Ba*sin(phi_BK) + alpha_BK*Bs*cos(phi_BK)...
       - omega_BK^2*Bs*sin(phi_BK) + 2*Bv*omega_BK*cos(phi_BK);

Dddx = Mddx + Da*cos(phi_DM) - alpha_DM*Ds*sin(phi_DM)...
       - omega_DM^2*Ds*cos(phi_DM) - 2*Dv*omega_DM*sin(phi_DM);
Dddy = Mddy + Da*sin(phi_DM) + alpha_DM*Ds*cos(phi_DM)...
       - omega_DM^2*Ds*sin(phi_DM) + 2*Dv*omega_DM*cos(phi_DM);

end