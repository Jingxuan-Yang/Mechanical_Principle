%输入
%Dx,Dy,Ddx,Ddy,Dddx,Dddy D点位置，速度，加速度
%Kx,Ky,Kdx,Kdy,Kddx,Kddy K点位置，速度，加速度
%phi,dphi,ddphi 导路角位置，速度，加速度
%LBC,LCD 杆长
%输出
%Ex,Ey,Edx.Edy.Eddx,Eddy--E点的位置，速度，加速度
%s,v,a------滑块沿导路的位置，速度，及加速度
function [Ex,Ey,Edx,Edy,Eddx,Eddy,s,v,a] = RRP(Dx,Dy,Ddx,Ddy,Dddx,Dddy,Kx,Ky,Kdx,Kdy,Kddx,Kddy,phi,dphi,ddphi,LDE,LEE)
A0 = (Dx - Kx)*sin(phi) - (Dy - Ky)*cos(phi);
theta_BC = asin((A0+LEE)/LDE) + phi;
Cx = Dx + LDE*cos(theta_BC);
Cy = Dy + LDE*sin(theta_BC);
%滑块位移
s = (Cx - Kx + LEE*sin(phi))/cos(phi);
%滑块位置方程
Ex = Kx + s*cos(phi);
Ey = Ky + s*sin(phi);
%滑块速度
Q1 = Kdx - Ddx - dphi*(s*sin(phi) + LEE*cos(phi));
Q2 = Kdy - Ddy + dphi*(s*cos(phi) - LEE*sin(phi));
Q3 = LDE*sin(theta_BC)*sin(phi) + LDE*cos(theta_BC)*cos(phi);
omega_BC = (-Q1*sin(phi) + Q2*cos(phi))/Q3;
v = -(Q1*LDE*cos(theta_BC) + Q2*LDE*sin(theta_BC))/Q3;
%D点的速度
Edx = Kdx + v*cos(phi) - s*dphi*sin(phi);
Edy = Kdy + v*sin(phi) + s*dphi*cos(phi);
%加速度，LBC杆的角加速度和滑块沿导程的加速度
Q4 = Kddx - Dddx + omega_BC*omega_BC*LDE*cos(theta_BC) - ddphi*(s*sin(phi) + LEE*cos(phi))-dphi^2*(s*cos(phi) - LEE*sin(phi)) - 2*v*dphi*sin(phi);
Q5 = Kddy - Dddy + omega_BC*omega_BC*LDE*sin(theta_BC) + ddphi*(s*cos(phi) - LEE*sin(phi))-dphi^2*(s*sin(phi) + LEE*sin(phi)) + 2*v*dphi*cos(phi);
alpha_BC = (-Q4*sin(phi) + Q5*cos(phi))/Q3;
a = (-Q4*LDE*cos(theta_BC) - Q5*LDE*sin(theta_BC))/Q3;
%E点加速
Eddx = Kddx + a*cos(phi) - s*ddphi*sin(phi) - s*dphi^2*cos(phi) - 2*v*dphi*sin(phi);
Eddy = Kddy + a*sin(phi) + s*ddphi*cos(phi) - s*dphi^2*sin(phi) + 2*v*dphi*cos(phi);
end