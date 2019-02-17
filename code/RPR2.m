%输入
%LCD---导杆的长度
%输入Cx,Cy,Cdx,Cdy,Cddx,Cddy,Bx,By,Bdx,Bdy,Bddx,Bddy----B,C两点的位置，速度，加速度
%输出
%Dx,Dy,Ddx,Ddy,Dddx,Dddy---D点的位置，速度，加速度
%theta_CD,omega_CD,alpha_CD----杆件CD的角度，角速度，角加速度
function [Dx,Dy,Ddx,Ddy,Dddx,Dddy,theta_CD,omega_CD,alpha_CD] = RPR(Cx,Cy,Cdx,Cdy,Cddx,Cddy,Bx,By,Bdx,Bdy,Bddx,Bddy,LCD)
%位置
A = Bx - Cx;
B = By - Cy;
C = 0;
s = sqrt(A*A+B*B-C*C);

if atan2(A,-B) < 0
theta_CD = atan2(A,-B) + 2*pi - atan2(s,C);
else
theta_CD = atan2(A,-B) - atan2(s,C);

Dx = Cx+s*cos(theta_CD)+(LCD-s)*cos(theta_CD);
Dy = Cy+s*sin(theta_CD)+(LCD-s)*sin(theta_CD);
%速度

G4 = (Bx-Cx)*cos(theta_CD)+(By-Cy)*sin(theta_CD);
omega_CD = ((Bdy-Cdy)*cos(theta_CD)-(Bdx-Cdx)*sin(theta_CD))/G4;
Ddx = Cdx-omega_CD*LCD*sin(theta_CD);
Ddy = Cdy+omega_CD*LCD*cos(theta_CD);
ds = ((Bdx-Cdx)*(Bx-Cx)+(Bdy-Cdy)*(Bdy-Cdy))/G4;
%加速度 

G5 = Bddx-Cddx+omega_CD*omega_CD*(Bx-Cx)+2*ds*omega_CD*sin(theta_CD);
G6 = Bddy-Cddy-omega_CD*omega_CD*(By-Cy)-2*ds*omega_CD*cos(theta_CD);
alpha_CD = (G6*cos(theta_CD)-G5*sin(theta_CD))/G4;
Dddx = Cddx-alpha_CD*LCD*sin(theta_CD)-omega_CD*omega_CD*LCD*cos(theta_CD);
Dddy = Cddy+alpha_CD*LCD*cos(theta_CD)-omega_CD*omega_CD*LCD*sin(theta_CD);

end
