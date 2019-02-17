%RPR杆组计算
%输入
%LCD---导杆的长度
%输入
%Cx,Cy,Cdx,Cdy,Cddx,Cddy,Bx,By,Bdx,Bdy,Bddx,Bddy----C,B两点的位置，速度，加速度
%输出
%Dx,Dy,Ddx,Ddy,Dddx,Dddy---D点的位置，速度，加速度
%theta_CD,omega_CD,alpha_CD----杆件CD的角度，角速度，角加速度
%设置参数
LCD = 579.55;LAB = 90.59;
Cx = 0;Cy = 0;Cdx = 0;Cdy = 0;Cddx = 0;Cddy =0;
omega_1 = 2*pi*49/60;
Ox = 0;Oy = 579.55;Odx = 0;Ody = 0;Oddx = 0 ;Oddy =0;
phi =0 ;dphi=0;ddphi=0;
LDE = 162.27;LEE = 0;
dr = pi/180;
deg = 0:6:360;
m = length(deg);
%初始化矩阵
%t = 0:2*pi/(5.13*60):2*pi/5.13;
%
theta_1 = ones(m,1); 
Bx = ones(m,1); By = ones(m,1); Bdx = ones(m,1); Bdy = ones(m,1);
Bddx = ones(m,1); Bddy = ones(m,1);
Dx = ones(m,1); Dy = ones(m,1); Ddx = ones(m,1); Ddy = ones(m,1);
Dddx = ones(m,1); Dddy = ones(m,1);
Ex = ones(m,1); Ey = ones(m,1); Edx = ones(m,1); Edy = ones(m,1);
Eddx = ones(m,1); Eddy = ones(m,1);
s = ones(1,m);v = ones(m,1);a = ones(m,1);
for n =1:m
theta_1(n) = deg(n)*dr;
Bx(n) = LAB*cos(theta_1(n));
By(n) = 350+LAB*sin(theta_1(n));
Bdx(n) = -LAB.*omega_1*sin(theta_1(n));
Bdy(n) = LAB*omega_1*cos(theta_1(n));
Bddx(n) = -LAB*omega_1*omega_1*cos(theta_1(n));
Bddy(n) = -LAB*omega_1*omega_1*sin(theta_1(n));
%
[Dx(n),Dy(n),Ddx(n),Ddy(n),Dddx(n),Dddy(n),theta_CD,omega_CD,alpha_CD] = RPR2(Cx,Cy,Cdx,Cdy,Cddx,Cddy,Bx(n),By(n),Bdx(n),Bdy(n),Bddx(n),Bddy(n),LCD);
%
[Ex(n),Ey(n),Edx(n),Edy(n),Eddx(n),Eddy(n),s(n),v(n),a(n)] = RRP2(Dx(n),Dy(n),Ddx(n),Ddy(n),Dddx(n),Dddy(n),Ox,Oy,Odx,Ody,Oddx,Oddy,phi,dphi,ddphi,LDE,LEE);
end


