%function: RR

%input
%x,y,dx,dy,ddx,ddy -----A: position, velocity and acceleration
%theta,omega,alpha----link: angle, angular velocity and angular acceleration
%LAB----length of the link AB

%output
%Bx,By,Bdx,Bdy,Bddx,Bddy-----B: position, velocity and acceleration

function [Bx,By,Bdx,Bdy,Bddx,Bddy] = RR(Ax,Ay,Adx,Ady,Addx,Addy,theta,omega,alpha,LAB)

%position
Bx = Ax + LAB*cos(theta);
By = Ay + LAB*sin(theta);

%velocity
Bdx = Adx - omega*LAB*sin(theta);
Bdy = Ady + omega*LAB*cos(theta);

%acceleration
Bddx = Addx - omega*omega*LAB*cos(theta) - alpha*LAB*sin(theta);
Bddy = Addy - omega*omega*LAB*sin(theta) + alpha*LAB*cos(theta);

end