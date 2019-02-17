
% Profile design for disc cam with oscillating follower
% date: 2018/6/25
% designer: XuanYuan_huan

% the roller radius is 15mm
% the push-travel motion angle is 53 degrees
% the farthest dwell angle is 0 degrees
% the return motion angle is 53 degrees
% the nearest dwell angle is 254 degrees
% push travel is fifth order polynomial motion
% travel distance is 15 degrees
% return travel is also fifth order polynomial motion
% cam rotates clockwise
% allowable push travel pressure angle is 40 degrees
% allowable return travel pressure angle is 50 degrees

clc; clear;

rd = 180/pi;  %rad. -> deg.
dr = pi/180;  %deg. -> rad.

r0 = 70;      %initial base circle radius
rr = 15;      %roller radius
rc = 8;       %cutter radius
psim = 15*dr; %travel distance
l = 122;      %length of following member bar
a = 180;      %length of AF
L_AC = 360.0; %length of AC

deltar0 = 1;  %base circle radius increase value

alpha1allow = 40*dr; %allowable angle for push travel
alpha2allow = 50*dr; %allowable angle for return travel

delta1 = 53*dr;    %push-travel motion angle
delta2 = 0*dr;       %farthest dwell angle
delta3 = 53*dr;    %motion angle for return travel
delta4 = 254*dr;   %nearest dwell angle

%accumulative angle
delta12 = delta1 + delta2;
delta13 = delta1 + delta2 + delta3;
delta14 = delta1 + delta2 + delta3 + delta4;

deltaDeg = 1;         %angle distance
n = 52;
omega1 = 2*pi*n/60;   %angular velocity of cam
deg = 0:deltaDeg:360; %degree of cam
N = length(deg);      %number of points

% initialize matrices
psi = ones(N,1); omega = ones(N,1); alpha = ones(N,1);x = ones(N,1); y = ones(N,1); 
dx = ones(N,1); dy = ones(N,1); ddx = ones(N,1); ddy = ones(N,1);
xr = ones(N,1); yr = ones(N,1); xc = ones(N,1); yc = ones(N,1);

while 1
    
    psi0 = acos((a^2+l^2-r0^2)/(2*a*l));
    rhomin = 1000;deltarhomin = 0;
    alpha1max = 0; deltaalpha1max = 0;
    alpha2max = 0; deltaalpha2max = 0;

    for n = 1:N

        %push travel
        rdeg = deg(n)*dr; % deg. to rad.
        
        if  rdeg <= delta1
            
            T1 = rdeg/delta1;
            
            psi(n) = psim*(10*T1^3 - 15*T1^4 + 6*T1^5);
           
            omega(n) = 30*psim*omega1*T1^2*(1 - 2*T1 + T1^2)/delta1; 
            dpsi = omega(n);
            
            alpha(n) = 60*psim*omega1^2*T1*(1 - 3*T1 + 2*T1^2)/delta1^2;
            ddpsi = alpha(n);
            
            %pressure angle for push travel
            alpha1 = abs(atan((30*l*psim*(T1^2 - 2*T1^3 + T1^4)/delta1 - a*cos(psi0 + psi(n)) + l)/...
                     (a*sin(psi0 + psi(n)))));
            
                %select max pressure angle for push travel
                if  alpha1 > alpha1max
                    alpha1max = alpha1;
                    deltaalpha1max = rdeg;
                end
                
        %farthest dwell angle
        elseif rdeg > delta1 && rdeg <= delta12
                psi(n) = psim; omega(n) = 0;
                dpsi = omega(n);

        %return travel
        elseif rdeg > delta12 && rdeg <= delta13
            
               degback = rdeg - delta12;             
               T2 = degback/delta3;
               
               psi(n) = psim*(1-(10*T2^3 - 15*T2^4 + 6*T2^5));
               
               omega(n) = - 30*psim*omega1*T2^2*(1 - 2*T2 +T2^2)/(delta3);
               dpsi = omega(n);
               
               alpha(n) = - 60*psim*omega1^2*T2*(1 - 3*T2 +2*T2^2)/(delta3)^2;
               ddpsi = alpha(n);

            %pressure angle for return travel
            alpha2 = abs(atan((l*psim*(-30*(T2^2 - 2*T2^3 + T2^4)/delta3) - a*cos(psi0 + psi(n)) + l)/...
                     (a*sin(psi0 + psi(n)))));
            
                %select max pressure angle for return travel
                if alpha2 > alpha2max
                   alpha2max = alpha2;
                   deltaalpha2max = rdeg;
                end

        %nearest dwell angle
        elseif rdeg > delta13 && rdeg <= delta14
               psi(n) = 0; omega(n) = 0;
               dpsi = omega(n);
        end
        
        %------calculating cam profile curve---------------------
        
        %theoretical cam profile
        x(n) = l*sin(psi0 + psi(n) + rdeg) - a*sin(rdeg);
        y(n) = a*cos(rdeg) - l*cos(psi0 + psi(n) + rdeg);
        
        dx(n) = l*cos(psi0 + psi(n) + rdeg) - a*cos(rdeg);
        dy(n) = l*sin(psi0 + psi(n) + rdeg) - a*sin(rdeg);

        %real cam profile
        stheta = dx(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
        ctheta = -dy(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
        
        %inner envelope profile
        xr(n) = x(n) + rr*ctheta;
        yr(n) = y(n) + rr*stheta;
        
        %cutter profile
        xc(n) = x(n) + (rr - rc)*ctheta;
        yc(n) = y(n) + (rr - rc)*stheta;
        
    end %w.r.t. for

        %If parameters are not suitable, adjust the radius of the base circle
        if alpha1max > alpha1allow || alpha2max > alpha2allow
            r0 = r0 + deltar0;
            continue;
        else
            break;        
        end
        
end %w.r.t. while

%print related parameters and draw cam profile
fprintf('base circle radius\n');
fprintf('%6.4f\n', r0);

fprintf('max angle for push travel, corresponding cam angle\n');
fprintf('%6.4f %6.4f\n',alpha1max*rd,deltaalpha1max*rd);

fprintf('max angle for return travel, corresponding cam angle\n');
fprintf('%6.4f %6.4f\n',alpha2max*rd,deltaalpha2max*rd);

%--------nominal profile points-------------------------------
fprintf('Results: nominal profile points \n' );
fprintf('n   x   y \n');

for i = 1:N
fprintf('%d\t %6.4f\t %6.4f \n', i,x(i),y(i));
end

%--------actual profile points---------------------------------
fprintf('Results: actual profile points \n' );
fprintf('n   x   y \n');

for i = 1:N
fprintf('%d\t %6.4f\t %6.4f \n', i,xr(i),yr(i));
end

%-------------nominal profile----------------------------------
figure(1)
hold on; grid on; axis equal;
title('Design of disc cam with oscillating follower');
xlabel('x/mm');
ylabel('y/mm');
plot(x,y + L_AC, 'r-');
ylim([270,480]);
%------------base circle------------------------------------
ct = linspace(0, 2*pi);
plot(r0*cos(ct), r0*sin(ct) + L_AC, 'g-'); 

%------------actual profile------------------------------------
plot(xr,yr + L_AC, 'b-');

%------------cutter profile------------------------------------
plot(xc,yc + L_AC, 'c-');

legend('nominal profile','base circle','actual profile','cutter profile');

%------------motion diagram------------------------------------
figure(2)
plot(deg,psi*rd,'r-');
title('Motion diagram');
xlabel('\phi/\circ');
ylabel('\psi/\circ');
%-------------end----------------------------------------------

