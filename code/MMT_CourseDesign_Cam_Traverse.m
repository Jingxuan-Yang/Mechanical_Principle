
% Profile design for disc cam with oscillating follower

% date: 2018/6/25
% designer: XuanYuan_huan

% the roller radius is 15mm
% the push-travel motion angle is 28.6 degrees
% the farthest dwell angle is 0 degrees
% the return motion angle is 28.6 degrees
% the nearest dwell angle is 302.8 degrees
% push travel is fifth order polynomial motion
% travel distance is 15 degrees
% return travel is also fifth order polynomial motion
% cam rotates clockwise
% push pressure angle is 40 degrees
% return pressure angle is 50 degrees

clc; clear;

rd = 180/pi; %rad. -> deg.
dr = pi/180; %deg. -> rad.

rr0 = 10:1:2000; %base circle radius
s = length(rr0);
rr = 15; %roller radius
rc = 8; %cutter radius
psim = 15*dr; %travel distance
l = 122; %length of following member bar
ra = 50:1:2000; %length of AF
t = length(ra);
L_AC = 360.0; %length of AC

deltar0 = 1; %base circle radius increase value

alpha1allow = 40*dr; %allowable angle for push travel
alpha2allow = 50*dr; %allowable angle for return travel

delta1 = 53*dr; %push-travel motion angle
delta2 = 0*dr; %farthest dwell angle
delta3 = 53*dr; %motion angle for return travel
delta4 = 254*dr; %nearest dwell angle

%accumulative angle
delta12 = delta1 + delta2;
delta13 = delta1 + delta2 + delta3;
delta14 = delta1 + delta2 + delta3 + delta4;

deltaDeg = 1; %angle distance
n = 52;
omega1 = 2*pi*n/60; %angular velocity of cam
deg = 0:deltaDeg:360; %degree of cam
N = length(deg); %number of points
B = s*t;%number of combinations
b = 1;

% initialize matrices
psi = ones(N,1);
omega = ones(N,1);
alpha = ones(N,1);
x = ones(N,1);
y = ones(N,1);
dx = ones(N,1);
dy = ones(N,1);
ddx = ones(N,1);
ddy = ones(N,1);
xr = ones(N,1);
yr = ones(N,1);
xc = ones(N,1);
yc = ones(N,1);

Alpha1max = ones(B,1);
Alpha2max = ones(B,1);
Alphamax = ones(B,1);

for p=1:s
    
    r0 = rr0(p);
    
    for q = 1:t
        
        a = ra(q);
        
        if (r0 + l > a && r0 + a > l && a + l > r0)
            
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
                    alpha1 = abs(atan((30*l*psim*(T1^2 - 2*T1^3 + T1^4)/delta1-a*cos(psi0 + psi(n)) + l)/...
                             (a*sin(psi0 + psi(n)))));

                        %select max pressure angle for push travel
                        if  alpha1 > alpha1max
                            alpha1max = alpha1;
                            deltaalpha1max = rdeg;
                            rp1 = r0;
                            ap1 = a;
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
                    alpha2 = abs(atan((l*psim*(-30*(T2^2 - 2*T2^3 + T2^4)/delta3)-a*cos(psi0 + psi(n)) + l)/...
                             (a*sin(psi0 + psi(n)))));

                        %select max pressure angle for return travel
                        if alpha2 > alpha2max
                           alpha2max = alpha2;
                           deltaalpha2max = rdeg;
                           rp2 = r0;
                           ap2 = a;
                        end

                %nearest dwell angle
                elseif rdeg > delta13 && rdeg <= delta14
                       psi(n) = 0; omega(n) = 0;
                       dpsi = omega(n);
                end
%                 %------calculating cam profile curve---------------------
%                 %theoretical cam profile
% 
%                 x(n) = l*sin(psi0 + psi(n) + rdeg) - a*sin(rdeg);
%                 y(n) = a*cos(rdeg) - l*cos(psi0 + psi(n) + rdeg);
% 
%                 dx(n) = l*cos(psi0 + psi(n) + rdeg) - a*cos(rdeg);
%                 dy(n) = l*sin(psi0 + psi(n) + rdeg) - a*sin(rdeg);
% 
%                 %real cam profile
%                 stheta = dx(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
%                 ctheta = -dy(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
% 
%                 %inner envelope profile
%                 xr(n) = x(n) + rr*ctheta;
%                 yr(n) = y(n) + rr*stheta;
% 
%                 %cutter profile
%                 xc(n) = x(n) + (rr - rc)*ctheta;
%                 yc(n) = y(n) + (rr - rc)*stheta;

            end %w.r.t. for

            Alpha1max(b) = alpha1max*rd;
            Alpha2max(b) = alpha2max*rd;
            b = b + 1;
        
        else
            
            Alpha1max(b) = 90;
            Alpha2max(b) = 90;
            b = b + 1;
            
        end %if
        
    end %q
    
end %p

[M1,P1] = min(Alpha1max);
[M2,P2] = min(Alpha2max);

for c = 1:B
    
    if Alpha1max(c) < alpha1allow*rd && Alpha2max(c) < alpha2allow*rd
        
        p = 0.661; %weight of push travel pressure angle
        q = 1 - p; %weight of return travel pressure angle
        
        %figure of merit for both pressure angle
        Alphamax(c) = p*Alpha1max(c) + q*Alpha2max(c);
        
    else
        
        Alphamax(c) = 90; %maximum angle adding 
        
    end
    
end

[M3,P3] = min(Alphamax);

%index for min alpha1max
if (P1/t ~= fix(floor(P1/t)))
    s1 = floor(P1/t) + 1;
    t1 = P1 - (s1 - 1)*t;
else
    s1 = floor(P1/t);
    t1 = t;
end

%index for min alpha2max
if (P2/t ~= fix(floor(P2/t)))
    s2 = floor(P2/t) + 1;
    t2 = P2 - (s2 - 1)*t;
else
    s2 = floor(P2/t);
    t2 = t;
end

%index for min figure of merit-Alphamin
if (P3/t ~= fix(floor(P3/t)))
    s3 = floor(P3/t) + 1;
    t3 = P3 - (s3 - 1)*t;
else
    s3 = floor(P3/t);
    t3 = t;
end

fprintf('max angle for push travel, corresponding r0 and a\n');
fprintf('%6.4f %6.0f %6.0f\n',M1,rr0(s1),ra(t1));

fprintf('max angle for return travel, corresponding r0 and a\n');
fprintf('%6.4f %6.0f %6.0f\n',M2,rr0(s2),ra(t2));

fprintf('max angle for push travel and return travel, corresponding r0 and a\n');
fprintf('%6.4f %6.4f %6.0f %6.0f\n',Alpha1max(P3),Alpha2max(P3),rr0(s3),ra(t3));

