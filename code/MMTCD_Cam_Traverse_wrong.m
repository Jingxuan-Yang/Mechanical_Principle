
% Profile design for disc cam with oscillating follower

% initial base circle radius is 20mm
% the roller radius is 15mm
% the push-travel motion angle is 29.21 degrees
% the farthest dwell angle is 0 degrees
% the return motion angle is 29.21 degrees
% the nearest dwell angle is 301.58 degrees
% push travel is fifth order polynomial motion
% travel distance is 15 degrees
% return travel is also fifth order polynomial motion
% cam rotates clockwise
% push pressure angle is 40 degrees
% return pressure angle is 50 degrees

clc; clear;

rd = 180/pi; %rad. -> deg.
dr = pi/180; %deg. -> rad.

rr = 15; %roller radius
rc = 8; %cutter radius
psim = 15*dr; %travel distance
l = 122; %length of following member bar

L_AC = 360.0; %length of AC

deltar0 = 1; %base circle radius increase value

alpha1allow = 40*dr; %allowable angle for push travel
alpha2allow = 50*dr; %allowable angle for return travel

delta1 = 28.6*dr; %push-travel motion angle
delta2 = 0*dr; %farthest dwell angle
delta3 = 28.6*dr; %motion angle for return travel
delta4 = 302.8*dr; %nearest dwell angle

%accumulative angle
delta12 = delta1 + delta2;
delta13 = delta1 + delta2 + delta3;
delta14 = delta1 + delta2 + delta3 + delta4;

nr = 52;
omega1 = 2*pi*nr/60; %angular velocity of cam

deltaDeg = 1; %angle distance
deg = 0:deltaDeg:360; %degree of cam
m = length(deg); %number of points

r0 = 80:1:120;
s = length(r0);

a = 120:1:200;
t = length(a);

N = m*s*t;
B = s*t;

n = 1;
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
Rp = ones(B,1);
Ra = ones(B,1);
 
 for p = 1:s
     
     rr0 = r0(p);
     
     for q = 1:t
         
         ra = a(q);

         if(rr0 + l > ra && rr0 + ra > l && ra + l > rr0)

            psi0 = acos((ra^2+l^2-rr0^2)/(2*ra*l)); 
            alpha1max = 0; deltaalpha1max = 0;
            alpha2max = 0; deltaalpha2max = 0;

            for k = 1:m
            %push travel
            rdeg = deg(k)*dr; % deg. to rad.

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
%                         rp1 = rr0;
%                         ap1 = ra;
                    end

            %farthest dwell angle
            elseif rdeg > delta1 && rdeg <= delta12
                
                    psi(n) = psim; 
                    omega(n) = 0;
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
%                        rp2 = rr0;
%                        ap2 = ra;
                    end

            %nearest dwell angle
            elseif rdeg > delta13 && rdeg <= delta14
                
                   psi(n) = 0; 
                   omega(n) = 0;
                   dpsi = omega(n);
            end

            n = n + 1;

            end %w.r.t for k
            
            Alpha1max(b) = alpha1max*rd;
            Alpha2max(b) = alpha2max*rd;
%             Rp(b) = rp;
%             Ra(b) = ra;

         else
             break;
             
         end
         
         b = b + 1;

     end
                  
 end
        
%print related parameters and draw cam profile
% fprintf('base circle radius\n');
% fprintf('%6.4f,%6.4f\n', rp,ap);

[M1,P1] = min(Alpha1max);
[M2,P2] = min(Alpha2max);

s1 = floor(P1/t);
t1 = P1 - s1*t;

s2 = floor(P2/t);
t2 = P2 - s2*t;

fprintf('max angle for push travel, corresponding cam angle\n');
fprintf('%6.4f %6f %6f\n',M1,r0(s1),a(t1));

fprintf('max angle for return travel, corresponding cam angle\n');
fprintf('%6.4f %6f %6f\n',M2,r0(s2),a(t2));
