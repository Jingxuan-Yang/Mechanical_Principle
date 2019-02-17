
% Profile design for offset disk cam

% initial base circle radius is 20mm
% the roller radius is 14mm
% the push-travel motion angle is 165 degrees
% the farthest dwell angle is 50 degrees
% the return motion angle is 100 degrees
% the nearest dwell angle is 45 degrees
% push travel is fifth order polynomial motion
% travel distance is 80mm
% return travel is cosine acceleration motion
% cam rotates clockwise
% center of offset circle is on the left side of the push rod
% push pressure angle is 30 degrees
% return pressure angle is 75 degrees
% radius of offset circle: radius of base circle = 1:4

clc; clear;

rd = 180/pi; %rad. -> deg.
dr = pi/180; %deg. -> rad.

r0 = 20; %base circle radius
rr = 14; %roller radius
h = 30; %travel distance
e = 1/4*r0; %offset distance
deltar0 = 1; %base circle radius increase value

alpha1allow = 30*dr; %allowable angle for push travel
alpha2allow = 75*dr; %allowable angle for return travel

rhominallow = 0.3*rr; %allowable minimum radius of curvature

delta1 = 165*dr; %push-travel motion angle
delta2 = 50*dr; %farthest dwell angle
delta3 = 100*dr; %motion angle for return travel
delta4 = 45*dr; %nearest dwell angle

%accumulative angle
delta12 = delta1 + delta2;
delta13 = delta1 + delta2 + delta3;
delta14 = delta1 + delta2 + delta3 + delta4;

deltaDeg = 3; %angle distance
omega1 = 1; %angular velocity of cam
deg = 0:deltaDeg:360; %degree of cam
N = length(deg); %number of points

% initialize matrices
s = ones(N,1);
v = ones(N,1);
a = ones(N,1);
x = ones(N,1);
y = ones(N,1);
dx = ones(N,1);
dy = ones(N,1);
ddx = ones(N,1);
ddy = ones(N,1);
xr = ones(N,1);
yr = ones(N,1);

% rhomin = 1000;deltarhomin = 0; 
% initialize min radius of curvature, corresponding cam angle

% alpha1max = 0; deltaalpha1max = 0; 
% initialize max angle for push travel, corresponding cam angle

% alpha2max = 0; deltaalpha2max = 0; 
% initialize max angle for return travel, corresponding cam angle

% x, y, xr, xr: theoretical and actual profile points, resp.

while 1
    
    e = 1/4*r0;
    s0 = sqrt(r0*r0 - e*e);
    rhomin = 1000;deltarhomin = 0;
    alpha1max = 0; deltaalpha1max = 0;
    alpha2max = 0; deltaalpha2max = 0;

    for n = 1:N

        %push travel
        rdeg = deg(n)*dr; % deg. to rad.
        if  rdeg <= delta1
            T = rdeg/delta1;
            s(n) = h*(10*T^3 - 15*T^4 + 6*T^5);
           
            v(n) = 30*h*T^2*(1 - 2*T + T^2)/delta1; %v=ds/ddelta
            ds = v(n);
            
            a(n) = 60*h*T*(1 - 3*T + 2*T^2)/delta1^2;
            dds = a(n);
            
            %pressure angle for push travel
            alpha1 = atan((abs(ds) - e)/(s(n) + s0));
            
                %select max pressure angle for push travel
                if  alpha1 > alpha1max
                    alpha1max = alpha1;
                    deltaalpha1max = rdeg;
                end
                
        %select min radius of curvature        
        %1st-derivative w.r.t. delta
        dx(n) = (ds - e)*sin(rdeg) + (s0 + s(n))*cos(rdeg);
        dy(n) = (ds - e)*cos(rdeg) - (s0 + s(n))*sin(rdeg);
        
        %2nd-derivative w.r.t. delta
        ddx(n) = (dds - s0 - s(n))*sin(rdeg) + (2*ds - e)*cos(rdeg);
        ddy(n) = (dds - s0 - s(n))*cos(rdeg) - (2*ds - e)*sin(rdeg);
        
        %radius of curvature
        rho = ((ddx(n)^2 + ddy(n)^2)^(3/2))/(dx(n)*ddy(n) - dy(n)*ddx(n)); 
        
        if rho < rhomin && rho > 0
            rhomin = rho;
            deltarhomin = rdeg;
        end
        
        %farthest dwell angle
        elseif rdeg > delta1 && rdeg <= delta12
                s(n) = h; v(n) = 0;
                ds = v(n);

        %return travel
        elseif rdeg > delta12 && rdeg <= delta13
               degback = rdeg - delta12;
               s(n) = 0.5*h*(1 + cos(pi*degback/delta3));
               
               v(n) = -0.5*pi*h*sin(pi*degback/delta3)/(delta3); 
               ds = v(n);
               
               a(n) = -0.5*pi^2*h*cos(pi*degback/delta3)/(delta3);
               dds = a(n);

            %pressure angle for return travel
            alpha2 = atan((abs(ds) + e)/(s(n) + s0));
            
                %select max pressure angle for return travel
                if alpha2 > alpha2max
                   alpha2max = alpha2;
                   deltaalpha2max = rdeg;
                end

        %nearest dwell angle
        elseif rdeg > delta13 && rdeg <= delta14
               s(n) = 0; v(n) = 0;
               ds = v(n);
        end
        
        %select min radius of curvature
        %1st-derivative w.r.t. delta
        dx(n) = (ds - e)*sin(rdeg) + (s0 + s(n))*cos(rdeg);
        dy(n) = (ds - e)*cos(rdeg) - (s0 + s(n))*sin(rdeg);
        
        %2nd-derivative w.r.t. delta
        ddx(n) = (dds - s0 - s(n))*sin(rdeg) + (2*ds - e)*cos(rdeg);
        ddy(n) = (dds - s0 - s(n))*cos(rdeg) - (2*ds - e)*sin(rdeg);
        
        %radius of curvature
        rho = ((ddx(n)^2 + ddy(n)^2)^(3/2))/(dx(n)*ddy(n) - dy(n)*ddx(n)); 
        
        if rho < rhomin && rho > 0
            rhomin = rho;
            deltarhomin = rdeg;
        end

        %------calculating cam profile curve---------------------
        %theoretical cam profile

        x(n) = (s0 + s(n))*sin(rdeg) + e*cos(rdeg);
        y(n) = (s0 + s(n))*cos(rdeg) - e*sin(rdeg);

        %real cam profile
        stheta = dx(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
        ctheta = -dy(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
        
        %inner envelope profile
        xr(n) = x(n) - rr*ctheta;
        yr(n) = y(n) - rr*stheta;
        
    end %w.r.t. for

        %If parameters are not suitable, adjust the radius of the base circle
        if alpha1max > alpha1allow || alpha2max > alpha2allow || rhomin - e < rhominallow
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

fprintf('min radius of curvature, corresponding cam angle\n');
fprintf('%6.4f %6.4f\n',rhomin,deltarhomin*rd);

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
title('Design of offset disk cam');
xlabel('x/mm');
ylabel('y/mm');
plot(x,y, 'r-');
ct = linspace(0, 2*pi);
plot(r0*cos(ct), r0*sin(ct), 'g-'); %base circle
plot(e*cos(ct), e*sin(ct), 'c-'); %offset circle

%------------actual profile------------------------------------
plot(xr,yr, 'b-');

%------------motion diagram------------------------------------
figure(2)
plot(deg,s,'r-');
title('Motion diagram');
xlabel('\phi/\circ');
ylabel('s/mm')
%--------------------------------------------------------------

