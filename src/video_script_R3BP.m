clc; close all; clear;

% First-Order Linear
%{
omega = 2.0152105515;
time_vec = 0:0.01:1/omega;
scaling_factor = 1.495978714e8;
w_p = 2.086453455;
w_v = 2.0152105515;
Ax = 206000/scaling_factor;
Ay = 665000/scaling_factor;
Az = 110000/scaling_factor;

t = 0:0.05:100;%2*pi;%100;
angle1 = 0:0.01:2*pi;
angle2 = 0:0.01:2*pi;

xpoint = 0.932385;
ypoint = 0;
zpoint = 0;

mu = -9.537e-4;
xsun = -mu;
ysun = 0;
zsun = 0;
xearth = 1 + xsun;
yearth = 0;
zearth = 0;

theta = 0 : 0.01 : 2*pi;
eradius = 6378/scaling_factor;
sradius = 696000/scaling_factor;

phi = 0;
m = 1;
psi = m*(pi/2) + phi;

% PLOTS 1
%{
x1 = -Ax*cos(angle1 + phi) + xpoint;
y1 = Ay*sin(angle1 + phi) + ypoint;
z1 = Az*sin(angle2 + 1*(pi/2) + phi) + zpoint;
x3 = -Ax*cos(angle1 + phi) + xpoint;
y3 = Ay*sin(angle1 + phi) + ypoint;
z3 = Az*sin(angle2 + 3*(pi/2) + phi) + zpoint;

subplot(3,1,1);
hold on;
plot(xpoint, ypoint, 'b*', 'DisplayName', 'L1 Lagrange Point')
plot(xearth, yearth, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'Earth')
plot(xsun, ysun, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'DisplayName', 'Sun')
plot(x1, y1, 'r-', 'DisplayName', 'M = 1')
plot(x3, y3, 'k-', 'DisplayName', 'M = 3')
xlabel('X')
ylabel('Y')
grid on
grid(gca, 'minor')
lgd = legend();
lgd.FontSize = 20;
hold off;

subplot(3,1,2);
hold on;
plot(xearth, zearth, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')
plot(xsun, zsun, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y')
plot(xpoint, zpoint, 'b*')
plot(x1, z1, 'r-', 'DisplayName', 'M = 1')
plot(x3, z3, 'k-', 'DisplayName', 'M = 3')
xlabel('X') 
ylabel('Z')
grid on
grid(gca, 'minor')
hold off;

subplot(3,1,3);
hold on;
plot(yearth, zearth, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')
plot(ysun, zsun, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y')
plot(ypoint, zpoint, 'b*')
plot(y1, z1, 'r-', 'DisplayName', 'M = 1')
plot(y3, z3, 'k-', 'DisplayName', 'M = 3')
xlabel('Y') 
ylabel('Z')
grid on
grid(gca, 'minor')
hold off;
%}

% PLOTS 2
%{
m = 1;
phi = 0;
psi = m*(pi/2) + phi;
x1 = -Ax*cos(w_p*t + phi) + xpoint;
y1 =  Ay*sin(w_p*t + phi) + ypoint;
z1 =  Az*sin(w_v*t + psi) + zpoint;

hold on;
[X,Y,Z] = sphere;
X2 = X * eradius;
Y2 = Y * eradius;
Z2 = Z * eradius;
XS = X * sradius;
YS = Y * sradius;
ZS = Z * sradius;
surf(X2 + xearth, Y2, Z2)
surf(XS + xsun, YS, ZS)
plot3(xpoint, ypoint, zpoint, 'b*')
plot3(x1, y1, z1, 'r', 'DisplayName', 'M = 1')
zlim([-0.01 0.01])
ylim([-0.01 0.01])
grid on;
xlabel('X') 
ylabel('Y') 
zlabel('Z') 
hold off;
%}

% PLOTS 3
%{
m = 1;
phi = 0;
psi = m*(pi/2) + phi;
x1 = -Ax*cos(w_p*t + phi) + xpoint;
y1 =  Ay*sin(w_p*t + phi) + ypoint;
z1 =  Az*sin(w_v*t + psi) + zpoint;

hold on;
plot3(xpoint, ypoint, zpoint, 'b*')
plot3(xearth, yearth, zearth, 'bo', 'MarkerSize', 5)
plot3(xsun, ysun, zsun, 'ro', 'MarkerSize', 10)
plot3(x1, y1, z1, 'r', 'DisplayName', 'M = 1')
% xlim([0.925 0.94])
grid on;
xlabel('X') 
ylabel('Y') 
zlabel('Z') 
hold off;
%}

% PLOTS 4
%{
m = 1;
phi = 0;
psi = m*(pi/2) + phi;
x1 = -Ax*cos(w_p*t + phi) + xpoint;
y1 =  Ay*sin(w_p*t + phi) + ypoint;
z1 =  Az*sin(w_v*t + psi) + zpoint;
subplot(3,1,1);
hold on;

plot(xpoint, ypoint, 'b*', 'DisplayName', 'L1 Lagrange Point')
plot(xearth, yearth, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'Earth')
plot(xsun, ysun, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'DisplayName', 'Sun')
plot(x1, y1, 'r-', 'DisplayName', 'M = 1')
xlabel('X')
ylabel('Y')
grid on
grid(gca, 'minor')
lgd = legend();
lgd.FontSize = 20;
hold off;

subplot(3,1,2);
hold on;
plot(xearth, zearth, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')
plot(xsun, zsun, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y')
plot(xpoint, zpoint, 'b*')
plot(x1, z1, 'r-', 'DisplayName', 'M = 1')
xlabel('X') 
ylabel('Z')
grid on
grid(gca, 'minor')
hold off;

subplot(3,1,3);
hold on;
plot(yearth, zearth, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')
plot(ysun, zsun, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'y')
plot(ypoint, zpoint, 'b*')
plot(y1, z1, 'r-', 'DisplayName', 'M = 1')
xlabel('Y')
ylabel('Z')
grid on
grid(gca, 'minor')
hold off;
%}
%}



%%%%%%% Lagrange Point Information for Earth-Sun System %%%%%%%
x_L1 = 1;
x_L2 = 1;
%%%%%%% Lagrange Point Information for Earth-Sun System %%%%%%%


% Third-Order Linear

msun = 1.98847e30;
mearth = 5.9736e24;
MuScale = mearth/(msun + mearth);

gamma = 1.001090475e-2;
LScale = 1.495978714e8;
x = 0.932385*LScale;
x_ = (x - 1 + MuScale + gamma)/gamma;

T = 1.53536256e7; % L1 orbital period in seconds for halo orbits

VScale = 29.784;
TScale = 3.147e7;


k = 3.2292680962;
c2 = 4.0610735668;
c3 = 3.0200105081;
c4 = 3.0305378797;
lam = sqrt((c2 + sqrt(9*c2^2 - 8*c2))/2);

d1 = ((3*lam^2)/k)*(k*(6*lam^2 - 1) - 2*lam); 
d2 = ((8*lam^2)/k)*(k*(11*lam^2 - 1) - 2*lam);
a21 = 3*c3*(k^2 - 2)/(4*(1 + 2*c2)); 
a22 = 3*c3/(4*(1 + 2*c2));

a23 = -((3*c3*lam)/(4*k*d1))*(3*(k^3)*lam - 6*k*(k-lam) + 4); 
a24 = -((3*c3*lam)/(4*k*d1))*(3*k*lam + 2);

b21 = -3*c3*lam*(3*k*lam - 4)/(2*d1);
b22 = -3*c3*lam/d1;

a31 = (-9*lam/(4*d2))*(4*c3*(k*a23 - b21) + k*c4*(4 + k^2)) ...
    + ((9*lam^2 + 1 - c2)/(2*d2))*(3*c3*(2*a23 - k*b21) + c4*(2 + 3*k^2));

d21 = -c3/(2*lam^2);
d31 = (3/(64*lam^2))*(4*c3*a24 + c4);
d32 = (3/(64*lam^2))*(4*c3*(a23 - d21) + c4*(4 + k^2));

a32 = (-9*lam/(4*d2))*(4*c3*(3*k*a24 - b22) + k*c4) ... 
    - 3*((9*lam^2 + 1 - c2)/(2*d2))*(c3*(k*b22 - d21 - 2*a24) - c4);

b31 = (3/(8*d2))*8*lam*(3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k^2)) ... 
     + (3/(8*d2))*((9*lam^2 + 1 + 2*c2)*(4*c3*(k*a23 - b21)) + k*c4*(4 + k^2));

b32 = ((9*lam)/d2)*(c3*(k*b22 + d21 + -2*a24) - c4) ... 
    + 3*((9*lam^2 + 1 + 2*c2)/(8*d2))*(4*c3*(k*a24 - b22) + k*c4); 


scaling_factor = 1.495978714e8;

% xpoint = 0.932385;%*scaling_factor;
xpoint = 0.5;%*scaling_factor;
ypoint = 0;
zpoint = 0;
mu = 3.040423398444176e-6;
xsun = -mu;
ysun = 0;
zsun = 0;
xearth = 1 + xsun;
yearth = 0;
zearth = 0;

phi = 0;
w_p = 2.086453455;
tau = 0:0.01:pi;%100;
tau_1 = w_p*tau + phi;

eradius = 6378/scaling_factor;
sradius = 696000/scaling_factor;
Ax = 206000/scaling_factor;   % Ax is in km
Ay = 665000/scaling_factor;   % Ay is in km
Az = 110000/scaling_factor;   % Az is in km
m = 1;
dm = 2 - m;
x1 =    -Ax*cos(tau_1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau_1) + (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau_1) + a21*Ax^2 + a22*Az^2 + xpoint;
y1 =     Ay*sin(tau_1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau_1) + (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau_1) + ypoint;
z1 =  dm*Az*cos(tau_1) + dm*d21*Ax*Az*cos(2*tau_1 - 3) + dm*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau_1) + zpoint;

hold on;
plot3(xpoint, ypoint, zpoint, 'b*')
[X,Y,Z] = sphere;
X2 = X * eradius;
Y2 = Y * eradius;
Z2 = Z * eradius;
XS = X * sradius;
YS = Y * sradius;
ZS = Z * sradius;
surf(X2 + xearth, Y2, Z2)
alpha 0.3
surf(XS + xsun, YS, ZS, 'FaceColor', 'y', 'EdgeColor', 'k')
plot3(x1, y1, z1, 'r', 'DisplayName', 'M = 1')

% Draw XY plane
% (1, 1, 0)
% (-1, 1, 0)
% (-1, -1, 0)
% (1, -1, 0)

p = patch([xearth xsun xsun xearth], [max(y1) max(y1) min(y1) min(y1)], [0 0 0 0]);
% p = patch([max(x1) min(x1) min(x1) max(x1)], [max(y1) max(y1) min(y1) min(y1)], [0 0 0 0]);
p.FaceAlpha = 0.1;
hold off

grid on;
xlabel('X') 
ylabel('Y') 
zlabel('Z') 
hold off;
%}




% BACKUP
% clc; close all; clear all;
% 
% omega = 2.0152105515;
% 
% time_vec = 0:0.01:1/omega;
% 
% % [x, y, z] = sphere(40);
% 
% scaling_factor = 1.497610041e6;
% w_p = 2.086453455;
% w_v = 2.0152105515;
% A_x = 206000/scaling_factor;
% A_z = 110000/scaling_factor;
% K = 3.2292680962;
% 
% m = 1;
% phi = 0;
% psi = m*(pi/2) + phi;
% 
% angle1 = 0:0.01:2*pi;
% angle2 = 0:0.01:2*pi;
% 
% xpoint = 0.932385;
% ypoint = 0;
% zpoint = 0;
% 
% mu = -9.537e-4;
% xsun = -mu;
% ysun = 0;
% zsun = 0;
% 
% xearth = 1 + xsun;
% yearth = 0;
% zearth = 0;
% 
% % xr = -A_x*cos(angle1) + xpoint;
% % yr = K*A_x*sin(angle1) + ypoint;
% % zr = A_z*sin(angle2) + zpoint;
% 
% x1 = -A_x*cos(angle1 + phi) + xpoint;
% y1 = K*A_x*sin(angle1 + phi) + ypoint;
% z1 = A_z*sin(angle2 + 1*(pi/2) + phi) + zpoint;
% 
% x3 = -A_x*cos(angle1 + phi) + xpoint;
% y3 = K*A_x*sin(angle1 + phi) + ypoint;
% z3 = A_z*sin(angle2 + 3*(pi/2) + phi) + zpoint;
% 
% gamma = 1.001090475e-2;
% % xbar_r = (xr - 1 + mu + gamma)/gamma;
% % ybar_r = yr/gamma;
% % zbar_r = zr/gamma;
% 
% xbar_1 = (x1 - 1 + mu + gamma)/gamma;
% ybar_1 = y1/gamma;
% zbar_1 = z1/gamma;
% 
% xbar_3 = (x3 - 1 + mu + gamma)/gamma;
% ybar_3 = y3/gamma;
% zbar_3 = z3/gamma;
% 
% % plot3(xr, yr, zr, 'b-*', 'DisplayName', 'Reference')
% % hold on;
% % plot3(x1, y1, z1, 'r-*', 'DisplayName', 'M = 1')
% % plot3(x3, y3, z3, 'k-*', 'DisplayName', 'M = 3')
% % hold off;
% % legend()
% % xlabel('X') 
% % ylabel('Y') 
% % zlabel('Z') 
% 
% 
% subplot(3,1,1);
% hold on;
% plot(xpoint, ypoint, 'b*')
% plot(xearth, yearth, 'bo', 'MarkerSize', 5)
% plot(xsun, ysun, 'ro', 'MarkerSize', 10)
% % plot(xr, yr, 'b-', 'DisplayName', 'Reference')
% % plot(x1, y1, 'r-', 'DisplayName', 'M = 1')
% % plot(x3, y3, 'k-', 'DisplayName', 'M = 3')
% plot(xbar_1, ybar_1, 'r-', 'DisplayName', 'M = 1')
% plot(xbar_3, ybar_3, 'k-', 'DisplayName', 'M = 3')
% xlabel('X') 
% ylabel('Y')
% % xlim([0.4, 1.4])
% legend()
% hold off;
% 
% subplot(3,1,2);
% hold on;
% plot(xearth, zearth, 'bo', 'MarkerSize', 5)
% plot(xsun, zsun, 'ro', 'MarkerSize', 10)
% plot(xpoint, zpoint, 'b*')
% % plot(xr, zr, 'b-')
% % plot(x1, z1, 'r-')
% % plot(x3, z3, 'k-')
% plot(xbar_1, zbar_1, 'r-', 'DisplayName', 'M = 1')
% plot(xbar_3, zbar_3, 'k-', 'DisplayName', 'M = 3')
% xlabel('X') 
% ylabel('Z') 
% hold off;
% 
% subplot(3,1,3);
% hold on;
% plot(yearth, zearth, 'bo', 'MarkerSize', 5)
% plot(ysun, zsun, 'ro', 'MarkerSize', 10)
% plot(ypoint, zpoint, 'b*')
% % plot(yr, zr, 'b-')
% % plot(y1, z1, 'r-')
% % plot(y3, z3, 'k-')
% plot(ybar_1, zbar_1, 'r-', 'DisplayName', 'M = 1')
% plot(ybar_3, zbar_3, 'k-', 'DisplayName', 'M = 3')
% xlabel('Y') 
% ylabel('Z') 
% hold off;
% 
% % BACKUP
% % x = -A_x*cos(time_vec.*w_p + phi);
% % y = K*A_x*sin(time_vec.*w_p + phi);
% % z = A_z*sin(time_vec.*w_v + psi);
% 
% % plot(y, z, 'b')
% % hold on;
% % plot(y_, z_, 'r')
% % hold off;
% 
