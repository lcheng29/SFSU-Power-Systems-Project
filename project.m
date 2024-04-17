% % First example on NR method
% 
% % Tolerance epsilon
% 
% eps = 1e-8;
% 
% % initial estimate
% x = 0.1;
% error = 1;
% n = 0;
% C = 3.993e-4;
% F = -1*x^3 + 0.0165*x^2;
% df = -3*x^2 + 0.033*x;
% 
% % C - f(x) = 0
% while abs(error) > eps
%     n = n + 1;
%     x_adj = (C - F)/(df);
%     x = x + x_adj;
%     error = x_adj;
% end
% 
% x
% x_adj
% error
% n

%% example power network with 2 buses

clear
clc
% equations
% 5 - x1^2 - x2 = 0 
% 3 + x1 - x2^2 = 0

x0 = [5 -5].'; % initial guess
eps = 1e-6; % tolerance
x = x0;
c = [5 3].';

f1 = @(x) x(1)^2 + x(2);
f2 = @(x) -x(1) + x(2)^2;

f = @(x) [f1(x) f2(x)].';

df11 = @(x) 2*x(1); df12 = @(x) 1;
df21 = @(x) -1; df22 = @(x) 2*x(2);

n = 0; % iterations

while norm((c-f(x)), 2) > eps

    J = [df11(x) df12(x); 
            df21(x) df22(x)];

    delta = inv(J) * (c-f(x));
    x = x + delta;
    n = n+1;
end

n, x, error = c - f(x)

%% Part 1 of Project: Three buses (slack, PV, PQ)
% slack bus
v1 = 1 % 1 per-unit (pu)
theta1 = 0

% PV bus
v2 = 1.05 % 1.05 pu
p2_sp = 2 % 2 pu

% PQ bus
p3_sp = 5 % 5 pu
q3_sp = 1 % 1 pu

%% system parameters (impedance, admittance, conductance, susceptance, etc) all in per unit (pu)

% susceptances (B)
b12 = 0.8034
b21 = b12

b13 = 1.0712
b31 = b13

b23 = 0.8034
b32 = b23

b33 = 0

% impedances (Z)

z12 = 0.0047 + 0.0474j
z13 = 0.0062 + 0.0632j
z23 = 0.0047 + 0.0474j

% admittances (Y)

y12 = 1/z12
y13 = 1/z13
y23 = 1/z23

% conductances (G)
g12 = real(y12)
g21 = g12

g13 = real(y23)
g31 = g13

g23 = real(y23)
g32 = g23

g33 = 0
%% initial guess

x0 = [0 0 1].'; % x0 = [theta2 @ (0), theta3 @ (0), v3 @ (0) = 1]

p_2 = 0 + v2*(v1*g21*cos(theta21) + v1*b21*sin(theta21) + v3*g23*cos(theta23) + v3*b23*sin(theta23)) 

p_3 = g33*v3^2 + v3*(v1*g31*cos(theta31) + v1*b31*sin(theta31) + v2*g32*cos(theta32) + v2*b32*sin(theta32)) 

q_3 = -b33*v3^2 + v3*(v1*g31*sin(theta31) - v1*b31*cos(theta31) + v2*g32*sin(theta32) - v2*b32*cos(theta32))

%% Jacobian matrix

d_p2_th2 = v2*v1 * (-g21 * sin(theta21))
d_p2_th3 = v2*v3 * (g23 * sin(theta23) - b23 * cos(theta23))
d_p2_v3  = v2 * (g23 * cos(theta23) + b23 * sin(theta23))

d_p3_th2 = v3 * v2 * (g32*sin(theta32) - b32*cos(theta32))
d_p3_th3 = v3 * v1 * (-g31*sin(theta31) + b31*cos(theta31)) + v3*v2(-g32*sin(theta32) + b32*cos(theta32))
d_p3_v3  = 2*g33*v3 + v1*(g31*cos(theta31) + b31*sin(theta31)) + v2*(g32*cos(theta32) + b32*sin(theta32))

d_q3_th2 = v3*v2 * (-g32*cos(theta32) - b32*sin(theta32))
d_q3_th3 = v3*v1 * (g31*cos(theta31) + b31*sin(theta31)) + v3*v2 * (g32*cos(theta32) + b32*sin(theta32))
d_q3_v3  = -2*b33*v3 + v1 * (g31*sin(theta31) - b31*cos(theta31)) + v2 * (g32*sin(theta32) - b32*cos(theta32))

% full Jacobian matrix
J = [d_p2_th2, d_p2_th3, d_p2_v3; d_p3_th2, d_p3_th3, d_p3_v3; d_q3_th2, d_q3_th3, d_q3_v3]












