%% 448 Final Project, 3-bus system

close all, clear, clc;

% given parameters of the 3-bus system

%Algorihtim for the Zbase
V_base = (345 * 10^3);
VA_base = (100*10^6);
Z_base = (V_base)^2/VA_base;

% Maximum acceptable tolerance
epsilon = 1e-5;

% Given parameters of each bus
V1 = 1;
theta1 = 0;
V2 = 1.05;
P2sp = 2;
P3sp = 5;
Q = 1;
c = [P2sp ; P3sp ; Q ];

% initial guesses
theta2 = 0;
theta3 = 0;
V3 = 1;
x = [theta2 ; theta3 ; V3];

%% Calculating circuit parameters

% define distance
D1 = 150;
D2 = 200;
% define the amount of busses
% defining Z_12 before normalization by the base
Z_12prebase = D1*(0.037 + (j*0.376));
Z_13prebase = D2*(0.037 + (j*0.376));
%defining the series impedence per unit for each Z_mk
Z_12 = Z_12prebase/Z_base;
Z_13 = Z_13prebase/Z_base;

% Busses are NOT equidistant !

% mutual impedances
Z_21 = Z_12; % bus connections 13 and 23  have same impedance
Z_23 = Z_12;
Z_32 = Z_23;

Z_31 = Z_13;

% self-impedances, ignore
Z_11 = 0;
Z_22 = 0;
Z_33 = 0;

% Shunt susceptence will be equal to about zero
Y_kg = 0;


% Self admittances
Y_11 = (1/Z_12) + (1/Z_13);
Y_22 = (1/Z_21) + (1/Z_23);
Y_33 = (1/Z_31) + (1/Z_32); % equal to Y_11

% Mutual Admittances
% Admittance between bus 1 and bus 2
Y_12 = -1/Z_12;
Y_21 = Y_12;

% Admittance between bus 2 and bus 3
Y_23 = -1/Z_23;
Y_32 = Y_23;

% Admittance between bus 1 and bus 3

Y_13 = -1/Z_13;
Y_31 = Y_13;

%% Admittance Matrix

Y_matrix = [Y_11 Y_12 Y_13; Y_21 Y_22 Y_23; Y_31 Y_32 Y_33];

G = real(Y_matrix); % Conductance

B = imag(Y_matrix); % Susceptance

%% Power matrix
P2 = @(x) G(2,2)*(V2^2)+V2*((V1*G(2,1)*cos(x(1)-theta1)+V1*B(2,1)*sin(x(1)-theta1))  +  (x(3)*G(2,3)*cos(x(1)-x(2))+x(3)*B(2,3)*sin(x(1)-x(2))));
P3 = @(x) G(3,3)*(x(3)^2)+x(3)*((V1*G(3,1)*cos(x(2)-theta1)+V1*B(3,1)*sin(x(2)-theta1))  +  (V2*G(3,2)*cos(x(2)-x(1))+V2*B(3,2)*sin(x(2)-x(1))));
Q3 = @(x) -B(3,3)*(x(3)^2)+x(3)*((V1*G(3,1)*sin(x(2)-theta1)-V1*B(3,1)*cos(x(2)-theta1))  +  (V2*G(3,2)*sin(x(2)-x(1))-V2*B(3,2)*cos(x(2)-x(1))));

funcMatrix = @(x) [P2(x); P3(x); Q3(x)];

%% Jacobian Matrix

%partial derivative of P2 in respect to theta2, theta3, and V3
J_11 = @(x) V2*V1*(-G(2,1) * sin(x(1) - theta1) + B(2,1) * cos(x(1) - theta1))   +   V2*x(2)*(-G(2,3)*sin(x(1) - x(3))+B(2,3)* cos(x(1) - x(2)));
J_12 = @(x) V2*x(3)*(G(2,3) * sin(x(1) - x(2))- B(2,3)*cos(x(1) - x(2)));
J_13 = @(x) V2*(G(2,3)*cos(x(1) - x(2))+ B(2,3)*sin(x(1) - x(2)));
%partial derivative of P3 in respect to theta2, theta3, V3
J_21 = @(x) x(3)*V2*(G(3,2)*sin(x(1) - x(2)) - B(3,2)*cos(x(1) - x(2)));
J_22 = @(x) x(3)*V1*(-G(3,1) * sin(theta1 - x(1)) + B(3,1)* cos(theta1 - x(2))) + x(3)*V2*(-G(3,2)*sin(x(1) - x(2))+ B(3,2)*cos(x(1) - x(2)));
J_23 = @(x) 2 * G(3,3) * x(3) + V1*(G(3,1) * cos(theta1 - x(2)) + B(3,1) * sin(theta1 - x(2))) + V2*(G(3,2)*cos(x(1) - x(2)) + B(3,2)*sin(x(1) - x(2)));
%partial derivative of Q3 in respect to theta2, theta3, V3
J_31 = @(x) x(3)*V2*(-G(3,2)*cos(x(1) - x(2)))- B(3,2)*sin(x(1) - x(2));
J_32 = @(x) x(3)*V1*(G(3,1) * cos(theta1 - x(2))+ B(3,1) * sin(theta1 - x(2))) + x(3)*V2*(G(3,2) * cos(x(1) - x(2)) + B(3,2) * sin(x(1) - x(2)));
J_33 = @(x) -2* B(3,3) * x(3) + V1*(G(3,1)*sin(theta1 - x(2)) - B(3,1) * cos(theta1 - x(2)))  +   V2*(G(3,2)*sin(x(1) - x(2)) - B(3,2)*cos(x(1) - x(2)));
%setting up the jacobian matrix

%% Iteration process

n = 0; % iterations
x_plot = x;

while norm((c-funcMatrix(x))) > epsilon % error matrix is normalized to a real value. 
    n = n+1;
    
    JacobMatrix = [J_11(x) J_12(x) J_13(x); J_21(x) J_22(x) J_23(x); J_31(x) J_32(x) J_33(x)];
    delta_x = inv(JacobMatrix)*(c-funcMatrix(x)); % delta to approach tolerance
    
    % limit delta phase angles between 0 and 2pi radians
    delta_x(1) = wrapTo2Pi(delta_x(1)); 
    delta_x(2) = wrapTo2Pi(delta_x(2));

    x = x + delta_x;

    % plot containing each iteration
    x_plot = [x_plot, x]
end
%% print result for iterations, guesses, etc.

% limit final x phase angles between 0 to 2pi

    x(1) = wrapTo2Pi(x(1)); 
    x(2) = wrapTo2Pi(x(2));

% print results
error = norm(c - funcMatrix(x))
n
x

% plot of unknowns versus iterations
stem(x_plot(1,:));
title("Plots of unknown variables and their change across each iteration")
ylabel("theta 2 (radians)");
xlabel("# of iterations");

figure;
stem(x_plot(2,:));
title("Plots of unknown variables and their change across each iteration")
ylabel("theta 3 (radians)");
xlabel("# of iterations");

figure;
stem(x_plot(3,:));
title("Plots of unknown variables and their change across each iteration")
ylabel("V3 (per-unit)");
xlabel("# of iterations");











