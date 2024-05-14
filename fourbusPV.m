%448 Final Project
clear;clc;
 
%Algorihtim for the Zbase
V1 = 1;
V_base = (345 * 10^3);
VA_base = (100*10^6);
Z_base = (V_base)^2/VA_base;
theta1 = 0;
V2 = 1.05;
V4 = 0.95;

%defining givens we have 
esp = 1e-3;
ne=1.1*eps;     % So I can enter my while loop
P2sp = 2;
P3sp = 5;
P4sp = 1.5;
Q = 1;
%vector to store power values
c = [P2sp ; P3sp ; Q ; P4sp];
%Vector to store angles and 
theta2 = 0;
theta3 = 0;
theta4 = 0;
V3 = 1;
x = [theta2 ; theta3 ; V3 ; theta4];


% define distance 
D = 150;
D2 = sqrt(2) * 150;

% define the amount of busses

n = 0;

% defining Z_12 before normalization by the base
Z_12prebase = D*(0.037 + (j*0.376));

Z_12 = Z_12prebase/Z_base;
Z_23 = Z_12;
Z_34 = Z_12;
Z_14 = Z_12;

Z_13 = 2* Z_12;
Z_24 = Z_13;

Z_42 = Z_24;

Z_21 = Z_12;
Z_31 = Z_13;
Z_41 = Z_14;

Z_32 = Z_23;

Z_43 = Z_34;




Y_kg = 0;
summyIn = 2.*(1./Z_12);
Y_12 = (-1)./Z_12;
Y_13 = (-1)./(Z_13);
Y_14 = (-1)./(Z_14);

Y_21 = (-1)./(Z_21);
Y_23 = (-1)./(Z_23);
Y_24 = (-1)./(Z_24);

Y_31 = (-1)./(Z_13);
Y_32 = (-1)./(Z_23);
Y_34 = (-1)./(Z_34);

Y_41 = (-1)./(Z_14);
Y_42 = (-1)./(Z_24);
Y_43 = (-1)./(Z_43);

Y_11 = (1/Z_12) + (1/Z_13) + (1/Z_14);
Y_22 = (1/Z_12) + (1/Z_23) + (1/Z_24); 
Y_33 = (1/Z_13) + (1/Z_23) + (1/Z_34);
Y_44 = (1/Z_14) + (1/Z_24) + (1/Z_34);

Y_kk = Y_kg + summyIn;

Y_admit = [Y_11 Y_12 Y_13 Y_14; Y_21 Y_22 Y_23 Y_24; Y_31 Y_32 Y_33 Y_34; Y_41 Y_42 Y_43 Y_44]


%m = 0;
G = real(Y_admit);
B = imag(Y_admit);

% Define the Jacobian matrix elements


% Partial derivatives of P2 with respect to theta2, theta3, V3, and theta4
j1= @(x) -V2 * (V1 * cos(theta1 - theta2) * B(2, 1) - V2 * cos(theta1 - theta2) * B(2, 3) - V3 * cos(theta2 - theta4) * B(2, 4) + V1 * sin(theta1 - theta2) * G(2, 1) + V2 * sin(theta2 - theta3) * G(2, 3) + V3 * sin(theta2 - theta4) * G(2, 4));
j2 =@(x) V2^2 * sin(theta2 - theta3) * G(2, 3);
j3 =@(x) V2 * (sin(theta2 - theta4) * G(2, 4) - cos(theta2 - theta4) * B(2, 4));
j4 =@(x) -V2 * (V3 * cos(theta2 - theta4) * B(2, 4) - V3 * sin(theta2 - theta4) * G(2, 4));

% Partial derivatives of P3 with respect to theta2, theta3, V3, and theta4
j5 =@(x) V3^2 * (G(3, 2) * cos(theta3 - theta2) + B(3, 2) * sin(theta3 - theta2)) + V3 * (V1 * cos(theta1 - theta2) * G(3, 1) - V1 * sin(theta1 - theta2) * B(3, 1) + V2 * cos(theta2 - theta3) * G(3, 2) - V2 * sin(theta2 - theta3) * B(3, 2) + V3 * cos(theta3 - theta4) * G(3, 4) + V3 * sin(theta3 - theta4) * B(3, 4));
j6 =@(x) V3^2 * (-G(3, 2) * cos(theta3 - theta2) + B(3, 2) * sin(theta3 - theta2)) + V3 * (-V1 * cos(theta1 - theta2) * G(3, 1) + V1 * sin(theta1 - theta2) * B(3, 1) - V2 * cos(theta2 - theta3) * G(3, 2) + V2 * sin(theta2 - theta3) * B(3, 2) + V3 * cos(theta3 - theta4) * G(3, 4) + V3 * sin(theta3 - theta4) * B(3, 4));
j7 =@(x) 2 * V3 * (G(3, 3) * V3 + B(3, 3) * V3) + (V1 * cos(theta1 - theta2) * G(3, 1) - V1 * sin(theta1 - theta2) * B(3, 1) + V2 * cos(theta2 - theta3) * G(3, 2) - V2 * sin(theta2 - theta3) * B(3, 2) + V3 * cos(theta3 - theta4) * G(3, 4) + V3 * sin(theta3 - theta4) * B(3, 4));
j8 =@(x) V3^2 * (G(3, 4) * cos(theta3 - theta4) + B(3, 4) * sin(theta3 - theta4)) + V3 * (G(3, 3) * V3 + B(3, 3) * V3);

% Partial derivatives of Q with respect to theta2, theta3, V3, and theta4
j9 =@(x) V3^2 * (-B(3, 2) * cos(theta3 - theta2) - G(3, 2) * sin(theta3 - theta2)) + V3 * (-V1 * sin(theta1 - theta2) * G(3, 1) - V1 * cos(theta1 - theta2) * B(3, 1) - V2 * sin(theta2 - theta3) * G(3, 2) - V2 * cos(theta2 - theta3) * B(3, 2) - V3 * sin(theta3 - theta4) * G(3, 4) + V3 * cos(theta3 - theta4) * B(3, 4));
j10=@(x) V3^2 * (B(3, 2) * cos(theta3 - theta2) + G(3, 2) * sin(theta3 - theta2)) + V3 * (V1 * sin(theta1 - theta2) * G(3, 1) - V1 * cos(theta1 - theta2) * B(3, 1) + V2 * sin(theta2 - theta3) * G(3, 2) - V2 * cos(theta2 - theta3) * B(3, 2) - V3 * sin(theta3 - theta4) * G(3, 4) + V3 * cos(theta3 - theta4) * B(3, 4));
j11 =@(x) 2 * V3 * (-G(3, 3) * V3 + B(3, 3) * V3) + (-V1 * sin(theta1 - theta2) * G(3, 1) - V1 * cos(theta1 - theta2) * B(3, 1) - V2 * sin(theta2 - theta3) * G(3, 2) - V2 * cos(theta2 - theta3) * B(3, 2) - V3 * sin(theta3 - theta4) * G(3, 4) + V3 * cos(theta3 - theta4) * B(3, 4));
j12 =@(x) V3^2 * (-G(3, 4) * cos(theta3 - theta4) - B(3, 4) * sin(theta3 - theta4)) + V3 * (-G(3, 3) * V3 + B(3, 3) * V3);

% Partial derivatives of P4 with respect to theta2, theta3, V3, and theta4
j13 =@(x) -V4 * (V1 * cos(theta1 - theta4) * B(4, 1) - V1 * sin(theta1 - theta4) * G(4, 1) + V2 * cos(theta2 - theta4) * G(4, 2) - V2 * sin(theta2 - theta4) * B(4, 2) + V3 * cos(theta3 - theta4) * G(4, 3) - V3 * sin(theta3 - theta4) * B(4, 3));
j14 =@(x) -V4 * (V3 * cos(theta3 - theta4) * G(4, 3) - V3 * sin(theta3 - theta4) * B(4, 3));
j15 =@(x) V4 * (sin(theta3 - theta4) * B(4, 3) - cos(theta3 - theta4) * G(4, 3));
j16 =@(x) -V4 * (V3 * cos(theta3 - theta4) * B(4, 3) - V3 * sin(theta3 - theta4) * G(4, 3));




P2 = @(x) ((G(2,2))*V2^2) + V2*(V1*(G(2,1))*cos(x(1) - theta1)+ V1*(B(2,1))*sin(x(1) - theta1) + V3*(G(2,3))*cos(x(2) - x(1)) + V3*(B(2,3))*sin(x(2) - x(1)) + V4*(G(2,4))*cos(x(4) - x(1)) + ...
V4*(B(2,4))* sin(x(4) - x(1)));

P3 = @(x) ((G(3,3))*x(3)^2) + x(3)*(V1*(G(3,1))*cos(x(2) - theta1)+ V1*(B(3,1))*sin(x(2) - theta1) + V2*(G(3,2))*cos(x(2) - x(1)) + V2*(B(3,2))*sin(x(2) - x(1)) + V4*(G(3,4))*cos(x(4) - x(2)) + ...
V4*(B(3,4))* sin(x(4) - x(2)));

Q3 = @(x) ((B(3,3))*x(3)^2) + x(3)*(V1*(G(3,1))*sin(x(2) - theta1)- V1*(B(3,1))*cos(x(2) - theta1) + V2*(G(3,2))*sin(x(2) - x(1)) - V2*(B(3,2))*cos(x(2) - x(1)) + V4*(G(3,4))*sin(x(4) - x(2)) - ...
V4*(B(3,4))* cos(x(4) - x(2))) ;
 
P4 = @(x) (G(4,4) * V4^2) + V4*(V1*(G(4,1))*cos(x(4) - theta1) +V1*(B(4,1))*sin(x(4) - theta1)+ V2*(G(4,2))* cos(x(4)-x(1)) + V2*(B(4,1))*sin(x(4)-x(1)) + x(3)*(G(4,2))*cos(x(4)-x(2)) + x(3)*(B(4,2))*sin(x(4)-x(2))) ;
 
funcMatrix2 = @(x) [P2(x); P3(x); Q3(x); P4(x)];


while norm(c-funcMatrix2(x)) > esp
JacobianMatrix = [j1(x) j2(x) j3(x) j4(x); j5(x) j6(x) j7(x) j8(x); j9(x) j10(x) j11(x) j12(x); j13(x) j14(x) j15(x) j16(x)];

delta_x = inv(JacobianMatrix)*(c-funcMatrix2(x)) % delta to approach tolerance
delta_x(1) = wrapTo2Pi(delta_x(1));
delta_x(2) = wrapTo2Pi(delta_x(2));
delta_x(4) = wrapTo2Pi(delta_x(4));
x = x + delta_x;

ne = norm(esp);
n = n+1
end
error = norm(c - funcMatrix2(x));
n
x
c
ne
