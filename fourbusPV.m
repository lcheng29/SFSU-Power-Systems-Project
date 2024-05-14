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
esp = 1e-5;
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

% define the amount of busses

n = 4;

% defining Z_12 before normalization by the base
Z_12prebase = D*(0.037 + (j*0.376));

Z_12 = Z_12prebase/Z_base;

Y_kg = 0;
summyIn = 2.*(1./Z_12);
Y_km = (-1)./Z_12;
Y_kk = Y_kg + summyIn;

Y_admit = [Y_kk Y_km Y_km Y_km; Y_km Y_kk Y_km Y_km; Y_km Y_km Y_kk Y_km; Y_km Y_km Y_km Y_kk];

%imaginary portion of Y_km for susceptence
G_km = real(Y_km);


%real portion of Y_km for conductance
B_km = imag(Y_km);

% real portion of Y_kk for susceptence
G_kk = real(Y_kk);

% real portion of Y_kk for conductance
B_kk = imag(Y_kk);

%m = 0;
G = real(Y_admit);
B = imag(Y_admit);



j1 =@(x) -V2*(V1*cos(theta1 - x(1))*(B_km) - V2*cos(x(1) - x(2))*(B_km) - V3*cos(x(1) - x(2))*(B_km) + V1*sin(theta1 - x(1))*(G_km) + V2*sin(x(1) - x(2))*(G_km) + V3*sin(x(1) - x(4))*(G_km));
 
 
j2 =@(x) -V2*(V2*cos(x(1) - x(2))*(B_km) - V2*sin(x(1) - x(2))*(G_km));
 
 
j3 =@(x) V2*(cos(x(1) - x(4))*(G_km) + sin(x(1) - x(4))*(B_km));
 
 
j4 =@(x) -V2*(V3*cos(x(1) - x(4))*(B_km) - V3*sin(x(1) - x(4))*(G_km));
 
 
j5 =@(x) V3*(V2*cos(x(1) - x(2))*(B_km) + V2*sin(x(1) - x(2))*(G_km));
 
 
j6 =@(x) -V3*(V1*cos(theta1 - x(2))*(B_km) + V2*cos(x(1) - x(2))*(B_km) - V3*cos(x(2) - x(4))*(B_km) + V1*sin(theta1 - x(2))*(G_km) + V2*sin(x(1) - x(2))*(G_km) + V3*sin(x(2) - x(4))*(G_km));
 
 
j7 =@(x) V3*(cos(x(2) - x(4))*(G_km) + sin(x(2) - x(4))*(B_km)) + 2*G_kk*V3 - V1*cos(theta1 - x(2))*(G_km) - V2*cos(x(1) - x(2))*(G_km) + V3*cos(x(2) - x(4))*(G_km) + V1*sin(theta1 - x(2))*(B_km) + ...
V2*sin(x(1) - x(2))*(B_km) + V3*sin(x(2) - x(4))*(B_km);
 
 
j8 =@(x) -V3*(V3*cos(x(2) - x(4))*(B_km) - V3*sin(x(2) - x(4))*(G_km));
 
 
j9 =@(x) V3*(V2*cos(x(1) - x(2))*(G_km) + V2*sin(x(1) - x(2))*(B_km));
 
 
j10 =@(x) -V3*(V1*cos(theta1 - x(2))*(G_km) + V2*cos(x(1) - x(2))*(G_km) - V3*cos(x(2) - x(3))*(G_km) + V1*sin(theta1 - x(3))*(B_km) + V2*sin(x(1) - x(2))*(B_km) + V3*sin(x(2) - x(4))*(B_km));
 
 
j11 =@(x) V3*(cos(x(2) - x(4))*(B_km) + sin(x(2) - x(4))*(G_km)) + 2*B_kk*V3 - V1*cos(theta1 - x(2))*(B_km) - V2*cos(x(1) - x(2))*(B_km) + V3*cos(x(2) - x(4))*(B_km) + V1*sin(theta1 - x(2))*(G_km) + V2*sin(x(1) - x(2))*(G_km) + V3*sin(x(2) - x(4))*(G_km);
 
 
j12 =@(x) -V3*(V3*cos(x(2) - x(4))*(G_km) - V3*sin(x(2) - x(4))*(B_km));
 
 
j13 =@(x) V4*(V2*cos(x(1) - x(4))*(B_km) + V2*sin(x(1) - x(4))*(G_km));
 
 
j14 =@(x) V4*(V3*cos(x(2) - x(4))*(B_km) + V3*sin(x(2) - x(4))*(G_km));
 
 
j15 =@(x) -V4*(cos(x(2) - x(4))*(G_km) - sin(x(2) - x(4))*(B_km));
 
 
j16 =@(x)-V4*(V1*cos(theta1 - x(4))*(B_km) + V2*cos(x(4) - x(1))*(B_km) + V3*cos(x(4) - x(2))*(B_km) + V1*sin(theta1 - x(4))*(G_km) + V2*sin(x(4) - x(1))*(G_km) + V3*sin(x(4) - x(2))*(G_km));

P2 = @(x) ((G(2,2))*V2^2) + V2*(V1*(G(2,1))*cos(x(1) - theta1)+ V1*(B(2,1))*sin(x(1) - theta1) + V2*(G(2,3))*cos(x(2) - x(1)) + V2*(B(2,3))*sin(x(2) - x(1)) + V3*(G(2,4))*cos(x(4) - x(1)) + ...
V3*(B(4,1))* sin(x(4) - x(1)));

P3 = @(x) ((G(3,3))*V3^2) + V3*(V1*(G(3,1))*cos(x(2) - theta1)+ V1*(B(3,1))*sin(x(2) - theta1) + V2*(G(3,2))*cos(x(2) - x(1)) + V2*(B(3,2))*sin(x(2) - x(1)) + V3*(G(3,4))*cos(x(4) - x(2)) + ...
V3*(B(3,4))* sin(x(4) - x(2)));

Q3 = @(x) ((B(3,3))*V3^2) + V3*(V1*(G(3,1))*sin(x(2) - theta1)- V1*(B(3,1))*cos(x(2) - theta1) + V2*(G(3,2))*sin(x(2) - x(1)) - V2*(B(3,2))*cos(x(2) - x(1)) + V3*(G(3,4))*sin(x(4) - x(2)) - ...
V3*(B(3,4))* cos(x(4) - x(2))) ;
 
P4 = @(x) (G(4,4) * V4^2) + V4*(V1*(G(4,1))*cos(x(4) - theta1) +V1*(B(4,1))*sin(x(4) - theta1)+ V2*(G(4,2))* cos(x(4)-x(1)) + V2*(B(4,2))*sin(x(4)-x(1)) + V3*(G(4,3))*cos(x(4)-x(2)) + V3*(B(4,3))*sin(x(4)-x(2))) ;
 
funcMatrix2 = @(x) [P2(x); P3(x); Q3(x); P4(x)];



while norm(c-funcMatrix2(x)) > esp
JacobMatrix = [j1(x) j2(x) j3(x) j4(x);
	           j5(x) j6(x) j7(x) j8(x);
			   j9(x) j10(x) j11(x) j12(x);
			   j13(x) j14(x) j15(x) j16(x)];

delta_x = inv(JacobMatrix)*(c-funcMatrix2(x)); % delta to approach tolerance
delta_x(1) = wrapTo2Pi(delta_x(1));
delta_x(2) = wrapTo2Pi(delta_x(2));
delta_x(4) = wrapTo2Pi(delta_x(4));
x = x + delta_x;

ne = norm(esp);
n = n+1;
end
error = norm(c - funcMatrix2(x))
n
x
c
ne
