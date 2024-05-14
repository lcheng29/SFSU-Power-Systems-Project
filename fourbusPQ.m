clear;clc;

%Algorihtim for the Zbase
V1 = 1;
V_base = (345 * 10^3);
VA_base = (100*10^6);
Z_base = (V_base)^2/VA_base;
theta1 = 0;
V2 = 1.05;

%defining givens we have
eps = 1e-5;

P2sp = 2;
P3sp = 5;
P4sp = 4.0;
Q3sp = 1;
Q4sp = 1.5;
%vector to store power values
c = [P2sp ; P3sp ; Q3sp ; P4sp; Q4sp];

%Vector to store angles and initial guesses
theta2 = 0;
theta3 = 0;
theta4 = 0;
V3 = 1;
V4 = 1;
x = [theta2 ; theta3 ; V3 ; theta4 ; V4];

% define distance
D = 150;
% define the amount of busses

Z_12prebase = D*(0.037 + (j*0.376));

% mutual impedances
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

%% Admittance matrix formation

% mutual admittances
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

% self-admittances
Y_11 = (1/Z_12) + (1/Z_13) + (1/Z_14);
Y_22 = (1/Z_12) + (1/Z_23) + (1/Z_24); 
Y_33 = (1/Z_13) + (1/Z_23) + (1/Z_34);
Y_44 = (1/Z_14) + (1/Z_24) + (1/Z_34);


Y_admit = [Y_11 Y_12 Y_13 Y_14; Y_21 Y_22 Y_23 Y_24; Y_31 Y_32 Y_33 Y_34; Y_41 Y_42 Y_43 Y_44];

G = real(Y_admit);
B = imag(Y_admit);


%% Power matrix

P2 = @(x) (G(1,1)*V2^2) + V2*(V1*(G(2,1))*cos(x(1) - theta1)+ V1*(B(2,1))*sin(x(1) - theta1) + ...
x(3)*(G(2,3))*cos(x(1) - x(2)) + x(3)*(B(2,3))*sin(x(1) - x(2)) + x(5)*(G(2,4))*cos(x(1) - x(4)) + ...
x(5)*(B(2,4))* sin(x(1) - x(4))); % real power of bus 2

P3 = @(x) ((G(3,3))*x(3)^2) + x(3)*(V1*G(3,1)*cos(x(2) - theta1)+ V1*(B(3,1))*sin(x(2) - theta1) + ...
V2*(G(3,2))*cos(x(2) - x(1)) + V2*(B(3,2))*sin(x(2) - x(1)) + x(5)*(G(3,4))*cos(x(2) - x(4)) + ...
x(5)*(B(3,4)* sin(x(2) - x(4)))); % real power of bus 3

Q3 = @(x) (-B(3,3)*x(3)^2) + x(3)*(V1*G(3,1)*sin(x(2) - theta1) - V1*B(3,1)*cos(x(2) - theta1) + ...
V2*G(3,2)*sin(x(2) - x(1)) - V2*B(3,2)*cos(x(2) - x(1)) + x(3)*G(3,4)*sin(x(2) - x(4)) - ...
x(3)*B(3,4)*cos(x(2) - x(4))); % reactive power of bus 3

P4 = @(x) (G(4,4) * x(5)^2) + x(5)*((V1*(G(4,1))*cos(x(4) - theta1)) +V1*(B(4,1))*sin(x(4) - theta1)+ V2*(G(4,2))*...
cos(x(4)-x(1)) + V2*(B(4,2))*sin(x(4)-x(1)) + x(3)*(G(4,3))*cos(x(4)-x(2)) + x(3) * B(4,3) * sin(x(4) - x(2)));
% real power of bus 4

Q4 =  @(x) (-B(4,4) * x(5)^2) + x(5)*((V1*G(4,1)*sin(x(4) - theta1))- (V1*B(4,1)*cos(x(4) - theta1)) + ...
(V2*G(4,2) *sin(x(4) - x(1))) - V2*B(4,2)*cos(x(4) - x(1)) + x(3)*G(4,3)*sin(x(4) - x(2)) - ...
x(3) * B(4,3) * cos(x(4) - x(2))); % reactive power of bus 4

funcMatrix3 = @(x) [P2(x), P3(x), Q3(x), P4(x), Q4(x)].';

% please make Jacobian matrix 
% partial deriv. of P2 based on theta2, theta3, V3, theta4, V4
% partial deriv. of P3 based on theta2, theta3, V3, theta4, V4
% partial deriv. of Q3 based on theta2, theta3, V3, theta4, V4
% partial deriv. of P4 based on theta2, theta3, V3, theta4, V4
% partial deriv. of Q4 based on theta2, theta3, V3, theta4, V4
%% Jacobian matrix

% partial deriv. of P2 based on theta2, theta3, V3, theta4, V4
j1 =@(x) -V2*(V1*cos(theta1 - x(1))*(B(1,2)) - V2*cos(x(1) - x(2))*(B(2,3)) - x(3)*cos(x(1) - x(4))*(B(2,4)) + V1*sin(theta1 - x(1))*(G(1,2)) + V2*sin(x(1) - x(2))*(G(2,3)) + x(3)*sin(x(1) - x(4))*(G(2,4)));
j2 =@(x) -V2*(V2*cos(x(1) - x(3))*(B(2,3)) - V2*sin(x(1) - x(2))*(G(2,3)));
j3 =@(x) V2*(cos(x(1) - x(4))*(G(2,4) + sin(x(1) - x(4))*(B(2,4))));
j4 =@(x) -V2*(x(3)*cos(x(1) - x(4))*(B(2,4)) - x(3)*sin(x(1) - x(4))*(G(2,4)));
j5 =@(x) 0;

% partial deriv. of P3 based on theta2, theta3, V3, theta4, V4
j6 =@(x) x(3)*(V2*cos(x(1) - x(2))*(B(2,3)) + V2*sin(x(1) - x(2))*(G(2,3)));
j7 =@(x) -x(3)*(V1*cos(theta1 - x(2))*(B(1,3)) + V2*cos(x(1) - x(2))*(B(2,3)) - x(3)*cos(x(2) - x(4))*(B(3,4)) + V1*sin(theta1 - x(2))*(G(1,3)) + V2*sin(x(1) - x(2))*(G(2,3)) + x(3)*sin(x(2) - x(4))*(G(3,4)));
j8 =@(x) x(3)*(cos(x(2) - x(4))*(G(3,4)) + sin(x(2) - x(4))*(B(3,4))) + 2*G(3,3)*x(3) - V1*cos(theta1 - x(2))*(G(1,3)) - V2*cos(x(1) - x(2))*(G(2,3)) + x(3)*cos(x(2) - x(4))*(G(3,4)) + V1*sin(theta1 - x(2))*(B(1,3)) + V2*sin(x(1) - x(2))*(B(2,3)) + x(3)*sin(x(2) - x(4))*(B(3,4));
j9 =@(x) -x(3)*(x(3)*cos(x(2) - x(4))*(B(3,4)) - x(3)*sin(x(2) - x(4))*(G(3,4)));
j10 =@(x) 0;

% partial deriv. of Q3 based on theta2, theta3, V3, theta4, V4
j11 =@(x) x(3)*(V2*cos(x(1) - x(2))*(G(2,3)) - V2*sin(x(1) - x(2))*(B(2,3)));
j12 =@(x) x(3)*(x(3)*cos(x(2) - x(4))*(G(3,4)) - V2*cos(x(1) - x(2))*(G(2,3)) - V1*cos(theta1 - x(2))*(G(1,3)) + V1*sin(theta1 - x(2))*(B(1,3)) + V2*sin(x(1) - x(2))*(B(2,3)) + x(3)*sin(x(2) - x(4))*(B(3,4)));
j13 =@(x) 2*B(3,3)*x(3) - x(3)*(cos(x(2) - x(4))*(B(3,4)) - sin(x(2) - x(4))*(G(3,4))) + V1*cos(theta1 - x(2))*(B(1,3)) + V2*cos(x(1) - x(2))*(B(2,3)) - x(3)*cos(x(2) - x(4))*(B(3,4)) + V1*sin(theta1 - x(2))*(G(1,3)) + V2*sin(x(1) - x(2))*(G(2,3)) + x(3)*sin(x(2) - x(4))*(G(3,4));
j14 =@(x) -x(3)*(x(3)*cos(x(2) - x(4))*(G(3,4)) + x(3)*sin(x(2) - x(4))*(B(3,4)));
j15 =@(x) 0;

% partial deriv. of P4 based on theta2, theta3, V3, theta4, V4
j16 =@(x) x(5)*(V2*cos(x(1) - x(4))*(B(2,4)) + V2*sin(x(1) - x(4))*(G(2,4)));
j17 =@(x) x(5)*(x(3)*cos(x(2) - x(4))*(B(3,4)) + x(3)*sin(x(2) - x(4))*(G(3,4)));
j18 =@(x) -x(5)*(cos(x(2) - x(4))*(G(3,4)) - sin(x(2) - x(4))*(B(3,4)));
j19 =@(x) -x(5)*(V1*cos(theta1 - x(4))*(B(3,4)) + V2*cos(x(1) - x(4))*(B(2,4)) + x(3)*cos(x(2) - x(4))*(B(3,4)) + V1*sin(theta1 - x(4))*(G(1,4)) + V2*sin(x(1) - x(4))*(G(2,4)) + x(3)*sin(x(2) - x(4))*(G(3,4)));
j20 =@(x) 2*G(4,4)*x(5) - V1*cos(theta1 - x(4))*(G(1,4)) - V2*cos(x(1) - x(4))*(G(2,4)) - x(3)*cos(x(2) - x(4))*(G(3,4)) + V1*sin(theta1 - x(4))*(B(1,4)) + V2*sin(x(1) - x(4))*(B(2,4)) + x(3)*sin(x(2) - x(4))*(B(3,4));

% partial deriv. of Q4 based on theta2, theta3, V3, theta4, V4
j21 =@(x) x(5)*(V2*cos(x(1) - x(4))*(G(2,4)) - V2*sin(x(1) - x(4))*(B(2,4)));
j22 =@(x) x(5)*(x(3)*cos(x(2) - x(4))*(G(3,4)) - x(3)*sin(x(2) - x(4))*(B(3,4)));
j23 =@(x) x(5)*(cos(x(2) - x(4))*(B(3,4)) + sin(x(2) - x(4))*(G(3,4)));
j24 =@(x) -x(5)*(V1*cos(theta1 - x(4))*(G(1,4)) + V2*cos(x(1) - x(4))*(G(2,4)) + x(3)*cos(x(2) - x(4))*(G(3,4)) - V1*sin(theta1 - x(4))*(B(1,4)) - V2*sin(x(1) - x(4))*(B(2,4)) - x(3)*sin(x(2) - x(4))*(B(3,4)));
j25 =@(x) 2*B(4,4)*x(5) + V1*cos(theta1 - x(4))*(B(1,4)) + V2*cos(x(1) - x(4))*(B(2,4)) + x(3)*cos(x(2) - x(4))*(B(3,4)) + V1*sin(theta1 - x(4))*(G(1,4)) + V2*sin(x(1) - x(4))*(G(2,4)) + x(3)*sin(x(2) - x(4))*(G(3,4));

n = 0;

%% while loop iterations

while norm(c-funcMatrix3(x)) > eps
    
    JacobMatrix = [j1(x) j2(x) j3(x) j4(x) j5(x);
    j6(x) j7(x) j8(x) j9(x) j10(x);
    j11(x) j12(x) j13(x) j14(x) j15(x);
    j16(x) j17(x) j18(x) j19(x) j20(x);
    j21(x) j22(x) j23(x) j24(x) j25(x)];
    
    
    
    delta_x = inv(JacobMatrix) * (c-funcMatrix3(x)); % delta to approach tolerance
    % delta_x(1) = wrapTo2Pi(delta_x(1));
    % delta_x(2) = wrapTo2Pi(delta_x(2));
    % delta_x(4) = wrapTo2Pi(delta_x(4));
    x = x + delta_x;
    
    n = n+1;

end

x(1) = wrapTo2Pi(x(1)); 
x(2) = wrapTo2Pi(x(2));
x(4) = wrapTo2Pi(x(4));
error = norm(c - funcMatrix3(x))
n
x

 
