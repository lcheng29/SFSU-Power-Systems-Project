syms G1 G2 G3 G4 B1 B2 B3 B4 V1 V2 V3 V4 theta1 theta2 theta3 theta4 G22 G33 G44 B33 G55 B55 G5 B5
P2 = ((G22)*V2^2) + V2*(V1*(G2-G1)*cos(theta2 - theta1)+ V1*(B2-B1)*sin(theta2 - theta1) + ...
V2*(G2-G3)*cos(theta2 - theta3) + V2*(B2-B3)*sin(theta2 - theta3) + V3*(G2-G4)*cos(theta2 - theta4) + ...
V3*(B2-B4)* sin(theta2 - theta4));

syms G1 G2 G3 G4 B1 B2 B3 B4 V1 V2 V3 V4 theta1 theta2 theta3 theta4 G22 G33 G44 B33
P3 = ((G33)*V3^2) + V3*(V1*(G3-G1)*cos(theta3 - theta1)+ V1*(B3-B1)*sin(theta3 - theta1) + ...
V2*(G3-G2)*cos(theta3 - theta2) + V2*(B3-B2)*sin(theta3 - theta2) + V3*(G3-G4)*cos(theta3 - theta4) + ...
V3*(B3-B4)* sin(theta3 - theta4));

syms G1 G2 G3 G4 B1 B2 B3 B4 V1 V2 V3 V4 theta1 theta2 theta3 theta4 G22 G33 G44 B33
Q3 = ((B33)*V3^2) + V3*(V1*(G3-G1)*sin(theta3 - theta1)- V1*(B3-B1)*cos(theta3 - theta1) + ...
V2*(G3-G2)*sin(theta3 - theta2) - V2*(B3-B2)*cos(theta3 - theta2) + V3*(G3-G4)*sin(theta3 - theta4) - ...
V3*(B3-B4)* cos(theta3 - theta4));

syms G1 G2 G3 G4 B1 B2 B3 B4 V1 V2 V3 V4 theta1 theta2 theta3 theta4 G22 G33 G44 B33
P4 = (G44 * V4^2) + V4*(V1*(G4 - G1)*cos(theta4 - theta1) +V1*(B4 - B1)*sin(theta4 - theta1)+ V2*(G4 - G2)*...
cos(theta4-theta2) + V2*(B4-B2)*sin(theta4-theta2) + V3*(G4-G3)*cos(theta4-theta3) + V3*(B4-B3)*sin(theta4-theta3));

syms G1 G2 G3 G4 B1 B2 B3 B4 V1 V2 V3 V4 theta1 theta2 theta3 theta4 G22 G33 G44 B33 B44
Q4 = ((B44)*V4^2) + V4*(V1*(G4-G1)*sin(theta4 - theta1)- V1*(B4-B1)*cos(theta4 - theta1) + ...
V2*(G4-G2)*sin(theta4 - theta2) - V2*(B4-B2)*cos(theta4 - theta2) + V3*(G4-G3)*sin(theta4 - theta3) - ...
V3*(B4-B3)* cos(theta4 - theta3));

f = [P2; P3 ; Q3 ; P4 ; Q4]; %making a 1x5 of known outputs
unknown = [theta2 theta3 V3 theta4 V4] %making a 5x1 matrixs for unknown/dynamic inputs
jjwatt = jacobian(f,unknown)
j1 =jjwatt(1,1)
j2 =jjwatt(1,2)
j3 =jjwatt(1,3)
j4 =jjwatt(1,4)
j5 =jjwatt(1,5)
j6 =jjwatt(2,1)
j7 =jjwatt(2,2)
j8 =jjwatt(2,3)
j9 =jjwatt(2,4)
j10 =jjwatt(2,5)
j11 =jjwatt(3,1)
j12 =jjwatt(3,2)
j13 =jjwatt(3,3)
j14 =jjwatt(3,4)
j15 =jjwatt(3,5)
j16 =jjwatt(4,1)
j17 =jjwatt(4,2)
j18 =jjwatt(4,3)
j19 =jjwatt(4,4)
j20 =jjwatt(4,5)
j21 =jjwatt(5,1)
j22 =jjwatt(5,2)
j23 =jjwatt(5,3)
j24 =jjwatt(5,4)
j25 =jjwatt(5,5)
