function R = Cz(theta)
% R3  Find the rotation matrix about the number one axis for a rotation of
% theta radians
%

R = [cos(theta) sin(theta) 0;
     -sin(theta) cos(theta)  0;
     0          0           1];