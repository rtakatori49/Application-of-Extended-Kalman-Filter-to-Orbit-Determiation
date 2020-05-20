function R = Cy(theta)
% R2  Find the rotation matrix about the number one axis for a rotation of
% theta radians
%

R = [cos(theta) 0 -sin(theta);
      0          1 0;
     sin(theta) 0 cos(theta)];