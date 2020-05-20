%% Two Body Problem
function [dstatedt] = two_body(t,state)
global mue
%% Position [km]
% Component
x = state(1);
y = state(2);
z = state(3);
% Vector
r = [x y z]';
% Magnitude
R = norm(r);

%% Velocity [km/s]
% Component
dx = state(4);
dy = state(5);
dz = state(6);
% Vector
v = [dx dy dz]';
% Magnitude
V = norm(v);

%% Acceleration [km/s^2]
% Component
ddx = -mue*x/R^3;
ddy = -mue*y/R^3;
ddz = -mue*z/R^3;
% Vector
a = [ddx ddy ddz]';

dstatedt = [v;a]; % Return state vector for next step
end

