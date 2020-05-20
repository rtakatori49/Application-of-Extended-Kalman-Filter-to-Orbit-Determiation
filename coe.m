%% Aero 215 HW3
% Ryo Takatori
% 10/27/17
% Introduction to Aerospace Design

function [a, E, H, inc, RAAN, omega, theta, coe_list] = coe(mu, r, v)
%% Given inertial and velocity vectors
R = norm(r); % Inertial postion vector [km]
V = norm(v); % Velocity vector [km/s]

%% Gravitational Constant and Semi-major axis (a)
epsilon = ((V.^2)/2)-(mu/R); % Equation for specific mechanical energy [MJ/kg]
a = -mu/(2*epsilon); % Equation for semi-major axis

%% Eccentricity (e)
ecc = (1/mu)*((((V^2)-(mu/R))*r)-(dot(r,v)*v)); % Equation for eccentricity
E = norm(ecc); % Magnitude of the eccentricity

%% Inclination (inc)
h = cross(r,v); % Cross product of r and v
H = norm(h); % Magnitude of h
inc = acosd(h(3)/H); % Equation for inclination angle [degree]

%% RAAN
k = [0,0,1]; % K hat vector
I = [1,0,0]; % I hat vector
n = cross(k,h); % Cross product of k hat and h
N = norm(n); % Magnitude of n
RAAN = acosd(dot(I,n)/N); % Equation for RAAN [degree]
if n(2) < 0 % if statement for checking quadrant ambiguity
    RAAN = 360-RAAN;
end

%% Argument of perigee (omega)
omega = acosd(dot(n,ecc)/(N*E)); % Equation for the argument of perigee [degree]
if ecc(3) < 0 % if statement for checking quadrant ambiguity
    omega = 360-omega;
end
    
%% True anomally (theta)
theta = acosd(dot(ecc,r)/(E*R)); % Equation for true anomally [degree]
if dot(r,v) < 0 % if statement for checking quadrant ambiguity
    theta = 360-theta;
end
coe_list = [a, E, H, inc, RAAN, omega, theta]';
end