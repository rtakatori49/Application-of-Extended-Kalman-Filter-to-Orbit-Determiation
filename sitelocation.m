%% Site Location
% Finds observation site using LST in ECI [km]
function [rsite, theta, phi] = sitelocation(date, lla)
global re
rp = 6357; % Radius at pole [km]
f = (re-rp)/re; % Oblateness factor
phi = deg2rad(lla(1)); % Geodetic latitude [rad]
long = lla(2); % Longitude [deg]
H = lla(3)/1000; % Elevation of site [km]
[theta] = LST(date, long); % LST [deg]
theta = deg2rad(theta); % LST [rad]
rsite = horzcat(((re/sqrt(1-(2*f-f^2)*sin(phi)^2))+H)*cos(phi)*[cos(theta) ...
    sin(theta)], (((re*(1-f)^2)/sqrt(1-(2*f-f^2)*sin(phi)^2))+H)*sin(phi))';
end