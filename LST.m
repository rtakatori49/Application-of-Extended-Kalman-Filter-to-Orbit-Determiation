%% Local Sidereal Time
% Time it takes for a star to return to its same position overhead
% Based on class notes
function [theta] = LST(date, long)
y = date(1);
m = date(2);
d = date(3);
ut = date(4) + (date(5)/60) + (date(6)/(60*60)); % UT [hr]
% Date calculation
j_0 = 367*y-floor((7*(y + floor((m + 9)/12)))/4) + floor((275*m)/9) + d + 1721013.5;

% Time calculation
j2k = 2451545.0;
jc = 36525;
t_0 = (j_0-j2k)/jc;

% Greenwich sidereal time at 0hr
theta_g_0 = 100.4606184 + 36000.77004*t_0 + 0.000387933*((t_0)^2) -(2.58*10^-8)*((t_0)^3);

% Make LST between 0 and 360
while theta_g_0 > 360
    theta_g_0 = theta_g_0 - 360;
end
while theta_g_0 < -360
    theta_g_0 = theta_g_0 + 360;
end
% Greenwich sidereal time at time
theta_g = theta_g_0 + 360.98564724*(ut/24);

% Local sidereal time
theta = theta_g + long;

% Make LST between 0 and 360
while theta > 360
    theta = theta - 360;
end
while theta < -360
    theta = theta + 360;
end
end

