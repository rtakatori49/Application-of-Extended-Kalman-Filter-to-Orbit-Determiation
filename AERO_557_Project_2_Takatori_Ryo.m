%% Final Project
% AERO 557
% Ryo Takatori
% 03/14/2020

% Housekeeping
clc, clear all, close all

% Common variables
global mue re options
mue = 398600; % Earth gravitational constant [km^3/s^2]
re = 6378; % Radius of Earth [km]
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Ode settings

%% Extended Kalman Filter Applied to Homework 2 Problem 2
disp('Extended Kalman Filter Applied to Homework 2 Problem 2:')
% Data
data = load('Data4Problem2Hw2_2020test.txt');
date_ut = data(:,1:3); % UT date
time_ut = data(:,4:6); % UT time
ut = data(:,1:6); % UT date + time
range = data(:,9); % Range [km]
az = data(:,7); % Azimuth [deg]
el = data(:,8); % Elevation [deg]
% Site
lat = 21.5748; % Latitude [deg]
long = -158.2706; % Longitude [deg]
alt = 300.2; % Altitude [m]
lla = [lat long alt];
% Error
error_range = 92.5/1000; % Range error [m]
error_az = deg2rad(0.0224); % Azimuth error [rad]
error_el = deg2rad(0.0139); % Elevation error [rad]
% Some initial epoch chosen
ut0 = [1995 1 29 2 38 0];
% True orbit taken from Vallado Text book
trueorbit = [5753.173e3,2673.361e3,3440.304e3,4.324207e3,-1.924299e3,-5.728216e3]/1000';
% Extended Kalman Filter
[x0, x, obs] = ekf(lla, ut, range, az, el, error_range, error_az, error_el, ut0, trueorbit);

%% Functions

%% Orbit Determination using Extended Kalman Filter
function [xekf0, x_store, rho] = ekf(lla, ut, range, az, el, erho, eaz, eel, ut0, trueorbit)
global options mue
N = length(range); % Sample size
N3 = floor(N/3); % Amount of observation sets available for gibbs
% Site, rho vector, juliandate, and time difference
for i = 1:N
    [obs(:,i), theta, phi] = sitelocation(ut(i,:), lla);
    [ra(i), dec(i)] = AzEl2RaDec(az(i),el(i),lla(1),lla(2),datestr(ut(i,:),...
        'yyyy/mm/dd HH:MM:SS'));
    rho(:,i) = range(i)*[cosd(ra(i))*cosd(dec(i)) sind(ra(i))*cosd(dec(i)) sind(dec(i))]';
    jd(i) = juliandate(ut(i,:));
    dt(i) = (jd(1)-jd(i))*24*60*60;
end
% R vector
r = rho + obs;
% Gibbs to get v2 vector
for i = 1:N3
    r2(:,i) = r(:,3*i-1);
    [v2(:,i)] = gibbs(r(:,3*i), r2(:,i), r(:,3*i-2),...
        jd(:,3*i), jd(:,3*i-1), jd(:,3*i-2));
    state(:,i) = [r2(:,i); v2(:,i)];
    tspan(i,:) = [0 dt(:,3*i-1)];
    [tnew{i}, statenew{i}] = ode45(@two_body, tspan(i,:), state(:,i), options); % ode45;
end
% Nominal RV vector from average of three vectors
for i = 1:N3
    for j = 1:3
        rnomvec(j,i) = statenew{i}(end,j);
        vnomvec(j,i) = statenew{i}(end,j+3);
    end
end
for i = 1:3
    rnom(i,1) = mean(rnomvec(i,:));
    vnom(i,1) = mean(vnomvec(i,:));
end
Xnom = [rnom;vnom]; %Xnom vector
% Some initial epoch to start
t0 = (juliandate(ut0)-jd(1))*24*60*60;
tspan0 = [0 t0];
[t0vec, x0vec] = ode45(@two_body, tspan0, Xnom, options); % ode45
x0 = x0vec(end,:);
trueobs = horzcat(horzcat(range,az),el); % True observation
% Weight matrix
W = ([(1/erho)^2 0 0;
    0 (1/eaz)^2 0;
    0 0 (1/eel)^2]);
n = length(trueobs(:,1)); % Number of observation

% Extended Kalman Filter
% Initialize
% Initial Covariance matrix
P = zeros(6);
for i=1:3
    P(i,i)=1e8;
end
for i=4:6
    P(i,i)=1e3;
end
x = x0';
Q = 0; % We set Q to zero for this problem, could be changed
R = inv(W); % R noise matrix
radobs = deg2rad(trueobs(:,2:3)); % Observation angles in radians
for i = 1:n
    x_old = x; % Store x as previous x
    phi = eye(6,6); % Initial STM is an identity matrix
    % Calculate time difference between each observation
    if i == 1
        t = -t0;
    else
        t = (jd(i)-jd(i-1))*24*60*60;
    end
    % Propgate prediction x
    tspankf = [0 t];
    [tvec, xvec] = ode45(@two_body, tspankf, x_old, options);
    x = xvec(end,:)';
    % Construct state matrix F
    F = sm(x(1:3));
    % Propagate STM phi
    [tvec,phivec] = ode45(@stm, tspankf, phi, options, F);
    philong = phivec(end,:);
    for j = 1:6
        phi(:,j) = philong(1+6*(j-1):6+6*(j-1));
    end
    % Predict covariance matrix
    P = phi*P*phi'+Q;
    % Calculate predicted observation angles
    [rhopre, rangepre, azpre, elpre, theta] = rv2razel(x(1:3), ut(i,:), lla);
    Csez2eci = Cz(-deg2rad(theta))*Cy(-((pi/2)-deg2rad(lla(2)))); % SEZ to ECI
    rhopre = Csez2eci'*rhopre; % Set from ECI to SEZ
    % Construct observation matrix H
    [H] = om(rhopre, lla, theta);
    % Calculate residual
    b = [trueobs(i,1) radobs(i,:)]' - [rangepre deg2rad(azpre) deg2rad(elpre)]';
    % Update Kalman gain
    K = P*H'*pinv(H*P*H'+R);
    % Calculate difference
    delx = K*(b);
    % Update state
    x = x + delx;
    % Update covariance
    P = (eye(6)-K*H)*P;
    % Store state for graph
    x_store(:,i) = x;
end

% Propagate back to initial epoch to get estimate epoch
tspanzero = [0 dt(end)];
[tzero, xzero] = ode45(@two_body, tspanzero, x, options);
[ttrue, xtrue] = ode45(@two_body, -tspanzero, trueorbit, options);
xekf0 = xzero(end,:)';

% Plots
figure
plot3(xzero(:,1),xzero(:,2),xzero(:,3))
hold on
plot3(xtrue(:,1),xtrue(:,2),xtrue(:,3))
plot3(x_store(1,:),x_store(2,:),x_store(3,:),'o');
title('Orbit for Observation Duration')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('EKF Orbit','True Orbit','EKF Points')
grid on

tspanfull = [0 24*60*60];
[tzero2, xzero2] = ode45(@two_body, tspanfull, x, options);
[ttrue2, xtrue2] = ode45(@two_body, tspanfull, trueorbit, options);
figure
plot3(xzero2(:,1),xzero2(:,2),xzero2(:,3))
hold on
plot3(xtrue2(:,1),xtrue2(:,2),xtrue2(:,3))
plot3(x_store(1,:),x_store(2,:),x_store(3,:),'o');
title('Entire Orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('EKF Orbit','True Orbit','EKF Points')
grid on

e = eig(P);
for i = 1:3
confidence(i) = sqrt(P(i,i))*1000;
end
fprintf('R vector: [%f %f %f] [km]\n',xekf0(1:3));
fprintf('R Magnitude: %f [km]\n',norm(xekf0(1:3)));
fprintf('V vector: [%f %f %f] [km/s]\n',xekf0(4:6));
fprintf('V Magnitude: %f [km/s]\n',norm(xekf0(4:6)));
fprintf('Confidence in R vector: [%f %f %f] [m]\n',confidence)
fprintf('Error ellipsoid in R vector: [%f %f %f] [m]\n',e(1:3)*1000)

% Comparison of states and orbital elements
[atrue, Etrue, Htrue, inctrue, RAANtrue, omegatrue, thetatrue, coe_list_true] = coe(mue, trueorbit(1:3), trueorbit(4:6));
[aEKF, EEKF, HEKF, incEKF, RAANEKF, omegaEKF, thetaEKF, coe_list_EKF] = coe(mue, xekf0(1:3), xekf0(4:6));
variablenames = {'EKF','True','Difference'};
rownames = {'Rx';'Ry';'Rz';'Vx';'Vy';'Vz';'R';'V';'a';'e';'h';'i';'RAAN';'w';'TA'};
list1 = [xekf0;norm(xekf0(1:3));norm(xekf0(4:6));coe_list_EKF];
list2 = [trueorbit';norm(trueorbit(1:3));norm(trueorbit(4:6));coe_list_true];
difference = abs(list1-list2);
table_1 = table(list1,list2,difference,'VariableNames',variablenames,'RowNames',rownames);
disp('Comparison of States and Orbital Elements:')
disp(table_1)
end

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

function [RA DEC] = AzEl2RaDec(Az,El,lat,lon,time)
% Programed by Darin C. Koblick 6/16/2009
%--------------------------------------------------------------------------
% Updated:                                                      Date:
% - quadrant check bug fix                                  1/22/2010
% - vectorized for speed                                    1/22/2010
%--------------------------------------------------------------------------
% External Function Call Sequence:
% [RA DEC] = AzEl2RaDec(0,0,0,-104,'1992/08/20 12:14:00')
%
% Worked Example pg. 262 Vallado
% [RA DEC] = AzEl2RaDec(210.8250667,23.8595052,39.007,-104.883,'1994/05/14 13:11:20.59856')
% [alpha_t,delta_t] = AzEl2RaDec(Beta,el,phi,lamda,'yyyy/mm/dd hh:mm:ss')
%
% Function Description:
%--------------------------------------------------------------------------
% AzEl2RaDec will take the Azimuth and Elevation in the local horizon
% reference frame, site latitude and longitude as well as a time in GMT
% and output the Right Ascension and Declination in the topocentric coordinate frame.
%
% Inputs:                                                       Format:
%--------------------------------------------------------------------------
% Local Azimuth Angle   (degrees)                               [N x 1]
% Local Elevation Angle (degrees)                               [N x 1]
% Lat (Site Latitude in degrees -90:90 -> S(-) N(+))            [N x 1]
% Lon (Site Longitude in degrees -180:180 W(-) E(+))            [N x 1]
% UTC (Coordinated Universal Time YYYY/MM/DD hh:mm:ss)          [N x 1]
%
% Outputs:                                                      Format:
%--------------------------------------------------------------------------
% Topocentric Right Ascension (Degrees)   [N x 1]
% Topocentric Declination Angle (Degrees)                       [N x 1]
%
%
% External Source References:
% Fundamentals of Astrodynamics and Applications 
% D. Vallado, Second Edition
% Example 3-5. Finding Local Siderial Time (pg. 192) 
% Algorithm 28: AzElToRaDec (pg. 259)
% -------------------------------------------------------------------------

%Example 3-5
[yyyy mm dd HH MM SS] = datevec(datenum(time,'yyyy/mm/dd HH:MM:SS'));
JD = juliandate(yyyy,mm,dd,HH,MM,SS);
T_UT1 = (JD-2451545)./36525;
ThetaGMST = 67310.54841 + (876600*3600 + 8640184.812866).*T_UT1 ...
+ .093104.*(T_UT1.^2) - (6.2*10^-6).*(T_UT1.^3);
ThetaGMST = mod((mod(ThetaGMST,86400.*(ThetaGMST./abs(ThetaGMST)))./240),360);
ThetaLST = ThetaGMST + lon;

%Algorithm 28
DEC = asind(sind(El).*sind(lat)+cosd(El).*cosd(lat).*cosd(Az));
LHA = atan2(-sind(Az).*cosd(El)./cosd(DEC), ...
    (sind(El)-sind(DEC).*sind(lat))./(cosd(DEC).*cosd(lat))).*(180/pi);
RA = mod(ThetaLST-LHA,360);

function jd = juliandate(year, month, day, hour, min, sec) 
YearDur = 365.25;
for i = length(month):-1:1
    if (month(i)<=2)
        year(i)=year(i)-1;
        month(i)=month(i)+12;
    end
end
A = floor(YearDur*(year+4716));
B = floor(30.6001*(month+1));
C = 2;
D = floor(year/100);
E = floor(floor(year/100)*.25);
F = day-1524.5;
G = (hour+(min/60)+sec/3600)/24;
jd =A+B+C-D+E+F+G;
end
end

%% Gibbs Method
function [v2] = gibbs(r1, r2, r3, jd1, jd2, jd3)
global mue
Z12 = cross(r1, r2);
Z23 = cross(r2, r3);
Z31 = cross(r3, r1);
r1mag = norm(r1);
r2mag = norm(r2);
r3mag = norm(r3);
% Coplanar check
alpha_cop = 90 - acosd(dot(Z23,r1)/(norm(Z23)*r1mag));
if alpha_cop > 3
    error('Not coplanar.')
end
% Check spacing
alpha12 = acosd(dot(r1,r2)/(r1mag*r2mag));
alpha23 = acosd(dot(r2,r3)/(r2mag*r3mag));
if alpha12 < 1 && alpha23 < 1
    fprintf('Spacing too small. Moving into Herrick-Gibbs.\n\n')
%     prompt1 = 'Please enter the Julian date corresponding to r1:\n';
%     jd1 = input(prompt1);
%     prompt2 = 'Please enter the Julian date corresponding to r3:\n';
%     jd2 = input(prompt2);
%     prompt3 = 'Please enter the Julian date corresponding to r3:\n';
%     jd3 = input(prompt3);
    [v2] = herrickgibbs(r1, r2, r3, jd1, jd2, jd3);
else
    N = r1mag*Z23 + r2mag*Z31 + r3mag*Z12;
    D = Z12 + Z23 + Z31;
    S = (r2mag-r3mag)*r1 + (r3mag-r1mag)*r2 + (r1mag-r2mag)*r3;
    B = cross(D, r2);
    Lg = sqrt(mue/(norm(N)*norm(D)));
    v2 = (Lg/r2mag)*B + Lg*S;
end
end

%% Herrick Gibbs
function [v2] = herrickgibbs(r1, r2, r3, jd1, jd2, jd3)
global mue
dt31 = (jd3 - jd1)*(24*60*60);
dt32 = (jd3 - jd2)*(24*60*60);
dt21 = (jd2 - jd1)*(24*60*60);
r1mag = norm(r1);
r2mag = norm(r2);
r3mag = norm(r3);
v2 = -dt32*((1/(dt21*dt31))+(mue/(12*r1mag^3)))*r1 +...
    (dt32-dt21)*((1/(dt21*dt32))+(mue/(12*r2mag^3)))*r2 +...
    dt21*((1/(dt32*dt31))+(mue/(12*r3mag^3)))*r3;
end

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

%% State Matrix of Two-Body Motion
function [F] = sm(rvec)
global mue
r = norm(rvec);
F = zeros(6);
F(1:3,4:6) = eye(3);
for i = 1:3
    for j = 1:3
        if i == j
            F(i+3,j) = (-mue/r^3)+((3*mue*rvec(i)^2)/r^5);
        end
        F(i+3,j) = (3*mue*rvec(i)*rvec(j))/r^5;
    end
end
end

%% State Transition Matrix
function [dstatedt] = stm(t, state, F)
global mue
phi = zeros(6);
for i = 1:6
    phi(:,i) = state(1+6*(i-1):6+6*(i-1));
end
dphi = F*phi;

dstatedt = [dphi(:,1); dphi(:,2); dphi(:,3); dphi(:,4); dphi(:,5); dphi(:,6)];
end

%% RV Vector to Range Azimuth Elevation
function [rho, range, az, el, theta] = rv2razel(r, date, lla)
[rsite, theta, phi] = sitelocation(date, lla); % Site (ECI)
    rho = r - rsite; % Line of sight (ECI)
    range = norm(rho); % Range [km]
    unitrho = rho/range; % Unit LoS (ECI)
    dec = asind(unitrho(3)); % Declination [deg]
    ra = acosd(unitrho(1)/cosd(dec)); % Right ascention [deg]
    % Check RA angle
    if rho(2) <= 0
        ra = 360 - ra;
    end
    % RA DEC to Az El
    [az, el] = RaDec2AzEl(ra,dec,lla(1),lla(2),datestr(date,'yyyy/mm/dd HH:MM:SS'));
end

function [Az El] = RaDec2AzEl(Ra,Dec,lat,lon,time)
% Programed by Darin C. Koblick 01/23/2010
%--------------------------------------------------------------------------
% External Function Call Sequence:
% [Az El] = RaDec2AzEl(0,0,0,-104,'1992/08/20 12:14:00')
%
% Worked Example: pg. 262 Vallado
%[Az El] = RaDec2AzEl(294.9891115,-20.8235624,39.007,-104.883,'1994/05/14 13:11:20.59856')
%[210.7514  23.9036] = RaDec2AzEl(294.9891115,-20.8235624,39.007,-104.883,'1994/05/14 13:11:20.59856')
%
% Worked Example: http://www.stargazing.net/kepler/altaz.html
% [Az El] = RaDec2AzEl(344.95,42.71667,52.5,-1.91667,'1997/03/14 19:00:00')
% [311.92258 22.40100] = RaDec2AzEl(344.95,42.71667,52.5,-1.91667,'1997/03/14 19:00:00')
%
% [Beta,el] = RaDec2AzEl(alpha_t,delta_t,phi,lamda,'yyyy/mm/dd hh:mm:ss')
%
% Function Description:
%--------------------------------------------------------------------------
% RaDec2AzEl will take the Right Ascension and Declination in the topocentric 
% reference frame, site latitude and longitude as well as a time in GMT
% and output the Azimuth and Elevation in the local horizon
% reference frame.
%
% Inputs:                                                       Format:
%--------------------------------------------------------------------------
% Topocentric Right Ascension (Degrees)                         [N x 1]
% Topocentric Declination Angle (Degrees)                       [N x 1]
% Lat (Site Latitude in degrees -90:90 -> S(-) N(+))            [N x 1]
% Lon (Site Longitude in degrees -180:180 W(-) E(+))            [N x 1]
% UTC (Coordinated Universal Time YYYY/MM/DD hh:mm:ss)          [N x 1]
%
% Outputs:                                                      Format:
%--------------------------------------------------------------------------
% Local Azimuth Angle   (degrees)                               [N x 1]
% Local Elevation Angle (degrees)                               [N x 1]
%
%
% External Source References:
% Fundamentals of Astrodynamics and Applications 
% D. Vallado, Second Edition
% Example 3-5. Finding Local Siderial Time (pg. 192) 
% Algorithm 28: AzElToRaDec (pg. 259)
% -------------------------------------------------------------------------

%Example 3-5
[yyyy mm dd HH MM SS] = datevec(datenum(time,'yyyy/mm/dd HH:MM:SS'));
JD = juliandate(yyyy,mm,dd,HH,MM,SS);
T_UT1 = (JD-2451545)./36525;
ThetaGMST = 67310.54841 + (876600*3600 + 8640184.812866).*T_UT1 ...
+ .093104.*(T_UT1.^2) - (6.2*10^-6).*(T_UT1.^3);
ThetaGMST = mod((mod(ThetaGMST,86400*(ThetaGMST./abs(ThetaGMST)))/240),360);
ThetaLST = ThetaGMST + lon;

%Equation 4-11 (Define Siderial Time LHA)
LHA = mod(ThetaLST - Ra,360);

%Equation 4-12 (Elevation Deg)
El = asind(sind(lat).*sind(Dec)+cosd(lat).*cosd(Dec).*cosd(LHA));

%Equation 4-13 / 4-14 (Adaptation) (Azimuth Deg)
%Az = mod(atand(-(sind(LHA).*cosd(Dec)./(cosd(lat).*sind(Dec) -
%sind(lat).*cosd(Dec).*cosd(LHA)))),360);
Az = mod(atan2(-sind(LHA).*cosd(Dec)./cosd(El),...
    (sind(Dec)-sind(El).*sind(lat))./(cosd(El).*cosd(lat))).*(180/pi),360);


function jd = juliandate(year, month, day, hour, min, sec) 
YearDur = 365.25;
for i = length(month):-1:1
    if (month(i)<=2)
        year(i)=year(i)-1;
        month(i)=month(i)+12;
    end
end
A = floor(YearDur*(year+4716));
B = floor(30.6001*(month+1));
C = 2;
D = floor(year/100);
E = floor(floor(year/100)*.25);
F = day-1524.5;
G = (hour+(min/60)+sec/3600)/24;
jd =A+B+C-D+E+F+G;
end
end

function R = Cz(theta)
% R3  Find the rotation matrix about the number one axis for a rotation of
% theta radians
%

R = [cos(theta) sin(theta) 0;
     -sin(theta) cos(theta)  0;
     0          0           1];
end

function R = Cy(theta)
% R2  Find the rotation matrix about the number one axis for a rotation of
% theta radians
%

R = [cos(theta) 0 -sin(theta);
      0          1 0;
     sin(theta) 0 cos(theta)];
end

%% Observation Matrix
function [H] = om(rhosez, lla, theta)
Csez2eci = Cz(-deg2rad(theta))*Cy(-((pi/2)-deg2rad(lla(2)))); % SEZ to ECI
rho = norm(rhosez);
rho12 = norm(rhosez(1:2));
dobsdrhosez = [rhosez(1)/rho rhosez(2)/rho rhosez(3)/rho;
    -rhosez(2)/rho12^2 rhosez(1)/rho12^2 0;...
    rhosez(1)*rhosez(3)/((rho^2)*rho12) rhosez(2)*rhosez(3)/((rho^2)*rho12) -rho12/(rho^2)];
H = [dobsdrhosez*Csez2eci' zeros(3)];
end

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
