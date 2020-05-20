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


