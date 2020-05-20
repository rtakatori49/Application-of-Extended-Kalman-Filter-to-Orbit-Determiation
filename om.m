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

