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

