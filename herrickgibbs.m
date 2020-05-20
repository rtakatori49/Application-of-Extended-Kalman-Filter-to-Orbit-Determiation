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

