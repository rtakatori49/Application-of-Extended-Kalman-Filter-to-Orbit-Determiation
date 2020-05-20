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