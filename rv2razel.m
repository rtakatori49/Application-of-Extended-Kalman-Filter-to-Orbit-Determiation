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

