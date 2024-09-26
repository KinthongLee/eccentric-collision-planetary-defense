function delta_r_min = Propagation(target,delta_v,min_distance,dt,Step,N_Step)


% Initial States
% create initial datetime object
years = year(dt);
months = month(dt);
days = day(dt);
hours = hour(dt);
minutes = minute(dt);
secs = second(dt);

Mjd0_UTC = Mjday(years, months, days, hours, minutes, secs);

% Get starting position & velocity from SPICE:
months = monthToString(months);
str = sprintf(' %s %g , %g %g:%g:%g', months, days, years, hours, minutes, secs);
et_start = cspice_str2et(str);
[Y, ~] = cspice_spkezr(target, et_start, 'J2000', 'NONE', 'SUN');
Y = Y.*1000; % m/s

% Generate delta_v to 径向
Y(4:6,1) = Y(4:6) + delta_v;
% ------------------ propagation-------------------------------------------
 Eph_eci = Ephemeris(Y, N_Step, Step, Mjd0_UTC);

et = et_start : Step : et_start + N_Step*Step;

%-------------------------------- Get Height value & time----------------------
[Y_EARTH, ~] = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');
Y_EARTH = Y_EARTH.*1000;
% [Y_APOPHIS, ~] = cspice_spkezr(target, et, 'J2000', 'NONE', 'SUN');
% Y_APOPHIS = Y_APOPHIS.*1000;

diff_r = zeros(length(Y_EARTH),1);
% diff_r2 = zeros(length(Y_EARTH),1);

for i = 1 : length(Y_EARTH)
    r = Eph_eci(i,2:4)' - Y_EARTH(1:3,i);
    % r2 = Y_APOPHIS(1:3,i) - Y_EARTH(1:3,i);
    diff_r(i) = norm( r );
    % diff_r2(i) = norm(r2);
end

% Get the minimum height

% delta_r_min = min(diff_r) - min(diff_r2);
delta_r_min = min(diff_r) - min_distance;

end