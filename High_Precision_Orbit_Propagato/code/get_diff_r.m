diff_r = zeros(size(Eph_eci,1),1);
% Define the start and stop times for the simulation
et_start = cspice_str2et('January 1, 2000 00:00:00');
% et_stop = cspice_str2et('December  26, 2022 00:00:00');
et_stop = cspice_str2et('April  13, 2029:23:0:0');

% Define the step size for the simulation
et_step = 3600; %s

% Generate the time vector for the simulation
et = et_start:et_step:et_stop;
[Y_EARTH, ~] = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');

for i = 1 : size(Eph_eci,1)
    r = Eph_eci(i,2:4)' - Y_EARTH(1:3,i).*1000;
    diff_r(i) = norm( r );
end
x = 1 : size(Eph_eci,1);
x = x / 365.25 / 24  ;
plot(x,diff_r)
ylabel('高度')
xlabel('年')
set(gca, 'YScale', 'log')