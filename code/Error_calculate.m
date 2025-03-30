addpath('C:\Users\windows\Desktop\High Precision Orbit Propagator\mice\src\mice')
addpath('C:\Users\windows\Desktop\High Precision Orbit Propagator\mice\lib')
addpath('C:\Users\windows\Desktop\High Precision Orbit Propagator\mice')
cspice_furnsh({...
    'C:\Users\windows\Desktop\High Precision Orbit Propagator\mice\de441_part-1.bsp',...
    'C:\Users\windows\Desktop\High Precision Orbit Propagator\mice\de441_part-2.bsp', ...
    'C:\Users\windows\Desktop\High Precision Orbit Propagator\mice\naif0007.txt',...
    'C:\Users\windows\Desktop\High Precision Orbit Propagator\mice\2099942.bsp'...
    });


% Define the start and stop times for the simulation
et_start = cspice_str2et('January 1, 2000 00:00:00');
% et_stop = cspice_str2et('December  26, 2022 00:00:00');
et_stop = cspice_str2et('April  13, 2029 23:00:00');

% Define the step size for the simulation
et_step = 3600; %s

% Generate the time vector for the simulation
et = et_start:et_step:et_stop;

% Preallocate arrays to store the results
positions = zeros(length(et),3);
velocities = zeros(length(et),3);

% Loop through the time vector and compute the positions and velocities
for i = 1:length(Eph_eci)
    [pos, ltime] = cspice_spkezr('2099942', et(i), 'J2000', 'NONE', 'SUN');
    positions(i, 1:3) = pos(1:3).*1000; % m/s
    velocities(i, 1:3) = pos(4:6).*1000; % km/s
end

diffx = zeros( length(Eph_eci) ,1);
diffy = zeros( length(Eph_eci) ,1);
diffz = zeros( length(Eph_eci) ,1);
r_real = zeros( length(Eph_eci) ,1);
r = zeros( length(Eph_eci) ,1);
diff_r = zeros( length(Eph_eci) ,1);
% 作误差曲线
for i = 1 : length(Eph_eci)
    diffx(i) = abs( positions(i,1) - Eph_eci(i,2) ) ./ positions(i,1) .* 100; % in %
    diffy(i) = abs( positions(i,2) - Eph_eci(i,3) ) ./ positions(i,2) .* 100; % in %
    diffz(i) = abs( positions(i,3) - Eph_eci(i,4) ) ./ positions(i,3) .* 100; % in %
    r_real(i) = sqrt(...
        (positions(i,1))^2 + ...
        (positions(i,2))^2 + ...
        (positions(i,3))^2 ...
        );
    r(i) = sqrt(  Eph_eci(i,2)^2 + Eph_eci(i,3)^2 + Eph_eci(i,4)^2  );
    diff_r(i) = abs(r_real(i)-r(i)) / r_real(i) * 100;
end
% figure(1)

x = 1:length(Eph_eci);
x = x / (365.25);
plot(x,diff_r)
title('与DE441星历表对比的相对误差')
ylabel('误差%')
xlabel('年')
set(gca, 'YScale', 'log')

% figure(2) 
% a1 = plot(x,diffx); m1 = 'diff_x'; hold on
% a2 = plot(x,diffy); m2 = 'diff_y'; hold on
% a3 = plot(x,diffz); m3 = 'diff_z'; 
% legend(m1,m2,m3)
% title('各分量相对DE441的相对误差')
% ylabel('误差%')
% xlabel('年')
% set(gca, 'YScale', 'log')

error_r = 0;
% 计算标准差 （平方差）
for i = 1 : size(r,1)
    error_r = error_r + ( r_real(i) - r(i) )^2;
end
error_r = sqrt( error_r / size(r,1) ); % m
disp('与DE441星历表的标准差 (m):')
disp(error_r)



