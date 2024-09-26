%get the Acceleration



a_earth = zeros( length(Eph_eci),1 );
a_mars = zeros( length(Eph_eci),1 );
a_J2 = zeros( length(Eph_eci),1 );
a_jupiter = zeros( length(Eph_eci),1 );
a_moon = zeros( length(Eph_eci),1 );
a_relative = zeros( length(Eph_eci),1 );
a_saturn = zeros( length(Eph_eci),1 );
a_solarP = zeros( length(Eph_eci),1 );
a_sun = zeros( length(Eph_eci),1 );
a_venus = zeros( length(Eph_eci),1 );
a_yarkovsky = zeros( length(Eph_eci),1 );


for i = 1 : size(Eph_eci,1)
    MJD_UTC = Mjd0_UTC+Eph_eci(i,1)/86400;
    [year,month, day, fd] = iauJd2cal( 2400000.5, MJD_UTC);
    [hour, minute, sec] = fd_to_hms(fd);
    [year, month, day, hour, minute, sec] = fix_seconds(year, month, day, hour, minute, sec);
    sec = round(sec,3);
    month = monthToString(month);

    % Define the start and stop times for the simulation
    str = sprintf(' %s %g , %g %g:%g:%g', month, day, year, hour, minute, sec);
    et= cspice_str2et(str);
    [r_Sun, ~] = cspice_spkezr('Sun', et, 'J2000', 'NONE', 'SUN');
    r_Sun = r_Sun(1:3).*1000; % m/s
    [r_Moon, ~] = cspice_spkezr('Moon', et, 'J2000', 'NONE', 'SUN');
    r_Moon = r_Moon(1:3).*1000; % m/s
    [r_Mercury, ~] = cspice_spkezr('Mercury', et, 'J2000', 'NONE', 'SUN');
    r_Mercury = r_Mercury(1:3).*1000; % m/s
    [r_Venus, ~] = cspice_spkezr('Venus', et, 'J2000', 'NONE', 'SUN');
    r_Venus = r_Venus(1:3).*1000; % m/s
    [r_Mars, ~] = cspice_spkezr('Mars Barycenter', et, 'J2000', 'NONE', 'SUN');
    r_Mars = r_Mars(1:3).*1000; % m/s
    [r_Jupiter, ~] = cspice_spkezr('Jupiter Barycenter', et, 'J2000', 'NONE', 'SUN');
    r_Jupiter = r_Jupiter(1:3).*1000; % m/s
    [r_Saturn, ~] = cspice_spkezr('Saturn Barycenter', et, 'J2000', 'NONE', 'SUN');
    r_Saturn = r_Saturn(1:3).*1000; % m/s
    [r_Uranus, ~] = cspice_spkezr('Uranus Barycenter', et, 'J2000', 'NONE', 'SUN');
    r_Uranus = r_Uranus(1:3).*1000; % m/s
    [r_Neptune, ~] = cspice_spkezr('Neptune Barycenter', et, 'J2000', 'NONE', 'SUN');
    r_Neptune = r_Neptune(1:3).*1000; % m/s
    [r_Pluto, ~] = cspice_spkezr('Pluto Barycenter', et, 'J2000', 'NONE', 'SUN');
    r_Pluto = r_Pluto(1:3).*1000; % m/s
    [r_Earth, ~] = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');
    r_Earth = r_Earth(1:3).*1000; % m/s


    Y(1:3,1) = Eph_eci(i,2:4);
    v(1:3,1) = Eph_eci(i,5:7);
    disp(i)
    a_earth(i,1) =norm( AccelPointMass(Y(1:3),r_Earth,398600441500000));
    a_mars(i,1) = norm(AccelPointMass(Y(1:3),r_Mars,42828375816000));
    a_venus(i,1) = norm(AccelPointMass(Y(1:3),r_Venus,324858592000000));
    a_sun(i,1) = norm(Accel_two_body(Y(1:3),132712440041279422464));
    a_moon(i,1) = norm(AccelPointMass(Y,r_Moon,4902800192171.3935546875));
    a_solarP(i,1) =  norm(AccelsRad(Y,149597870700));
    a_yarkovsky(i,1) = norm(AccelYarkovsky(Y,149597870700,v));
    a_relative(i,1) = norm(Relativity(Y,v));
    a_jupiter(i,1) = norm(AccelPointMass(Y(1:3),r_Jupiter,126712764100000000));
    a_saturn(i,1) = norm(AccelPointMass(Y(1:3),r_Saturn,37940584841800000));
    a_J2(i,1) = norm(AccelJ2(Y(1:3),r_Earth,398600441500000));
end
% plot graph

x = 1 : size(Eph_eci,1);
x = x / 365.25   ;
% a1 = plot(x,a_earth); m1 = 'Earth';hold on
% a2 = plot(x,a_sun); m2 = 'Sun';hold on
% a3 = plot(x,a_moon); m3 = 'Moon';hold on
% a4 = plot(x,a_solarP); m4 = 'SolarP';hold on
% a5 = plot(x, a_yarkovsky); m5 = 'Yarkovsky';hold on
% a6 = plot( x,a_relative); m6 = 'Relative';hold on
% a7 = plot(x,a_venus); m7 = 'Venus';hold on
% a8 = plot(x,a_mars,'k'); m8 = 'Mars';hold on
% a9 = plot(x,a_jupiter); m9 = 'Jupiter';hold on
a10 = plot(x,a_J2,'B'); m10 = 'J2';hold on
% legend(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
set(gca, 'YScale', 'log')
ylabel('加速度m/s^2')
xlabel('year')
