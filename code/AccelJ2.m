function a = AccelJ2(Y,r_Earth,GM)

J2 = 1.08262668e-3;
Re = 6371000; %m
r_ECI = Y- r_Earth;
r = norm(r_ECI);
a_ECI = - GM *...
    (...
    ( 3 * J2) / (r^3 *2) * (Re/r)^2 *...
    [...
    (1 - 5*r_ECI(3)^2/r^2)*r_ECI(1);...
    (1 - 5*r_ECI(3)^2/r^2)*r_ECI(2);...
    (3 - 5*r_ECI(3)^2/r^2)*r_ECI(3)...
    ]...
    );
a = a_ECI;



      
end 