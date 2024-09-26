
function a = AccelDrag(dens, r, v, T, Area, mass, CD, Omega)

% Earth angular velocity vector [rad/s]
omega = [0.0, 0.0, Omega]';

% Position and velocity in true-of-date system
r_tod = T * r;
v_tod = T * v;
  
% Velocity relative to the Earth's atmosphere
v_rel = v_tod - cross(omega, r_tod);
v_abs = norm(v_rel);

% Acceleration
a_tod = -0.5*CD*(Area/mass)*dens*v_abs*v_rel;

a = T' * a_tod;

