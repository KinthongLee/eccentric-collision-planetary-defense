function a = Accel_two_body(r, GM)

% Relative position vector of satellite w.r.t. point mass 
d = r ;

% Acceleration 
a = -GM * ( d/(norm(d)^3) );
end